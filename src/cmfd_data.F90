module cmfd_data

!==============================================================================
! CMFD_DATA -- This module processes the cmfd tally object to generate
! parameters for CMFD calculation.
!==============================================================================

  use, intrinsic :: ISO_FORTRAN_ENV

  use constants

  implicit none
  private
  public :: set_up_cmfd, neutron_balance

contains

!==============================================================================
! SET_UP_CMFD configures cmfd object for a CMFD eigenvalue calculation
!==============================================================================

  subroutine set_up_cmfd()

    use cmfd_header,         only: allocate_cmfd
    use constants,           only: CMFD_NOACCEL
    use global,              only: cmfd, cmfd_coremap, cmfd_downscatter

    ! Check for core map and set it up
    if ((cmfd_coremap) .and. (cmfd%mat_dim == CMFD_NOACCEL)) call set_coremap()

    ! Calculate all cross sections based on reaction rates from last batch
    call compute_xs()

    ! DEBUG: read in Ds, cross-sections when requested
    call read_in_reference_parameters()
    
    ! Compute effective downscatter cross section
    if (cmfd_downscatter) call compute_effective_downscatter()

    ! Check neutron balance
    call neutron_balance()
    
    ! Calculate dtilde
    call compute_dtilde()

    ! Calculate dhat
    call compute_dhat()

  end subroutine set_up_cmfd

!===============================================================================
! COMPUTE_XS takes tallies and computes macroscopic cross sections
!===============================================================================

  subroutine compute_xs()

    use constants,    only: FILTER_MESH, FILTER_ENERGYIN, FILTER_ENERGYOUT,     &
                            FILTER_SURFACE, IN_RIGHT, OUT_RIGHT, IN_FRONT,      &
                            OUT_FRONT, IN_TOP, OUT_TOP, CMFD_NOACCEL, ZERO,     &
                            ONE, TINY_BIT
    use error,        only: fatal_error
    use global,       only: cmfd, n_cmfd_tallies, cmfd_tallies, meshes,&
                            matching_bins, keff
    use mesh,         only: mesh_indices_to_bin
    use mesh_header,  only: StructuredMesh
    use string,       only: to_str
    use tally_header, only: TallyObject
    use output,       only: write_message

    integer :: nx            ! number of mesh cells in x direction
    integer :: ny            ! number of mesh cells in y direction
    integer :: nz            ! number of mesh cells in z direction
    integer :: ng            ! number of energy groups
    integer :: i             ! iteration counter for x
    integer :: j             ! iteration counter for y
    integer :: k             ! iteration counter for z
    integer :: g             ! iteration counter for g
    integer :: h             ! iteration counter for outgoing groups
    integer :: q             ! iteration counter for quadrature (1,2,3,4)
    integer :: nq            ! number of quadratures: 2 for CMFD, 4 for LOO
    integer :: ital          ! tally object index
    integer :: ijk(3)        ! indices for mesh cell
    integer :: score_index   ! index to pull from tally object
    integer :: i_mesh        ! index in meshes array
    integer :: i_filter_mesh ! index for mesh filter
    integer :: i_filter_ein  ! index for incoming energy filter
    integer :: i_filter_eout ! index for outgoing energy filter
    integer :: i_filter_surf ! index for surface filter
    integer :: cnt           ! counter for non-zero entripes in openmc_src
    integer :: n             ! number of realizations
    real(8) :: flux          ! temp variable for flux
    type(TallyObject),    pointer :: t => null() ! pointer for tally object
    type(StructuredMesh), pointer :: m => null() ! pointer for mesh object
    ! for debugging
    integer :: ii
    integer :: jj
    integer :: kk

    ! Extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! Associate tallies and mesh
    t => cmfd_tallies(1)
    i_mesh = t % filters(t % find_filter(FILTER_MESH)) % int_bins(1)
    m => meshes(i_mesh)

    ! Set mesh widths
    cmfd % hxyz(1,:,:,:) = m % width(1) ! set x width
    cmfd % hxyz(2,:,:,:) = m % width(2) ! set y width
    cmfd % hxyz(3,:,:,:) = m % width(3) ! set z width

    ! Save a copy of the fission source from the previous batch,
    ! because LOO needs a copy of the fission source before this batch
    ! of MC is performed.
    if (sum(cmfd % openmc_src) > ZERO) then
       do k = 1, nz
          do j = 1, ny
             do i = 1, nx
                do h = 1, ng
                   cmfd % openmc_src_old(h, i, j, k) = &
                          cmfd % openmc_src(h, i, j, k)
                end do
             end do
          end do
       end do
    else 
       cmfd % openmc_src_old = ONE
    end if

    ! reset parameters before computation
    cmfd % openmc_src = ZERO
    cmfd % openmc_total_src = ZERO
    cmfd % flux = ZERO
    cmfd % keff_bal = ZERO
    cmfd % totalxs = ZERO
    cmfd % p1scattxs = ZERO
    cmfd % diffcof = ZERO
    cmfd % scattxs = ZERO
    cmfd % nfissxs = ZERO
    cmfd % current = ZERO
    cmfd % quad_current = ZERO

   ! Begin loop around tallies
   TAL: do ital = 1, n_cmfd_tallies

     ! Associate tallies and mesh
     t => cmfd_tallies(ital)
     i_mesh = t % filters(t % find_filter(FILTER_MESH)) % int_bins(1)
     m => meshes(i_mesh)

     ! debug: pull out number of realizations
     !n = dble(t % n_realizations)
     n = 1

     i_filter_mesh = t % find_filter(FILTER_MESH)
     i_filter_ein  = t % find_filter(FILTER_ENERGYIN)
     i_filter_eout = t % find_filter(FILTER_ENERGYOUT)
     i_filter_surf = t % find_filter(FILTER_SURFACE)

     ! Begin loop around space
     ZLOOP: do k = 1,nz

       YLOOP: do j = 1,ny

          XLOOP: do i = 1,nx

            ! Check for active mesh cell
            if (allocated(cmfd%coremap)) then
              if (cmfd%coremap(i,j,k) == CMFD_NOACCEL) then
                cycle
              end if
            end if

            ! Loop around energy groups
            OUTGROUP: do h = 1,ng

              ! Start tally 1: reaction types that involve a single energy level: 
              ! flux, total reaction, p1 scattering.
              TALLY: if (ital == 1) then

                ! Reset all bins to 1
                matching_bins(1:t%n_filters) = 1

                ! Set ijk as mesh indices
                ijk = (/ i, j, k /)

                ! Get bin number for mesh indices
                matching_bins(i_filter_mesh) = mesh_indices_to_bin(m,ijk)

                ! Apply energy in filter
                if (i_filter_ein > 0) then
                  matching_bins(i_filter_ein) = ng - h + 1
                end if

                ! Calculate score index from bins
                score_index = sum((matching_bins(1:t%n_filters) - 1) * t%stride) + 1

                ! Get flux
                flux = t % results(1,score_index) % sum / n
                cmfd % flux(h,i,j,k) = flux 

                ! Detect zero flux, abort if located
                if ((flux - ZERO) < TINY_BIT) then
                  do ii= 1,nx
                    do jj = 1,ny
                      do kk = 1,nz
                         call write_message("flux at (" // to_str(ii) //&
                             &to_str(jj) // to_str(kk) // &
                             &" 1) group 1 = "//&
                             &to_str(cmfd % flux(ii, jj, kk, 1)), 5)
                      enddo
                    enddo
                  enddo
                  call fatal_error('Detected zero flux without coremap overlay &
                       &at: (' // to_str(i) // ',' // to_str(j) // ',' // &
                       &to_str(k) // ') in group ' // to_str(h))
                end if

                ! Get total rr and convert to total xs
                cmfd % totalxs(h,i,j,k) = t % results(2,score_index) % sum / flux / n

                ! Get p1 scatter rr and convert to p1 scatter xs
                cmfd % p1scattxs(h,i,j,k) = t % results(3,score_index) % sum / flux / n

                ! Calculate diffusion coefficient
                cmfd % diffcof(h,i,j,k) = ONE/(3.0_8*(cmfd % totalxs(h,i,j,k) - &
                     cmfd % p1scattxs(h,i,j,k)))

              else if (ital == 2) then

                ! Begin loop to get energy out tallies, that is,
                ! events (scattering and fission) that leads to a
                ! neutron changing energy from group h to g.
                INGROUP: do g = 1, ng
                  ! Reset all bins to 1
                  matching_bins(1:t%n_filters) = 1

                  ! Set ijk as mesh indices
                  ijk = (/ i, j, k /)

                  ! Get bin number for mesh indices
                  matching_bins(i_filter_mesh) = mesh_indices_to_bin(m,ijk)

                  if (i_filter_ein > 0) then
                    ! Apply energy in filter
                    matching_bins(i_filter_ein) = ng - h + 1

                    ! Set energy out bin
                    matching_bins(i_filter_eout) = ng - g + 1
                  end if

                  ! Calculate score index from bins
                  score_index = sum((matching_bins(1:t%n_filters) - 1) * t%stride) + 1

                  ! Get scattering
                  cmfd % scattxs(h,g,i,j,k) = t % results(1,score_index) % sum /&
                       cmfd % flux(h,i,j,k) / n

                  ! Get nu-fission
                  cmfd % nfissxs(h,g,i,j,k) = t % results(2,score_index) % sum /&
                       cmfd % flux(h,i,j,k) / n
 
                  ! Bank fission source: any fission event from group
                  ! h (outgroup) to g (ingroup) is added to group g's
                  ! fission source counter openmc_src
                  cmfd % openmc_src(g,i,j,k) = cmfd % openmc_src(g,i,j,k) + &
                       t % results(2,score_index) % sum / n

                  cmfd % keff_bal = cmfd % keff_bal + &
                       t % results(2,score_index) % sum / &
                       dble(t % n_realizations)

                  ! Bank total source (fission + scattering) for LOO:
                  ! any scattering event (1) or fission event (2) from
                  ! group h (outgroup) to g (ingroup) is added to
                  ! group g's total source counter openmc_total_src
                  cmfd % openmc_total_src(g,i,j,k) = &
                       cmfd % openmc_total_src(g,i,j,k) + &
                       t % results(1,score_index) % sum + &
                       t % results(2,score_index) % sum / keff
                end do INGROUP

             else if (ital == 3) then

                ! Set number of quadrature to be 2 for CMFD, 4 for LOO
                 nq = 2

                ! Initialize and filter for energy
                matching_bins(1:t%n_filters) = 1
                if (i_filter_ein > 0) then
                  matching_bins(i_filter_ein) = ng - h + 1
                end if

                ! Left surface: from left neighbor
                matching_bins(i_filter_mesh) = mesh_indices_to_bin(m, &
                     (/ i-1, j, k /) + 1, .true.)
                do q = 1, nq
                   matching_bins(i_filter_surf) = q
                   score_index = sum((matching_bins(1:t%n_filters) - 1) &
                        * t % stride) + 1
                   cmfd % current(q,h,i,j,k) = &
                        t % results(1,score_index) % sum / n
                end do

                ! Right surface
                matching_bins(i_filter_mesh) = mesh_indices_to_bin(m, &
                     (/ i, j, k /) + 1, .true.)
                do q = 1, nq
                   matching_bins(i_filter_surf) = q
                   score_index = sum((matching_bins(1:t%n_filters) - 1) &
                        * t % stride) + 1
                   cmfd % current(q+nq,h,i,j,k) = &
                        t % results(1,score_index) % sum / n
                end do

                ! Back surface
                matching_bins(i_filter_mesh) = mesh_indices_to_bin(m, &
                     (/ i, j-1, k /) + 1, .true.)
                do q = nq + 1, 2 * nq
                   matching_bins(i_filter_surf) = q
                   score_index = sum((matching_bins(1:t%n_filters) - 1) &
                        * t % stride) + 1
                   cmfd % current(q+nq,h,i,j,k) = &
                        t % results(1,score_index) % sum / n
                end do

                ! Front surface
                matching_bins(i_filter_mesh) = mesh_indices_to_bin(m, &
                     (/ i, j, k /) + 1, .true.)
                do q = nq + 1, 2 * nq
                   matching_bins(i_filter_surf) = q
                   score_index = sum((matching_bins(1:t%n_filters) - 1) &
                        * t % stride) + 1
                   cmfd % current(q+2*nq,h,i,j,k) = &
                        t % results(1,score_index) % sum / n
                end do

                ! Bottom surface
                matching_bins(i_filter_mesh) = mesh_indices_to_bin(m, &
                     (/ i, j, k-1 /) + 1, .true.)
                do q = 2 * nq + 1, 3 * nq
                   matching_bins(i_filter_surf) = q
                   score_index = sum((matching_bins(1:t%n_filters) - 1) &
                        * t % stride) + 1
                   cmfd % current(q+2*nq,h,i,j,k) = &
                        t % results(1,score_index) % sum / n
                end do

                ! Top surface
                matching_bins(i_filter_mesh) = mesh_indices_to_bin(m, &
                     (/ i, j, k /) + 1, .true.)
                do q = 2 * nq + 1, 3 * nq
                   matching_bins(i_filter_surf) = q
                   score_index = sum((matching_bins(1:t%n_filters) - 1) &
                        * t % stride) + 1
                   cmfd % current(q+3*nq,h,i,j,k) = &
                        t % results(1,score_index) % sum / n
                end do

                ! Set number of quadrature to be 4 for LOO quad currents
                 nq = 4

                ! Initialize and filter for energy
                matching_bins(1:t%n_filters) = 1
                if (i_filter_ein > 0) then
                  matching_bins(i_filter_ein) = ng - h + 1
                end if

                ! Left surface: from left neighbor
                matching_bins(i_filter_mesh) = mesh_indices_to_bin(m, &
                     (/ i-1, j, k /) + 1, .true.)
                do q = 1, nq
                   matching_bins(i_filter_surf) = 6 + q
                   score_index = sum((matching_bins(1:t%n_filters) - 1) &
                        * t % stride) + 1
                   cmfd % quad_current(q,h,i,j,k) = &
                        t % results(1,score_index) % sum / n
                end do

                ! Right surface
                matching_bins(i_filter_mesh) = mesh_indices_to_bin(m, &
                     (/ i, j, k /) + 1, .true.)
                do q = 1, nq
                   matching_bins(i_filter_surf) = 6 + q
                   score_index = sum((matching_bins(1:t%n_filters) - 1) &
                        * t % stride) + 1
                   cmfd % quad_current(q+nq,h,i,j,k) = &
                        t % results(1,score_index) % sum / n
                end do

                ! Back surface
                matching_bins(i_filter_mesh) = mesh_indices_to_bin(m, &
                     (/ i, j-1, k /) + 1, .true.)
                do q = nq + 1, 2 * nq
                   matching_bins(i_filter_surf) = 6 + q
                   score_index = sum((matching_bins(1:t%n_filters) - 1) &
                        * t % stride) + 1
                   cmfd % quad_current(q+nq,h,i,j,k) = &
                        t % results(1,score_index) % sum / n
                end do 

                ! Front surface
                matching_bins(i_filter_mesh) = mesh_indices_to_bin(m, &
                     (/ i, j, k /) + 1, .true.)
                do q = nq + 1, 2 * nq
                   matching_bins(i_filter_surf) = 6 + q
                   score_index = sum((matching_bins(1:t%n_filters) - 1) &
                        * t % stride) + 1
                   cmfd % quad_current(q+2*nq,h,i,j,k) = &
                        t % results(1,score_index) % sum / n
                end do

                ! FIXME: add top & bottom for 3D implementation
              end if TALLY

            end do OUTGROUP

          end do XLOOP

        end do YLOOP

      end do ZLOOP

    end do TAL

    ! Normalize openmc source distribution such that the average is 1.0
    cnt = ZERO

    do k = 1,nz

       do j = 1,ny

          do i = 1,nx

            do h = 1,ng
            
                if (cmfd % openmc_src(h, i, j, k) > ZERO) cnt = cnt + 1

            end do

          end do

       end do 

    end do

    cmfd % openmc_src = cmfd % openmc_src / sum(cmfd % openmc_src) * cnt

    ! Nullify all pointers
    if (associated(t)) nullify(t)
    if (associated(m)) nullify(m)

  end subroutine compute_xs

!===============================================================================
! READ_IN_REFERENCE_PARAMETERS overrite the cross sections by compute_xs()
!===============================================================================

  subroutine read_in_reference_parameters()

    use constants,    only: FILTER_MESH, FILTER_ENERGYIN, FILTER_ENERGYOUT,     &
         FILTER_SURFACE, IN_RIGHT, OUT_RIGHT, IN_FRONT,      &
         OUT_FRONT, IN_TOP, OUT_TOP, CMFD_NOACCEL, ZERO,     &
         ONE, TINY_BIT
    use error,        only: fatal_error
    use global,       only: cmfd, k_generation, overall_gen
    use string,       only: to_str
    use output,       only: write_message

    integer :: nx            ! number of mesh cells in x direction
    integer :: ny            ! number of mesh cells in y direction
    integer :: nz            ! number of mesh cells in z direction
    integer :: ng            ! number of energy groups
    integer :: i             ! iteration counter for x
    integer :: j             ! iteration counter for y
    integer :: k             ! iteration counter for z
    integer :: h             ! iteration counter for energy group
    integer :: g             ! iteration counter for energy group       
    integer :: s             ! iteration counter for surfaces
    integer :: index, index2

    !real, dimension(60) :: current
    !real, dimension(80) :: quad_current
    !real(8), dimension(5) :: openmc_src_old
    real(8), allocatable :: totalxs(:), d(:), scattxs(:), nfissxs(:)
    real(8), allocatable :: flux(:)

    ! Extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    allocate(totalxs(ng))
    allocate(d(ng))
    allocate(scattxs(ng*ng))
    allocate(nfissxs(ng*ng))
    allocate(flux(nx*ng))

    !k_generation(overall_gen) =  1.6010261

    ! note: openmc_src_old, flux are volume-integrated, current, quad_currents are area-integrated
    if (ng == 1) then 
       totalxs = (/0.927057/)
       scattxs = (/0.8996023/)
       nfissxs = (/0.0371365/)
       d = (/1.19/)
       flux(:) = 1.0 
    else if (ng == 2) then
       totalxs = (/0.64336155, 1.94585169/)
       scattxs = (/0.608241, 0.02568, 0.0, 1.853706/)
       nfissxs = (/0.005233, 0.0, 0.151707, 0.0/)
       d = (/1.44, 0.284/)
       flux = (/ 3.92672896,   1.02983295,   8.8257216 ,   2.45199824,&
        13.40635594,   3.73209005,  17.66224929,   4.91811123,&
        21.5082574 ,   5.98930568,  24.85776969,   6.92183052,&
        27.6353359 ,   7.69520216,  29.77977404,   8.29268241,&
        31.23987571,   8.6992619 ,  31.97863565,   8.90487998,&
        31.98218281,   8.90611406,  31.24706938,   8.70134445,&
        29.79362647,   8.29655209,  27.65130049,   7.70007488,&
        24.87142261,   6.92575664,  21.52007272,   5.99232399,&
        17.67133909,   4.92064943,  13.41131647,   3.73341871,&
         8.8293219 ,   2.4530031 ,   3.92759676,   1.0300203 /)
    else 
       return
    end if

    ! Begin loop around space
    ZLOOP: do k = 1,nz

       YLOOP: do j = 1,ny

          XLOOP: do i = 1,nx

             ! Check for active mesh cell
             if (allocated(cmfd%coremap)) then
                if (cmfd%coremap(i,j,k) == CMFD_NOACCEL) then
                   cycle
                end if
             end if

             ! Loop around energy groups
             OUTGROUP: do h = 1,ng

                index = h;
                cmfd % diffcof(h,i,j,k) = d(index);

                ! the rest of the cross-sections:
                !cmfd % totalxs(h, i, j, k) = totalxs(index);
                !do g = 1,ng
                !   index2 = (h - 1) * ng + g;                
                !   cmfd % scattxs(h, g, i, j, k) = scattxs(index2);
                !   cmfd % nfissxs(h, g, i, j, k) = nfissxs(index2);
                !end do

                !cmfd % flux(h, i, j, k) = flux( (i - 1) * ng + h);
                !cmfd % openmc_src_old(h, i, j, k) = openmc_src_old(index);
                !do s = 1, 12
                !   index2 = (index - 1) * 12 + s;
                !   cmfd % current(s, h, i, j, k) = current(index2);
                !end do
                !do s = 1, 16
                !   index2 = (index - 1) * 16 + s;
                !   cmfd % quad_current(s, h, i, j, k) = quad_current(index2)
                !end do
             end do OUTGROUP

          end do XLOOP

       end do YLOOP

    end do ZLOOP

    cmfd % openmc_src_old = cmfd % openmc_src_old / sum (cmfd % openmc_src_old) * (nx * ny * nz)

  end subroutine read_in_reference_parameters

!===============================================================================
! SET_COREMAP is a routine that sets the core mapping information
!===============================================================================

  subroutine set_coremap()

    use constants,  only: CMFD_NOACCEL
    use global,     only: cmfd

    integer :: counter=1 ! counter for unique fuel assemblies
    integer :: nx        ! number of mesh cells in x direction
    integer :: ny        ! number of mesh cells in y direction
    integer :: nz        ! number of mesh cells in z direction
    integer :: i         ! iteration counter for x
    integer :: j         ! iteration counter for y
    integer :: k         ! iteration counter for z

    ! Extract spatial indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)

    ! Count how many fuel assemblies exist
    cmfd % mat_dim = sum(cmfd % coremap - 1)

    ! Allocate indexmap
    if (.not. allocated(cmfd % indexmap)) &
         allocate(cmfd % indexmap(cmfd % mat_dim,3))

    ! Begin loops over spatial indices
    ZLOOP: do k = 1, nz

      YLOOP: do j = 1, ny

        XLOOP: do i = 1, nx

          ! Check for reflector
          if (cmfd % coremap(i,j,k) == 1) then

            ! reset value to CMFD no acceleration constant
            cmfd % coremap(i,j,k) = CMFD_NOACCEL

          else

            ! Must be a fuel --> give unique id number
            cmfd % coremap(i,j,k) = counter
            cmfd % indexmap(counter,1) = i
            cmfd % indexmap(counter,2) = j
            cmfd % indexmap(counter,3) = k
            counter = counter + 1

          end if

        end do XLOOP

      end do YLOOP

    end do ZLOOP

  end subroutine set_coremap

!===============================================================================
! NEUTRON_BALANCE computes the RMS neutron balance over the CMFD mesh
!===============================================================================

  subroutine neutron_balance()

    use constants,    only: ONE, ZERO, CMFD_NOACCEL, CMFD_NORES
    use global,       only: cmfd, keff, k_generation, current_batch, &
         cmfd_rebalance, overall_gen

    integer :: nx           ! number of mesh cells in x direction
    integer :: ny           ! number of mesh cells in y direction
    integer :: nz           ! number of mesh cells in z direction
    integer :: ng           ! number of energy groups
    integer :: i            ! iteration counter for x
    integer :: j            ! iteration counter for y
    integer :: k            ! iteration counter for z
    integer :: g            ! iteration counter for g
    integer :: h            ! iteration counter for outgoing groups
    integer :: l            ! iteration counter for leakage
    integer :: cnt          ! number of locations to count neutron balance
    real(8), allocatable :: leakage(:)    ! leakage term in neutron balance
    real(8) :: interactions ! total number of interactions in balance
    real(8) :: scattering   ! scattering term in neutron balance
    real(8), allocatable :: fission(:)      ! fission term in neutron balance
    real(8) :: construction, destruction
    real(8) :: res          ! residual of neutron balance for whole geometry
    real(8) :: rms          ! RMS of the residual
    real(8) :: flux1             ! group 1 volume int flux
    real(8) :: flux2             ! group 2 volume int flux
    real(8) :: sigt1             ! group 1 total xs
    real(8) :: sigt2_old         ! group 2 total xs
    real(8) :: sigt2             ! group 2 total xs
    real(8) :: sigs11            ! scattering transfer 1 --> 1
    real(8) :: sigs21            ! scattering transfer 2 --> 1
    real(8) :: sigs12            ! scattering transfer 1 --> 2
    real(8) :: sigs22            ! scattering transfer 2 --> 2
    real(8) :: siga1             ! group 1 abs xs
    real(8) :: siga2             ! group 2 abs xs

    ! Extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    if (.not. allocated(leakage)) allocate(leakage(ng))
    if (.not. allocated(fission)) allocate(fission(ng))

    ! Allocate res dataspace
    if (.not. allocated(cmfd % resnb)) allocate(cmfd % resnb(ng,nx,ny,nz))

    ! Reset rms and cnt
    rms = ZERO
    cnt = 0

    construction = 0
    destruction = 0
    
    ! Begin loop around space and energy groups
    ZLOOP: do k = 1, nz

      YLOOP: do j = 1, ny

        XLOOP: do i = 1, nx

           leakage = ZERO
           fission = ZERO
           
           GROUPG: do g = 1, ng

            ! Check for active mesh
            if (allocated(cmfd % coremap)) then
              if (cmfd % coremap(i,j,k) == CMFD_NOACCEL) then
                cmfd % resnb(g,i,j,k) = CMFD_NORES
                cycle
              end if
            end if

            ! Get leakage
            LEAK: do l = 1, 3

              leakage(g) = leakage(g) + ((cmfd % current(4*l,g,i,j,k) - &
                   cmfd % current(4*l-1,g,i,j,k))) - &
                   ((cmfd % current(4*l-2,g,i,j,k) - &
                   cmfd % current(4*l-3,g,i,j,k)))

            end do LEAK

            GROUPH: do h = 1, ng

              fission(g) = fission(g) + cmfd % nfissxs(h,g,i,j,k) * &
                   cmfd % flux(h,i,j,k)

            end do GROUPH

          end do GROUPG

          if ((cmfd_rebalance) .and. (ng == 2)) then
             flux1 = cmfd % flux(1,i,j,k)
             flux2 = cmfd % flux(2,i,j,k)
             sigt1 = cmfd % totalxs(1,i,j,k)
             sigt2_old = cmfd % totalxs(2,i,j,k)
             sigs11 = cmfd % scattxs(1,1,i,j,k)
             sigs21 = cmfd % scattxs(2,1,i,j,k)
             sigs12 = cmfd % scattxs(1,2,i,j,k)
             sigs22 = cmfd % scattxs(2,2,i,j,k)
             siga1 = sigt1 - sigs11 - sigs12
             siga2 = sigt2 - sigs22 - sigs21

             sigs12 = ((ONE / keff) * fission(1) + sigs21 * flux2 - leakage(1) &
                  & - siga1 * flux1) / flux1
             cmfd % scattxs(1,2,i,j,k) = sigs12
             sigt1 = siga1 + sigs11 + sigs12
             cmfd % totalxs(1, i, j, k) = sigt1
             
             siga2 = (sigs12 * flux1 + (ONE / keff) * fission(2) - leakage(2) &
                  & - sigs21 * flux2 ) / flux2
             sigt2 =  siga2 + sigs22 + sigs21
             cmfd % totalxs(2, i, j, k) = sigt2
          end if

          GROUPG2: do g = 1, ng

             ! Interactions
             interactions = cmfd % totalxs(g,i,j,k) * cmfd % flux(g,i,j,k)

             ! Get scattering and fission
             scattering = ZERO
             GROUPH3: do h = 1, ng

                scattering = scattering + cmfd % scattxs(h,g,i,j,k) * &
                     cmfd % flux(h,i,j,k)
             end do GROUPH3
             
             ! Compute residual
             res = leakage(g) + interactions - scattering - &
                  & (ONE /   k_generation(overall_gen)) * fission(g)

             ! DEBUG
             if (.false.) then
                write(OUTPUT_UNIT, FMT='(" res_bef=",ES14.6,ES14.6,ES14.6,ES14.6,ES14.6)') &
                     res, leakage(g), interactions- scattering, fission(g), &
                     k_generation(overall_gen)
             end if

             ! Normalize by flux
             res = res / cmfd % flux(g,i,j,k)
             
             ! Bank res in cmfd object
             cmfd % resnb(g,i,j,k) = res
             
             ! Take square for RMS calculation
             rms = rms + res ** 2
             cnt = cnt + 1
             ! write(OUTPUT_UNIT, FMT='("res_aft =", ES9.2  )') res

             ! accumulate into total counters
             destruction = destruction + leakage(g) + interactions - scattering
             construction = construction + fission(g)

          end do GROUPG2

        end do XLOOP

      end do YLOOP

    end do ZLOOP

    ! DEBUG: print out k constructed from given NDA parameters
    !write(OUTPUT_UNIT, FMT='(" NDA parameters produce k = ",ES14.6)') &
    !     construction / destruction

    ! Calculate RMS and record in vector for this batch
    cmfd % balance(current_batch) = sqrt(ONE/dble(cnt)*rms)

  end subroutine neutron_balance

!===============================================================================
! COMPUTE_DTILDE precomputes the diffusion coupling coefficient
!===============================================================================

  subroutine compute_dtilde()

    use constants,  only: CMFD_NOACCEL, ZERO_FLUX, TINY_BIT
    use global,     only: cmfd, cmfd_coremap

    integer :: nx           ! maximum number of cells in x direction
    integer :: ny           ! maximum number of cells in y direction
    integer :: nz           ! maximum number of cells in z direction
    integer :: ng           ! maximum number of energy groups
    integer :: nxyz(3,2)    ! single vector containing boundary locations
    integer :: i            ! iteration counter for x
    integer :: j            ! iteration counter for y
    integer :: k            ! iteration counter for z
    integer :: g            ! iteration counter for groups
    integer :: l            ! iteration counter for leakages
    integer :: xyz_idx      ! index for determining if x,y or z leakage
    integer :: dir_idx      ! index for determining - or + face of cell
    integer :: shift_idx    ! parameter to shift index by +1 or -1
    integer :: neig_idx(3)  ! spatial indices of neighbour
    integer :: bound(6)     ! vector containing indices for boudary check
    real(8) :: albedo(6)    ! albedo vector with global boundaries
    real(8) :: cell_dc      ! diffusion coef of current cell
    real(8) :: cell_hxyz(3) ! cell dimensions of current ijk cell
    real(8) :: neig_dc      ! diffusion coefficient of neighbor cell
    real(8) :: neig_hxyz(3) ! cell dimensions of neighbor cell
    real(8) :: dtilde       ! finite difference coupling parameter
    real(8) :: ref_albedo   ! albedo to reflector

    ! Get maximum of spatial and group indices
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! Create single vector of these indices for boundary calculation
    nxyz(1,:) = (/1,nx/)
    nxyz(2,:) = (/1,ny/)
    nxyz(3,:) = (/1,nz/)

    ! Get boundary condition information
    albedo = cmfd%albedo

    ! Loop over group and spatial indices
    ZLOOP: do k = 1, nz

      YLOOP: do j = 1, ny

        XLOOP: do i = 1, nx

          GROUP: do g = 1, ng

            ! Check for active mesh cell
            if (allocated(cmfd%coremap)) then
              if (cmfd%coremap(i,j,k) == CMFD_NOACCEL) cycle
            end if

            ! Get cell data
            cell_dc = cmfd % diffcof(g,i,j,k)
            cell_hxyz = cmfd % hxyz(:,i,j,k)

            ! Setup of vector to identify boundary conditions
            bound = (/i,i,j,j,k,k/)

            ! Begin loop around sides of cell for leakage
            LEAK: do l = 1, 6

              ! Define xyz and +/- indices
              xyz_idx = int(ceiling(real(l)/real(2)))  ! x=1, y=2, z=3
              dir_idx = 2 - mod(l,2) ! -=1, +=2
              shift_idx = -2*mod(l,2) + 1 ! shift neig by -1 or +1

              ! Check if at a boundary
              if (bound(l) == nxyz(xyz_idx,dir_idx)) then

                ! Compute dtilde with albedo boundary condition
                dtilde = (2*cell_dc*(1-albedo(l)))/(4*cell_dc*(1+albedo(l)) + &
                     (1-albedo(l))*cell_hxyz(xyz_idx))

                ! Check for zero flux
                if (abs(albedo(l) - ZERO_FLUX) < TINY_BIT) dtilde = 2*cell_dc / &
                     cell_hxyz(xyz_idx)

              else  ! not a boundary

                ! Compute neighboring cell indices
                neig_idx = (/i,j,k/)                ! begin with i,j,k
                neig_idx(xyz_idx) = shift_idx + neig_idx(xyz_idx)

                ! Get neigbor cell data
                neig_dc = cmfd%diffcof(g,neig_idx(1),neig_idx(2),neig_idx(3))
                neig_hxyz = cmfd%hxyz(:,neig_idx(1),neig_idx(2),neig_idx(3))

                ! Check for fuel-reflector interface
                if (cmfd_coremap) then

                  if (cmfd % coremap(neig_idx(1),neig_idx(2),neig_idx(3)) == &
                       CMFD_NOACCEL .and. cmfd % coremap(i,j,k) /= CMFD_NOACCEL) then

                    ! Get albedo
                    ref_albedo = get_reflector_albedo(l,g,i,j,k)

                    ! Compute dtilde
                    dtilde = (2*cell_dc*(1-ref_albedo))/(4*cell_dc*(1+ &
                         ref_albedo)+(1-ref_albedo)*cell_hxyz(xyz_idx))

                  else ! Not next to a reflector or no core map

                    ! Compute dtilde to neighbor cell
                    dtilde = (2*cell_dc*neig_dc)/(neig_hxyz(xyz_idx)*cell_dc + &
                         cell_hxyz(xyz_idx)*neig_dc)

                  end if

                else ! no core map

                  ! Compute dtilde to neighbor cell
                  dtilde = (2*cell_dc*neig_dc)/(neig_hxyz(xyz_idx)*cell_dc + &
                       cell_hxyz(xyz_idx)*neig_dc)

               end if

              end if

              ! Record dtilde in cmfd object
              cmfd%dtilde(l,g,i,j,k) = dtilde

            end do LEAK

          end do GROUP

        end do XLOOP

      end do YLOOP

    end do ZLOOP

  end subroutine compute_dtilde

!===============================================================================
! COMPUTE_DHAT computes the nonlinear coupling coefficient
!===============================================================================

  subroutine compute_dhat()

    use constants,  only: CMFD_NOACCEL, ZERO
    use global,     only: cmfd, cmfd_coremap, dhat_reset
    use output,     only: write_message
    use string,     only: to_str

    integer :: nx             ! maximum number of cells in x direction
    integer :: ny             ! maximum number of cells in y direction
    integer :: nz             ! maximum number of cells in z direction
    integer :: ng             ! maximum number of energy groups
    integer :: nxyz(3,2)      ! single vector containing boundary locations
    integer :: i              ! iteration counter for x
    integer :: j              ! iteration counter for y
    integer :: k              ! iteration counter for z
    integer :: g              ! iteration counter for groups
    integer :: l              ! iteration counter for leakages
    integer :: xyz_idx        ! index for determining if x,y or z leakage
    integer :: dir_idx        ! index for determining - or + face of cell
    integer :: shift_idx      ! parameter to shift index by +1 or -1
    integer :: neig_idx(3)    ! spatial indices of neighbour
    integer :: bound(6)       ! vector containing indices for boudary check
    real(8) :: cell_dtilde(6) ! cell dtilde for each face
    real(8) :: cell_flux      ! flux in current cell
    real(8) :: current(12)    ! area integrated cell current at each face
    real(8) :: net_current    ! net current on a face
    real(8) :: neig_flux      ! flux in neighbor cell
    real(8) :: dhat           ! dhat equivalence parameter

    ! Get maximum of spatial and group indices
    nx = cmfd%indices(1)
    ny = cmfd%indices(2)
    nz = cmfd%indices(3)
    ng = cmfd%indices(4)

    ! Create single vector of these indices for boundary calculation
    nxyz(1,:) = (/1,nx/)
    nxyz(2,:) = (/1,ny/)
    nxyz(3,:) = (/1,nz/)

    ! Geting loop over group and spatial indices
    ZLOOP:  do k = 1,nz

      YLOOP: do j = 1,ny

        XLOOP: do i = 1,nx

          GROUP: do g = 1,ng

            ! Check for active mesh cell
            if (allocated(cmfd%coremap)) then
              if (cmfd%coremap(i,j,k) == CMFD_NOACCEL) then
                cycle
              end if
            end if

            ! Get cell data
            cell_dtilde = cmfd%dtilde(:,g,i,j,k)
            cell_flux = cmfd%flux(g,i,j,k)/product(cmfd%hxyz(:,i,j,k))
            current = cmfd % current(:,g,i,j,k)

            ! Setup of vector to identify boundary conditions
            bound = (/i,i,j,j,k,k/)

            ! Begin loop around sides of cell for leakage
            LEAK: do l = 1,6

              ! Define xyz and +/- indices
              xyz_idx = int(ceiling(real(l)/real(2)))  ! x=1, y=2, z=3
              dir_idx = 2 - mod(l,2) ! -=1, +=2
              shift_idx = -2*mod(l,2) +1          ! shift neig by -1 or +1

              ! Calculate net current on l face (divided by surf area)
              net_current = (current(2*l) - current(2*l-1)) / &
                   product(cmfd%hxyz(:,i,j,k)) * cmfd%hxyz(xyz_idx,i,j,k)

              ! Check if at a boundary
              if (bound(l) == nxyz(xyz_idx,dir_idx)) then

                ! Compute dhat
                dhat = (net_current - shift_idx*cell_dtilde(l)*cell_flux) / &
                     cell_flux

              else  ! not a boundary

                ! Compute neighboring cell indices
                neig_idx = (/i,j,k/)                ! begin with i,j,k
                neig_idx(xyz_idx) = shift_idx + neig_idx(xyz_idx)

                ! Get neigbor flux
                neig_flux = cmfd%flux(g,neig_idx(1),neig_idx(2),neig_idx(3)) / &
                     product(cmfd%hxyz(:,neig_idx(1),neig_idx(2),neig_idx(3)))

                ! Check for fuel-reflector interface
                if (cmfd_coremap) then

                  if (cmfd % coremap(neig_idx(1),neig_idx(2),neig_idx(3)) == &
                       CMFD_NOACCEL .and. cmfd % coremap(i,j,k) /= CMFD_NOACCEL) then

                    ! compute dhat
                    dhat = (net_current - shift_idx*cell_dtilde(l)*cell_flux) /&
                         cell_flux

                  else ! not a fuel-reflector interface

                    ! Compute dhat
                    dhat = (net_current + shift_idx*cell_dtilde(l)* &
                         (neig_flux - cell_flux))/(neig_flux + cell_flux)

                  end if

                else ! not for fuel-reflector case

                  ! Compute dhat
                  dhat = (net_current + shift_idx*cell_dtilde(l)* &
                       (neig_flux - cell_flux))/(neig_flux + cell_flux)

                end if

              end if

              ! record dhat in cmfd object
              cmfd%dhat(l,g,i,j,k) = dhat

              ! check for dhat reset
              if (dhat_reset) then
                cmfd%dhat(l,g,i,j,k) = ZERO
              end if

            end do LEAK

          end do GROUP

        end do XLOOP

      end do YLOOP

    end do ZLOOP

    ! write that dhats are zero
    if (dhat_reset) then
      call write_message('Dhats reset to zero.', 1)
    end if

  end subroutine compute_dhat

!===============================================================================
! GET_REFLECTOR_ALBEDO is a function that calculates the albedo to the reflector
!===============================================================================

  function get_reflector_albedo(l, g, i, j, k)

    use constants,  only: ONE
    use global,     only: cmfd

    real(8) :: get_reflector_albedo ! reflector albedo
    integer, intent(in) :: i ! iteration counter for x
    integer, intent(in) :: j ! iteration counter for y
    integer, intent(in) :: k ! iteration counter for z
    integer, intent(in) :: g ! iteration counter for groups
    integer, intent(in) :: l ! iteration counter for leakages

    integer :: shift_idx   ! parameter to shift index by +1 or -1
    real(8) :: current(12) ! partial currents for all faces of mesh cell
    real(8) :: albedo      ! the albedo

    ! Get partial currents from object
    current = cmfd % current(:,g,i,j,k)

    ! Define xyz and +/- indices
    shift_idx = -2*mod(l,2) + 1          ! shift neig by -1 or +1

    ! Calculate albedo
    if ((shift_idx ==  1 .and. current(2*l  ) < 1.0e-10_8) .or. &
        (shift_idx == -1 .and. current(2*l-1) < 1.0e-10_8)) then
      albedo = ONE
    else
      albedo = (current(2*l-1)/current(2*l))**(shift_idx)
    end if

    ! Assign to function variable
    get_reflector_albedo = albedo

  end function get_reflector_albedo

!===============================================================================
! COMPUTE_EFFECTIVE_DOWNSCATTER changes downscatter rate for zero upscatter
!===============================================================================

  subroutine compute_effective_downscatter()

    use constants, only: ZERO, CMFD_NOACCEL
    use global,    only: cmfd

    integer :: nx                ! number of mesh cells in x direction
    integer :: ny                ! number of mesh cells in y direction
    integer :: nz                ! number of mesh cells in z direction
    integer :: ng                ! number of energy groups
    integer :: i                 ! iteration counter for x
    integer :: j                 ! iteration counter for y
    integer :: k                 ! iteration counter for z
    real(8) :: flux1             ! group 1 volume int flux
    real(8) :: flux2             ! group 2 volume int flux
    real(8) :: sigt1             ! group 1 total xs
    real(8) :: sigt2             ! group 2 total xs
    real(8) :: sigs11            ! scattering transfer 1 --> 1
    real(8) :: sigs21            ! scattering transfer 2 --> 1
    real(8) :: sigs12            ! scattering transfer 1 --> 2
    real(8) :: sigs22            ! scattering transfer 2 --> 2
    real(8) :: siga1             ! group 1 abs xs
    real(8) :: siga2             ! group 2 abs xs
    real(8) :: sigs12_eff        ! effective downscatter xs

    ! Extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! Return if not two groups
    if (ng /= 2) return

    ! Begin loop around space and energy groups
    ZLOOP: do k = 1, nz

      YLOOP: do j = 1, ny

        XLOOP: do i = 1, nx

          ! Check for active mesh
          if (allocated(cmfd%coremap)) then
            if (cmfd%coremap(i,j,k) == CMFD_NOACCEL) cycle
          end if

          ! Extract cross sections and flux from object
          flux1 = cmfd % flux(1,i,j,k)
          flux2 = cmfd % flux(2,i,j,k)
          sigt1 = cmfd % totalxs(1,i,j,k)
          sigt2 = cmfd % totalxs(2,i,j,k)
          sigs11 = cmfd % scattxs(1,1,i,j,k)
          sigs21 = cmfd % scattxs(2,1,i,j,k)
          sigs12 = cmfd % scattxs(1,2,i,j,k)
          sigs22 = cmfd % scattxs(2,2,i,j,k)

          ! Compute absorption xs
          siga1 = sigt1 - sigs11 - sigs12
          siga2 = sigt2 - sigs22 - sigs21

          ! Compute effective downscatter xs
          sigs12_eff = sigs12 - sigs21*flux2/flux1

          ! Recompute total cross sections (use effective and no upscattering)
          sigt1 = siga1 + sigs11 + sigs12_eff
          sigt2 = siga2 + sigs22

          ! Record total xs
          cmfd % totalxs(1,i,j,k) = sigt1
          cmfd % totalxs(2,i,j,k) = sigt2

          ! Record effective downscatter xs
          cmfd % scattxs(1,2,i,j,k) = sigs12_eff

          ! Zero out upscatter cross section
          cmfd % scattxs(2,1,i,j,k) = ZERO

        end do XLOOP

      end do YLOOP

    end do ZLOOP

  end subroutine compute_effective_downscatter

end module cmfd_data
