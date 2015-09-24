module cmfd_execute

!==============================================================================
! CMFD_EXECUTE -- This module is the highest level cmfd module that controls the
! cross section generation, diffusion calculation, and source re-weighting
!==============================================================================
  use loo_pass_data,      only: pass_data_into_loo

  use, intrinsic :: ISO_FORTRAN_ENV
  use global

  implicit none
  private
  public :: execute_cmfd, cmfd_init_batch

contains

!==============================================================================
! EXECUTE_CMFD runs the CMFD calculation
!==============================================================================

  subroutine execute_cmfd()

    use cmfd_data,              only: set_up_cmfd
    use cmfd_solver,            only: cmfd_solver_execute
    use error,                  only: warning, fatal_error

    ! CMFD single processor on master
    if (master) then

      ! Start cmfd timer
      call time_cmfd % start()

      ! Create cmfd data from OpenMC tallies
      call set_up_cmfd()

      ! Run loo routine here if it is requested
      if (loo_run) call pass_data_into_loo()

      ! Call solver
      call cmfd_solver_execute()

      ! Save k-effective
      cmfd % k_cmfd(current_batch) = cmfd % keff
      cmfd % k_loo(current_batch) = cmfd % loo_keff

      ! check to perform adjoint on last batch
      if (current_batch == n_batches .and. cmfd_run_adjoint) then
        call cmfd_solver_execute(adjoint=.true.)
      end if

    end if

    ! calculate fission source
    if (master) call calc_fission_source()

    ! print fission sources to file
    if (master) call print_fission_sources()

    ! calculate weight factors
    if (cmfd_feedback) call cmfd_reweight(.true.)

    ! stop cmfd timer
    if (master) call time_cmfd % stop()

  end subroutine execute_cmfd

!==============================================================================
! CMFD_INIT_BATCH handles cmfd options at the start of every batch
!==============================================================================

  subroutine cmfd_init_batch()

    use global,            only: cmfd_begin, loo_tally, cmfd_on, &
                                 cmfd_reset, cmfd_run,            &
                                 current_batch

    ! Check to activate CMFD diffusion and possible feedback this
    ! guarantees that when cmfd begins at least one batch of tallies
    ! are accumulated
    !
    ! The flag loo_tally is turned on one iteration becfore
    ! acceleration starts to. It would not trigger the actual CMFD
    ! kernel. Instead it would parse the MC data and stores some
    ! values.
    !
    ! FIXME: after loo flag is implemented, the following cmfd_run
    ! should be replaced by loo_run because we only care about storing
    ! an old copy of the tallies for LOO.
    if (cmfd_run .and. cmfd_begin == current_batch + 1) then
       loo_tally = .true.
    else
       loo_tally = .false.
    end if
    if (cmfd_run .and. cmfd_begin == current_batch) then
       cmfd_on = .true.
       loo_tally = .false.
    end if

    ! If this is a restart run and we are just replaying batches leave
    if (restart_run .and. current_batch <= restart_batch) return

    ! Check to reset tallies
    if (cmfd_run .and. cmfd_reset % contains(current_batch)) then
      call cmfd_tally_reset()
    end if

  end subroutine cmfd_init_batch

!===============================================================================
! CALC_FISSION_SOURCE calculates the cmfd fission source
!===============================================================================

  subroutine calc_fission_source()

    use constants,  only: CMFD_NOACCEL, ZERO, TWO
    use global,     only: cmfd, cmfd_coremap, master, entropy_on, current_batch
    use string,     only: to_str
    use error,            only: fatal_error

#ifdef MPI
    use global,     only: mpi_err
    use message_passing
#endif

    integer :: nx      ! maximum number of cells in x direction
    integer :: ny      ! maximum number of cells in y direction
    integer :: nz      ! maximum number of cells in z direction
    integer :: ng      ! maximum number of energy groups
    integer :: n       ! total size
    integer :: np      ! number of cells in 3D 
    integer :: i       ! iteration counter for x
    integer :: j       ! iteration counter for y
    integer :: k       ! iteration counter for z
    integer :: g       ! iteration counter for groups
    integer :: idx     ! index in vector
    real(8) :: hxyz(3) ! cell dimensions of current ijk cell
    real(8) :: vol     ! volume of cell
    real(8), allocatable :: source(:,:,:,:)  ! tmp source array for entropy
    real(8), allocatable :: loo_src(:,:,:,:)
    real(8) :: avg     ! avg source if it is flat
    real(8) :: temp   ! temporary counter
    real(8) :: temp1   ! temporary counter
    real(8) :: temp2   ! temporary counter
    real(8) :: temp3   ! temporary counter
    real(8) :: temp4   ! temporary counter 
    real(8) :: siga1
    real(8) :: siga2
    real(8) :: nusigf1
    real(8) :: nusigf2
    real(8) :: sigs12
    ! Get maximum of spatial and group indices
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)
    n  = ng*nx*ny*nz
    np = nx * ny * nz

    avg = 1.0 

    ! Allocate cmfd source if not already allocated and allocate buffer
    if (.not. allocated(cmfd % cmfd_src)) then
    !   allocate(cmfd % cmfd_src(ng,nx,ny,nz))
       call fatal_error("cmfd_src should be allocated by the time &
            & it reaches calc_fission_source")
    end if

    ! Reset cmfd source to 0
    cmfd % cmfd_src = ZERO

    ! Only perform for master
    if (master) then

      ! Loop around indices to map to cmfd object
      ZLOOP: do k = 1, nz

        YLOOP: do j = 1, ny

          XLOOP: do i = 1, nx

            GROUP: do g = 1, ng

              ! Check for core map
              if (cmfd_coremap) then
                if (cmfd % coremap(i,j,k) == CMFD_NOACCEL) then
                  cycle
                end if
              end if

              ! Get dimensions of cell
              hxyz = cmfd % hxyz(:,i,j,k)

              ! Calculate volume
              vol = hxyz(1)*hxyz(2)*hxyz(3)

              ! Get first index
              idx = get_matrix_idx(1,i,j,k,ng,nx,ny)

              ! Compute fission source
              cmfd % cmfd_src(g,i,j,k) = sum(cmfd % nfissxs(:,g,i,j,k) * &
                     cmfd % phi(idx:idx + (ng - 1)))*vol

            end do GROUP

          end do XLOOP

        end do YLOOP

      end do ZLOOP

      ! Normalize source such that it sums to 1.0
      cmfd % cmfd_src = cmfd % cmfd_src/sum(cmfd % cmfd_src)
      cmfd % loo_src = cmfd % loo_src / sum(cmfd % loo_src)

      ! Compute entropy
      if (entropy_on) then

        ! Allocate tmp array
         if (.not.allocated(source)) allocate(source(ng,nx,ny,nz))
         if (.not.allocated(loo_src)) allocate(loo_src(ng,nx,ny,nz))

        ! Initialize the source
        source = ZERO
        loo_src = ZERO

        ! Compute log
        where (cmfd % cmfd_src > ZERO)
           source = cmfd % cmfd_src*log(cmfd % cmfd_src)/log(TWO)
           loo_src = cmfd % loo_src*log(cmfd % loo_src)/log(TWO)
        end where

        ! Sum that source
        cmfd % entropy(current_batch) = -sum(source)
        cmfd % loo_entropy(current_batch) = -sum(loo_src)

        ! Deallocate tmp array
        if (allocated(source)) deallocate(source)
        if (allocated(loo_src)) deallocate(loo_src)

      end if

      ! Normalize source so average is 1.0
      !cmfd % cmfd_src = cmfd % cmfd_src/sum(cmfd % cmfd_src) * cmfd % norm
      !cmfd % loo_src = cmfd % loo_src/sum(cmfd % loo_src) * cmfd % norm

    end if

!#ifdef MPI
!    ! Broadcast full source to all procs
!    call MPI_BCAST(cmfd % cmfd_src, n, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
!    call MPI_BCAST(cmfd % loo_src, n, MPI_REAL8, 0, MPI_COMM_WORLD, mpi_err)
!#endif

  end subroutine calc_fission_source


!===============================================================================
! PRINT_FISSION_SOURCES print to file the mesh-cell homogenized
! fission sources generated by MC, CMFD and LOO
! ===============================================================================

  subroutine print_fission_sources()

    use constants,  only: CMFD_NOACCEL, ZERO, TWO
    use global,     only: cmfd, cmfd_coremap, master, cmfd_begin, current_batch

    integer :: nx      ! maximum number of cells in x direction
    integer :: ny      ! maximum number of cells in y direction
    integer :: nz      ! maximum number of cells in z direction
    integer :: ng      ! maximum number of energy groups
    integer :: i       ! iteration counter for x
    integer :: j       ! iteration counter for y
    integer :: k       ! iteration counter for z
    integer :: g       ! iteration counter for groups
    logical :: exist

    ! Only perform for master
    if (master) then
       ! Get maximum of spatial and group indices
       nx = cmfd % indices(1)
       ny = cmfd % indices(2)
       nz = cmfd % indices(3)
       ng = cmfd % indices(4)

       ! Open files
       inquire(file = "fs.dat", exist = exist)

       if (exist) then
          ! In two cases we replace the file: either this is the first of
          ! a restart fun (which starts at restart_batch + 1), or that
          ! this is the first of the acceleration run (which starts at cmfd_begin)
          if (((.not. restart_run) .and. (current_batch == cmfd_begin)) &
               .or. (restart_run .and. (current_batch == restart_batch + 1))) then
             open(unit = 2, file = "fs.dat", status = "replace", action = "write")
             ! If it is not one of the conditions where we replace, then we
             ! consider append to the file if acceleration is on.
          elseif (current_batch > cmfd_begin) then
             open(unit = 2, file = "fs.dat", status = "old", position = "append", &
                  action = "write")
          end if
       else
          open(unit = 2, file = "fs.dat", status = "new", action = "write")
       end if

       ! Loop around indices to map to cmfd object
       ZLOOP: do k = 1, nz

          YLOOP: do j = 1, ny

             XLOOP: do i = 1, nx

                GROUP: do g = 1, ng

                   ! Check for core map
                   if (cmfd_coremap) then
                      if (cmfd % coremap(i,j,k) == CMFD_NOACCEL) then
                         cycle
                      end if
                   end if

                   write(2, *) current_batch, g, i, j, k, &
                        cmfd % openmc_src_old(g, i, j, k), &
                        cmfd % cmfd_src(g, i, j, k), &
                        cmfd % loo_src(g, i, j, k), &
                        entropy_p(g, i, j, k), &
                        cmfd % openmc_src(g, i, j, k), &
                        entropy_s_old(g, i, j, k)
                end do GROUP
             end do XLOOP
          end do YLOOP
       end do ZLOOP

    ! Close file
    close(2)
    end if


end subroutine print_fission_sources

!===============================================================================
! CMFD_REWEIGHT performs weighting of particles in the source bank
!===============================================================================

  subroutine cmfd_reweight(new_weights)

    use constants,   only: ZERO, ONE
    use error,       only: warning, fatal_error
    use global,      only: meshes, source_bank, work, n_user_meshes, cmfd, &
                           master
    use mesh_header, only: StructuredMesh
    use mesh,        only: count_bank_sites, get_mesh_indices
    use search,      only: binary_search
    use string,      only: to_str

#ifdef MPI
    use global,      only: mpi_err
    use message_passing
#endif

    logical, intent(in) :: new_weights ! calcualte new weights

    integer :: nx       ! maximum number of cells in x direction
    integer :: ny       ! maximum number of cells in y direction
    integer :: nz       ! maximum number of cells in z direction
    integer :: ng       ! maximum number of energy groups
    integer :: i        ! iteration counter
    integer :: ijk(3)   ! spatial bin location
    integer :: e_bin    ! energy bin of source particle
    integer :: n_groups ! number of energy groups
    logical :: outside  ! any source sites outside mesh
    logical :: in_mesh  ! source site is inside mesh

    type(StructuredMesh), pointer :: m ! point to mesh

    ! Associate pointer
    m => meshes(n_user_meshes + 1)

    ! Get maximum of spatial and group indices
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    ! allocate arrays in cmfd object (can take out later extend to multigroup)
    if (.not.allocated(cmfd%sourcecounts)) then
      allocate(cmfd%sourcecounts(ng,nx,ny,nz))
      cmfd % sourcecounts = 0
    end if
    if (.not.allocated(cmfd % weightfactors)) then
      allocate(cmfd % weightfactors(ng,nx,ny,nz))
      cmfd % weightfactors = ONE
    end if

    ! Compute new weight factors
    if (new_weights) then

      ! Set weight factors to a default 1.0
      cmfd%weightfactors = ONE

      ! Count bank sites in mesh and reverse due to egrid structure
      call count_bank_sites(m, source_bank, cmfd%sourcecounts, cmfd % egrid, &
           sites_outside=outside, size_bank=work)
      cmfd % sourcecounts = cmfd%sourcecounts(ng:1:-1,:,:,:)

      ! Check for sites outside of the mesh
      if (master .and. outside) then
        call fatal_error("Source sites outside of the CMFD mesh!")
      end if

      ! Have master compute weight factors (watch for 0s)
      if (master) then
         if (loo_run) then
            where(cmfd % loo_src > ZERO .and. cmfd % sourcecounts > ZERO)
               cmfd % weightfactors = cmfd % loo_src/sum(cmfd % loo_src)* &
                    sum(cmfd % sourcecounts) / cmfd % sourcecounts
            end where
         else
            where(cmfd % cmfd_src > ZERO .and. cmfd % sourcecounts > ZERO)
               cmfd % weightfactors = cmfd % cmfd_src/sum(cmfd % cmfd_src)* &
                    sum(cmfd % sourcecounts) / cmfd % sourcecounts
            end where
         end if
      end if

      if (.not. cmfd_feedback) return

      ! Broadcast weight factors to all procs
#ifdef MPI
      call MPI_BCAST(cmfd % weightfactors, ng*nx*ny*nz, MPI_REAL8, 0, &
           MPI_COMM_WORLD, mpi_err)
#endif
    end if

    ! begin loop over source bank
    do i = 1, int(work,4)

      ! Determine spatial bin
      call get_mesh_indices(m, source_bank(i) % xyz, ijk, in_mesh)

      ! Determine energy bin
      n_groups = size(cmfd % egrid) - 1
      if (source_bank(i) % E < cmfd % egrid(1)) then
        e_bin = 1
        if (master) call warning('Source pt below energy grid')
      elseif (source_bank(i) % E > cmfd % egrid(n_groups + 1)) then
        e_bin = n_groups
        if (master) call warning('Source pt above energy grid')
      else
        e_bin = binary_search(cmfd % egrid, n_groups + 1, source_bank(i) % E)
      end if

      ! Reverese energy bin (lowest grp is highest energy bin)
      e_bin = n_groups - e_bin + 1

      ! Check for outside of mesh
      if (.not. in_mesh) then
        call fatal_error('Source site found outside of CMFD mesh')
      end if

      ! Reweight particle
      source_bank(i) % wgt = source_bank(i) % wgt * &
             cmfd % weightfactors(e_bin, ijk(1), ijk(2), ijk(3))

    end do

  end subroutine cmfd_reweight

!===============================================================================
! GET_MATRIX_IDX takes (x,y,z,g) indices and computes location in matrix
!===============================================================================

  function get_matrix_idx(g, i, j, k, ng, nx, ny) result (matidx)

    use global, only: cmfd, cmfd_coremap

    integer :: matidx ! the index location in matrix
    integer, intent(in) :: i  ! current x index
    integer, intent(in) :: j  ! current y index
    integer, intent(in) :: k  ! current z index
    integer, intent(in) :: g  ! current group index
    integer, intent(in) :: nx ! maximum number of cells in x direction
    integer, intent(in) :: ny ! maximum number of cells in y direction
    integer, intent(in) :: ng ! maximum number of energy groups

    ! Check if coremap is used
    if (cmfd_coremap) then

      ! Get idx from core map
      matidx = ng*(cmfd % coremap(i,j,k)) - (ng - g)

    else

      ! Compute index
      matidx = g + ng*(i - 1) + ng*nx*(j - 1) + ng*nx*ny*(k - 1)

    end if

  end function get_matrix_idx

!===============================================================================
! CMFD_TALLY_RESET resets all cmfd tallies
!===============================================================================

  subroutine cmfd_tally_reset()

    use global,  only: n_cmfd_tallies, cmfd_tallies
    use output,  only: write_message
    use tally,   only: reset_result

    real(8), allocatable :: openmc_src(:,:,:,:)

    integer :: n ! loop counter
    integer :: i
    integer :: j
    integer :: k
    integer :: h
    integer :: nx
    integer :: ny
    integer :: nz
    integer :: ng

    ! Extract spatial and energy indices from object
    nx = cmfd % indices(1)
    ny = cmfd % indices(2)
    nz = cmfd % indices(3)
    ng = cmfd % indices(4)

    allocate(openmc_src(ng,nx,ny,nz))
    openmc_src = ZERO

    ! Print message
    call write_message("CMFD tallies reset", 7)

    ! Save old fission source for LOO
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             do h = 1, ng
                openmc_src(h, i, j, k) = &
                     cmfd % openmc_src(h, i, j, k)
             end do
          end do
       end do
    end do

    ! Begin loop around CMFD tallies
    do n = 1, n_cmfd_tallies

      ! Reset that tally
      cmfd_tallies(n) % n_realizations = 0
      call reset_result(cmfd_tallies(n) % results)

    end do

    go to 100
    ! Disable this feature for now; FIXME: look into starting a batch
    ! right after a flush with  some fission source that's not just zero.
    ! Restore old fission source for LOO
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             do h = 1, ng
                cmfd % openmc_src(h, i, j, k) = &
                     openmc_src(h, i, j, k)
             end do
          end do
       end do
    end do
100 continue

    deallocate(openmc_src)
  end subroutine cmfd_tally_reset

end module cmfd_execute
