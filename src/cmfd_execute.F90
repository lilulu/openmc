module cmfd_execute

!==============================================================================
! CMFD_EXECUTE -- This module is the highest level cmfd module that controls the
! cross section generation, diffusion calculation, and source re-weighting
!==============================================================================
  use, intrinsic :: ISO_FORTRAN_ENV
  use global

  implicit none
  private
  public :: execute_acceleration, cmfd_init_batch, print_fission_sources

contains

!==============================================================================
! EXECUTE_ACCELERATION runs the CMFD calculation
!==============================================================================

  subroutine execute_acceleration()

    use cmfd_data,              only: set_up_cmfd
    use cmfd_solver,            only: cmfd_solver_execute
    use loo_solver,             only: loo_solver_execute
    use error,                  only: warning, fatal_error

    ! CMFD single processor on master
    if (master) then

      ! Start cmfd timer
      call time_cmfd % start()

      ! Create acceleration data from OpenMC tallies
      call set_up_cmfd()

      ! Call solver
      if (cmfd_on) then 
         call cmfd_solver_execute()
         
         ! Save k-effective
         cmfd % k_cmfd(current_batch) = cmfd % keff
         
         ! check to perform adjoint on last batch
         if (current_batch == n_batches .and. cmfd_run_adjoint) then
            call cmfd_solver_execute(adjoint=.true.)
         end if
      end if
      
      ! Run loo routine here if it is requested
      if (loo_on) then 
         call loo_solver_execute()      
         cmfd % k_loo(current_batch) = cmfd % loo_keff
      end if
      
      ! calculate fission source
      call calc_fission_source()

      ! print fission sources to file
      call print_fission_sources()
    end if
    
    ! calculate weight factors
    if (cmfd_feedback) call cmfd_reweight(.true.)

    ! stop cmfd timer
    if (master) call time_cmfd % stop()

  end subroutine execute_acceleration

!==============================================================================
! CMFD_INIT_BATCH handles cmfd options at the start of every batch
!==============================================================================

  subroutine cmfd_init_batch()

    use global,            only: cmfd_begin, loo_tally, cmfd_on, &
                                 loo_on, cmfd_reset, cmfd_run, &
                                 current_batch

    ! Check to activate CMFD diffusion and possible feedback this
    ! guarantees that when cmfd begins at least one batch of tallies
    ! are accumulated
    !
    ! The flag loo_tally is turned on one iteration before
    ! acceleration starts. It would not trigger the actual CMFD
    ! kernel. Instead it would parse the MC data and stores some
    ! values.
    if (loo_run .and. cmfd_begin == current_batch + 1) then
       loo_tally = .true.
    else if (loo_run .and. cmfd_begin == current_batch) then
       loo_on = .true.
       loo_tally = .false.
    end if
    if (cmfd_run .and. cmfd_begin == current_batch) then
       cmfd_on = .true.
    end if

    ! debug: special feature to turn off acceleration at active batches
    !if ((loo_on) .and. (current_batch == n_inactive + 1)) then 
    !   loo_on = .false.
    !end if
    

    ! If this is a restart run and we are just replaying batches leave
    if (restart_run .and. current_batch <= restart_batch) return

    ! Check to reset tallies when tally reset is requested
    if ((cmfd_run .or. loo_run) .and. cmfd_reset % contains(current_batch)) then
       call cmfd_tally_reset()
    end if

    ! Check to reset tallies during CMFD with moving window FIXME:
    ! this might not be necessary because the parameters were
    ! over-written anyway.
    if ((cmfd_run .or. loo_run) .and. cmfd_flush_every .and. &
         current_batch > cmfd_begin) then
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
    integer :: cnt_cmfd, cnt_loo ! counter for number of non-zero entries
    real(8) :: hxyz(3) ! cell dimensions of current ijk cell
    real(8) :: vol     ! volume of cell
    real(8), allocatable :: source(:,:,:,:)  ! tmp source array for entropy
    real(8), allocatable :: loo_src(:,:,:,:)
    real(8) :: avg     ! avg source if it is flat
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
    cnt_cmfd = 0
    cnt_loo = 0

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
              if (cmfd_on) then 
                 cmfd % cmfd_src(g,i,j,k) = sum(cmfd % nfissxs(:,g,i,j,k) * &
                      cmfd % phi(idx:idx + (ng - 1)))*vol
                 
                 if (cmfd % cmfd_src(g, i, j, k) > ZERO) cnt_cmfd = cnt_cmfd + 1
              end if

              if (loo_on) then 
                 if (cmfd % loo_src(g, i, j, k) > ZERO) cnt_loo = cnt_loo + 1
              end if

            end do GROUP

          end do XLOOP

        end do YLOOP

      end do ZLOOP

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
           cmfd % cmfd_src = cmfd % cmfd_src/sum(cmfd % cmfd_src)
           source = cmfd % cmfd_src*log(cmfd % cmfd_src)/log(TWO)
        end where
        where (cmfd % loo_src > ZERO)
           cmfd % loo_src = cmfd % loo_src / sum(cmfd % loo_src)
           loo_src = cmfd % loo_src*log(cmfd % loo_src)/log(TWO)
        end where

        ! Sum that source
        cmfd % entropy(current_batch) = -sum(source)
        cmfd % loo_entropy(current_batch) = -sum(loo_src)

        ! Deallocate tmp array
        if (allocated(source)) deallocate(source)
        if (allocated(loo_src)) deallocate(loo_src)

      end if

      ! Normalize source such that it average to 1.0
      if (sum(cmfd % cmfd_src) > 1e-5) then 
         cmfd % cmfd_src = cmfd % cmfd_src / sum(cmfd % cmfd_src) * cnt_cmfd
      end if
      if (sum(cmfd % loo_src) > 1e-5) then 
         cmfd % loo_src = cmfd % loo_src / sum(cmfd % loo_src) * cnt_loo
      end if
      
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
    use global,     only: cmfd, cmfd_coremap, master, current_batch

    integer :: nx      ! maximum number of cells in x direction
    integer :: ny      ! maximum number of cells in y direction
    integer :: nz      ! maximum number of cells in z direction
    integer :: ng      ! maximum number of energy groups
    integer :: i       ! iteration counter for x
    integer :: j       ! iteration counter for y
    integer :: k       ! iteration counter for z
    integer :: g       ! iteration counter for groups
    logical :: exist
    type(StructuredMesh), pointer :: m => null()

    ! Only perform for master
    if (master) then
       ! Open files
       inquire(file = "fs.dat", exist = exist)

       if (exist) then
          ! In one case we replace the file: this is the first run of a fresh run
          if ((.not. restart_run) .and. (current_batch == 1)) then
             open(unit = 2, file = "fs.dat", status = "replace", action = "write")
             ! If it is not one of the conditions where we replace, then we
             ! consider append to the file if acceleration is on.
          elseif (current_batch > 0) then
             open(unit = 2, file = "fs.dat", status = "old", position = "append", &
                  action = "write")
          end if
       else
          open(unit = 2, file = "fs.dat", status = "new", action = "write")
       end if

       if (cmfd_on .or. loo_on) then 
       ! Debug: for pin mesh, might need to turn off printing here
       !if (.false.) then
          ! Get maximum of spatial and group indices
          nx = cmfd % indices(1)
          ny = cmfd % indices(2)
          nz = cmfd % indices(3)
          ng = cmfd % indices(4)

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
       else
          m => entropy_mesh
          nx = m % dimension(1)
          ny = m % dimension(2)
          nz = m % dimension(3)
          ng = 1
          ! Loop around indices to map to entropy meshes
          do k = 1, nz

             do j = 1, ny

                do i = 1, nx

                   do g = 1, ng
                      ! TODO: temporarily put place holder such that
                      ! entropy_p is at the same location as in the
                      ! accelerated case such that we can use the same
                      ! script to process files
                      write(2, *) current_batch, g, i, j, k, &
                           i, j, k, entropy_p(g, i, j, k)          
                   end do
                end do
             end do
          end do
       end if
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
                           master, cmfd_begin, current_batch
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

    ! no reason to proceed if this is the LOO run that turns on before
    ! actual run
    if ((loo_run) .and. (cmfd_begin == current_batch + 1)) return

    ! do not feedback cmfd at first batch
    if ((cmfd_run .or. loo_run) .and. (current_batch == 1)) return

    ! FIXME: temporarily no feedback on and after 128 batches for -a128
    ! 
    !if (current_batch > n_inactive) return

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
      cmfd % sourcecounts = ZERO
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
         if (loo_on) then
            where(cmfd % loo_src > ZERO .and. cmfd % sourcecounts > ZERO)
               cmfd % weightfactors = cmfd % loo_src/sum(cmfd % loo_src)* &
                    sum(cmfd % sourcecounts) / cmfd % sourcecounts
            end where
         elseif (cmfd_on) then
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

    use global,  only: n_cmfd_tallies, cmfd_tallies, cmfd_flush_every
    use output,  only: write_message
    use tally,   only: reset_result

    integer :: n ! loop counter

    ! Print message. notice printing message is disalbed for rolling
    ! windows case because it's too distracting to have it printed out
    ! every batch
    if (.not. cmfd_flush_every) then 
       call write_message("CMFD tallies reset", 7)
    end if

    ! Begin loop around CMFD tallies
    do n = 1, n_cmfd_tallies

      ! Reset that tally
      cmfd_tallies(n) % n_realizations = 0
      call reset_result(cmfd_tallies(n) % results)

    end do

  end subroutine cmfd_tally_reset

end module cmfd_execute
