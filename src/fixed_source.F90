module fixed_source

#ifdef MPI
  use message_passing
#endif

  use constants,       only: ZERO, MAX_LINE_LEN
  use global
  use output,          only: write_message, header
  use particle_header, only: Particle
  use random_lcg,      only: set_particle_seed
  use source,          only: sample_external_source, copy_source_attributes
  use state_point,     only: write_state_point
  use string,          only: to_str
  use tally,           only: synchronize_tallies, setup_active_usertallies
  use trigger,         only: check_triggers
  use tracking,        only: transport

  implicit none

contains

  subroutine run_fixedsource()

    integer(8)     :: i ! index over histories in single cycle
    type(Particle) :: p

    if (master) call header("FIXED SOURCE TRANSPORT SIMULATION", level=1)

    ! Allocate particle and dummy source site
!$omp parallel
    allocate(source_site)
!$omp end parallel

    ! Turn timer and tallies on
    tallies_on = .true.
!$omp parallel
    call setup_active_usertallies()
!$omp end parallel
    call time_active % start()

    ! ==========================================================================
    ! LOOP OVER BATCHES
    BATCH_LOOP: do current_batch = 1, n_max_batches

      ! In a restart run, skip any batches that have already been simulated
      if (restart_run .and. current_batch <= restart_batch) then
        if (current_batch > n_inactive) n_realizations = n_realizations + 1
        cycle BATCH_LOOP
      end if

      call initialize_batch()

      ! Start timer for transport
      call time_transport % start()

      ! =======================================================================
      ! LOOP OVER PARTICLES
!$omp parallel do schedule(static) firstprivate(p)
      PARTICLE_LOOP: do i = 1, work

        ! Set unique particle ID
        p % id = (current_batch - 1)*n_particles + work_index(rank) + i

        ! set particle trace
        trace = .false.
        if (current_batch == trace_batch .and. current_gen == trace_gen .and. &
             work_index(rank) + i == trace_particle) trace = .true.

        ! set random number seed
        call set_particle_seed(p % id)

        ! grab source particle from bank
        call sample_source_particle(p)

        ! transport particle
        call transport(p)

      end do PARTICLE_LOOP
!$omp end parallel do

      ! Accumulate time for transport
      call time_transport % stop()

      call finalize_batch()

      if (satisfy_triggers) exit BATCH_LOOP

    end do BATCH_LOOP

    call time_active % stop()

    ! ==========================================================================
    ! END OF RUN WRAPUP

    if (master) call header("SIMULATION FINISHED", level=1)

  end subroutine run_fixedsource

!===============================================================================
! INITIALIZE_BATCH
!===============================================================================

  subroutine initialize_batch()

    call write_message("Simulating batch " // trim(to_str(current_batch)) &
         &// "...", 1)

    ! Reset total starting particle weight used for normalizing tallies
    total_weight = ZERO

  end subroutine initialize_batch

!===============================================================================
! FINALIZE_BATCH
!===============================================================================

  subroutine finalize_batch()

    ! Collect and accumulate tallies
    call time_tallies % start()
    call synchronize_tallies()
    call time_tallies % stop()

    ! Check_triggers
    if (master) call check_triggers()
#ifdef MPI
    call MPI_BCAST(satisfy_triggers, 1, MPI_LOGICAL, 0, &
         MPI_COMM_WORLD, mpi_err)
#endif
    if (satisfy_triggers .or. &
         (trigger_on .and. current_batch == n_max_batches)) then
      call statepoint_batch % add(current_batch)
    end if

    ! Write out state point if it's been specified for this batch
    if (statepoint_batch % contains(current_batch)) then
      call write_state_point()
    end if

  end subroutine finalize_batch

!===============================================================================
! SAMPLE_SOURCE_PARTICLE
!===============================================================================

  subroutine sample_source_particle(p)

    type(Particle), intent(inout) :: p

    ! Set particle
    call p % initialize()

    ! Sample the external source distribution
    call sample_external_source(source_site)

    ! Copy source attributes to the particle
    call copy_source_attributes(p, source_site)

    ! Determine whether to create track file
    if (write_all_tracks) p % write_track = .true.

  end subroutine sample_source_particle

end module fixed_source
