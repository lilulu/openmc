module loo_solver

  use global,         only: cmfd, keff, current_batch, entropy_s, entropy_p, &
       entropy_s_old, overall_gen, gen_per_batch, k_generation
  use iso_c_binding,  only: c_int, c_double, c_loc
  use, intrinsic :: ISO_FORTRAN_ENV

  implicit none
  private
  public :: loo_solver_execute

  ! C interface
  interface
     real (c_double) function new_loo(indices, k, albedo, hxyz, flux, &
          src_old, totalxs, nfissxs, scattxs, p1scattxs, &
          current, quad_current, loo_src) &
          bind (C)
       use iso_c_binding
       type (c_ptr), value :: indices
       type (c_ptr), value :: k
       type (c_ptr), value :: albedo
       type (c_ptr), value :: hxyz
       type (c_ptr), value :: flux
       type (c_ptr), value :: src_old
       type (c_ptr), value :: totalxs
       type (c_ptr), value :: nfissxs
       type (c_ptr), value :: scattxs
       type (c_ptr), value :: p1scattxs
       type (c_ptr), value :: current
       type (c_ptr), value :: quad_current
       type (c_ptr), value :: loo_src
     end function new_loo
  end interface

contains

!===============================================================================
! LOO_SOLVER_EXECUTE access the parameters from cmfd object and calls C++ function
!===============================================================================

  subroutine loo_solver_execute()
    ! format for the following declaration:
    ! fortran data type (c_wrapper_type), target :: internal values
    integer :: i
    integer (c_int), target :: indices(4)
    real (c_double) :: loo_k
    real (c_double), target :: k
    real (c_double), target :: albedo(6)
    real (c_double), allocatable, target :: hxyz(:,:,:,:)
    real (c_double), allocatable, target :: flux(:,:,:,:)
    real (c_double), allocatable, target :: src_old(:,:,:,:)
    real (c_double), allocatable, target :: totalxs(:,:,:,:)
    real (c_double), allocatable, target :: nfissxs(:,:,:,:,:)
    real (c_double), allocatable, target :: scattxs(:,:,:,:,:)
    real (c_double), allocatable, target :: p1scattxs(:,:,:,:)
    real (c_double), allocatable, target :: current(:,:,:,:,:)
    real (c_double), allocatable, target :: quad_current(:,:,:,:,:)
    real (c_double), allocatable, target :: loo_src(:,:,:,:)

    indices = cmfd % indices
    k = k_generation(overall_gen)
    albedo = cmfd % albedo
    hxyz = cmfd % hxyz
    flux = cmfd % flux
    !src_old = cmfd % openmc_src
    src_old = cmfd % openmc_src_old
    !src_old = entropy_s_old
    totalxs = cmfd % totalxs
    nfissxs = cmfd % nfissxs
    scattxs = cmfd % scattxs
    p1scattxs = cmfd % p1scattxs
    current = cmfd % current
    quad_current = cmfd % quad_current
    loo_src = cmfd % loo_src

    ! note: k value is not over-written, loo_src is;
    loo_k = new_loo(c_loc(indices), c_loc(k), c_loc(albedo), &
         c_loc(hxyz), c_loc(flux(1,1,1,1)), c_loc(src_old(1,1,1,1)), &
         c_loc(totalxs(1,1,1,1)), c_loc(nfissxs(1,1,1,1,1)), &
         c_loc(scattxs(1,1,1,1,1)), c_loc(p1scattxs(1,1,1,1)), &
         c_loc(current(1,1,1,1,1)), c_loc(quad_current(1,1,1,1,1)), &
         c_loc(loo_src(1,1,1,1)))

    cmfd % loo_src = loo_src
    cmfd % loo_keff = loo_k

    if (.false.) then
       write(*,*)'old mc fs'
       write(*,FMT='(F8.5, ", ", F8.5, ", ", F8.5, ", ", F8.5, ", ", F8.5)') &
            src_old(1,1,1,1), src_old(1,2,1,1), src_old(1,3,1,1), &
            src_old(1,4,1,1), src_old(1,5,1,1)

       !write(*,*)'removal'
       !write(*,FMT='(F8.5, ", ", F8.5, ", ", F8.5, ", ", F8.5, ", ", F8.5)') &
       !     totalxs(1,1,1,1) - scattxs(1,1,1,1,1), &
       !     totalxs(1,2,1,1) - scattxs(1,1,2,1,1), &
       !     totalxs(1,3,1,1) - scattxs(1,1,3,1,1), &
       !     totalxs(1,4,1,1) - scattxs(1,1,4,1,1), &
       !     totalxs(1,5,1,1) - scattxs(1,1,5,1,1)

       !write(*,*)'dtilde'
       !do i = 1, 5
       !   write(*,FMT='(  F8.5, ", ", F8.5, ", ", F8.5, ", ", F8.5, ", & 
       !        &", F8.5, ", ", F8.5)') &
       !        cmfd % dtilde(1,1,i,1,1), cmfd % dtilde(2,1,i,1,1), &
       !        cmfd % dtilde(3,1,i,1,1), cmfd % dtilde(4,1,i,1,1), &
       !        cmfd % dtilde(5,1,i,1,1), cmfd % dtilde(6,1,i,1,1)
       !end do

       !write(*,*)'dhat'
       !do i = 1, 5
       !   write(*,FMT='(  F8.5, ", ", F8.5, ", ", F8.5, ", ", F8.5, ", & 
       !        &", F8.5, ", ", F8.5)') &
       !        cmfd % dhat(1,1,i,1,1), cmfd % dhat(2,1,i,1,1), &
       !        cmfd % dhat(3,1,i,1,1), cmfd % dhat(4,1,i,1,1), &
       !        cmfd % dhat(5,1,i,1,1), cmfd % dhat(6,1,i,1,1)
       !end do
       
       !write(*,FMT='(  F8.5, ", ", F8.5, ", ", F8.5, ", ", F8.5, ", & 
       !     &", F8.5, ", ", F8.5)') &
       !     cmfd % albedo(1), cmfd % albedo(2), cmfd % albedo(3), cmfd % albedo(4), &
       !     cmfd % albedo(5), cmfd % albedo(6)

       write(*,*)'current'
       do i = 1, 5
          write(*,FMT='(  F8.5, ", ", F8.5, ", ", F8.5, ", ", F8.5)') &
               cmfd%current(1,1,i,1,1), cmfd%current(2,1,i,1,1), &
               cmfd%current(3,1,i,1,1), cmfd%current(4,1,i,1,1)
       end do

       write(*,*)'cmfd'
       write(*,FMT='(F8.5, ", ", F8.5, ", ", F8.5, ", ", F8.5, ", ", F8.5)') &
            5.0*cmfd%cmfd_src(1,1,1,1),  5.0*cmfd%cmfd_src(1,2,1,1), &
            5.0*cmfd%cmfd_src(1,3,1,1), 5.0*cmfd%cmfd_src(1,4,1,1), &
            5.0*cmfd%cmfd_src(1,5,1,1)

       write(*,*)'loo'
       write(*,FMT='(F8.5, ", ", F8.5, ", ", F8.5, ", ", F8.5, ", ", F8.5)') &
            loo_src(1,1,1,1), loo_src(1,2,1,1), &
            loo_src(1,3,1,1), loo_src(1,4,1,1), &
            loo_src(1,5,1,1)
    end if
 
  end subroutine loo_solver_execute
end module loo_solver
