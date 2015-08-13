module loo_pass_data

  use global,         only: cmfd, keff, current_batch
  use iso_c_binding,  only: c_int, c_double, c_loc
  use, intrinsic :: ISO_FORTRAN_ENV

  implicit none
  private
  public :: pass_data_into_loo

  ! C interface
  interface
     function new_loo(indices, k, albedo, hxyz, flux, &
          src_old, &
          totalxs, nfissxs, scattxs, &
          current, quad_current, loo_src) &
          result(rms) bind (C)
       use iso_c_binding
       real (c_double) :: rms
       type (c_ptr), value :: indices
       type (c_ptr), value :: k
       type (c_ptr), value :: albedo
       type (c_ptr), value :: hxyz
       type (c_ptr), value :: flux
       type (c_ptr), value :: src_old
       type (c_ptr), value :: totalxs
       type (c_ptr), value :: nfissxs
       type (c_ptr), value :: scattxs
       type (c_ptr), value :: current
       type (c_ptr), value :: quad_current
       type (c_ptr), value :: loo_src
     end function new_loo
  end interface

contains

!===============================================================================
! PASS_DATA_INTO_LOO access the parameters from cmfd object and calls C++ function
!===============================================================================

  subroutine pass_data_into_loo()
    ! format for the following declaration:
    ! fortran data type (c_wrapper_type), target :: internal values
    integer (c_int), target:: indices(4)
    real (c_double), target:: k
    real (c_double), target:: albedo(6)
    real (c_double), allocatable, target:: hxyz(:,:,:,:)
    real (c_double), allocatable, target:: flux(:,:,:,:)
    real (c_double), allocatable, target:: src_old(:,:,:,:)
    real (c_double), allocatable, target:: totalxs(:,:,:,:)
    real (c_double), allocatable, target:: nfissxs(:,:,:,:,:)
    real (c_double), allocatable, target:: scattxs(:,:,:,:,:)
    real (c_double), allocatable, target:: current(:,:,:,:,:)
    real (c_double), allocatable, target:: quad_current(:,:,:,:,:)
    real (c_double), allocatable, target:: loo_src(:,:,:,:)
    real (c_double) rms

    indices = cmfd % indices
    k = keff
    albedo = cmfd % albedo
    hxyz = cmfd % hxyz
    flux = cmfd % flux
    src_old = cmfd % openmc_src_old
    totalxs = cmfd % totalxs
    nfissxs = cmfd % nfissxs
    scattxs = cmfd % scattxs
    current = cmfd % current
    quad_current = cmfd % quad_current
    loo_src = cmfd % loo_src

    rms = new_loo(c_loc(indices), c_loc(k), c_loc(albedo), &
         c_loc(hxyz), c_loc(flux(1,1,1,1)), &
         c_loc(src_old(1,1,1,1)), &
         c_loc(totalxs(1,1,1,1)), &
         c_loc(nfissxs(1,1,1,1,1)), c_loc(scattxs(1,1,1,1,1)), &
         c_loc(current(1,1,1,1,1)), c_loc(quad_current(1,1,1,1,1)), &
         c_loc(loo_src(1,1,1,1)))
    cmfd % src_cmp_loo(current_batch) = rms
    cmfd % loo_src = loo_src

  end subroutine pass_data_into_loo
end module loo_pass_data
