module loo_pass_data

  use constants,      only: CMFD_NOACCEL, ZERO
  use global,         only: cmfd, cmfd_coremap, keff
  use iso_c_binding,  only: c_int, c_double, c_loc

  implicit none
  private
  public :: pass_data_into_loo

  ! C interface
  interface
     subroutine new_loo(indices, keff, albedo, hxyz, flux, total_src_old, &
          totalxs, nfissxs, scattxs, &
          current, quad_current) bind (C, name='new_loo')
       use iso_c_binding
       type (c_ptr), value :: indices
       type (c_ptr), value :: keff
       type (c_ptr), value :: albedo
       type (c_ptr), value :: hxyz
       type (c_ptr), value :: flux
       type (c_ptr), value :: total_src_old
       type (c_ptr), value :: totalxs
       type (c_ptr), value :: nfissxs
       type (c_ptr), value :: scattxs
       type (c_ptr), value :: current
       type (c_ptr), value :: quad_current
     end subroutine new_loo
  end interface

contains

!===============================================================================
! PASS_DATA_INTO_LOO access the parameters from cmfd object and calls C++ function
!===============================================================================

  subroutine pass_data_into_loo()
    ! fortran data type (c_wrapper_type), target :: internal values
    integer (c_int), target :: indices(4)
    real (c_double), target :: k
    real (c_double), target:: albedo(6)
    real (c_double), allocatable, target:: hxyz(:,:,:,:)
    real (c_double), allocatable, target:: flux(:,:,:,:)
    real (c_double), allocatable, target:: total_src_old(:,:,:,:)
    real (c_double), allocatable, target:: totalxs(:,:,:,:)
    real (c_double), allocatable, target:: nfissxs(:,:,:,:,:)
    real (c_double), allocatable, target:: scattxs(:,:,:,:,:)
    real (c_double), allocatable, target:: current(:,:,:,:,:)
    real (c_double), allocatable, target:: quad_current(:,:,:,:,:)

    indices = cmfd % indices
    k = keff
    albedo = cmfd % albedo
    hxyz = cmfd % hxyz
    flux = cmfd % flux
    total_src_old = cmfd % openmc_total_src_old
    totalxs = cmfd % totalxs
    nfissxs = cmfd % nfissxs
    scattxs = cmfd % scattxs
    current = cmfd % current
    quad_current = cmfd % quad_current
    call new_loo(c_loc(indices), c_loc(k), c_loc(albedo), c_loc(hxyz), &
         c_loc(flux(1,1,1,1)), c_loc(total_src_old(1,1,1,1)), &
         c_loc(totalxs(1,1,1,1)), &
         c_loc(nfissxs(1,1,1,1,1)), c_loc(scattxs(1,1,1,1,1)), &
         c_loc(current(1,1,1,1,1)), c_loc(quad_current(1,1,1,1,1)))

  end subroutine pass_data_into_loo
end module loo_pass_data
