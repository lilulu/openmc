module loo_pass_data

  use constants,      only: CMFD_NOACCEL, ZERO
  use global,         only: cmfd, cmfd_coremap
  use iso_c_binding,  only: c_int, c_double, c_loc
  
  implicit none
  private
  public :: pass_data_into_loo

  ! C interface
  interface
     subroutine new_loo(indices, flux, totalxs, nfissxs, scattxs, current) &
          &bind (C, name='new_loo')
       use iso_c_binding
       type (c_ptr), value :: indices
       type (c_ptr), value :: flux
       type (c_ptr), value :: totalxs
       type (c_ptr), value :: nfissxs
       type (c_ptr), value :: scattxs
       type (c_ptr), value :: current
     end subroutine new_loo
  end interface

contains

!===============================================================================
! PASS_DATA_INTO_LOO access the parameters from cmfd object and calls C++ function
!===============================================================================

  subroutine pass_data_into_loo()
    ! fortran data type (c_wrapper_type), target :: internal values 
    integer (c_int), target :: indices(4)
    real (c_double), allocatable, target:: flux(:,:,:,:)
    real (c_double), allocatable, target:: totalxs(:,:,:,:)
    real (c_double), allocatable, target:: nfissxs(:,:,:,:,:)
    real (c_double), allocatable, target:: scattxs(:,:,:,:,:)
    real (c_double), allocatable, target:: current(:,:,:,:,:)

    indices = cmfd % indices
    flux = cmfd % flux
    totalxs = cmfd % totalxs
    nfissxs = cmfd % nfissxs
    scattxs = cmfd % scattxs
    current = cmfd % current
    call new_loo(c_loc(indices), c_loc(flux(1,1,1,1)),&
         & c_loc(totalxs(1,1,1,1)), c_loc(nfissxs(1,1,1,1,1)),&
         & c_loc(scattxs(1,1,1,1,1)), c_loc(current(1,1,1,1,1)))

  end subroutine pass_data_into_loo
end module loo_pass_data
