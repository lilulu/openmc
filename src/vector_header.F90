module vector_header

  use constants,  only: ZERO

  implicit none
  private

# ifdef PETSC
#  include <finclude/petsc.h90>
# endif

  type, public :: Vector 
    integer :: n        ! number of rows/cols in matrix
    real(8), pointer :: val(:) ! matrix value vector
#  ifdef PETSC
    Vec :: petsc_vec
#  endif
   contains
     procedure :: create       => vector_create
     procedure :: destroy      => vector_destroy
     procedure :: add_value    => vector_add_value
     procedure :: setup_petsc  => vector_setup_petsc
     procedure :: write_petsc_binary => vector_write_petsc_binary
  end type Vector

contains

!===============================================================================
! VECTOR_CREATE allocates and initializes a vector 
!===============================================================================

  subroutine vector_create(self, n)

    integer       :: n
    class(Vector) :: self

    ! preallocate vector
    if (.not.associated(self % val)) allocate(self % val(n))

    ! set n
    self % n = n

    ! initialize to zero
    self % val = ZERO

  end subroutine vector_create

!===============================================================================
! VECTOR_DESTROY deallocates all space associated with a vector
!===============================================================================

  subroutine vector_destroy(self)

    class(Vector) :: self

    if (associated(self % val)) deallocate(self % val)

  end subroutine vector_destroy

!===============================================================================
! VECTOR_ADD_VALUE adds a value to the vector
!===============================================================================

  subroutine vector_add_value(self, idx, val)

    integer       :: idx
    real(8)       :: val
    class(Vector) :: self

    self % val(idx) = val

  end subroutine vector_add_value

!===============================================================================
! VECTOR_SETUP_PETSC
!===============================================================================

  subroutine vector_setup_petsc(self)

    class(Vector) :: self

    integer :: petsc_err

    ! link to petsc
#ifdef PETSC
    call VecCreateSeqWithArray(PETSC_COMM_WORLD, 1, self % n, self % val, &
         self % petsc_vec, petsc_err) 
#endif

  end subroutine vector_setup_petsc

!===============================================================================
! VECTOR_WRITE_PETSC_BINARY
!===============================================================================

  subroutine vector_write_petsc_binary(self, filename)

    character(*) :: filename
    class(Vector) :: self

    integer :: petsc_err
#ifdef PETSC
    PetscViewer :: viewer

    call PetscViewerBinaryOpen(PETSC_COMM_WORLD, trim(filename), &
         FILE_MODE_WRITE, viewer, petsc_err)
    call VecView(self % petsc_vec, viewer, petsc_err)
    call PetscViewerDestroy(viewer, petsc_err)
#endif

  end subroutine vector_write_petsc_binary

end module vector_header