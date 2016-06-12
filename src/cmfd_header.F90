module cmfd_header 

  use constants,  only: CMFD_NOACCEL, ZERO, ONE

  implicit none
  private
  public :: allocate_cmfd, deallocate_cmfd

  type, public :: cmfd_type

    ! Indices for problem
    integer :: indices(4)

    ! Albedo boundary condition
    real(8) :: albedo(6)

    ! Core overlay map
    integer, allocatable :: coremap(:,:,:)
    integer, allocatable :: indexmap(:,:)
    integer :: mat_dim = CMFD_NOACCEL 

    ! Energy grid
    real(8), allocatable :: egrid(:)

    ! Batch index and reaction rates for moving/rolling windows right now
    integer :: idx
    real(8), allocatable :: flux_rate(:,:,:,:,:)
    real(8), allocatable :: total_rate(:,:,:,:,:)
    real(8), allocatable :: p1scatt_rate(:,:,:,:,:)
    real(8), allocatable :: scatt_rate(:,:,:,:,:,:)
    real(8), allocatable :: nfiss_rate(:,:,:,:,:,:)
    real(8), allocatable :: current_rate(:,:,:,:,:,:)
    real(8), allocatable :: quad_current_rate(:,:,:,:,:,:)
    real(8), allocatable :: openmc_src_rate(:,:,:,:,:)
    
    ! Cross sections
    real(8), allocatable :: totalxs(:,:,:,:)
    real(8), allocatable :: p1scattxs(:,:,:,:)
    real(8), allocatable :: scattxs(:,:,:,:,:)
    real(8), allocatable :: nfissxs(:,:,:,:,:)

    ! Diffusion coefficient
    real(8), allocatable :: diffcof(:,:,:,:)

    ! Current
    real(8), allocatable :: current(:,:,:,:,:)
    real(8), allocatable :: quad_current(:,:,:,:,:)

    ! Flux
    real(8), allocatable :: flux(:,:,:,:)

    ! Coupling coefficients and equivalence parameters
    real(8), allocatable :: dtilde(:,:,:,:,:)
    real(8), allocatable :: dhat(:,:,:,:,:)

    ! Dimensions of mesh cells ([hu,hv,hw],xloc,yloc,zloc)
    real(8), allocatable :: hxyz(:,:,:,:)

    ! Fission source distributions
    real(8), allocatable :: cmfd_src(:,:,:,:)
    real(8), allocatable :: loo_src(:,:,:,:)
    real(8), allocatable :: openmc_src(:,:,:,:)
    real(8), allocatable :: openmc_src_old(:,:,:,:)

    ! Total source distributions (fission + scattering)
    real(8), allocatable :: openmc_total_src(:,:,:,:)

    ! Source sites in each mesh box
    real(8), allocatable :: sourcecounts(:,:,:,:)

    ! Weight adjustment factors 
    real(8), allocatable :: weightfactors(:,:,:,:)

    ! Eigenvector/eigenvalue from cmfd run
    real(8), allocatable :: phi(:)
    real(8) :: keff = ZERO

    ! Eigenvector/eigenvalue from adjoint run
    real(8), allocatable :: adj_phi(:)
    real(8) :: adj_keff = ZERO

    ! Eigenvalue for loo run
    real(8) :: loo_keff = ZERO
    
    ! Residual for neutron balance
    real(8), allocatable :: resnb(:,:,:,:)

    ! Openmc source normalization factor
    real(8) :: norm = ONE

    ! Shannon entropy from cmfd and loo fission source
    real(8), allocatable :: entropy(:)
    real(8), allocatable :: loo_entropy(:)

    ! RMS of neutron balance equations
    real(8), allocatable :: balance(:)

    ! Dominance ratio
    real(8), allocatable :: dom(:)

    ! List of CMFD and LOO k
    real(8), allocatable :: k_cmfd(:)
    real(8), allocatable :: k_loo(:)

    ! Balance keff
    real(8) :: keff_bal

  end type cmfd_type

contains

!==============================================================================
! ALLOCATE_CMFD allocates all data in of cmfd type
!==============================================================================

  subroutine allocate_cmfd(this, n_batches, n_save)

    integer, intent(in)            :: n_batches ! number of batches in calc
    integer, intent(in)            :: n_save    ! number of batches to save
    type(cmfd_type), intent(inout) :: this      ! cmfd instance

    integer :: nx  ! number of mesh cells in x direction
    integer :: ny  ! number of mesh cells in y direction
    integer :: nz  ! number of mesh cells in z direction
    integer :: ng  ! number of energy groups

   ! Extract spatial and energy indices from object
    nx = this % indices(1)
    ny = this % indices(2)
    nz = this % indices(3)
    ng = this % indices(4)

    ! Allocate reaction rate counters for moving/rolling window
    if (.not. allocated(this % openmc_src_rate))allocate(this % openmc_src_rate(ng,nx,ny,nz,n_save))
    if (.not. allocated(this % flux_rate))   allocate(this % flux_rate(ng,nx,ny,nz,n_save))
    if (.not. allocated(this % total_rate))  allocate(this % total_rate(ng,nx,ny,nz,n_save))
    if (.not. allocated(this % p1scatt_rate))allocate(this % p1scatt_rate(ng,nx,ny,nz,n_save))
    if (.not. allocated(this % scatt_rate))  allocate(this % scatt_rate(ng,ng,nx,ny,nz,n_save))
    if (.not. allocated(this % nfiss_rate))  allocate(this % nfiss_rate(ng,ng,nx,ny,nz,n_save))
    if (.not. allocated(this % current_rate))allocate(this % current_rate(12,ng,nx,ny,nz,n_save))
    if (.not. allocated(this % quad_current_rate)) &
         allocate(this % quad_current_rate(16,ng,nx,ny,nz,n_save))
    
    ! Allocate flux, cross sections and diffusion coefficient
    if (.not. allocated(this % flux))       allocate(this % flux(ng,nx,ny,nz))
    if (.not. allocated(this % totalxs))    allocate(this % totalxs(ng,nx,ny,nz))
    if (.not. allocated(this % p1scattxs))  allocate(this % p1scattxs(ng,nx,ny,nz))
    if (.not. allocated(this % scattxs))    allocate(this % scattxs(ng,ng,nx,ny,nz))
    if (.not. allocated(this % nfissxs))    allocate(this % nfissxs(ng,ng,nx,ny,nz))
    if (.not. allocated(this % diffcof))    allocate(this % diffcof(ng,nx,ny,nz))

    ! Allocate dtilde and dhat 
    if (.not. allocated(this % dtilde))     allocate(this % dtilde(6,ng,nx,ny,nz))
    if (.not. allocated(this % dhat))       allocate(this % dhat(6,ng,nx,ny,nz))

    ! Allocate dimensions for each box (here for general case)
    if (.not. allocated(this % hxyz))       allocate(this % hxyz(3,nx,ny,nz))

    ! Allocate surface currents
    ! FIXME: for now 12 for CMFD, 16 for 2D LOO
    if (.not. allocated(this % current))    allocate(this % current(12,ng,nx,ny,nz))
    if (.not. allocated(this % quad_current)) allocate(this % quad_current(16,ng,nx,ny,nz))

    ! Allocate source distributions
    if (.not. allocated(this % cmfd_src)) &
         allocate(this % cmfd_src(ng,nx,ny,nz))
    if (.not. allocated(this % loo_src)) &
         allocate(this % loo_src(ng,nx,ny,nz))
    if (.not. allocated(this % openmc_src)) &
         allocate(this % openmc_src(ng,nx,ny,nz))
    if (.not. allocated(this % openmc_src_old)) &
         allocate(this % openmc_src_old(ng,nx,ny,nz))
    if (.not. allocated(this % openmc_total_src)) &
         allocate(this % openmc_total_src(ng,nx,ny,nz))

    ! Allocate source weight modification vars
    if (.not. allocated(this % sourcecounts)) allocate(this % sourcecounts(ng,nx,ny,nz))
    if (.not. allocated(this % weightfactors)) allocate(this % weightfactors(ng,nx,ny,nz))

    ! Allocate batchwise parameters
    if (.not. allocated(this % entropy)) allocate(this % entropy(n_batches))
    if (.not. allocated(this % loo_entropy)) allocate(this % loo_entropy(n_batches))
    if (.not. allocated(this % balance)) allocate(this % balance(n_batches))
    if (.not. allocated(this % dom)) allocate(this % dom(n_batches))
    if (.not. allocated(this % k_cmfd)) allocate(this % k_cmfd(n_batches))
    if (.not. allocated(this % k_loo)) allocate(this % k_loo(n_batches))

    ! Set everthing to 0 except weight multiply factors if feedback isnt on
    this % flux          = ZERO
    this % totalxs       = ZERO
    this % p1scattxs     = ZERO
    this % scattxs       = ZERO
    this % nfissxs       = ZERO
    this % openmc_src_rate= ZERO
    this % flux_rate     = ZERO
    this % total_rate    = ZERO
    this % p1scatt_rate  = ZERO
    this % scatt_rate    = ZERO
    this % nfiss_rate    = ZERO
    this % diffcof       = ZERO
    this % dtilde        = ZERO
    this % dhat          = ZERO
    this % hxyz          = ZERO
    this % current       = ZERO
    this % current_rate  = ZERO
    this % quad_current  = ZERO
    this % quad_current_rate  = ZERO
    this % cmfd_src      = ZERO
    this % loo_src       = ZERO
    this % openmc_src    = ZERO
    this % openmc_src_old &
                         = ZERO
    this % openmc_total_src &
                         = ZERO
    this % sourcecounts  = ZERO
    this % weightfactors = ONE
    this % balance       = ZERO
    this % dom           = ZERO
    this % k_cmfd        = ZERO
    this % k_loo         = ZERO
    this % entropy       = ZERO
    this % loo_entropy   = ZERO

    ! Set starting index to 1
    this % idx           = 1

  end subroutine allocate_cmfd

!===============================================================================
! DEALLOCATE_CMFD frees all memory of cmfd type 
!===============================================================================

  subroutine deallocate_cmfd(this)

    type(cmfd_type), intent(inout) :: this ! cmfd instance

    if (allocated(this % egrid))         deallocate(this % egrid)
    if (allocated(this % totalxs))       deallocate(this % totalxs)
    if (allocated(this % p1scattxs))     deallocate(this % p1scattxs)
    if (allocated(this % scattxs))       deallocate(this % scattxs)
    if (allocated(this % nfissxs))       deallocate(this % nfissxs)
    if (allocated(this % diffcof))       deallocate(this % diffcof)
    if (allocated(this % current))       deallocate(this % current)
    if (allocated(this % quad_current))  deallocate(this % quad_current)
    if (allocated(this % flux))          deallocate(this % flux)
    if (allocated(this % total_rate))    deallocate(this % total_rate)
    if (allocated(this % p1scatt_rate))  deallocate(this % p1scatt_rate)
    if (allocated(this % scatt_rate))    deallocate(this % scatt_rate)
    if (allocated(this % nfiss_rate))    deallocate(this % nfiss_rate)
    if (allocated(this % current_rate))  deallocate(this % current_rate)
    if (allocated(this % quad_current_rate))  deallocate(this % quad_current_rate)
    if (allocated(this % openmc_src_rate))deallocate(this % openmc_src_rate)    
    if (allocated(this % flux_rate))     deallocate(this % flux_rate)    
    if (allocated(this % dtilde))        deallocate(this % dtilde)
    if (allocated(this % dhat))          deallocate(this % dhat)
    if (allocated(this % hxyz))          deallocate(this % hxyz)
    if (allocated(this % coremap))       deallocate(this % coremap)
    if (allocated(this % indexmap))      deallocate(this % indexmap)
    if (allocated(this % phi))           deallocate(this % phi)
    if (allocated(this % sourcecounts))  deallocate(this % sourcecounts)
    if (allocated(this % weightfactors)) deallocate(this % weightfactors)
    if (allocated(this % cmfd_src))      deallocate(this % cmfd_src)
    if (allocated(this % loo_src))       deallocate(this % loo_src)
    if (allocated(this % openmc_src))    deallocate(this % openmc_src)
    if (allocated(this % openmc_src_old)) &
                                         deallocate(this % openmc_src_old)
    if (allocated(this % openmc_total_src)) &
                                         deallocate(this % openmc_total_src)
    if (allocated(this % balance))       deallocate(this % balance)
    if (allocated(this % dom))           deallocate(this % dom)
    if (allocated(this % k_cmfd))        deallocate(this % k_cmfd)
    if (allocated(this % k_loo))        deallocate(this % k_loo)
    if (allocated(this % entropy))       deallocate(this % entropy)
    if (allocated(this % loo_entropy))   deallocate(this % loo_entropy)
    if (allocated(this % resnb))         deallocate(this % resnb)

  end subroutine deallocate_cmfd

end module cmfd_header
