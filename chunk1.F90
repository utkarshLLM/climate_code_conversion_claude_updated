module micro_mg3_0
    !---------------------------------------------------------------------------------
    ! Purpose:
    !   MG microphysics version 3.0 - Update of MG microphysics with
    !                                 prognostic hail OR graupel.
    !
    ! Author: Andrew Gettelman, Hugh Morrison
    !
    ! Version 3 history: Sep 2016: development begun for hail, graupel 
    !
    ! Version 2 history: Sep 2011: Development begun.
    !                    Feb 2013: Added of prognostic precipitation.
    !                    Aug 2015: Published and released version
    ! Contributions from:  Sean Santos, Peter Caldwell, Xiaohong Liu and Steve Ghan
    !
    ! invoked in CAM by specifying -microphys=mg3
    !
    ! References: 
    !
    !           Gettelman, A. and H. Morrison, Advanced Two-Moment Microphysics for Global Models. 
    !
    !           Part I: Off line tests and comparisons with other schemes. 
    !
    !           J. Climate, 28, 1268-1287. doi: 10.1175/JCLI-D-14-00102.1, 2015. 
    !
    !
    !
    !           Gettelman, A., H. Morrison, S. Santos, P. Bogenschutz and P. H. Caldwell 
    !
    !           Advanced Two-Moment Microphysics for Global Models. 
    !
    !           Part II: Global model solutions and Aerosol-Cloud Interactions. 
    !
    !           J. Climate, 28, 1288-1307. doi:10.1175/JCLI-D-14-00103.1 , 2015. 
    !
    ! for questions contact Hugh Morrison, Andrew Gettelman
    ! e-mail: morrison@ucar.edu, andrew@ucar.edu
    !---------------------------------------------------------------------------------
    !
    ! NOTE: Modified to allow other microphysics packages (e.g. CARMA) to do ice
    ! microphysics in cooperation with the MG liquid microphysics. This is
    ! controlled by the do_cldice variable.
    !
    ! If do_cldice is false, then MG microphysics should not update CLDICE or
    ! NUMICE; it is assumed that the other microphysics scheme will have updated
    ! CLDICE and NUMICE. The other microphysics should handle the following
    ! processes that would have been done by MG:
    !   - Detrainment (liquid and ice)
    !   - Homogeneous ice nucleation
    !   - Heterogeneous ice nucleation
    !   - Bergeron process
    !   - Melting of ice
    !   - Freezing of cloud drops
    !   - Autoconversion (ice -> snow)
    !   - Growth/Sublimation of ice
    !   - Sedimentation of ice
    !
    ! This option has not been updated since the introduction of prognostic
    ! precipitation, and probably should be adjusted to cover snow as well.
    !
    !---------------------------------------------------------------------------------
    ! Version 3.O based on micro_mg2_0.F90 and WRF3.8.1 module_mp_morr_two_moment.F
    !---------------------------------------------------------------------------------
    ! Based on micro_mg (restructuring of former cldwat2m_micro)
    ! Author: Andrew Gettelman, Hugh Morrison.
    ! Contributions from: Xiaohong Liu and Steve Ghan
    ! December 2005-May 2010
    ! Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
    !                 Gettelman et al., 2010 J. Geophys. Res. - Atmospheres (G2010)
    ! for questions contact Hugh Morrison, Andrew Gettelman
    ! e-mail: morrison@ucar.edu, andrew@ucar.edu
    !---------------------------------------------------------------------------------
    ! Code comments added by HM, 093011
    ! General code structure:
    !
    ! Code is divided into two main subroutines:
    !   subroutine micro_mg_init --> initializes microphysics routine, should be called
    !                                  once at start of simulation
    !   subroutine micro_mg_tend --> main microphysics routine to be called each time step
    !                                this also calls several smaller subroutines to calculate
    !                                microphysical processes and other utilities
    !
    ! List of external functions:
    !   qsat_water --> for calculating saturation vapor pressure with respect to liquid water
    !   qsat_ice --> for calculating saturation vapor pressure with respect to ice
    !   gamma   --> standard mathematical gamma function
    ! .........................................................................
    ! List of inputs through use statement in fortran90:
    ! Variable Name                      Description                Units
    ! .........................................................................
    ! gravit          acceleration due to gravity                    m s-2
    ! rair            dry air gas constant for air                  J kg-1 K-1
    ! tmelt           temperature of melting point for water          K
    ! cpair           specific heat at constant pressure for dry air J kg-1 K-1
    ! rh2o            gas constant for water vapor                  J kg-1 K-1
    ! latvap          latent heat of vaporization                   J kg-1
    ! latice          latent heat of fusion                         J kg-1
    ! qsat_water      external function for calculating liquid water
    !                 saturation vapor pressure/humidity              -
    ! qsat_ice        external function for calculating ice
    !                 saturation vapor pressure/humidity              pa
    ! rhmini          relative humidity threshold parameter for
    !                 nucleating ice                                  -
    ! .........................................................................
    ! NOTE: List of all inputs/outputs passed through the call/subroutine statement
    !       for micro_mg_tend is given below at the start of subroutine micro_mg_tend.
    !---------------------------------------------------------------------------------
    
    ! Procedures required:
    ! 1) An implementation of the gamma function (if not intrinsic).
    ! 2) saturation vapor pressure and specific humidity over water
    ! 3) svp over ice
    
    #ifndef HAVE_GAMMA_INTRINSICS
    use shr_spfn_mod, only: gamma => shr_spfn_gamma
    #endif
    
    use wv_sat_methods, only: &
         qsat_water => wv_sat_qsat_water, &
         qsat_ice => wv_sat_qsat_ice
    
    ! Parameters from the utilities module.
    use micro_mg_utils, only: &
         r8, &
         pi, &
         omsm, &
         qsmall, &
         mincld, &
         rhosn, &
         rhoi, &
         rhow, &
         rhows, &
         ac, bc, &
         ai, bi, &
         aj, bj, &
         ar, br, &
         as, bs, &
         ag, bg, &
         ah, bh, &
         rhog,rhoh, &
         mi0, &
         rising_factorial
    
    implicit none
    private
    save
    
    public :: &
         micro_mg_init, &
         micro_mg_get_cols, &
         micro_mg_tend
    
    ! Switches for specification rather than prediction of droplet and crystal number
    ! note: number will be adjusted as needed to keep mean size within bounds,
    ! even when specified droplet or ice number is used
    !
    ! If constant cloud ice number is set (nicons = .true.),
    ! then all microphysical processes except mass transfer due to ice nucleation
    ! (mnuccd) are based on the fixed cloud ice number. Calculation of
    ! mnuccd follows from the prognosed ice crystal number ni.
    
    logical :: nccons ! nccons = .true. to specify constant cloud droplet number
    logical :: nicons ! nicons = .true. to specify constant cloud ice number
    logical :: ngcons = .false. ! ngcons = .true. to specify constant graupel number
    
    ! specified ice and droplet number concentrations
    ! note: these are local in-cloud values, not grid-mean
    real(r8) :: ncnst  ! droplet num concentration when nccons=.true. (m-3)
    real(r8) :: ninst  ! ice num concentration when nicons=.true. (m-3)
    real(r8) :: ngnst = 0.1e6_r8  ! graupel num concentration when ngcons=.true. (m-3)
    
    !=========================================================
    ! Private module parameters
    !=========================================================
    
    !Range of cloudsat reflectivities (dBz) for analytic simulator
    real(r8), parameter :: csmin = -30._r8
    real(r8), parameter :: csmax = 26._r8
    real(r8), parameter :: mindbz = -99._r8
    real(r8), parameter :: minrefl = 1.26e-10_r8    ! minrefl = 10._r8**(mindbz/10._r8)
    
    ! autoconversion size threshold for cloud ice to snow (m)
    real(r8) :: dcs
    
    ! minimum mass of new crystal due to freezing of cloud droplets done
    ! externally (kg)
    real(r8), parameter :: mi0l_min = 4._r8/3._r8*pi*rhow*(4.e-6_r8)**3
    
    ! Ice number sublimation parameter. Assume some decrease in ice number with sublimation if non-zero. Else, no decrease in number with sublimation. 
      real(r8), parameter :: sublim_factor =0.0_r8      !number sublimation factor. 
    
    !=========================================================
    ! Constants set in initialization
    !=========================================================
    
    ! Set using arguments to micro_mg_init
    real(r8) :: g           ! gravity
    real(r8) :: r           ! dry air gas constant
    real(r8) :: rv          ! water vapor gas constant
    real(r8) :: cpp         ! specific heat of dry air
    real(r8) :: tmelt       ! freezing point of water (K)
    
    ! latent heats of:
    real(r8) :: xxlv        ! vaporization
    real(r8) :: xlf         ! freezing
    real(r8) :: xxls        ! sublimation
    
    real(r8) :: rhmini      ! Minimum rh for ice cloud fraction > 0.
    
    ! flags
    logical :: microp_uniform
    logical :: do_cldice
    logical :: use_hetfrz_classnuc
    logical :: do_hail
    logical :: do_graupel
    
    real(r8) :: rhosu       ! typical 850mn air density
    
    real(r8) :: icenuct     ! ice nucleation temperature: currently -5 degrees C
    
    real(r8) :: snowmelt    ! what temp to melt all snow: currently 2 degrees C
    real(r8) :: rainfrze    ! what temp to freeze all rain: currently -5 degrees C
    
    ! additional constants to help speed up code
    real(r8) :: gamma_br_plus1
    real(r8) :: gamma_br_plus4
    real(r8) :: gamma_bs_plus1
    real(r8) :: gamma_bs_plus4
    real(r8) :: gamma_bi_plus1
    real(r8) :: gamma_bi_plus4
    real(r8) :: gamma_bj_plus1
    real(r8) :: gamma_bj_plus4
    real(r8) :: xxlv_squared
    real(r8) :: xxls_squared
    
    character(len=16)  :: micro_mg_precip_frac_method  ! type of precipitation fraction method
    real(r8)           :: micro_mg_berg_eff_factor     ! berg efficiency factor
    
    logical  :: allow_sed_supersat ! Allow supersaturated conditions after sedimentation loop
    logical  :: do_sb_physics ! do SB 2001 autoconversion or accretion physics
    
    !===============================================================================
    contains
    !===============================================================================
    
    subroutine micro_mg_init( &
         kind, gravit, rair, rh2o, cpair,    &
         tmelt_in, latvap, latice,           &
         rhmini_in, micro_mg_dcs,            &
         micro_mg_do_hail_in,micro_mg_do_graupel_in, &
         microp_uniform_in, do_cldice_in, use_hetfrz_classnuc_in, &
         micro_mg_precip_frac_method_in, micro_mg_berg_eff_factor_in, &
         allow_sed_supersat_in, do_sb_physics_in, &
         nccons_in, nicons_in, ncnst_in, ninst_in, ngcons_in, ngnst_in, errstring)
    
      use micro_mg_utils, only: micro_mg_utils_init
    
      !-----------------------------------------------------------------------
      !
      ! Purpose:
      ! initialize constants for MG microphysics
      !
      ! Author: Andrew Gettelman Dec 2005
      !
      !-----------------------------------------------------------------------
    
      integer,  intent(in)  :: kind         ! Kind used for reals
      real(r8), intent(in)  :: gravit
      real(r8), intent(in)  :: rair
      real(r8), intent(in)  :: rh2o
      real(r8), intent(in)  :: cpair
      real(r8), intent(in)  :: tmelt_in     ! Freezing point of water (K)
      real(r8), intent(in)  :: latvap
      real(r8), intent(in)  :: latice
      real(r8), intent(in)  :: rhmini_in    ! Minimum rh for ice cloud fraction > 0.
      real(r8), intent(in)  :: micro_mg_dcs
    
    !MG3 dense precipitating ice. Note, only 1 can be true, or both false.
      logical,  intent(in)  :: micro_mg_do_graupel_in    ! .true. = configure with graupel
                                                       ! .false. = no graupel (hail possible)
      logical,  intent(in)  :: micro_mg_do_hail_in    ! .true. = configure with hail
                                                       ! .false. = no hail (graupel possible)
      logical,  intent(in)  :: microp_uniform_in    ! .true. = configure uniform for sub-columns
                                                ! .false. = use w/o sub-columns (standard)
      logical,  intent(in)  :: do_cldice_in     ! .true. = do all processes (standard)
                                                ! .false. = skip all processes affecting
                                                !           cloud ice
      logical,  intent(in)  :: use_hetfrz_classnuc_in ! use heterogeneous freezing
    
      character(len=16),intent(in)  :: micro_mg_precip_frac_method_in  ! type of precipitation fraction method
      real(r8),         intent(in)  :: micro_mg_berg_eff_factor_in     ! berg efficiency factor
      logical,  intent(in)  ::  allow_sed_supersat_in ! allow supersaturated conditions after sedimentation loop
      logical,  intent(in)  ::  do_sb_physics_in ! do SB autoconversion and accretion physics
    
      logical, intent(in)   :: nccons_in
      logical, intent(in)   :: nicons_in
      real(r8), intent(in)  :: ncnst_in
      real(r8), intent(in)  :: ninst_in
    
      logical, intent(in)   :: ngcons_in
      real(r8), intent(in)  :: ngnst_in
    
      character(128), intent(out) :: errstring    ! Output status (non-blank for error return)
    
      !-----------------------------------------------------------------------
    
      dcs = micro_mg_dcs
    
      ! Initialize subordinate utilities module.
      call micro_mg_utils_init(kind, rair, rh2o, cpair, tmelt_in, latvap, latice, &
           dcs, errstring)
    
      if (trim(errstring) /= "") return
    
      ! declarations for MG code (transforms variable names)
    
      g= gravit                 ! gravity
      r= rair                   ! dry air gas constant: note units(phys_constants are in J/K/kmol)
      rv= rh2o                  ! water vapor gas constant
      cpp = cpair               ! specific heat of dry air
      tmelt = tmelt_in
      rhmini = rhmini_in
      micro_mg_precip_frac_method = micro_mg_precip_frac_method_in
      micro_mg_berg_eff_factor    = micro_mg_berg_eff_factor_in
      allow_sed_supersat          = allow_sed_supersat_in
      do_sb_physics               = do_sb_physics_in
    
      nccons = nccons_in
      nicons = nicons_in
      ncnst  = ncnst_in
      ninst  = ninst_in
      ngcons = ngcons_in
      ngnst  = ngnst_in
    
      ! latent heats
    
      xxlv = latvap         ! latent heat vaporization
      xlf  = latice         ! latent heat freezing
      xxls = xxlv + xlf     ! latent heat of sublimation
    
      ! flags
      microp_uniform = microp_uniform_in
      do_cldice  = do_cldice_in
      use_hetfrz_classnuc = use_hetfrz_classnuc_in
      do_hail = micro_mg_do_hail_in
      do_graupel = micro_mg_do_graupel_in
    
      ! typical air density at 850 mb
    
      rhosu = 85000._r8/(rair * tmelt)
    
      ! Maximum temperature at which snow is allowed to exist
      snowmelt = tmelt + 2._r8
      ! Minimum temperature at which rain is allowed to exist
      rainfrze = tmelt - 40._r8
    
      ! Ice nucleation temperature
      icenuct  = tmelt - 5._r8
    
      ! Define constants to help speed up code (this limits calls to gamma function)
      gamma_br_plus1=gamma(1._r8+br)
      gamma_br_plus4=gamma(4._r8+br)
      gamma_bs_plus1=gamma(1._r8+bs)
      gamma_bs_plus4=gamma(4._r8+bs)
      gamma_bi_plus1=gamma(1._r8+bi)
      gamma_bi_plus4=gamma(4._r8+bi)
      gamma_bj_plus1=gamma(1._r8+bj)
      gamma_bj_plus4=gamma(4._r8+bj)
    
      xxlv_squared=xxlv**2
      xxls_squared=xxls**2
    
    end subroutine micro_mg_init

end module micro_mg3_0