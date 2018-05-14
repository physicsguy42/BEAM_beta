!--------------------------------------------------------------------
!
!  Author: Lindsey S. Chambers
!  This is a ray-tracing indirect Monte Carlo radiative transfer code.
!  It is meant to be used for modeling the I/F of Saturn's rings.
!  Subroutines are in alphabetical order following the main program.
!--------------------------------------------------------------------
! Changes by D. Olson
!  Changed seed for ran1 so that it will be different everytime code is run
!  modified input so the a planar albedo (R) is supplied by user
!  July 8, 2015 added open statement for dummy input filename in input()
!  Aug 3, 2015 modified obslong expression in input()
!  Aug 17, 2015 reset idum seed in init_pht_pos() and scater_phot()
!  April 2016 modified albedo() to load apf values from a 
!  tabulated file instead of calculating it each time the 
!  program is run.
!  Dec 2, 2016 fixed bug in flux_to_obs so that phase angle is updated 
!  on current direction of photon.
!
! July 19, 2017
!     corrected typo in flux_to_obs in Minnert scattering law. 
!     Was only giving Lambert result no matter the value of k -DMO
!
! Aug 24, 2017
!     Did some clean up and eliminated some redundant phase function options.
!
! September 7, 2017
!      Revised Europa phase function. Noticed Lindsey's quadratic model 
!      was giving unrealistic results and was negaive at high phase angles. 
!      It is now a power law like Callisto.
!
! March 6,2018
!      Added Shadowed Lambert point particle phase function option. 
!      This is mainly for consistancy checks with the facet scattering.
!      
!--------------------------------------------------------------------


! Source code can be compiled by typing
! > gfortran -O3 mc.F90 -o ../bin/mcx
!
! and run by entering,
! > ./mcx < mc-input.in > output

! Where mc-input.in is a file containing all the input parameters for mcx 


module variables      
  integer, parameter :: nscatter_max = 20 ! MAX NBR OF SCATTERINGS
  integer, parameter :: npmax = 1000000 ! MAX NBR OF PARTICLES
  integer, parameter :: nbins = 1000 ! MAX NBR BINS FOR INTEGRATION
  integer, parameter :: nxbox = 6 ! NBR OF X GRID BOXES FOR CODE SPEED-UP
  integer, parameter :: nybox = 6 ! NBR OF Y GRID BOXES FOR CODE SPEED-UP
  integer, parameter :: nzbox = 6 ! NBR OF Z GRID BOXES FOR CODE SPEED-UP
  integer, parameter :: nint = 10 ! NBR OF INTERVALS FOR UNIFORM GRID SETUP
  integer, parameter :: npts = 10000 ! NBR OF TEST VALUES TO GENERATE WHEN
  ! TESTING MU AND PHI
  integer, parameter :: nhistbins = 10 ! NBR OF BINS TO USE FOR TESTING 
  ! MU AND PHI
  integer, parameter :: nmu0 = 2001 ! NBR OF MU0 TO USE FOR FLUX WEIGHTING
  integer, parameter :: nmu = 2001 ! NBR OF MU TO USE FOR FLUX WEIGHTING
  integer, parameter :: nphi = 2001 ! NBR OF PHI TO USE FOR FLUX WEIGHTING
  integer :: nphot! NBR OF PHOTONS
  integer :: nphotsat ! NBR OF SATURNSHINE PHOTONS REQUESTED
  integer :: nphotpass, nphotpassSat, nphotpassSolar ! NBR OF PHOTONS THAT PASS THROUGH RINGS WITH NO SCATTERING
  integer :: npart ! NBR OF PARTICLES
  integer :: idum ! DUMMY VARIABLE FOR RANDOM NBR GENERATOR
  integer :: phasefunc ! NBR ID OF PARTICLE PHASE FUNCTION FOR RINGS
  integer :: satphasefunc ! NBR ID OF PARTICLE PHASE FUNCTION FOR SATURN
  integer :: nsatlat ! NBR SATSHINE LATITUDE BINS
  integer :: nsatlong ! NBR SATSHINE LONGITUDE BINS
  integer :: satswitch ! SWTICH TO TURN SATURNSHINE ON/OFF
  integer :: solarswitch ! SWITCH TO TURN SOLAR PHOTONS ON/OFF - DEBUGGING
  integer :: iphotbin ! COUNTER
  integer :: idebug ! SWITCH TO TURN WRITING DEBUG FILE ON/OFF
  integer :: npartbox(nxbox,nybox,nzbox) ! NRB OF PARTICLES IN EACH GRID 
  ! BOX FOR CODE SPEED-UP
  integer :: whichpartbox(nxbox,nybox,nzbox,npmax/(nxbox*nybox*nzbox)*3) 
  ! LIST OF WHICH PARTICLES (IDED BY THEIR NBR) ARE IN EACH GRID BOX FOR 
  ! CODE SPEED-UP
  integer :: ncells(nint) ! NBR OF CELLS IN EACH INTERVAL OF MU_0 FOR 
  ! UNIFORM GRID SETUP
  integer :: iter ! COUNTER
  integer :: nhist ! NBR OF BINS TO USE WHEN TESTING MU AND PHI
  
  real*8, parameter :: pi = 0.314159265358979D+01
  real*8, parameter :: satrad = 6.033d9 ! RADIUS OF SATURN IN CM
  real*8 :: sunlat,sunlong,obslat,obslong ! LATITUDE AND LONGITUDE OF SUN 
  ! AND OBSERVER; 0 LAT IS EQUATOR, 0 LONG IS THROUGH RING CALC PATCH
  real*8 :: alpha,sha ! PHASE ANGLES AND SOLAR HOUR ANGLE
  ! SHA IS MEASURED FROM MIDNIGHT, COUNTERCLOCKWISE
  real*8 :: raddeg ! CONVERSION FACTOR TO GO FROM DEGREES TO RADIANS
  real*8 :: flux(nscatter_max) ! ARRAY FOR FLUX AT EACH SCATTERING - DIRECT 
  ! SUN
  real*8 :: fluxsat(nscatter_max) ! ARRAY FOR FLUX AT EACH SCATTERING -
  !SATSHINE
  real*8 :: e_sun(3),e_obs(3) ! DIRECTION VECTORS FOR SUN AND OBSERVER IN
  ! RING-WAKE FRAME
  real*8 :: e_sunprime(3),e_obsprime(3) ! DIRECTION VECTORS FOR SUN AND 
  ! OBSERVER IN SATURN FRAME
  real*8 :: xmin,xmax,ymin,ymax,zmin,zmax,area ! DIMENSIONS OF CALC AREA
  real*8 :: p_part(3,npmax),rad_part(npmax) ! PARTICLE POSITION AND RADIUS 
  ! AT EACH SCATTERING
  real*8 :: k ! EXPONENT FOR MINNEART PHASE FUNC
  real*8 :: steep ! STEEPNESS PARAMETER FOR SHADOW-MINNAERT PHASE FUNC
  real*8 :: expo ! EXPONENT FOR POWER-LAW PHASE FUNC
  real*8 :: power_alpha(nbins+1),power_rand(nbins+1) ! ARRAYS CONTAINING 
  ! TABLE LOOK-UP INFO (RANDOM NUMBER AND ALPHA) FOR POWER-LAW PHASE FUNC
  real*8 :: ring_dist ! DISTANCE FROM SATURN CENTER TO RING CALCULATION AREA
  real*8 :: ngridphot(10000) ! ARRAY FOR NBR OF PHOTONS EXPECTED FROM 
  ! EACH GRID CELL ON SATURN
  real*8 :: snorm(3,10000) ! NORMALIZED DIRECTION VECTOR FOR SATSHINE PHOTONS
  ! IN RING-WAKE FRAME
  real*8 :: snormprime(3,10000) ! NORMALIZED DIRECTION VECTOR FOR SATSHINE 
  ! PHOTONS IN SATURN FRAME
  real*8 :: ngridphottot ! NBR OF SATURNSHINE PHOTONS SUMMED FROM GRID CELL
  real*8 :: xminbox,xmaxbox,yminbox,ymaxbox,zminbox,zmaxbox
  real*8 :: dxbox,dybox,dzbox ! EXTENT OF EACH GRID BOX FOR CODE SPEED-UP
  real*8 :: xboxmid(nxbox),yboxmid(nybox),zboxmid(nzbox)
  real*8 :: rbox(3,nxbox,nybox,nzbox) ! MIDPTS OF EACH GRID BOX FOR 
  ! CODE SPEED-UP
  real*8 :: phiwake ! WAKE ORIENTATION ANGLE
  real*8 :: cosphiwake, sinphiwake ! COS AND SIN OF PHIWAKE
  real*8 :: sunlongwake, obslongwake ! LONGITUDES IN RING-WAKE FRAME
  real*8 :: rcross(3) ! VECTOR FROM CENTER OF SATURN ALONG X'-AXIS IN 
  ! SATURN FRAME
  real*8 :: rcrossxy ! DISTANCE FROM SATURN CENTER TO PLACE IN X'-Y' PLANE
  ! WHERE E_SUNPRIME VECTOR CROSSES RING PLANE
  real*8 :: lambda ! DISTANCE FROM POINT ON SATURN TO X'-AXIS IN 
  ! SATURN FRAME
  real*8 :: mu_m ! PLACEHOLDER MU VALUE FOR SETTING UP SHADOW-MINNAERT PHASE
  ! FUNC
  real*8 :: intlength(4) ! LENGTHS OF INTERVALS IN EACH DIRECTION FOR 
  ! UNIFORM GRID SETUP
  real*8 :: limbound(4,2) ! LIMITS OF BOUNDING RECTANGLE FOR UNIFORM GRID
  real*8 :: gridcoords(nint+1,4) ! STORES COORDINATES OF GRID PTS FOR 
  ! UNIFORM GRID SETUP
  real*8 :: mu0np
  real*8 :: normconst,meantmp,binphilimit
  real*8 :: moments,histcounts(nhistbins,nhistbins)
  real*8 :: genpts(npts,2) ! ARRAY STORING MU AND PHI VALUES
  real*8 :: varin,reduced,func_to_reduce
  real*8 :: cdbleprime ! SMOOTH PARTICLE NORMALIZATION PARAMETER
  real*8 :: apf(nmu0) ! PLANAR ALBEDO
  real*8 :: curlypi ! BOND ALBEDO
  ! LOOKUP TABLES
  real*8 :: weight(nscatter_max) ! WEIGHT OF PHOTON AFTER SCATTERING
  integer :: time_temp
  real*8 :: frac_temp
  logical :: iflag  
  integer :: numrejects  ! number of scattering events blocked by <n,esun> criteria  
  integer :: numObstructed ! number of cases where flux is blocked by a particle 
  character*20 :: ftestfile   ! file that records the flux from each mpi process for testing
  character*20 :: photcount
  integer :: intnlPhotons  ! count internal photons
  integer :: numfluxSat, numfluxSolar  ! number of photons used to calculate flux
  real*8 :: c_S(11) ! normalization constants for lambert
  integer :: ishad
  

  type cellstruct ! CELL INFO FOR UNIFORM GRID
     integer :: mu_count
     integer :: phi_count
     integer :: gamma_count
     integer :: goodness
  end type cellstruct

  type(cellstruct) :: cell ! STORES INFO ON STARTER CELL FOR UNIFORM GRID
  type(cellstruct) :: cellarr(nint,10000)
!  type(cellstruct),allocatable :: cellarr(:,:) ! ARRAY STORING INFO ON 
! ALL CELLS FOR UNIFORM GRID, SIZE OF THE ARRAY WILL BE DETERMINE LATER 
! DURING 'DRY RUN' IN INIT_MINNAERT

end module variables

!--------------------------------------------------------------------

use variables
implicit none


call initialize ()
call input ()

if (satswitch == 1) call satshine()
if (solarswitch == 1) call ray_tracing()   
call output ()

end

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! SUBROUTINES, LISTED IN ALPHABETICAL ORDER, ARE BELOW
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!--------------------------------------------------------------------
! ALBEDO():
! CALCULATES THE WEIGHTING FACTORS FOR REFLECTION LAWS THAT NEED
! NUMERICAL INTEGRATION
! *** CURRENTLY ONLY WRITTEN FOR MINNAERT AND LAMBERT ***
! 
! April 20, 2016 - removed calculation of Apf in favor of loading 
! apf values from tabulated files
!    we only provide apf tables for k=1,1.05, and 1.15
!
! August 2, 2017 - fixed typo that was referencing wround array for 
!                  k=1.05 and k=1.15
!--------------------------------------------------------------------

subroutine albedo()
  use variables
  implicit none

  integer :: imu0

  real*8 :: phasenorm ! SEE FUNCTION
  real*8 :: smooth! SMOOTH = CALL TO PHASENORM WITH STEEP = 0
  real*8 :: dS ! spacing of sampled shadowing parameter
    
! give full directory for apf file  
  character*100 :: apf_file1 ! apf table with k = 1.0
  character*100 :: apf_file2 ! apf table with k = 1.05
  character*100 :: apf_file3 ! apf table with k = 1.15 
  integer, parameter :: nSteep = 21
  real*8 :: apf0(nmu0,nSteep),apf1(nmu0,nSteep),apf2(nmu0,nSteep)
  character*300 junk0
  character*300 junk1
  character*300 junk2
  integer :: iSteep  ! index of steepness parameter
  real*8 :: k0,k1,k2
  real*8 :: mu0_min, mu0_max
  real*8 :: dmu0
  character*80 :: Path, albedoPath
  integer :: islash
  
  call get_command_argument(0,Path)
    
  islash = index(Path,'/',BACK=.True.)
  
  albedoPath = Path(1:islash)//'../supportfiles/'
    
  apf_file1 = trim(albedoPath)//'apf_table_k1.dat'
  apf_file2 = trim(albedoPath)//'apf_table_k105.dat'
  apf_file3 = trim(albedoPath)//'apf_table_k115.dat'
  
  write(*,*)'Calculating planar albedo (APF)...'

! FIRST CALCULATE THE NORMALIZATION PARAMETERS: CDBLEPRIME

! IF SOLVING FOR SHADOWED LAMBERT, SET K=1 AND FOLLOW REST OF 
! ROUTINE FOR SOLVING SHADOWED MINNAERT 

  if (phasefunc.eq.7) then
     k = 1.0d0
  endif
  
  mu0_min = 0.d0
  mu0_max = 1.d0

  dS = 0.1d0
  iSteep = int(steep/dS) + 1  ! calculate index od steepness parameter
  dmu0 = (mu0_max-mu0_min)/dble(nmu0-1)

  k0=1.0d0
  k1=1.05d0
  k2=1.15d0
  
! USE NORMALIZATION FACTOR CDBLEPRIME ONLY FOR SPHERICAL PHASE
! FUNCTIONS

  if (phasefunc.eq.3.or.phasefunc.eq.4) then
     smooth = phasenorm(0.d0)
     cdbleprime = 2.d0 / smooth
  else
     cdbleprime = 1.d0
  endif
  
  open(unit=50,file=apf_file1, status='old',action = 'read')  
  open(unit=51,file=apf_file2, status='old',action = 'read')  
  open(unit=52,file=apf_file3, status='old',action = 'read')  

! load k=1 planar flux
  rewind(50)
  
  read(50,*) junk0
  do imu0=1,nmu0
    read(50,*) apf0(imu0,:)
  enddo

! load k=1.5 planar flux
  rewind(51)
  
  read(51,*) junk1
  do imu0=1,nmu0
    read(51,*) apf1(imu0,:)
  enddo

! load k=2 planar flux
  rewind(52)
  
  read(52,*) junk2
  do imu0=1,nmu0
    read(52,*) apf2(imu0,:)
  enddo
  
  close(50)
  close(51)
  close(52)
  
  if ( k .eq. k0) then
     do imu0 = 1,nmu0
        apf(imu0) = apf0(imu0,iSteep)
     enddo
     
  elseif(k .eq. k1) then
     do imu0 = 1,nmu0
        apf(imu0) = apf1(imu0,iSteep)
     enddo

  elseif (  k .eq. k2) then ! it's highly unlikely we will use a k>= 1.15
     do imu0 = 1,nmu0           ! but just in case we include this one
        apf(imu0) = apf2(imu0,iSteep)
     enddo
  else
     write(*,*) 'k value for Apf not found'
     stop
  endif  

end subroutine albedo

!--------------------------------------------------------------------
! FLUX_TO_OBS(P,E,INOW,NSCATTER):
! CALCULATE THE FLUX RECEIVED AT THE OBSERVER WHEN PHOTON
! SCATTERS OFF PARTICLE INOW.  IF VIEW IS BLOCKED, FLUX IS ZERO.
!
! July 19, 2017
!     corrected typo in Minnert scattering law. Was only giving Lambert result 
!     no matter the value of k -DMO
! April 30, 2018
!      Fixed photon blocking conditions so that we get the correct 
!      Flux when using particle phase functions.
!--------------------------------------------------------------------

subroutine flux_to_obs (p,e,inow,nscatter,iphot)
  use variables
  implicit none

  integer :: i,j
  integer :: inext,inow,nscatter ! COUNTERS, INOW WILL HOLD INDEX OF 
                                 ! PARTICLE IN QUESTION

  real*8 :: p(3),e(3),tnext,p_copy(3) ! PHOTON DIRECTION AND POSITION
  real*8 :: n(3),n2 ! UNIT NORMAL VECTOR, AND DOT PRODUCT OF N * N
  real*8 :: mu,mu0! MU = COS(EMERGENT ANGLE),MU0 = COS(INCIDENT ANGLE) 
  real*8 :: sf ! S_F = FLUX FROM PHASE FUNC
!  real*8 :: facto,ki,cn ! VARIABLES IN POWER-LAW PHASE FUNC, SEE DONES ET AL 1993
  real*8 :: facto(20),ki,cn ! VARIABLES IN POWER-LAW PHASE FUNC, SEE DONES ET AL 1993
  real*8 :: x
  real*8 :: talpha2 ! TAN(ALPHA/2)
  real*8 :: sta2 ! SQRT[TAN(ALPHA/2)]
  real*8 :: dmu0
  real*8 :: calpha, alphanew  ! for computing current phase angle
  integer :: imu1
  logical :: facto_init = .false.
  integer :: iphot
  real*8 :: Palpha

  save facto, facto_init ! PRESERVES THESE VARIABLES FROM ONE CALL TO THE NEXT
  
! index for Normalization constant c_S sampled at increments of 0.2
  
  ishad = int(steep/0.2d0+1.0d0)
  
! Count number of photons that are used for flux calculations  
  if (nscatter == 1) then
     if (satswitch == 1)  then
      numfluxSat = numfluxSat + 1
     else
      numfluxSolar = numfluxSolar + 1
     endif    
  endif   
  
  dmu0 = 1.0d0/dble(nmu0-1) !grid scacing in mu0, used to get correct value 
                            !of apf for weighting in case of facet 
                            ! (aka surface element) scattering. 

  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write (30,*) 'inside flux_to_obs'
     close(30)
  endif

! MAKE A COPY OF THIS PHOTON SO ITS POSITION DOESN'T GET CHANGED
 
  p_copy = p

! CALCULATE COMPONENTS OF UNIT NORMAL VECTOR FROM CURRENT PARTICLE
! POSITION TO CENTER OF PHOTON

  n(1:3) = p(1:3) - p_part(1:3,inow)
  n2 = dot_product(n,n)
  n = n / sqrt(n2)

! IGNORE THIS SCATTERING EVENT IF OBSERVER IS HIDDEN BY SCATTERING PARTICLE ITSELF
! ONLY IN FACET SCATTERING CASE

! For surface element scattering, if mu < 0 the scattering point is not visible 
! to the observer. So we do not calculate flux for any scattering event 
! for any of these cases.

! For particle phase functions the scattering particle is treated like a point
! (althought scattering is done on the surface rather then the particle's 
! center) so the observer will still recieve the flux. 

  mu = dot_product(n,e_obs)  
  if(phasefunc == 6 .or. phasefunc == 7) then
    if (mu < 0) then
      numrejects = numrejects +1
      return
     endif    
  endif   
  mu0 = dot_product(-e,n)
  
! IGNORE THIS SCATTERING EVENT IF OBSERVER IS HIDDEN BY ANOTHER PARTICLE
! We use which_part_next() with e=e_obs to determine if there is 
! a particle in the way of the observer seeing the flux.

  call which_part_next (p_copy,e_obs,inext,tnext,nscatter)

  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write (30,*) 'back from which_part_next, inside flux_to_obs'
     close(30)
  endif

! inext is the index of the particle that would be in the way 
! of the observer seeing the flux from the photon. A negative value 
! for inext means there is no ring particle in the way. If inext is 
! positive there is a obstruction between photon and observer. The index 
! of particle the photon is "sitting on" is inow. It is posible to get 
! inext equal to inow. However, this is the same as the check above for 
! "self-blocking" (mu < 0) so we filter out that case. 

  if ((inext > 0) .and. (inext .ne. inow)) then
     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write (30,*) 'scattering pt hidden from observer'
        close(30)
     endif
     numObstructed = numObstructed + 1
     return
  endif


  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write (30,*) 'n,e_obs,mu:',n,e_obs,mu
     close(30)
  endif
  
! update phase angle for multiple scatterings    
  calpha = dot_product(e_obs, -e)
  
  alphanew = acos(calpha)
  
  talpha2 = tan(alphanew/2.d0)
  sta2 = sqrt(talpha2)
  if(nscatter > 1) then
     if(phasefunc == 6 .or. phasefunc == 7) then
! find value in lookup table closest to calculated mu0       
       imu1 = int(mu0/dmu0+0.5d0)
! only update weight with apf for multiple scattering     
       weight(nscatter) = weight(nscatter-1)*apf(imu1)
     else
       weight(nscatter) = weight(nscatter-1) ! set all the weights equal
                                             ! for point phase functions
     endif                                        
  endif         
  
! CALCULATE FLUX TO OBSERVER, BASED ON PARTICLE PHASE FUNC

  if (phasefunc == 3) then
     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*) 'seen by observer'
        write(30,*) 'Callisto exponent = ',expo
        close(30)
     endif
     ki = 0.0d0 
! SETTING UP FOR FACTORIAL TERM THAT COMES INTO KI CALCULATION BELOW.
     if (.not. facto_init) then
        do i = 1, 20
           facto(i) = 1.0d0
           do j = 1, 2 * (i - 1) + 1
              facto(i) = facto(i) * j
           enddo
        enddo
        facto_init = .true.
     endif
     do i = 1, 20 ! KI OVERFLOWS BETWEEN 50 AND 75, DONES ET AL 93 USE 20
        ki = ki + pi ** (2.0d0 * (i - 1) + expo + 2.0d0) * (-1) ** (i - 1) &
             & / (2.0d0 * (i - 1) + expo + 2.0d0) / facto(i)
     enddo
     cn = 2.0d0 / ki
     sf = cn / (4.0d0 * pi) * (pi - acos(dot_product(-e,e_obs))) ** expo &
          & * weight(nscatter)
          
     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*) 'ki,cn,sf:',ki,cn,sf
        close(30)
     endif
          
     if (satswitch == 1) then
        fluxsat(nscatter) = fluxsat(nscatter) + sf
        if (idebug == 1) then
           open (30, file='debug.out', access='append')        
           write(30,*) 'nscatter:',nscatter
           write(30,*) 'fluxsat(nscatter)',fluxsat(nscatter)
           close(30)
        endif
     else
        flux(nscatter) = flux(nscatter) + sf        
        if (idebug == 1) then
           open (30, file='debug.out', access='append')
           write(30,*) 'nscatter:',nscatter
           write(30,*) 'flux(nscatter)',flux(nscatter)
           close(30)
        endif
     end if

     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write (30,*) 'idebug: ',idebug
        write(30,*) 'done with Callisto'
        close(30)
     endif

  else if (phasefunc == 4) then
     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*) 'seen by observer'
        close(30)
     endif
     
     ki = 0.0d0 
! SETTING UP FOR FACTORIAL TERM THAT COMES INTO KI CALCULATION BELOW.
     if (.not. facto_init) then
        do i = 1, 20
           facto(i) = 1.0d0
           do j = 1, 2 * (i - 1) + 1
              facto(i) = facto(i) * j
           enddo
        enddo
        facto_init = .true.
     endif
     do i = 1, 20 ! KI OVERFLOWS BETWEEN 50 AND 75, DONES ET AL 93 USE 20
        ki = ki + pi ** (2.0d0 * (i - 1) + expo + 2.0d0) * (-1) ** (i - 1) &
             & / (2.0d0 * (i - 1) + expo + 2.0d0) / facto(i)
     enddo
     cn = 2.0d0 / ki
     sf = cn / (4.0d0 * pi) * (pi - acos(dot_product(-e,e_obs))) ** expo &
          & * weight(nscatter)

     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*) 'ki,cn,sf:',ki,cn,sf
        close(30)
     endif


     if (satswitch == 1) then
        fluxsat(nscatter) = fluxsat(nscatter) + sf
     else
        flux(nscatter) = flux(nscatter) + sf
     end if

 else if (phasefunc == 6) then
     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*) 'seen by observer'
        write(30,*) 'steep, alpha, mu0, mu, k =',steep,alpha,mu0,mu,k
        close(30)
     endif
     sf = exp(-1.d0*steep*sta2)*(1.0d0+k)*(mu0**(k-1.0d0)*mu**k) &
          & /(2.0d0*pi)*weight(nscatter)
     
     if (satswitch == 1) then
        fluxsat(nscatter) = fluxsat(nscatter) + sf
     else
        flux(nscatter) = flux(nscatter) + sf                   
     end if

  else if (phasefunc == 7) then
     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write (30,*) 'seen by observer'
        close(30)
     endif
     sf = mu / pi * exp(-1.d0 * steep * sta2)*weight(nscatter)

     if (satswitch == 1) then
        fluxsat(nscatter) = fluxsat(nscatter) + sf
     else
        flux(nscatter) = flux(nscatter) + sf

     end if

  else if (phasefunc == 8) then
     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write (30,*) 'seen by observer'
        close(30)
     endif
     Palpha = (8.0/(3*pi))*(sin(alphanew)+(pi-alphanew)*cos(alphanew))
     sf = (c_S(ishad)/(4*pi))*exp(-1.d0 * steep * sta2)*Palpha*weight(nscatter)
     
     if (satswitch == 1) then
        fluxsat(nscatter) = fluxsat(nscatter) + sf
     else
        flux(nscatter) = flux(nscatter) + sf
     endif
  endif
  
     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*)'nscatter=',nscatter
        write(30,*)'phasefunc=',phasefunc
        write(30,*)'sf = ',sf
        write(30,*)'flux(nscatter) = ',flux(nscatter)
        write(30,*)'fluxsat(nscatter) before norm = ',fluxsat(nscatter)
        close(30)
     endif

  return
end subroutine flux_to_obs

!--------------------------------------------------------------------
! INIT_GRID_SPEEDUP():
! INITIAL SETUP OF GRID OF RING PARTICLES - USED TO SPEED UP CODE
!--------------------------------------------------------------------

subroutine init_grid_speedup()
  use variables
  implicit none

  integer :: ix,iy,iz,m,ip
  integer :: npartboxtot

  write(*,*)'Initialize grid set-up for ring partcles..'

  dxbox = (xmax - xmin) / dble(nxbox)
  dybox = (ymax - ymin) / dble(nybox)
  dzbox = (zmax - zmin) / dble(nzbox)

  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write(30,*) 'Setting up grid for code speed-up'
     write(30,*) 'xmin, xmax = ',xmin,xmax
     write(30,*) 'ymin, ymax = ',ymin,ymax
     write(30,*) 'zmin, zmax = ',zmin,zmax
     write(30,*) 'nxbox, nybox, nzbox = ',nxbox,nybox,nzbox
     write(30,*) 'dxbox, dybox, dzbox = ',dxbox,dybox,dzbox
     close(30)
  endif

  npartboxtot = 0

  xminbox = xmin

  do ix = 1, nxbox
     xmaxbox = xminbox + dxbox
     xboxmid(ix) = xminbox + (xmaxbox - xminbox) / 2.0d0

     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*) 'xminbox, xmaxbox, xboxmid = ',xminbox,xmaxbox,xboxmid(ix)
        close(30)
     endif

     xminbox = xmaxbox
  enddo

  yminbox = ymin

  do iy = 1, nybox
     ymaxbox = yminbox + dybox
     yboxmid(iy) = yminbox + (ymaxbox - yminbox) / 2.0d0
     
     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*) 'yminbox, ymaxbox, yboxmid = ',yminbox,ymaxbox,yboxmid(iy)
        close(30)
     endif

     yminbox = ymaxbox
  enddo

  zminbox = zmin

  do iz = 1, nzbox
     zmaxbox = zminbox + dzbox
     zboxmid(iz) = zminbox + (zmaxbox - zminbox) / 2.0d0
     
     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*) 'zminbox, zmaxbox, zboxmid = ',zminbox,zmaxbox,zboxmid(iz)
        close(30)
     endif

     zminbox = zmaxbox
  enddo

! GRID BOXES ARE NUMBER 1 TO NXBOX*NYBOX*NZBOX, GOING ROW BY ROW; SO IF 
! NXBOX*NYBOX*NZBOX = 216, BOXES 1-6 ARE THE FIRST ROW (XMIN TO XMAX BUT ALL 
! AT YMIN AND ZMIN)

  do iz = 1,nzbox
     do iy = 1, nybox
        do ix = 1, nxbox
           rbox(1,ix,iy,iz) = xboxmid(ix)
           rbox(2,ix,iy,iz) = yboxmid(iy)
           rbox(3,ix,iy,iz) = zboxmid(iz)

           if (idebug == 1) then
              open (30, file='debug.out', access='append')
              write(30,*) 'ix, iy, iz, rbox(1:3,ix,iy,iz):',&
                   & ix,iy,iz,rbox(1:3,ix,iy,iz)
              close(30)
           endif

        enddo
     enddo
  enddo

! DETERMINING WHICH BOX EACH PARTICLE GOES IN AND HOW MANY PARTICLES ARE 
! IN EACH BOX

  do ip = 1, npart
     do ix = 1, nxbox
        xminbox = xmin + dxbox * dble(ix - 1)
        xmaxbox = xminbox + dxbox
        if (p_part(1,ip).ge.xminbox.and.p_part(1,ip).lt.xmaxbox) then
           do iy = 1, nybox
              yminbox = ymin + dybox * dble(iy - 1)
              ymaxbox = yminbox + dybox
              if (p_part(2,ip).ge.yminbox.and.p_part(2,ip).lt.ymaxbox) then
                 do iz = 1, nzbox
                    zminbox = zmin + dzbox * dble(iz - 1)
                    zmaxbox = zminbox + dzbox
                    if (p_part(3,ip).ge.zminbox.and. &
                         & p_part(3,ip).lt.zmaxbox) then
                       npartbox(ix,iy,iz) = npartbox(ix,iy,iz) + 1
                       whichpartbox(ix,iy,iz,npartbox(ix,iy,iz)) = ip
                    endif
                 enddo
              endif
           enddo
        endif
     enddo
  enddo

  if (idebug == 1) then
     open (30, file='debug.out', access='append')
    
     do ix = 1, nxbox
        do iy = 1, nybox
           do iz = 1, nzbox
              write(30,*) 'npartbox(',ix,',',iy,',',iz,') =',npartbox(ix,iy,iz)
              npartboxtot = npartboxtot + npartbox(ix,iy,iz)
              do m = 1, npartbox(iy,iy,iz)
                 write(30,*) 'whichpartbox(',ix,',',iy,',',iz,',',m,') =',&
                      & whichpartbox(ix,iy,iz,m)
              enddo
           enddo
        enddo
     enddo

     write(30,*) 'npartboxtot = ',npartboxtot

     close(30)
  endif
  
  return
end subroutine init_grid_speedup 

!--------------------------------------------------------------------
! INIT_PHOT_POS(P,E):
! INITIAL PHOTON POSITION
!--------------------------------------------------------------------

subroutine init_phot_pos (p,e)
  use variables
  implicit none

  real*8 :: p(3), e(3) ! PHOTON POSITION AND DIRECTION
  real :: ran1 ! RANDOM NUMBER

! INITIALLY POSITION PHOTON JUST ABOVE OR JUST BELOW PARTICLE FIELD
! AND RANDOMLY WITHIN X AND Y RANGE
! INITIAL DIRECTION OF PHOTON IS FROM THE SUN FOR SOLAR PHOTONS; 
! FROM SATURN IF USING SATSHINE

  p(1) = xmin + (xmax - xmin) * dble( ran1(idum) )
  p(2) = ymin + (ymax - ymin) * dble( ran1(idum) )

  if (satswitch == 0) then

     if (e_sun(3) > 0) then
        p(3) = zmax + (zmax - zmin) * 0.01d0
     else
        p(3) = zmin - (zmax - zmin) * 0.01d0
     end if
     
     e = -e_sun

     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write (30,*) 'Photon initial position: ',p(1:3)
        write (30,*) 'Photon initial direction: ',e(1:3)
        close(30)
     endif

  else if (satswitch == 1) then
     
     if (snorm(3,iphotbin) < 0) then
        p(3) = zmax + (zmax - zmin) * 0.01d0
     else
        p(3) = zmin - (zmax - zmin) * 0.01d0
     end if

     e(1) = snorm(1,iphotbin)
     e(2) = snorm(2,iphotbin)
     e(3) = snorm(3,iphotbin)

     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write (30,*) 'Satshine photon initial position: ',p(1:3)
        write (30,*) 'Satshine photon initial direction: ',e(1:3)
        close(30)
     endif

  end if

  return
end subroutine init_phot_pos

!--------------------------------------------------------------------
! INIT_POWER():
! INITIALIZE FOR POWER-LAW PHASE FUNCTION TABLE LOOK-UP
!--------------------------------------------------------------------

subroutine init_power ()
  use variables
  implicit none

  integer :: i

  real*8 :: simpsons ! SIMPSON'S RULE FUNCTION
  real*8 :: power_den ! DENOMINATOR OF POWER-LAW INTEGRATION
  real*8 :: power_num ! NUMERATOR OF POWER-LAW INTEGRATION
  
! USING SIMPON'S RULE TO NUMERICALLY INTEGRATE
! SETTING UP TABLE TO USE WITH EQN 9 FROM SALO AND KARJALAIMEN (2003)

  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write (30,*) 'inside init_power'
     close(30)
  endif

  do i = 1,nbins+1
     power_alpha(i) = ( pi * (i - 1) ) / dble(nbins)
  enddo

! DENOMINATOR (NORMALIZATION FACTOR) FROM EQN 9 IN SK2003
 
  power_den = simpsons(0.0d0,pi)

  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write (30,*) 'power_den = ',power_den
     close(30)
  endif

  do i = 1,nbins+1
     power_num = simpsons(0.0d0,power_alpha(i))
     power_rand(i) = power_num / power_den

     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write (30,*) 'power_num, power_rand(i) = ',power_num,power_rand(i)
        close(30)
  endif

enddo
  
end subroutine init_power

!--------------------------------------------------------------------
! INITIALIZE():
! INITIALIZATION
!--------------------------------------------------------------------

subroutine initialize ()
  use variables
  implicit none

  write(*,*)'Initializing variables...'

  alpha = 0.0d0
  apf = 0.d0
  area = 0.0d0
  binphilimit = 0.0d0
  cdbleprime = 0.d0
  curlypi = 0.d0
  dxbox = 0.0d0
  dybox = 0.0d0
  e_obs = 0.0d0
  e_obsprime = 0.0d0
  e_sun = 0.0d0
  e_sunprime = 0.0d0
  expo = 0.0d0
  flux = 0.0d0
  fluxsat = 0.0d0
  func_to_reduce = 0.d0
  genpts = 0.d0
  gridcoords = 0.d0
  histcounts = 0.d0
  idebug = 0
 ! idum = -1234
  intlength = 0.d0
  iter = 0
  k = 0.0d0
  lambda = 0.0d0
  limbound = 0.d0
  meantmp = 0.d0
  moments = 0.0d0
  mu_m = 0.d0
  mu0np = 0.d0
  ncells = 0
  nhist = 0
  normconst = 0.d0
  npart = 0
  nphot = 0
  nphotpass = 0
  nphotpassSat = 0
  nphotpassSolar = 0
  nphotsat = 1 ! SET TO 1 TO AVOID DIVISION BY ZERO IN CASE NOT USING SATSHINE
  ngridphot = 0.0d0
  ngridphottot = 0.0d0
  nsatlat = 0
  nsatlong = 0
  obslat = 0.0d0
  obslong = 0.0d0
  obslongwake = 0.0d0
  p_part= 0.0d0
  phasefunc = 0
  phiwake = 0.0d0
  power_alpha = 0.0d0
  power_rand = 0.0d0
  rad_part = 0.0d0
  raddeg = 0.0d0
  reduced = 0.d0
  rcross = 0.0d0
  rcrossxy = 0.0d0
  ring_dist = 0.0d0
  satphasefunc = 0
  satswitch = 0
  steep = 0.0d0
  sha = 0.0d0
  solarswitch = 1
  sunlat = 0.0d0
  sunlong = 0.0d0
  sunlongwake = 0.0d0
  varin = 0.d0
  weight = 0.d0
  xmax = 0.0d0
  xmin = 0.0d0
  xmaxbox = 0.0d0
  xminbox = 0.0d0
  ymax = 0.0d0
  ymin = 0.0d0
  ymaxbox = 0.0d0
  yminbox = 0.0d0
  zmax = 0.0d0
  zmin = 0.0d0
  iflag = .true.
  numrejects = 0
  numObstructed = 0
  intnlPhotons = 0
  numfluxSat = 0
  numfluxSolar = 0
  
! to ensure we get a unique seed everytime code is executed 
! we base it on system time.
  
  if (iflag .eqv. .true. ) then
	  call system_clock(time_temp)
	  frac_temp = real(time_temp)/10000
	  idum = -1*int(10000*(frac_temp-real(int(frac_temp))))
  else
          idum = -12345
  endif
  

  return
end subroutine initialize

!--------------------------------------------------------------------
! INPUT():
! INPUT
!--------------------------------------------------------------------

subroutine input ()
  use variables
  implicit none

  integer :: i
  character*200 :: infile,ctemp
  logical :: test
  integer :: status
  
  real*8 :: tempRatio
  
  write(*,*) 'input variables'

  write(*,*)'Write debug file? 1 = yes, 0 = no'
  read(*,*)idebug
  
  write(*,*) 'debug = ', idebug

  if (idebug == 1) then
     open (30, file='debug.out',status='replace')
     write(*,*)'Writing debug file'
     close(30)
  endif

  write(*,*)'Enter nbr of photons from sun:'
  read(*,*) nphot
  
  write(*,*) 'nphot = ', nphot


  write(*,*)'Enter sun latitude:'
  read(*,*) sunlat

  write(*,*) 'sunlat = ', sunlat


  write(*,*)'Enter observer latitude:'
  read(*,*) obslat
  
  write(*,*)'Enter phase angle:'
  read(*,*) alpha

  write(*,*)'Enter solar hour angle in degress:'
  read(*,*) sha
  
  sunlong = 180.0 - sha

! CONVERT DEGREES TO RADIANS

  raddeg = pi / 180.0d0 

  sunlat = sunlat * raddeg
  sunlong = sunlong * raddeg
  obslat = obslat * raddeg
  alpha = alpha * raddeg
  sha = sha * raddeg
  
! Normalization constant for lambert phase function at
! S=0.0 to 2.0 at increments of 0.2. Taken from table 2 Cuzzietal2017 
! "Rough surfaces: Is the dark stuff just shadow?"

  c_S = (/1.00,1.16,1.35,1.56,1.80,2.08,2.38,2.73,3.12,3.57,4.06/)
  
! If alpha is less then abs(obslat-sunlat) or 
! greater then 180-(abslat+sunlat) there will be message to check the 
! three angles because in implies that the cosine of obslong-sunlong is 
! greater then 1 or less then -1.


  tempRatio = (cos(alpha) - sin(obslat) * sin(sunlat)) / &
       & (cos(obslat) * cos(sunlat))
       
  if ( abs(tempRatio) .gt. 1 ) then
      write(*,*) '**unable to find the azimuthal angle**'
      write(*,*) '**between incident and emitted ray**'
      write(*,*) '**check values for obslat,sunlat and phase angle***'
      stop
      
 ! **** need to terminate run here**** 
  endif
        
  obslong = acos( (cos(alpha) - sin(obslat) * sin(sunlat)) / &
       & (cos(obslat) * cos(sunlat)) ) + sunlong

  write(*,*)'sunlat, sunlong, obslat, obslong in Saturn frame:'
  write(*,*)sunlat/raddeg,sunlong/raddeg,obslat/raddeg,obslong/raddeg
  
  write(*,*)'Enter wake orientation angle in degrees:'
  read(*,*) phiwake
  
  phiwake = phiwake * raddeg

  sunlongwake = sunlong + phiwake
  obslongwake = obslong + phiwake

  write(*,*)'sunlat, sunlong, obslat, obslong in ring/wake frame:'
  write(*,*)sunlat/raddeg,sunlongwake/raddeg,obslat/raddeg,obslongwake/raddeg

! VECTOR POINTING TOWARDS SUN IN THE SATURN FRAME

    e_sunprime(1) = cos(sunlat) * cos(sunlong)
    e_sunprime(2) = cos(sunlat) * sin(sunlong)
    e_sunprime(3) = sin(sunlat)

! VECTOR POINTING TOWARDS SUN IN THE RING-WAKE FRAME

    e_sun(1) = cos(sunlat) * cos(sunlongwake)
    e_sun(2) = cos(sunlat) * sin(sunlongwake)
    e_sun(3) = sin(sunlat)

! VECTOR POINTING TOWARDS OBSERVER IN THE SATURN FRAME

    e_obsprime(1) = cos(obslat) * cos(obslong)
    e_obsprime(2) = cos(obslat) * sin(obslong)
    e_obsprime(3) = sin(obslat)

! VECTOR POINTING TOWARDS OBSERVER IN THE RING-WAKE FRAME

    e_obs(1) = cos(obslat) * cos(obslongwake)
    e_obs(2) = cos(obslat) * sin(obslongwake)
    e_obs(3) = sin(obslat)

  do
     write (*,*) 'Enter name of file containing particle positions'
     read (*,'(a)') infile
     inquire (file=infile, exist=test)
     if (test) exit
     write (*,'(/,a)') 'File not found. Check the file name and path.'
  end do
  
  open(10,file=infile,status='old',iostat=status)

  read(10,*)ctemp ! CONTAINS FILE PATH AND NAME
  read(10,*)xmax,ymax
  xmin = -xmax
  ymin = -ymax
  read(10,*)npart
  write(*,*)'Number of particles: ',npart

  do i = 1, npart
     read(10,*)rad_part(i),p_part(1:3,i)
  enddo

  close(10)

! CHECK THE DIMENSIONS OF THE PARTICLE FIELD IN THE Z DIRECTION
! (X, Y DIRRECTIONS LESS CRITICAL SINCE WILL USE GHOST PARTICLES)

  zmin = 9.9d99
  zmax = -9.9d99

  do i = 1, npart
     zmax = max(zmax, p_part(3,i) + rad_part(i))
     zmin = min(zmin, p_part(3,i) - rad_part(i))
  enddo

  write(*,*)'xmin,xmax: ',xmin,xmax
  write(*,*)'ymin,ymax: ',ymin,ymax
  write(*,*)'zmin,zmax: ',zmin,zmax

! CALCULATE THE AREA OF PARTICLE FIELD
  
  area = (xmax - xmin) * (ymax - ymin)

  call init_grid_speedup()

! CHOOSE PHASE FUNCTION AND CALCULATE FLUX TO OBSERVER

  write(*,*)'Enter ID nbr for particle phase function :'
  write(*,*)'3 = Power-law spherical phase function for Callisto (n=3.301)'
  write(*,*)'4 = Power-law spherical phase function for Europa'
  write(*,*)'6 = Shadowed-Minnaert surface reflection law - John table'
  write(*,*)'7 = Shadowed-Lambert surface reflection law'
  write(*,*)'8 = Shadowed-Lambert point particle phase function'  
  write(*,*)'these number codes are adopted for historical reasons'

  read(*,*)phasefunc

! FOR CALLISTO PHASE FUNCTION, SET THE EXPONENT N = 3.301, BEST FIT FOR 
! POWER-LAW TO CALLISTO FROM DONES ET AL. 1993, AND INITIALIZE LOOK-UP
! TABLE FOR NEW SCATTERING DIRECTION

  if (phasefunc == 3) then
     expo = 3.301

    if (idebug == 1) then
       open (30, file='debug.out', access='append')
       write (30,*)'calling init_power'
       close(30)
    endif

     call init_power ()

    if (idebug == 1) then
       open (30, file='debug.out', access='append')
       write (30,*)'back from call to init_power'
       close(30)
    endif

  end if

  if (phasefunc == 4) then
     expo = 2.400
     call init_power ()
  end if

  if (phasefunc == 6) then
     write(*,*)'Enter Minnaert parameter k to use:'
     read(*,*)k
     write(*,*)'Enter steepness parameter to use:'
     read(*,*)steep

  end if

  if (phasefunc == 7) then
     write(*,*)'Enter steepness parameter to use:'
     read(*,*)steep

  end if

  if (phasefunc == 8) then
     write(*,*)'Enter steepness parameter to use:'
     read(*,*)steep
     call init_power()
  end if

! SOLAR PHOTON SWITCH - DEFAULT IS ON - CAN TURN OFF TO SEE JUST EFFECTS
! OF SATURNSHINE AND/OR FOR DEBUGGING

  write(*,*) 'Solar photons on? 0 = off, 1 = on, default = 1'
  read(*,*) solarswitch
  
! FOR SATURNSHINE

  write(*,*) 'Saturnshine? 0 = off, 1 = on, default = 0'
  read(*,*) satswitch
  
  if (satswitch == 1) then
     write(*,*) 'Enter nbr of Saturnshine photons:'
     read(*,*) nphotsat

     write(*,*) 'Enter nbr of lat/long grids on facing Saturn hemisphere:'
     write(*,*) '(min should be 10 x 10)'
     read(*,*) nsatlat, nsatlong

     write(*,*) 'Enter distance of ring calculation area from Saturn center in km:'
     read(*,*) ring_dist

     ring_dist = ring_dist * 1.0d5   ! convert to cm
  end if
  
  if (phasefunc == 6 .or. phasefunc == 7) then
    call albedo()
  endif  

  return
end subroutine input

!--------------------------------------------------------------------
! INSIDE_RING_PLANE(P,E,FLAG):
! INSIDE RING PLANE - RETURNS TRUE IF PHOTON STILL INSIDE RING PLANE
! RETURNS FALSE IF OUTSIDE RING PLANE
!--------------------------------------------------------------------

subroutine inside_ring_plane (p,e,flag)
  use variables
  implicit none

  real*8 :: p(3),e(3) ! PHOTON POSITION AND DIRECTION

  logical :: flag

  flag = .true.
  if (p(3) > zmax.and.e(3) >= 0) flag = .false.
  if (p(3) < zmin.and.e(3) <= 0) flag = .false.

  return
end subroutine inside_ring_plane

!--------------------------------------------------------------------
! MAKE_GHOST_PHOTON(P,E):
! MAKE GHOST PHOTON - IF PHOTON MOVES OFF THE EDGE OF THE PARTICLE
! FIELD IN THE X OR Y DIRECTION BUT IS STILL IN RING PLANE, RESTART
! PHOTON FROM OPPOSITE EDGE OF PARTICLE FIELD
!--------------------------------------------------------------------

subroutine make_ghost_photon (p,e)
  use variables
  implicit none

  real*8 :: p(3),e(3),t,tnext ! PHOTON POSITION AND DIRECTION, AND TIME
                              ! TO NEXT COLLISION
  real*8 :: tnextx,tnexty

  character*4 :: edge ! WHICH EDGE PHOTON IS EXITING

  tnext = 9.9d99
  tnextx = 9.9d99
  tnexty = 9.9d99

! FIND TIME AT WHICH PHOTON WILL CROSS OUT OF PARTICLE FIELD IN X DIRECTION

   if (idebug == 1) then
       open (30, file='debug.out', access='append')
       write (30,*)'Inside make_ghost_phot:'
       write (30,*)'p=',p
       write (30,*)'e=',e
       write (30,*)'xmin,xmax=',xmin,xmax
       write (30,*)'ymin,ymax=',ymin,ymax
       close(30)
    endif

  if (e(1) > 0.0d0 ) then
     t = (xmax - p(1)) / e(1)
     if (t < tnext) then
        tnext = t
        tnextx = tnext
        edge = 'xmax'
     end if
  else if (e(1) < 0.d0) then
      t = (xmin  -  p(1)) / e(1)
      if (t < tnext) then
         tnext = t
         tnextx = tnext
         edge = 'xmin'
      end if
   else 
      edge = 'notx'
   end if

   if (idebug == 1) then
       open (30, file='debug.out', access='append')
       write (30,*)'x edge crossed: ',edge
       close(30)
    endif
   
   ! FIND TIME AT WHICH PHOTON WILL CROSS OUT OF PARTICLE FIELD IN Y DIRECTION
   
   if (e(2) > 0.d0) then
      t = (ymax  -  p(2)) / e(2)
      if (t < tnext) then
         tnext = t
         tnexty = tnext
         edge = 'ymax'
      end if
   else if (e(2) < 0.d0) then
      t = (ymin  -  p(2)) / e(2)
      if (t < tnext) then
         tnext = t
         tnexty = tnext
         edge = 'ymin'
      end if
   else
      edge = 'noty'
    end if

    if (idebug == 1) then
       open (30, file='debug.out', access='append')
       write (30,*)'y edge crossed: ',edge
       close(30)
    endif

    if (tnextx < tnexty) then
       tnext = tnextx
       if (e(1) > 0.d0) then
          edge = 'xmax'
       else if (e(1) < 0.d0) then
          edge = 'xmin'
       endif
    else if (tnextx > tnexty) then
       tnext = tnexty
       if (e(2) > 0.d0) then
          edge = 'ymax'
       else if (e(2) < 0.d0) then
          edge = 'ymin'
       endif
    endif

    if (idebug == 1) then
       open (30, file='debug.out', access='append')
       write (30,*)'soonest edge crossed: ',edge
       close(30)
    endif

! MOVE PHOTON TO THE EDGE OF THE PARTICLE FIELD

    p = p + e * tnext

! RESTART THE PHOTON FROM THE OPPOSITE EDGE OF THE PARTICLE FIELD

    if (idebug == 1) then
       open (30, file='debug.out', access='append')
       write (30,*)'photon pos before edge: ',p(1:3)
       close(30)
    endif

    if (edge == 'xmax') p(1) = xmin
    if (edge == 'xmin') p(1) = xmax
    if (edge == 'ymax') p(2) = ymin
    if (edge == 'ymin') p(2) = ymax

    if (idebug == 1) then
       open (30, file='debug.out', access='append')
       write (30,*)'photon pos after edge: ',p(1:3)
       close(30)
    endif

    return
end subroutine make_ghost_photon

!--------------------------------------------------------------------
! OUTPUT():
! OUTPUT
!--------------------------------------------------------------------

subroutine output ()
  use variables
  implicit none

  integer :: i,n,j
  character*80 :: fluxfile
  real*8 :: i_f ! I/F
  real*8 :: i_f_sat ! Satshine I/F
  real*8 :: i_f_sol ! Solar I/F  
  real*8 :: fac ! NORMALIZATION FACTOR
  real*8 :: fluxtot(nscatter_max) ! FLUX+FLUXSAT
  real*8 :: albedo,albedo_n ! ALBEDO AND ALBEDO TO POWER N
  real*8 :: albedo0
  
  open(20, file='flux.out', status='replace')
  write(20,*)k,steep,sunlat,obslat,nphot
  write(20,*)'i   solar flux    satshine flux      total flux'
  write(20,*)'-----------------------------------------------'

  do i = 1,nscatter_max
  
! WEIGHT SATSHINE FLUX AGAINST FLUX FROM DIRECT SOLAR PHOTONS
  
     fluxsat(i) = fluxsat(i) * ngridphottot / dble(nphotsat)
     write(20,*)i,flux(i),fluxsat(i),flux(i)+fluxsat(i)
  enddo
  close(20)
  
! I/F is calculated by computing the a weigthed sum of the flux across all orders
! for facet scattering the flux at order n is weighted by the curly-R (planar flux albedo
! for a smooth surface) raised to the nth power. For particle phase fuctions scattering 
! (such as Callisto) the flux of order n is weighted by the spherical albedo to the 
! nth power. Once the weighted sum is calculated it is multiplied by the conversion 
! factor appropriate for the viewing geometry. This is done in the following lines 
! of code.  
  
! CALCULATE FLUX NORMALIZATION FACTOR
  
  fac = pi * abs(sin(sunlat) / sin(obslat)) / dble(nphot)  ! same normalization for saturn?
  
  ! WRITE OUT TABLE OF I/F FOR SEVERAL ALBEDO VALUES
  
  if (phasefunc == 8) then
     albedo0 = 0.01d0/c_S(ishad)
  else
     albedo0 = 0.01d0
  endif   
  
  open(40, file='if.out', status='replace')
  
  if(phasefunc == 6 .or. phasefunc == 7) then
      write (40,*) 'Curly R     I/F(total)  I/F(solar)   I/F(saturn)'
  else
      write (40,*) 'albedo      I/F(total)  I/F(solar)   I/F(saturn)'
  endif    
  write (40,*) '---------------------------------------------------'
  do i = 1, 100
     albedo = albedo0 * dble(i)
     albedo_n = albedo
     i_f = 0.d0
     i_f_sat = 0.d0
     i_f_sol = 0.d0        
     do j = 1, nscatter_max
        i_f = i_f  +  albedo_n * (flux(j) + fluxsat(j))
        i_f_sol = i_f_sol  +  albedo_n * flux(j)
        i_f_sat = i_f_sat  +  albedo_n * fluxsat(j)
        albedo_n = albedo_n * albedo
     end do
     i_f = i_f * fac
     i_f_sol = i_f_sol * fac
     i_f_sat = i_f_sat * fac
     write (40,'(2x,f5.3,3x,ES11.4E2,3x,ES11.4E2,3x,ES11.4E2)') albedo,i_f,i_f_sol,i_f_sat
     
  end do
  
  close(40)
     
  open(60, file='nphotpass.out', status='replace')

  write(60,*)'nbr of photons that pass through rings without scattering off of any ring particles:'
  write(60,*) 'total: ', nphotpass
  write(60,*) 'nbr of noninteracting Satshine photons : ', nphotpassSat  
  write(60,*) 'nbr of noninteracting Solar photons : ', nphotpassSolar    
  write(60,*) 'nbr of Satshine photons used for flux: ', numfluxSat
  write(60,*) 'nbr of Solar photons used for flux: ', numfluxSolar  
  write(60,*) 'number of self blocked scattering events : ', numrejects
  write(60,*) 'number of scattering events hidden by another particle : ',numObstructed

  close(60)
  return
end subroutine output

!--------------------------------------------------------------------
! RAY_TRACING():
! RAY TRACING
!--------------------------------------------------------------------

subroutine ray_tracing ()
  use variables
  implicit none

  integer :: inext, nscatter ! COUNTERS
  integer :: iphot
  real*8 :: p(3), e(3), tnext ! P = PHOTON POSITION VECTOR
                              ! E = PHOTON DIRECTION VECTOR
                              ! TNEXT = TIME TO NEXT PARTICLE FOR 
                              ! SCATTERING
  real*8 :: dp(3),ddp
  integer :: iscatter

  inext = 0
  tnext = 0.0d0

  write(*,*)'Begin ray tracing of SOLAR photons...'

  do iphot = 1, nphot

     if(mod(iphot,1000).eq.0) then
        write(*,*) 'Solar iphot =',iphot
     endif

     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write (30,*)
        write (30,*) 'photon ',iphot
        close(30)
     endif
 
     call init_phot_pos(p,e)
     nscatter = 0
     
! initialize all weights to 1.0     
     do iscatter =1,nscatter_max
        weight(iscatter) = 1.0d0
     enddo   
     
     do while(nscatter.le.nscatter_max)
! DETERMINE WHICH PARTICLE THE PHOTON WILL HIT NEXT

        call which_part_next (p,e,inext,tnext,nscatter)
        
        if (idebug == 1) then
           open (30, file='debug.out', access='append')
           write (30,*)'back from which_part_next, inside ray_tracing'
           write (30,*) '   hit:',inext
           close(30)
        endif

! IF NO MORE PARTICLE HITS, GO TO NEXT PHOTON

        if (inext < 0) exit

! IF SCATTERING OCCURES, COUNT IT AND MOVE PHOTON TO COLLSION POINT
        
        nscatter = nscatter + 1
        p = p + tnext * e
        
        if (idebug == 1) then
           open (30, file='debug.out', access='append')
           write (30,*) 'Photon position just before scatter: ',p(1:3)
           write (30,*) 'Particle position: ',p_part(1:3,inext)
           write (30,*) 'Photon direction just before scatter: ',e(1:3)
           close(30)
        endif

! CALCULATE FLUX SEEN BY OBSERVER

        call flux_to_obs (p,e,inext,nscatter,iphot)

        if (idebug == 1 ) then
           open (30, file='debug.out', access='append')
           write(30,*) 'Back from flux_to_obs, inside ray_tracing'
           close(30)
        endif

!        hit = .false. ! RESET HIT FLAG FOR NEXT SCATTERING EVENT

! SCATTER PHOTON OFF NEXT PARTICLE
        
        if (idebug == 1) then
           open (30, file='debug.out', access='append')
           write (30,*) 'Direction before scatter: ',e(1:3)
           close(30)
        endif

        call scatter_phot (p,e,inext,tnext,nscatter)
        
        if (idebug == 1) then
           open (30, file='debug.out', access='append')
           write (30,*) 'Direction after  scatter: ',e(1:3)
           close(30)
        endif
     
     enddo
     intnlPhotons = intnlPhotons + 1
  enddo
    
  return
end subroutine ray_tracing

!--------------------------------------------------------------------
! RUF_MINNAERT(COSI,COSE,PHI):
! QUICK AND DIRTY ROUTINE TO CALCULATE MU, PHI PAIRS
! GIVEN MU0 USING THE MINNAERT SCATTERING LAW.  MU, PHI VALUES ARE
! CALCULATED ON A GRID USING A LOOK-UP TABLE WITH DIMENSION NDIM
! IN EACH DIRECTION, PLUS THE MU0 DIRECTION.  THE LOOK-UP TABLE IS
! STORED IN A 3D ARRAY CALLED PROB.  THE TABLE IS CALCULATED ONLY
! ONCE ON THE FIRST CALL TO THE SUBROUTINE.
! MU AND PHI ARE CALCULATED USING UNIFORM RANDOM VARIBALES WITHIN
! EACH GRID CELL IN THE LOOK-UP TABLE.
! N.B. IF NDIM IS LARGER THAN 128, ADD MORE ENTRIES TO THE ARRAY MM
! (ADDITIONAL POWERS OF 2 IN REVERSE ORDER)
! WRITTEN BY JOHN E. CHAMBERS (2010)
!--------------------------------------------------------------------

subroutine ruf_minnaert (cosi,cose,phi)
  use variables
  implicit none

    integer,parameter::NDIM = 32
    integer,parameter::NMM = 14

    integer::i,imu0,imu,iphi,itest,idum_local,entry,imax
    integer::mm(NMM)=(/8192,4096,2048,1024,512,256,128,64,32,16,8,4,2,1/)
    
    real*8::cosi,cose,phi,calpha,x,totprob,temp,lo,hi,mid
    real*8::mu0_mid(NDIM),mu_mid(NDIM),phi_mid(NDIM)
    real*8::mu0_lo(NDIM), mu_lo(NDIM), phi_lo(NDIM)
    real*8::mu0_hi(NDIM), mu_hi(NDIM), phi_hi(NDIM)
    real*8::romm0s(NDIM),romms(NDIM),cphi(NDIM)
    real*8::prob(NDIM,NDIM,NDIM)
    real*8::talpha2 ! TAN(ALPHA/2)
    real*8::sta2 ! SQRT[TAN{ALPHA/2)]
    real::ran1
    real*8 :: alpha1
    
    integer :: time_ruf
    real*8 :: frac_ruf
    
    

    logical::firstflag = .true.

    save mu0_lo,mu0_hi,mu_lo,mu_hi,phi_lo,phi_hi
    save prob,idum_local,imax,firstflag

! IF FIRST CALL TO SUBROUTINE, SET UP THE LOOK-UP TABLE
    if (firstflag) then
      do i = 1, NDIM
        mid = (dble(i)  -  0.5d0) / dble(NDIM)
        lo  = (dble(i)  -  1.d0)  / dble(NDIM)
        hi  =  dble(i)  / dble(NDIM)
        mu0_lo(i) = lo
        mu_lo(i)  = lo
        phi_lo(i) = lo * PI
        mu0_hi(i) = hi
        mu_hi(i)  = hi
        phi_hi(i) = hi * PI
        mu0_mid(i) = mid
        mu_mid(i)  = mid
        phi_mid(i) = mid * PI
        romm0s(i)  = sqrt(1.d0  -  mu0_mid(i) * mu0_mid(i))
        romms(i)   = sqrt(1.d0  -  mu_mid(i)  * mu_mid(i))
        cphi(i)    = cos(phi_mid(i))
      end do

      do imu0 = 1, NDIM
        totprob = 0.d0

        do imu = 1, NDIM
          do iphi = 1, NDIM
            calpha = romm0s(imu0) * romms(imu) * cphi(iphi) &
                   + mu_mid(imu) * mu0_mid(imu0)
            alpha1 = acos(calpha)
            talpha2 = tan(alpha1/2.d0)
            sta2 = sqrt(talpha2)
            
            totprob = totprob  +  mu_mid(imu)**k * exp(-steep*sta2)
            prob(imu0,imu,iphi) = totprob
          end do
        end do

! NORMALIZE THE PROBABILTIES FOR EACH VALUE OF MU0
        temp = 1.d0 / totprob
        prob(imu0,1:NDIM,1:NDIM) = prob(imu0,1:NDIM,1:NDIM) * temp
      end do

      idum_local = -1

!      call system_clock(time_ruf)
!      frac_ruf = real(time_temp)/10000
!      idum_local = -1*int(1000*(frac_ruf-real(int(frac_ruf))))

      imax = int(2.d0 * log10(dble(NDIM)) / log10(2.d0))
      firstflag = .false.
    end if

! STANDARD CALL TO SUBROUTINE FROM THIS POINT ON
!
! DETERMINE WHICH MU0 SLICE OF LOOK-UP TABLE TO USE
    imu0 = int(cosi * NDIM  +  1.d0)
    imu0 = min(imu0, NDIM)

! CHOOSE A RANDOM NUMBER
    x = dble(ran1(idum_local))

! FIND THE CORRESPONDING ENTRY IN THE LOOK-UP TABLE
    entry = 1
    do i = 1, NMM
      itest = entry  +  mm(i)
      if (itest > NDIM * NDIM) cycle
      imu  =     (itest - 2) / NDIM  +  1
      iphi = mod((itest - 2), NDIM)  +  1
      if (x > prob(imu0,imu,iphi)) entry = itest
    end do

! DETERMINE THE CORRESPONDING VALUES OF MU AND PHI
    imu  =     (entry - 1) / NDIM  +  1
    iphi = mod((entry - 1), NDIM)  +  1

    cose  = mu_lo(imu)   + (mu_hi(imu)   - mu_lo(imu))   * ran1(idum_local)
    phi = phi_lo(iphi) + (phi_hi(iphi) - phi_lo(iphi)) * ran1(idum_local) 

! PHI CAN BE POSITIVE OR NEGATIVE
    if (ran1(idum_local) < 0.5) phi = -phi

    return
end subroutine ruf_minnaert

!--------------------------------------------------------------------
! SATSHINE():
! SATURNSHINE - SETS UP GRID ON FACING SATURN SIDE, DETERMINES 
! WHICH GRID CELLS LIT AND VISIBLE TO RING CALC PATCH, CALCULATES
! BRIGHTNESS OF EACH GRID CELL, CALCULATES PROBABILITY OF PHOTON
! COMING FROM EACH GRID CELL
!--------------------------------------------------------------------

subroutine satshine ()
  use variables
  implicit none

  integer :: i,j ! COUNTERS
  integer :: isattab ! COUNTERS
  integer :: ilat, ilong, nkeep ! COUNTERS
  integer :: inext, nscatter ! COUNTERS
  
  logical :: flagvis
  logical :: flaglit

  real*8 :: gridw, gridh ! WIDTH AND HEIGHT OF EACH GRID CELL
  real*8 :: satlat, satlong ! LAT AND LONG OF MIDPOIT OF EACH
  ! SATURN GRID CELL
  real*8 :: mu0sat ! COS OF INCIDENT ANGLE OF SOLAR PHOTONS ONTO SATURN
  real*8 :: musat ! COS OF EMERGENT ANGLE OF SOLAR PHOTONS 
  ! SCATTERING FROM SATURN
  real*8 :: s(3) ! VECTOR FROM SATURN GRID CELL TO RING CALC PATCH IN 
  ! RING-WAKE FRAME
  real*8 :: sprime(3) ! VECTOR FROM SATURN GRID CELL TO RING CALC PATCH
  ! IN SATURN FRAME
  real*8 :: cos_iring ! COS OF INCIDENT ANGLE ONTO RING CALC PATCH OF PHOTONS
  ! FROM SATURN
  real*8 :: domega ! SOLID ANGLE THAT PT ON SATURN CASTS OUT TO RING CALC PATCH
  real*8 :: nprime(3) ! NORMAL VECTOR OF SATURN GRID PT
  real*8 :: ndotesun ! DOT PRODUCT N' AND E_SUN'
  real*8 :: ndots ! DOT PRODUCT N' AND S'
  real*8 :: s2, s2root ! DOT PRODUCT OF S' WITH ITSELF - FOR NORMALIZATION
  real*8 :: ifsat ! I/F BRIGHTNESS OF SATURN GRID CELL
  real*8 :: a, b ! COEFFICIENTS FOR BARKSTROM'S LAW, SEE DONES ET AL 1993
  real*8 :: alphasat ! PHASE ANGLE BETWEEN SOLAR INCIDENT RAY ONTO SATURN AND 
  ! EMERGENT RAY TO RING CALC PATCH
  real*8 :: alphasattab(6) ! TABLE OF PHASE ANGLES FOR USE WITH BARKSTROM'S LAW
  real*8 :: atab(6), btab(6) ! A AND B COEFFICIENTS FOR BARKSTROM'S LAW; DATA
  ! FROM DONES ET AL 1993
  real*8 :: alpha0, alpha1, a0, a1, b0, b1 ! VARIABLE TO DO LINEAR 
  ! INTERPOLATION OF TABLE
  real*8 :: ngridphotnorm(10000) ! ARRAY STORING NBR OF PHOTONS FROM EACH 
  ! GRID PT ON SATURN NORMALIZED
  real*8 :: ngridphotcum(10000) ! ARRAY STORING CUMMULATIVE NBR OF PHOTONS
  ! FROM SATURN GRID PTS
  real :: ran1 ! RANDOM NBR
  real*8 :: x ! RANDOM NBR
  real*8 :: p(3), e(3), tnext ! P = PHOTON POSITION VECTOR
                              ! E = PHOTON DIRECTION VECTOR
                              ! TNEXT = TIME TO NEXT PARTICLE FOR 
                              ! SCATTERING
  real*8 :: dp(3),ddp
  real*8 :: mk,ma ! MINNAERT K AND GEOMETRIC ALBEDO FOR SATURN
  integer :: iscatter
  
  

! CHOOSE PHASE FUNCTION FOR SATURN

  write (*,*) 'What phase func for Saturn?'
  write (*,*) '1 = Barkstrom Law with tabulated A & B'
  write (*,*) '2 = Barkstrom Law with user-chosen A & B'
  write (*,*) '3 = Minnaert'
  read (*,*) satphasefunc

  if (satphasefunc == 1) then

  ! SETTING UP TABLE LOOK-UP AND INTERPOLATION FROM DONES ET AL 1993
  ! VALUES OF A AND B COEFFICIENTS FROM BARKSTROM'S LAW

  alphasattab = (/0, 30, 60, 90, 120, 150/)
  do i = 1, 6
     alphasattab(i) = alphasattab(i) * pi / 180.0d0
  enddo

     write (*,*) 'Which brightness profile table for Saturnshine?'
     write (*,*) '1 = blue belts'
     write (*,*) '2 = blue zones'
     write (*,*) '3 = red belts'
     write (*,*) '4 = red zones'
     read (*,*) isattab

     if (isattab == 1) then
        atab = (/0.6736, 0.6197, 0.5726, 0.5948, 0.7869, 2.1323/)
        btab = (/1.1339, 1.1254, 1.1552, 1.1970, 1.2385, 1.4539/)
     else if (isattab == 2) then
        atab = (/0.5950, 0.5659, 0.5406, 0.5290, 0.6015, 1.5886/)
        btab = (/1.0948, 1.0938, 1.1372, 1.1606, 1.1636, 1.3727/)
     else if (isattab == 3) then
        atab = (/1.6854, 1.5946, 1.4529, 1.3386, 1.3653, 2.2539/)
        btab = (/1.4812, 1.4779, 1.4586, 1.4224, 1.3637, 1.3385/)
     else if (isattab == 4) then
        atab = (/1.6911, 1.5995, 1.4567, 1.3436, 1.3680, 2.2021/)
        btab = (/1.4825, 1.4794, 1.4592, 1.4221, 1.3618, 1.3341/)
     endif
     
     gridw = pi / nsatlong
     gridh = pi / nsatlat

     nkeep = 0

! STARTING FROM -90 LAT, -90 LONG FINDING MIDPT OF EACH GRID CELL

     do ilat = 1, nsatlat
        satlat = -pi / 2.0d0 + gridh * ( dble(ilat) - 0.5 )
        
        do ilong = 1, nsatlong
           satlong = -pi / 2.0d0 + gridw * ( dble(ilong) - 0.5 )
           
           nprime(1) = cos(satlat) * cos(satlong)
           nprime(2) = cos(satlat) * sin(satlong)
           nprime(3) = sin(satlat)
           
           sprime(1) = ring_dist - satrad * nprime(1)
           sprime(2) = -satrad * nprime(2)
           sprime(3) = -satrad * nprime(3)
           
           ndotesun = dot_product(nprime,e_sunprime)
           ndots = dot_product(nprime,sprime)
           
           lambda = -satrad * nprime(3) / e_sunprime(3)
           rcross(1) = satrad * nprime(1) + lambda * e_sunprime(1)
           rcross(2) = satrad * nprime(2) + lambda * e_sunprime(2)
           rcrossxy = sqrt ( rcross(1) * rcross(1) + rcross(2) * rcross(2) )

           ! CHECK TO SEE IF THIS GRID PT IS LIT, IF NOT FIND NEW LAT/LONG PT
           ! CHECKS TO SEE IF LAT/LONG PT IS IN RING SHADOW - ONLY FOR B AND
           ! A RING DOES IT BLOCK SUNLIGHT (B RING: 92000 - 117580 KM; 
           ! A RING: 122170 - 136780 KM

           flagvis = .false.
           flaglit = .true.
           
           if (ndotesun.lt.0.0d0) then
              flaglit = .false.
           else
              if (nprime(3)*e_sunprime(3).lt.0) then
                 if ((rcrossxy.gt.92000.d5.and.rcrossxy.lt.117580.d5).or. &
                      & (rcrossxy.gt.122170.d5.and.rcrossxy.lt.136780.d5)) then
                    flaglit = .false.
                 endif
              endif
           endif

!           if (.not.flaglit) then
!              ! RECORD ONLY THOSE GRID CELLS IN SHADOW BOTH FROM 
!              ! RING AND AT LIMB
!              open (50, file='shadow.out', access='append')
!              write (50,*) satlat/raddeg, satlong/raddeg
!              close (50)
!           endif

           if (flaglit.and.idebug == 1) then
              open (30, file='debug.out', access='append')
              write(30,*) 'Lit'
              close(30)
           endif
           
           ! CHECK TO SEE IF THIS GRID PT IS VISIBLE TO RING CALC PATCH
           ! IF NOT FIND NEW LAT/LONG PT
           
           if (ndots.gt.0.0d0) flagvis = .true.
           if (flagvis.and.idebug == 1) then
              open (30, file='debug.out', access='append')
              write(30,*) 'Visible'
              close(30)
           endif

           if (flagvis.and.flaglit.and.idebug == 1) then
              open (30, file='debug.out', access='append')
              write(30,*) 'Both'
              close(30)
           endif
           
           ! IF GRID PT IS BOTH LIT AND VISIBLE TO RING CALC PATCH, 
           ! CALCULATE MU0 AND MU
           
           if (flaglit.and.flagvis) then
              nkeep = nkeep + 1
              
              s2 = dot_product(sprime,sprime)
              s2root = sqrt(s2)
              
              snormprime(1,nkeep) = sprime(1)/s2root
              snormprime(2,nkeep) = sprime(2)/s2root
              snormprime(3,nkeep) = sprime(3)/s2root
              
              musat = ndots/s2root
              mu0sat = ndotesun
              
              alphasat = acos ( dot_product(e_sunprime,sprime)/s2root )
              
              if (alphasat > alphasattab(6)) then
                 
                 ! FUDGE FACTOR FOR WHEN OUTSIDE BOUNDS OF DONES TABLE FOR 
                 ! BARKSTROM'S LAW
                 alphasat = alphasattab(6)
                 
              end if
              
              ! LINEAR INTERPOLATION OF DONES ET AL 1993 DATA FOR BARKSTROM LAW
              
              do i = 2,6
                 if (alphasat < alphasattab(i).and.&
                      & alphasat > alphasattab(i-1)) then
                    
                    alpha0 = alphasattab(i-1)
                    alpha1 = alphasattab(i)
                    a0 = atab(i-1)
                    a1 = atab(i)
                    b0 = btab(i-1)
                    b1 = btab(i)
                    
                    a = (alphasat - alpha0) / (alpha1 - alpha0) * &
                         & (a1 - a0) + a0
                    b = (alphasat - alpha0) / (alpha1 - alpha0) * &
                         & (b1 - b0) + b0
                    
                 end if
                 
              enddo
              
              do i = 1,6
                 if (alphasat == alphasattab(i)) then
                    a = atab(i)
                    b = btab(i)
                 end if
              enddo
              
              ! CALCULATE I/F FOR EACH PT ON SATURN
              
              ifsat = a / musat * (mu0sat * musat / (mu0sat +  musat) ) ** b

! TURNING OFF N. HEMISPHERE AS HACK FOR RING SHADOWING EFFECT
!
!              if (satlat.gt.0) then
!                 ifsat = a / musat * (mu0sat * musat / &
!                      & (mu0sat +  musat) ) ** b
!              else
!                 ifsat = 0.0d0
!              end if

              ! CALCULATE NBR OF PHOTONS COMING FROM EACH GRID PT ON SATURN 
              ! RELATIVE TO THE NBR OF PHOTONS COMING DIRECTLY FROM THE SUN
              
              cos_iring = satrad * sin(satlat) / sqrt( satrad * satrad + & 
                & ring_dist * ring_dist - 2.0d0 * ring_dist * satrad * &
                & cos(satlong) * cos(satlat) )
              
              domega = gridw * satrad * gridh * satrad * cos(satlat) * &
                   & musat / s2
              
              ngridphot(nkeep) = dble(nphot) * ifsat / pi * abs(cos_iring) &
                   & * domega / abs(sin(sunlat))
              
              ngridphottot = ngridphottot + ngridphot(nkeep)
              
              if (idebug == 1) then
                 open (30, file='debug.out', access='append')
                 write(30,*) 'nkeep:', nkeep
                 write(30,*) 'ifsat:', ifsat
                 write(30,*) 'musat, mu0sat:',musat,mu0sat
                 write(30,*) 'a,b:', a,b
                 write(30,*) 'cos_iring:', cos_iring
                 write(30,*) 'domega:', domega
                 write(30,*) 'ngridphot(nkeep):',ngridphot(nkeep)
                 write(30,*) 'nphot:',nphot
                 close(30)
              endif
                 
           endif
        enddo
     enddo
        
  endif

  if (satphasefunc == 2) then

     write (*,*) 'Barkstrom A ~ Saturn albedo'
     read (*,*) a
     write (*,*) 'Barkstrom B ~ amount of backscatter'
     read (*,*) b

     gridw = pi / nsatlong
     gridh = pi / nsatlat

!     open (50, file='shadow.out')
!     write (50,*) gridw/raddeg, gridh/raddeg
!     close (50)

     nkeep = 0

! STARTING FROM -90 LAT, -90 LONG FINDING MIDPT OF EACH GRID CELL

     do ilat = 1, nsatlat
        satlat = -pi / 2.0d0 + gridh * ( dble(ilat) - 0.5 )
        
        do ilong = 1, nsatlong
           satlong = -pi / 2.0d0 + gridw * ( dble(ilong) - 0.5 )
           
           nprime(1) = cos(satlat) * cos(satlong)
           nprime(2) = cos(satlat) * sin(satlong)
           nprime(3) = sin(satlat)
           
           sprime(1) = ring_dist - satrad * nprime(1)
           sprime(2) = -satrad * nprime(2)
           sprime(3) = -satrad * nprime(3)
           
           ndotesun = dot_product(nprime,e_sunprime)
           ndots = dot_product(nprime,sprime)
           
           lambda = -satrad * nprime(3) / e_sunprime(3)
           rcross(1) = satrad * nprime(1) + lambda * e_sunprime(1)
           rcross(2) = satrad * nprime(2) + lambda * e_sunprime(2)
           rcrossxy = sqrt ( rcross(1) * rcross(1) + rcross(2) * rcross(2) )

           ! CHECK TO SEE IF THIS GRID PT IS LIT, IF NOT FIND NEW LAT/LONG PT
           ! CHECKS TO SEE IF LAT/LONG PT IS IN RING SHADOW - ONLY FOR B AND
           ! A RING DOES IT BLOCK SUNLIGHT (B RING: 92000 - 117580 KM; 
           ! A RING: 122170 - 136780 KM

           flagvis = .false.
           flaglit = .true.
           
           if (ndotesun.lt.0.0d0) then
              flaglit = .false.
           else
              if (nprime(3)*e_sunprime(3).lt.0) then
                 if ((rcrossxy.gt.92000.d5.and.rcrossxy.lt.117580.d5).or. &
                      & (rcrossxy.gt.122170.d5.and.rcrossxy.lt.136780.d5)) then
                    flaglit = .false.
                 endif
              endif
           endif

!           if (.not.flaglit) then
!              ! RECORD ONLY THOSE GRID CELLS IN SHADOW BOTH FROM 
!              ! RING AND AT LIMB
!              open (50, file='shadow.out', access='append')
!              write (50,*) satlat/raddeg, satlong/raddeg
!              close (50)
!           endif

           if (flaglit.and.idebug == 1) then
              open (30, file='debug.out', access='append')
              write(30,*) 'Lit'
              close(30)
           endif
           
           ! CHECK TO SEE IF THIS GRID PT IS VISIBLE TO RING CALC PATCH
           ! IF NOT FIND NEW LAT/LONG PT
           
           if (ndots.gt.0.0d0) flagvis = .true.
           if (flagvis.and.idebug == 1) then
              open (30, file='debug.out', access='append')
              write(30,*) 'Visible'
              close(30)
           endif
           if (flagvis.and.flaglit.and.idebug == 1) then
              open (30, file='debug.out', access='append')
              write(30,*) 'Both'
              close(30)
           endif
           
           ! IF GRID PT IS BOTH LIT AND VISIBLE TO RING CALC PATCH, 
           ! CALCULATE MU0 AND MU
           
           if (flaglit.and.flagvis) then
              nkeep = nkeep + 1
              
              s2 = dot_product(sprime,sprime)
              s2root = sqrt(s2)
              
              snormprime(1,nkeep) = sprime(1)/s2root
              snormprime(2,nkeep) = sprime(2)/s2root
              snormprime(3,nkeep) = sprime(3)/s2root
              
              musat = ndots/s2root
              mu0sat = ndotesun
              
              alphasat = acos ( dot_product(e_sunprime,sprime)/s2root )
              
              ! CALCULATE I/F FOR EACH PT ON SATURN
              
              ifsat = a / musat * (mu0sat * musat / (mu0sat +  musat) ) ** b

! TURNING OFF N. HEMISPHERE AS HACK FOR RING SHADOWING EFFECT
!
!              if (satlat.gt.0) then
!                 ifsat = a / musat * (mu0sat * musat / &
!                      & (mu0sat +  musat) ) ** b
!              else
!                 ifsat = 0.0d0
!              end if
              
              ! CALCULATE NBR OF PHOTONS COMING FROM EACH GRID PT ON SATURN 
              ! RELATIVE TO THE NBR OF PHOTONS COMING DIRECTLY FROM THE SUN
              
              cos_iring = satrad * sin(satlat) / sqrt( satrad * satrad + & 
                & ring_dist * ring_dist - 2.0d0 * ring_dist * satrad * &
                & cos(satlong) * cos(satlat) )
              
              domega = gridw * satrad * gridh * satrad * cos(satlat) * &
                   & musat / s2
              
              ngridphot(nkeep) = dble(nphot) * ifsat / pi * abs(cos_iring) &
                   & * domega / abs(sin(sunlat))
              
              ngridphottot = ngridphottot + ngridphot(nkeep)
              
              if (idebug == 1) then
                 open (30, file='debug.out', access='append')
                 write(30,*) 'nkeep:', nkeep
                 write(30,*) 'ifsat:', ifsat
                 write(30,*) 'musat, mu0sat:',musat,mu0sat
                 write(30,*) 'a,b:', a,b
                 write(30,*) 'cos_iring:', cos_iring
                 write(30,*) 'domega:', domega
                 write(30,*) 'ngridphot(nkeep):',ngridphot(nkeep)
                 write(30,*) 'nphot:',nphot
                 close(30)
              endif              
           endif
        enddo
     enddo
    
  endif

  if (satphasefunc == 3) then

     write (*,*) 'Minnaert k?'
     read (*,*) mk
     write (*,*) 'Minnaert albedo for Saturn'
     read (*,*) ma

     gridw = pi / nsatlong
     gridh = pi / nsatlat

!     open (50, file='shadow.out')
!     write (50,*) gridw/raddeg, gridh/raddeg
!     close (50)

     nkeep = 0

! STARTING FROM -90 LAT, -90 LONG FINDING MIDPT OF EACH GRID CELL

     do ilat = 1, nsatlat
        satlat = -pi / 2.0d0 + gridh * ( dble(ilat) - 0.5 )
        
        do ilong = 1, nsatlong
           satlong = -pi / 2.0d0 + gridw * ( dble(ilong) - 0.5 )
           
           nprime(1) = cos(satlat) * cos(satlong)
           nprime(2) = cos(satlat) * sin(satlong)
           nprime(3) = sin(satlat)
           
           sprime(1) = ring_dist - satrad * nprime(1)
           sprime(2) = -satrad * nprime(2)
           sprime(3) = -satrad * nprime(3)
           
           ndotesun = dot_product(nprime,e_sunprime)
           ndots = dot_product(nprime,sprime)

           lambda = -satrad * nprime(3) / e_sunprime(3)
           rcross(1) = satrad * nprime(1) + lambda * e_sunprime(1)
           rcross(2) = satrad * nprime(2) + lambda * e_sunprime(2)
           rcrossxy = sqrt ( rcross(1) * rcross(1) + rcross(2) * rcross(2) )

           ! CHECK TO SEE IF THIS GRID PT IS LIT, IF NOT FIND NEW LAT/LONG PT
           ! CHECKS TO SEE IF LAT/LONG PT IS IN RING SHADOW - ONLY FOR B AND
           ! A RING DOES IT BLOCK SUNLIGHT (B RING: 92000 - 117580 KM; 
           ! A RING: 122170 - 136780 KM

           flagvis = .false.
           flaglit = .true.
           
           if (ndotesun.lt.0.0d0) then
              flaglit = .false.
           else
              if (nprime(3)*e_sunprime(3).lt.0) then
                 if ((rcrossxy.gt.92000.d5.and.rcrossxy.lt.117580.d5).or. &
                      & (rcrossxy.gt.122170.d5.and.rcrossxy.lt.136780.d5)) then
                    flaglit = .false.
                 endif
              endif
           endif

           if (.not.flaglit) then
              ! RECORD ONLY THOSE GRID CELLS IN SHADOW BOTH FROM 
              ! RING AND AT LIMB
!              open (50, file='shadow.out', access='append')
!              write (50,*) satlat/raddeg, satlong/raddeg
!              close (50)
           endif

           if (flaglit.and.idebug == 1) then
              open (30, file='debug.out', access='append')
              write(30,*) 'Lit'
              close(30)
           endif
           
           ! CHECK TO SEE IF THIS GRID PT IS VISIBLE TO RING CALC PATCH
           ! IF NOT FIND NEW LAT/LONG PT
           
           if (ndots.gt.0.0d0) flagvis = .true.
           if (flagvis.and.idebug == 1) then
              open (30, file='debug.out', access='append')
              write(30,*) 'Visible'
              close(30)
           endif
           if (flagvis.and.flaglit.and.idebug == 1) then
              open (30, file='debug.out', access='append')
              write(30,*) 'Both'
              close(30)
           endif
           
           ! IF GRID PT IS BOTH LIT AND VISIBLE TO RING CALC PATCH, 
           ! CALCULATE MU0 AND MU
           
           if (flaglit.and.flagvis) then
              nkeep = nkeep + 1
              
              s2 = dot_product(sprime,sprime)
              s2root = sqrt(s2)
              
              snormprime(1,nkeep) = sprime(1)/s2root
              snormprime(2,nkeep) = sprime(2)/s2root
              snormprime(3,nkeep) = sprime(3)/s2root
              
              musat = ndots/s2root
              mu0sat = ndotesun
              
              alphasat = acos ( dot_product(e_sunprime,sprime)/s2root )
              
              ! CALCULATE I/F FOR EACH PT ON SATURN
              
              ifsat = (mk + 0.5) * ma * mu0sat ** mk * musat ** (mk - 1)

! TURNING OFF N. HEMISPHERE AS HACK FOR RING SHADOWING EFFECT
!
!              if (satlat.gt.0) then
!              ifsat = (mk + 0.5) * ma * mu0sat ** mk &
!              & * musat ** (mk - 1)
!              else
!                 ifsat = 0.0d0
!              end if
              
              ! CALCULATE NBR OF PHOTONS COMING FROM EACH GRID PT ON SATURN 
              ! RELATIVE TO THE NBR OF PHOTONS COMING DIRECTLY FROM THE SUN
              
              cos_iring = satrad * sin(satlat) / sqrt( satrad * satrad + & 
                & ring_dist * ring_dist - 2.0d0 * ring_dist * satrad * &
                & cos(satlong) * cos(satlat) )
              
              domega = gridw * satrad * gridh * satrad * cos(satlat) * &
                   & musat / s2
              
              ngridphot(nkeep) = dble(nphot) * ifsat / pi * abs(cos_iring) &
                   & * domega / abs(sin(sunlat))
              
              ngridphottot = ngridphottot + ngridphot(nkeep)
              
              if (idebug == 1) then
                 open (30, file='debug.out', access='append')
                 write(30,*) 'nkeep:', nkeep
                 write(30,*) 'ifsat:', ifsat
                 write(30,*) 'musat, mu0sat:',musat,mu0sat
                 write(30,*) 'cos_iring:', cos_iring
                 write(30,*) 'domega:', domega
                 write(30,*) 'ngridphot(nkeep):',ngridphot(nkeep)
                 write(30,*) 'nphot:',nphot
                 close(30)
              endif
              
           endif
        enddo
     enddo
    
  endif
 
  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write(30,*) 'Total sat photons:', ngridphottot
     close(30)
  endif

  ! CALCULATE ARRAYS FOR NBR OF PHOTONS FROM EACH GRID PT ON SATURN 
  ! NORMALIZED AND CUMMULATIVE

  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write(30,*)'i,ngridphot(i),ngridphotnorm(i),ngridphotcum(i):'
     close(30)
  endif

  cosphiwake = cos (phiwake)
  sinphiwake = sin (phiwake)

  do i = 1, nkeep
     ngridphotnorm(i) = ngridphot(i)/ngridphottot
     
     if (i == 1) then
        ngridphotcum(i) = ngridphotnorm(i)
     else
        ngridphotcum(i) = ngridphotcum(i-1) + ngridphotnorm(i)
     endif
     
     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*) i,ngridphot(i),ngridphotnorm(i),ngridphotcum(i)
        close(30)
     endif

     ! NORMALIZE SNORM VECTOR INTO RING-WAKE FRAME TO USE WITH REST OF 
     ! PROGRAM
     
     snorm(1,i) = snormprime(1,i) * cosphiwake + snormprime(2,i) * sinphiwake
     snorm(2,i) = snormprime(2,i) * cosphiwake - snormprime(1,i) * sinphiwake
     snorm(3,i) = snormprime(3,i)
     
  enddo

  ! FIRE PHOTONS FROM SATURN INTO RING PATCH

  inext = 0
  tnext = 0.0d0

  write(*,*)'Begin ray tracing of SATSHINE photons...'

  do i = 1, nphotsat

     x = dble( ran1(idum) )

     ! DETERMINE WHICH GRID PT PHOTON SHOULD COME FROM

     iphotbin = 0

     do j = 1, nkeep
        
        if ( x.gt.ngridphotcum(j) ) then
           iphotbin = j
        endif

     enddo

     iphotbin = iphotbin + 1

     ! START RAY TRACING FOR SATSHINE PHOTONS
     
     if(mod(i,1000).eq.0) then
        write(*,*) 'Satshine iphot =',i
     endif

     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write (30,*)
        write (30,*) 'satshine photon and iphotbin ',i,iphotbin
        close(30)
     endif

     call init_phot_pos(p,e)
     nscatter = 0
! initialize all weights to 1.0     
     do iscatter =1,nscatter_max
        weight(iscatter) = 1.0d0
     enddo   
       
     do while(nscatter.le.nscatter_max)
! DETERMINE WHICH PARTICLE THE PHOTON WILL HIT NEXT

        call which_part_next (p,e,inext,tnext,nscatter)
        if (idebug == 1) then
           open (30, file='debug.out', access='append')
           write (30,*)'back from which_part_next, inside satshine'
           write (30,*) '   hit:',inext
           close(30)
        endif

! IF NO MORE PARTICLE HITS, GO TO NEXT PHOTON

        if (inext < 0) exit

! IF SCATTERING OCCURES, COUNT IT AND MOVE PHOTON TO COLLSION POINT
        
        nscatter = nscatter + 1
        p = p + tnext * e

        if (idebug == 1) then
           open (30, file='debug.out', access='append')
           write (30,*) 'Photon position just before scatter: ',p(1:3)
           write (30,*) 'Particle position: ',p_part(1:3,inext)
           write (30,*) 'Photon direction just before scatter: ',e(1:3)
           close(30)
        endif

! CALCULATE FLUX SEEN BY OBSERVER

        call flux_to_obs (p,e,inext,nscatter,i)

!        hit = .false. !RESET HIT FLAG FOR NEXT SCATTERING EVENT

! SCATTER PHOTON OFF NEXT PARTICLE

        if (idebug == 1) then
           open (30, file='debug.out', access='append')
           write(30,*) 'Back from flux_to_obs, inside satshine'
           write (30,*) 'Direction before scatter: ',e(1:3)
           close(30)
        endif

       call scatter_phot (p,e,inext,tnext,nscatter)

       if (idebug == 1) then
          open (30, file='debug.out', access='append')
          write (30,*) 'Direction after  scatter: ',e(1:3)
          close(30)
       endif

     enddo     

  enddo

! TURNING SATSHINE SWITCH OFF AFTER STORING PHOTON AND RAY TRACING INFO ON 
! SATSHINE PHOTONS, SO REST OF PROGRAM RUNS WITH JUST DIRECT SOLAR PHOTONS

  satswitch = 0
  
  return
end subroutine satshine

!--------------------------------------------------------------------
! SCATTER_PHOT(P,E,INEXT,TNEXT):
! SCATTER PHOTON - DETERMINE NEW DIRECTION OF PHOTON 
! AFTER SCATTERING FROM PARTICLE
!--------------------------------------------------------------------

subroutine scatter_phot (p,e,inext,tnext,nscatter)
  use variables
  implicit none
  
  integer :: i,j
  integer :: inext
  
  real*8 :: p(3),e(3),tnext ! PHOTON POSITION AND DIRECTION, AND TIME
  ! TO NEXT COLLISION
  real*8 :: phi,cose,sine,cosphi,sinphi ! PHI = AZIMUTH OF NEW PHOTON 
  ! DIRECTION
  ! COSE, SINE = COS AND SIN OF EMERGENT ANGLE
  ! COSPHI, SINPHI = COS AND SIN OF AZIMUTHAL ANGLE
  real*8 :: eold(3),eold2 ! COPY OF PRE-SCATTERED PHOTON DIRECTION
  real*8 :: psi,cospsi,sinpsi ! PSI = AZIMUTH OF INCOMING PHOTON DIRECTION
  real*8 :: n(3),n2 ! UNIT NORMAL VECTOR
  real*8 :: root,root2 
  real :: ran1 ! RANDOM NUMBER
  real*8 :: temp
  real*8 :: cosi,sini ! COS AND SIN OF INCIDENT ANGLE
  real*8 :: alp ! ALPHA FROM POWER-LAW TABLE INTERPOLATION
  real*8 :: enew(3),enew2 ! COPY OF E POST-SCATTERING
  integer :: nscatter
  real*8 :: cosalp, sinalp

! PHASEFUNC = 7: USING LAMBERT SURFACE ELEMENT REFLECTION LAW

  if (phasefunc == 7) then
     phi = 2.d0 * pi * dble( ran1(idum) )
     cosphi = cos(phi)
     sinphi = sin(phi)
     sine = sqrt( dble( ran1(idum) ) )
     cose = sqrt(1.d0  -  sine * sine)
     
! COMPONENTS OF NORMAL FORM PARTICLE CENTER TO PHOTON

     n(1:3) = p(1:3)  -  p_part(1:3,inext)
     n2 = dot_product(n,n)
     n = n / sqrt(n2)

! NEW PHOTON DIRECTION

     if (n(3) * n(3) >= 1.d0) then
        e(1) = sine * cosphi
        e(2) = sine * sinphi
        e(3) = n(3) * cose
     else
        root2 = 1.d0  -  n(3) * n(3)
        root = sqrt(root2)
        e(3) = n(3) * cose  +  sine * cosphi * root
        e(1) = (n(1)*(cose - n(3)*e(3))  -  &
             & sine * sinphi * n(2) * root) / root2
        e(2) = (n(2)*(cose - n(3)*e(3))  +  &
             & sine * sinphi * n(1) * root) / root2
     end if
     
  end if
 
! PHASEFUNC = 3: USING CALLISTO POWER-LAW SPHERICAL PHASE FUNCTION
! PHASEFUNC = 4: USING EUROPA POWER-LAW SPHERICAL PHASE FUNCTION
! PHASEFUNC = 8: USING SHADOWED LAMBERT PARTICLE PHASE FUNCTION

  if (phasefunc == 3.or.phasefunc == 4.or.phasefunc==8) then
     
! COMPONENTS OF INCOMING RAY
     
     eold = -e 
     eold2 = dot_product(eold,eold)
     eold = eold / sqrt(eold2)

     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*)'eold = ',eold(1:3)
        close(30)
     endif
     
200  continue
     psi = 2.d0 * pi * dble( ran1(idum) )
     cospsi = cos(psi)
     sinpsi = sin(psi)
 
     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*)'psi, cospsi, sinpsi: ',psi,cospsi,sinpsi
        close(30)
     endif

! COMPONENTS OF NORMAL FORM PARTICLE CENTER TO PHOTON
     
     n(1:3) = p(1:3)  -  p_part(1:3,inext)
     n2 = dot_product(n,n)
     n = n / sqrt(n2)
     
     cosi = dot_product(-e,n)
     sini = sqrt (1.0d0 - cosi * cosi)
 
     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*) 'n = ',n
        write(30,*) 'cosi, sini: ',cosi,sini
        close(30)
     endif

     temp = dble( ran1(idum) )

     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*) 'temp = ',temp
        close(30)
     endif

     do i = 1, nbins+1

!        if (idebug == 1) then
!           open (30, file='debug.out', access='append')
!           write(30,*) 'i, power_rand(i) = ',i,power_rand(i)
!           close(30)
!       endif

        if (power_rand(i) > temp) then
          
! INTERPOLATE BETWEEN i-1 AND i VALUES OF POWER_ALPHA AND POWER_RAND

           alp = ( power_alpha(i) * (temp - power_rand(i-1)) - &
                & power_alpha(i-1) * (temp - power_rand(i)) ) / &
                & ( power_rand(i) - power_rand(i-1) )


                      
           if (idebug == 1) then
              open (30, file='debug.out', access='append')
              write(30,*)'alp, cosalp, sinalp: ',alp,cosalp,sinalp
              close(30)
           endif
          
! NEW PHOTON DIRECTION adapting Yu A. Shreider (1966) pgs 151-153 

! for these equations we really need scattering angle not phase angle 
		  cosalp = cos(pi-alp)
		  sinalp = sin(pi-alp)


           if (eold(3) * eold(3) >= 1.d0) then
              e(1) = sinalp * cospsi
              e(2) = sinalp * sinpsi
              e(3) = eold(3) * cosalp
           else
              root2 = 1.d0  -  eold(3) * eold(3)
              root = sqrt(root2)
              e(3) = eold(3) * cosalp  +  sinalp * cospsi * root
              e(1) = (eold(1)*(cosalp - eold(3)*e(3))  -  &
                   & sinalp * sinpsi * eold(2) * root) / root2
              e(2) = (eold(2)*(cosalp - eold(3)*e(3))  +  & 
                   & sinalp * sinpsi * eold(1) * root) / root2             
           end if

           if (idebug == 1) then
              open (30, file='debug.out', access='append')
              write(30,*)'eold(3)^2 = ',eold(3) * eold(3)
              write(30,*)'new direction e = ',e(1:3)
              close(30)
           endif
        
           enew = e
           enew2 = dot_product(enew,enew)
           enew = enew / sqrt(enew2)
           if (idebug == 1) then
              open (30, file='debug.out', access='append')
              write(30,*)'enew = ',enew(1:3)
              write(30,*)'enew dot n = ',dot_product(enew,n)
              close(30)
           endif

           if (dot_product(enew,n) < 0) then
              go to 200
           end if
           
           exit
           
        end if
     enddo
          
  end if

! PHASEFUNC = 6: USING SHADOW-MINNAERT SURFACE ELEMENT REFLECTION LAW - 
! JOHN'S TABLE

  if (phasefunc == 6) then

     n(1:3) = p(1:3)  -  p_part(1:3,inext)
     n2 = dot_product(n,n)
     n = n / sqrt(n2)
     
     cosi = dot_product(-e,n)
     sini = sqrt (1.0d0 - cosi * cosi)

     call ruf_minnaert(cosi,cose,phi)
     
     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*)'final Shadow-Minnaert: cosi, cose, phi = ',cosi,cose,phi
        close(30)
     endif

     cosphi = cos(phi)
     sinphi = sin(phi)
     sine = sqrt(1.d0  -  cose * cose)

! COMPONENTS OF NORMAL FORM PARTICLE CENTER TO PHOTON

     n(1:3) = p(1:3)  -  p_part(1:3,inext)
     n2 = dot_product(n,n)
     n = n / sqrt(n2)

! NEW PHOTON DIRECTION

     if (n(3) * n(3) >= 1.d0) then
        e(1) = sine * cosphi
        e(2) = sine * sinphi
        e(3) = n(3) * cose
     else
        root2 = 1.d0  -  n(3) * n(3)
        root = sqrt(root2)
        e(3) = n(3) * cose  +  sine * cosphi * root
        e(1) = (n(1)*(cose - n(3)*e(3))  -  &
             & sine * sinphi * n(2) * root) / root2
        e(2) = (n(2)*(cose - n(3)*e(3))  +  & 
             & sine * sinphi * n(1) * root) / root2
     end if
     
  end if

  return
end subroutine scatter_phot

!--------------------------------------------------------------------
! SHAD_MINNAERT(COSI,COSE,PHI):
! SHADOW-MINNAERT PHASE FUNC CALCULATION - THIS SUBROUTINE WILL 
! GENERATE LOTS OF RANDOM VALUES OF MU (ALSO CALLED COSE) AND PHI
! USING UNIGRID_GENERATE AND THEN COMPARE THE MEAN OF THEM TO 
! AUTOMATICALLY-DETERMINED TRUE MEAN (COMPUTED FROM DISTRIBUTION)
!--------------------------------------------------------------------

subroutine shad_minnaert (cosi,cose,phi)
  use variables
  implicit none

  integer :: i,j,itot,i_mu,i_phi ! COUNTERS
  integer :: complete, complete_last

  real*8 :: cosi,cose,phi
  real*8 :: muint(2)
  real*8 :: sum_mu,sum_phi,mean_mu,mean_phi ! USED TO CALCULATE THE MEAN
  ! OF MU AND PHI VALUES GENERATED WITH UNIGRID_GENERATE
  real*8 :: x,y,dx,dy ! FOR 2D INTEGRATION
  real*8 :: phasefunc_scaled ! FUNCTION THAT CALCULATES SHADOW-MINNAERT
  ! PHASE FUNC
  real*8 :: realmoment_mu, realmoment_phi ! MEAN VALUES OF MU AND PHI AT MU_0

  nhist = nhistbins

  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write (30,*) 'Beginning test of numerical output of uniform grid &
          & generation for cosi =',cosi 
     write(30,*) 'Generating ',npts,' test points...'
     close(30)
  endif

! STORE OUTPUT OF UNIGRID_GENERATE

  itot = 0
  i_mu = 0
  i_phi = 0
  complete = 0
  muint(1) = 0.d0
  muint(2) = 1.d0 
  mu0np = cosi

  do i = 1,npts

     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*) 'Before call to unigrid_generate'
        write(30,*) 'cosi,cose,phi=',cosi,cose,phi
        close(30)
     endif

     call unigrid_generate(cosi,cose,phi)

     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*) 'After call to unigrid_generate'
        close(30)
     endif

     genpts(i,1) = cose
     genpts(i,2) = phi
     itot = itot + iter
     complete_last = complete
     complete = ceiling(i*1d1/npts)
     if(complete.gt.complete_last.and.idebug.eq.1) then
        write(30,*) complete*10,'%'
     endif
  enddo

  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write (30,*) 'generated ',npts,'points in ',itot,'total iterations'
     write (30,*) 'mean nbr iterations = ',itot*1d0/npts
     close(30)
  endif

! TESTING - DETERMINE MOMENTS AND COMPARE TO ANALYTICALLY-COMPUTED MOMENTS

  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write (30,*) 'computing moments of test values'
     close(30)
  endif

  sum_mu = 0.d0
  sum_phi = 0.d0
  mean_mu = 0.d0
  mean_phi = 0.d0

  do i=1,npts
     sum_mu = sum_mu + genpts(i,1)
     sum_phi = sum_phi + genpts(i,2)
  enddo

  mean_mu = sum_mu/dble(npts)
  mean_phi = sum_phi/dble(npts)

  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write (30,*) 'computing normalized constant & actual moments'
     close(30)
  endif

! 2D INTEGRATION OF PHASEFUNC

  normconst = 0.d0

  dx = (1.d0 - 0.d0) / dble(96)
  dy = (pi + pi) / dble(96)

  do i = 1,96
     do j = 1,96
        x = 0.d0 + dx * (dble(i) - 0.5d0)
        y = -pi + dy * (dble(j) - 0.5d0)
        normconst = normconst + phasefunc_scaled(mu0np,x,y)
     enddo
  enddo

  normconst = normconst * dx * dy

! 2D ITEGRATION TO GET MEAN VALUES OF MU AND PHI AT MU_0

  realmoment_mu = 0.d0
  realmoment_phi = 0.d0

  do i = 1,96
     do j = 1,96
        x = 0.d0 + dx * (dble(i) - 0.5d0)
        y = -pi + dy * (dble(j) - 0.5d0)
        realmoment_mu = realmoment_mu + cose*phasefunc_scaled(mu0np,x,y)
        realmoment_phi = realmoment_phi + phi*phasefunc_scaled(mu0np,x,y)
     enddo
  enddo
  
  realmoment_mu = realmoment_mu * dx * dy
  realmoment_phi = realmoment_phi * dx * dy

  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write (30,*) 'mean values of mu and phi at mu_0=',cosi,&
          & ' are ',realmoment_mu,realmoment_phi
     write(30,*) 'mean values of test points generated are:',mean_mu,mean_phi
     write(30,*) 'Mean testing complete.  &
          & Beginning density-of-points histogram...'
     close(30)
  endif

  do i = 0,npts-1
     i_mu = ceiling(genpts(i,1)*nhistbins)
     i_phi = ceiling( (genpts(i,2)+pi) / (2.d0*pi) * nhistbins )
     histcounts(i_mu,i_phi) = histcounts(i_mu,i_phi) + 1.d0
  enddo

  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write (30,*) 'counts = ',histcounts
     close(30)
  endif

  histcounts = histcounts * 1.d0 / (dble(npts) * 2.d0 * pi / (&
       & size(histcounts) * size(histcounts)))

! THE PROBABILITY OF LANDNG IN BIN X,Y MUST BE THE TOTAL COUNTS IN
! THAT BIN DIVIDED BY THE TOTAL COUNTS IN ALL BINS.  FUTHERMORE, THE
! PROBABILITY DENSTY IN A BIN IS THE PROBABILITY (COUNTS/TOTAL) OVER
! THE AREA COVERED BY THE BIN (2*PI/BINS^2). THUS, THE VALUE HERE IS LIKE
! THE "NORMALIZED CONSTANT" IN THIS MANNER IN THAT IT TURNS THE
! HISTOGRAM INTO A SORT OF PROBABILITY DENSITY ARRAY, SCALED SO THAT
! IT MAY BE COMPARED TO PLOTS OF THE ACTUAL FUNCTION

  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write (30,*) '...done.'
     close(30)
  endif

  return
end subroutine shad_minnaert

!--------------------------------------------------------------------
! UNIGRID_GENERATE(COSI,COSE,PHI):
! UNIFORM GRID GENERATOR - INPUTS COSI AND CALCULATES COSE AND PHI
!--------------------------------------------------------------------

subroutine unigrid_generate(cosi,cose,phi)

  use variables
  implicit none

  integer*8 :: i_mu0,i_cell ! 64-BIT INTEGER COUNTERS

  real :: ran1
  real*8 :: cosi,cose,phi ! MU_0,MU,PHI
  real*8 :: phasefunc_scaled ! FUNCTION THAT CALCULATES SHADOW-MINNAERT
  ! PHASE FUNC

  logical :: notfound

  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write(30,*) 'Inside unigrid'
     write(30,*) 'cosi,cose,phi=',cosi,cose,phi
     close(30)
  endif

  cose = 0.d0
  phi = 0.d0
  notfound = .true.

  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write(30,*) 'Before i_mu0 ceiling'
     close(30)
  endif

  i_mu0 = ceiling(cosi*nint)
  if (i_mu0 == 0) i_mu0 = 1
  
  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write(30,*) 'After i_mu0 ceiling, i_mu0=',i_mu0
     close(30)
  endif

  do while(notfound)

     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*) 'Inside do'
        close(30)
     endif

     i_cell = ceiling(ran1(idum) * ncells(i_mu0))
     
     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*) 'i_cell=',i_cell
        close(30)
     endif

     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*) 'gridcoords mu=',&
             & gridcoords(cellarr(i_mu0,i_cell)%mu_count,2)
        write(30,*) 'gridcoords phi=',&
             & gridcoords(cellarr(i_mu0,i_cell)%phi_count,3)
        write(30,*) 'intlentgth 2,3=',intlength(2),intlength(3)
        write(30,*) 'ran1(idum)=',ran1(idum)
        close(30)
     endif

     cose = gridcoords(cellarr(i_mu0,i_cell)%mu_count,2) + ran1(idum) &
          & * intlength(2) 

     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*) 'cose=',cose
        close(30)
     endif

     phi = gridcoords(cellarr(i_mu0,i_cell)%phi_count,3) + ran1(idum) &
          & * intlength(3)

     if (idebug == 1) then
        open (30, file='debug.out', access='append')
        write(30,*) 'cose,phi=',cose,phi
        close(30)
     endif

     if(cellarr(i_mu0,i_cell)%goodness.eq.1) then
        notfound = .false.
     else if (gridcoords(cellarr(i_mu0,i_cell)%gamma_count,4) + &
          & ran1(idum) * intlength(4).le.phasefunc_scaled(cosi,cose,phi)) then
        notfound = .false.
     endif
     iter = iter + 1
  enddo
  
  return
end subroutine unigrid_generate

!--------------------------------------------------------------------
! WHICH_PART_NEXT(P,E,INEXT,TNEXT,NSCATTER):
! DETERMINE WHICH PARTICLE PHOTON HITS NEXT
!--------------------------------------------------------------------

subroutine which_part_next (p,e,inext,tnext,nscatter)
  use variables
  implicit none

  integer :: inext, i ! INEXT IS INDEX OF PARTICLE PHOTON HITS NEXT 
  integer :: ix, iy, iz, ip 
  integer :: iparticle ! INDEX OF PARTICLE
  integer :: nscatter ! COUNTER - NBR OF SCATTERINGS

  real*8 :: tnext ! TIME OF NEXT SCATTERING
  real*8 :: dp(3) ! DISTANCE BETWEEN PHOTON AND PARTICLE IN QUESTION
  real*8 :: ddp ! DOT PRODUCT OF E AND DP - TO MAKE PARTICLE NOT IN
                ! OPPOSITE DIRECTION OF PHOTON TRAJECTORY
  real*8 :: a,b,c ! QUADRATIC VARIABLES, DEFINED BELOW
  real*8 :: temp 
  real*8 :: t ! TIME TO NEXT PARTICLE
  real*8 :: p(3),e(3) ! PHOTON POSTION AND DIRECTION
  real*8 :: deltacrit ! CRITIAL DISTANCE FROM THE PHOTON RAY THAT PARTICLE
  ! IN GRID BOX COULD BE = DIAGONAL OF BOX + RADIUS OF BIGGEST PARTICLE
  real*8 :: rmax ! MAXIMUM PARTICLE RADIUS
  real*8 :: pminusr(3) ! p - rbox
  real*8 :: dpminusr ! DOT PRODUCT OF E AND PMINUSR

  logical :: flag ! TO SEE IF PHOTON IS LEAVING PARTICLE FIELD

  inext = -1
  tnext = 9.9d99

  rmax = maxval( rad_part )
  deltacrit = sqrt( (dxbox/2.0d0)*(dxbox/2.0d0) + (dybox/2.0d0)* & 
       & (dybox/2.0d0) + (dzbox/2.0d0)*(dzbox/2.0d0) ) + rmax

  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write (30,*)'rmax, deltacrit: ',rmax,deltacrit
     write (30,*)'dxbox,dybox,dzbox:',dxbox,dybox,dzbox
     close(30)
  endif

100 continue
  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write (30,*)'photon position (top of which_part_next): ',p(1:3)
     close(30)
  endif

  do ix = 1, nxbox
     do iy = 1, nybox
        do iz = 1, nzbox

! DETERMINE IF PHOTON COMES WITHIN DELTACRIT OF THIS GRID BOX
! SIMILAR TO SALO AND KARJALAINEN 2003, EQN 4

           pminusr(1:3) = p(1:3) - rbox(1:3,ix,iy,iz)  
           dpminusr = dot_product (e,pminusr)

           a = dot_product (e,e)
           b = 2.0d0 * dpminusr
           c = dot_product(pminusr,pminusr) - & 
                & deltacrit*deltacrit
           temp = b * b - 4.0d0 * a * c
           
       if (idebug == 1) then
          open (30, file='debug.out', access='append')
          write (30,*)'rbox:',rbox(1:3,ix,iy,iz)
          write (30,*)'p:',p(1:3)
          write (30,*)'e:',e(1:3)
          write (30,*)'pminusr:',pminusr(1:3)
          close(30)
       endif

! IF PHOTON PASSES THROUGH THIS GRID BOX
     
           if (temp >= 0.0d0) then

              if (idebug == 1) then
                 open (30, file='debug.out', access='append')
                 write (30,*)'photon goes through this box'
                 write (30,*)'box ix,iy,iz: ',ix,iy,iz
                 close(30)
              endif

              do ip = 1, npartbox(ix,iy,iz) 
                 iparticle = whichpartbox(ix,iy,iz,ip)
                 dp(1:3) = p(1:3) - p_part(1:3,iparticle)
                 ddp = dot_product (e,dp)

! IF PHOTON IS MOVING AWAY FROM PARTICLE, IGNORE THIS PARTICLE
! AND GO ON TO NEXT ONE

                 if (ddp > 0) cycle

! DETERMINE IF PHOTON CAN HIT THIS PARTICLE ON ITS CURRENT TRAJECTORY
! TAKEN FROM SALO AND KARJALAINEN 2003, EQN 4

                 a = dot_product (e,e)
                 b = 2.0d0 * ddp
                 c = dot_product(dp,dp) - & 
                      & rad_part(iparticle) * rad_part(iparticle)
                 temp = b * b - 4.0d0 * a * c
                 
! IF PHOTON HITS PARTICLE, DETERMINE WHEN

                 if (temp >= 0) then
                    t = -(b + sqrt(temp))/ (2.0d0 * a)

! IF THIS IS THE EARLIEST COLLISION SO FAR, SAVE TIME AND PARTICLE INDEX
              
                    if (t < tnext) then
                       inext = iparticle
                       tnext = t
                    endif
                 endif

              enddo

           endif
        
        enddo
     enddo
  enddo

  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write(30,*)'which_part, done cycling particles'
     write(30,*)'which_part, inext = ',inext
     write(30,*)'which_part,  tnext = ',tnext
     close(30)
  endif

! IF PHOTON HITS A PARTICLE, RETURN TO RAY TRACING

  if (inext > 0) return

! IF CALCULATING THE PHOTON PATH BETWEEN PARTICLE AND OBSERVER TO LOOK FOR 
! OBSTACLES (CALLED FROM FLUX_TO_OBS), THEN HIT HAS ALREADY OCCURRED.  
! DO NOT MAKE GHOST PHOTON.  RETURN TO FLUX_TO_OBS

! CHECK WHETHER PHOTON HAS LEFT THE RING PLANE
! IF SO, GO TO NEXT PHOTON

  call inside_ring_plane (p,e,flag)
  if (.not.flag) then
     if (nscatter == 0) then
         nphotpass = nphotpass + 1
         if (satswitch == 1) then
             nphotpassSat = nphotpassSat +1
         elseif (satswitch == 0) then    
             nphotpassSolar = nphotpassSolar +1
         endif    
     endif
     return
  endif

! PHOTON HAS CROSSED EDGE OF PARTICLE FIELD IN X OR Y DIRECTION 
! BUT IS STILL IN RING PLANE
! SHOOT PHOTON IN FROM THE OTHER SIDE OF THE PARTICLE FIELD
! AND CHECK WHICH PARTICLE IT WILL HIT NEXT

  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write (30,*) 'making ghost photon'
     close(30)
  endif

  call make_ghost_photon (p,e)
  go to 100

end subroutine which_part_next


!--------------------------------------------------------------------
! EVAL_FUNC(X):
! FUNCTION TO EVALUATE FOR SIMPSON'S RULE
!--------------------------------------------------------------------

function eval_func (x)

  use variables
  implicit none

  real*8 :: x,eval_func,y
  real*8 :: stx2,talphx2 
  

  if ((phasefunc == 3) .or. (phasefunc == 4)) then

! POWER-LAW CALLISTO TAKEN FROM DONES ET AL 1993 EQN 12

! POWER LAW TO FIT EUROPA DATA FROM DONES ET AL 1993
! FIG 4 BY DMO

    
     eval_func = (pi - x)**expo * sin(x)

  else if (phasefunc == 8) then
  
     talphx2 = tan(x/2.d0)
     stx2 = sqrt(talphx2)
    
     eval_func = exp(-1.0d0*steep*stx2)*((sin(x)+(pi-x)*cos(x))*sin(x))
  
  endif

!  if (idebug == 1) then
!     open (30, file='debug.out', access='append')
!     write (30,*) 'phasefunc,expo,x,eval_func = ',phasefunc,expo,x,eval_func
!     close(30)
!  endif

end function eval_func

!--------------------------------------------------------------------
! PHASEFUNC_SCALED(COSI,COSE,PHI):
! SCALED FUNCTION TO EVALUATE SHADOW-MINNAERT PHASE FUNC
!--------------------------------------------------------------------

function phasefunc_scaled(cosi,cose,phi)

  use variables
  implicit none

  real*8 :: cosi,cose,phi ! MU_0,MU,PHI
  real*8 :: phasefunc_un, phasefunc_scaled
  real*8 :: temp1,temp2

  phasefunc_scaled = 0.d0

! AVOIDING DIVISION BY ZERO
  if (cosi.eq.0.d0) then
     temp1 = 1.d-10
  else
     temp1 = cosi
  endif

  temp2 = max(mu_m,cosi)
  phasefunc_scaled = phasefunc_un(temp1,cose,phi) / &
       & phasefunc_un(temp1,temp2,0.d0)

  return
end function phasefunc_scaled

!--------------------------------------------------------------------
! PHASEFUNC_UN(COSI,COSE,PHI):
! SCALED FUNCTION TO EVALUATE SHADOW-MINNAERT PHASE FUNC
!--------------------------------------------------------------------

function phasefunc_un(cosi,cose,phi)

  use variables
  implicit none

  real*8 :: cosi,cose,phi ! MU_0,MU,PHI
  real*8 :: phasefunc_un
  real*8 :: talpha2 ! TAN(ALPHA/2)
  real*8 :: sta2 ! SQRT[TAN(ALPHA/2)]

  phasefunc_un = 0.d0

  talpha2 = tan(alpha/2.d0)
  sta2 = sqrt(talpha2)

!  alpha = acos( sqrt((1.d0-cosi*cosi) * (1.d0-cose*cose)) * cos(phi) &
!       & + cosi*cose )
  phasefunc_un = exp(-steep * sta2) * (cosi*cose)**(k)

  return
end function phasefunc_un

!--------------------------------------------------------------------
! PHASENORM(SLOCAL):
! FUNCTION FOR CALCULATING PHASE FUNCTION NORMALIZATION FACTORS 
!--------------------------------------------------------------------

function phasenorm(slocal)

  use variables
  implicit none

  integer, parameter :: na = 100 ! NBR OF ALPHA BINS TO INTEGRATE OVER
  integer, parameter :: nl = 100 ! NBR OF LAMBDA BINS TO INEGRATE OVER
  integer :: ia, il ! COUNTERS

  real*8 :: a, l ! ALPHA, LAMBDA
  real*8 :: talpha2 ! TAN(ALPHA/2)
  real*8 :: sta2 ! SQRT[TAN(ALPHA/2)]
  real*8 :: amin, amax, da ! LIMITS OF ALPHA INTEGRATION
  real*8 :: lmin, lmax, dl ! LIMITS OF LAMBDA INTEGRATION
  real*8 :: int1, int2 ! INTEGRALS OVER ALPHA AND LAMBDA

  real*8 :: slocal ! LOCAL STEEP
  real*8 :: phasenorm
  
  amin = 0.d0
  amax = pi
  da = (amax - amin) / dble(na)
  lmax = pi / 2.d0
  int1 = 0.d0
  
  do ia = 1, na
     a = amin + (dble(ia) - 0.5) / dble(na) * (amax - amin)
     talpha2 = tan(a/2.d0)
     sta2 = sqrt(talpha2)
     lmin = -pi / 2.d0 + a
     dl = (lmax - lmin) / dble(nl)
     int2 = 0.d0
     
     do il = 1, nl
        l = lmin + (dble(il) - 0.5) / dble(nl) * (lmax - lmin)
        int2 = int2 + (cos(l))**k * (cos(a-l))**k * dl
     end do
     
     int1 = int1 + int2 * exp(-1.d0*sta2*slocal) * sin(a) * da
     
  enddo

  phasenorm = int1

  return
end function phasenorm

!--------------------------------------------------------------------
! RAN1(IDUM):
! RANDOM NBR GENERATOR
!--------------------------------------------------------------------
! ROUTINE FROM "NUMERICAL RECIPES IN FORTRAN", W.H. PRESS ET AL. (1992)
!
! RETURNS A UNIFORM RANDOM DEVIATE BETWEEN 0.0 AND 1.0 (EXCLUSIVE OF
! ENDPOIT VALUES). CALL WITH IDUM, A NEGATIVE INTEGER TO INITIALIZE,
! THEREAFTER DO NOT ALTER IDUM BETWEEN SUCCESSIVE DEVIATES IN A
! SEQUENCE. RNMX SHOULD APPROXIMATE THE LARGEST FLOATING VALUE THAT IS
! LESS THAN 1.
!--------------------------------------------------------------------

function ran1(idum)
  integer idum,ia,im,iq,ir,ntab,ndiv
  real ran1,am,eps,rnmx
  parameter (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836, &
       ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-7,rnmx=1.-eps)

  integer j,k,iv(ntab),iy
  save iv,iy
  data iv/ntab*0/,iy/0/

  if(idum.le.0.or.iy.eq.0)then
     idum=max(-idum,1)
     do 11 j=ntab+8,1,-1
        k=idum/iq
        idum=ia*(idum-k*iq)-ir*k
        if(idum.lt.0)idum=idum+im
        if(j.le.ntab)iv(j)=idum
11   continue
     iy=iv(1)
  endif
  k=idum/iq
  idum=ia*(idum-k*iq)-ir*k
  if(idum.lt.0)idum=idum+im
  j=1+iy/ndiv
  iy=iv(j)
  iv(j)=idum
  ran1=min(am*iy,rnmx)
  return
end function ran1
  
!--------------------------------------------------------------------
! SIMPSONS(A,B):
! SIMPSON'S RULE NUMERICAL INTEGRATOR
!--------------------------------------------------------------------

function simpsons (a,b)

  use variables
  implicit none

  integer :: i

  real*8 :: simpsons ! RETURNS EVALUATION OF INTEGRAL
  real*8 :: f, eval_func ! FUNCTION TO EVALUATE
  real*8 :: a,b ! INTEGRATION LIMITS
  real*8 :: h ! WIDTH OF INTEGRATION BINS
  real*8 :: simp_inte ! SUM FROM INTEGRATION

! SEE KREYSZIG (1999) PG 873 EQN 7

!  simpsons = 0.d0
!  return

  h = (b - a) / dble(nbins)

  if (idebug == 1) then
     open (30, file='debug.out', access='append')
     write (30,*) 'a,b,nbins = ',a,b,nbins
     close(30)
  endif

  simp_inte = 0.d0

  do i=1,nbins+1

     f = eval_func( h * dble(i - 1) )

     if ((i-1) == 0) then
        simp_inte = simp_inte + f

     else if( mod((i-1),2) == 0) then
        simp_inte = simp_inte + 2.0d0 * f

     else if( mod((i-1),2) == 1.and.i < nbins+1) then
        simp_inte = simp_inte + 4.0d0 * f

     else if(i == nbins+1) then
        simp_inte = simp_inte + f

     end if

  enddo

  simpsons = simp_inte * h / 3.0d0

  return
end function simpsons

