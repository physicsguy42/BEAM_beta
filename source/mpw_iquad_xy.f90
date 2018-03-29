program make_particles_wakes
!
!  Author: Lindsey S. Chambers
!  This code will generate a field of particles with wakes to be read 
!  in and used by the MC code.  Wake properties courtesy Josh Colwell.
!

!  update by D.OLSON
!     Instead of selecting new random position when an overlap occurs 
!     Code now moves overlapping paticle by sum of radii of overlapping particles.
!   Sept23,2014
!     Fixed binning so that largest binned particle has a radius equal to smax
!   Oct 20, 2014
!     Rewrite of calcultion of parameters. Made number of bins a fixed parameter 
!     and number of paricle is determined just mainly from gemetric parameters(smin, smax, number of bins).
!     Particle field size is now dtermined from phyical parameters with no explicit 
!     dependance on number of paricles.
!   Nov18,2015
!     Allow for additional factor (Lfactor) in the input file to expand ring patch 
!     in the x and y direction also modify npart by this factor squared.
!   Dec 14,2015 
!     changed seed so that a different 4 digit seed is used each time code is run.
!   Jan 6, 2015
!     Fixed overlaping fix so that particles don't get pushed out of 
!     calculation patch
!
!   Aug 15, 2017
!      Fixed error that made particles move far outside the vertical gaussian envelope.
!   
!   Oct 17, 2017 (Update by P. Estrada)
!      New code that can do size distributions with amin down to 0.4 cm and set up the
!      particle target in a reasonable amount of time. New code is called:
!
!      mpw_iquad_xy.f90
!
!      In order to make these small amin cases (which have millions of particles) run
!      in reasonable time, I had to break up the box into "quadrants" (32 for the xy code). 
!       This is done by defining box edges at -xmax/2,xmax/2,-ymax/2,ymax/2. So initially 
!      the code lays down all particles but does not allow particles to overlap these axes 
!      at all. This is the tradeoff in order to be able to successfully do this, but the 
!      target is large and the bulk of particles very small so it shouldn't really be 
!      an issue. Particles are allowed to "touch" at axes but those particles belong to 
!      different quadrants. Once particles are laid down, then the subroutine QUADRANT where 
!      each particles quadrant is assigned in the array iquad(). The subroutine REARRANGE is 
!      called and places all particles in the 2D arrays xtmp,ytmp,ztmp,rtmp by quadrant. 
!      These can then be run through the overlapping particle search routine and overlaps 
!      are removed by moving the particle randomly within its quadrant. The REARRANGE is
!      called again and particles are placed back in their original order according to
!      particle size bins.
! 
!      NOTE: only set up right now for non-wakes case.****************
! 

!
! to compile the original code we type,
!
! >gfortran -O3 make_particles_wakes.F90 -o ../bin/wakex
!
! for the new code with nmax=5000000, type using the -mcmodel=medium (or large) flag
! Note: this is not needed for nmax=4000000. 
! 
! Beaware that the code will not compile with -mcmodel=large on MacOS systems.
! This is probably due to the fact that -mcmodel is not really fully implemented 
! in gfortran.
!
! >gfortran -O3 -mcmodel=medium mpw_iquad_xy.f90 -o ../bin/wakex
!
! and to run type,
!
! ../bin/wakex < ../input/part-in.in > output_particles


  implicit none

! on MacOS maximum nmax = 3888552 with -mcmodel=medium, if we use -mcmodel = large 
! with MacOS gfortran we get a linkage error

  integer, parameter :: nmax = 3500000, nquad = 32 
  integer :: i,j,icount,ibin,itot,ioverlap,iovercount,k ! COUNTERS
  integer :: idum ! DUMMY SEED FOR RANDOM NBR GENERATOR
  integer :: idist ! SWITCH FOR VERTICAL PARTICLE DISTRIBUTION
  integer :: iwake ! SWITCH FOR WAKES OR NO WAKES
  integer :: npart  ! TOTAL NBR OF PARTICLES
  integer :: npart2,npart3,npart4,npart5 ! MULTIPLES OF NPART
  integer :: npartw,npartg ! NBR PARTICLES IN WAKE AND GAP REGIONS
  integer :: nbins ! NBR OF SIZE BINS
  integer :: ns(nmax) ! NBR OF PARTICLES IN EACH SIZE BIN
  integer :: ntotg,ntotw ! CHECK OF TOTAL NBR OF PARTICLES AS SUM OF NS ARRAY
  integer :: indexi(nmax),indexj(nmax),iquad(nmax)

  character (len=100) :: partdat  ! NAME OF PARTICLE FILE, END WITH .PARTDAT
  character (len=100) :: partpath ! DIRECTORY PATH FOR PARTICLE FILE
  character (len=200) :: partfile ! CONCAT PARTPATH AND PARTDAT

  logical, dimension(4) :: xslquad ! set to TRUE depending where x(sl) is
                                    ! Used in overlap routine
  real*8, parameter :: pi = 4.0d0*DATAN(1.0d0)
  real*8 :: dw,dg,tauw,taug ! VOLUME FILLING FACTOR, OPTICAL DEPTH OF 
  ! WAKE AND GAP
  real*8 :: a,b ! TEMP VARIABLES 
  real*8 :: x(nmax),y(nmax),z(nmax) ! PARTICLE POSITION COORDIATES
                                  ! X = RADIAL, Y = TANGENTIAL, Z = VERTICAL
  real*8 :: s(nmax),smax,smin,ds,r(nmax) ! PARTICLE RADIUS: S = ARRAY FOR 
                                       ! PARTICLE RADIUS OF EACH SIZE BIN;
                                       ! R = RADIUS OF EACH PHYSICAL PARTICLE
  real*8 :: q,constg,constw,const,oneq ! SIZE DISTRIBUTION CONSTANTS
  real*8 :: l,h ! SIZE OF CALCULATION SPACE
  real*8 :: binratio ! BIN SIZE RATIO
  real*8 :: slow ! MIN RADIUS FOR LARGEST SIZE BIN
  real*8 :: ntop ! NBR PARTICLES IN LARGEST SIZE BIN (NEED AT LEAST 1)
  real*8 :: s1,s2 ! LIMITS FOR EACH SIZE BIN
  real*8 :: dx,dy,dz,dist ! CALCULATING DISTANCE BETWEEN PARTICLES
  real :: ran1 ! RANDOM NUMBER
  real*8 :: sumrad, sumrad2, sumrad3 ! SUM OF PARTICLE RADII, SUM SQUARED, SUM 
                                   ! CUBED, AS CALUCALTED FROM SIZE BINS
  real*8 :: sumr ! SUM OF PHYSICAL PARTICLE RADII
  real*8 :: hw, sw ! HW IS HEIGH TO WIDTH OF WAKE RATIO, SW IS SEPARATION TO 
                 ! WIDTH OF WAKE RATIO - VALUES FROM JOSH COLWELL
  real*8 :: sh ! SH IS SEPARATION TO HEIGHT OF WAKE RAIO
  real*8 :: phi,sinphi,cosphi ! VIEWING ANGLE W.R.T WAKE ALIGNMENT MEASURE 
  ! IN RING PLANE
  real*8 :: xmax,xmin,ymax,ymin,zmax ! DEFINE CENTER GRID CELL
  real*8 :: r1(nmax),r2(nmax),r3(nmax),r4(nmax),r5(nmax),rnew(nmax*5)
  real*8 :: x1(nmax),x2(nmax),x3(nmax),x4(nmax),x5(nmax)
  real*8 :: y1(nmax),y2(nmax),y3(nmax),y4(nmax),y5(nmax) 
  real*8 :: z1(nmax),z2(nmax),z3(nmax),z4(nmax),z5(nmax) ! X,Y,Z FOR GRID 
  ! SETTING UP FOR ROTATION
  real*8 :: xnew(nmax*5),ynew(nmax*5),znew(nmax*5) ! X,Y,Z OF ALL GRID CELLS
  real*8 :: xrot(nmax*5),yrot(nmax*5),zrot(nmax*5) ! ROTATED X,Y,Z
  real*8 :: h_scale ! SCALE HEIGHT, USED FOR GAUSSIAN PARTICLE DIST
  real :: gasdev, gasdev2 ! GAUSSIAN RANDOM NUMBER
  integer :: counter
  real*8 :: Lsq,l1,lcu,l2,numer1,denom1,numer2,denom2
  real*8 :: dw1, sfac  ! SFAC allows for particles to overlap a grid boundary by
                        ! a factor. SFAC=1 means no overlap
  real*8 :: sl, wl, xsl
  integer :: time_temp
  real*8 :: frac_temp
  real*8 :: Lfactor, lfac ! ENLARGEMENT FACTOR, wake case REDUCTION factor
  character*15 output2
  real*8 :: x2new, y2new
  real*8 :: x1only,y1only
  integer :: counter2
  real*8, dimension(2*nmax/nquad,nquad) :: xtmp,ytmp,ztmp,rtmp
  real*8, dimension(nquad), parameter :: zq = (/ 1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,&
      &  1.0d0,1.0d0,1.0d0,1.0d0,1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,&
      &  -1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0/)
  integer, dimension(2*nmax/nquad,nquad) :: indx
  integer, dimension(2*nmax/nquad) :: jquad

! INITIALIZATION
  i = 0; j = 0; idum = 0
  npart = 0; npart2 = 0; npart3 = 0; npart4 = 0; npart5 = 0
  npartw = 0; npartg = 0; nbins = 0; 
  dw = 0.0d0; dg = 0.0d0; a = 0.0d0; b = 0.0d0
  hw = 0.0d0; sh = 0.0d0; sw = 0.0d0
  phi = 0.0d0; tauw = 0.0d0; taug = 0.0d0
  smin = 0.0d0; smax = 0.0d0; ds = 0.0d0
  constg = 0.0d0; constw = 0.d0; q = 0.0d0
  l = 0.0d0; h = 0.0d0; xsl = 0.0d0;
  binratio = 0.0d0; slow = 0.0d0; ntop = 0.0d0; s1 = 0.0d0; s2 = 0.0d0
  dx = 0.0d0; dy = 0.0d0; dz = 0.0d0; dist = 0.0d0
  sumrad = 0.0d0; sumrad2 = 0.0d0; sumrad3 = 0.0d0;sumr = 0.0d0
  xmax = 0.0d0; xmin = 0.0d0; ymax = 0.0d0; ymin = 0.0d0; zmax = 0.0d0
  x1 = 0.0d0; y1 = 0.0d0; z1 = 0.0d0; x2 = 0.0d0; y2 = 0.0d0; z2 = 0.0d0
  x3 = 0.0d0; y3 = 0.0d0; z3 = 0.0d0; x4 = 0.0d0; y4 = 0.0d0; z4 = 0.0d0
  x5 = 0.0d0; y5 = 0.0d0; z5 = 0.0d0; s = 0.0d0
  r1 = 0.0d0; r2 = 0.0d0; r3 = 0.0d0; r4 = 0.0d0; r5 = 0.0d0
  rnew = 0.0d0; xnew = 0.0d0; ynew = 0.0d0; znew = 0.0d0
  xrot = 0.0d0; yrot = 0.0d0; zrot = 0.0d0; iquad = 0
  idist = 1;h_scale=0.d0;iwake=1;Lfactor = 1.0d0
  ioverlap = 0; iovercount = 0
  xslquad = .FALSE.
  

! replaced hardwired directory path with call to intrinsic 
! get current working directory function

  call getcwd(partpath) 

! Code is now generalized for all q  
  q = 3.00d0
  oneq = 1.0d0 - q
  sfac = 1.0d0
  lfac = 1.0d0
!  idum = -5328

! set random seed 
  call system_clock(time_temp)
  frac_temp = real(time_temp,8)/10000
  idum = -1*int(10000*(frac_temp-real(int(frac_temp),8)))
  
!  output2 = 'testfile.out'
!  open(12,file=output2,status='unknown')
  

! END INITIALIZATION

  write(*,*),'Enter name for particle file ending with .partdat:'
  read*,partdat
  
  write(*,*) 'partdat = ', partdat
  partfile = trim(partpath) // '/' // partdat
  
  write(*,*) 'filename = ', partfile

  write(*,*),'Include wakes? no wakes=0, wakes=1:'
  read*,iwake
  
  write(*,*) 'iwake= ', iwake
  
  write(*,*),'number of bins : '
  read*,nbins
  
  if(iwake.eq.1) then

     write(*,*),'Enter tau of wake:'
     read*,tauw     

     write(*,*),'Enter tau of gap:'
     read*,taug
     write(*,*) 'taug = ',taug

     write(*,*),'Enter volume filling factor (D) in wakes:'
     read*,dw

     write(*,*) 'Dw = ',dw

     write(*,*),'Enter H/W and S/W ratio for wakes:'
     read(*,*) hw, sw
     
     write(*,*) 'hw = ',hw

     sh = sw/hw

     write(*,*),'Enter orientation angle of wakes:'
     read*,phi

     phi = phi * pi / 180.0d0

     write(*,*),'Enter min and max particle radii in cm:'
     read*,smin,smax

     write(*,*) 'Vertical distribution: 0=uniform, 1=Gaussian:'
     read*,idist
     
     write(*,*) 'Enter an expansion factor for the length of ring patch:'
     read*,Lfactor
     
! binratio
     binratio = (smax/smin)**(1.0d0/nbins)
     
! bin boundary values for largest size bin
     
     s2 = smax*sqrt(binratio)
     s1 = s2/binratio
     
! the initial values of the bin boudaries are set so that the gemetric 
! center of the first bin is smax
     
! number of paricles in the gap is fixed by the constraint that 
! the lagest size bin contains 1 particle

     npartg=(smin**(oneq)-smax**(oneq))/(s1**(oneq)-s2**(oneq))  ! q=3 hardwired
     npartw= ((tauw/taug)*(smin**(oneq)-smax**(oneq)))/(sw*(s1**(oneq)-s2**(oneq)))

     dg = dw*taug/tauw
     
! new normalization

     if (q.eq.3) then
      constw = 3.0d0*dw/(4.0d0*pi*(smax-smin))     
      constg = 3.0d0*dg/(4.0d0*pi*(smax-smin))          
      h = 4.0d0*taug*(smax-smin)/(3.0d0*dg*log(smax/smin))    
      l = sqrt(2*(1+sw)/(h*constg*(s1**(oneq)-s2**(oneq))*sw))
     else
      constw = 3.0d0*dw*(4.0d0-q)/(4.0d0*pi*(smax**(4.0d0-q)-smin**(4.0d0-q)))
      constg = 3.0d0*dg*(4.0d0-q)/(4.0d0*pi*(smax**(4.0d0-q)-smin**(4.0d0-q)))
      h = (4.0d0/3.0d0)*(taug/dg)*(3.0d0-q)*(smax**(4.0d0-q)-smin**(4.0d0-q))/&
           &   ((4.0d0-q)*(smax**(3.0d0-q)-smin**(3.0d0-q)))
      l = DSQRT(-oneq*(1.0d0+sw)/(sw*h*constg*(s1**(oneq)-s2**(oneq))))
     endif

     write(*,*) 'w= ', h/hw
          
     l=Lfactor*l/lfac
     
     npartg = int(Lfactor**2)*npartg/int(lfac**2)
     npartw = int(Lfactor**2)*npartw/int(lfac**2)     
     
     npart = npartg + npartw
     
     sl = l*sw/(1+sw)
     wl = l/(1+sw)
     
     write(*,*) 'smin, smax = ',smin,smax
     write(*,*) 'nbins = ', nbins
     write(*,*)  'binratio = ', binratio
     write(*,*) 's1, s2 = ',s1,s2
     write(*,*) 'sw = ',sw          
     write(*,*) 'l = ',l     
     write(*,*) 'taug = ',taug
     write(*,*) 'npartw = ',npartw
     write(*,*) 'npartg = ',npartg
     write(*,*) 'constw = ',constw
     write(*,*) 'constg = ',constg   
     write(*,*) 'l/(s+w) = ', l/(sl + wl)
               
     write(*,*)'l (in cm):', l
     write(*,*)'s,w,h:',sl,wl,h
     write(*,*)'taug, tauw:',taug,tauw
     write(*,*)'Dw:',dw
     write(*,*)'Dg:',dg
     write(*,*)'npartg, npartw, npart:',npartg,npartw,npart     

  else if(iwake.eq.0) then

     write(*,*)'Enter D and tau:'
     read*,dg,taug

     write(*,*),'Enter min and max particle radii in cm:'
     read*,smin,smax

     write(*,*) 'Vertical distribution: 0=uniform, 1=Gaussian:'
     read*,idist

     write(*,*) 'Enter an expansion factor for the length of ring patch:'
     read*,Lfactor
     
     
! bin ratio
     binratio = (smax/smin)**(1.0d0/nbins)
     
! bin boundary values for largest size bin
     
     s2 = smax*sqrt(binratio)
     
     s1 = s2/binratio
     
     if (q.eq.3.0d0) then
      const = 3.0d0*dg/(4.0d0*pi*(smax-smin))      
      h = 4.0d0*taug*(smax-smin)/(3.0d0*dg*log(smax/smin))    
      l = sqrt(2.0d0/(h*const*(s1**(oneq)-s2**(oneq)))) ! q=3 hardwired     
     else
      const = 3.0d0*dg*(4.0d0-q)/(4.0d0*pi*(smax**(4.0d0-q)-smin**(4.0d0-q)))
      h = (4.0d0/3.0d0)*(taug/dg)*(3.0d0-q)*(smax**(4.0d0-q)-smin**(4.0d0-q))/&
           &   ((4.0d0-q)*(smax**(3.0d0-q)-smin**(3.0d0-q)))
      l = DSQRT(-oneq/(h*const*(s1**(oneq) - s2**(oneq))))
     endif
     l = Lfactor*l/lfac ! make l larger
     
     
     npart = (smin**(oneq)-smax**(oneq))/(s1**(oneq)-s2**(oneq))  
     
     npart = int(Lfactor**2)*npart/int(lfac**2)

     
     write(*,*) 'number of bins :', nbins,binratio,smin,smax
     write(*,*) 's1, s2 : ',s1,s2
     write(*,*)'l (in cm):', l
     write(*,*)'h:',h
     write(*,*)'tau:',taug
     write(*,*)'d:',dg
     write(*,*)'npart:',npart     

  endif

  xmax = l / 2.0d0
  xmin = -l / 2.0d0
  ymax = l / 2.0d0
  ymin = -l / 2.0d0
  zmax = h 

  if(iwake.eq.1) then
   xsl = sl - l/2.0d0
   if(xsl.lt.-xmax/2.0d0) xslquad(1) = .TRUE.
   if(xsl.gt.-xmax/2.0d0 .and. xsl.lt.0.0d0) xslquad(2) = .TRUE.
   if(xsl.gt.0.0d0 .and. xsl.lt.xmax/2.0d0) xslquad(3) = .TRUE.
   if(xsl.gt.xmax/2.0d0) xslquad(4) = .TRUE.
   write(*,*) 'xsl = ',xsl
  endif

  write(*,*),'Writing to ',trim(partfile)
  open(11,file=trim(partfile),status='replace')
  write(11,*)trim(partfile)
  write(11,*)xmax,ymax
  write(11,*)npart

 100 continue
  ntotg = 0; ntotw = 0
  x = 0.0d0; y = 0.0d0; z = 0.0d0; r = 0; s = 0
  ns = 0; icount = 0; ibin = 0
  s2 = smax*sqrt(binratio)    
  s1 = s2/binratio

  if (smax.eq.smin) go to 40
  if (iwake.eq.0) go to 70

! ------GAPS-------

! NORMALIZATION OF SIZE DIST TO DETERMINE CONST 
!
!  constg = npartg * oneq / (smax**oneq - smin**oneq)

! DETERMINING NUMBER OF SIZE BINS SUCH THAT BIN SIZE RATIO IS EQUAL ACROSS
! ALL BINS AND AT LEAST ONE PARTICLE POPULATES LARGEST BIN FOR GAPS
! MAY BE MORE THAN ONE OF LARGEST PARTICLE IN LARGEST BIN FOR WAKES

10 continue 

   
! POPULATE SIZE BINS for gap

  do i = 1,nbins+1
     
     ntop = l*sl*h*constg / oneq * (s2**oneq - s1**oneq) ! number of particles in current bin
     ns(i) = int(ntop + 0.5)
     
     ntotg = ntotg + ns(i)
     s(i) = sqrt(s1*s2)  ! geometric "middle" of bin
     s2 = s1
     s1 = s2/binratio
     
  enddo

!  write(*,*)'Actual nbr of particles used in gap:',ntotg
 
! POPULATE PARTICLE FIELD

  ibin = 1

! LAYING DOWN FIRST PARTICLES 
  do i = 1,int(Lfactor**2) 
   x(i) = -l / 2.0d0 + sl * real(ran1(idum),8)
! Checks to make sure x is away from xmax/2
   if(x(i).gt.-xmax/2.0d0-s(ibin)/sfac .and. x(i).le.-xmax/2.0d0) then
    x(i) = -xmax/2.0d0 - s(ibin)/sfac - real(ran1(idum),8)*(xmax/2.0d0 - s(ibin)/sfac)
   elseif(x(i).ge.-xmax/2.0d0 .and. x(i).lt.-xmax/2.0d0-s(ibin)/sfac) then
    x(i) = MAX(-xmax/2.0d0 + s(ibin)/sfac + real(ran1(idum),8)*(xmax/2.0d0 - DABS(xsl) - &
        &      2.0d0*s(ibin)/sfac),-xmax/2.0d0 + s(ibin)/sfac)
   elseif(x(i).gt.xmax/2.0d0-s(ibin)/sfac .and. x(i).le.xmax/2.0d0) then
    x(i) = xmax/2.0d0 - s(ibin)/sfac - real(ran1(idum),8)*(xmax/2.0d0 - 2.0d0*s(ibin)/sfac)
   elseif(x(i).ge.xmax/2.0d0 .and. x(i).lt.xmax/2.0d0+s(ibin)/sfac) then
    x(i) = MAX(xmax/2.0d0 + s(ibin)/sfac + real(ran1(idum),8)*(xsl - xmax/2.0d0 - s(ibin)/sfac),&
        &    xmax/2.0d0 + s(ibin)/sfac)
   endif
! Makes sure x is away from 0
   x(i) = SIGN(1.0d0,x(i))*MAX(DABS(x(i)),s(ibin)/sfac)
   y(i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * l
! Checks to make sure y is away from ymax/2 
   if(DABS(y(i)).ge.ymax/2.0d0 .and. DABS(y(i)).lt.(ymax/2.0d0 + s(ibin)/sfac)) then
    y(i) = SIGN(1.0d0,y(i))*(ymax/2.0d0 + s(ibin)/sfac + real(ran1(idum),8)*&
      &       (ymax/2.0d0 - s(ibin)/sfac))
   elseif(DABS(y(i)).le.ymax/2.0d0 .and. DABS(y(i)).gt.(ymax/2.0d0 - s(ibin)/sfac)) then
    y(i) = SIGN(1.0d0,y(i))*(ymax/2.0d0 - s(ibin)/sfac - real(ran1(idum),8)*&
       &       (ymax/2.0d0 - 2.0d0*s(ibin)/sfac))
   endif
! Makes sure y is away from 0
   y(i) = SIGN(1.0d0,y(i))*MAX(DABS(y(i)),s(ibin)/sfac)

   if(idist.eq.0) then
      z(i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * h
   else if (idist.eq.1) then
      z(i) = real(gasdev2(idum),8) * h/2
   end if
! Makes sure z is away from 0
   z(i) = SIGN(1.0d0,z(i))*MAX(DABS(z(i)),s(ibin)/sfac)
   r(i) = s(ibin)
   
   icount = 1

   if(icount.ge.ns(ibin)) then
      ibin = ibin + 1
      icount = 0
   endif
 
  end do
! LAYING DOWN SUBSEQUENT PARTICLES - added a "+1" here since the first particle
! was done above...

  do i = int(Lfactor**2) + 1,npartg
  
20   continue
     x(i) = -l / 2.0d0 + sl * real(ran1(idum),8)
! Checks to make sure x is away from xmax/2
     if(x(i).gt.-xmax/2.0d0-s(ibin)/sfac .and. x(i).le.-xmax/2.0d0) then
      x(i) = -xmax/2.0d0 - s(ibin)/sfac - real(ran1(idum),8)*(xmax/2.0d0 - s(ibin)/sfac)
     elseif(x(i).ge.-xmax/2.0d0 .and. x(i).lt.-xmax/2.0d0-s(ibin)/sfac) then
      x(i) = MAX(-xmax/2.0d0 + s(ibin)/sfac + real(ran1(idum),8)*(xmax/2.0d0 - DABS(xsl) - &
          &      2.0d0*s(ibin)/sfac),-xmax/2.0d0 + s(ibin)/sfac)
     elseif(x(i).gt.xmax/2.0d0-s(ibin)/sfac .and. x(i).le.xmax/2.0d0) then
      x(i) = xmax/2.0d0 - s(ibin)/sfac - real(ran1(idum),8)*(xmax/2.0d0 - 2.0d0*s(ibin)/sfac)
     elseif(x(i).ge.xmax/2.0d0 .and. x(i).lt.xmax/2.0d0+s(ibin)/sfac) then
      x(i) = MAX(xmax/2.0d0 + s(ibin)/sfac + real(ran1(idum),8)*(xsl - xmax/2.0d0 - s(ibin)/sfac),&
           &    xmax/2.0d0 + s(ibin)/sfac)
     endif
! Makes sure x is away from 0
     x(i) = SIGN(1.0d0,x(i))*MAX(DABS(x(i)),s(ibin)/sfac)
     y(i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * l
! Checks to make sure y is away from ymax/2 
     if(DABS(y(i)).gt.ymax/2.0d0 .and. DABS(y(i)).lt.(ymax/2.0d0 + s(ibin)/sfac)) then
      y(i) = SIGN(1.0d0,y(i))*(ymax/2.0d0 + s(ibin)/sfac + real(ran1(idum),8)*&
         &       (ymax/2.0d0 - s(ibin)/sfac))
     elseif(DABS(y(i)).lt.ymax/2.0d0 .and. DABS(y(i)).gt.(ymax/2.0d0 - s(ibin)/sfac)) then
      y(i) = SIGN(1.0d0,y(i))*(ymax/2.0d0 - s(ibin)/sfac - real(ran1(idum),8)*&
         &       (ymax/2.0d0 - 2.0d0*s(ibin)/sfac))
     endif
! Makes sure y is away from 0
     y(i) = SIGN(1.0d0,y(i))*MAX(DABS(y(i)),s(ibin)/sfac)

     if(idist.eq.0) then
        z(i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * h
     else if (idist.eq.1) then
        z(i) = real(gasdev2(idum),8) * h/2
     end if
     z(i) = SIGN(1.0d0,z(i))*MAX(DABS(z(i)),s(ibin)/sfac)
     r(i) = s(ibin)

     icount = icount + 1

     if(icount.ge.ns(ibin)) then
        ibin = ibin + 1
        icount = 0
     endif
 
  enddo

! ------WAKES-------

! NORMALIZATION OF SIZE DIST TO DETERMINE CONST 

70 continue

  
! Added by D.OLSON on Aug 20, 2014
  if( iwake == 0) then
    do i = 1,nbins+1     
     ntop = (l**2)*h*const / oneq * (s2**oneq - s1**oneq) ! number of particles in current bin
     ns(i) = int(ntop + 0.5)
     ntotg = ntotg + ns(i)
     s(i) = sqrt(s1*s2)  ! geometric "middle" of bin
     
     s2 = s1
     s1 = s2/binratio
    enddo
   
  endif
  

! POPULATE SIZE BINS
! reset to largest size bin
  s2 = smax*sqrt(binratio)
  s1 = s2/binratio! check against first particle, if overlaps, recalculate position
  if (iwake == 1) then
  do i = 1,nbins+1
      
     ns(i) = ns(i)*(tauw/taug)*(1/sw) ! multiply by scale factor to get number of paricles in bin  
     s(i) = sqrt(s1*s2)               ! assume that the number of particles in each bin is scaled 
     s2 = s1                          ! like the total number of particles between gap and wake
     s1 = s2/binratio
  enddo
  endif

!  write(*,*)'Actual nbr of particles used in wake:',ntotw
 
! POPULATE PARTICLE FIELD

  ibin = 1

! LAYING DOWN FIRST PARTICLE
  do i=1,int(Lfactor**2)
   if (iwake.eq.1) then
     x(npartg+i) = (sl - l/2.0d0) + wl * real(ran1(idum),8)
! Checks to make sure x is away from xmax/2 WAKES
     if(x(npartg+i).gt.-xmax/2.0d0-s(ibin)/sfac .and. x(npartg+i).le.-xmax/2.0d0) then
      x(npartg+i) = MIN(-xmax/2.0d0 - s(ibin)/sfac - real(ran1(idum),8)*(DABS(xsl) - xmax/2.0d0 -&
          &   s(ibin)/sfac),-xmax/2.0d0-s(ibin)/sfac)
     elseif(x(npartg+i).ge.-xmax/2.0d0 .and. x(npartg+i).lt.-xmax/2.0d0-s(ibin)/sfac) then
      x(npartg+i) = -xmax/2.0d0 + s(ibin)/sfac + real(ran1(idum),8)*(xmax/2.0d0 - 2.0d0*s(ibin)/sfac)
     elseif(x(npartg+i).gt.xmax/2.0d0-s(ibin)/sfac .and. x(npartg+i).le.xmax/2.0d0) then
      x(npartg+i) = MIN(xmax/2.0d0 - s(ibin)/sfac - real(ran1(idum),8)*(xmax/2.0d0 - xsl - 2.0d0*&
          &   s(ibin)/sfac),xmax/2.0d0-s(ibin)/sfac)
     elseif(x(npartg+i).ge.xmax/2.0d0 .and. x(npartg+i).lt.xmax/2.0d0+s(ibin)/sfac) then
      x(npartg+i) = xmax/2.0d0 + s(ibin)/sfac + real(ran1(idum),8)*(xmax/2.0d0 - s(ibin)/sfac)
     endif
   else if (iwake.eq.0) then
     npartg =0
     npartw = npart
     x(npartg+i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * l
     if(DABS(x(npartg+i)).ge.xmax/2.0d0 .and. DABS(x(npartg+i)).lt.&
        &       (xmax/2.0d0 + s(ibin)/sfac)) then
      x(npartg+i) = SIGN(1.0d0,x(npartg+i))*(xmax/2.0d0 + s(ibin)/sfac + real(ran1(idum),8)*&
        &       (xmax/2.0d0 - s(ibin)/sfac))
     elseif(DABS(x(npartg+i)).le.xmax/2.0d0 .and. DABS(x(npartg+i)).gt.&
        &       (xmax/2.0d0 - s(ibin)/sfac)) then
      x(npartg+i) = SIGN(1.0d0,x(npartg+i))*(xmax/2.0d0 - s(ibin)/sfac - real(ran1(idum),8)*&
        &       (xmax/2.0d0 - 2.0d0*s(ibin)/sfac))
     endif
   endif

! Make sure x is away from 0 
  x(npartg+i) = SIGN(1.0d0,x(npartg+i))*MAX(DABS(x(npartg+i)),s(ibin)/sfac)
  y(npartg+i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * l
  if(DABS(y(npartg+i)).ge.ymax/2.0d0 .and. DABS(y(npartg+i)).lt.&
   &       (ymax/2.0d0 + s(ibin)/sfac)) then
     y(npartg+i) = SIGN(1.0d0,y(npartg+i))*(ymax/2.0d0 + s(ibin)/sfac + real(ran1(idum),8)*&
   &       (ymax/2.0d0 - s(ibin)/sfac))
  elseif(DABS(y(npartg+i)).le.ymax/2.0d0 .and. DABS(y(npartg+i)).gt.&
   &       (ymax/2.0d0 - s(ibin)/sfac)) then
     y(npartg+i) = SIGN(1.0d0,y(npartg+i))*(ymax/2.0d0 - s(ibin)/sfac - real(ran1(idum),8)*&
   &       (ymax/2.0d0 - 2.0d0*s(ibin)/sfac))
  endif
! Make sure y is away from 0 
  y(npartg+i) = SIGN(1.0d0,y(npartg+i))*MAX(DABS(y(npartg+i)),s(ibin)/sfac)
  if(idist.eq.0) then
     z(npartg+i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * h
  else if (idist.eq.1) then
     z(npartg+i) = real(gasdev2(idum),8) * h/2
  end if
! Make sure z is away from 0 
  z(npartg+i) = SIGN(1.0d0,z(npartg+i))*MAX(DABS(z(npartg+i)),s(ibin)/sfac)
  r(npartg+i) = s(ibin)
  
  end do

  icount = Lfactor**2

  if(icount.ge.ns(ibin)) then
     ibin = ibin + 1
     icount = 0
  endif
 
! LAYING DOWN SUBSEQUENT PARTICLES
!
! Particles are just laid down randomly within the cube defined
! by -xmax < x < xmax, -ymax < y < ymax, and -zmax < z < zmax
!
! Note for this version if the particle code, particles are not
! allowed to be the axes.In other words the entire particle must
! be completely within a quadrant.
!
  do i = npartg+int(Lfactor**2) + 1,npartg+npartw
30   continue
     if (iwake.eq.1) then
        x(i) = (sl - l/2.0d0) + wl * real(ran1(idum),8)
! Checks to make sure x is away from xmax/2 WAKES
        if(x(i).gt.-xmax/2.0d0-s(ibin)/sfac .and. x(i).le.-xmax/2.0d0) then
         x(i) = MIN(-xmax/2.0d0 - s(ibin)/sfac - real(ran1(idum),8)*(DABS(xsl) - xmax/2.0d0 -&
             &   s(ibin)/sfac),-xmax/2.0d0-s(ibin)/sfac)
        elseif(x(i).ge.-xmax/2.0d0 .and. x(i).lt.-xmax/2.0d0-s(ibin)/sfac) then
         x(i) = -xmax/2.0d0 + s(ibin)/sfac + real(ran1(idum),8)*(xmax/2.0d0 - 2.0d0*s(ibin)/sfac)
        elseif(x(i).gt.xmax/2.0d0-s(ibin)/sfac .and. x(i).le.xmax/2.0d0) then
         x(i) = MIN(xmax/2.0d0 - s(ibin)/sfac - real(ran1(idum),8)*(xmax/2.0d0 - xsl - 2.0d0*&
             &   s(ibin)/sfac),xmax/2.0d0-s(ibin)/sfac)
        elseif(x(i).ge.xmax/2.0d0 .and. x(i).lt.xmax/2.0d0+s(ibin)/sfac) then
         x(i) = xmax/2.0d0 + s(ibin)/sfac + real(ran1(idum),8)*(xmax/2.0d0 - s(ibin)/sfac)
        endif
     else if (iwake.eq.0) then
        x(i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * l
! Checks to make sure x is away from xmax/2
        if(DABS(x(i)).ge.xmax/2.0d0 .and. DABS(x(i)).lt.(xmax/2.0d0 + s(ibin)/sfac)) then
         x(i) = SIGN(1.0d0,x(i))*(xmax/2.0d0 + s(ibin)/sfac + real(ran1(idum),8)*&
             &       (xmax/2.0d0 - s(ibin)/sfac))
        elseif(DABS(x(i)).le.xmax/2.0d0 .and. DABS(x(i)).gt.(xmax/2.0d0 - s(ibin)/sfac)) then
         x(i) = SIGN(1.0d0,x(i))*(xmax/2.0d0 - s(ibin)/sfac - real(ran1(idum),8)*&
             &       (xmax/2.0d0 - 2.0d0*s(ibin)/sfac))
        endif
     endif
! Makes sure x,y are away from 0
     x(i) = SIGN(1.0d0,x(i))*MAX(DABS(x(i)),s(ibin)/sfac)
     y(i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * l
! Checks to make sure y is away from ymax/2 
     if(DABS(y(i)).ge.ymax/2.0d0 .and. DABS(y(i)).lt.(ymax/2.0d0 + s(ibin)/sfac)) then
        y(i) = SIGN(1.0d0,y(i))*(ymax/2.0d0 + s(ibin)/sfac + real(ran1(idum),8)*&
       &       (ymax/2.0d0 - s(ibin)/sfac))
     elseif(DABS(y(i)).le.ymax/2.0d0 .and. DABS(y(i)).gt.(ymax/2.0d0 - s(ibin)/sfac)) then
        y(i) = SIGN(1.0d0,y(i))*(ymax/2.0d0 - s(ibin)/sfac - real(ran1(idum),8)*&
       &       (ymax/2.0d0 - 2.0d0*s(ibin)/sfac))
     endif
     y(i) = SIGN(1.0d0,y(i))*MAX(DABS(y(i)),s(ibin)/sfac)
     if(idist.eq.0) then
        z(i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * h
     else if (idist.eq.1) then
        z(i) = real(gasdev2(idum),8) * h/2
     end if
! Makes sure z is away from 0
     z(i) = SIGN(1.0d0,z(i))*MAX(DABS(z(i)),s(ibin)/sfac)
     r(i) = s(ibin)
 
     icount = icount + 1
     if(icount.ge.ns(ibin)) then
        ibin = ibin + 1
        icount = 0
     endif
 
  enddo
! Identify which quadrant a particle is in. Returns iquad.
  iquad = 0
  CALL quadrant(xmax,ymax,x,y,z,iquad,npartg,npartw,nmax)
!
! Temporarily rearrange particles into the arrays xtmp,ytmp,
! and ztmp by quadrant for check.
  jquad = 0; indx = 0
  xtmp = 0; ytmp = 0; ztmp = 0; rtmp = 0
  CALL rearrange(x,y,z,r,xtmp,ytmp,ztmp,rtmp,iquad,jquad,indx,npartg,&
     &      npartw,nmax,nquad,1)
!
! Implement overlap routine. For each quadrant, reassign overlapping particles
! until none remain. Constraint remains that the particles must remain in their
! original quadrant.
  counter2 = 0
  do k = 1,nquad
   ioverlap = 0
   do i = 1,jquad(k)
    counter = 0
37  do j = 1,i-1
     dx = xtmp(i,k) - xtmp(j,k)
     dy = ytmp(i,k) - ytmp(j,k)
     dz = ztmp(i,k) - ztmp(j,k)
     dist = dsqrt(dx*dx + dy*dy + dz*dz)
     if(dist.lt.rtmp(i,k)+rtmp(j,k)) then     
      ioverlap = ioverlap + 1
      counter = counter + 1
!     Target is unresolvable. Produce a new target
      if(counter.gt.10000) then
       write(6,*) 'unresolvable target. Generating new target...'
       GOTO 100
      endif
! Move in x
      if(k.eq.1 .or. k.eq.3 .or. k.eq.14 .or. k.eq.16 .or. k.eq.17 .or. k.eq.19&
           .or. k.eq.30 .or. k.eq.32) then
       if(xslquad(4)) then
        if(xtmp(i,k).le.xsl) then
         xtmp(i,k) = MAX(xmax/2.0d0 + rtmp(i,k)/sfac + real(ran1(idum),8)*(xsl - xmax/2.0d0 - &
           &         rtmp(i,k)/sfac),xmax/2.0d0 + rtmp(i,k)/sfac)
        else
         xtmp(i,k) = MAX(xmax/2.0d0 + rtmp(i,k)/sfac,xsl) + real(ran1(idum),8)*(xmax - xsl - &
           &         rtmp(i,k)/sfac) 
        endif 
       else 
        xtmp(i,k) = xmax/2.0d0 + rtmp(i,k)/sfac + real(ran1(idum),8)*(xmax/2.0d0 - rtmp(i,k)/sfac)
       endif
      elseif(k.eq.2 .or. k.eq.4 .or. k.eq.13 .or. k.eq.15 .or. k.eq.18 .or. k.eq.20&
           .or. k.eq.29 .or. k.eq.31) then
       if(xslquad(3)) then
        if(xtmp(i,k).le.xsl) then
         xtmp(i,k) = MIN(xmax/2.0d0 - rtmp(i,k)/sfac,xsl) - real(ran1(idum),8)*(xsl - 2.0d0*rtmp(i,k)/sfac)
        else
         xtmp(i,k) = MIN(xmax/2.0d0 - rtmp(i,k)/sfac - real(ran1(idum),8)*(xmax/2.0d0 - xsl - 2.0d0*&
             &   rtmp(i,k)/sfac),xmax/2.0d0 - rtmp(i,k)/sfac)
        endif
       else
        xtmp(i,k) = xmax/2.0d0 - rtmp(i,k)/sfac - real(ran1(idum),8)*(xmax/2.0d0 - 2.0d0*rtmp(i,k)/sfac)
       endif
      elseif(k.eq.5 .or. k.eq.7 .or. k.eq.10 .or. k.eq.12 .or. k.eq.21 .or. k.eq.23&
           .or. k.eq.26 .or. k.eq.28) then
       if(xslquad(2)) then
        if(xtmp(i,k).le.xsl) then
         xtmp(i,k) = MAX(-xmax/2.0d0 + rtmp(i,k)/sfac + real(ran1(idum),8)*(xmax/2.0d0 - DABS(xsl) - &
            &      2.0d0*rtmp(i,k)/sfac),-xmax/2.0d0 + rtmp(i,k)/sfac)
        else
         xtmp(i,k) = MAX(-xmax/2.0d0 + rtmp(i,k)/sfac,xsl) + real(ran1(idum),8)*(DABS(xsl) - &
            &      2.0d0*rtmp(i,k)/sfac)
        endif
       else
        xtmp(i,k) = -xmax/2.0d0 + rtmp(i,k)/sfac + real(ran1(idum),8)*(xmax/2.0d0 - 2.0d0*rtmp(i,k)/sfac)
       endif
      else
       if(xslquad(1)) then
        if(xtmp(i,k).le.xsl) then
         xtmp(i,k) = MIN(-xmax/2.0d0 - rtmp(i,k)/sfac,xsl) - real(ran1(idum),8)*(DABS(xmax)- DABS(xsl)&
            &     - rtmp(i,k)/sfac)
        else
         xtmp(i,k) = MIN(-xmax/2.0d0 - rtmp(i,k)/sfac - real(ran1(idum),8)*(DABS(xsl) - xmax/2.0d0 -&
             &   rtmp(i,k)/sfac),-xmax/2.0d0 - rtmp(i,k)/sfac)
        endif
       else
        xtmp(i,k) = -xmax/2.0d0 - rtmp(i,k)/sfac - real(ran1(idum),8)*(xmax/2.0d0 - rtmp(i,k)/sfac)
       endif
      endif
! Move in y
      if(k.eq.1 .or. k.eq.2 .or. k.eq.5 .or. k.eq.6 .or. k.eq.17 .or. k.eq.18&
           .or. k.eq.21 .or. k.eq.22) then
       ytmp(i,k) = ymax/2.0d0 + rtmp(i,k)/sfac + real(ran1(idum),8)*(ymax/2.0d0 - rtmp(i,k)/sfac)
      elseif(k.eq.3 .or. k.eq.4 .or. k.eq.7 .or. k.eq.8 .or. k.eq.19 .or. k.eq.20&
           .or. k.eq.23 .or. k.eq.24) then
       ytmp(i,k) = ymax/2.0d0 - rtmp(i,k)/sfac - real(ran1(idum),8)*(ymax/2.0d0 - 2.0d0*rtmp(i,k)/sfac)
      elseif(k.eq.9 .or. k.eq.10 .or. k.eq.13 .or. k.eq.14 .or. k.eq.25 .or. k.eq.26&
           .or. k.eq.29 .or. k.eq.30) then
       ytmp(i,k) = -ymax/2.0d0 + rtmp(i,k)/sfac + real(ran1(idum),8)*(ymax/2.0d0 - 2.0d0*rtmp(i,k)/sfac)
      else
       ytmp(i,k) = -ymax/2.0d0 - rtmp(i,k)/sfac - real(ran1(idum),8)*(ymax/2.0d0 - rtmp(i,k)/sfac)
      endif
!
! Keeping this here but it is probably fine to just move in x and y and preserve the
! vertical structure, especially if it is gaussian.
!
! move in z
      if(idist.eq.0) then
       ztmp(i,k) = zq(k)*DABS((0.5d0 - 1.0d0 * real(ran1(idum),8)) * h)
      else if (idist.eq.1) then
       ztmp(i,k) = zq(k)*DABS(real(gasdev2(idum),8) * h/2)
      end if
      ztmp(i,k) = SIGN(1.0d0,ztmp(i,k))*MAX(DABS(ztmp(i,k)),rtmp(i,k)/sfac)
      GOTO 37
     endif
    end do
 
   end do

   write(6,*) 'overlapping particles for quad',k,' =',ioverlap 
  end do
 
  CALL rearrange(x1,y1,z1,r1,xtmp,ytmp,ztmp,rtmp,iquad,jquad,indx,npartg,&
     &      npartw,nmax,nquad,2)
 
  go to 60

! ---------------------------------------------------------
! --- FOR SPECIAL CASE WHEN SMAX = SMIN                ----
! ---------------------------------------------------------

40 continue

  s = smax
  nbins = 1
  ibin = 1
  
  do i = 1,npart
50   continue
     x(i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * l
     y(i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * l
     if(idist.eq.0) then
        z(i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * h
     else if (idist.eq.1) then
        z(i) = real(gasdev2(idum),8) * h/2
     end if
     r(i) = s(ibin)
     
     counter = 0
55   do j = 1, i-1
        dx = x(i) - x(j)
        dy = y(i) - y(j)
        dz = z(i) - z(j)
        dist = sqrt(dx * dx + dy * dy + dz * dz)
        
!        if (dist.ge.r(i)+r(j)) then
!           go to 50
!        endif

! Move overlapping particles by sum of radii along vector connecting their centers, if possible

         x2new = x(i) + ((r(i)+r(j))/dist-1)*dx
         y2new = y(i) + ((r(i)+r(j))/dist-1)*dy
         
         y1only = y(i) + ((r(i)+r(j))/dist-1)*(dx+dy)  ! for only incremting y and z
         x1only = x(i) + ((r(i)+r(j))/dist-1)*(dx+dy)  ! for only incremting x and z
        
       
         if (dist.le.r(i)+r(j)) then
           if((x2new .le. xmax .and. x2new .ge. xmin) .and. &   ! Don't allow 
              (y2new .le. ymax .and. y2new .ge. ymin)) then     ! particle to get 
                x(i) = x(i) + ((r(i)+r(j))/dist-1)*dx           ! moved off the field
                y(i) = y(i) + ((r(i)+r(j))/dist-1)*dy
                z(i) = z(i) + ((r(i)+r(j))/dist-1)*dz   
                counter = counter + 1
                if (counter .lt. 10) then
                   go to 55
                endif   
           elseif((x2new .gt. xmax .or. x2new .lt. xmin) .and. &
                 (y1only .le. ymax .and. y1only .ge. ymin)) then 
                y(i) = y(i) + ((r(i)+r(j))/dist-1)*(dx+dy)
                z(i) = z(i) + ((r(i)+r(j))/dist-1)*dz
           elseif((y2new .gt. ymax .or. y2new .lt. ymin) .and. &
                 (x1only .le. xmax .and. x1only .ge. xmin)) then
                x(i) = x(i) + ((r(i)+r(j))/dist-1)*(dx+dy)
                z(i) = z(i) + ((r(i)+r(j))/dist-1)*dz
           else
                counter2 = counter2 + 1
                if (counter2 .lt. 10) then
                   goto 50
                endif                
           endif  
         endif

     enddo
     
!     write(12,*) r(i),x(i),y(i),z(i)
 
     r1(i) = r(i)
     x1(i) = x(i)
     y1(i) = y(i)
     z1(i) = z(i)
    
  enddo

!---------------------------------------
! ROTATION - FOR BOTH SIZE DIST AND NON
!---------------------------------------

60 continue

! SETTING UP GRID - MAKING COPIES OF ORIGINAL BOX

  do i = 1,npart
!     write(12,*) SNGL(r(i)),SNGL(x(i)),SNGL(y(i)),SNGL(z(i)),iquad(i)

     x2(i) = x1(i) + l
     y2(i) = y1(i)
     z2(i) = z1(i)
     r2(i) = r1(i)
 

     x3(i) = x1(i)
     y3(i) = y1(i) - l
     z3(i) = z1(i)
     r3(i) = r1(i)
 

     x4(i) = x1(i) - l
     y4(i) = y1(i)
     z4(i) = z1(i)
     r4(i) = r1(i)
 

     x5(i) = x1(i)
     y5(i) = y1(i) + l
     z5(i) = z1(i)
     r5(i) = r1(i)


  enddo

! SAVING ALL X, Y, Z VALUES IN SINGLE ARRAY

  npart2 = 2*npart
  npart3 = 3*npart
  npart4 = 4*npart
  npart5 = 5*npart

  do i = 1,npart
     xnew(i) = x1(i)
     ynew(i) = y1(i)
     znew(i) = z1(i)
     rnew(i) = r1(i)
  enddo

  do i = npart+1,npart2
     xnew(i) = x2(i-npart)
     ynew(i) = y2(i-npart)
     znew(i) = z2(i-npart)
     rnew(i) = r2(i-npart)
  enddo

  do i = npart2+1,npart3
     xnew(i) = x3(i-npart2)
     ynew(i) = y3(i-npart2)
     znew(i) = z3(i-npart2)
     rnew(i) = r3(i-npart2)
  enddo

  do i = npart3+1,npart4
     xnew(i) = x4(i-npart3)
     ynew(i) = y4(i-npart3)
     znew(i) = z4(i-npart3)
     rnew(i) = r4(i-npart3)
  enddo

  do i = npart4+1,npart5
     xnew(i) = x5(i-npart4)
     ynew(i) = y5(i-npart4)
     znew(i) = z5(i-npart4)
     rnew(i) = r5(i-npart4)
  enddo

! ROTATING ALL PARTICLES

  write(*,*)'rotating particles for wake orientation...'

  cosphi = cos(phi)
  sinphi = sin(phi)

  do i = 1,npart5
     xrot(i) = xnew(i) * cosphi - ynew(i) * sinphi
     yrot(i) = xnew(i) * sinphi + ynew(i) * cosphi
     zrot(i) = znew(i)
  enddo

! KEEPING ONLY THOSE PARTICLES IN ORIGINAL BOX DEFINED BY XMAX, XMIN
! YMAX, YMIN

  itot = 0

  do i = 1,npart5

     if (xrot(i).ge.xmin.and.xrot(i).le.xmax.and.&
          & yrot(i).ge.ymin.and.yrot(i).le.ymax ) then
!        write(11,*)rnew(i)/100.0,xrot(i)/100.0,yrot(i)/100.0,zrot(i)/100.0
        write(11,*)rnew(i),xrot(i),yrot(i),zrot(i)
        itot = itot + 1
!        write(*,*)'laying down particle:',itot
!     else
!        write(12,*) rnew(i),xrot(i),xmin,xmax,yrot(i),ymin,ymax,zrot(i)
    endif

  enddo

  close(11)
  close(12)

end program make_particles_wakes
!***************************************************************
! 
!
  subroutine quadrant(xmax,ymax,x,y,z,iquad,npartg,npartw,nmax)
!
! Routine to determine which quadrant (octant, or more...) that
! a particle is in, assigned as 1-32 (1-16 are z > 0, 17-32 are z < 0)
! Returns values in the array iquad.
!***************************************************************
  implicit none
  integer,intent(in) :: npartg,npartw,nmax
  real*8,intent(in) :: xmax,ymax
  real*8,dimension(nmax),intent(in) :: x,y,z
  integer,dimension(nmax),intent(out) :: iquad
  integer :: i

  do i = 1,npartg+npartw
   if(x(i).gt.xmax/2.0d0 .and. y(i).gt.ymax/2.0d0) then
    if(z(i).gt.0.0d0 .and. iquad(i).eq.0) iquad(i) = 1
    if(z(i).lt.0.0D0 .and. iquad(i).eq.0) iquad(i) = 17
   elseif((x(i).gt.0.0d0 .and. x(i).lt.xmax/2.0d0) .and. y(i).gt.ymax/2.0d0) then   
    if(z(i).gt.0.0d0 .and. iquad(i).eq.0) iquad(i) = 2
    if(z(i).lt.0.0D0 .and. iquad(i).eq.0) iquad(i) = 18
   elseif(x(i).gt.xmax/2.0d0 .and. (y(i).gt.0.0d0 .and. y(i).lt.ymax/2.0d0)) then
    if(z(i).gt.0.0d0 .and. iquad(i).eq.0) iquad(i) = 3
    if(z(i).lt.0.0D0 .and. iquad(i).eq.0) iquad(i) = 19
   elseif((x(i).gt.0.0d0 .and. x(i).lt.xmax/2.0d0) .and. (y(i).gt.0.0d0 .and.&
       &   y(i).lt.ymax/2.0d0)) then
    if(z(i).gt.0.0d0 .and. iquad(i).eq.0) iquad(i) = 4
    if(z(i).lt.0.0D0 .and. iquad(i).eq.0) iquad(i) = 20
   elseif((x(i).lt.0.0d0 .and. x(i).gt.-xmax/2.0d0) .and. y(i).gt.ymax/2.0d0) then
    if(z(i).gt.0.0d0 .and. iquad(i).eq.0) iquad(i) = 5
    if(z(i).lt.0.0D0 .and. iquad(i).eq.0) iquad(i) = 21
   elseif(x(i).lt.-xmax/2.0d0 .and. y(i).gt.ymax/2.0d0) then
    if(z(i).gt.0.0d0 .and. iquad(i).eq.0) iquad(i) = 6
    if(z(i).lt.0.0D0 .and. iquad(i).eq.0) iquad(i) = 22
   elseif((x(i).lt.0.0d0 .and. x(i).gt.-xmax/2.0d0) .and. (y(i).gt.0.0d0 .and.&
       &   y(i).lt.ymax/2.0d0)) then
    if(z(i).gt.0.0d0 .and. iquad(i).eq.0) iquad(i) = 7
    if(z(i).lt.0.0D0 .and. iquad(i).eq.0) iquad(i) = 23
   elseif(x(i).lt.-xmax/2.0d0 .and. (y(i).gt.0.0d0 .and. y(i).lt.ymax/2.0d0)) then
    if(z(i).gt.0.0d0 .and. iquad(i).eq.0) iquad(i) = 8
    if(z(i).lt.0.0D0 .and. iquad(i).eq.0) iquad(i) = 24
   elseif(x(i).lt.-xmax/2.0d0 .and. (y(i).lt.0.0d0 .and. y(i).gt.-ymax/2.0d0)) then
    if(z(i).gt.0.0d0 .and. iquad(i).eq.0) iquad(i) = 9
    if(z(i).lt.0.0D0 .and. iquad(i).eq.0) iquad(i) = 25
   elseif((x(i).lt.0.0d0 .and. x(i).gt.-xmax/2.0d0) .and. (y(i).lt.0.0d0 .and.&
       &     y(i).gt.-ymax/2.0d0)) then
    if(z(i).gt.0.0d0 .and. iquad(i).eq.0) iquad(i) = 10
    if(z(i).lt.0.0D0 .and. iquad(i).eq.0) iquad(i) = 26
   elseif(x(i).lt.-xmax/2.0d0 .and. y(i).lt.-ymax/2.0d0) then
    if(z(i).gt.0.0d0 .and. iquad(i).eq.0) iquad(i) = 11
    if(z(i).lt.0.0D0 .and. iquad(i).eq.0) iquad(i) = 27
   elseif((x(i).lt.0.0d0 .and. x(i).gt.-xmax/2.0d0) .and. y(i).lt.-ymax/2.0d0) then
    if(z(i).gt.0.0d0 .and. iquad(i).eq.0) iquad(i) = 12
    if(z(i).lt.0.0D0 .and. iquad(i).eq.0) iquad(i) = 28
   elseif((x(i).gt.0.0d0 .and. x(i).lt.xmax/2.0d0) .and. (y(i).lt.0.0d0 .and.&
       &    y(i).gt.-ymax/2.0d0)) then
    if(z(i).gt.0.0d0 .and. iquad(i).eq.0) iquad(i) = 13
    if(z(i).lt.0.0D0 .and. iquad(i).eq.0) iquad(i) = 29
   elseif(x(i).gt.xmax/2.0d0 .and. (y(i).lt.0.0d0 .and. y(i).gt.-ymax/2.0d0)) then
    if(z(i).gt.0.0d0 .and. iquad(i).eq.0) iquad(i) = 14
    if(z(i).lt.0.0D0 .and. iquad(i).eq.0) iquad(i) = 30
   elseif((x(i).gt.0.0d0 .and. x(i).lt.xmax/2.0d0) .and. y(i).lt.-ymax/2.0d0) then
    if(z(i).gt.0.0d0 .and. iquad(i).eq.0) iquad(i) = 15
    if(z(i).lt.0.0D0 .and. iquad(i).eq.0) iquad(i) = 31
   else
    if(z(i).gt.0.0d0 .and. iquad(i).eq.0) iquad(i) = 16
    if(z(i).lt.0.0D0 .and. iquad(i).eq.0) iquad(i) = 32
   endif
  end do  
  return
 end subroutine quadrant
!***************************************************************
! 
!
  subroutine rearrange(x,y,z,r,xtmp,ytmp,ztmp,rtmp,iquad,jquad,&
     &    indx,npartg,npartw,nmax,nquad,icase)
!
! Rearranges the particles by quadrant for ease of manipulation.
! Returns the 2D arrays xtmp,ytmp,ztmp,rtmp, and jquad.
!***************************************************************
  implicit none
  integer,intent(in) :: icase,npartg,npartw,nmax,nquad
  integer,dimension(nmax),intent(in) :: iquad
  real*8,dimension(nmax),intent(inout) :: x,y,z,r
  integer,dimension(2*nmax/nquad),intent(inout) :: jquad
  integer,dimension(2*nmax/nquad,nquad),intent(inout) :: indx
  real*8,dimension(2*nmax/nquad,nquad),intent(inout) :: xtmp,ytmp,ztmp,rtmp
  integer :: i,j,k

  if(icase.eq.1) then
! Assign particles into their quadrants
   do k = 1,npartg+npartw
    i = iquad(k)
    jquad(i) = jquad(i) + 1
    j = jquad(i)
    xtmp(j,i) = x(k)
    ytmp(j,i) = y(k)
    ztmp(j,i) = z(k)
    rtmp(j,i) = r(k)
    indx(j,i) = k
   end do
  else
! Return particles to their original array positions
   do k = 1,nquad
    do i = 1,jquad(k)
     x(indx(i,k)) = xtmp(i,k)
     y(indx(i,k)) = ytmp(i,k)
     z(indx(i,k)) = ztmp(i,k)
     r(indx(i,k)) = rtmp(i,k)
    end do
   end do
  endif

!  do  j = 1,jquad(1)
!   write(6,*) j,rtmp(j,1)
!  end do
!  stop
 
  return
 end subroutine rearrange 
!***********************************************************
!       Auxillary routines from NUMERICAL RECIPES
!
!**************************************************************
  function gasdev(idum)
!**************************************************************
    integer idum
    real gasdev
    integer iset
    real fac,gset,rsq,v1,v2,ran1
    save iset,gset
    data iset/0/
!
    if (iset == 0) then
 1    v1 = 2.0 * ran1(idum) - 1.0
      v2 = 2.0 * ran1(idum) - 1.0
      rsq = v1**2 + v2**2
      if (rsq >= 1.or.rsq == 0) goto 1
      fac = sqrt(-2.0 * log(rsq)/rsq)
      gset = v1 * fac
      gasdev = v2 * fac
      iset = 1
    else
      gasdev = gset
      iset = 0
    end if
!
    return
    end
!**************************************************************
       FUNCTION RAN1(IDUM)
!**************************************************************
!       creates uniformly distributed pseudorandom numbers
!       rnd(1)

      DIMENSION R(97)
      PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
      PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
      PARAMETER (M3=243000,IA3=4561,IC3=51349)
!      DATA IFF /0/
!      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN

      save
      IF (IDUM.LT.0) THEN
        IFF=1
        IX1=MOD(IC1-IDUM,M1)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IX1,M2)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX3=MOD(IX1,M3)
        DO 11 J=1,97
          IX1=MOD(IA1*IX1+IC1,M1)
          IX2=MOD(IA2*IX2+IC2,M2)
          R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
11      CONTINUE
        IDUM=1
      ENDIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      J=1+(97*IX3)/M3
      IF(J.GT.97.OR.J.LT.1) then
        j=1
        endif
      RAN1=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      RETURN
      END
      
!**************************************************************
  function gasdev2(idum)
!**************************************************************
! exclude random numbers greater then 2.0 

! *****why do we still sometimes get gasdev2 > 2? More testing ****
      integer idum
      
      real    gasdev2, temp
      
      
   11 temp = gasdev(idum)
      
      if(abs(temp) > 2.0) then
          
          go to 11
      else
      
          gasdev2 = temp
     end if     
     
     return
     
     end
      
      
      
