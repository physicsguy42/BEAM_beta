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

! to compile type,
!
! >gfortran -O3 make_particles_wakes.F90 -o ../bin/wakex
! 
! and to run type,
!
! ../bin/wakex < ../input/part-in.in > output_particles


  implicit none

  integer, parameter :: nmax = 4000000 
  integer :: i,j,icount,ibin,itot ! COUNTERS
  integer :: idum ! DUMMY SEED FOR RANDOM NBR GENERATOR
  integer :: idist ! SWITCH FOR VERTICAL PARTICLE DISTRIBUTION
  integer :: iwake ! SWITCH FOR WAKES OR NO WAKES
  integer :: npart  ! TOTAL NBR OF PARTICLES
  integer :: npart2,npart3,npart4,npart5 ! MULTIPLES OF NPART
  integer :: npartw,npartg ! NBR PARTICLES IN WAKE AND GAP REGIONS
  integer :: nbins ! NBR OF SIZE BINS
  integer :: ns(nmax) ! NBR OF PARTICLES IN EACH SIZE BIN
  integer :: ntotg,ntotw ! CHECK OF TOTAL NBR OF PARTICLES AS SUM OF NS ARRAY

  character (len=100) :: partdat  ! NAME OF PARTICLE FILE, END WITH .PARTDAT
  character (len=100) :: partpath ! DIRECTORY PATH FOR PARTICLE FILE
  character (len=200) :: partfile ! CONCAT PARTPATH AND PARTDAT

  real*8, parameter :: pi = 3.14159
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
  real*8 :: xmax,xmin,ymax,ymin ! DEFINE CENTER GRID CELL
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
  real*8 :: dw1
  real*8 :: sl, wl
  integer :: time_temp
  real*8 :: frac_temp
  real*8 :: Lfactor
  character*15 output2
  real*8 :: x2new, y2new
  real*8 :: x1only,y1only
  integer :: counter2
  integer :: c_overlaps 
  

! INITIALIZATION
  i = 0; j = 0; icount = 0; ibin = 0; idum = 0
  npart = 0; npart2 = 0; npart3 = 0; npart4 = 0; npart5 = 0
  npartw = 0; npartg = 0; nbins = 0; ntotg = 0; ntotw = 0
  dw = 0.0d0; dg = 0.0d0; a = 0.0d0; b = 0.0d0
  hw = 0.0d0; sh = 0.0d0; sw = 0.0d0
  phi = 0.0d0; tauw = 0.0d0; taug = 0.0d0
  smin = 0.0d0; smax = 0.0d0; ds = 0.0d0
  constg = 0.0d0; constw = 0.d0; q = 0.0d0
  l = 0.0d0; h = 0.0d0
  binratio = 0.0d0; slow = 0.0d0; ntop = 0.0d0; s1 = 0.0d0; s2 = 0.0d0
  dx = 0.0d0; dy = 0.0d0; dz = 0.0d0; dist = 0.0d0
  sumrad = 0.0d0; sumrad2 = 0.0d0; sumrad3 = 0.0d0;sumr = 0.0d0
  xmax = 0.0d0; xmin = 0.0d0; ymax = 0.0d0; ymin = 0.0d0
  ns = 0; x = 0.0d0; y = 0.0d0; z = 0.0d0
  x1 = 0.0d0; y1 = 0.0d0; z1 = 0.0d0; x2 = 0.0d0; y2 = 0.0d0; z2 = 0.0d0
  x3 = 0.0d0; y3 = 0.0d0; z3 = 0.0d0; x4 = 0.0d0; y4 = 0.0d0; z4 = 0.0d0
  x5 = 0.0d0; y5 = 0.0d0; z5 = 0.0d0; s = 0.0d0
  r = 0.0d0; r1 = 0.0d0; r2 = 0.0d0; r3 = 0.0d0; r4 = 0.0d0; r5 = 0.0d0
  rnew = 0.0d0; xnew = 0.0d0; ynew = 0.0d0; znew = 0.0d0
  xrot = 0.0d0; yrot = 0.0d0; zrot = 0.0d0
  idist = 1;h_scale=0.d0;iwake=1;Lfactor = 1.0d0

! replaced hardwired directory path with call to intrinsic 
! get current working directory function

  call getcwd(partpath) 
  
  q = 3.0d0
  oneq = 1.0d0 - q
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

     npartg=(smin**(-2)-smax**(-2))/(s1**(-2)-s2**(-2))  ! q=3 hardwired
     npartw= ((tauw/taug)*(smin**(-2)-smax**(-2)))/(sw*(s1**(-2)-s2**(-2)))

     dg = dw*taug/tauw
     
! new normalization

     constw = 3*dw/(4*pi*(smax-smin))
     
     constg = 3*dg/(4*pi*(smax-smin))
     
          
     h = 4*taug*(smax-smin)/(3*dg*log(smax/smin))    

     write(*,*) 'w= ', h/hw
     
     l = sqrt(2*(1+sw)/(h*constg*(s1**(-2)-s2**(-2))*sw))
     
     l=Lfactor*l
     
     npartg = int(Lfactor**2)*npartg
     npartw = int(Lfactor**2)*npartw     
     
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
     
     const = 3*dg/(4*pi*(smax-smin)) 
     
     h = 4*taug*(smax-smin)/(3*dg*log(smax/smin))
     
     l = sqrt(2/(h*const*(s1**(-2)-s2**(-2)))) ! q=3 hardwired
     
     l = Lfactor*l ! make l larger
     
     
     npart = (smin**(-2)-smax**(-2))/(s1**(-2)-s2**(-2)) ! q=3 hardwired
     
     npart = int(Lfactor**2)*npart

     
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

  write(*,*),'Writing to ',trim(partfile)
  open(11,file=trim(partfile),status='replace')
  write(11,*)trim(partfile)
  write(11,*)xmax,ymax
  write(11,*)npart

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
  y(i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * l
  if(idist.eq.0) then
     z(i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * h
  else if (idist.eq.1) then
     z(i) = real(gasdev2(idum),8) * h/2
  end if
  r(i) = s(1)
   
  icount = 1

  if(icount.ge.ns(ibin)) then
     ibin = ibin + 1
     icount = 0
  endif

!  write(*,*)'Laying down first gap particle...'
!  write(11,*)r(1),x(1),y(1),z(1)

  r1(i) = r(i)
  x1(i) = x(i)
  y1(i) = y(i)
  z1(i) = z(i)
  end do
! LAYING DOWN SUBSEQUENT PARTICLES
  c_overlaps = 0
  do i = int(Lfactor**2),npartg
  
20   continue
     x(i) = -l / 2.0d0 + sl * real(ran1(idum),8)
     y(i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * l
     if(idist.eq.0) then
        z(i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * h
     else if (idist.eq.1) then
        z(i) = real(gasdev2(idum),8) * h/2
        

     end if
     r(i) = s(ibin)
     
! CHECKING FOR PARTICLE OVERLAP
       counter = 0
       counter2 = 0
  25   do j = 1, i-1
        dx = x(i) - x(j)
        dy = y(i) - y(j)
        dz = z(i) - z(j)
        dist = sqrt(dx * dx + dy * dy + dz * dz)
        
!        if (dist.le.r(i)+r(j)) then
!           go to 20
!        endif

! Move overlapping particles by sum of radii along vector connecting their centers, if possible

         x2new = x(i) + ((r(i)+r(j))/dist-1)*dx
         y2new = y(i) + ((r(i)+r(j))/dist-1)*dy
         
         y1only = y(i) + ((r(i)+r(j))/dist-1)*(dx+dy)  ! for only incremting y and z
         x1only = x(i) + ((r(i)+r(j))/dist-1)*(dx+dy)  ! for only incremting x and z
        
       
         if (dist.le.r(i)+r(j)) then
           c_overlaps = c_overlaps + 1
           if((x2new .le. xmax .and. x2new .ge. xmin) .and. &   ! Don't allow 
              (y2new .le. ymax .and. y2new .ge. ymin)) then     ! particle to get 
                x(i) = x(i) + ((r(i)+r(j))/dist-1)*dx           ! moved off the field
                y(i) = y(i) + ((r(i)+r(j))/dist-1)*dy
                z(i) = z(i) + ((r(i)+r(j))/dist-1)*dz   
                counter = counter + 1
                if (counter .lt. 10) then
                   go to 25
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
                   goto 20
                endif                
           endif  
         endif

     enddo
     


     icount = icount + 1

     if(icount.ge.ns(ibin)) then
        write(*,*) 'overlap fract: ',&
     &    s(ibin),ns(ibin),c_overlaps,float(c_overlaps)/float(ns(ibin))
        c_overlaps = 0
        ibin = ibin + 1
        icount = 0
     endif

!     write(*,*)'Laying down gap particle ',i
!     write(11,*)r(i),x(i),y(i),z(i)

     r1(i) = r(i)
     x1(i) = x(i)
     y1(i) = y(i)
     z1(i) = z(i)

  enddo

! ------WAKES-------

! NORMALIZATION OF SIZE DIST TO DETERMINE CONST 

70 continue

  
! Added by D.OLSON on Aug 20, 2014
  if( iwake == 0) then
    do i = 1,nbins+1     
     ntop = (l**2)*h*const / oneq * (s2**oneq - s1**oneq) ! number of particles in current bin
!     write(*,*) 'i, ntop = ', i, int(ntop + 0.5)
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
  s1 = s2/binratio
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
  else if (iwake.eq.0) then
     npartg =0
     npartw = npart
     x(npartg+i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * l
  endif
  y(npartg+i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * l
  if(idist.eq.0) then
     z(npartg+i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * h
  else if (idist.eq.1) then
     z(npartg+i) = real(gasdev2(idum),8) * h/2
  end if
  r(npartg+i) = s(1)

  r1(npartg+i) = r(npartg+i)
  x1(npartg+i) = x(npartg+i)
  y1(npartg+i) = y(npartg+i)
  z1(npartg+i) = z(npartg+i)
!  write(12,*) 'i,x/L, y/L =',i, x(i)/l,y(i)/l 
  
  end do

  icount = Lfactor**2

  if(icount.ge.ns(ibin)) then
     ibin = ibin + 1
     icount = 0
  endif

! LAYING DOWN SUBSEQUENT PARTICLES
!  write(*,*) 'laying down particles ****'
  do i = npartg+int(Lfactor**2) + 1,npartg+npartw
30   continue
     if (iwake.eq.1) then
        x(i) = (sl - l/2.0d0) + wl * real(ran1(idum),8)
     else if (iwake.eq.0) then
        x(i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * l
     endif
     y(i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * l
     if(idist.eq.0) then
        z(i) = (0.5d0 - 1.0d0 * real(ran1(idum),8)) * h
     else if (idist.eq.1) then
        z(i) = real(gasdev2(idum),8) * h/2
     end if
     r(i) = s(ibin)
!     write(*,*) 'particle number,radius = ',i,r(i)
!     write(12,*) 'i,x/L, y/L =',i, x(i)/l,y(i)/l 
     
! CHECKING FOR PARTICLE OVERLAP
     counter = 0
 35  do j = 1, i-1
        dx = x(i) - x(j)
        dy = y(i) - y(j)
        dz = z(i) - z(j)
        dist = sqrt(dx * dx + dy * dy + dz * dz)
        
! Move overlapping particles by sum of radii along vector connecting their centers, if possible

         x2new = x(i) + ((r(i)+r(j))/dist-1)*dx
         y2new = y(i) + ((r(i)+r(j))/dist-1)*dy
         
         y1only = y(i) + ((r(i)+r(j))/dist-1)*(dx+dy)  ! for only incremting y and z
         x1only = x(i) + ((r(i)+r(j))/dist-1)*(dx+dy)  ! for only incremting x and z
        
       
         if (dist.le.r(i)+r(j)) then
           c_overlaps = c_overlaps + 1
!           write(*,*) 'ibin,ns(ibin),r(i),r(j): ',ibin,ns(ibin),r(i),r(j)    
           if((x2new .le. xmax .and. x2new .ge. xmin) .and. &   ! Don't allow 
              (y2new .le. ymax .and. y2new .ge. ymin)) then     ! particle to get 
                x(i) = x(i) + ((r(i)+r(j))/dist-1)*dx           ! moved off the field
                y(i) = y(i) + ((r(i)+r(j))/dist-1)*dy
                z(i) = z(i) + ((r(i)+r(j))/dist-1)*dz   
                counter = counter + 1
                if (counter .lt. 10) then
                   go to 35
                else
                   exit
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
                   goto 30
                endif                
           endif  
         endif
     enddo

     icount = icount + 1

     if(icount.ge.ns(ibin)) then
!         write(*,*) 'overlap fract: ',&
!      &    s(ibin),ns(ibin),c_overlaps,float(c_overlaps)/float(ns(ibin))
        c_overlaps = 0     
        ibin = ibin + 1
        icount = 0
     endif

     r1(i) = r(i)
     x1(i) = x(i)
     y1(i) = y(i)
     z1(i) = z(i)

  enddo

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
!     write(12,*) r(i),x(i),y(i),z(i)

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
        write(11,*)rnew(i),xrot(i),yrot(i),zrot(i)
        itot = itot + 1
        if (modulo(i,100) == 0) then
           write(*,*)'laying down particle:',itot
        endif   
!     else
!        write(12,*) rnew(i),xrot(i),xmin,xmax,yrot(i),ymin,ymax,zrot(i)
    endif

  enddo

  close(11)
!  close(12)

end program make_particles_wakes

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
      
      
      
