! program generates a table of  planar flux albedo (Apf) for a variety of values of mu0 and shadowing parameter S
! This program is only valid for shadowed Minert surfaces

! this version uses simpson's rule to compute the integral

! compile under the gnu compiler collection use the command,
! > gfortran -O3 planarflux_simpsons.f90 -o planarflux_simpsons

! to generate the table type,
! > ./planarflux_simpsons > apf_table.dat


program planarflux_simpsons

implicit none

real*8 :: mu,mu0, phi,S,pi,alpha,calpha,rmu,rmu0,muk,mu0k
integer :: i,j,l,m,n,p,q
integer, parameter :: NMAX = 2000
integer, parameter :: MUMAX = 200
integer, parameter :: SMAX = 20
real*8 :: cphi(MUMAX+1)  ! store values of cos(phi)
real*8 :: apf(NMAX+1,SMAX+1)
real*8 :: mu_min,mu_max
real*8 :: phi_min,phi_max
real*8 :: S1,S2
real*8 :: curlyr
real*8 :: dS,dmu,dmu0,dphi ! increments
real*8 :: intp ! accumulator for integral
real*8 :: talpha2, sqrt_ta
integer :: wphi(MUMAX+1), wmu(MUMAX+1)
character*10 :: label(NMAX+1)
character*3 :: str1
real*8 :: k

! Set some parameters

    pi = 4*atan(1.0d0)  ! tan(pi/4) = 1.0
    
    phi_min=-pi
    phi_max=pi
    mu_min=0.0d0  
    mu_max=1.0d0   
    S1=0.0d0
    S2=2.0d0
    dphi = (phi_max-phi_min)/MUMAX
    dmu =(mu_max-mu_min)/MUMAX
    dmu0 =(mu_max-mu_min)/NMAX
    dS = (S2-S1)/SMAX
    mu=0.0d0
    mu0=0.0d0
    
    
    k=1.15d0
    curlyr = 1.0d0
    
! store cos(phi) values as an array
    do i = 1, MUMAX+1
       phi = phi_min + dphi * (real(i)-1)  
       cphi(i) = cos(phi)
    enddo

! Compute weights for simpson's rule

! Weights for mu integral
    do p = 1,MUMAX+1
       wmu(p) = 4-2*mod(p,2)
    enddo    
    wmu(1) = 1
    wmu(MUMAX+1) = 1
    
! Weights for phi integral     
    do q = 1,MUMAX+1
       wphi(q) = 4-2*mod(q,2)
    enddo    
    wphi(1) = 1
    wphi(MUMAX+1) = 1

! loop over roughness parameter from 0.0 to 2.0
    do m=1,SMAX+1

        S = S1 + dS*(m-1)
        write(str1,'(F3.1)') S
        label(m) = 'Apf(S='//str1//')'
    
! loop over mu0

      do n=1,NMAX+1 
    
        mu0 = mu_min + (n-1)*dmu0
        rmu0 = (1.0d0-mu0**2)**0.5d0
        mu0k = mu0**(k-1)
        intp = 0.0d0
                
        do p=1,MUMAX+1
           
           mu=mu_min+(real(p)-1.0d0)*dmu
!           write(*,*) 'p,mu = ',p,mu
           rmu = (1.0d0-mu**2)**0.5d0
           muk = mu**k
           
           
           do q =1,MUMAX+1
           
              calpha = cphi(q)*rmu0*rmu + mu*mu0
! Due to floating point arithmetic it is possible to get abs(calpha) slightly larger then 1               
              if (calpha .gt. 1.0d0 ) calpha = 1.0d0
              if (calpha .lt. -1.0d0 ) calpha = -1.0d0              
              alpha = acos(calpha)
              talpha2 = tan(alpha/2.0d0)
              sqrt_ta = sqrt(talpha2)
              intp = intp + real(wphi(q)*wmu(p))*(muk*exp(-1.0d0*S*sqrt_ta)*dmu*dphi)/(2.0d0*pi)
              
           enddo
       enddo
              
       apf(n,m)= (1.0d0+k)*mu0k*curlyr*intp/9.0d0 ! Divide each integral by 3 for simpson's rule
                                                  ! weights are actually, 1/3, 2/3, 4/3 
              
      enddo
   enddo
   

   write(*, '(*(a10,1x))') (label(j) ,j=1,SMAX+1)
   
   do i =1, NMAX+1
     
       write(*, '(*(F7.5,4X))') (real(apf(i,j)) ,j=1,SMAX+1)
           
   enddo
        
end program planarflux_simpsons

    
    