program integrate

  use refractivity
  real(8) :: hmin, hmax, betamin,betamax, dh 
  real(8) :: beta, hh, steep, inclined
  integer :: nh, nbeta
  integer(8) :: c1, c2
  real(8) :: cr

!!$  print*,'Enter hmin,hmax,betamin,betamax'
!!$  read(*,*) hmin,hmax,betamin,betamax
!!$  print*,'You entered ',hmin,hmax,betamin,betamax
!!$
!!$  print*,'Enter nh, nbeta'
!!$  read(*,*) nh,nbeta
!!$  print*,'You entered ',nh,nbeta

!!$  dh = (hmax-hmin)/nh
!!$  dbeta = (betamax-betamin)/nbeta
!!$
!!$  do j=0,nh
!!$     h=hmin+j*dh
!!$     do i=0,nbeta
!!$        beta=betamin+i*dbeta
!!$        inclined = inclined_approx(beta,h)
!!$        steep = steep_approx(beta,h)
!!$        print*,'beta,h,inclined,steep = ',beta,h,inclined,steep
!!$     enddo
!!$  enddo

!!$  print*,'Enter h, beta'
!!$  read(*,*) h, beta
  hh=10.d0
  beta=50.d0
  call system_clock(c1,cr)
  do i=1,1000000
     inclined = inclined_approx(beta,hh)
  enddo
  call system_clock(c2,cr)
  print*,'Mean elapsed time for inclined (ns)= ',(c2-c1)/cr/1000000.*1e9
  call system_clock(c1,cr)
  do i=1,1000000
     steep=steep_approx(beta,hh)
  enddo
  call system_clock(c2,cr)
  print*,'Mean elapsed time for steep (ns)= ',(c2-c1)/cr/1000000.*1e9
  print*,'beta,h,inclined,steep = ',beta,hh,inclined,steep

end program integrate
