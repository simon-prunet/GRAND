program exact
  use refractivity
  use splines
  implicit none

  real(8) :: beta, hh
  real(8), allocatable, dimension(:) :: tabh, tabcos, tabcos2, tabh2
  integer, parameter :: nh=6, ncos=11
  real(8), parameter :: coszenmin=0.001d0, coszenmax=0.999d0, hmin=0.1d0, hmax=30.d0
  real(8), allocatable, dimension(:) :: lspace, coeffs
  real(8), allocatable, dimension(:,:) :: res, res2, app2
  integer :: i

  tabcos = logspace(coszenmin,coszenmax,ncos)
  tabh = logspace(hmin,hmax,nh)
  allocate (tabcos2(ncos-1))
  do i=1,ncos-1
    tabcos2(i)=(tabcos(i)+tabcos(i+1))/2.0d0
  enddo
  allocate (tabh2(nh-1))
  do i=1,nh-1
    tabh2(i) = (tabh(i)+tabh(i+1))/2.0d0
  enddo

  res = compute_table(tabcos,tabh,R_earth)
  res2 = compute_table(tabcos2,tabh2,R_earth)
  
  print*,res
  print*,res2
  coeffs = compute_coeffs(tabcos,tabh,res,0.d0,3)
  !print*,coeffs
  !app2 = spline_table(coeffs,tabcos2,tabh2)

  !print*,app2-res2

  ! hh=10.0d0
  ! beta=50.0d0
  ! cb = cos(beta*pi/180.d0)

  ! print*,'cb = ',cb
  ! sA = s(hh,R_earth,cb)
  ! res = compute_exact_integral(sA,R_earth,cb)

  ! print*, 'exact integral for h=',hh,' beta=',beta,' is ',res

  ! lspace = logspace(1.d0,1000.d0,4)
  ! print*,'lspace = ',lspace

end program exact
