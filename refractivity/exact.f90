program exact
  use refractivity
  use splines
  implicit none

  real(8), allocatable, dimension(:) :: tabh, tabcos, tabcos2, tabh2
  integer, parameter :: nh=51, ncos=101
  real(8), parameter :: coszenmin=0.001d0, coszenmax=0.999d0, hmin=0.1d0, hmax=30.d0
  real(8), allocatable, dimension(:) :: coeffs
  real(8), allocatable, dimension(:,:) :: res, res2, app2
  integer :: i

  tabcos = logspace(coszenmin,coszenmax,ncos)
  tabh = logspace(hmin,hmax,nh)

  !print*, tabcos
  !print*, tabh

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
  
  !print*,res2
  coeffs = compute_coeffs(tabh,tabcos,res,0.d0)
  !print*,coeffs
  !print*,tx
  !print*,ty
  print*, tabh(1),tabcos(1),eval(coeffs,tabh(1),tabcos(1))
  print*, tabh2(1),tabcos2(1),eval(coeffs,tabh2(1),tabcos2(1))
  !print*, tabcos(1),tabh(1),res(1,1)

  app2 = grid_eval(coeffs,tabh2,tabcos2)
  !print*, app2
  !print*, res
  !print*,stdev(app2-res2)

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
