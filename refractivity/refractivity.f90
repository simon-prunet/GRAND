module refractivity

use fadeeva, only : fadeeva_erfc

real(kind=8), parameter :: R_earth=6370.949d0, k=0.1218d0, Rs=315.d0, pi=4.d0*atan(1.d0),sqtpi=sqrt(4.d0*atan(1.d0))

contains


function s(h,R,coszen)
  real(8), intent(in) :: h, R, coszen
  real(8) :: s
  s = -R*coszen + sqrt(h*h + 2*h*R + coszen**2 * R*R)
  return
end function s

function h(s,R,coszen)
  real(8), intent(in) :: s, R, coszen
  real(8) :: h
  h = -R + sqrt(s**2 + R*R + 2*R*s*coszen)
  return
end function h

subroutine compute_misc(xA,xB,xG_earth,coszeng,hA,hB)
  real(8), intent(in), dimension(3) :: xA, xB
  real(8), intent(out), dimension(3) :: xG_earth
  real(8), intent(out) :: coszeng, hA, hB
  real(8), dimension(3) :: xA_earth, xB_earth, xAB
  real(8) :: nxAB, nxAAB, nxA_earth, nxB_earth, l

  ! Compute coordinates at earth center
  xA_earth = xA
  xB_earth = xB
  xA_earth(2) = xA_earth(2) + R_earth
  xB_earth(2) = xB_earth(2) + R_earth
  ! Compute Intersection of AB segment with earth surface: point G
  xAB = xB-xA
  nxAB = norm2(xAB)
  nxAAB = dot_product(xA_earth,xAB)
  nxA_earth = norm2(xA_earth)
  nxB_earth = norm2(xB_earth)

  l = (-nxAAB + sqrt(nxAAB**2 - nxAB**2*(nxA_earth**2-R_earth**2)))/nxAB**2
  xG_earth = xA_earth + l*xAB
  coszeng = dot_product(xG_earth,xAB)/(R_earth*nxAB)
  hA = nxA_earth - R_earth
  hB = nxB_earth - R_earth
  return
end subroutine compute_misc

real(8) function compute_exact_integral(sA,R,coszen)
  
  real(8), intent(in) :: sA, R, coszen
  real(8) :: ds, hh, ss
  integer, parameter :: npoints=300000

  ds = sA/real(npoints,kind=8)
  ss = 0.d0
  compute_exact_integral = 1.0
  do i=1,npoints-1
    ss = ds*i
    hh = h(ss,R,coszen)
    compute_exact_integral = compute_exact_integral + 2.0*exp(-k*hh)
  enddo
  ! Deal with last point separately as there is no factor 2 for that last point
  ss = ss + ds
  hh = h(ss,R,coszen)
  compute_exact_integral = compute_exact_integral + exp(-k*hh)
  ! Normalize by ds, Rs
  compute_exact_integral = compute_exact_integral * (Rs*ds/2.d0)
  return

end function compute_exact_integral

function logspace(xmin,xmax,npoints) result (res)

  real(8), intent(in) :: xmin, xmax
  integer, intent(in) :: npoints
  real(8), allocatable, dimension (:) :: res
  real(8) :: lxmin, dlx

  allocate(res(npoints))
  dlx = (log(xmax)-log(xmin))/real(npoints-1,8)
  res(1)=xmin
  lxmin = log(xmin)
  do i=1,npoints-1
    res(i+1) = exp(lxmin + i*dlx)   
  enddo    
end function logspace

function compute_table(tabcos,tabh, R) result(res)

  real(8), dimension(:), intent(in) :: tabcos, tabh
  real(8), intent(in) :: R
  integer :: ncos,nh
  real(8), allocatable, dimension(:,:) :: res 
  real(8) :: sA

  ncos = size(tabcos)
  nh = size(tabh)
  print*, 'ncos, nh = ',ncos,nh
  
  allocate(res(ncos,nh))

  do j=1,nh
    do i=1,ncos
      sA = s(tabh(j),R,tabcos(i))
      res(i,j) = compute_exact_integral(sA,R,tabcos(i))
    enddo
  enddo

end function compute_table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TEST APPROXIMATIONS (ASSUMES ANTENNA AT GROUND LEVEL) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(kind=8) function IK1(u,z,w)
real(kind=8), intent(in) :: u,z,w
real(kind=8) :: shw2,shw22,chw2,chw22,shw,chw,beta0,beta1,sq2z
  ! Compute trigonometric function of w/2 and of w
  shw2 = sinh(w/2.d0) 
  shw22 = shw2*shw2 ! sh(w/2)**2
  chw22 = 1.d0 + shw22 !ch(w/2)**2
  chw2 = sqrt(chw22) ! ch(w/2)
  chw = chw22 + shw22 ! ch(w)
  shw = 2.d0*chw2*shw2 ! sh(w)

  ! Now compute beta0 and beta1 terms
  beta0 = chw/shw -0.5d0/shw2
  beta1 = 1.d0/shw - (chw*chw)/(shw*shw*shw)+1.d0/(8.d0*shw2*shw2*shw2)-3.d0/(16.d0*shw2)

  sq2z = sqrt(2.d0*z)

  IK1 = sqtpi/sq2z * (1.+3.d0/(8.d0*z))*fadeeva_erfc(sq2z*shw2)*exp(-z+u) + &
       (beta0/z + beta1/z**2)*exp(-z*chw+u)
  
  return
  
end function IK1

function inclined_approx(beta,hx)
  real(kind=8), intent(in) :: beta, hx
  real(kind=8) :: sb,cb,t1,t2,kR,betarad,cosbeta,sx,inclined_approx

  
  betarad = beta * pi/180.d0 ! Convert to radians
  cb = cos(betarad)
  sx = s(hx,R_earth,cb)
  sb = sqrt(1.d0 - cb**2)

  t1 = asinh(cb/sb)
  t2 = asinh((sx+R_earth*cb)/(R_earth*sb))
  kR = k*R_earth
  
  inclined_approx = Rs*R_earth*sb*( IK1(kR,kR*sb,t1)-IK1(kR,kR*sb,t2) )
  return

end function inclined_approx

function steep_approx(beta,hx)
  real(8), intent(in) :: beta, hx
  real(8) :: sx, betarad, cb, steep_approx
  real(8) :: SB, SB2, SB3, SB5, SB7, SB9
  real(8) :: KR, KR2, KR3, KR4
  real(8) :: HK, HK2, HK3, EHK, EMHK
  real(8) :: t1, t3, t5, t7, t9

  betarad = beta *pi/180.d0
  cb = cos(betarad)
  sx = s(hx,R_earth,cb)

  SB = 1.d0/cb ! Secant
  SB2 = SB*SB
  SB3 = SB2*SB
  SB5 = SB3*SB2
  SB7 = SB5*SB2
  SB9 = SB7*SB2
  
  KR = k*R_earth
  KR2 = KR*KR
  KR3 = KR2*KR
  KR4 = KR3*KR
  
  HK = hx * k 
  HK2 = HK*HK
  HK3 = HK2*HK
  EHK = exp(HK)
  EMHK = 1.d0/EHK
  t1 = 8.d0*KR3*(-1.d0-HK-KR + EHK*(1.d0+KR))*SB
  t3 = 4.d0*KR* ( 6.d0+HK3+3.d0*HK2*(1.d0+KR) + 2.d0*KR*(3.d0+KR) - 2.d0*EHK*(3.d0+KR*(3.d0+KR)) + 2.d0*HK*(3.d0+KR*(3.d0+KR)))*SB3
  t5 = 3.d0*(-5.d0*(24.d0+HK*(24.d0+HK*(12.d0+HK*(4.d0+HK)))) - 8.d0*(6.d0+HK*(6.d0+HK*(3.d0+HK)))*KR & 
       - 4.d0*(2.d0+HK*(2.d0+HK))*KR2 + 8.d0*EHK*(15.d0+KR*(6.d0+KR)))*SB5
  t7 = 10.d0*(5.d0*(24.d0+HK*(24.d0+HK*(12.d0+HK*(4.d0+HK)))) + 2.d0*(6.d0+HK*(6.d0+HK*(3.d0+HK)))*KR - 12.d0*EHK*(10.d0+KR))*SB7
  t9 = 35.d0*(-24.d0+24.d0*EHK - HK*(24.d0+HK*(12.d0+HK*(4.d0+HK))))*SB9
  steep_approx = Rs*EMHK/(8.d0*k*KR4)*(t1+t3+t5+t7+t9)
  return

end function steep_approx

end module refractivity
