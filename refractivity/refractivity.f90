module refractivity

use fadeeva, only : fadeeva_erfc

real(kind=8), parameter :: R_earth=6370.949d0, k=0.1218d0, Rs=315.d0, pi=4.d0*atan(1.d0),sqtpi=sqrt(4.d0*atan(1.d0))

contains

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

function s(h,R,zen)
  real(8), intent(in) :: h, zen
  real(8) :: s
  s = -R*cos(zen) + sqrt(h*h + 2*h*R + cos(zen)**2 * R*R)
  return
end function s

function inclined_approx(beta,hx)
  real(kind=8), intent(in) :: beta, hx
  real(kind=8) :: sb,cb,t1,t2,kR,betarad,sx,inclined_approx

  
  betarad = beta * pi/180.d0 ! Convert to radians
  sx = s(hx,betarad)
  sb = sin(betarad)
  cb = sqrt(1.d0-sb**2) ! cos(beta)

  t1 = asinh(cb/sb)
  t2 = asinh((sx+R*cb)/(R*sb))
  kR = k*R
  
  inclined_approx = Rs*R*sb*( IK1(kR,kR*sb,t1)-IK1(kR,kR*sb,t2) )
  return

end function inclined_approx

function steep_approx(beta,hx)
  real(8), intent(in) :: beta, hx
  real(8) :: sx, betarad, steep_approx
  real(8) :: SB, SB2, SB3, SB5, SB7, SB9
  real(8) :: KR, KR2, KR3, KR4
  real(8) :: HK, HK2, HK3, EHK, EMHK
  real(8) :: t1, t3, t5, t7, t9

  betarad = beta *pi/180.d0
  sx = s(hx,betarad)

  SB = 1.d0/cos(betarad)
  SB2 = SB*SB
  SB3 = SB2*SB
  SB5 = SB3*SB2
  SB7 = SB5*SB2
  SB9 = SB7*SB2
  
  KR = k*R
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
  t5 = 3.d0*(-5.d0*(24.d0+HK*(24.d0+HK*(12.d0+HK*(4.d0+HK)))) - 8.d0*(6.d0+HK*(6.d0+HK*(3.d0+HK)))*KR - 4.d0*(2.d0+HK*(2.d0+HK))*KR2 + 8.d0*EHK*(15.d0+KR*(6.d0+KR)))*SB5
  t7 = 10.d0*(5.d0*(24.d0+HK*(24.d0+HK*(12.d0+HK*(4.d0+HK)))) + 2.d0*(6.d0+HK*(6.d0+HK*(3.d0+HK)))*KR - 12.d0*EHK*(10.d0+KR))*SB7
  t9 = 35.d0*(-24.d0+24.d0*EHK - HK*(24.d0+HK*(12.d0+HK*(4.d0+HK))))*SB9
  steep_approx = Rs*EMHK/(8.d0*k*KR4)*(t1+t3+t5+t7+t9)
  return

end function steep_approx

end module refractivity
