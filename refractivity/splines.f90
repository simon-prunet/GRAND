module splines

implicit none
integer, parameter :: kx = 5, ky = 5
real(8), allocatable, DIMENSION(:) :: tx, ty
integer :: nx, ny


contains

function compute_coeffs(x,y,z,s) result(c)

   real(8), dimension(:), intent(in) :: x,y
   real(8), dimension(:,:), intent(in) :: z
   real(8), intent(in) :: s
   real(8), dimension(:), allocatable :: c


   integer :: mx,my,nxest,nyest,lwrk,kwrk
   integer, allocatable, dimension(:) :: iwrk
   real(8), allocatable, dimension(:) :: wrk, packed
   real(8) :: xb,xe,yb,ye,fp
   integer :: ier

   mx = size(x)
   my = size(y)
 

   nxest = mx+kx+1    
   nyest = my+ky+1
   allocate(tx(nxest))
   allocate(ty(nyest))
   allocate(c((nxest-kx-1)*(nyest-ky-1)))

   lwrk = (4+nxest*(my+2*kx+5) + nyest*(2*ky+5) + mx*(kx+1) + my*(ky+1) + max(my,nxest) )*2
   kwrk = (3+mx+my+nxest+nyest)*2
   allocate(wrk(lwrk),iwrk(kwrk))
   xb = minval(x)
   xe = maxval(x)
   yb = minval(y)
   ye = maxval(y)

   packed = pack(z,.true.)
   call regrid(0,mx,x,my,y,packed,xb,xe,yb,ye,kx,ky,s,nxest,nyest,nx,tx,ny,ty,c,fp,wrk,lwrk,iwrk,kwrk,ier)
   deallocate(wrk,iwrk)

end function compute_coeffs

real(8) function eval(c,x,y)
   
   real(8), DIMENSION(:), INTENT(IN) :: c
   real(8) :: x, y
   real(8), DIMENSION(1) :: xx, yy, zz
   real(8), ALLOCATABLE, DIMENSION(:) :: wrk
   integer :: lwrk, ier
   lwrk = kx+ky+3
   allocate(wrk(lwrk))
   xx(1) = x
   yy(1) = y
   call bispeu(tx,nx,ty,ny,c,kx,ky,xx,yy,zz,1,wrk,lwrk,ier)
   eval = zz(1)   
end function eval

function grid_eval(c,x,y) result(res)

   real(8), DIMENSION(:), INTENT(IN) :: c, x, y
   real(8), DIMENSION(:,:), ALLOCATABLE :: res

   real(8), DIMENSION(:), ALLOCATABLE :: wrk
   integer, DIMENSION(:), ALLOCATABLE :: iwrk
   integer :: lwrk, kwrk, numx, numy, ier

   numx = size(x)
   numy = size(y)
   lwrk = numx*(kx+1)+numy*(ky+1)
   kwrk = numx+numy
 
   allocate(res(numx,numy))
   allocate(wrk(lwrk),iwrk(kwrk))
   call bispev(tx,nx,ty,ny,c,kx,ky,x,numx,y,numy,res,wrk,lwrk,iwrk,kwrk,ier)
   deallocate(wrk,iwrk)
end function grid_eval

end module splines

