module splines

implicit none
contains

function compute_coeffs(x,y,z,s,k) result(res)

   real(8), dimension(:), intent(in) :: x,y
   real(8), dimension(:,:), intent(in) :: z
   real(8), intent(in) :: s
   integer, intent(in) :: k
   real(8), dimension(:), allocatable :: res

   integer :: mx,my,kx,ky,nxest,nyest,nx,ny,lwrk,kwrk
   integer, allocatable, dimension(:) :: iwrk
   real(8), allocatable, dimension(:) :: wrk, tx, ty,packed
   real(8) :: xb,xe,yb,ye,fp
   integer :: ier

   mx = size(x)
   my = size(y)
   kx = k
   ky = k

   nxest = mx+kx+1    
   nyest = my+ky+1
   allocate(tx(nxest))
   allocate(ty(nyest))
   allocate(res(nxest*nyest))

   lwrk = (4+nxest*(my+2*kx+5) + nyest*(2*ky+5) + mx*(kx+1) + my*(ky+1) + max(my,nxest) )*2
   kwrk = (3+mx+my+nxest+nyest)*2
   xb = minval(x)
   xe = maxval(x)
   yb = minval(y)
   ye = maxval(y)

   packed = pack(z,.true.)

   call regrid(0,mx,x,my,y,packed,xb,xe,yb,ye,kx,ky,s,nxest,nyest,nx,tx,ny,ty,res,fp,wrk,lwrk,iwrk,kwrk,ier)

end function compute_coeffs

end module splines

