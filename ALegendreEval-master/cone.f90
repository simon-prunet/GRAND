module cone

use alegendreeval

contains 

subroutine alegendre_array_p_eval(dnu,dmu,t_array,res)
implicit none

double precision, intent(in):: dnu, dmu
double precision, intent(in), dimension(:) :: t_array
double precision, intent(inout), dimension(:) :: res

double precision :: alpha, alphader, vallogp, vallogq, valp, valq
integer :: i,n

n = size(t_array)

do i=1,n
   call alegendre_eval(dnu,dmu,t_array(i),alpha,alphader,vallogp,vallogq,valp,valq)
   res(i)=valp
enddo
return
end subroutine alegendre_array_p_eval
   
subroutine alegendre_array_proots(dnu,dmu,roots)
! Routine computing all roots of a given associated legendre function
double precision, intent(in) :: dnu, dmu
double precision, intent(inout), dimension(:) :: roots

integer :: i, n
double precision :: rt
n = size(roots)

do i=1,n
   call alegendre_proot(dnu,dmu,i,rt)
   roots(i)=rt
enddo
return
end subroutine alegendre_array_proots

subroutine alegendre_array_jacobi(dmu,n,t_array,wht_array)

double precision, intent(in):: dmu
integer, intent(in):: n
double precision, dimension(:), intent(inout):: t_array
double precision, dimension(:), intent(inout):: wht_array

integer :: j,nn
double precision :: t, wht

nn = ceiling(n/2.0)
do j=1,nn
   call alegendre_jacobi(dmu,n,j,t,wht)
   print*,"t=,wht=",t,wht
   t_array(j) = t
   wht_array(j) = wht
enddo
return

end subroutine alegendre_array_jacobi

end module cone
