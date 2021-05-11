program testsmart

  double precision :: r1(3), r2(3)
  double precision :: n1, dr1r2, neff
  integer*8         c1,c2
  double precision  cr

  print*,'r1 coordinates ?'
  read(*,*) r1(1), r1(2), r1(3)
  print*,'r2 coordinates ?'
  read(*,*) r2(1), r2(2), r2(3)

  call system_clock(c1,cr)
  do i=1,1000000
    call smartexponentialneff(r1, r2, n1, dr1r2, neff)
  enddo
  call system_clock(c2,cr)
  print*,'Mean elapsed time for inclined (ns)= ',(c2-c1)/cr/1000000.*1e9

  print*,'neff =',neff

end program testsmart
