!-*- f90 -*-
subroutine setvals(a,b,c)
  implicit none
  double precision a,b,c
  call dummy()
  a = 2.41445D0
  b = 6.443D-3
!  c = a / b
  return
end subroutine setvals
subroutine setvals4(a,b,c)
  implicit none
  real a,b,c
  call dummy()
  a = 2.41445E0
  b = 6.443E-3
!  c = a / b
  return
end subroutine setvals4
