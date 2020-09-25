!-*- f90 -*-
!
! implemented by R. Bader following numerical recipes
!
!@A
!
SUBROUTINE polint(xa,ya,n,x,y,dy)
  !-------------------------------------------------------------------------------
  !  Wert einer Polynom-Interpolation mit dem Neville-Algorithmus
  !
  !  Erlaeuterung der globalen Variablen:
  !  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! >  XA  : DOUBLE XA(N). Argumentwerte, fuer die die Funktionswerte
  !          bekannt sind.
  ! >  YA  : DOUBLE YA(N). Funktionswerte.
  ! >  N   : INTEGER. Zahl der Funktionswerte.
  ! >  X   : DOUBLE. Argumentwert, fuer den die Interpolierte auszuwerten ist.
  !    Y > : DOUBLE. Wert des Interpolationspolynoms.
  !    DY> : DOUBLE. Fehlerabschaetzung.
  !-------------------------------------------------------------------------------
  !@Z
  IMPLICIT NONE
  !  .. parameters ..
  !  ..
  INTEGER          NMAX
  PARAMETER        (NMAX=10)
  !  .. globals ..
  !  ..
  INTEGER          n
  DOUBLE PRECISION dy,x,y,xa(n),ya(n)
  !  .. locals ..
  !  ..
  INTEGER          i,m,ns
  DOUBLE PRECISION den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
  !  .. executable statements ..
  !  ..

  if (NMAX.lt.n) STOP 'POLINT: NMAX too small'
  ns=1
  dif=abs(x-xa(1))
  do i=1,n
     dift=abs(x-xa(i))
     if (dift.lt.dif) then
        ns=i
        dif=dift
     endif
     c(i)=ya(i)
     d(i)=ya(i)
  end do
  y=ya(ns)
  ns=ns-1
  do m=1,n-1
     do i=1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        if (den.eq.0.d0) stop 'POLINT: IDENTICAL XA OCCURRED'
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
     end do
     if (2*ns.lt.n-m)then
        dy=c(ns+1)
     else
        dy=d(ns)
        ns=ns-1
     endif
     y=y+dy
  end do
  return
end SUBROUTINE polint


