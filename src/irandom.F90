integer FUNCTION irandom ()
  implicit none
  INTEGER  :: x,y,z,w
  common / random / x, y, z, w
  integer :: k
  ! The  KISS (Keep It Simple Stupid) random number generator. Combines:
  ! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
  ! (2) A 3-shift shift-register generator, period 2^32-1,
  ! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
  !  Overall period>2^123;  Default seeds x,y,z,w.

  x = 69069 * x + 1327217885
  y = mieorishft1 (mieorishft2 (mieorishft3 (y, 13), - 17), 5)
  z = 18000 * iand (z, 65535) + ishft (z, - 16)
  w = 30903 * iand (w, 65535) + ishft (w, - 16)
  irandom = x + y + ishft (z, 16) + w
contains
  integer function mieorishft1(k,n)
    integer :: k,n
    mieorishft1 = ieor (k, ishft (k, n) )
  end function mieorishft1
  integer function mieorishft2(k,n)
    integer :: k,n
    mieorishft2 = ieor (k, ishft (k, n) )
  end function mieorishft2
  integer function mieorishft3(k,n)
    integer :: k,n
    mieorishft3 = ieor (k, ishft (k, n) )
  end function mieorishft3
END FUNCTION irandom
integer FUNCTION irandomet (ix,iy,iz,iw)
  implicit none
  INTEGER :: ix, iy, iz, iw
  INTEGER  :: x, y, z, w
  common / random / x, y, z, w
  x = ix
  y = iy
  z = iz
  w = iw
  irandomet = 1
end function irandomet
  blockdata marsaglia
    INTEGER  :: x,y,z,w
    common /random/ x, y, z, w
    data x,y,z,w /123456789, 362436069, 521288629, 916191069/
  end blockdata
