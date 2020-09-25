!-*- f90 -*-
! This is needed for evaluation of weighted performance
! and should be called by only one MPI process and the
! OpenMP master thread
!
!
!  Usage: (1) initialization - also allocates storage
!         (2) call after measurement for each vector length
!         (3) finalize - calculates suitable average,
!             also deallocates storage
!
!
module rinf1_wp
  use rinf1_admin
  implicit none
  integer :: rinf1_wp_nval, rinf1_wp_min, rinf1_wp_max, rinf1_wp_perf_size
  character(1) :: rinf1_wp_mode
  integer , parameter:: CACHE_CRIT = 10000
#ifndef ALLOCATE
  integer, dimension(nnmax) :: rinf1_wp_veclen
  double precision, dimension(nnmax) :: rinf1_wp_perf
#else
  integer, allocatable :: rinf1_wp_veclen(:)
  double precision, allocatable :: rinf1_wp_perf(:)
#endif
contains
  subroutine rinf1_wp_init(mode,nmax)
    character(1), intent(in) :: mode
    integer, intent(in) :: nmax
    integer :: all_stat
    rinf1_wp_mode = mode
    rinf1_wp_perf_size = nmax
    rinf1_wp_nval = 0
#ifdef ALLOCATE
    allocate (rinf1_wp_perf(nmax), rinf1_wp_veclen(nmax), stat=all_stat)
    if (all_stat /= 0) then
       write(rinf1_out,'(''rinf1_wp: Tried to allocate '',i12,'' words. '')') 2*nmax
       call rinf1_abort('rinf1_wp: Could not allocate auxiliaries')
    end if
#endif
    return
  end subroutine rinf1_wp_init
  subroutine rinf1_wp_addval(veclen, mflops)
    integer, intent(in) :: veclen
    double precision, intent(in) :: mflops
    rinf1_wp_nval = rinf1_wp_nval + 1
    if (rinf1_wp_nval > rinf1_wp_perf_size) then
       call rinf1_abort('rinf1_wp: Gone beyond size of auxiliaries')
    else
       rinf1_wp_veclen(rinf1_wp_nval) = veclen
       rinf1_wp_perf(rinf1_wp_nval) = mflops
    end if
    return
  end subroutine rinf1_wp_addval
  double precision function rinf1_wp_finalize()
    implicit none
    integer :: i, imin, imax
    integer :: vmin
    double precision :: result, w, d
    result = 0
    select case(rinf1_wp_mode)
    case('c')
       imax = rinf1_wp_nval
       imin = rinf1_wp_nval
       vmin = CACHE_CRIT                                                        
       do i=rinf1_wp_nval-1, 1, -1
          imin = i
          if (rinf1_wp_veclen(i) .le. vmin) exit
       end do
       if (mp_id == 0) then
          Write(rinf1_out,1)  'Linear weight starting with veclen=', rinf1_wp_veclen(imin)
          Write(rinf1_fout,1) 'Linear weight starting with veclen=', rinf1_wp_veclen(imin)
       endif
1      Format(/,A,i5)
       if (imin < imax) then
          d = 5.0d-1/( dble(rinf1_wp_veclen(imax)*1.) - dble(rinf1_wp_veclen(imin)*1.) )
          result = rinf1_wp_perf(imin)*( dble(rinf1_wp_veclen(imin+1)*1.) - dble(rinf1_wp_veclen(imin)*1.) ) + &
               &   rinf1_wp_perf(imax)*( dble(rinf1_wp_veclen(imax)*1.)   - dble(rinf1_wp_veclen(imax-1)*1.))
          do i=imin+1,imax-1
             result = result + &
               &   rinf1_wp_perf(i)   *( dble(rinf1_wp_veclen(i+1)*1.)    - dble(rinf1_wp_veclen(i-1)*1.) )
          end do
       else
	    d=1
          result=rinf1_wp_perf(imax)
       end if
       rinf1_wp_finalize = result * d
    case('L')
       imax = rinf1_wp_nval
       imin = rinf1_wp_nval
       vmin = CACHE_CRIT                    
       do i=rinf1_wp_nval-1, 1, -1
          imin = i
          if (rinf1_wp_veclen(i) .le. vmin) exit
       end do
       if (mp_id == 0) then
          write(rinf1_out,1)  'Logarithmic weight starting with veclen=', rinf1_wp_veclen(imin)
          write(rinf1_fout,1) 'Logarithmic weight starting with veclen=', rinf1_wp_veclen(imin)
       endif
       if (imin < imax) then
          d = 5.0d-1/( ALOG(rinf1_wp_veclen(imax)*1.) - Alog(rinf1_wp_veclen(imin)*1.) )
          result = rinf1_wp_perf(imin)*( alog(rinf1_wp_veclen(imin+1)*1.) - Alog(rinf1_wp_veclen(imin)*1.) ) + &
               &   rinf1_wp_perf(imax)*( alog(rinf1_wp_veclen(imax)*1.)   - alog(rinf1_wp_veclen(imax-1)*1.))
          do i=imin+1,imax-1
             result = result + &
               &   rinf1_wp_perf(i)   *( alog(rinf1_wp_veclen(i+1)*1.)    - alog(rinf1_wp_veclen(i-1)*1.) )
          end do
       else
	    d=1
          result=rinf1_wp_perf(imax)
       end if
       rinf1_wp_finalize = result * d
    case('l')
       imax = rinf1_wp_nval
       imin = 1
       if (mp_id == 0) then
          write(rinf1_out,1) 'Logarithmic weight starting with veclen=', rinf1_wp_veclen(imin)
          write(rinf1_fout,1) 'Logarithmic weight starting with veclen=', rinf1_wp_veclen(imin)
       endif
       if (imin < imax) then
          d = 5.0d-1/( ALOG(rinf1_wp_veclen(imax)*1.) - Alog(rinf1_wp_veclen(imin)*1.) )
          result = rinf1_wp_perf(imin)*( alog(rinf1_wp_veclen(imin+1)*1.) - Alog(rinf1_wp_veclen(imin)*1.) ) + &
               &   rinf1_wp_perf(imax)*( alog(rinf1_wp_veclen(imax)*1.)   - alog(rinf1_wp_veclen(imax-1)*1.))
          do i=imin+1,imax-1
             result = result + &
               &   rinf1_wp_perf(i)   *( alog(rinf1_wp_veclen(i+1)*1.)    - alog(rinf1_wp_veclen(i-1)*1.) )
          end do
       else
          d=1
          result=rinf1_wp_perf(imax)
       end if
       rinf1_wp_finalize = result * d
    case('s')
       imax = rinf1_wp_nval
       imin = rinf1_wp_nval
       vmin = CACHE_CRIT 
       do i=rinf1_wp_nval-1, 1, -1
          if (rinf1_wp_veclen(i) .le. vmin) exit
          imin = i
       end do
       if (mp_id == 0) then
          write(rinf1_out,1)  'Simple average starting with veclen=',rinf1_wp_veclen(imin)
          write(rinf1_fout,1)  'Simple average starting with veclen=',rinf1_wp_veclen(imin)
       endif
       if (imin < imax) then
          d = 1.0d0/(imax-imin+1)
          result=0
          do I=imin,imax
             result = result +  rinf1_wp_perf(i)
          end do
       else
          d=1
          result =rinf1_wp_perf(imax)
       end if
       rinf1_wp_finalize = result * d
    case default
       imax = rinf1_wp_nval
       imin = 1
       if (mp_id == 0) then
          Write(rinf1_out,1) 'Linear weight starting with veclen=', rinf1_wp_veclen(imin)
          Write(rinf1_fout,1) 'Linear weight starting with veclen=', rinf1_wp_veclen(imin)
       endif
       d = 5.0d-1/dble(rinf1_wp_veclen(imax) - rinf1_wp_veclen(imin))
       result = rinf1_wp_perf(imin)*dble(rinf1_wp_veclen(imin+1)-rinf1_wp_veclen(imin)) + &
            & rinf1_wp_perf(imax)*dble(rinf1_wp_veclen(imax)-rinf1_wp_veclen(imax-1))
       do i=imin+1,imax-1
          result = result + rinf1_wp_perf(i)*dble(rinf1_wp_veclen(i+1)-rinf1_wp_veclen(i-1))
       end do
       rinf1_wp_finalize = result * d
    end select
    rinf1_wp_nval = 0
#ifdef ALLOCATE
    deallocate (rinf1_wp_perf, rinf1_wp_veclen, stat=all_stat)
    if (all_stat /= 0) then
       call rinf1_abort('rinf1_wp: Could not deallocate auxiliaries')
    end if
#endif
    return
  end function rinf1_wp_finalize
end module rinf1_wp
