!-*- f90 -*-
subroutine select_lengths(ind_length, nvec, error, nvecmax, num_lengths)
  implicit none
  ! > ind_length   : integer. Counting index.
  !   nvec        >: integer. Calculated vector length.
  ! > error        : external subroutine with string argument.
  ! > nvecmax      : integer, optional. Maximal vector length.
  ! > num_lengths  : integer, optional. Number of vector lengths.
  !
  ! If nvecmax and num_lengths are present, an internally stored array
  ! containing the vector lengths for each counting index are computed.
  !
  ! In any case, the present vector length to be used is calculated.
  !
  integer :: ndim
  parameter (ndim=1000)
  !
  integer, intent(in) :: ind_length
  integer, intent(in), optional :: nvecmax, num_lengths
  integer, intent(out) :: nvec
  external error
  ! locals
  real (kind=8) :: s,ss
  integer :: i, ibs, il, lsum, n, nrem, num_aux
  integer, save :: num = 0 
  logical :: firstcall
  character(len=80) errstr
  ! internal administration
  integer, save :: nvec_array(ndim)
  integer :: l

  if (present (num_lengths) .and. present(nvecmax)) then
     num = num_lengths
     if (num_lengths >= ndim) then
        write(errstr,fmt=10) num_lengths
10      format("Please increase NDIM in select_lengths to ",i10)
        call error(errstr)
     end if
     if (num_lengths <= 20 ) then
        write(errstr,fmt=11) num_lengths
11      format("Please increase the number of vector lengths to more than 20. Now: ",i10 )
        call error(errstr)
     endif
     if (nvecmax <= 1000 ) then
        write(errstr,fmt=12) nvecmax
12      format("Please increase the length of nvecmax to more than 1000. Now: ",i10 )
        call error(errstr)
     endif


     s = Dlog10(dble(nvecmax))

     nvec_array(1) = 1
     nvec_array(2) = 2
     nvec_array(3) = 3
     nvec_array(4) = 4
     nvec_array(5) = 5
     nvec_array(6) = 8
     nvec_array(7) = 15
     nvec_array(8) = 31
     do I=9,num_lengths
        ss=((i*1./num_lengths))**0.9
        nvec_array(i) = nint(10.**(s*ss))
        if (nvec_array(i) <= nvec_array(i-1)+33  ) then 
            nvec_array(i) =  nvec_array(i-1)+33
        endif
     enddo

  end if
  if (ind_length > num) then
     call error('Call of select_lengths with ind_length > num_lengths.')
  end if
  nvec = nvec_array(ind_length)
  return
end subroutine
