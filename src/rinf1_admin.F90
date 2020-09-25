!-*- f90 -*-
! Administrative variables for rinf1
module rinf1_admin
!$  use omp_lib
  implicit none
  interface
     subroutine select_lengths(ind_length, nvec, error, nvecmax, num_lengths)
       integer, intent(in) :: ind_length
       integer, intent(in), optional :: nvecmax, num_lengths
       integer, intent(out) :: nvec
       external error
     end subroutine select_lengths
   end interface
!no$  include 'omp.inc'
  !
  integer :: nitmax, num_lengths
  logical :: rinf1_check = .FALSE.
  logical :: rinf1_check_read = .TRUE.
  integer :: rinf1_checkdim
  double precision :: fconst
  logical :: rinf1_debug = .false.
  parameter (rinf1_checkdim = 128)
  !
  ! NUM_LENGTHS the number of vectors for each case
  ! Adjust so that the final case exceeds cache if available
  ! Note: SMPs may accumulate a lot of cache
  ! rinf1_check: perform correctness checks if true
  ! rinf1_check_read: default is reading. See the subroutine
  !                   rinf1_check_init for details
#ifndef RINF1_MAXLEN
#define RINF1_MAXLEN 1000000
#endif
  integer, parameter :: nnmax = RINF1_MAXLEN
  ! NNMAX is the maximum vector length
  ! NITMAX is the number of iterations (repeats) to circumvent
  !   timer granularity
  !   1.0E8*(CLOCK TICK)*(Mflop/s)
  !   e. g.  if the expected performance is 1 Mflop/s  
  !       =   1000 if TICK is 1.0E-5 sec
  !       = 100000 if TICK is 1.0E-3 sec 
  !   NITMAX=10000 to 50000 for testing
  !   NITMAX=100000 to 1000000 for measuring
  !   In this version of rinf1 only used to measure the offset
  !   (call of dummy)
  parameter(nitmax = 10000000)
  parameter(fconst = -1.0067D0)
  !
  ! input and output data structures
  !
  type rinf1_input
     !   icase    : Which loop type to use
     !   istride  : Which stride to make in loop
     !   num_lengths : How many vectors with varying length to run
     !   niter(3) : Number of iterations for 3 vector lengths 
     !              (up to 3 runs are kept for extrapolation)
     !   titer(3) : Times for the given number of iterations
     !   nit      : Number of iterations for present run
     !   nvec     : Vector length for present run
     integer :: icase, istride, num_lengths, nvec
     double precision :: titer(3)
     integer :: nit, niter(3)
     integer :: nitisone
     logical :: realloc = .true.
     logical :: first_touch = .true.
     character(len=1) :: mode = 'L'
  end type rinf1_input
  type rinf1_output
     !   label    : Characterize performed loop and stride
     !   nops     : Number of operations by data type
     !              1: INTEGER*4
     !              2: INTEGER*8
     !              3: REAL*4
     !              4: REAL*8
     !              5: COMPLEX*8
     !              6: COMPLEX*16
     !   nbytes   : Number of Bytes loaded (1) viz. stored (2).
     !   time     : Time required for 1 iteration.
     !   bw       : Memory bandwidth achieved (in MB/s)
     !   bwdev    : Maximum deviation in %
     !   perf     : Performance in (formal) MFlop/s
     !   perfdev  : Maximum deviation in %
     !   testsum  : for testing correctness
     character(70) label
     integer, dimension(6) :: nops
     integer(kind=8) , dimension(2) :: nbytes
     double precision :: time
     double precision :: bw,bwdev,bwmax,bwmin
     double precision :: perf,perfdev,perfmax,perfmin
     double precision :: testsum
  end type rinf1_output
  double precision, allocatable :: allperfs(:), allbws(:)
  !
  ! global timing variables
  !  t0  time overhead for subroutine call
  !  t1, t2 auxiliaries for begin and end of measurement
  !  tconst target constant for each call of doall
  double precision :: t0, t1, t2, tconst
  !
  ! one-dimensional types and data
  !   please take note:
  !   storage requirement is 140 x nnmax
  !   this means you can have nnmax up to around 13 Million before
  !   needing to switch to 64 Bit
  double precision :: s
  integer, parameter :: i_8 = selected_int_kind(10)
#ifndef ALLOCATE
  double precision, dimension(nnmax) :: a, b, c, d, e, f
  double complex, dimension(nnmax) :: ca,cb,cc,cd
  real, dimension(nnmax) :: sa,sb,sc,sd
  integer(kind=i_8), dimension(nnmax) :: ia,ib,ic,id
  integer(kind=i_8) :: kkaux,kkaux1
  common / arrays / b, c, d, a, e, f, ia, ib, ic, id 
  common / timings / t0
  common / b / kkaux,kkaux1
#ifdef FORTRAN90
  logical, dimension(nnmax) :: mask
#endif
#else
  double precision, dimension (:), allocatable :: a, b, c, d, e, f
  double complex, dimension(:), allocatable :: ca,cb,cc,cd
  real, dimension(:), allocatable :: sa,sb,sc,sd
  integer(kind=i_8), dimension(:), allocatable :: ia,ib,ic,id
#ifdef FORTRAN90
  logical, dimension(:), allocatable :: mask
#endif
#endif
  !
  ! check data
  !
  double precision, dimension(:,:), allocatable :: arr_c
  integer, parameter :: num_tests = 37
  !
  ! file I/O
  !
  integer :: rinf1_err, rinf1_out, rinf1_fout, rinf1_checkio, rinf1_aux, len, len2
  character(80) :: rinf1_config, resfilename, rinf1_check_data
  character(80) :: internalio, internalio2
  !save rinf1_err, rinf1_out, rinf1_fout, rinf_config
  !
  ! MPI and OpenMP tracking
  !
  integer :: all_stat
  integer :: mp_err=0, mp_id=0, mp_processes=1
  integer :: omp_id=0, omp_threads=1, omp_chunk=512
contains
  !
  ! io_init initializes the io variables
  !
  subroutine rinf1_io_init()
    implicit none
    rinf1_out = 6
    rinf1_fout = 12
    rinf1_err = 13
    rinf1_aux = 20
    rinf1_checkio = 21
    rinf1_config = 'rinf1.conf'
    rinf1_check_data = 'rinf1check.dat'
    return
  end subroutine rinf1_io_init
  !
  ! rinf1_init performs setup and reads configuration
  !
  subroutine rinf1_init()
    implicit none
#ifdef USE_MPI
    include 'mpif.h'
#endif
    integer :: i
    double precision :: storage
    double precision :: MPI_WTIM
    external MPI_WTIM
    mp_processes = 1
    mp_id = 0
#ifdef USE_MPI
    call mpi_init(mp_err)
    call mpi_comm_rank(mpi_comm_world,mp_id,mp_err)
    call mpi_comm_size(mpi_comm_world,mp_processes,mp_err)
    if (mp_processes>999999)  call rinf1_abort('More than 9999 MPI processes')
#endif
    if (mp_id == 0 .and. rinf1_debug) then
       rinf1_debug=.true.
    else
       rinf1_debug=.false.
    endif
    ! allocate arrays for collecting performance and
    ! bandwidth results
    allocate(allperfs(mp_processes),allbws(mp_processes),stat=all_stat)
    if (all_stat /= 0) then
       call rinf1_abort('rinf1_init: Could not allocate perf arrays')
    end if
    omp_id = 0
    omp_threads = 1
!$  omp_id = omp_get_thread_num()
!$  omp_threads = omp_get_max_threads()

   if (mp_id == 0) then
       write(6,*) 'Processes: ', mp_processes
       write(6,*) 'Threads:   ', omp_threads
   endif

    ! fixme: construct file name from OpenMP and MPI ID
    write(internalio,'(I3.3,"_",I6.6,"_",I11.11)') omp_threads, mp_processes,RINF1_MAXLEN
    len=len_trim(internalio)
    resfilename = 'rinf1_' // internalio(1:len) // '.res'
    if (rinf1_debug) then
       write(6,fmt=*) ' internalio = ',internalio(1:len)
       write(6,fmt=*) ' resfilename = ',resfilename
    end if
    if (mp_id == 0 ) then
       open(rinf1_fout, file=resfilename, form='FORMATTED', status='UNKNOWN')
       open(rinf1_err, file='rinf1.err', form='FORMATTED', status='UNKNOWN')
    endif
    !
    ! read configuration file
    !
    tconst = 0.5d0
    open(rinf1_aux, file=rinf1_config, form='FORMATTED', status='OLD', err=100)
!   skip header
    read(rinf1_aux,*) 
    read(rinf1_aux, fmt=*, end=20, err=100) num_lengths, tconst
!   skip header
    read(rinf1_aux,*) 
    if (num_lengths == 0) then
       write(6,'(''rinf1_init: Initializing correctness checks'')')
       call rinf1_initchecks()
    end if
    !
    ! compute timing overhead
    !
    t1 = MPI_WTIME()
    do i=1,nitmax
       call dummy()
    end do
    t2 = MPI_WTIME()
    t0 = (t2 - t1)/dble(nitmax)
    if (rinf1_debug) then
       write(6,fmt='(''rinf1_init: Overhead is '',1pe12.5, '' seconds.'')') t0
    end if
    !

    storage = 160.0d0*dble(nnmax) / (1024.0d0**2)
    if (mp_id == 0 ) write(rinf1_out, fmt=1) storage
1   Format(' Memory requirement approx.', f10.1, ' Mbyte')

#ifdef ALLOCATE
    if (mp_id == 0 ) write(6,*) 'Allocating ...  '
    allocate(a(nnmax),b(nnmax),c(nnmax),d(nnmax),e(nnmax),f(nnmax),stat=all_stat)
    if (all_stat /= 0) then
       call rinf1_abort('rinf1_init: Allocating REAL*8 Arrays failed')
    end if
    allocate(ca(nnmax),cb(nnmax),cc(nnmax),cd(nnmax),stat=all_stat)
    if (all_stat /= 0) then
       call rinf1_abort('rinf1_init: Allocating DOUBLE COMPLEX Arrays failed')
    end if
    allocate(sa(nnmax),sb(nnmax),sc(nnmax),sd(nnmax),stat=all_stat)
    if (all_stat /= 0) then
       call rinf1_abort('rinf1_init: Allocating REAL*4 Arrays failed')
    end if
    allocate(ia(nnmax),ib(nnmax),ic(nnmax),id(nnmax),stat=all_stat)
    if (all_stat /= 0) then
       call rinf1_abort('rinf1_init: Allocating INTEGER*4 Arrays failed')
    end if
#ifdef FORTRAN90
    allocate(mask(nnmax),stat=all_stat)
    if (all_stat /= 0) then
       call rinf1_abort('rinf1_init: Allocating LOGICAL Array failed')
    end if 
#endif

#endif
    if (mp_id == 0) then
       write(rinf1_out, fmt=*) 'rinf1_init: Starting benchmark. Please, wait...'
    end if
    return
20  continue
    call rinf1_abort('rinf1_init: Could not read num_lengths.')
    !
    ! error exit
    !
100 continue
    call rinf1_abort('Error reading Config File')
  end subroutine rinf1_init
  !
  ! initialization of array data
  ! to be performed before every call
  ! 
  !
  subroutine rinf1_arrays(input)
    implicit none
    type (rinf1_input) input
    integer :: i,i1
    double precision :: dbli
!
!  possibly reallocate
!
#ifdef ALLOCATE
    if (input%realloc) then
!!!!       if (mp_id==0) wrste(6,*) ' Reallocating arrays with size ',input%nvec
       if (allocated(a))  deallocate(a,b,c,d,e,f)
       if (allocated(ca)) deallocate(ca,cb,cc,cd)
       if (allocated(sa)) deallocate(sa,sb,sc,sd)
       if (allocated(ia)) deallocate(ia,ib,ic,id)
       i = input%nvec
       allocate(a(i),b(i),c(i),d(i),e(i),f(i),stat=all_stat)
       if (all_stat /= 0) then
          call rinf1_abort('rinf1_arrays: Allocating REAL*8 Arrays failed')
       end if
       allocate(ca(i),cb(i),cc(i),cd(i),stat=all_stat)
       if (all_stat /= 0) then
          call rinf1_abort('rinf1_arrays: Allocating DOUBLE COMPLEX Arrays failed')
       end if
       allocate(sa(i),sb(i),sc(i),sd(i),stat=all_stat)
       if (all_stat /= 0) then
          call rinf1_abort('rinf1_arrays: Allocating REAL*4 Arrays failed')
       end if
       allocate(ia(i),ib(i),ic(i),id(i),stat=all_stat)
       if (all_stat /= 0) then
          call rinf1_abort('rinf1_arrays: Allocating INTEGER*4 Arrays failed')
       end if
    end if
#endif

!initialize now
    omp_chunk=input%nvec/omp_threads
    if (input%nvec > 65536*omp_threads) omp_chunk=65536
!$omp parallel if (input%first_touch) 
!no$omp do private(dbli, i1) schedule(static,omp_chunk)
!$omp do private(dbli, i1)
    do i=1,input%nvec
       dbli = dble(i)
       a(i) = 5.842D0 * dbli
       b(i) = 1.39675D0 * cos(dbli)**2 + 4.333D-2
       c(i) = 4.5693D0 / dbli
       d(i) = 9.8124D0 * sqrt(dbli)
       e(i) = -4.815 * dbli
       f(i) = fconst * dbli
       sa(i) = real(a(i))
       sb(i) = real(b(i))
       sc(i) = real(c(i))
       sd(i) = real(d(i))
       i1 = b(i) * dbli
       ib(i) = mod(i1,nnmax) + 1
       ib(i) = mod(i1,nnmax) + 1
       i1 = a(i)
       ic(i) = mod(i1,nnmax) + 1
       ia(i) = (mod(ib(i),60_i_8)*nnmax)/60_i_8 + 1_i_8
       ca(i) = cmplx(a(i),b(i))
       cb(i) = cmplx(b(i),c(i))
       cc(i) = cmplx(c(i),a(i))
       cd(i) = -cb(i)*cc(i)
       if (input%realloc) then
          ia(i)=mod(ia(i),0_i_8+input%nvec)+1_i_8          
          ib(i)=mod(ib(i),0_i_8+input%nvec)+1_i_8         
       endif
    end do
!$omp end do
!$omp end parallel 
!  if (mp_id == 0 )write(6,*) 'first touch done'
  end subroutine rinf1_arrays

  subroutine setzero(a,n)
    implicit none
    integer :: n
    double precision, dimension(n) :: a
    integer :: i
    do i=1,n
      a(i) = 0
    end do
  end subroutine setzero
  subroutine sparse_arrays1(inp)
    implicit none
    type (rinf1_input) inp
    integer :: i1,i
    a = 0.0d0
    i1 = 1333
    if (inp%realloc) then
       do i=1,inp%nvec
       i1 = mod(i1 + i/17 + 47,nnmax)+1
       ia(i) = mod(i1,inp%nvec+1)+1
       if (mod(i,2).eq.0) then
          b(i) = 1
          c(i) = -1
       else
          b(i) = -1
          c(i) = 1
       end if
       end do
    else
       do i=1,nnmax
       i1 = mod(i1 + i/17 + 47,nnmax)+1
       ia(i) = mod(i1,inp%nvec+1)+1
       if (mod(i,2).eq.0) then
          b(i) = 1
          c(i) = -1
       else
          b(i) = -1
          c(i) = 1
       end if
       end do
    endif
  end subroutine sparse_arrays1
  subroutine sparse_arrays2(inp)
    implicit none
    type (rinf1_input) inp
    integer :: i,ioff
    A=0.
    B=1.
    C=1.
    do I=1,nnmax
       ioff=i*5.1243
       ioff = mod(ioff,inp%nvec)+1
       ia(I)=ioff
       ioff=ioff+77+i/15
       ioff = mod(ioff,inp%nvec)+1
       IB(I)=ioff
       ioff=ioff+99+i/19
       ioff = mod(ioff,inp%nvec)+1
       IC(I)=ioff
       if (mod(i,1).eq.0) then
          B(I)=1.
       else
          b(I)=-1.
       endif
    enddo
  end subroutine sparse_arrays2
  !
  ! rinf1_finalize completes measurement process 
  !
  subroutine rinf1_finalize()
    close(rinf1_aux)
    if (rinf1_check) close(rinf1_checkio)
    if (mp_id == 0) then
       write(6,fmt=*) 'Benchmark rinf1 completed successfully'
       close(rinf1_fout)
       close(rinf1_err)
    end if
#ifdef USE_MPI
    call mpi_finalize(mp_err)
#endif
    close(rinf1_aux)
    if (rinf1_check) close(rinf1_checkio)
    stop 
  end subroutine rinf1_finalize
  !
  ! rinf1_abort in case of errors
  !
  subroutine rinf1_abort(error)
    implicit none
#ifdef USE_MPI
    include 'mpif.h'
#endif
    character (len=*) error
    close(rinf1_aux)
    if (rinf1_check) close(rinf1_checkio)
    write(rinf1_out,fmt=*) error
    write(rinf1_err,fmt=*) error
#ifdef USE_MPI
    call mpi_abort(MPI_COMM_WORLD,mp_err)
#endif
    stop 'Aborting rinf1.'
  end subroutine rinf1_abort
  subroutine rinf1_initchecks()
    implicit none
    logical ex
    rinf1_check = .TRUE.
    num_lengths = 1
    allocate(arr_c(rinf1_checkdim,num_tests),stat=all_stat)
    if (all_stat /= 0) then
       call rinf1_abort('rinf1_initchecks: Could not allocate checking array')
    end if
    ! if data exist, then read, otherwise write
    inquire(file=rinf1_check_data,exist=ex)
    if (.not. ex) then
       rinf1_check_read=.FALSE.
       open(rinf1_checkio,file=rinf1_check_data,form='FORMATTED',status='NEW')
    else
       open(rinf1_checkio,file=rinf1_check_data,form='FORMATTED',status='OLD')
       read(rinf1_checkio, fmt='(3(1pd25.16))') arr_c
    end if
    return
  end subroutine rinf1_initchecks
  double precision function rinf1_do_checking(icase)
    implicit none
    integer :: i,icase
    double precision dev_mx, v_mx
! Note: not all data are completely covered in all cases
    if (icase == 4) then
       a(1) = s
       a(2:nnmax) = 0
    else if (icase == 6) then
       a = b(1:nnmax)
    else if (11 .le. icase .and. icase .le. 14 .or. icase == 20) then
       a = dble(ia(1:nnmax))
    else if (icase == 8 .or. icase == 16) then
       a = c(1:nnmax)
    else if (icase == 19) then
       a = dble(ca(1:nnmax))
    else if (icase == 22) then
       a = f(1:nnmax)
    else if (icase == 27) then
       a = sa(1:nnmax)
    end if
    if (rinf1_check_read) then
       dev_mx = 0
       v_mx = 0
       do i=1,rinf1_checkdim
          dev_mx = max(dev_mx, abs(a(i) - arr_c(i,icase)))
          v_mx = max(v_mx, abs(arr_c(i,icase)))
       end do
       if (v_mx .ne. 0) then
          rinf1_do_checking = dev_mx/v_mx
       else
          rinf1_do_checking = dev_mx
       end if
    else
       arr_c(1:rinf1_checkdim,icase) = a(1:rinf1_checkdim)
       rinf1_do_checking = -1.0
    end if
    return
  end function rinf1_do_checking
end module rinf1_admin












