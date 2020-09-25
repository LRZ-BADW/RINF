!-*- f90 -*-
subroutine doall(inp, res)
  !
  ! doall contains all loops in Fortran 77 style. 
  !
  ! > inp    : RINF1_INPUT. 
  !   res   >: RINF1_OUTPUT.
  ! See Module rinf1_admin for explanation of these types.
  !==========================================================================
  use rinf1_admin
  implicit none
#ifdef USE_MPI
  include 'mpif.h'
#endif
  !
  ! globals
  !
  type (rinf1_input) :: inp
  type (rinf1_output) :: res
  !
  ! locals
  !
  integer :: i, istr, istrrt, j, it, llen, lrt, ndiff, ioff, ioff1, ioff2, ioff3
  integer (kind=i_8) k1,k2
  integer :: irandom, irandomet
  external :: irandom, irandomet
  double precision :: s1, s2, s3, s4, s5, s6
  double precision :: MPI_WTIM
  external MPI_WTIM

  ! the following must stay in to prevent over-optimization of loop 38
  ! ipo is not allowed!
  common / b / k1,k2
  !
  ! executable statements
  !
  if (rinf1_debug) then
     write(6,fmt=*) ' Calling doall with '
     write(6,fmt=*) ' case number       ',inp%icase
     write(6,fmt=*) ' vector length     ',inp%nvec
     write(6,fmt=*) ' stride            ',inp%istride
  end if
  ndiff = 0
  do i=1,6
     res%nops(i) = 0
  end do
  do i=1,2
     res%nbytes(i) = 0
  end do
  llen = inp%nvec

#ifdef STRIDED
  inp%istride = 0
  do it=1 STRIDED,1
     inp%istride = inp%istride + 1
  end do
#endif

  istr = inp%istride
  !
  select case(inp%icase)
     !
     case(1)
        res%label = 'REAL_8: DYADS: a=b*c'
        res%nops(4) = 1
        res%nbytes(1) = 2 * 8
        res%nbytes(2) = 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_1
OPTION_1
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              a(i) = b(i) * c(i)
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(2)
        res%label = 'REAL_8: TRIADS: a=b*c+d'
        res%nops(4) = 2
        res%nbytes(1) = 3 * 8
        res%nbytes(2) = 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!no$OMP do schedule(static,omp_chunk)
!$OMP do simd
#ifdef OPTION_2
OPTION_2
#endif
#ifdef OPTION_2A
OPTION_2A
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
              do i=1,llen STRIDED
#endif
              a(i) = d(i) + b(i) * c(i) 
           end do
!$OMP end do simd
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(3)
        res%label = 'REAL_8: 4-OPS: a=b*c+d*e+f'
        res%nops(4) = 4
        res%nbytes(1) = 5 * 8
        res%nbytes(2) = 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel  private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_3
OPTION_3
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              a(i) = b(i) * c(i) + d(i) * e(i) + f(i)
           end do
!$OMP end do
        end do
!$OMP end parallel 
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
     !
     case(4)
        res%label = 'REAL_8: DOT SCALAR_PRODUCT: s+=b*c'
        res%nops(4) = 1
        res%nbytes(1) = 2 * 8
        res%nbytes(2) = 0 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
        s = 0.0d0
!$OMP parallel private(i,it) reduction(+:s)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_4
OPTION_4
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              s = s + b(i) * c(i)
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        res%testsum = s
        !
        !
     case(5)
        res%label = 'REAL_8_INT_8: RANDOM_GATHER: a=c(ib)'
        res%nops(4) = 1
        res%nbytes(1) = 1 * 8 + 1 * 8
        res%nbytes(2) = 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_5
OPTION_5
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              a(i) = c(ib(i))
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(6)
        res%label = 'REAL_8_INT_8: RANDOM_SCATTER: b(ib)=c'
        res%nops(4) = 1
        res%nbytes(1) = 1 * 8
        res%nbytes(2) = 1 * 8 + 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_6
OPTION_6
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              b(ib(i)) = c(i)
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(7)
        res%label = 'REAL_8: FIRST ORDER RECURRENCE: a(i)=b(i)*a(i-1)+d(i)'
        res%nops(4) = 2
        res%nbytes(1) = 2 * 8
        res%nbytes(2) = 1 * 8
! Load of a not counted since available from previous iteration
        ndiff = 1
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
        do it=1,inp%nit
           call dummy()
           !  OpenMP parallelization not possible without change to code
#ifdef OPTION_7
OPTION_7
#endif
#ifndef STRIDED
           do i=2,llen,istr
#else
           do i=2,llen STRIDED
#endif
!             write(6,*) 'i,a,b,d: ',i,a(i-1),b(i),d(i)
              a(i) = b(i)*a(i-1)+d(i)
           end do
        end do
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(8)
        res%label = 'REAL_8_INT_8: CHARGE ASSIGNMENT: a(ia) = a(ia)+const'
        res%nops(4) = 1
        res%nbytes(1) = 1 * 8 + 1 * 8
        res%nbytes(2) = 0 * 8
! store of c not counted since written back from load
        ndiff = 1
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
        s = 1.7349d0
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP  do
#ifdef OPTION_8
OPTION_8
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              c(ia(i)) = c(ia(i)) + s
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(9)
        res%label = 'REAL_8: DAXPY: a=const*b+c'
        res%nops(4) = 2
        res%nbytes(1) = 2 * 8
        res%nbytes(2) = 1 * 8
        s = 1.7349d0
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_9
OPTION_9
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              a(i) = s * b(i) + c(i)
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(10)
        res%label = 'REAL_8_INT_8: INDIRECT DAXPY: a(ic)=const*b(ib)+c(ic)'
        res%nops(4) = 2
        res%nbytes(1) = 2 * 8 + 2 * 8
        res%nbytes(2) = 1 * 8
        s = 1.7349d0
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_10
OPTION_10
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              a(ic(i)) = s * b(ib(i)) + c(ic(i))
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
     !
     case(11)
        res%label = 'INT_8: ADD: ia=ib+ic'
        res%nops(1) = 1
        res%nbytes(1) = 2 * 8
        res%nbytes(2) = 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_11
OPTION_11
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              ia(i) = ib(i) + ic(i)
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
     !
     case(12)
        res%label = 'INT_8: MULT: ia=ib*ic'
        res%nops(1) = 1
        res%nbytes(1) = 2 * 8
        res%nbytes(2) = 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_12
OPTION_12
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              ia(i) = ib(i) * ic(i)
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
     !
     case(13)
        res%label = 'INT_8: TRIADS: ia=ib*ic+id'
        res%nops(1) = 2
        res%nbytes(1) = 3 * 8
        res%nbytes(2) = 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_13
OPTION_13
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              ia(i) = ib(i) * ic(i) + id(i)
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
     !
     case(14)
        res%label = 'ALL_TYPES: VARIOUS INTRINSICS: int,dble,cmplx,mod'
        res%nops(1) = 2
        res%nops(4) = 2 
        res%nbytes(1) = 1 * 4 + 2 * 8 + 1 * 16
        res%nbytes(2) = 2 * 4 + 2 * 8 + 1 * 16
! b not loaded, ia not loaded
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_14
OPTION_14
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              ia(i) = int(a(i))
              b(i)  = dble(ib(i))
              c(i)  = dble(cc(i))
              ca(i) = cmplx(b(i),c(i))
              ib(i) = mod(ia(i),10_i_8)
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(15)
        res%label = 'REAL_8: TRIGO_EXP: a=acos(sin(exp(b)))'
        res%nops(4) = 120
        !  THIS IS ONLY TRUE FOR CRAY T90 BUT DO NOT CHANGE THIS NUMBER 
        res%nbytes(1) = 1 * 8
        res%nbytes(2) = 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_15
OPTION_15
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              a(i) = acos(sin(exp(b(i))))
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
     !
     case(16)
        res%label = 'REAL_8: NORM: c=sqrt(a*a+b*b)'
        res%nops(4) = 4
        res%nbytes(1) = 2 * 8
        res%nbytes(2) = 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_16
OPTION_16
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              c(i) = sqrt(a(i)*a(i) + b(i)*b(i))
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
     !
     case(17)
        res%label = 'REAL_8: DIVIDE: a=b/c'
        res%nops(4) = 1
        res%nbytes(1) = 2 * 8
        res%nbytes(2) = 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_17
OPTION_17
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              a(i) = b(i)/c(i)
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
     !
     case(18)
        res%label = 'REAL_8: RANDOM_GEN: random_number(a)'
        res%nops(4) = 4
        !  THIS IS ONLY TRUE FOR CRAY T90 BUT DO NOT CHANGE THIS NUMBER 
        res%nbytes(1) = 0
        res%nbytes(2) = 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              call random_number(a(i))
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(19)
        res%label = 'COMPLEX_16: TRIADS: ca=cb*cc+cd'
        res%nops(6) = 2
        res%nbytes(1) = 3 * 16
        res%nbytes(2) = 1 * 16
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_19
OPTION_19
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              ca(i) = cb(i) * cc(i) + cd(i)
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(20)
        ! intrinsic maxval has been moved to Fortran 90 specific code
        res%label = 'REAL_8_INT_8: MAXVAL: maxval()'
        res%nops(1) = 1
        res%nops(4) = 1
        res%nbytes(1) = 1 * 8 + 1 * 8
        res%nbytes(2) = 0
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
           ib(1) = ia(1)
           b(1) = a(1)
!$OMP do
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              if (ib(1) .lt. ia(i)) ib(1) = ia(i)
              if (b(1) .lt. a(i)) b(1) = a(i)
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(21)
        res%label = 'REAL_8: LOOP WITH IFS: b = c +/- d'
        res%nops(4) = 1
        res%nbytes(1) = 3 * 8
        res%nbytes(2) = 1 * 8
        s = acos(-1.0d0)/100.0d0
        do i=1,llen,istr
           a(i) = sin(s*mod(i,1000))
        end do
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_21
OPTION_21
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              if (a(i) >= 0.5) then
                 b(i) = c(i) + d(i)
              else if (a(i) >= 0.0 .and. a(i) < 0.5) then
                 b(i) = c(i) - d(i)
              else if (a(i) >= -0.5 .and. a(i) < 0.0) then
                 b(i) = c(i) + e(i)
              else
                 b(i) = c(i) - e(i)
              end if
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(22)
        res%label = 'REAL_8: MESSY LOOP: all floats '
        res%nops(4) = 40
        res%nbytes(1) = 5 * 8
        res%nbytes(2) = 1 * 8
        call setzero(f,nnmax)
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
        do it=1,inp%nit
           call dummy()
#ifdef OPTION_22
OPTION_22
#endif
#ifndef STRIDED
           do i=2,inp%nvec-1,inp%istride
#else
           do i=2,llen-1 STRIDED
#endif
!             write(6, *) 'i = ',i
             s1 = a(I)*b(I-1)
             s2 = a(I)**1.01d0 - c(I)**0.99d0
             s3 = s1*s2 + b(I+1)
             d(I) = s1+s2*s3 + 1.
             s1 = e(I)*c(I) + b(I+1)
             s2 = abs(s3 * s1 - e(I))
             s1 = (a(I) + b(I))*s2
             e(i) = a(i) + c(I)*b(I)*I + sqrt(s1**1.02) -1.
             S4=min(a(i),b(i),d(i))
             s5=max(a(i),b(i),e(i))
             if (s5.ge.0) then
                S6=s4**1.003 + s5**1.005
             else
                s6=1.
             endif
             if (s4.ge.0) s4=max(1.D0,sqrt(s4))
             d(I) = D(i) + s6 - s4
!             write(6, *) 'before line: ',f(i-1),a(i+1),b(i+1),c(i),d(i)
!             write(6, *) '...        : ',s1,s2,s5,s6
             f(i) = f(i)+f(i-1)*(a(I+1)*s1+b(i+1)*s2+c(i)*s5+d(i)*s6)
!             write(6, *) 'after line '
          end do
      end do
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        ndiff = 2
        !
     case(23)
        res%label = 'REAL_8: ADD_UNROLL_4: a=b+c'
        res%nops(4) = 1
        res%nbytes(1) = 2 * 8
        res%nbytes(2) = 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_23
OPTION_23
#endif
           do i=1,inp%nvec-3,4
              a(i) = b(i) + c(i)
              a(i+1) = b(i+1) + c(i+1)
              a(i+2) = b(i+2) + c(i+2)
              a(i+3) = b(i+3) + c(i+3)
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        ndiff = inp%nvec - i + 1 
        !
     case(24)
        res%label = 'REAL_8: SCALE_UNROLL_4: a=s*c'
        res%nops(4) = 1
        res%nbytes(1) = 1 * 8
        res%nbytes(2) = 1 * 8
        s = 1.0001d0
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_24
OPTION_24
#endif
           do i=1,inp%nvec-3,4
              a(i) = s * c(i)
              a(i+1) = s * c(i+1)
              a(i+2) = s * c(i+2)
              a(i+3) = s * c(i+3)
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        ndiff = inp%nvec - i + 1 
        !
     case(25)
        res%label = 'REAL_8: 4_OP_UNROLL_4: a=b*c-d*e'
        res%nops(4) = 3
        res%nbytes(1) = 4 * 8
        res%nbytes(2) = 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_25
OPTION_25
#endif
           do i=1,inp%nvec-3,4
              a(i) = b(i)*c(i) - d(i)*e(i)
              a(i+1) = b(i+1)*c(i+1) - d(i+1)*e(i+1)
              a(i+2) = b(i+2)*c(i+2) - d(i+2)*e(i+2)
              a(i+3) = b(i+3)*c(i+3) - d(i+3)*e(i+3)
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        ndiff = inp%nvec - i + 1 
        !
     case(26)
        res%label = 'REAL_8: RECURSION_UNROLL_4:  '
        res%nops(4) = 2
        res%nbytes(1) = 2 * 8
        res%nbytes(2) = 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel  private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_26
OPTION_26
#endif
           do i=2,inp%nvec-3,4
              a(i) = a(i-1)*b(i) + c(i)
              a(i+1) = a(i)*b(i+1) + c(i+1)
              a(i+2) = a(i+1)*b(i+2) + c(i+2)
              a(i+3) = a(i+2)*b(i+3) + c(i+3)
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        ndiff = inp%nvec - i + 2
        !
     case(27)
        res%label = 'REAL_4: DYADS: a=b*c'
        res%nops(3) = 1
        res%nbytes(1) = 2 * 4
        res%nbytes(2) = 1 * 4
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_27
OPTION_27
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              sa(i) = sb(i) * sc(i)
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(28)
        res%label = 'REAL_8: SPARSE_JAGGED_FORMAT: v=M*v'
        res%nops(4) =   100 * 2
        res%nbytes(1) = 100 * (1 * 4 + 3 * 8)
        res%nbytes(2) = 100 * 8
        call sparse_arrays1(inp)
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(ioff,i,it)
        do it=1,inp%nit
           ioff = it*5.1234
           call dummy()
           do j=1,100
              ioff = mod(ioff+111+j,nnmax)+1
              ioff = max(0,min(nnmax-inp%nvec,ioff-10))
#ifdef OPTION_28
OPTION_28
#endif
!$OMP do
#ifndef STRIDED
              do i=1,inp%nvec
#else
              do i=1,llen
#endif
                 a(i) = a(i) + b(ioff+i)*c(ia(ioff+i))
              end do
!$OMP end do
           end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(29)
        res%label = 'REAL_8: SPARSE_JAGGED_FORMAT_UNROLL_4 : v=M*v'
        it = mod(inp%istride,4) 
        if (inp%istride < 4) it = it - 4
        inp%istride = inp%istride - it
        res%nops(4) = inp%istride * 2
        res%nbytes(1) = inp%istride * (1 * 4 + 3 * 8)
        res%nbytes(2) = inp%istride * 8
        call sparse_arrays1(inp)
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(it,i,j,ioff,ioff1,ioff2,ioff3)
        do it=1,inp%nit
           ioff = it*5.1234
           ioff1 = ioff
           ioff2 = ioff
           ioff3 = ioff
           call dummy()
!$OMP do
           do j=1,inp%istride-3,4
              ioff = mod(ioff+111+j,nnmax)+1
              ioff = max(0,min(nnmax-inp%nvec,ioff-10))
              ioff1 = mod(ioff1+112+j,nnmax)+1
              ioff1 = max(0,min(nnmax-inp%nvec,ioff1-10))
              ioff2 = mod(ioff2+113+j,nnmax)+1
              ioff2 = max(0,min(nnmax-inp%nvec,ioff2-10))
              ioff3 = mod(ioff3+114+j,nnmax)+1
              ioff3 = max(0,min(nnmax-inp%nvec,ioff3-10))
#ifdef OPTION_29
OPTION_29
#endif
#ifndef STRIDED
              do i=1,inp%nvec
#else
              do i=1,llen
#endif
                 a(i) = a(i) + b(ioff+i)*c(ia(ioff+i))  &
                      + b(ioff1+i) * c(ia(ioff1+i)) &
                      + b(ioff2+i) * c(ia(ioff2+i)) &
                      + b(ioff3+i) * c(ia(ioff3+i))
              end do
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(30)
        res%label = 'REAL_8: SPARSE_FEM: v=M*v'
        res%nops(4) = inp%istride * 2
        res%nbytes(1) = inp%istride * (3 * 4 + 3 * 8)
        res%nbytes(2) = inp%istride * 8
        call sparse_arrays2(inp)
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(ioff,i,j)
        do it=1,inp%nit
           ioff = mod(47,inp%nvec)+1
           ioff  = max(0,min(nnmax-inp%nvec,ioff-10))
           call dummy()
! No repeated indices, so the following is allowed
!DIR$   IVDEP
!OCL  NOVREC
!$OMP do
#ifdef OPTION_30
OPTION_30
#endif
#ifndef STRIDED
           do j=1,inp%nvec
#else
           do j=1,llen
#endif
              i = j + ioff
              a(ia(i)) = a(ia(i)) + b(ib(i))*c(ic(i))
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(31)
        res%label = 'REAL_8: Square Root: a=SQRT(b)'
        res%nops(4) = 1
        res%nbytes(1) = 1 * 8
        res%nbytes(2) = 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_31
OPTION_31
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              a(i) = sqrt(b(i)) 
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(32)
        res%label = 'REAL_8_INT_8: RANDOM_GATHER_UNROLL4: a=c(ib)'
        res%nops(4) = 1
        res%nbytes(1) = 1 * 8 + 1 * 8
        res%nbytes(2) = 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_32
OPTION_32
#endif
           do i=1,inp%nvec-3,4
              a(i) = c(ib(i))
              a(i+1) = c(ib(i+1))
              a(i+2) = c(ib(i+2))
              a(i+3) = c(ib(i+3))
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        ndiff = inp%nvec - i + 2
        !
     case(33)
        res%label = 'REAL_8_INT_8: RANDOM_SCATTER_UNROLL4: a=c(ib)'
        res%nops(4) = 1
        res%nbytes(1) = 1 * 8 + 1 * 8
        res%nbytes(2) = 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_33
OPTION_33
#endif
           do i=1,inp%nvec-3,4
              b(ib(i)) = c(i)
              b(ib(i+1)) = c(i+1)
              b(ib(i+2)) = c(i+2)
              b(ib(i+3)) = c(i+3)
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        ndiff = inp%nvec - i + 2
        !
        !
     case(34)
        res%label = 'REAL_8_INT_8: RANDOM_GATHER_NOPERM: b=c(ib)'
        res%nops(4) = 1
        res%nbytes(1) = 1 * 8 + 1 * 8
        res%nbytes(2) = 1 * 8
        do i=1,llen
           ib(i)=i
        end do
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_34
OPTION_34
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              b(i) = c(ib(i))
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
     case(35)
        res%label = 'REAL_8_INT_8: CONDITIONAL_INDEX: if a>b ib(i+)=index'
        res%nops(4) = 1
        res%nbytes(1) = 2 * 8
        res%nbytes(2) = 0.5 * 8
        do i=1,llen
           a(i) = b(i) + (-1)**i 
        end do
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
           j = 0
!$OMP do
#ifdef OPTION_35
OPTION_35
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              if (a(i)>b(i)) then
                 j = j+1
                 ib(j) = i
              end if
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(36)
        res%label = 'REAL_8: CONDITIONAL_INDEX: if a>0 b=a'
        res%nops(4) = 1
        res%nbytes(1) = 2 * 8
        res%nbytes(2) = 0
        do i=1,llen
           a(i) =  (-1)**i 
        end do
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_36
OPTION_36
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
           do i=1,llen STRIDED
#endif
              if (a(i)>0.0) then
                 b(i) = a(i)
              end if
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(37)
        res%label = 'REAL_8: COMPLEX IF LOOP:  s = sum th(b) x th(a)'
        res%nops(4) = 2
        res%nbytes(1) = 2 * 8
        res%nbytes(2) = 0
        do i=1,llen
           a(i) = (-1)**i 
           b(i) = - a(i)
        end do
        lrt = int(sqrt(real(llen)+0.1))
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
        s = 0.0d0
        d = 0.0d0
!$OMP parallel private(i,j,it) reduction (+:s)
        do it=1,inp%nit
           call dummy()
!$OMP do 
#ifdef OPTION_37
OPTION_37
#endif
           do i=1,lrt,istr
              do j=1,lrt,istr
                 if (a(i*j) > 0.0d0 .and. b(i*j) > 0.0d0) then
                    d(i + istr*(j-1)) =  b(i*j)*a(i*j)
                 end if
              end do
           end do
!$OMP end do
           s = s + sum(d)
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(38)
        res%label = 'INT_8: weighted inverse Latency : Random read'
        res%nops(4) = 1 
! counted as single op for perf. evaluation
        res%nbytes(1) = 1 * 8
        res%nbytes(2) = 0
!        llen = 200
!        i = irandomet(123456789,362436069,521288629,916191069)
!        do i=1,llen 
!           k1 = irandom()
!           if (k1 <= 0) k1 = -k1
!           ia(i) = mod(int(k1),llen)+1
!test           ia(i) = abs(dble(k1)*(dble(llen)/dble(huge(ia(1)))))
!test           write(6,*) ' ia = ', ia(i)
!        end do
! try to enforce chain length llen
        ib(2:llen) = 0
        ib(1) = 1
        k2 = 1
        do i=1,llen
           k1 = abs(irandom())
           k1 = mod(int(k1),llen) + 1
           do while (ib(k1) /= 0 .and. i /= llen)
              k1 = k1 + 1
              if (k1 > llen) k1 = k1 - llen
           end do
           ib(k1) = 1
           ia(k2) = k1
           if (i == llen) ia(k2) = 1
           k2 = k1
        end do
        s1 = 0.0d0
        s2 = 0.0d0
        s3 = 0.0d0
        k1 = ia(1)  
        do i=2,llen 
          s1 = s1 + dble(ia(k1) - k1)
          s2 = s2 + abs(dble(ia(k1) - k1))
          s3 = s3 + abs(dble(ia(k1) - k1)**2)
          k1 = ia(k1)
        end do
        s1 = s1 / dble(llen)
        s2 = s2 / dble(llen)
        s3 = s3 / dble(llen)
        s3 = sqrt(dble(llen)/dble(llen-1)*(s3 - s2**2))
!        stop
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,j,k1,k2,it) 
        do it=1,inp%nit
           call dummy()
           k1 = ia(1)
!$OMP do 
#ifdef OPTION_38
OPTION_38
#endif
           do i=2,llen
              k2 = ia(k1)
              k1 = k2
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(39)
        res%label = 'REAL_4: TRIADS: sa=sb*sc+sd'
        res%nops(4) = 2
        res%nbytes(1) = 3 * 4
        res%nbytes(2) = 1 * 4
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_39
OPTION_39
#endif
#ifdef OPTION_39A
OPTION_39A
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
              do i=1,llen STRIDED
#endif
              sa(i) = sd(i) + sb(i) * sc(i)
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(40)
        res%label = 'REAL_8: SNORM: s=s+a*a'
        res%nops(4) = 2
        res%nbytes(1) = 1 * 8
        res%nbytes(2) = 0 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
        s = 0.0d0
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_40
OPTION_40
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
              do i=1,llen STRIDED
#endif
              s = s + a(i) * a(i)
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(41)
        res%label = 'REAL_8: SCALE: a=k*a + s'
        res%nops(4) = 2
        res%nbytes(1) = 1 * 8
        res%nbytes(2) = 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
        s1 = 1.0d0
        s = 1.0d0
!$OMP parallel private(i,it)
        do it=1,inp%nit
           call dummy()
!$OMP do
#ifdef OPTION_40
OPTION_40
#endif
#ifndef STRIDED
           do i=1,llen,istr
#else
              do i=1,llen STRIDED
#endif
              a(i) = s + s1 * a(i)
           end do
!$OMP end do
        end do
!$OMP end parallel
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !


     case default
        write(rinf1_out,fmt=*) 'doall: Case ',inp%icase,' not implemented'
        call rinf1_abort('       Please add case in doall.F')
  end select
  !
  ! compute timing results etc  
  !
  res%time =  (t2 - t1)/dble(inp%nit) - t0
  res%perf = 0
  do i=1,6
     res%nops(i) = (res%nops(i) * (inp%nvec - ndiff)) / inp%istride
     res%perf = res%perf + res%nops(i)
  end do
  res%perf = res%perf / (1.0D6 * res%time)

  if (rinf1_debug) then
      write(6,'(2x, A,I12,4(es12.5,1x))') 'nops,perf,time/it,time :', sum(res%nops(1:6)),res%perf, res%time,  res%time*dble(inp%nit)
      write(6,*)
  endif


  do i=1,2
     res%nbytes(i) = (res%nbytes(i) * (inp%nvec - ndiff)) / inp%istride
  end do
  res%bw = (res%nbytes(1) + res%nbytes(2))/(1.0D6 * res%time)

  if (inp%icase == 38) then
     open(unit=98,file='latency.res',position='APPEND',status='UNKNOWN')
     write(98,'(i10,1x,4(es12.5,1x))') llen,s1,s2,s3,res%time
     close(98)
  end if

  return
end subroutine doall




