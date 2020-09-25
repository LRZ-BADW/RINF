!-*- f90 -*-
subroutine doall90(inp, res)
  !
  ! doall contains all loops in Fortran 90 style. 
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
  integer :: i, j, it, llen, lrt, ndiff, ioff, ioff1, ioff2, ioff3, istr
  double precision :: s1, s2, s3, s4, s5, s6
  double precision :: MPI_WTIM
  external MPI_WTIM
  !
  ! executable statements
  !
  if (rinf1_debug) then
     write(6,fmt=*) 'Calling doall90 with '
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
  do it=1 STRIDED
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
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_1
OPTION90_1
#endif
           a(1:llen:istr) = b(1:llen:istr) * c(1:llen:istr)
        end do
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
!$omp parallel private(it)
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_2
OPTION90_2
#endif
!$omp workshare
           a(1:llen:istr) = b(1:llen:istr) * c(1:llen:istr) + d(1:llen:istr) 
!$omp end workshare
        end do
!$omp end parallel
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
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_3
OPTION90_3
#endif
           a(1:llen:istr) = b(1:llen:istr) * c(1:llen:istr) + d(1:llen:istr)*e(1:llen:istr) + f(1:llen:istr)
        end do
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
     !
     case(4)
        res%label = 'REAL_8: SCALAR_PRODUCT: s+=b*c'
        res%nops(4) = 1
        res%nbytes(1) = 2 * 8
        res%nbytes(2) = 0 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_4
OPTION90_4
#endif
           s = dot_product(b(1:llen:istr),c(1:llen:istr))
        end do
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        res%testsum = s
        !
        !
        ! cases (5),(6),(7),(8) removed. not possible with array syntax
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
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_9
OPTION90_9
#endif
           a(1:llen:istr) = s * b(1:llen:istr) + c(1:llen:istr)  
        end do
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
        ! case (10) removed. not possible with array syntax
     case(11)
        res%label = 'INT_8: ADD: ia=ib+ic'
        res%nops(1) = 1
        res%nbytes(1) = 2 * 8
        res%nbytes(2) = 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_11
OPTION90_11
#endif
           ia(1:llen:istr) = ib(1:llen:istr) + ic(1:llen:istr)
        end do
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
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_12
OPTION90_12
#endif
           ia(1:llen:istr) = ib(1:llen:istr) * ic(1:llen:istr)
        end do
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
     !
     case(13)
        res%label = 'INT_8: TRIAD: ia=ib*ic+id'
        res%nops(1) = 2
        res%nbytes(1) = 3 * 8
        res%nbytes(2) = 1 * 8
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_13
OPTION90_13
#endif
           ia(1:llen:istr) = ib(1:llen:istr) * ic(1:llen:istr) + id(1:llen:istr)
        end do
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
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_14
OPTION90_14
#endif
           ia(1:llen:istr) = int(a(1:llen:istr))
           b(1:llen:istr) = dble(ib(1:llen:istr))
           c(1:llen:istr) = dble(cc(1:llen:istr))
           ca(1:llen:istr) = cmplx(b(1:llen:istr),c(1:llen:istr))
           ib(1:llen:istr) = mod(ia(1:llen:istr),10_8)
        end do
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
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_15
OPTION90_15
#endif
           a(1:llen:istr) = acos(sin(exp(b(1:llen:istr))))
        end do
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
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_16
OPTION90_16
#endif
           c(1:llen:istr) = sqrt(a(1:llen:istr)*a(1:llen:istr)+b(1:llen:istr)*b(1:llen:istr))
        end do
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
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_17
OPTION90_17
#endif
           a(1:llen:istr) = b(1:llen:istr)/c(1:llen:istr)
        end do
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
        do it=1,inp%nit
           call dummy()
           call random_number(a(1:llen:istr))
        end do
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
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_19
OPTION90_19
#endif
           ca(1:llen:istr) = cb(1:llen:istr)*cc(1:llen:istr)+cd(1:llen:istr)
        end do
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(20)
        res%label = 'REAL_8_INT_8: MAXVAL: maxval()'
        res%nops(1) = 1
        res%nops(4) = 1
        res%nbytes(1) = 1 * 8 + 1 * 8
        res%nbytes(2) = 0
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
        do it=1,inp%nit
           call dummy()
           ib(1) = maxval(ia(1:llen:istr))
           b(1) = maxval(a(1:llen:istr))
        end do
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
        do i=1,inp%nvec,inp%istride
           a(i) = sin(s*mod(i,1000))
        end do
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_21
OPTION90_21
#endif
           where(a(1:llen:istr) >= 0.5)
               b(1:llen:istr) = c(1:llen:istr) + d(1:llen:istr)
           elsewhere(a(1:llen:istr) >= 0.0 .and. a(1:llen:istr) < 0.5 )
               b(1:llen:istr) = c(1:llen:istr) - d(1:llen:istr)
           elsewhere(a(1:llen:istr) >= -0.5 .and. a(1:llen:istr) < 0.0 )
               b(1:llen:istr) = c(1:llen:istr) + e(1:llen:istr)
           elsewhere(a(1:llen:istr) < -0.5)
               b(1:llen:istr) = c(1:llen:istr) - e(1:llen:istr)
           endwhere
        end do
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
        ! removed case(22)
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
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_23
OPTION90_23
#endif
           a(1:llen-3:4) = b(1:llen-3:4) + c(1:llen-3:4)
           a(2:llen-2:4) = b(2:llen-2:4) + c(2:llen-2:4)
           a(3:llen-1:4) = b(3:llen-1:4) + c(3:llen-1:4)
           a(4:llen:4) = b(4:llen:4) + c(4:llen:4)
        end do
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
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_24
OPTION90_24
#endif
           a(1:llen-3:4) = s*c(1:llen-3:4)
           a(2:llen-2:4) = s*c(2:llen-2:4)
           a(3:llen-1:4) = s*c(3:llen-1:4)
           a(4:llen:4) = s*c(4:llen:4)          
        end do
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
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_25
OPTION90_25
#endif
           a(1:llen-3:4) = b(1:llen-3:4) * c(1:llen-3:4) -  d(1:llen-3:4) * e(1:llen-3:4)
           a(2:llen-2:4) = b(2:llen-2:4) * c(2:llen-2:4) -  d(2:llen-2:4) * e(2:llen-2:4)
           a(3:llen-1:4) = b(3:llen-1:4) * c(3:llen-1:4) -  d(3:llen-1:4) * e(3:llen-1:4)
           a(4:llen:4) = b(4:llen:4) * c(4:llen:4) - d(4:llen:4) * e(4:llen:4)
        end do
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        ndiff = inp%nvec - i + 1 
        ! removed  case(26)
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
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_27
OPTION90_27
#endif
           sa(1:llen:istr) = sb(1:llen:istr) * sc(1:llen:istr)
        end do
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(28)
        res%label = 'REAL_8: SPARSE_JAGGED_FORMAT: v=M*v'
        res%nops(4) = inp%istride * 2
        res%nbytes(1) = inp%istride * (1 * 4 + 3 * 8)
        res%nbytes(2) = inp%istride * 8
        call sparse_arrays1(inp)
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
!$OMP parallel private(ioff,i,it)
        do it=1,inp%nit
           ioff = it*5.1234
           call dummy()
!$OMP do
           do j=1,inp%istride
              ioff = mod(ioff+111+j,nnmax)+1
              ioff = max(0,min(nnmax-inp%nvec,ioff-10))
#ifdef OPTION_28
OPTION_28
#endif
#ifndef STRIDED
              do i=1,inp%nvec
#else
              do i=1,llen
#endif
                 a(i) = a(i) + b(ioff+i)*c(ia(ioff+i))
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
        res%nbytes(1) = inp%istride * (1 * 4 + 3 * 8)
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
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_31
OPTION90_31
#endif
           a(1:llen:istr) = sqrt(b(1:llen:istr))
        end do
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
        res%label = 'REAL_8_INT_8: PACK: c into b'
        res%nops(4) = 1
        res%nbytes(1) = 1 * 8 + 1 * 8
        res%nbytes(2) = 1 * 8
        mask = .false.
        do i=1,llen,istr
           mask(i) = .true.
        end do
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t1 = MPI_WTIME()
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_34
OPTION90_34
#endif
           b(1:llen:istr) = pack(c(1:llen:istr),mask(1:llen:istr))
        end do
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
           do i=1,inp%nvec,inp%istride
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
        res%label = 'REAL_8: CONDITIONAL_INDEX: where a>0 b=a'
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
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_36
OPTION90_36
#endif
           where(a(1:llen:istr)>0.0) b(1:llen:istr)=a(1:llen:istr)
        end do
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case(37)
        res%label = 'REAL_8: FORALL LOOP:  s = sum th(b) x th(a)'
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
        do it=1,inp%nit
           call dummy()
#ifdef OPTION90_37
OPTION90_37
#endif
           forall (i=1:lrt:istr,j=1:lrt:istr,a(i*j)>0.0d0 .and. b(i*j)>0.0d0)
              d(i + istr*(j-1)) = b(i*j)*a(i*j)
           end forall
           s = s + sum(d)
        end do
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t2 = MPI_WTIME()
        !
     case default
        write(rinf1_out,fmt=*) 'doall90: Case ',inp%icase,' not implemented'
        t1 = 0
        t2 = 1
        res%nops = 1
        res%nbytes = 0
  end select
  !
  ! compute timing results etc  
  !
  res%time =  (t2 - t1)/dble(inp%nit) - t0
	write(6,*) 'res%time' , res%time

  res%perf = 0
  do i=1,6
     res%nops(i) = res%nops(i) * ((inp%nvec - ndiff) / inp%istride)
     res%perf = res%perf + res%nops(i)
  end do
!test  write(6,*) ' dat = ',res%perf, inp%nvec, inp%istride, ndiff
  res%perf = res%perf / (1.0D6 * res%time)
  do i=1,2
     res%nbytes(i) = res%nbytes(i) * ((inp%nvec - ndiff) / inp%istride)
  end do
  res%bw = (res%nbytes(1) + res%nbytes(2))/(1.0D6 * res%time)

  return
end subroutine doall90




