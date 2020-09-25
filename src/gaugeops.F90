!-*- f90 -*-
program gaugeops
  !
  ! measure operation performance for all available
  ! Fortran operations
  !
  implicit none
  integer*8 :: itbase
  double precision :: delta
  parameter (delta = 2.0d0)
  parameter (itbase = 100000000)
  !
  logical :: firsttry
  integer*8 :: i,iter,irep,j
  real :: a4,b4,c4
  double precision :: a,b,c
  double precision :: t1,t2,tdeduct,toper,tdeduct4
  double precision :: sum
  double precision :: MPI_WTIM
  external MPI_WTIM
  !
  ! executable statements
  !
  irep = 20
  sum = 0.0d0
  do j=1,irep
     t1 = MPI_WTIME()
     do i=1,itbase
        call setvals(a,b,c)
     end do
     t2 = MPI_WTIME()
     tdeduct = (t2 - t1)
     sum = sum + tdeduct
     write(6,'(''gaugeops: Overhead '',D10.3,'' Acc. Mean '',D10.3)') tdeduct, sum/dble(j)
  end do
  tdeduct = sum / dble(irep)
  iter = itbase
  ! Overhead for REAL*4
  sum = 0.0d0
  do j=1,irep
     t1 = MPI_WTIME()
     do i=1,iter
        call setvals4(a4,b4,c4)
     end do
     t2 = MPI_WTIME()
     tdeduct4 = t2 - t1
     sum = sum + tdeduct4
     write(6,'(''gaugeops: Overhead '',D10.3,'' Acc. Mean '',D10.3)') tdeduct, sum/dble(j)
  end do
  tdeduct4 = sum / dble(irep)
  !
  ! addition
  !
  toper = 0
  do j=1,irep
     t1 = MPI_WTIME()
     do i=1,iter
        call setvals(a,b,c)
        c = a + b
     end do
     t2 = MPI_WTIME()
     toper = toper + t2 - t1
  end do
  toper = toper/dble(irep) - tdeduct
  write(6,fmt='(''gaugeops: Addition took '',F10.3,'' secs, '',F10.3,'' MFlop/s '')') toper, dble(iter)/toper*1.0D-6
  !
  ! addition (real*4)
  !
  toper = 0
  do j=1,irep
     t1 = MPI_WTIME()
     do i=1,iter
        call setvals4(a4,b4,c4)
        c4 = a4 + b4
     end do
     t2 = MPI_WTIME()
     toper = toper + t2 - t1
  end do
  toper = toper/dble(irep) - tdeduct4
  write(6,fmt='(''gaugeops: Addition (4) took '',F10.3,'' secs, '',F10.3,'' MFlop/s '')') toper, dble(iter)/toper*1.0D-6
  !
  ! multiplication
  !
  toper = 0
  do j=1,irep
     t1 = MPI_WTIME()
     do i=1,iter
        call setvals(a,b,c)
        c = a * b
     end do
     t2 = MPI_WTIME()
     toper = toper + t2 - t1
  end do
  toper = toper/dble(irep) - tdeduct
  write(6,fmt='(''gaugeops: Multiply took '',F10.3,'' secs, '',F10.3,'' MFlop/s '')') toper, dble(iter)/toper*1.0D-6
  !
  ! multiplication (real*4)
  !
  toper = 0
  do j=1,irep
     t1 = MPI_WTIME()
     do i=1,iter
        call setvals4(a4,b4,c4)
        c4 = a4 * b4
     end do
     t2 = MPI_WTIME()
     toper = toper + t2 - t1
  end do
  toper = toper / dble(irep) - tdeduct4
  write(6,fmt='(''gaugeops: Multiply (4) took '',F10.3,'' secs, '',F10.3,'' MFlop/s '')') toper, dble(iter)/toper*1.0D-6
  !
  ! multiply-add
  !
  toper = 0
  do j=1,irep
     t1 = MPI_WTIME()
     do i=1,iter
        call setvals(a,b,c)
        c = a * b + c
     end do
     t2 = MPI_WTIME()
     toper = toper + t2 - t1
  end do
  toper = toper/dble(irep) - tdeduct   
  write(6,fmt='(''gaugeops: Multiply-Add took '',F10.3,'' secs, '',F10.3,'' MFlop/s '')') toper, 2*dble(iter)/toper*1.0D-6
  !
  ! multiply-add (real*4)
  !
  toper = 0
  do j=1,irep
     t1 = MPI_WTIME()
     do i=1,iter
        call setvals4(a4,b4,c4)
        c4 = a4 * b4 + c4
     end do
     t2 = MPI_WTIME()
     toper = toper + t2 - t1
  end do
  toper = toper/dble(irep) - tdeduct4
  write(6,fmt='(''gaugeops: Multiply-Add (4) took '',F10.3,'' secs, '',F10.3,'' MFlop/s '')') toper, 2*dble(iter)/toper*1.0D-6
  !
  ! division
  !
  toper = 0
  do j=1,irep
     t1 = MPI_WTIME()
     do i=1,iter
        call setvals(a,b,c)
        c = a / b
     end do
     t2 = MPI_WTIME()
     toper = toper + t2 - t1
  end do
  toper = toper / dble(irep) - tdeduct
  write(6,fmt='(''gaugeops: Division took '',F10.3,'' secs, '',F10.3,'' MFlop/s '')') toper, dble(iter)/toper*1.0D-6
  !
  ! square root
  !
  toper = 0
  do j=1,irep
     t1 = MPI_WTIME()
     do i=1,iter
        call setvals(a,b,c)
        c = sqrt(a)
     end do
     t2 = MPI_WTIME()
     toper = toper + t2 - t1
  end do
  toper = toper / dble(irep) - tdeduct
  write(6,fmt='(''gaugeops: Sq.Root took '',F10.3,'' secs, '',F10.3,'' MFlop/s '')') toper, dble(iter)/toper*1.0D-6

  stop
end program gaugeops


