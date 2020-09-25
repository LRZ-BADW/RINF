!-*- f90 -*-
program rinf1
  use rinf1_admin
  use rinf1_wp
  implicit none
#ifdef USE_MPI
  include 'mpif.h'
#endif
  type (rinf1_input) :: input
  type (rinf1_output) :: results
  integer :: ind,ival,nc,nval
  integer, allocatable :: nit_send(:)
  double precision :: xiter(3)
  double precision :: average,tval,xnew,dt,tnew,t_all,t_doall,t_doall0
  character(6) :: xxx
  !
  double precision :: MPI_WTIM
  external MPI_WTIM
  !
  call rinf1_io_init()
  call rinf1_init()
  !
  ! Now start loop over all instances
  !
  nc = 0
  t_all = MPI_WTIME()
  t_doall = 0.0

  do while(.true.)
     input%icase=1 
     input%istride=1
     input%realloc=.TRUE.
     !input%realloc='L'
     read(rinf1_aux, fmt=*, end=200, err=100) input%icase, input%istride, input%realloc, input%mode
     nc = nc + 1
     !
     ! First call:
     !  get label for output
     !  check values 
     !  calibrate first iteration number
     !  
     call select_lengths(1,input%nvec,rinf1_abort,nvecmax=nnmax,num_lengths=num_lengths)

     results%time = 0.0d0
     input%nit = 1

     do while (results%time * input%nit < 1.0 )
        input%nit = input%nit * 10
        call rinf1_arrays(input)
#ifdef FORTRAN77
        call doall(input,results)
#endif
#ifdef FORTRAN90
        call doall90(input,results)
#endif
#ifdef M_C
        results%label=''
        call doallc(input,results)
#endif
        if (rinf1_debug) write(6,fmt=*) 'rinf1 Callibration: Time ',results%time * input%nit,' with nit = ',input%nit
     end do

     if (mp_id == 0) then
        write(rinf1_out,fmt=12) input%icase,results%label
        write(rinf1_fout,fmt=12) input%icase,results%label
12      format(/,/,'(',I2,')  ',A70)     
        write(rinf1_out,fmt=13) 
        write(rinf1_fout,fmt=13)
13      format(/, &
             '      veclen',&
             '      stride',&
             '        iter',&
             '  time/iter.',&
             '  avg. Perf.',&
             '  min. Perf.',&
             '     avg. BW',&
             '     min. BW',&
             '      maxdev',&
                      /,36x,& 
             '           s',&
             '     MFlop/s',&
             '     MFlop/s',&
             '     MByte/s',&
             '     MByte/s',&
             '           %' )
     end if
     input%nit = int(tconst / results%time)
     tval = tconst
     input%niter(1) = input%nvec
     input%titer(1) = results%time
     !
     if (rinf1_debug) then
        write(rinf1_out,'(''rinf1:    # Vectors   '',I5)') num_lengths
        write(rinf1_out,'(''rinf1:    Line '',I5,'' from rinf1.conf '')') nc
        write(rinf1_out,'(''rinf1:    Loop Type   '',I5)') input%icase
        write(rinf1_out,'(''rinf1:    Loop Stride '',I5)') input%istride
     end if

     if (mp_id == 0) then
        call rinf1_wp_init(input%mode,nnmax)
     end if
     !
     ! various vector lengths
     !
     do ind=1,num_lengths

        if (rinf1_debug .and. nc == 1) then
           write(rinf1_out,fmt=*) 
           write(rinf1_out,fmt=*) 'rinf1: Length for run ',ind,' is ',input%nvec
        end if
 
        ! extrapolate time for single iteration to new vector length
        ! and recalculate iteration count for constant total time
        !
        if (ind <= 2) then
           nval = ind
        else if (ind >= 3) then
           nval = 3
        end if
        xiter = input%niter
        xnew = input%nvec
        call polint(xiter,input%titer,nval,xnew,tnew,dt)
        if (tnew < 0.0d0) tnew = input%titer(3)
        input%nit = tval/tnew
        if (input%nit < 5) input%nit = 5
        ! prevent too many longish 1-iterations
        if (input%nit /= 1) then
           input%nitisone = 0
        else
           input%nitisone = input%nitisone + 1
        end if
        ! actually we use the results from node 0 and scatter those
#ifdef USE_MPI
        allocate(nit_send(1:mp_processes),stat=all_stat)
        nit_send = input%nit
        call mpi_scatter(nit_send,1,mpi_integer,input%nit,1,mpi_integer,0,mpi_comm_world,mp_err)
        deallocate(nit_send)
#endif
        if (rinf1_debug) then
           write(6,fmt=*) ' Iterations on proc ',mp_id,': ',input%nit
        end if
        !
        ! initialize arrays and measure
        !
        call rinf1_arrays(input)



#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
#endif
        t_doall0 = MPI_WTIME()
#ifdef FORTRAN77
        call doall(input,results)
#endif
#ifdef FORTRAN90
        call doall90(input,results)
#endif
#ifdef M_C
        results%label=''
        call doallc(input,results)
#endif
#ifdef USE_MPI
        call mpi_barrier(mpi_comm_world,mp_err)
        t_doall = t_doall + MPI_WTIME() - t_doall0
#endif

        !
        ! Write out results
        !
        ! fixme: MPI case, calculate deviations
#ifdef USE_MPI

        results%perf=results%perf/dble(omp_threads)
        results%bw =results%bw   /dble(omp_threads)

        call mpi_gather(results%perf,1,mpi_double_precision,&
             allperfs,1,mpi_double_precision,0,mpi_comm_world,mp_err)
        call mpi_gather(results%bw,1,mpi_double_precision,&
             allbws,1,mpi_double_precision,0,mpi_comm_world,mp_err)
#else
        allperfs = results%perf/dble(omp_threads)
        allbws   = results%bw/dble(omp_threads)
#endif
        !
        if (mp_id == 0) then

           results%perf = sum(allperfs(1:mp_processes))/dble(mp_processes)
           results%perfmin = minval(allperfs(1:mp_processes)) 

           call rinf1_wp_addval(input%nvec,results%perfmin)

           if (results%perf > 1.0D-10) then
              results%perfdev = maxval(abs(allperfs(1:mp_processes)-results%perf))/results%perf*100.0D0
           else
              results%perfdev = 0.0D0
           end if
           results%bw = sum(allbws(1:mp_processes))/dble(mp_processes)
           results%bwmin = minval(allbws(1:mp_processes))
           write(rinf1_out ,fmt=112) input%nvec,input%istride,input%nit,&
                results%time, &
                results%perf,results%perfmin,&
                results%bw,  results%bwmin,&
                results%perfdev
           write(rinf1_fout,fmt=112) input%nvec,input%istride,input%nit,&
                results%time, &
                results%perf,results%perfmin,&
                results%bw,  results%bwmin,&
                results%perfdev
112        format(3(1x,I11),5(1X,1PE11.4),3(1X,EN11.1),1X,F11.2)
        end if

        !
        ! Update for next step
        !
        if (ind.lt.num_lengths) then
           call select_lengths(ind+1,input%nvec,rinf1_abort)
        end if
        if (nval < 3) then
           input%niter(nval + 1) = input%nvec
           input%titer(nval + 1) = results%time
        else
           do ival=1,2
              input%niter(ival) = input%niter(ival+1)
              input%titer(ival) = input%titer(ival+1)
           end do
           input%niter(3) = input%nvec
           input%titer(3) = results%time
        end if

     end do

     if (mp_id == 0) then
!       average over alle loops
        average = rinf1_wp_finalize()
        write(rinf1_out,fmt=123)  input%icase,input%istride,average
        write(rinf1_fout,fmt=123) input%icase,input%istride,average
        write(rinf1_out,fmt=124)  input%icase,input%istride,average*omp_threads
        write(rinf1_fout,fmt=124) input%icase,input%istride,average*omp_threads
        write(rinf1_out,fmt=125)  input%icase,input%istride,average*omp_threads*mp_processes
        write(rinf1_fout,fmt=125) input%icase,input%istride,average*omp_threads*mp_processes

123     format('===>> Case: ',i2,'  Stride:',i2,'  Weighted Performance Value PER THREAD (MFlop/s): ',1pe12.5)
124     format('===>> Case: ',i2,'  Stride:',i2,'  Weighted Performance Value PER TASK   (MFlop/s): ',1pe12.5)
125     format('===>> Case: ',i2,'  Stride:',i2,'  Weighted Performance Value PER JOB    (MFlop/s): ',1pe12.5,/)
     end if

  end do

100 continue
200 continue



  t_all = MPI_WTIME() - t_all
  if (mp_id == 0) then
     write(rinf1_out,'(''rinf1 total time: '',F12.3)') t_all
     write(rinf1_out,'(''rinf1 loop time : '',F12.3)') t_doall
     write(rinf1_out,'(''rinf1 init time : '',F12.3)') t_all - t_doall
  end if
  call rinf1_finalize()
end program rinf1










