MPICC = mpicc
F90PP = cpp -P
MPIF90 = mpif90
F90_MPI_LOADER= mpif90
OMPFLAG= -fopenmp
OBJGEN=-c
F90_OPT3 = -O3 -march=native -Ofast

DEFS=-DALLOCATE -DUSE_MPI -DRINF1_MAXLEN=$${SIZE:-2000000} -DSTRIDED=' '
F_DEFS = -DFORTRAN77  $(DEFS) 
C_DEFS = -DM_C        $(DEFS) 

USEFILE = doall.o
MOD_RINF1 = rinf1_admin.mod rinf1_wp.mod
OBJ_RINF1 = rinf1.o rinf1_admin.o rinf1_wp.o polint.o select_lengths.o dummy.o irandom.o  $(USEFILE) 


all: rinf1.small rinf1.large

rinf1.small:
	rm -f *.o *mod *.f90 rinf1.exe; touch *.F90; export   SIZE=2000000;   make -f $(MAKEFILE_LIST) exe 
	mv rinf1.exe rinf1.small.exe

rinf1.large:
	rm -f *.o *mod *.f90 rinf1.exe; touch *.F90; export SIZE=400000000; make -f $(MAKEFILE_LIST) exe
	mv rinf1.exe rinf1.large.exe


test: rinf1.small	
	cp rinf1.test.conf rinf1.conf; OMP_NUM_THREADS=1   mpiexec -n 1   rinf1.small.exe

run bench:  all
	cp rinf1.small.conf rinf1.conf; OMP_NUM_THREADS=1   mpiexec -n 28  rinf1.small.exe
	cp rinf1.large.conf rinf1.conf; OMP_NUM_THREADS=28  mpiexec -n 1   rinf1.large.exe
	
exe : rinf1.exe


%.o       : %.c
	$(MPICC)   $(OBJGEN) $(C_DEFS)  $(CC_OPT5)   $(OMPFLAG)  $<

%.o %.mod : %.F90
#C-style preprocessor needed
	$(F90PP)   $(F_DEFS) $*.F90 >$*.f90
	$(MPIF90)  $(OBJGEN) $(F_DEFS) $(F90FREE) $(F90_OPT3) $(OMPFLAG) $*.f90
	mv $*.f90 $*.f90.1

time.o : time.c time.body1.c time.body.c
	$(MPICC)   $(OBJGEN) $(C_DEFS)  $(CC_OPT0)   $<

rinf1.exe : $(MOD_RINF1) $(OBJ_RINF1) 
	$(F90_MPI_LOADER) $(OMPFLAG) -o $@  $(OBJ_RINF1) $(MPI_LIBS)  $(EXTRA)


# administrative stuff
install : rinf1.exe
	test -d $(XDIR) || mkdir $(XDIR)
	cp -p rinf1.exe rinf1.conf $(XDIR)

clean:
	rm -f *.o *mod *.f90 *.f90.1 *.exe
	rm -f *.*out *.*err *.res *.jpg

config :
	rinf1_setup
	cp rinf1.large.conf rinf1.conf
	$(MAKE) clean
