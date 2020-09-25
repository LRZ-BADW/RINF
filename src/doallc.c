#include "doallc.h"

#ifdef USE_MPI
#include "mpi.h"
int mp_err;
#endif

#ifdef _OPENMP
#include "omp.h"
#endif



void doallc_(Rinf1_input *inp, Rinf1_output *res) {
    int i, j, it, llen, lrt, ndiff, ioff, ioff1, ioff2, ioff3, istr;
    double s, s1, s2, s3, s4, s5, s6;
    long long int k2;
    long long int * k1, k3;

    ndiff = 0;
    for (i=0;i<6;i++) {
	res->nops[i] = 0;
    }
    res->nbytes[0] = 0;
    res->nbytes[1] = 0;
    llen = inp->nvec;
#ifndef STRIDED
    istr = inp->istride;
#else
#define istr C_STRIDE
#endif
/*    printf (" Now accessing arrays \n"); */
    a = arrays_.a;
    b = arrays_.b;
    c = arrays_.c;
    d = arrays_.d;
    e = arrays_.e;
    f = arrays_.f;
    ia = arrays_.ia;
    ib = arrays_.ib;
    ic = arrays_.ic;
    id = arrays_.id;
    k1 = &b_.k1;
/*    printf ("C: a(2) = %f \n",a[1]);  */
/*    printf ("C: llen = %d \n",llen);  */

/*  Here come all the test kernels */
    switch (inp->icase) {
	case 1:
	    sprintf(res->label, "REAL_8: DYADS: a=b*c");
	    res->nops[3] = 1;
	    res->nbytes[0] = 2 * 8;
	    res->nbytes[1] = 1 * 8;
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t1 = MPI_Wtime();
#pragma omp parallel
	    for (it=0;it<inp->nit;it++) {
		dummy_();
#pragma omp for
		for (i=0;i<llen;i+=istr) {
		    a[i] = b[i] * c[i];
		}
	    }
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t2 = MPI_Wtime();
	    break;
	case 2:
	    sprintf(res->label, "REAL_8: TRIADS: a=b*c+d");
	    res->nops[3] = 2;
	    res->nbytes[0] = 3 * 8;
	    res->nbytes[1] = 1 * 8;
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t1 = MPI_Wtime();
#pragma omp parallel
	    for (it=0;it<inp->nit;it++) {
/*		printf ("C: it = %d \n",it); */
		dummy_();
#pragma omp for
		for (i=0;i<llen;i+=istr) {
/*		    printf ("C: i = %d \n",i); */
		    a[i] = b[i] * c[i] + d[i];
		}
	    }
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t2 = MPI_Wtime();
	    break;
	case 3:
	    sprintf(res->label, "REAL_8: 4-OPS: a=b*c+d*e+f");
	    res->nops[3] = 4;
	    res->nbytes[0] = 5 * 8;
	    res->nbytes[1] = 1 * 8;
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t1 = MPI_Wtime();
#pragma omp parallel
	    for (it=0;it<inp->nit;it++) {
		dummy_();
#pragma omp for
		for (i=0;i<llen;i+=istr) {
		    a[i] = b[i] * c[i] + d[i]* e[i] + f[i];
		}
	    }
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t2 = MPI_Wtime();
	    break;
	case 4:
	    sprintf(res->label, "REAL_8: SCALAR_PRODUCT: s+=b*c");
	    res->nops[3] = 1;
	    res->nbytes[0] = 2 * 8;
	    res->nbytes[1] = 0 * 8;
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t1 = MPI_Wtime();
	    s = 0.0;
#pragma omp parallel
	    for (it=0;it<inp->nit;it++) {
		dummy_();
#pragma omp for
		for (i=0;i<llen;i+=istr) {
		    s = s + b[i] * c[i];
		}
	    }
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t2 = MPI_Wtime();
	    res->testsum = s;
	    break;
	case 5:
	    sprintf(res->label, "REAL_8_INT_8: RANDOM_GATHER: a=c(ib)");
	    res->nops[3] = 1;
	    res->nbytes[0] = 1 * 8 + 1 * 8;
	    res->nbytes[1] = 1 * 8;
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t1 = MPI_Wtime();
#pragma omp parallel
	    for (it=0;it<inp->nit;it++) {
		dummy_();
#pragma omp for
		for (i=0;i<llen;i+=istr) {
		    a[i] = c[ib[i]];
		}
	    }
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t2 = MPI_Wtime();
	    break;
	case 6:
	    sprintf(res->label, "REAL_8_INT_8: RANDOM_SCATTER: b(ib)=c");
	    res->nops[3] = 1;
	    res->nbytes[0] = 1 * 8;
	    res->nbytes[1] = 1 * 8 + 1 * 8;
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t1 = MPI_Wtime();
#pragma omp parallel
	    for (it=0;it<inp->nit;it++) {
		dummy_();
#pragma omp for
		for (i=0;i<llen;i+=istr) {
		    b[ib[i]] = c[i];
		}
	    }
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t2 = MPI_Wtime();
	    break;
	case 7:
	    sprintf(res->label, "REAL_8: FIRST ORDER RECURRENCE: a(i)=b(i)*a(i-1)+d(i)");
	    res->nops[3] = 2;
	    res->nbytes[0] = 2 * 8;
	    res->nbytes[1] = 1 * 8;
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t1 = MPI_Wtime();
	    for (it=0;it<inp->nit;it++) {
		dummy_();
		for (i=0;i<llen;i+=istr) {
		    a[i] = b[i]*a[i-1]+d[i];
		}
	    }
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t2 = MPI_Wtime();
	    break;
	case 8:
	    sprintf(res->label, "REAL_8_INT_8: CHARGE ASSIGNMENT: a(ia) = a(ia)+const");
	    res->nops[3] = 1;
	    res->nbytes[0] = 1 * 8 + 1 * 8;
	    res->nbytes[1] = 0 * 8;
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t1 = MPI_Wtime();
	    s = 1.7349;
#pragma omp parallel
	    for (it=0;it<inp->nit;it++) {
		dummy_();
#pragma omp for
		for (i=0;i<llen;i+=istr) {
		    c[ia[i]] += s;
		}
	    }
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t2 = MPI_Wtime();
	    break;
	case 9:
	    sprintf(res->label, "REAL_8: DAXPY: a=const*b+c");
	    res->nops[3] = 2;
	    res->nbytes[0] = 2 * 8;
	    res->nbytes[1] = 1 * 8;
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t1 = MPI_Wtime();
	    s = 1.7349;
#pragma omp parallel
	    for (it=0;it<inp->nit;it++) {
		dummy_();
#pragma omp for
		for (i=0;i<llen;i+=istr) {
		    a[i] = s * b[i] + c[i];
		}
	    }
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t2 = MPI_Wtime();
	    break;
	case 13:
	    sprintf(res->label, "INTEGER_8: TRIADS: ia=ib*ic+id");
	    res->nops[3] = 2;
	    res->nbytes[0] = 3 * 8;
	    res->nbytes[1] = 1 * 8;
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t1 = MPI_Wtime();
#pragma omp parallel
	    for (it=0;it<inp->nit;it++) {
		dummy_();
#pragma omp for
		for (i=0;i<llen;i+=istr) {
		    ia[i] = ib[i] * ic[i] + id[i];
		}
	    }
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t2 = MPI_Wtime();
	    break;
	case 38:
	    sprintf(res->label, "INT_8: weighted inverse Latency : Random read");
	    res->nops[3] = 1;
	    res->nbytes[0] = 1 * 8;
	    res->nbytes[1] = 0;
	    for (i=0;i<llen;i++) {
		k2 = (long long int) irandom_();
		if (k2 <= 0) { k2 = -k2; }
		ia[i] = k2%llen;
/*		printf(" i = %d  ia = %d\n",i,ia[i]); */
	    }
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t1 = MPI_Wtime();
	    for (it=0;it<inp->nit;it++) {
		dummy_();
		k1 = ia;
		for (i=1;i<llen;i++) {
		    k2 = ia[*k1];
		    *k1 = k2;
		}
	    }
#ifdef USE_MPI
	    mp_err = mpi_barrier(MPI_COMM_WORLD);
#endif
	    t2 = MPI_Wtime();
	    break;
	default:
	    t1 = 0.0;
	    t2 = -1.0;
	    sprintf(res->label," Case %d not implemented",inp->icase);
	    break;
    } /* switch (inp->icase) */
    res->time = (t2 - t1)/((double) inp->nit) - timings_.t0;
/*    printf(" nit is %d \n",inp->nit);
    printf(" t0 is  %e15.12 \n",timings_.t0);
    printf(" Time is %e t1, t2: %e %e \n",res->time,t1,t2); */
    res->perf = 0;
    for (i=0;i<6;i++) {
	res->nops[i] = res->nops[i] * ((inp->nvec - ndiff) / inp->istride);
	res->perf = res->perf + res->nops[i];
    }
    res->perf = res->perf / (1.0e6 * res->time);
    res->nbytes[0] = res->nbytes[0] * ((inp->nvec - ndiff) / inp->istride);
    res->nbytes[1] = res->nbytes[1] * ((inp->nvec - ndiff) / inp->istride);
    res->bw = (res->nbytes[0] + res->nbytes[1])/(1.0e6 * res->time);
/*    printf(" BW is %e \n",res->bw);  */
}
