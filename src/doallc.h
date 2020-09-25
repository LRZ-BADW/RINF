#include <stdio.h>

/* global variables and definitions*/


extern struct arraydata {
    double a[RINF1_MAXLEN];
    double b[RINF1_MAXLEN];
    double c[RINF1_MAXLEN];
    double d[RINF1_MAXLEN];
    double e[RINF1_MAXLEN];
    double f[RINF1_MAXLEN];
    long long int ia[RINF1_MAXLEN];
    long long int ib[RINF1_MAXLEN];
    long long int ic[RINF1_MAXLEN];
    long long int id[RINF1_MAXLEN];
} arrays_;
extern struct aux {
    long long int k1;
} b_;
extern struct timedata {
    double t0;
} timings_;


double t1,t2;

double *a,*b,*c,*d,*e,*f;
long long int *ia, *ib, *ic, *id;

typedef struct rinf1_input {
    int icase;
    int istride;
    int num_lengths;
    int nvec;
    double titer[3];
    int nit;
    int niter[3];
    int nitisone;
} Rinf1_input;

typedef struct rinf1_output {
    char label[70];
    int nops[6];
    int nbytes[2];
    double  time;
    double  bw;
    double bwdev;
    double perf;
    double perfdev;
    double testsum;
} Rinf1_output;


/* prototypes */

void doallc_(Rinf1_input *inp, Rinf1_output *res);
void rinf1_abort_(char *err);
double dwalltime00();
void dummy_();
int irandom_();
