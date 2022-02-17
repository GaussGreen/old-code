/*<NP+>*/

/*<I4+> default to INTEGER*4*/





double   sqrt(double);

double   atan2(double,double);

double   cabs(complexf);

complexf clog(complexf);



void     e1psh(char*);

void     e1pop(char*);

void     e1sti(int,int);

void     e1str(int,float);

void     e1std(int,double);

void     e1stc(int,complexf);

void     e1stz(int,complex);

void     e1stl(int,char*);

void     e1csm(int,double*);

void     e1usr(char*);

void     e1pos(int,int*,int*);

void     e1mes(int,int,char*);

void	 e1chk(double);

int      n1rty(int);

int      n1rcd(int);

int      n1rnof(int);



float    amach(int);

double   dmach(int);

void     umach(int,int*);

char     achar(int);

int      iachar(char*,unsigned);



	/* Chapter 1 */



void     page(int,int*);



void     l2nrg(int,float*,int,float*,int,float[],int[]);

void     lftrg(int,float*,int,float*,int,int[]);

void     dlftrg(int,double*,int,double*,int,int[]);

void     l2trg(int,float*,int,float*,int,int[],float[]);

void     dl2trg(int,double*,int,double*,int,int[],double[]);

void     lfsrg(int,float*,int,int[],float*,int,float*);

void     dlfsrg(int,double*,int,int[],double*,int,double*);



void     lslrg(int,float*,int,float*,int,float*);

void     dlslrg(int,double*,int,double*,int,double*);

void     lsarg(int,float*,int,float*,int,float*);

void     dlsarg(int,double*,int,double*,int,double*);

void     lslcg(int,complexf*,int,complexf[],int,complexf[]);

void     dlslcg(int,complex*,int,complex[],int,complex[]);



void     l2lrb(int*, float*, int*, int*, int*, float*, int*,float*, float*, int*,float*);

void     dl2lrb(int*, double*, int*, int*, int*, double*, int*,double*, double*, int*,double*);

void     l2crb(int*, float*, int*, int*, int*, float*,int*, int*, float*,float*);

void     dl2crb(int*, double*, int*, int*, int*, double*,int*, int*, double*,double*);

void     l2trb(int*, float*, int*, int*, int*, float*, int*, int*, float*);

void     dl2trb(int*, double*, int*, int*, int*, double*, int*, int*, double*);

void     lfsrb(int*, float*, int*, int*,int*, int*, float*, int*, float*);

void     dlfsrb(int*, double*, int*, int*,int*, int*, double*, int*, double*);





        /* Chapter 2 */



void     evcrg(int,float*,int,complexf[],complexf*,int);

void     devcrg(int,double*,int,complex[],complex*,int);





	/* Chapter 3 */

void     c1sor(int,float*,float*,float*,float*,int,int*);

void     dc1sor(int,double*,double*,double*,double*,int,int*);

void     c2int(int*,float*,float*,float*,float*,int*);

void     dc2int(int*,double*,double*,double*,double*,int*);

void     c2dec(int*,float*,float*,int*,float*,int*,float*,float*,float*,int*);

void     dc2dec(int*,double*,double*,int*,double,int*,double*,double*,double*,int*);

float    csval(float*,int*,float*,float*);

double   dcsval(double*,int*,double*,double*);

float    csder(int*,float*,int*,float*,float*);

double   dcsder(int*,double*,int*,double*,double*);

float    csitg(float*,float*,int*,float*,float*);

double   dcsitg(double*,double*,int*,double*,double*);

float    ppder(int,float,int,int,float*,float*);

double   dppder(int,double,int,int,double*,double*);

void     p3der(int,int,float*,float,int*);

void     dp3der(int,int,double*,double,int*);

float    ppitg(float,float,int,int,float*,float*);

double   dppitg(double,double,int,int,double*,double*);

float    b2der(int*, float*, int*, float*, int*, float*,float*, float*, float*);

double   db2der(int*, double*, int*, double*, int*, double*,double*, double*,double*);

float    b3der(int*, float*, int*, float*,int*, float*,float*, float*, float*);

double   db3der(int*, double*, int*, double*,int*, double*,double*, double*,double*);

void     b4der(float*, int*, float*, int*, int*);

void     db4der(double*, int*, double*, int*, int*);

void     b2int(int*, float*, float*, int*, float*, float*,float*, float*, float*,int*);

void     db2int(int*, double*, double*, int*, double*, double*,double*, double*, double*,int*);

void     b3int(int*, float*, int*);

void     db3int(int*, double*, int*);

void     b4int(float*, int*, float*, int*, float*,float*, float*);

void     db4int(double*, int*, double*, int*, double*,double*, double*);

void     b5int(int*, float*, float*, int*, float*, float*,float*, int*, float*,int*);

void     db5int(int*, double*, double*, int*, double*, double*,double*, int*, double*,int*);

float    b2itg(float*, float*, int*, float*, int*,float*, float*, float*, float*,float*);

double   db2itg(double*, double*, int*, double*, int*,double*, double*, double*, double*,double*);

float    b3itg(float*, float*, int*, float*, int*,float*, float*, float*, float*,float*);

double   db3itg(double*, double*, int*, double*, int*,double*, double*, double*, double*,double*);

float    b4itg(float*, int*, float*, int*,float*, float*, float*, int*);

double   db4itg(double*, int*, double*, int*,double*, double*, double*, int*);

void     b5itg(float*, int*, float*, int*, int*);

void     db5itg(double*, int*, double*, int*, int*);

void     b2nak(int*, float*, int*, float*, float*,int*);

void     db2nak(int*, double*, int*, double*, double*,int*);

void     c1not(char*, char*, int*, float*,int*, float*);

void     dc1not(char*, char*, int*, double*,int*, double*);

void     b2opk(int*, float*, int*, float*,int*, float*, int*);

void     db2opk(int*, double*, int*, double*,int*, double*, int*);

void     b3opk(int*, float*, int*, float*,float*, float*, float*,float*, float*, float*, float*, float*,int*, float*, int*);

void     db3opk(int*, double*, int*, double*,double*, double*, double*,double*, double*, double*, double*, double*,int*, double*, int*);

void     b4opk(int*, int*, float*, float*,float*, float*, float*, float*, float*, float*,int*);

void     db4opk(int*, int*, double*, double*,double*, double*, double*, double*, double*, double*,int*);

float    b22dr(int*, int*, float*, float*,int*, int*, float*,float*, int*, int*, float*, float*);

double   db22dr(int*, int*, double*, double*,int*, int*, double*,double*, int*, int*, double*, double*);

void     b32dr(char*, int*, int*);

void     db32dr(char*, int*, int*);

float    b22ig(float*, float*, float*, float*, int*, int*,float*, float*, int*,int*, float*, float*);

double   db22ig(double*, double*, double*, double*, int*, int*,double*, double*, int*,int*, double*, double*);

float    b32ig(float*, float*, float*, float*, int*, int*,float*, float*, int*,int*, float*, float*, float*, float*, float*, float*); 

double   db32ig(double*, double*, double*, double*, int*, int*,double*, double*, int*,int*, double*, double*, double*, double*, double*, double*); 

void     b22in(int*, float*, int*, float*, float*, int*,int*, int*, float*, float*, float*,float*, int*);

void     db22in(int*, double*, int*, double*, double*, int*,int*, int*, double*, double*, double*,double*, int*);

void     b32in(char*, int*, float*, int*);

void     db32in(char*, int*, double*, int*);

void     b42in(char*, int*, int*, float*, float*,int*, int*, float*, float*, float*,int*, float*, float*, float*, int*);

void     db42in(char*, int*, int*, double*, double*,int*, int*, double*, double*, double*,int*, double*, double*, double*, int*);

void     b2ls2(int*, float*, int*, float*,float*, int*, int*, int*,float*, float*, int*, int*,float*, float*, float*, float*); 

void     db2ls2(int*, double*, int*, double*,double*, int*, int*, int*,double*, double*, int*, int*,double*, double*, double*, double*); 

void     b3ls2(int*, float*, int*, float*,int*, float*, char*, int,int*);

void     db3ls2(int*, double*, int*, double*,int*, double*, char*, int,int*);

void     b4ls2(int*, float*, int*, float*,float*, int*, int*, int*, float*, float*, int*, int*,float*, float*, float*, float*,float*, float*, float*, float*);

void     db4ls2(int*, double*, int*, double*,double*, int*, int*, int*, double*, double*, int*, int*,double*, double*, double*, double*,double*, double*, double*, double*);

void     b5ls2(int*, float*, float*, float*,int*, float*, int*, float*,float*, float*, int*);

void     db5ls2(int*, double*, double*, double*,int*, double*, int*, double*,double*, double*, int*);

void     b6ls2(int*, int*, int*, float*,float*, float*);

void     db6ls2(int*, int*, int*, double*,double*, double*);

void     b7ls2(int*);

void     db7ls2(int*);

void     b2lsq(int*, float*, float*, float*, int*,float*, int*, float*, float*, float*,float*, float*, int*);

void     db2lsq(int*, double*, double*, double*, int*,double*, int*, double*, double*, double*,double*, double*, int*);

void     b3lsq(int*, int*, float*, float*, float *);

void     db3lsq(int*, int*, double*, double*, double *);

void     b4lsq(int*, float*, float*, float*, int*,float*, int*, float*, float*, float*);

void     db4lsq(int*, double*, double*, double*, int*,double*, int*, double*, double*, double*);

void     b5lsq(float*, int*, int*);

void     db5lsq(double*, int*, int*);

void     b6lsq(float*, int*, int*, float*);

void     db6lsq(double*, int*, int*, double*);

void     b2vls(int*, float*, float*, float*, int*, int*,float*, float*, float*, float*, int*, float*);

void     db2vls(int*, double*, double*, double*, int*, int*,double*, double*, double*, double*, int*, double*);

void     b3vls(int*, float*, float*, float*, int*, int*,float*, float*, float*, float*, float*, float*,float*, float*, float*, float*);

void     db3vls(int*, double*, double*, double*, int*, int*,double*, double*, double*, double*, double*, double*,double*, double*, double*, double*);

void     b4vls(int*, int*, float*, float*, float*, int*,float*, int*, float*, float*, float*, float*, float*, float*);

void     db4vls(int*, int*, double*, double*, double*, int*,double*, int*, double*, double*, double*, double*, double*, double*);

void     b5vls(float*, float*, int*, int*, float*, int*, float*);

void     db5vls(double*, double*, int*, int*, double*, int*, double*);

float    b6vls(int*, float*, float*, float*, int*, float*,int*, float*, float*, float*, float*);

double   db6vls(int*, double*, double*, double*, int*, double*,int*, double*, double*, double*, double*);

void     b7vls(int*, float*, float*, float*, float*, float*, float*, int*, float*);

void     db7vls(int*, double*, double*, double*, double*, double*, double*, int*, double*);

void     b8vls(float*, int*, int*, float*);

void     db8vls(double*, int*, int*, double*);

void     b2cpp(int*, float*, int*, float*, int*, float*,float*, float*);

void     db2cpp(int*, double*, int*, double*, int*, double*,double*, double*);

void     b3cpp(int*, float*, int*, float*, int*, float*,float*, float*, float*, float*, float*);

void     db3cpp(int*, double*, int*, double*, int*, double*,double*, double*, double*, double*, double*);

void     f2lsq(float (*f)(),   int*, int*, int*,float*,  float*,  int*,float*,  float*,  float*,  float*);

void     df2lsq(double (*f)(), int*, int*, int*,double*, double*, int*,double*, double*, double*, double*);

void     c2scv(int*, float[], float[], int*, float[], float*,float[], float*, float[], int[]);

void     dc2scv(int*, double[], double[], int*, double[], double*,double[], double*, double[], int[]);

void     c3scv(float[], float*, float[], float[], float*,int*, float[],float*, float*, float*); 

void     dc3scv(double[], double*, double[], double[], double*,int*, double[],double*, double*, double*); 

void     c4scv(float[], float*, float[], int*,float*, float*, float*, float*,float*, float[], float[], float*,float*, float*, float[], float[]);

void     dc4scv(double[], double*, double[], int*,double*, double*, double*, double*,double*, double[], double[], double*,double*, double*, double[], double[]);

void     c5scv(float[], float*, float[], int*, float*,float*, float[], float*, float[], float[]);

void     dc5scv(double[], double*, double[], int*, double*,double*, double[], double*, double[], double[]);

void     c2smh(int*, float[], float[], float[], float*, float[],float*, float[], int[]);

void     dc2smh(int*, double[], double[], double[], double*, double[],double*, double[], int[]);

void     c3smh(int*, float[], float[], float*, float[], float*,float[], float[], float[], float[], float[], float[], int[]);

void     dc3smh(int*, double[], double[], double*, double[], double*,double[], double[], double[], double[], double[], double[], int[]);

float    c4smh(int*, float*);

double   dc4smh(int*, double*);



                              

                              

	/* Chapter 4 */



void     qdng(float(*f)(),float,float,float,float,float*,float*);

void     dqdng(double(*f)(),double,double,double,double,double*,double*);

void     q3ng(float(*f)(),float,float,float,float,float*,float*,int*,int*);

void     dq3ng(float(*f)(),double,double,double,double,double*,double*,int*,int*);

void     q4ng(float*,float*,float*);

void     dq4ng(double*,double*,double*);



void     qdag(float(*f)(),float,float,float,float,int,float*,float*);

void     dqdag(double(*f)(),double,double,double,double,int,double*,double*);

void     q2ag(float(*f)(),float,float,float,float,int,float*,float*,int,int*,int*,float*,float*,float*,float*,int[]);

void     dq2ag(double(*f)(),double,double,double,double,int,double*,double*,int,int*,int*,double*,double*,double*,double*,int[]);

void     g2rcf(int*, float[], float[], int*, float[],float[], float[], float[]);

void     dg2rcf(int*, double[], double[], int*, double[],double[], double[], double[]);

float    g3rcf(float*, int*, float[], float[]);

double   dg3rcf(double*, int*, double[], double[]);

void     g4rcf(int*, float[], float[], float[], float[], float[]);

void     dg4rcf(int*, double[], double[], double[], double[], double[]);

void     g2rul(int*, int*, float*, float*, int*, float[], float[], float[], float[]);

void     dg2rul(int*, int*, double*, double*, int*, double[], double[], double[], double[]);

void     reccf(int*, int*, float*, float*, float[], float[]);

void     dreccf(int*, int*, double*, double*, double[], double[]);





	/* Chapter 5*/



void     ivprk(int*,int,float(*f)(),float*,float,float,float*,float*);

void     divprk(int*,int,double(*f)(),double*,double,double,double*,double*);

void     i2prk(int*,int,float(*f)(),float*,float,float,float*,float*,float(*f)(),float*);

void     di2prk(int*,int,double(*f)(),double*,double,double,double*,double*,double(*f)(),double*);

void     i3prk(int,float*,float*,float*,float*);

void     di3prk(int,double*,double*,double*,double*);



	/* Level 1 BLAS */



int      isamax(int,float*,int);

int      idamax(int,double*,int);

#ifndef INLINE_BLAS

void     sscal(int,float,float*,int);

void     dscal(int,double,double*,int);

void     saxpy(int,float,float*,int,float*,int);

void     daxpy(int,double,double*,int,double*,int);

void     sswap(int,float*,int,float*,int);

void     dswap(int,double*,int,double*,int);

void     iset(int,int,int*,int);

void     sset(int,float,float*,int);

void     dset(int,double,double*,int);

void     icopy(int,int*,int,int*,int);

void     scopy(int,float*,int,float*,int);

void     dcopy(int,double*,int,double*,int);

void     sadd(int,float,float*,int);

void     dadd(int,double,double*,int);

#endif

float	 sasum(int,float*,int);

double   dasum(int,double*,int);

float    snrm2(int,float*,int);

double   dnrm2(int,double*,int);

float    sdot(int,float*,int,float*,int);

double   ddot(int,double*,int,double*,int);



	/* Level 9 BLAS */



void     sger(int,int,float,float*,int,float*,int,float*,int);

void     dger(int,int,double,double*,int,double*,int,double*,int);

void     strsv(char*,char*,char*,int,float[],int,float[],int);

void     dtrsv(char*,char*,char*,int,double[],int,double[],int);



	/* Chapter 9 */



void     trnrr(int,int,float*,int,int,int,float*,int);

void     dtrnrr(int,int,double*,int,int,int,double*,int);

void     crgrg(int,float*,int,float*,int);

void     dcrgrg(int,double*,int,double*,int);



	/* Utilities*/



void     svrgp(int,float*,float*,int*);

void     dsvrgp(int,double*,double*,int*);

void     svrgn(int,float*,float*);

void     dsvrgn(int,double*,double*);

void     svign(int,int*,int*);



void     crbrb(int*, float*, int*, int*, int*,float*, int*, int*, int*);

void     dcrbrb(int*, double*, int*, int*, int*,double*, int*, int*, int*);

void     nr1rb(int*, float*, int*, int*, int*, float*);

void     dnr1rb(int*, double*, int*, int*, int*, double*);

void     stbsv(char*, char*, char*,int*, int*, float*, int*, float*, int*);

void     dstbsv(char*, char*, char*,int*, int*, double*, int*, double*, int*);





	/* SFUN/LIBRARY*/



int      inits(float*,int,float);

int      initds(double*,int,double);

float    csevl(float,float*,int);

double   dcsevl(double,double*,int);



float    erf(float);

double   derf(double);

float    erfc(float);

double   derfc(double);

float    erfi(float);

double   derfi(double);

float    erfci(float);

double   derfci(double);



float    beta(float,float);

double   dbeta(double,double);

float    betai(float,float,float);

double   dbetai(double,double,double);

float    gamma(float);

double   dgamma(double);

float    gami(float,float);

double   dgami(double,double);

float    albeta(float,float);

double   dlbeta(double,double);

float    alngam(float);

double   dlngam(double);

float    gamit(float,float);

double   dgamit(double,double);



float    bsj0(float);

double   dbsj0(double);

float    bsj1(float);

double   dbsj1(double);

float    bsi0(float);

double   dbsi0(double);

float    bsi1(float);

double   dbsi1(double);

float    bsy0(float);

double   dbsy0(double);

float    bsy1(float);

double   dbsy1(double);

float    bsk0(float);

double   dbsk0(double);

float    bsk1(float);

double   dbsk1(double);



float    alnrel(float);

float    bsi0e(float);

double   dbsi0e(double);

float    bsi1e(float);

double   dbsi1e(double);

float    bsk0e(float);

double   dbsk0e(double);

float    bsk1e(float);

double   dbsk1e(double);



float    anorin(float);

double   dnorin(double);

float    anordf(float);

double   dnordf(double);

float    chidf(float,float);

double   dchidf(double,double);

float    chiin(float,float);

double   dchiin(double,double);

float    fdf(float,float,float);

double   dfdf(double,double,double);

float    fin(float,float,float);

double   dfin(double,double,double);

float    tdf(float,float);

double   dtdf(double,double);

float    tin(float,float);

double   dtin(double,double);

float    gamdf(float,float);

double   dgamdf(double,double);

float    bindf(int,int,float);

double   dbindf(int,int,double);

float    hypdf(int,int,int,int);

double   dhypdf(int,int,int,int);

float    poidf(int,float);

double   dpoidf(int,double);

float    betin(float,float,float);

double   dbetin(double,double,double);

float    betdf(float,float,float);

double   dbetdf(double,double,double);



void     r9gaml(float*,float*);

float    r9lgmc(float);

float    r9lgic(float,float,float);

float    r9lgit(float,float,float);

float    r9gmit(float,float,float,float,float);

void     algams(float,float*,float*);

float    gamr(float);

double   d9lgmc(double);



double   dble(complex);

double   dimag(complex);

complexf cmplx(float,float);

complexf dcmplx(double,double);



float    carg(complexf);

double   zarg(complex);

complexf clngam(complexf);

complexf clnrel(complexf);

complexf c9lgmc(complexf);



	/* AT Tests */



void     ststn(int,int*,int,CHAR_UNSGN,float,int*,float,int*,int);

void	 rititm(int*,int*);

