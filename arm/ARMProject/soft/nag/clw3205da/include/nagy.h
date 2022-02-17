#ifndef NAGY
#define NAGY
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagy.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library y Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2171 (Feb 1998).
 */


#include <stdio.h>
#define TRUNCATE(str, pos)  if (strlen(str) > pos) str[pos] = '\0';

#ifdef NAG_PROTO
extern const char *global_resfile;
extern const char *global_datafile;
#else
extern char *global_resfile;
extern char *global_datafile;
#endif

#ifdef NAG_PROTO
extern void  y07aac(const char *outfile);
#else
extern void y07aac();
#endif

#ifdef NAG_PROTO
extern void  y07abc(const char *outfile);
#else
extern void y07abc();
#endif

#ifdef NAG_PROTO
extern void  y07acc(const char *outfile);
#else
extern void y07acc();
#endif

#ifdef NAG_PROTO
extern void  y07adc(FILE *fp_in);
#else
extern void y07adc();
#endif

#ifdef NAG_PROTO
extern void  y07aec(FILE *fp_in, const char *outfile);
#else
extern void y07aec();
#endif

#ifdef NAG_PROTO
extern void  y90aac( char NAG_HUGE *nag, Boolean NAG_HUGE *tnag);
#else
extern void y90aac();
#endif

#ifdef NAG_PROTO
extern void   y90abc(char NAG_HUGE *prog, Boolean NAG_HUGE *tprog);
#else
extern void  y90abc();
#endif

#ifdef NAG_PROTO
extern void   y90adc(char NAG_HUGE *nag, Integer NAG_HUGE *xnag);
#else
extern void  y90adc();
#endif

#ifdef NAG_PROTO
extern void   y90aec(char NAG_HUGE *prog, Integer NAG_HUGE *xprog);
#else
extern void  y90aec();
#endif

#ifdef NAG_PROTO
extern void  y90ccc(Integer type, Integer dtype, Integer m, Integer n, Complex NAG_HUGE *a, 
            Integer ia, Complex NAG_HUGE *d, Complex NAG_HUGE *diag, double cond, Complex scale,
            Complex NAG_HUGE *detman, Integer NAG_HUGE *detexp, Integer dist, Integer NAG_HUGE *seed,
            char pvrow, char pvcol, Integer NAG_HUGE *ipvrow, Integer NAG_HUGE *ipvcol, 
            Complex NAG_HUGE *u, Integer iu, Complex NAG_HUGE *v, Integer iv, Complex NAG_HUGE *work1, 
            Integer iwork1, Complex NAG_HUGE *work2, Integer iwork2);
#else
extern void y90ccc();
#endif

#ifdef NAG_PROTO
extern void  y90cdc(char side, char init, Integer m, Integer n, Complex NAG_HUGE *a, 
            Integer tda, Complex NAG_HUGE *x, Complex NAG_HUGE *d, Integer NAG_HUGE *seed);
#else
extern void y90cdc();
#endif

#ifdef NAG_PROTO
extern void  y90cfc(Integer ttype, char uplo, Integer m, Integer n, Complex NAG_HUGE *a,
            Integer ia, Complex NAG_HUGE *diag, Complex NAG_HUGE *odiag, Complex NAG_HUGE *detman, 
            Integer NAG_HUGE *detexp, Integer dist, Integer NAG_HUGE *seed);
#else
extern void y90cfc();
#endif

#ifdef NAG_PROTO
extern void  y90cgc(char det, Integer vtype, Integer n, Complex NAG_HUGE *v, Complex NAG_HUGE *vbound,
            double cond, Complex scale, Complex NAG_HUGE *detman, Integer NAG_HUGE *detexp, 
            Integer dist, Integer NAG_HUGE *seed);
#else
extern void y90cgc();
#endif

#ifdef NAG_PROTO
extern void  y90ctc(char *matrix, char *jobs, Integer n, Integer njord,
		   Integer jord[], Complex eig[], Complex a[],
		   Integer tda, double cond, Complex x[],
		   Integer tdx, Complex y[], Integer tdy,
		   Complex q[], Integer tdq, Complex vec[],
		   Complex wk1[], Integer tdwk1, Integer seed[]);
#else
extern void y90ctc();
#endif


#ifdef NAG_PROTO
extern void  y90dgc(char type, Integer nlose, Complex NAG_HUGE *cmat, Integer icmat, Integer m, 
            Integer n);
#else
extern void y90dgc();
#endif

#ifdef NAG_PROTO
extern void  y90dhc(char matrix, Integer m, Integer n, Complex NAG_HUGE *a, Integer nda,
            Complex NAG_HUGE *b, Integer ndb);
#else
extern void y90dhc();
#endif

#ifdef NAG_PROTO
extern void  y90djc(char transa, char transb, Integer m, Integer n, Integer k,
            Complex alpha, Complex NAG_HUGE *a, Integer lda, Complex NAG_HUGE *b, Integer ldb,
            Complex beta, Complex NAG_HUGE *c, Integer ldc);
#else
extern void y90djc();
#endif

#ifdef NAG_PROTO
extern void  y90dkc(char matrix, Complex NAG_HUGE *a, Integer nda, Integer n);
#else
extern void y90dkc();
#endif

#ifdef NAG_PROTO
extern void  y90dmc(Complex NAG_HUGE *vman, Integer NAG_HUGE *vexp, Integer scale);
#else
extern void y90dmc();
#endif

#ifdef NAG_PROTO
extern void  y90dnc(char matra, char matrb, Integer m, Integer n, Integer mn, 
            Complex alpha, Complex NAG_HUGE *a, Integer ia, Complex NAG_HUGE *b, Integer ib, 
            Complex beta, Complex NAG_HUGE *c, Integer ic);
#else
extern void y90dnc();
#endif

#ifdef NAG_PROTO
extern Complex  y90ebc(Integer idist, Integer NAG_HUGE *seed);
#else
extern Complex y90ebc();
#endif

#ifdef NAG_PROTO
extern void  y90hac(const char *outfile, Integer inag, char *nag, Boolean tnag, Boolean wnag, Integer xnag);
#else
     extern void y90hac();
#endif

#ifdef NAG_PROTO
extern void  y90hbc(const char *outfile, char *prog, Boolean *tprog, Integer xprog);
#else
     extern void y90hbc();
#endif

#ifdef NAG_PROTO
extern void  y90hbf(const char *outfile, char *prog,  Boolean tprog,
		   Integer xprog);
#else
extern void y90hbf();
#endif


#ifdef NAG_PROTO
extern void  y90hcc(const char *outfile, Integer itask, char *task, Boolean ttask,
	    Boolean wtask, Integer xtask, char **nag, Integer nnag);
#else
     extern void y90hcc();
#endif

#ifdef NAG_PROTO
extern void  y90hdc(const char *outfile, Integer itest, char tytest, Boolean test,
	    Boolean warn, char *ttest);
#else
     extern void y90hdc();
#endif

#ifdef NAG_PROTO
extern void  y90hdf(const char *outfile, Integer itest,  Boolean test, Boolean warn,
	    char *ttest);
#else
extern void y90hdf();
#endif

#ifdef NAG_PROTO
extern void  y90pac(const char *outfile, char NAG_HUGE *line);
#else
extern void y90pac();
#endif

#ifdef NAG_PROTO
extern void  y90pbc(const char *outfile, char matrix, char NAG_HUGE *line, Complex NAG_HUGE *a, Integer tda, Integer m,
            Integer n);
#else
extern void y90pbc();
#endif

#ifdef NAG_PROTO
extern void  y90pcc(const char *outfile, char matrix, char NAG_HUGE *line, char NAG_HUGE *mname1, char NAG_HUGE *mname2, Complex NAG_HUGE *a,
            Integer tda, Complex NAG_HUGE *b, Integer tdb, Integer m, Integer n);
#else
extern void y90pcc();
#endif

#ifdef NAG_PROTO
extern void  y90pdc(const char *outfile, char matrix, char NAG_HUGE *line, double NAG_HUGE *a, Integer tda, Integer m,
            Integer n);
#else
extern void y90pdc();
#endif

#ifdef NAG_PROTO
extern void  y90pdf(const char *outfile, char *matrix, char *line,  Integer m, Integer n,
            double a[],  Integer tda);
#else
extern void y90pdf();
#endif

#ifdef NAG_PROTO
extern void  y90pef(const char *outfile, char *matrix, char *line, char *mname1, char *mname2,
            Integer m, Integer n,  double a[],  Integer tda,
            double b[],  Integer tdb);
#else
extern void y90pef();
#endif


#ifdef NAG_PROTO
extern void  y90pec(const char *outfile, char matrix, char NAG_HUGE *line, char NAG_HUGE *mname1, char NAG_HUGE *mname2, double NAG_HUGE *a,
            Integer tda, double NAG_HUGE *b, Integer tdb, Integer m, Integer n);
#else
extern void y90pec();
#endif

#ifdef NAG_PROTO
extern void  y90pfc(const char *outfile, char NAG_HUGE *line, Complex cvar);
#else
extern void y90pfc();
#endif

#ifdef NAG_PROTO
extern void  y90pgc(const char *outfile, char NAG_HUGE *line, char NAG_HUGE *strvar);
#else
extern void y90pgc();
#endif

#ifdef NAG_PROTO
extern void  y90phc(const char *outfile, char NAG_HUGE *line, Integer ivar);
#else
extern void y90phc();
#endif

#ifdef NAG_PROTO
extern void  y90pjc(const char *outfile, char NAG_HUGE *line, Boolean bvar);
#else
extern void y90pjc();
#endif

#ifdef NAG_PROTO
extern void  y90pkc(const char *outfile, char NAG_HUGE *line, double dvar);
#else
extern void y90pkc();
#endif

#ifdef NAG_PROTO
extern void  y90plc(const char *outfile, char NAG_HUGE *line, Complex NAG_HUGE *vec, Integer nvec);
#else
extern void y90plc();
#endif

#ifdef NAG_PROTO
extern void  y90pmc(const char *outfile, char NAG_HUGE *line, char NAG_HUGE *vname1, char NAG_HUGE *vname2, Complex NAG_HUGE *vec1,
            Complex NAG_HUGE *vec2, Integer nvec);
#else
extern void y90pmc();
#endif

#ifdef NAG_PROTO
extern void  y90pnc(const char *outfile, char NAG_HUGE *line, Integer NAG_HUGE *vec, Integer nvec);
#else
extern void y90pnc();
#endif

#ifdef NAG_PROTO
extern void  y90pnf(const char *outfile, char *line,  Integer nveci,
		   Integer veci[], Integer iveci);
#else
extern void y90pnf();
#endif

#ifdef NAG_PROTO
extern void  y90ppc(const char *outfile, char NAG_HUGE *line, char NAG_HUGE *vname1, char NAG_HUGE *vname2, Integer NAG_HUGE *vec1,
            Integer NAG_HUGE *vec2, Integer nvec);
#else
extern void y90ppc();
#endif

#ifdef NAG_PROTO
extern void  y90ppf(const char *outfile, char *line, char *vname1, char *vname2,
            Integer nveci, Integer veci1[], Integer iveci1, Integer veci2[],
            Integer iveci2);
#else
extern void y90ppf();
#endif




#ifdef NAG_PROTO
extern void  y90prc(const char *outfile, char NAG_HUGE *line, double NAG_HUGE *vecr, Integer nvecr);
#else
extern void y90prc();
#endif

#ifdef NAG_PROTO
extern void  y90prf(const char *outfile, char *line,  Integer nvecr,  double vecr[],
            Integer ivecr);
#else
extern void y90prf();
#endif


#ifdef NAG_PROTO
extern void  y90psc(const char *outfile, char NAG_HUGE *line, char NAG_HUGE *vname1, char NAG_HUGE *vname2, double NAG_HUGE *vecr1,
            double NAG_HUGE *vecr2, Integer nvecr);
#else
extern void y90psc();
#endif

#ifdef NAG_PROTO
extern void  y90psf(const char *outfile, char *line, char *vname1, char *vname2,  Integer nvecr,
	    double vecr1[],  Integer ivecr1,  double vecr2[],
	    Integer ivecr2);
#else
     extern void y90psf();
#endif

#ifdef NAG_PROTO
extern void  y90ptc(const char *outfile, char NAG_HUGE *line, Boolean NAG_HUGE *vecl, Integer nvecl);
#else
extern void y90ptc();
#endif

#ifdef NAG_PROTO
extern void  y90puc(const char *outfile, char NAG_HUGE *line, char NAG_HUGE *vname1, char NAG_HUGE *vname2, Boolean NAG_HUGE *vecl1,
            Boolean NAG_HUGE *vecl2,
            Integer nvecl);
#else
extern void y90puc();
#endif

#ifdef NAG_PROTO
extern void  y90pvc(const char *outfile, char matrix, char NAG_HUGE *line,  Integer m, Integer n, Integer NAG_HUGE a[],
            Integer ia);
#else
extern void y90pvc();
#endif

#ifdef NAG_PROTO
extern void  y90pwc(const char *outfile, char matrix, char NAG_HUGE *line, char NAG_HUGE *mname1, char NAG_HUGE *mname2,
             Integer m, Integer n, Integer NAG_HUGE a[], Integer ia, Integer NAG_HUGE b[],
             Integer ib);
#else
extern void y90pwc();
#endif

#ifdef NAG_PROTO
extern void  y90rac(char matrix,  Integer dtype, Integer ttype, Integer m,
             Integer n, Integer kl, Integer ku,  double NAG_HUGE a[],
             Integer ia,  double NAG_HUGE d[], double NAG_HUGE diag[],
             double NAG_HUGE odiag[], double cond, double scale,
             double NAG_HUGE *detman,  Integer NAG_HUGE *detexp, Integer dist,
             Integer NAG_HUGE seed[]);
#else
extern void y90rac();
#endif

#ifdef NAG_PROTO
extern void  y90raf(char matrix,  Integer dtype, Integer ttype, Integer m,
             Integer n, Integer kl, Integer ku,  double NAG_HUGE a[],
             Integer ia,  double NAG_HUGE d[], double NAG_HUGE diag[],
             double NAG_HUGE odiag[], double cond, double scale,
             double NAG_HUGE *detman,  Integer NAG_HUGE *detexp, Integer dist,
             Integer NAG_HUGE seed[]);
#else
extern void y90raf();
#endif

#ifdef NAG_PROTO
extern void  y90rbc(Integer m, Integer n, Integer kl, Integer ku, Integer NAG_HUGE seed[],
             double NAG_HUGE d[], double NAG_HUGE a[],  Integer tda);
#else
extern void y90rbc();
#endif

#ifdef NAG_PROTO
extern void  y90rbf(Integer m, Integer n, Integer kl, Integer ku, Integer NAG_HUGE seed[],
             double NAG_HUGE d[], double NAG_HUGE a[],  Integer tda);
#else
extern void y90rbf();
#endif

#ifdef NAG_PROTO
extern void   y90rcc(Integer type, Integer dtype, Integer m, Integer n, double NAG_HUGE *a,
             Integer ia, double NAG_HUGE *d, double NAG_HUGE *diag, double cond, double scale,
             double NAG_HUGE *detman, Integer NAG_HUGE *detexp, Integer dist, Integer NAG_HUGE *seed,
             char pvrow, char pvcol, Integer NAG_HUGE *ipvrow, Integer NAG_HUGE *ipvcol, 
             double NAG_HUGE *u, Integer iu, double NAG_HUGE *v, Integer iv, double NAG_HUGE *work1,
             Integer iwork1, double NAG_HUGE *work2, Integer iwork2);
#else
extern void  y90rcc();
#endif

#ifdef NAG_PROTO
extern void  y90rdf(char *side, char *init, Integer m, Integer n,
		   double a[], Integer tda, double x[], double d[],
		   Integer seed[]);
#else
extern void y90rdf();
#endif

#ifdef NAG_PROTO
extern void  y90rdc( char side, char init, Integer m, Integer n, double NAG_HUGE *a,
            Integer lda, double NAG_HUGE *x, double NAG_HUGE *d, Integer NAG_HUGE *seed);
#else
extern void y90rdc();
#endif

#ifdef NAG_PROTO
extern void  y90rec( char uplo, Integer n, Integer kb, Integer NAG_HUGE *seed, double NAG_HUGE *d,
            double NAG_HUGE *a, Integer lda);
#else
extern void y90rec();
#endif

#ifdef NAG_PROTO
extern void  y90ref( char uplo, Integer n, Integer kb, Integer NAG_HUGE *seed, double NAG_HUGE *d,
            double NAG_HUGE *a, Integer lda);
#else
extern void y90ref();
#endif

#ifdef NAG_PROTO
extern void  y90rfc( Integer ttype,  char uplo,  Integer m, Integer n,  double NAG_HUGE *a,
            Integer ia,  double NAG_HUGE *diag,  double NAG_HUGE *odiag,  double NAG_HUGE *detman,
            Integer NAG_HUGE *detexp,  Integer dist,  Integer NAG_HUGE *seed);
#else
extern void y90rfc();
#endif

#ifdef NAG_PROTO
extern void  y90rgc( char det, Integer vtype, Integer n, double NAG_HUGE *v, double NAG_HUGE *vbound,
            double cond, double scale, double NAG_HUGE *detman, Integer NAG_HUGE *detexp, Integer dist,
            Integer NAG_HUGE *seed);
#else
extern void y90rgc();
#endif

#ifdef NAG_PROTO
extern void  y90rgf(char *det, Integer vtype, Integer n, double v[],
		   double vbound[], double cond, double scale,
		   double *detman, Integer *detexp, Integer dist,
		   Integer seed[]);
#else
extern void y90rgf();
#endif


#ifdef NAG_PROTO
extern void  y90rhc(Integer dtype, Integer type, Integer n, Integer kl, Integer NAG_HUGE *nrow,
            double NAG_HUGE *a, double NAG_HUGE *d, double NAG_HUGE *diag, double NAG_HUGE *odiag, double cond, 
            double scale, double NAG_HUGE *detman,  Integer NAG_HUGE *detexp, Integer dist,
            Integer NAG_HUGE *seed, Integer NAG_HUGE *irow, double NAG_HUGE *work1, double NAG_HUGE *work2,
            Integer iwork2);
#else
extern void y90rhc();
#endif

#ifdef NAG_PROTO
extern void  y90rjc(    Integer n, Integer k, Integer l, Integer dtype, double condv,
            double NAG_HUGE *diag, double NAG_HUGE *rrx, double NAG_HUGE *rix, Integer ncj, Integer NAG_HUGE *icj,
            double NAG_HUGE *a, Integer ia, double NAG_HUGE *vrx, Integer ivrx, double NAG_HUGE *vrix,
            Integer ivrix, Integer NAG_HUGE *seed, double NAG_HUGE *x, double NAG_HUGE *wk1, Integer iwk1,
            double NAG_HUGE *wk2, Integer iwk2, double NAG_HUGE *wk3, Integer iwk3);
#else
extern void y90rjc();
#endif


#ifdef NAG_PROTO
extern void  y90rtc(char *matrix, char *jobs, Integer n, Integer njord,
		   Integer jord[], double eigr[], double eigi[],
		   double a[], Integer tda, double cond, double x[],
		   Integer tdx, double yh[], Integer tdyh, double q[],
		   Integer tdq, double vec[], Integer ivec[],
		   double wk1[], Integer tdwk1, Integer seed[]);
#else
extern void y90rtc();
#endif

#ifdef NAG_PROTO
extern void  y90rtx(Integer n, Integer njord, Integer jord[], double eigr[],
		   double eigi[], double a[], Integer tda,
		   double x[], Integer tdx, double y[], Integer tdy,
		   double wk1[], Integer tdwk1);
#else
extern void y90rtx();
#endif


#ifdef NAG_PROTO
extern void  y90ruc(Integer dtype, Integer n,  double *dens,  Integer *nnz,
             double *a,  Integer *irow, Integer *icol, Integer idima,
             double *lambda, double *lambnd, double *cond,
             double scale,  Integer dist, Integer *seed, Integer *iwork,
             double *work);
#else
extern void y90ruc();
#endif


#ifdef NAG_PROTO
extern void  y90ruw(Integer n, Integer p, Integer *nnz,  double *a,
             Integer *irow, Integer *icol, Integer idima, Integer *istr,
             Integer *istc,  double *ap);
#else
 extern void y90ruw();
#endif


#ifdef NAG_PROTO
extern void  y90rux(Integer n, Integer p,  double *a,  Integer *irow,
             Integer *icol, Integer idima, Integer *istr, Integer *istc,
             double *ap);
#else
 extern void y90rux();
#endif

#ifdef NAG_PROTO
extern void  y90ruy(Integer n, Integer *nnz, Integer p, Integer q,  double c,
             double s, double *a,  Integer *irow, Integer *icol,
             Integer idima, Integer *istr, Integer *istc,  double *ap,
             double *aq);
#else
extern void y90ruy();
#endif

#ifdef NAG_PROTO
extern void  y90ruz(Integer n,  double dens,  Integer *nnz,  double *a,
             Integer *irow, Integer *icol, Integer idima, Integer *seed,
             double *work,  Integer *iwork);
#else
extern void y90ruz();
#endif



#ifdef NAG_PROTO
extern void  y90rvf(Integer dtype,  char *diag,  Integer n,  double *dens,
             Integer *nnz,  double a[],  Integer irow[], Integer icol[],
             Integer idima,  double lambda[], double lambnd[],
             double cond, double scale,  Integer dist, Integer seed[],
             Integer iwork[],  double work[]);
#else
extern void y90rvf();
#endif



#ifdef NAG_PROTO
extern void  y90rvw(Integer n,  double dens,  Integer *nnz,  double a[],
             Integer irow[], Integer icol[], Integer idima, Integer seed[],
             double work[],  Integer iwork[]);
#else
extern void y90rvw();
#endif


#ifdef NAG_PROTO
extern void  y90rvx(Integer n, Integer *nnz, Integer p, Integer q,  double c,
             double s,  Boolean pre,  double a[],  Integer irow[],
             Integer icol[], Integer idima, Integer istr[], Integer istc[],
             double ap[], double aq[]);
#else
extern void y90rvx();
#endif


#ifdef NAG_PROTO
extern void  y90rvy(Integer n, Integer p,  char *rc,  double a[],
             Integer irow[], Integer icol[], Integer idima, Integer istr[],
             Integer istc[],  double ap[]);
#else
extern void y90rvy();
#endif

#ifdef NAG_PROTO
extern void  y90rvz(Integer n, Integer p,  char *rc,  Integer *nnz,
             double a[],  Integer irow[], Integer icol[], Integer idima,
             Integer istr[], Integer istc[],  double ap[]);
#else
extern void y90rvz();
#endif


#ifdef NAG_PROTO
extern void  y90sbc( char conv, char matrix, Integer m, Integer n, Integer kl,
            Integer ku, double NAG_HUGE *a, Integer ia, double NAG_HUGE *b, Integer ib);
#else
extern void y90sbc();
#endif

#ifdef NAG_PROTO
extern void  y90scc( char matrix, Integer m, Integer mn, Integer n, Integer kl,
            Integer ku, double NAG_HUGE *a, Integer lda, double NAG_HUGE *b, Integer ldb, double NAG_HUGE *c,
            Integer ldc);
#else
extern void y90scc();
#endif

#ifdef NAG_PROTO
extern double  y90sdc( char norm, char matrix, Integer m, Integer n, Integer kl,
              Integer ku, double NAG_HUGE *a, Integer lda);
#else
extern double y90sdc();
#endif

#ifdef NAG_PROTO
extern void  y90sec( Integer m, Integer mn, Integer n, Integer kl, Integer ku,
            double NAG_HUGE *l, Integer il, double NAG_HUGE *u, Integer iu, double NAG_HUGE *a, Integer ia);
#else
extern void y90sec();
#endif

#ifdef NAG_PROTO
extern void  y90sfc( Integer m, Integer n, Integer kl, Integer ku, double NAG_HUGE *a,
            Integer ia, double NAG_HUGE *b, Integer ib);
#else
extern void y90sfc();
#endif

#ifdef NAG_PROTO
extern void  y90sgc( char type, Integer nlose, double NAG_HUGE *rmat, Integer irmat, Integer m,
            Integer n);
#else
extern void y90sgc();
#endif

#ifdef NAG_PROTO
extern void  y90shc( char matrix, Integer m, Integer n, double NAG_HUGE *a, Integer nda,
            double NAG_HUGE *b, Integer ndb);
#else
extern void y90shc();
#endif

#ifdef NAG_PROTO
extern void  y90shf(char *matrix, Integer m, Integer n, double a[],
		   Integer nda, double b[], Integer ndb);
#else
extern void y90shf();
#endif


#ifdef NAG_PROTO
extern void  y90sjc( char transa, char transb, Integer m, Integer n, Integer k,
            double alpha, double NAG_HUGE *a, Integer lda, double NAG_HUGE *b, Integer ldb, double beta,
            double NAG_HUGE *c, Integer ldc);
#else
extern void y90sjc();
#endif

#ifdef NAG_PROTO
extern void  y90skc( char matrix, double NAG_HUGE *a, Integer nda,
            Integer n);
#else
extern void y90skc();
#endif

#ifdef NAG_PROTO
extern void  y90smc( double NAG_HUGE *vman, Integer NAG_HUGE *vexp,
            Integer scale);
#else
extern void y90smc();
#endif

#ifdef NAG_PROTO
extern void  y90snc( char matra, char matrb, Integer m, Integer n, Integer mn,
            double alpha, double NAG_HUGE *a, Integer ia, double NAG_HUGE *b, Integer ib,
            double beta, double NAG_HUGE *c, Integer ic);
#else
extern void y90snc();
#endif

#ifdef NAG_PROTO
extern void  y90spc(Integer select, Integer m, Integer n, Integer NAG_HUGE *irow, double NAG_HUGE *l, 
            double NAG_HUGE *d, double NAG_HUGE *a, Integer ia, double NAG_HUGE *b, Integer ib, 
            double NAG_HUGE *work, Integer iwork);
#else
extern void y90spc();
#endif

#ifdef NAG_PROTO
extern void  y90sqc(Integer n, Integer NAG_HUGE *irow, double NAG_HUGE *l, double NAG_HUGE *d, double NAG_HUGE *a);
#else
extern void y90sqc();
#endif

#ifdef NAG_PROTO
extern double  y90src(char norm, Integer n, Integer NAG_HUGE *irow, double NAG_HUGE *a,double NAG_HUGE *work);
#else
extern double y90src();
#endif

#ifdef NAG_PROTO
extern double  y90tac( Integer NAG_HUGE *seed);
#else
extern double y90tac();
#endif

#ifdef NAG_PROTO
extern double  y90tbc( Integer dist, Integer NAG_HUGE *seed);
#else
extern double y90tbc();
#endif

#ifdef NAG_PROTO
extern void  y90tcc( char order, double NAG_HUGE *v, Integer n);
#else
extern void y90tcc();
#endif

#ifdef NAG_PROTO
extern void  y90tdc( Integer n, Integer k, Integer l, double NAG_HUGE *d, Integer NAG_HUGE *intger,
            double NAG_HUGE *a, Integer ia, double NAG_HUGE *w, Integer iw);
#else
extern void y90tdc();
#endif

#ifdef NAG_PROTO
extern void  y90tec( Integer n, Integer k, Integer l, double NAG_HUGE d[], Integer NAG_HUGE intger[],
            double NAG_HUGE a[], Integer ia, double NAG_HUGE w1[], Integer iw1, double NAG_HUGE w2[],
            Integer iw2);
#else
extern void y90tec();
#endif

#ifdef NAG_PROTO
extern void  y90tfc( Integer n, double NAG_HUGE *rr, double NAG_HUGE *ri, Integer NAG_HUGE *iord1, Integer NAG_HUGE *iord2);
#else
extern void y90tfc();
#endif

#ifdef NAG_PROTO
extern double  y90tgc(char norm, char matrix,  Integer m, Integer n,
              double NAG_HUGE a[],  Integer lda);
#else
extern double y90tgc();
#endif

#ifdef NAG_PROTO
extern Boolean  y90wac(char ca, char cb);
#else
extern Boolean y90wac();
#endif

#ifdef NAG_PROTO
extern void  y90wec(char NAG_HUGE *srname, Integer info);
#else
extern void y90wec();
#endif

#ifdef NAG_PROTO
extern void  y90zfc(char prtype, char matype, char select,  Integer NAG_HUGE *n,
             double NAG_HUGE *x);
#else
extern void y90zfc();
#endif

#ifdef NAG_PROTO
extern void  y99aac(char NAG_HUGE *msg, int code, char NAG_HUGE *name);
#else
extern void y99aac();
#endif

#ifdef NAG_PROTO
extern char *  NAG_HUGE y99abc(int errorcode, char NAG_HUGE errormessage[], NagError NAG_HUGE *fail);
#else
extern char *y99abc();
#endif

#ifdef __cplusplus
}
#endif
#endif /* not NAGY */
