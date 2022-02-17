#ifndef NAGF06
#define NAGF06
#ifdef __cplusplus
extern "C"
{
#endif


/* <nagf06.h>
 *
 * Copyright 1996 Numerical Algorithms Group
 *
 * Include file for NAG C Library f06 Chapter
 *
 * Mark 4, revised, 1996.
 * Mark 5 revised. IER-2153 (Feb 1998).
 */

#ifdef NAG_PROTO
extern void f06aac(double NAG_HUGE *a, double NAG_HUGE *b, double NAG_HUGE *c, double NAG_HUGE *s);
#else
extern void f06aac();
#endif

#ifdef NAG_PROTO
extern void f06aaz(char NAG_HUGE *srname, Integer info);
#else
extern void f06aaz();
#endif

#ifdef NAG_PROTO
extern void f06bac(double NAG_HUGE *a, double NAG_HUGE *b, double NAG_HUGE *c, double NAG_HUGE *s);
#else
extern void f06bac();
#endif

#ifdef NAG_PROTO
extern void f06bcc(double t, double NAG_HUGE *c, double NAG_HUGE *s);
#else
extern void f06bcc();
#endif

#ifdef NAG_PROTO
extern double f06blc(double a, double b, Boolean NAG_HUGE *failed);
#else
extern double f06blc();
#endif

#ifdef NAG_PROTO
extern double f06bmc(double scale, double ssq);
#else
extern double f06bmc();
#endif

#ifdef NAG_PROTO
extern double f06bnc(double a, double b);
#else
extern double f06bnc();
#endif

#ifdef NAG_PROTO
extern void f06cac(Complex NAG_HUGE *a, Complex NAG_HUGE *b, double NAG_HUGE *c, Complex NAG_HUGE *s);
#else
extern void f06cac();
#endif

#ifdef NAG_PROTO
extern void f06ccc(Complex t, double NAG_HUGE *c, Complex NAG_HUGE *s);
#else
extern void f06ccc();
#endif

#ifdef NAG_PROTO
extern Complex f06clc(Complex a, Complex b, Boolean NAG_HUGE *fail);
#else
extern Complex f06clc();
#endif

#ifdef NAG_PROTO
extern void f06fcz(Integer n,  double d[],  Integer incd,  double x[],
             Integer incx);
#else
extern void f06fcz();
#endif

#ifdef NAG_PROTO
extern void f06dbc(Integer n, Integer alpha, Integer NAG_HUGE x[], Integer incx);
#else
extern void f06dbc();
#endif

#ifdef NAG_PROTO
extern void f06dfc(Integer n, Integer NAG_HUGE *x, Integer incx, Integer NAG_HUGE *y, Integer incy);
#else
extern void f06dfc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP double NAG_CALL f06eac(Integer n, double NAG_HUGE *x, Integer incx, double NAG_HUGE *y, Integer incy);
#else
extern double f06eac();
#endif

#ifdef NAG_PROTO
extern void f06ecc(Integer n, double alpha, double NAG_HUGE *x, Integer incx, double NAG_HUGE *y,
            Integer incy);
#else
extern void f06ecc();
#endif

#ifdef NAG_PROTO
extern void f06edc(Integer n, double alpha, double NAG_HUGE *x, Integer incx);
#else
extern void f06edc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06efc(Integer n, double NAG_HUGE *x, Integer incx, double NAG_HUGE *y, Integer incy);
#else
extern void f06efc();
#endif

#ifdef NAG_PROTO
extern void f06egc(Integer n,double NAG_HUGE *x, Integer incx, double NAG_HUGE *y, Integer incy);
#else
extern void f06egc();
#endif

#ifdef NAG_PROTO
extern double f06ejc(Integer n,double NAG_HUGE *x, Integer incx);
#else
extern double f06ejc();
#endif

#ifdef NAG_PROTO
extern double f06ekc(Integer n, double NAG_HUGE *x, Integer incx);
#else
extern double f06ekc();
#endif

#ifdef NAG_PROTO
extern void f06epc(Integer n, double NAG_HUGE *x, Integer incx, double NAG_HUGE *y, Integer incy,
            double c, double s);
#else
extern void f06epc();
#endif

#ifdef NAG_PROTO
extern void f06fbc(Integer n, double cnst, double NAG_HUGE *x, Integer incx);
#else
extern void f06fbc();
#endif

#ifdef NAG_PROTO
extern void f06fpc(Integer n,  double x[],  Integer incx,  double y[],
             Integer incy,  double c, double s);
#else
extern void f06fpc();
#endif

#ifdef NAG_PROTO
extern void f06fcc(Integer n, double NAG_HUGE *d, Integer incd, double NAG_HUGE *x, Integer incx);
#else
extern void f06fcc();
#endif

#ifdef NAG_PROTO
extern void f06fdc(Integer n, double alpha, double NAG_HUGE *x, Integer incx,
            double NAG_HUGE *y, Integer incy);
#else
extern void f06fdc();
#endif

#ifdef NAG_PROTO
extern void f06fgc(Integer n, double NAG_HUGE *x, Integer incx);
#else
extern void f06fgc();
#endif

#ifdef NAG_PROTO
extern void f06fjc(Integer n, double NAG_HUGE *x, Integer incx, double NAG_HUGE *scale, double NAG_HUGE *sumsq);
#else
extern void f06fjc();
#endif

#ifdef NAG_PROTO
extern void f06flc(Integer n, double NAG_HUGE x[], Integer incx, double NAG_HUGE *xmax,
            double NAG_HUGE *xmin);
#else
extern void f06flc();
#endif

#ifdef NAG_PROTO
extern void f06fqc(PivotType pivot, SequenceDirection direct, Integer n,
            double NAG_HUGE *alpha, double NAG_HUGE *x, Integer incx, double NAG_HUGE *c, double NAG_HUGE *s);
#else
extern void f06fqc();
#endif

#ifdef NAG_PROTO
extern void f06frc(Integer n, double NAG_HUGE *alpha, double NAG_HUGE *x, Integer incx,
            double tol, double NAG_HUGE *zeta);
#else
extern void f06frc();
#endif

#ifdef NAG_PROTO
extern void f06fsc(Integer n, double NAG_HUGE *alpha, double NAG_HUGE *x, Integer incx, double tol,
            double NAG_HUGE *z1);
#else
extern void f06fsc();
#endif

#ifdef NAG_PROTO
extern void f06ftc(Integer n, double NAG_HUGE *delta, double NAG_HUGE *y, Integer incy,
            double zeta, double NAG_HUGE *z, Integer incz);
#else
extern void f06ftc();
#endif

#ifdef NAG_PROTO
extern void f06fuc(Integer n, double NAG_HUGE *z, Integer incz, double z1, double NAG_HUGE *alpha,
            double NAG_HUGE *x, Integer incx);
#else
extern void f06fuc();
#endif

#ifdef NAG_PROTO
/* Double Complex */ int f06gaf(Complex *ret_val, Integer n, Complex *x,
Integer incx, Complex *y, Integer incy);
#else
extern /* Double Complex */ int f06gaf();
#endif

#ifdef NAG_PROTO
extern Complex f06gbc(Integer n, Complex NAG_HUGE *x, Integer incx, Complex NAG_HUGE *y, Integer incy);
#else
extern Complex f06gbc();
#endif

#ifdef NAG_PROTO
extern void f06gcc(Integer n, Complex alpha, Complex NAG_HUGE *x, Integer incx, Complex NAG_HUGE *y,
            Integer incy);
#else
extern void f06gcc();
#endif

#ifdef NAG_PROTO
extern void f06gdc(Integer n, Complex alpha, Complex NAG_HUGE *x, Integer incx);
#else
extern void f06gdc();
#endif

#ifdef NAG_PROTO
extern void f06gfc(Integer n, Complex NAG_HUGE *x, Integer incx, Complex NAG_HUGE *y, Integer incy);
#else
extern void f06gfc();
#endif

#ifdef NAG_PROTO
extern void f06ggc(Integer n, Complex NAG_HUGE *x, Integer incx, Complex NAG_HUGE *y, Integer incy);
#else
extern void f06ggc();
#endif

#ifdef NAG_PROTO
extern void f06hbc(Integer n, Complex alpha, Complex NAG_HUGE *x, Integer incx);
#else
extern void f06hbc();
#endif

#ifdef NAG_PROTO
extern void f06hgc(Integer n, Complex NAG_HUGE *x, Integer incx);
#else
extern void f06hgc();
#endif

#ifdef NAG_PROTO
extern void f06hpc(Integer n, Complex NAG_HUGE *x, Integer incx, Complex NAG_HUGE *y, Integer incy, Complex c,
            Complex s);
#else
extern void f06hpc();
#endif

#ifdef NAG_PROTO
extern void f06hqc(PivotType pivot, SequenceDirection direct, Integer n, Complex NAG_HUGE *alpha,
            Complex NAG_HUGE *x, Integer incx, double NAG_HUGE *c, Complex NAG_HUGE *s);
#else
extern void f06hqc();
#endif

#ifdef NAG_PROTO
extern void f06hrc(Integer n, Complex NAG_HUGE *alpha, Complex NAG_HUGE *x, Integer incx, double tol,
            Complex NAG_HUGE *theta);
#else
extern void f06hrc();
#endif

#ifdef NAG_PROTO
extern void f06jdc(Integer n, double alpha, Complex NAG_HUGE *x, Integer incx);
#else
extern void f06jdc();
#endif

#ifdef NAG_PROTO
extern double f06jjc(Integer n, Complex NAG_HUGE *x, Integer incx);
#else
extern double f06jjc();
#endif

#ifdef NAG_PROTO
extern double f06jkc(Integer n, Complex NAG_HUGE *x, Integer incx);
#else
extern double f06jkc();
#endif

#ifdef NAG_PROTO
extern Integer f06jlc(Integer n, double NAG_HUGE *x, Integer incx);
#else
extern Integer f06jlc();
#endif

#ifdef NAG_PROTO
extern Integer f06jmc(Integer n, Complex NAG_HUGE *x, Integer incx);
#else
extern Integer f06jmc();
#endif

#ifdef NAG_PROTO
extern void f06kfc(Integer n, double NAG_HUGE *x, Integer incx, Complex NAG_HUGE *y, Integer incy);
#else
extern void f06kfc();
#endif

#ifdef NAG_PROTO
extern void f06kjc(Integer n, Complex NAG_HUGE *x, Integer incx, double NAG_HUGE *scale, double NAG_HUGE *sumsq);
#else
extern void f06kjc();
#endif

#ifdef NAG_PROTO
extern Integer f06klc(Integer n, double NAG_HUGE *x, Integer incx, double tol);
#else
extern Integer f06klc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06pac(MatrixTranspose Trans, Integer m, Integer n, 
            double alpha, const double NAG_HUGE a[], Integer tda, const double NAG_HUGE x[],
            Integer incx, double beta, double NAG_HUGE y[], Integer incy);
#else
extern void f06pac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06pbc(MatrixTranspose trans, Integer m, Integer n, Integer kl, 
             Integer ku, double alpha, const double NAG_HUGE a[], Integer tda,
             const double NAG_HUGE x[], Integer incx, double beta, double NAG_HUGE y[], 
             Integer incy);
#else
extern void f06pbc();
#endif

#ifdef NAG_PROTO
extern void f06pcf(char *uplo, Integer n, double alpha, double a[],
Integer lda, double x[], Integer incx, double beta,
double y[], Integer incy);
#else
extern void f06pcf();
#endif




#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06pcc(MatrixTriangle UpperLower, Integer n, double alpha, 
            const double NAG_HUGE a[], Integer tda, const double NAG_HUGE x[], 
            Integer incx, double beta, double NAG_HUGE y[], Integer incy);
#else
extern void f06pcc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06pdc(MatrixTriangle uplo, Integer n, Integer k, double alpha, 
            const double NAG_HUGE a[], Integer tda, const double NAG_HUGE x[], Integer incx,
            double beta, double NAG_HUGE y[], Integer incy);
#else
extern void f06pdc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06pec(MatrixTriangle UpperLower, Integer n, double alpha, 
            const double NAG_HUGE ap[], const double NAG_HUGE x[], Integer incx, 
            double beta, double NAG_HUGE y[], Integer incy);
#else
extern void f06pec();
#endif

#ifdef NAG_PROTO
extern void f06pef(MatrixTriangle UpperLower, Integer n, double alpha, double NAG_HUGE *ap,
            double NAG_HUGE *x, Integer incx, 
            double beta, double NAG_HUGE *y, Integer incy);
     
#else
     
extern void f06pef();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06pfc(MatrixTriangle UpperLower, MatrixTranspose Trans, 
            MatrixUnitTriangular TriangularMatrix, Integer n,
            const double NAG_HUGE a[], Integer tda, double NAG_HUGE x[], Integer incx);
#else
extern void f06pfc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06pgc(MatrixTriangle uplo, MatrixTranspose trans, 
            MatrixUnitTriangular diag, Integer n, Integer k, 
            const double NAG_HUGE a[], Integer tda, double NAG_HUGE x[], Integer incx);
#else
extern void f06pgc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06phc(MatrixTriangle UpperLower, MatrixTranspose Trans,
            MatrixUnitTriangular Diag, Integer n, const double NAG_HUGE ap[], 
            double NAG_HUGE x[], Integer incx);
#else
extern void f06phc();
#endif

#ifdef NAG_PROTO
extern void f06phf(MatrixTriangle UpperLower, MatrixTranspose Trans,
            MatrixUnitTriangular Diag, Integer n, double NAG_HUGE *ap, double NAG_HUGE *x,
            Integer incx);
#else
extern void f06phf();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06pjc(MatrixTriangle UpperLower, MatrixTranspose Trans, 
             MatrixUnitTriangular TriangularMatrix, Integer n, const double NAG_HUGE a[],
             Integer tda, double NAG_HUGE x[], Integer incx);
#else
extern void f06pjc();
#endif

#ifdef NAG_PROTO

extern NAG_DLL_EXPIMP void NAG_CALL f06pkc(MatrixTriangle uplo, MatrixTranspose trans,
            MatrixUnitTriangular diag, Integer n, Integer k,
            const double NAG_HUGE a[], Integer tda, double NAG_HUGE x[], Integer incx);
#else
extern void f06pkc();
#endif

#ifdef NAG_PROTO

extern void f06pkf(MatrixTriangle uplo, MatrixTranspose trans,
            MatrixUnitTriangular diag, Integer n, Integer k,
            double NAG_HUGE *a, Integer lda, double NAG_HUGE *x, Integer incx);
#else
extern void f06pkf();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06plc(MatrixTriangle UpperLower, MatrixTranspose Trans,
            MatrixUnitTriangular Diag, Integer n, const double NAG_HUGE ap[], 
            double NAG_HUGE x[], Integer incx);
#else
extern void f06plc();
#endif

#ifdef NAG_PROTO
extern void f06plf(MatrixTriangle UpperLower, MatrixTranspose Trans,
            MatrixUnitTriangular Diag, Integer n, double NAG_HUGE *ap, double NAG_HUGE *x,
            Integer incx);
#else
extern void f06plf();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06pmc(Integer m, Integer n, double alpha, const double NAG_HUGE x[],
            Integer incx, const double NAG_HUGE y[], Integer incy, 
            double NAG_HUGE a[], Integer tda);
#else
extern void f06pmc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06ppc(MatrixTriangle uplo, Integer n, double alpha, const double NAG_HUGE x[],
            Integer incx, double NAG_HUGE a[], Integer tda);
#else
extern void f06ppc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06pqc(MatrixTriangle uplo, Integer n, double alpha,
            const double NAG_HUGE x[], Integer incx, double NAG_HUGE ap[]);
#else
extern void f06pqc();
#endif

#ifdef NAG_PROTO
extern void f06pqf(MatrixTriangle uplo, Integer n, double alpha,
            double NAG_HUGE x[], Integer incx, double NAG_HUGE ap[]);
#else
extern void f06pqf();
#endif

#ifdef NAG_PROTO
extern void f06prf(char *uplo, Integer n, double alpha, double x[],
Integer incx, double y[], Integer incy, double a[],
Integer lda);
#else
extern void f06prf();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06prc(MatrixTriangle UpperLower, Integer n, double alpha, 
            const double NAG_HUGE x[], Integer incx, const double NAG_HUGE y[], Integer incy,
            double NAG_HUGE a[], Integer tda);
#else
extern void f06prc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06psc(MatrixTriangle UpperLower, Integer n, double alpha, 
            const double NAG_HUGE x[], Integer incx, const double NAG_HUGE y[], 
            Integer incy, double NAG_HUGE ap[]);
#else
extern void f06psc();
#endif

#ifdef NAG_PROTO
extern void f06psf(MatrixTriangle UpperLower, Integer n,  double alpha, double x[],
	    Integer incx,  double y[],  Integer incy,
	    double ap[]);
#else
     extern void f06psf();
#endif

#ifdef NAG_PROTO
extern void f06qfc(MatrixType matrix, Integer m, Integer n, double NAG_HUGE *a, Integer tda,
            double NAG_HUGE *b, Integer tdb);
#else
extern void f06qfc();
#endif

#ifdef NAG_PROTO
extern double f06qgc(NormType norm, MatrixType matrix, Integer m, Integer n,
              double NAG_HUGE *a, Integer tda);
#else
extern double f06qgc();
#endif

#ifdef NAG_PROTO
extern void f06qhc(MatrixType matrix, Integer m, Integer n, double constant,
            double diag, double NAG_HUGE *a, Integer tda);
#else
extern void f06qhc();
#endif

#ifdef NAG_PROTO
extern void f06qjc(OperationSide side, MatrixTranspose trans, Integer n, Integer NAG_HUGE *perm, Integer k,
            double NAG_HUGE *b, Integer ldb);
#else
extern void f06qjc();
#endif

#ifdef NAG_PROTO
extern void f06qkc(OperationSide Side, MatrixTranspose Trans, Integer n, 
            double NAG_HUGE *perm, Integer k, double NAG_HUGE *b, Integer tdb);
#else
extern void f06qkc();
#endif

#ifdef NAG_PROTO
extern void f06qkf(OperationSide side, MatrixTranspose trans, Integer n,
	    double NAG_HUGE perm[], Integer k,  double NAG_HUGE b[],  Integer ldb);
#else
extern void f06qkf();
#endif

#ifdef NAG_PROTO
extern void f06qnc(OperationSide side, Integer n, Integer k1, Integer k2,
            double NAG_HUGE s[], double NAG_HUGE a[], Integer tda);
#else
extern void f06qnc();
#endif

#ifdef NAG_PROTO
extern void f06qnz(char NAG_HUGE *side,  Integer n, Integer k1, Integer k2,
            double NAG_HUGE s[], double NAG_HUGE a[],  Integer lda);
#else
extern void f06qnz();
#endif

#ifdef NAG_PROTO
extern void f06qqc(Integer n, double alpha, double NAG_HUGE *x, Integer incx,
            double NAG_HUGE *a, Integer tda, double NAG_HUGE *c, double NAG_HUGE *s);
#else
extern void f06qqc();
#endif

#ifdef NAG_PROTO
extern void f06qrc(OperationSide side, Integer n, Integer k1, Integer k2, double NAG_HUGE *c,
            double NAG_HUGE *s, double NAG_HUGE *a, Integer tda);
#else
extern void f06qrc();
#endif

#ifdef NAG_PROTO
extern void f06qsc(OperationSide side,  Integer n, Integer k1, Integer k2,
            double NAG_HUGE c[], double NAG_HUGE s[], double NAG_HUGE a[],  Integer tda);
#else
extern void f06qsc();
#endif

#ifdef NAG_PROTO
extern void f06qtc(OperationSide side, Integer n, Integer k1, Integer k2, double NAG_HUGE *c,
            double NAG_HUGE *s, double NAG_HUGE *a, Integer tda);
#else
extern void f06qtc();
#endif

#ifdef NAG_PROTO
extern void f06qvc(OperationSide side, Integer n, Integer k1, Integer k2,
            double NAG_HUGE c[], double NAG_HUGE s[], double NAG_HUGE a[], Integer tda);
#else
extern void f06qvc();
#endif

#ifdef NAG_PROTO
extern void f06qwc(OperationSide side,  Integer n, Integer k1, Integer k2,
            double NAG_HUGE c[], double NAG_HUGE s[], double NAG_HUGE a[], Integer lda);
#else
extern void f06qwc();
#endif

#ifdef NAG_PROTO
extern void f06qxc(OperationSide side, PivotType pivot, SequenceDirection direct,
            Integer m, Integer  n, Integer k1, Integer k2, double NAG_HUGE *c,
            double NAG_HUGE *s, double NAG_HUGE *a, Integer tda);
#else
extern void f06qxc();
#endif

#ifdef NAG_PROTO
extern void f06qxf(OperationSide side, PivotType pivot, SequenceDirection direct,
	    Integer m, Integer n, Integer k1, Integer k2, double NAG_HUGE c[],
	    double NAG_HUGE s[], double NAG_HUGE a[], Integer lda);
#else
extern void f06qxf();
#endif

#ifdef NAG_PROTO
extern void f06qzz(char NAG_HUGE *hess,  Integer n, Integer k1, Integer k2,
            double NAG_HUGE c[], double NAG_HUGE s[], double NAG_HUGE a[],  Integer tda);
#else
extern void f06qzz();
#endif

#ifdef NAG_PROTO
extern double f06rac(NormType norm,  Integer m, Integer n,  double a[],
Integer tda,  double work[]);
#else
double f06rac();
#endif

#ifdef NAG_PROTO
extern double f06rcf(char *norm, char *uplo, Integer n, double a[],
Integer lda, double work[]);
#else
extern double f06rcf();
#endif

#ifdef NAG_PROTO
extern double f06rjc(char NAG_HUGE *norm, char NAG_HUGE *uplo, char NAG_HUGE *diag,  Integer m, Integer n,
              double NAG_HUGE a[],  Integer tda,  double NAG_HUGE work[]);
#else
extern double f06rjc();
#endif

#ifdef NAG_PROTO
extern double f06rmf(char *norm, Integer n, double a[], Integer lda,
double work[]);
#else
extern double f06rmf();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06sac(MatrixTranspose Trans, Integer m, Integer n, 
            Complex alpha, const Complex NAG_HUGE a[], Integer tda, const Complex NAG_HUGE x[],
            Integer incx, Complex beta, Complex NAG_HUGE y[], Integer incy);
#else
extern void f06sac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06sbc(MatrixTranspose trans, Integer m, Integer n, Integer kl, 
             Integer ku, Complex alpha, const Complex NAG_HUGE a[], Integer tda,
             const Complex NAG_HUGE x[], Integer incx, Complex beta, 
             Complex NAG_HUGE y[], Integer incy);
#else
extern void f06sbc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06scc(MatrixTriangle UpperLower, Integer n, Complex alpha,
            const Complex NAG_HUGE a[], Integer tda, const Complex NAG_HUGE x[], 
            Integer incx, Complex beta, Complex NAG_HUGE y[], Integer incy);
#else
extern void f06scc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06sdc(MatrixTriangle uplo, Integer n, Integer k, Complex alpha, 
            const Complex NAG_HUGE a[], Integer tda, const Complex NAG_HUGE x[], Integer incx, 
            Complex beta, Complex NAG_HUGE y[], Integer incy);
#else
extern void f06sdc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06sec(MatrixTriangle UpperLower, Integer n, Complex alpha, 
            const Complex NAG_HUGE ap[], const Complex NAG_HUGE x[], Integer incx, 
            Complex beta, Complex NAG_HUGE y[], Integer incy);
#else
extern void f06sec();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06sfc(MatrixTriangle UpperLower, MatrixTranspose Trans, 
             MatrixUnitTriangular TriangularMatrix, Integer n, 
             const Complex NAG_HUGE a[], Integer tda, Complex NAG_HUGE x[], Integer incx);
#else
extern void f06sfc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06sgc(MatrixTriangle uplo, MatrixTranspose trans, 
            MatrixUnitTriangular diag, Integer n, Integer k, 
            const Complex NAG_HUGE a[], Integer tda, Complex NAG_HUGE x[], Integer incx);
#else
extern void f06sgc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06shc(MatrixTriangle UpperLower, MatrixTranspose Trans,
            MatrixUnitTriangular Diag, Integer n, const Complex NAG_HUGE ap[],
            Complex NAG_HUGE x[], Integer incx);
#else
extern void f06shc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06sjc(MatrixTriangle UpperLower, MatrixTranspose Trans, 
            MatrixUnitTriangular TriangularMatrix, Integer n, 
            const Complex NAG_HUGE a[], Integer tda, Complex NAG_HUGE x[], Integer incx);
#else
extern void f06sjc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06skc(MatrixTriangle uplo, MatrixTranspose trans,
            MatrixUnitTriangular diag, Integer n, Integer k,
            const Complex NAG_HUGE a[], Integer tda, Complex NAG_HUGE x[], Integer incx);
#else
extern void f06skc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06slc(MatrixTriangle UpperLower, MatrixTranspose Trans,
            MatrixUnitTriangular Diag, Integer n, const Complex NAG_HUGE ap[],
            Complex NAG_HUGE x[], Integer incx);
#else
extern void f06slc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06smc(Integer m, Integer n, Complex alpha, const Complex NAG_HUGE x[], 
            Integer incx, const Complex NAG_HUGE y[], Integer incy, 
            Complex NAG_HUGE a[], Integer tda);
#else
extern void f06smc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06snc(Integer m, Integer n, Complex alpha, const Complex NAG_HUGE x[], 
            Integer incx, const Complex NAG_HUGE y[], Integer incy, Complex NAG_HUGE a[],
            Integer tda);
#else
extern void f06snc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06spc(MatrixTriangle uplo, Integer n, double alpha, 
            const Complex NAG_HUGE x[], Integer incx, Complex NAG_HUGE a[], Integer tda);
#else
extern void f06spc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06sqc(MatrixTriangle uplo, Integer n, double alpha,
            const Complex NAG_HUGE x[], Integer incx, Complex NAG_HUGE ap[]);
#else
extern void f06sqc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06src(MatrixTriangle UpperLower, Integer n, Complex alpha, 
            const Complex NAG_HUGE x[], Integer incx, const Complex NAG_HUGE y[],
            Integer incy, Complex NAG_HUGE a[], Integer tda);
#else
extern void f06src();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06ssc(MatrixTriangle UpperLower, Integer n, Complex alpha, 
            const Complex NAG_HUGE x[], Integer incx, const Complex NAG_HUGE y[], 
            Integer incy, Complex NAG_HUGE ap[]);
#else
extern void f06ssc();
#endif

#ifdef NAG_PROTO
extern void f06tfc(MatrixType matrix, Integer m, Integer n, Complex NAG_HUGE *a, Integer tda,
            Complex NAG_HUGE *b, Integer tdb);
#else
extern void f06tfc();
#endif

#ifdef NAG_PROTO
extern void f06thc(MatrixType matrix, Integer m, Integer n, Complex constant, 
            Complex diag, Complex NAG_HUGE *a, Integer tda);
#else
extern void f06thc();
#endif

#ifdef NAG_PROTO
extern void f06ttc(OperationSide side, Integer n, Integer k1, Integer k2, double NAG_HUGE *c,
            Complex NAG_HUGE *s, Complex NAG_HUGE *a, Integer tda);
#else
extern void f06ttc();
#endif

#ifdef NAG_PROTO
extern double f06uaf(char *norm, Integer m, Integer n, Complex *a,
	       Integer lda, double *work);
#else
extern double f06uaf();
#endif

#ifdef NAG_PROTO
extern double f06umf(char *norm, Integer n, Complex *a, Integer lda,
double *work);
#else
extern double f06umf();
#endif

#ifdef NAG_PROTO
extern void f06txc(OperationSide side, PivotType pivot, SequenceDirection direct,
            Integer m, Integer n, Integer k1, Integer k2, double NAG_HUGE *c,
            Complex NAG_HUGE *s, Complex NAG_HUGE *a, Integer tda);
#else
extern void f06txc();
#endif

#ifdef NAG_PROTO
extern double f06vgc(NormType norm, MatrixType matrix, Integer m, Integer n, 
              Complex NAG_HUGE *a, Integer tda);
#else
extern double f06vgc();
#endif

#ifdef NAG_PROTO
extern void f06vkc(OperationSide side, MatrixTranspose trans, Integer n, double NAG_HUGE *perm,
            Integer k, Complex NAG_HUGE *b, Integer tdb);
#else
extern void f06vkc();
#endif

#ifdef NAG_PROTO
extern void f06vxc(OperationSide side, PivotType pivot, SequenceDirection direct,
            Integer m, Integer n, Integer k1, Integer k2, double NAG_HUGE *c, double NAG_HUGE *s,
            Complex NAG_HUGE *a, Integer tda);
#else
extern void f06vxc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06yac(MatrixTranspose transa, MatrixTranspose transb,
            Integer m, Integer n,
            Integer k, double alpha, const double NAG_HUGE a[], Integer tda,
            const double NAG_HUGE b[], Integer tdb, double beta,
            double NAG_HUGE c[], Integer tdc);
#else
extern void f06yac();
#endif

#ifdef NAG_PROTO
extern void f06yaf(char NAG_HUGE *transa, char NAG_HUGE *transb, Integer m, Integer n, Integer k, 
            double alpha, double NAG_HUGE *a, Integer lda, double NAG_HUGE *b,
            Integer ldb, double beta, double NAG_HUGE *c, Integer ldc);
#else
extern void f06yaf();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06ycc(OperationSide side, MatrixTriangle uplo, Integer m, Integer n,
            double alpha, const double NAG_HUGE a[], Integer tda, const double NAG_HUGE b[],
            Integer tdb, double beta, double NAG_HUGE c[], Integer tdc);
#else
extern void f06ycc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06yfc(OperationSide side, MatrixTriangle uplo, MatrixTranspose transa,
            MatrixUnitTriangular diag, Integer m, Integer n, double alpha, 
            const double NAG_HUGE a[], Integer tda, double NAG_HUGE b[], Integer tdb);
#else
extern void f06yfc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06yjc(OperationSide side, MatrixTriangle uplo, MatrixTranspose transa,
            MatrixUnitTriangular diag, Integer m, Integer n, double alpha, 
            const double NAG_HUGE a[], Integer tda, double NAG_HUGE b[], Integer tdb);
#else
extern void f06yjc();
#endif

#ifdef NAG_PROTO
extern void f06yjf(char NAG_HUGE *side, char NAG_HUGE *uplo, char NAG_HUGE *transa, char NAG_HUGE *diag, Integer m,
            Integer n, double alpha, double NAG_HUGE *a, 
            Integer lda, double NAG_HUGE *b, Integer ldb);
#else
extern void f06yjf();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06ypc(MatrixTriangle uplo, MatrixTranspose trans, Integer n, Integer k,
            double alpha, const double NAG_HUGE a[], Integer tda,
            double beta, double NAG_HUGE c[], Integer tdc);
#else
extern void f06ypc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06yrc(MatrixTriangle uplo, MatrixTranspose trans, Integer n, Integer k,
            double alpha, const double NAG_HUGE a[], Integer tda, const double b [],
            Integer tdb, double beta, double NAG_HUGE c[], Integer tdc);
#else
extern void f06yrc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06zac(MatrixTranspose transa, MatrixTranspose transb,
            Integer m, Integer n,
            Integer k, Complex alpha, const Complex NAG_HUGE a[], Integer tda,
            const Complex NAG_HUGE b[], Integer tdb, Complex beta,
            Complex NAG_HUGE c[], Integer tdc);
#else
extern void f06zac();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06zcc(OperationSide side, MatrixTriangle uplo, Integer m, Integer n,
            Complex alpha, const Complex NAG_HUGE a[], Integer tda, const Complex NAG_HUGE b[],
            Integer tdb, Complex beta, Complex NAG_HUGE c[], Integer tdc);
#else
extern void f06zcc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06zfc(OperationSide side, MatrixTriangle uplo, MatrixTranspose transa,
            MatrixUnitTriangular diag, Integer m, Integer n, Complex alpha, 
            const Complex NAG_HUGE a[], Integer tda, Complex NAG_HUGE b[], Integer tdb);
#else
extern void f06zfc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06zjc(OperationSide side, MatrixTriangle uplo, MatrixTranspose transa,
            MatrixUnitTriangular diag, Integer m, Integer n, Complex alpha, 
            const Complex NAG_HUGE a[], Integer tda, Complex NAG_HUGE b[], Integer tdb);
#else
extern void f06zjc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06zpc(MatrixTriangle uplo, MatrixTranspose trans, Integer n, Integer k,
            double alpha, const Complex NAG_HUGE a[], Integer tda,
            double beta, Complex NAG_HUGE c[], Integer tdc);
#else
extern void f06zpc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06zrc(MatrixTriangle uplo, MatrixTranspose trans, Integer n, Integer k,
            Complex alpha, const Complex NAG_HUGE a[], Integer tda, 
            const Complex b [], Integer tdb, double beta, Complex NAG_HUGE c[], 
            Integer tdc);
#else
extern void f06zrc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06ztc(OperationSide side, MatrixTriangle uplo, Integer m, Integer n,
            Complex alpha, const Complex NAG_HUGE a[], Integer tda, const Complex NAG_HUGE b[],
            Integer tdb, Complex beta, Complex NAG_HUGE c[], Integer tdc);
#else
extern void f06ztc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06zuc(MatrixTriangle uplo, MatrixTranspose trans, Integer n, Integer k,
            Complex alpha, const Complex NAG_HUGE a[], Integer tda,
            Complex beta, Complex NAG_HUGE c[], Integer tdc);
#else
extern void f06zuc();
#endif

#ifdef NAG_PROTO
extern NAG_DLL_EXPIMP void NAG_CALL f06zwc(MatrixTriangle uplo, MatrixTranspose trans, Integer n, Integer k,
            Complex alpha, const Complex NAG_HUGE a[], Integer tda, 
            const Complex b [], Integer tdb, Complex beta, Complex NAG_HUGE c[], 
            Integer tdc);
#else
extern void f06zwc();
#endif


#ifdef __cplusplus
}
#endif
#endif /* not NAGF06 */
