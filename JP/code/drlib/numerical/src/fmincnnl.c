#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK l_min_con_nonlin(void (*fcn)(Mint, Mint, Mint,
                                Mfloat[], Mint[], Mfloat*, Mfloat[]),
                                Mint m, Mint me, Mint n, Mint ibtype,
                                Mfloat xlb[], Mfloat xub[],
                                va_list argptr);
static void l_n2onf(void (*fcns)(Mint, Mint, Mint, Mfloat[], Mint[],
                    Mfloat*, Mfloat[]), Mint *m, Mint *me, Mint *n,
                    Mfloat xguess[], Mint *ibtype, Mfloat xlb[],
                    Mfloat xub[], Mfloat xscale[], Mint *iprint,
                    Mint *maxitn, Mfloat x[], Mfloat *fvalue,
                    Mfloat wk[], Mint *lwk, Mint iwk[], Mint *liwk,
                    Mfloat conwk[], Mfloat *acc);
static void l_n2ong(void (*fcns)(Mint, Mint, Mint, Mfloat[], Mint[],
                    Mfloat*, Mfloat[]), void (*grad)(Mint, Mint,
                    Mint, Mint, Mfloat[], Mint[], Mfloat, Mfloat[],
                    Mfloat[], Mfloat[]), Mint *m, Mint *me, Mint *n,
                    Mfloat xguess[], Mint *ibtype, Mfloat xlb[],
                    Mfloat xub[], Mint *iprint, Mint *maxitn,
                    Mfloat x[], Mfloat *fvalue, Mfloat wk[],
                    Mint *lwk, Mint iwk[], Mint *liwk, Mfloat *acc);
static void l_n3onf(void (*fcns)(Mint, Mint, Mint, Mfloat[], Mint[],
                    Mfloat*, Mfloat[]), Mint *m, Mint *me,
                    Mint *mmax, Mint *n, Mint *nmax, Mint *mnn2,
                    Mfloat x[], Mfloat xscale[], Mfloat *f,
                    Mfloat g[], Mfloat df[], Mfloat dg[], Mint *lddg,
                    Mfloat u[], Mfloat xl[], Mfloat xu[], Mfloat *c,
                    Mint *ldc, Mfloat d[], Mfloat conwk[],
                    Mfloat *acc, Mfloat *scbou, Mint *maxfun,
                    Mint *maxit, Mint *iprint, Mint *mode,
                    Mint *ifail, Mfloat wa[], Mint *lwa, Mint kwa[],
                    Mint *lkwa, Mint active[], Mint *lactiv,
                    Mint *llise, Mint *lql);
static void l_n3ong(void (*fcns)(Mint, Mint, Mint, Mfloat[], Mint[],
                    Mfloat*, Mfloat[]), void (*grad)(Mint, Mint,
                    Mint, Mint, Mfloat[], Mint[], Mfloat, Mfloat[],
                    Mfloat[], Mfloat[]), Mint *m, Mint *me,
                    Mint *mmax, Mint *n, Mint *nmax, Mint *mnn2,
                    Mfloat x[], Mfloat *f, Mfloat g[], Mfloat df[],
                    Mfloat dg[], Mint *lddg, Mfloat u[], Mfloat xl[],
                    Mfloat xu[], Mfloat *c, Mint *ldc, Mfloat d[],
                    Mfloat *acc, Mfloat *scbou, Mint *maxfun,
                    Mint *maxit, Mint *iprint, Mint *mode,
                    Mint *ifail, Mfloat wa[], Mint *lwa, Mint kwa[],
                    Mint *lkwa, Mint active[], Mint *lactiv,
                    Mint *llise, Mint *lql);
static void l_n4onf(void (*fcns)(Mint, Mint, Mint, Mfloat[], Mint[],
                    Mfloat*, Mfloat[]), Mint *mmax, Mint *n,
                    Mint *nmax, Mfloat x[], Mfloat xs[], Mfloat g[],
                    Mfloat df[], Mfloat dg[], Mint *lddg, Mfloat u[],
                    Mfloat xl[], Mfloat xu[], Mfloat *dcl,
                    Mint *lddcl, Mfloat cd[], Mfloat cwk[],
                    Mfloat vmu[], Mfloat del[], Mfloat dla[],
                    Mfloat dclf[], Mfloat bdel[], Mfloat eta[],
                    Mfloat xold[], Mfloat dlaold[], Mfloat v[],
                    Mfloat w[], Mfloat vmuold[], Mfloat dphi[],
                    Mfloat rpen[], Mfloat scg[], Mfloat *fbest,
                    Mfloat dfbest[], Mfloat gbest[], Mfloat *dgbest,
                    Mfloat wa[], Mint *lwa, Mint *mnn2, Mint *mo1,
                    Mint *nfunc, Mint *ngrad, Mint *iter, Mint *nql,
                    Mint *iline, Mint *iflise, Mint *nopt, Mint iw[],
                    Mint *liw, Mfloat *phi, Mfloat *dfdel,
                    Mfloat *dbd, Mfloat *alpham, Mfloat *alphao,
                    Mfloat *scf, Mfloat *prd, Mint active[],
                    Mint l7[]);
static void l_n4ong(void (*fcns)(Mint, Mint, Mint, Mfloat[], Mint[],
                    Mfloat*, Mfloat[]), void (*grad)(Mint, Mint,
                    Mint, Mint, Mfloat[], Mint[], Mfloat, Mfloat[],
                    Mfloat[], Mfloat[]), Mint *mmax, Mint *n,
                    Mint *nmax, Mfloat x[], Mfloat g[], Mfloat df[],
                    Mfloat dg[], Mint *lddg, Mfloat u[], Mfloat xl[],
                    Mfloat xu[], Mfloat *dcl, Mint *lddcl,
                    Mfloat cd[], Mfloat vmu[], Mfloat del[],
                    Mfloat dla[], Mfloat dclf[], Mfloat bdel[],
                    Mfloat eta[], Mfloat xold[], Mfloat dlaold[],
                    Mfloat v[], Mfloat w[], Mfloat vmuold[],
                    Mfloat dphi[], Mfloat rpen[], Mfloat scg[],
                    Mfloat *fbest, Mfloat dfbest[], Mfloat gbest[],
                    Mfloat *dgbest, Mfloat wa[], Mint *lwa,
                    Mint *mnn2, Mint *mo1, Mint *nfunc, Mint *ngrad,
                    Mint *iter, Mint *nql, Mint *iline, Mint *iflise,
                    Mint *nopt, Mint iw[], Mint *liw, Mfloat *phi,
                    Mfloat *dfdel, Mfloat *dbd, Mfloat *alpham,
                    Mfloat *alphao, Mfloat *scf, Mfloat *prd,
                    Mint active[], Mint l7[]); 
static void l_n5onf(void (*fcn)(Mint, Mint, Mint, Mfloat[], Mint[],
                    Mfloat*, Mfloat[]), Mint *m, Mint *me,
                    Mint *mmax, Mint *n, Mfloat xc[],
                    Mfloat xscale[], Mint active[], Mfloat *fc,
                    Mfloat gc[], Mfloat df[], Mfloat dg[],
                    Mfloat work[]);
static void l_n5ong(Mint *mode, Mint *m, Mint *me, Mint *n,
                    Mint *mnn, Mint *nmnn, Mfloat *acc, Mfloat r[],
                    Mfloat *f, Mfloat df[], Mfloat g[], Mfloat dg[],
                    Mint *lddg, Mfloat v[], Mfloat u[], Mfloat x[],
                    Mfloat xl[], Mfloat xu[], Mfloat *phi, 
                    Mfloat dphi[], Mint active[], Mfloat wa[],
                    Mint *lwa); 
static void l_n6ong(Mint *m, Mint *me, Mint *mmax, Mint *n,
                    Mint *nmax, Mint *mnn, Mfloat *c, Mint *ldc,
                    Mfloat d[], Mfloat *a, Mint *lda, Mfloat b[],
                    Mfloat xl[], Mfloat xu[], Mfloat x[], Mfloat u[],
                    Mint *ifail, Mint *iprint, Mfloat war[],
                    Mint *lwar, Mint iwar[], Mint *liwar); 
static void l_n7ong(Mint *n, Mfloat *cl, Mint *ldcl, Mfloat d[],
                    Mfloat u[], Mfloat v[]);
static void l_n8ong(Mfloat *alpha, Mfloat *alpham, Mfloat *phi,
                    Mfloat *dphi, Mfloat *amue, Mfloat *imsl_beta,
                    Mint *iline, Mint *maxfun, Mint *ifail,
                    Mint *iprint, Mfloat wa[], Mint *lwa, Mint kwa[],
                    Mint *lkwa, Mint lowa[], Mint *llowa);
static void l_n9ong(Mint *n, Mint *m, Mint *meq, Mint *mmax,
                    Mint *mn, Mint *mnn, Mint *nmax, Mint *lql,
                    Mfloat *a, Mint *lda, Mfloat b[], Mfloat grad[],
                    Mfloat *g, Mint *ldg, Mfloat xl[], Mfloat xu[],
                    Mfloat x[], Mint *nact, Mint iact[], Mint *info,
                    Mfloat *diag, Mfloat w[], Mint *lw);
static Mfloat l_a1ot(Mint n, Mfloat sx[], Mint incx, Mfloat sy[],
                     Mint incy);
static Mint   l_ismax(Mint n, Mfloat sx[], Mint incx);
static Mfloat l_ssum(Mint *n, Mfloat sx[], Mint *incx);
static void   l_shprod(Mint *n, Mfloat sx[], Mint *incx, Mfloat sy[],
                       Mint *incy, Mfloat sz[], Mint *incz);
#else
static VA_LIST_HACK      l_min_con_nonlin();
static void l_n2onf();
static void l_n2ong();
static void l_n3onf();
static void l_n3ong();
static void l_n4onf();
static void l_n4ong();
static void l_n5onf();
static void l_n5ong();
static void l_n6ong();
static void l_n7ong();
static void l_n8ong();
static void l_n9ong();
static Mfloat l_a1ot();
static void l_shprod();
static Mint l_ismax();
static Mfloat l_ssum();
#endif

static Mfloat       *lv_value;

#ifdef ANSI
#if defined(COMPUTER_HP97C)
Mfloat *imsl_f_min_con_nonlin(void (*fcn) (Mint, Mint, Mint,
                              Mfloat*, Mint*, Mfloat*, Mfloat*),
                              Mint m, Mint me, Mint n, Mint ibtype,
                              Mfloat xlb[], Mfloat xub[], ...)
#else
Mfloat *imsl_f_min_con_nonlin(void (*fcn) (Mint, Mint, Mint,
                              Mfloat[], Mint[], Mfloat*, Mfloat[]),
                              Mint m, Mint me, Mint n, Mint ibtype,
                              Mfloat xlb[], Mfloat xub[], ...)
#endif
#else
Mfloat *imsl_f_min_con_nonlin(fcn, m, me, n, ibtype, xlb, xub,
                              va_alist)
    void     (*fcn) ();
    Mint     m;
    Mint     me; 
    Mint     n;
    Mint     ibtype;
    Mfloat   xlb[];
    Mfloat   xub[];
    va_dcl
#endif
{
    va_list argptr;

    VA_START(argptr, xub);
    E1PSH("imsl_f_min_con_nonlin", "imsl_d_min_con_nonlin");
    lv_value = NULL;
    IMSL_CALL(l_min_con_nonlin(fcn, m, me, n, ibtype, xlb, xub,
                               argptr));
    va_end(argptr);
    E1POP("imsl_f_min_con_nonlin", "imsl_d_min_con_nonlin"); 
    return lv_value;
}


#ifdef ANSI
#if defined(COMPUTER_HP97C)
static VA_LIST_HACK l_min_con_nonlin(void (*fcn) (Mint, Mint, Mint,
                                Mfloat*, Mint*, Mfloat*, Mfloat*),
                                Mint m, Mint me, Mint n, Mint ibtype,
                                Mfloat xlb[], Mfloat xub[],
                                va_list argptr)
#else
static VA_LIST_HACK l_min_con_nonlin(void (*fcn) (Mint, Mint, Mint,
                                Mfloat[], Mint[], Mfloat*, Mfloat[]),
                                Mint m, Mint me, Mint n, Mint ibtype,
                                Mfloat xlb[], Mfloat xub[],
                                va_list argptr)
#endif
#else
static VA_LIST_HACK l_min_con_nonlin(fcn, m, me, n, ibtype, xlb, xub,
                                argptr)
    void     (*fcn) ();
    Mint     m;
    Mint     me;
    Mint     n;
    Mint     ibtype;
    Mfloat   xlb[];
    Mfloat   xub[];
    va_list       argptr;
#endif
{
    Mint          i, code;
    Mint          arg_number    = 7;
    Mint          itmax         = 200;
    Mint          iprint        = 0;
    Mfloat        *xguess       = NULL;
    Mfloat        *xguess_float = NULL;
    Mint          user_xguess   = 0;
    Mfloat        *xscale_float = NULL;
    Mfloat        *xscale       = NULL;
    Mint          user_xscale   = 0;
    Mfloat        *x_float      = NULL;
    Mfloat        return_user   = 0;
    Mfloat        *xlb_float    = NULL;
    Mfloat        *xub_float    = NULL;
    Mfloat        err_rel;
    Mfloat        *obj          = NULL;
    Mfloat        *wk           = NULL;
    Mint          lwk;
    Mfloat        *conwk        = NULL;
    Mint          *iwk          = NULL;
    Mint          liwk;
#ifdef ANSI
    void        (*grad)(Mint, Mint, Mint, Mint, Mfloat[], Mint[], Mfloat, Mfloat[],
                        Mfloat[], Mfloat[]);
#else
    void        (*grad)();
#endif

    Mint          user_gradient = 0;
    Mint          obj_user      = 0;
    Mint          mmax;
    Mint          mnmax;

    err_rel  = sqrt(imsl_amach(4));

    code = 1;
    while (code > 0) {
        code = va_arg(argptr, Mint);
        arg_number++;
        switch(code) {
            case IMSL_XGUESS:
                arg_number++;
                xguess = va_arg(argptr, Mfloat*);
                user_xguess = 1;
                break;
            case IMSL_GRADIENT:
                arg_number++;
#ifdef ANSI
                grad = (void (*)(Mint, Mint, Mint, Mint, 
			Mfloat[], Mint[], Mfloat, Mfloat[], 
			Mfloat[], Mfloat[])) va_arg(argptr, void*);
#else
                grad = (void (*)())va_arg(argptr, void*);
#endif
                user_gradient = 1;
                break;
            case IMSL_ERR_REL:
                arg_number++; 
                err_rel = (Mfloat) va_arg(argptr, Mdouble);
                break;
            case IMSL_ERR_REL_ADR:
                arg_number++; 
                err_rel = *(va_arg(argptr, Mfloat *));
                break;
            case IMSL_XSCALE:
                arg_number++;
                xscale = va_arg(argptr, Mfloat*);
                user_xscale = 1;
                break; 
            case IMSL_PRINT:
                arg_number++;
                iprint = va_arg(argptr, Mint);
                break;
            case IMSL_ITMAX:
                arg_number++; 
                itmax = va_arg(argptr, Mint);
                break;
            case IMSL_RETURN_USER:
                arg_number++;
                lv_value = va_arg(argptr, Mfloat*);
                return_user = 1;
                break;
            case IMSL_OBJ:
                arg_number++;
                obj = va_arg(argptr, Mfloat*);
                obj_user = 1;
                break;
            case 0:
                break;
            default:
                imsl_e1sti (1, code);
                imsl_e1sti (2, arg_number);
                imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
                break;
        }
    }

    if (imsl_n1rty(0)) goto RETURN;

    if (n <= 0) {
        imsl_e1sti(1, n);
        imsl_ermes(IMSL_TERMINAL, IMSL_N_MUST_BE_POSITIVE);
        goto FREE_SPACE;
    }

    if (m < 0) {
        imsl_e1sti(1, m);
        imsl_ermes(IMSL_TERMINAL, IMSL_NONNEGATIVE_CONSTRAINTS);
        goto FREE_SPACE;
    } 

    if (ibtype < 0 || ibtype >= 4) {
        imsl_e1sti(1, ibtype);
        imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_IBTYPE_VALUE);
        goto FREE_SPACE;
    } 

    if (user_xscale) {
        for (i = 0; i < n; i++) {
             if (xscale[i] <= F_ZERO) {
                 imsl_e1sti(1, i);
                 imsl_e1str(1, xscale[i]);
                 imsl_ermes(IMSL_TERMINAL, IMSL_XSCALE_DIAGONAL_LT_ZERO);
                 goto FREE_SPACE;
             }
        }
    }

    mmax  = imsl_i_max(1, m);
    mnmax = imsl_i_max(m, n);
    liwk  = 19 + mnmax;
    lwk   = n*(3*n+38+mmax) + 6*mmax + 6*m + 72;

    wk = (Mfloat *) imsl_malloc (lwk*sizeof(*wk));
    conwk = (Mfloat *) imsl_malloc (mmax*sizeof(*conwk));
    iwk = (Mint *) imsl_malloc (liwk*sizeof(*iwk));
    x_float = (Mfloat *) imsl_malloc (n*sizeof(*x_float));
    xlb_float = (Mfloat *) imsl_malloc (n*sizeof(*xlb_float));
    xub_float = (Mfloat *) imsl_malloc (n*sizeof(*xub_float));
    xguess_float = (Mfloat *) imsl_malloc (n*sizeof(*xguess_float)); 
    xscale_float = (Mfloat *) imsl_malloc (n*sizeof(*xscale_float));

    if (!obj_user)  obj = (Mfloat *) imsl_malloc (1*sizeof(*obj));

    if (wk == NULL || conwk == NULL || iwk == NULL) {
        imsl_e1stl(1, "n");
        imsl_e1sti(1, n);
        imsl_e1stl(2, "m");
        imsl_e1sti(2, m);
        imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
        goto FREE_SPACE;
    }

    if (x_float == NULL || xlb_float == NULL || xub_float == NULL ||
        xguess_float == NULL || xscale_float == NULL || obj == NULL) {
          imsl_e1sti(1, n);
          imsl_e1stl(1, "n");
          imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
          goto FREE_SPACE;
    }

    if (user_xguess) 
       for (i=0; i<n; i++)  xguess_float[i] = (Mfloat) xguess[i];
    else
       for (i=0; i<n; i++)  xguess_float[i] = F_ZERO;

    if (user_xscale) 
       for (i=0; i<n; i++)  xscale_float[i] = (Mfloat) xscale[i]; 
    else
       for (i=0; i<n; i++)  xscale_float[i] = F_ONE; 
    
    for (i=0; i<n; i++)  xlb_float[i] = (Mfloat) xlb[i];
    for (i=0; i<n; i++)  xub_float[i] = (Mfloat) xub[i];

    if (lv_value == NULL) {
       lv_value = (Mfloat *) imsl_malloc (n*sizeof(*lv_value));
       if (lv_value == NULL) {
          imsl_e1sti(1, n);
          imsl_e1stl(1, "n");
          imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
          goto FREE_SPACE;
       }
    }

    if (user_gradient)
       l_n2ong(fcn, grad, &m, &me, &n, xguess_float, &ibtype,
               xlb_float, xub_float, &iprint, &itmax, x_float, obj,
               wk, &lwk, iwk, &liwk, &err_rel); 
    else
       l_n2onf(fcn, &m, &me, &n, xguess_float, &ibtype, xlb_float,
               xub_float, xscale_float, &iprint, &itmax, x_float,
               obj, wk, &lwk, iwk, &liwk, conwk, &err_rel);

    for (i=0; i<n; i++)  lv_value[i] = x_float[i];

FREE_SPACE:
    if (wk != NULL)                imsl_free (wk);
    if (iwk != NULL)               imsl_free (iwk);
    if (conwk != NULL)             imsl_free (conwk);
    if (x_float != NULL)           imsl_free (x_float);
    if (xguess_float != NULL)      imsl_free (xguess_float);
    if (xscale_float != NULL)      imsl_free (xscale_float);
    if (xlb_float != NULL)         imsl_free (xlb_float);
    if (xub_float != NULL)         imsl_free (xub_float);
    if (obj != NULL && !obj_user)  imsl_free (obj);

RETURN:
    if (imsl_n1rty(0) > 4) {
        if (!return_user && lv_value != NULL) 
            imsl_free(lv_value);
        lv_value = NULL;
    }
    return (argptr);
}

/* -----------------------------------------------------------------------
    IMSL Name:  N2ONF/DN2ONF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 11, 1988

    Purpose:    Solve a general nonlinear programming problem using the
                successive quadratic programming algorithm and a finite
                difference gradient.

    Usage:      CALL N2ONF (FCNS, M, ME, N, XGUESS, IBTYPE, XLB, XUB,
                            XSCALE, IPRINT, MAXITN, X, FVALUE, WK, LWK,
                            IWK, LIWK, CONWK)

    Arguments:  See NCONF.

    Remarks:    See NCONF.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_n2onf(void (*fcns)(Mint, Mint, Mint, Mfloat*, Mint*,
                    Mfloat*, Mfloat*), Mint *m, Mint *me, Mint *n,
                    Mfloat xguess[], Mint *ibtype, Mfloat xlb[],
                    Mfloat xub[], Mfloat xscale[], Mint *iprint,
                    Mint *maxitn, Mfloat x[], Mfloat *fvalue, Mfloat wk[],
                    Mint *lwk, Mint iwk[], Mint *liwk, Mfloat conwk[], 
                    Mfloat *acc)
#else
static void l_n2onf(void (*fcns)(Mint, Mint, Mint, Mfloat[], Mint[],
                    Mfloat*, Mfloat[]), Mint *m, Mint *me, Mint *n,
                    Mfloat xguess[], Mint *ibtype, Mfloat xlb[],
                    Mfloat xub[], Mfloat xscale[], Mint *iprint,
                    Mint *maxitn, Mfloat x[], Mfloat *fvalue, Mfloat wk[],
                    Mint *lwk, Mint iwk[], Mint *liwk, Mfloat conwk[], 
                    Mfloat *acc)
#endif
#else
static void l_n2onf(fcns, m, me, n, xguess, ibtype, xlb, xub, xscale, iprint,
                    maxitn, x, fvalue, wk, lwk, iwk, liwk, conwk, acc)
	void            (*fcns) ();
	Mint            *m, *me, *n;
	Mfloat           xguess[];
	Mint            *ibtype;
	Mfloat           xlb[], xub[], xscale[];
	Mint            *iprint, *maxitn;
	Mfloat           x[], *fvalue, wk[];
	Mint            *lwk, iwk[], *liwk;
	Mfloat           conwk[], *acc;
#endif
{
        FILE            *iout;
	Mint            active[1000], llise, lql;
	Mint            icd, idcl, idf, idg, ifail, ig, iu, kwk, lactiv, 
                        lwk1, m1, maxfun, mnn2, mode, n1;
	Mfloat          scbou;
        Mint            IMSLTRUE = 1, IMSLFALSE = 0;


	imsl_e1psh("N2ONF ");

	/* CHECK BOUND CONDITIONS */
	if (*ibtype == 1) {
	        sset(*n, F_ZERO, xlb, 1);
	} else if (*ibtype == 2) {
	        sset(*n, F_ZERO, xub, 1);
	} else if (*ibtype == 3) {
	        sset(*n - 1, xlb[0], &xlb[1], 1);
		sset(*n - 1, xub[0], &xub[1], 1);
	}

	scopy(*n, xguess, 1, x, 1);
	scbou = 1.0e3;
	maxfun = 5;
	mode = 0;
	ifail = 0;
	lql = IMSLFALSE;
	if (*n < 20)
	    lql = IMSLTRUE;
	llise = IMSLTRUE;
	imsl_umach(2, &iout);

	/* INITIAL ADDRESSES IN WK */
	n1 = *n + 1;
	m1 = imsl_i_max(1, *m);
	mnn2 = *m + n1 + n1;
	ig = 1;
	idf = ig + m1;
	idg = idf + *n;
	iu = idg + m1 * n1;
	idcl = iu + mnn2;
	icd = idcl + n1 * n1;
	kwk = icd + n1;
	lwk1 = *lwk - kwk + 1;
	lactiv = 200;

	l_n3onf(fcns, m, me, &m1, n, &n1, &mnn2, x, xscale, fvalue, &wk[ig - 1],
                &wk[idf - 1], &wk[idg - 1], &m1, &wk[iu - 1], xlb, xub,
                &wk[idcl - 1], &n1, &wk[icd - 1], conwk, acc, &scbou, &maxfun,
                maxitn, iprint, &mode, &ifail, &wk[kwk - 1], &lwk1, iwk, liwk,
                active, &lactiv, &llise, &lql);
	imsl_e1pop("N2ONF ");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  N2ONG/DN2ONG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    March 11, 1988

    Purpose:    Solve a general nonlinear programming problem using the
                successive quadratic programming algorithm and a
                user-supplied gradient.

    Usage:      CALL N2ONG (FCNS, GRAD, M, ME, N, XGUESS, IBTYPE, XLB,
                            XUB, IPRINT, MAXITN, X, FVALUE, WK, LWK,
                            IWK, LIWK)

    Arguments:  See NCONG.

    Remarks:    See NCONG.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_n2ong(void (*fcns)(Mint, Mint, Mint, Mfloat*, Mint*,
                    Mfloat*, Mfloat*), void (*grad)(Mint, Mint, Mint, Mint,
                    Mfloat*, Mint*, Mfloat, Mfloat*, Mfloat*, Mfloat*),
                    Mint *m, Mint *me, Mint *n, Mfloat xguess[], Mint *ibtype,
                    Mfloat xlb[], Mfloat xub[], Mint *iprint, Mint *maxitn,
                    Mfloat x[], Mfloat *fvalue, Mfloat wk[], Mint *lwk,
                    Mint iwk[], Mint *liwk, Mfloat *acc)
#else
static void l_n2ong(void (*fcns)(Mint, Mint, Mint, Mfloat[], Mint[],
                    Mfloat*, Mfloat[]), void (*grad)(Mint, Mint, Mint, Mint,
                    Mfloat[], Mint[], Mfloat, Mfloat[], Mfloat[], Mfloat[]),
                    Mint *m, Mint *me, Mint *n, Mfloat xguess[], Mint *ibtype,
                    Mfloat xlb[], Mfloat xub[], Mint *iprint, Mint *maxitn,
                    Mfloat x[], Mfloat *fvalue, Mfloat wk[], Mint *lwk,
                    Mint iwk[], Mint *liwk, Mfloat *acc)
#endif
#else
static void l_n2ong(fcns, grad, m, me, n, xguess, ibtype, xlb, xub, iprint,
                    maxitn, x, fvalue, wk, lwk, iwk, liwk, acc)
	void            (*fcns) (), (*grad) ();
	Mint            *m, *me, *n;
	Mfloat           xguess[];
	Mint            *ibtype;
	Mfloat           xlb[], xub[];
	Mint            *iprint, *maxitn;
	Mfloat           x[], *fvalue, wk[];
	Mint            *lwk, iwk[], *liwk;
        Mfloat          *acc;
#endif
{
	Mint            active[1000], llise, lql;
	Mint            icd, idcl, idf, idg, ifail, ig, iu, kwk,
	                lactiv, lwk1, m1, maxfun, mnn2, mode, n1;
	Mfloat          scbou;
        Mint            IMSLTRUE = 1;


	imsl_e1psh("N2ONG ");

	/* CHECK BOUND CONDITIONS */
	if (*ibtype == 1) {
	        sset(*n, F_ZERO, xlb, 1);
	} else if (*ibtype == 2) {
		sset(*n, F_ZERO, xub, 1);
	} else if (*ibtype == 3) {
		sset(*n - 1, xlb[0], &xlb[1], 1);
		sset(*n - 1, xub[0], &xub[1], 1);
	}

	scopy(*n, xguess, 1, x, 1);
	scbou = 1.0e3;
	maxfun = 5;
	mode = 0;
	ifail = 0;
	lql = IMSLTRUE;
	llise = IMSLTRUE;

	/* INITIAL ADDRESSES IN WK */
	n1 = *n + 1;
	m1 = imsl_i_max(1, *m);
	mnn2 = *m + n1 + n1;
	ig = 1;
	idf = ig + m1;
	idg = idf + *n;
	iu = idg + m1 * n1;
	idcl = iu + mnn2;
	icd = idcl + n1 * n1;
	kwk = icd + n1;
	lwk1 = *lwk - kwk + 1;
	lactiv = 200;

	l_n3ong(fcns, grad, m, me, &m1, n, &n1, &mnn2, x, fvalue, &wk[ig - 1],
                &wk[idf - 1], &wk[idg - 1], &m1, &wk[iu - 1], xlb, xub,
                &wk[idcl - 1], &n1, &wk[icd - 1], acc, &scbou, &maxfun, maxitn,
                iprint, &mode, &ifail, &wk[kwk - 1], &lwk1, iwk, liwk, active,
                &lactiv, &llise, &lql);

	imsl_e1pop("N2ONG ");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  N3ONF/DN3ONF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    December 2, 1985

    Purpose:    Update the Cholesky factorization.

    Usage:      CALL N3ONF (FCNS, M, ME, MMAX, N, NMAX, MNN2, X,
                            XSCALE, F, G, DF, DG, LDDG, U, XL, XU, C,
                            LDC, D, CONWK, ACC, SCBOU, MAXFUN, MAXIT,
                            IPRINT, MODE, IFAIL, WA, LWA, KWA, LKWA,
                            ACTIVE, LACTIV, LLISE, LQL)

    Arguments:
       FCNS   - User-supplied SUBROUTINE to evaluate the functions at
                a given poMint.  The usage is
                CALL FCNS (M, ME, MMAX, N, X, ACTIVE, F, G), where
                M      - Total number of constraMints.  (Input)
                ME     - Number of equality constraMints.  (Input)
                MMAX   - Maximum of (1, M).  (Input)
                N      - Number of variables.  (Input)
                X      - The poMint at which the function is evaluated.
                         (Input)
                         X should not be changed by FCNS.
                ACTIVE - Logical vector of length MMAX indicating the
                         active constraMints.  (Input)
                F      - The computed function value at the poMint X.
                         (Output)
                G      - Vector of length MMAX containing the values of
                         constraMints at poMint X.  (Output)
                FCNS must be declared EXTERNAL in the calling program.
       M      - Number of constraMints.  (Input)
       ME     - Number of equality constraMints.  (Input)
       MMAX   - Leading dimension of DG.  (Input)
                MMAX must be at least MAX(1,M).
       N      - Number of variables.  (Input)
       NMAX   - Leading dimension of C.  (Input)
                NMAX must be at least MAX(2,N).
       MNN2   - Length of U where MNN2 = M + N + N + 2.  (Input)
       X      - Vector of length N containing the initial guesses to the
                solution on input and the solution on output.
                (Input/Output)
       F      - Scalar containing the objective function value.  (Output)
       G      - Vector of length MMAX containing constraMint values.
                (Output)
       DF     - Vector of length NMAX containing the gradient of the
                of the objective function.  (Output)
       DG     - Array of dimension MMAX by MMAX containing the gradient
                of the constraMints.  (Output)
       LDDG   - Leading dimension of DG exactly as specified in the
                dimension statement of the calling program.  (Input)
       U      - Vector of length MNN2 containing the multipliers of the
                nonlinear constraMints and the bounds.  (Output)
       XL     - Vector of length N containing the lower bounds for the
                variables.  (Input)
       XU     - Vector of length N containing the upper bounds for the
                variables.  (Input)
       C      - Array of dimension NMAX by NMAX containing an the final
                approximation to the Hessian.  (Output)
       LDC    - Leading dimension of C exactly as specified in the
                dimension statement of the calling program.  (Input)
       D      - Vector of length NMAX containing the diagonal elements of
                the Hessian.  (Output)
       CONWK  - Work vector of length M.
       ACC    - Final accuracy.  (Input)
       SCBOU  - Scalar containing the scaling variable for the problem
                function.  (Input)
       MAXFUN - Scalar containing the maximum allowable function calls
                during the line search.  (Input)
       MAXIT  - Scalar containing the maximum allowable iterations.
                (Input)
       IPRINT - Specification for the desired print level.  (Input)
       MODE   - Desired solving version for the NLPQL1 algorithm.
                (Input)
       IFAIL  - Scalar containing error message information.  (Output)
       WA     - Work vector of length LWA.
       LWA    - Length of WA where LWA = 2*N*(N+16) + 4*MMAX + 5*M + 66.
                (Input)
       KWA    - Work vector of length LKWA.
       LKWA   - Length of KWA where LKWA = 19.  (Input)
       ACTIVE - Logical vector of length LACTIV containing information
                as to which constraMints are active.  (Output)
       LACTIV - Length of ACTIVE where LACTIV must be at least 200.
                (Input)
       LLISE  - Logical scalar determining if the line search algorithm
                is to be used if LLISE = TRUE.  (Input)
                Otherwise a user supplied line search is used.
       LQL    - Logical scalar determing whether the qudratic programming
                subproblem is to be solved with a positive definite
                Quasi-Newton matrix if LQL = TRUE.  (Input)
                Otherwise a Cholesky-decomposition is performed and
                updated.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
static struct t_n10nf{
	Mfloat           f1, acc1, scbou1, dbdfac, zefac, rpeno, rpens,
	                rpenu, zefacu, delta, imsl_beta, amue, alm;
	Mint             mm1, me1, maxfu1, maxit1, iprin1, mode1, ifail1;
	Mint            llise1, lql1, lmer;
}               n10nf;
static struct t_n11nf {
	Mint             n1, lact, no1, mnn, nmnn;
}               n11nf;
/* LACTIV is not used here, but leave the calling sequence intact.*/

#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_n3onf(void (*fcns)(Mint, Mint, Mint, Mfloat*, Mint*,
                    Mfloat*, Mfloat*), Mint *m, Mint *me, Mint *mmax, Mint *n,
                    Mint *nmax, Mint *mnn2, Mfloat x[], Mfloat xscale[], 
                    Mfloat *f, Mfloat g[], Mfloat df[], Mfloat dg[], Mint *lddg,
                    Mfloat u[], Mfloat xl[], Mfloat xu[], Mfloat *c, Mint *ldc,
                    Mfloat d[], Mfloat conwk[], Mfloat *acc, Mfloat *scbou,
                    Mint *maxfun, Mint *maxit, Mint *iprint, Mint *mode,
                    Mint *ifail, Mfloat wa[], Mint *lwa, Mint kwa[], Mint *lkwa,
                    Mint active[], Mint *lactiv, Mint *llise, Mint *lql)
#else
static void l_n3onf(void (*fcns)(Mint, Mint, Mint, Mfloat[], Mint[],
                    Mfloat*, Mfloat[]), Mint *m, Mint *me, Mint *mmax, Mint *n,
                    Mint *nmax, Mint *mnn2, Mfloat x[], Mfloat xscale[], 
                    Mfloat *f, Mfloat g[], Mfloat df[], Mfloat dg[], Mint *lddg,
                    Mfloat u[], Mfloat xl[], Mfloat xu[], Mfloat *c, Mint *ldc,
                    Mfloat d[], Mfloat conwk[], Mfloat *acc, Mfloat *scbou,
                    Mint *maxfun, Mint *maxit, Mint *iprint, Mint *mode,
                    Mint *ifail, Mfloat wa[], Mint *lwa, Mint kwa[], Mint *lkwa,
                    Mint active[], Mint *lactiv, Mint *llise, Mint *lql)
#endif
#else
static void l_n3onf(fcns, m, me, mmax, n, nmax, mnn2, x, xscale, f, g, df, dg, 
                    lddg, u, xl, xu, c, ldc, d, conwk, acc, scbou, maxfun,
                    maxit, iprint, mode, ifail, wa, lwa, kwa, lkwa, active,
                    lactiv, llise, lql)
	void            (*fcns) ();
	Mint            *m, *me, *mmax, *n, *nmax, *mnn2;
	Mfloat           x[], xscale[], *f, g[], df[], dg[];
	Mint            *lddg;
	Mfloat           u[], xl[], xu[], *c;
	Mint            *ldc;
	Mfloat           d[], conwk[], *acc, *scbou;
	Mint            *maxfun, *maxit, *iprint, *mode, *ifail;
	Mfloat           wa[];
	Mint            *lwa, kwa[], *lkwa;
	Mint            active[];
	Mint            *lactiv;
	Mint           *llise, *lql;
#endif
{
#define DG(I_,J_)	(dg+(I_)*(*lddg)+(J_))
#define C(I_,J_)	(c+(I_)*(*ldc)+(J_))
	Mint            rest;
        FILE            *iout;
	Mint             i, i1, i2, i3, i4, i5, i6, ibdel, iclf, idel, idfb,
	                idgb, idla, idlaol, idphi, ieta, ifb, igb,
	                ipen, iscf, iscg, iv, ivmu, ivmuol, iw,
	                iwaql, ixold, lkwa0, lwaql, m1,maxit2, mo1;
        Mint            IMSLFALSE = 0, IMSLTRUE = 1;

	/*
	 * PROVIDE CONSTANT PARAMETERS FOR CALLING N4ONF
	 */
	imsl_e1psh("N3ONF ");
	n10nf.dbdfac = 2.0e-1;
	n10nf.zefac = F_ONE;
	if (!*lql)
		n10nf.zefac = sqrt(n10nf.zefac);
	n10nf.rpeno = F_TEN;
	n10nf.rpens = 1.0e-2;
	n10nf.zefacu = 1.0e6;
	if (!*lql)
		n10nf.zefacu = sqrt(n10nf.zefacu);
	n10nf.rpenu = 1.0e9;
	n10nf.delta = 9.0e-1;
	n10nf.imsl_beta = 1.0e-1;
	n10nf.amue = 1.0e-4;
	n10nf.alm = F_ONE;
	imsl_umach(2, &iout);
	/* PRINT PARAMETERS */
	if (*iprint == 0)
		goto L_10;
	if (*ifail < 0)
		goto L_10;
	fprintf(stdout, "\n\n\n\n   ---------------------------------------------------------------------------\n     START OF THE SEQUENTIAL QUADRATIC PROGRAMMING ALGORITHM\n   ---------------------------------------------------------------------------\n");
	fprintf(stdout, "\n     PARAMETERS:\n        MODE =%3ld\n        ACC =%13.4e\n        SCBOU =%13.4e\n        MAXFUN =%3ld\n        MAXIT =%5ld\n        IPRINT =%4ld\n",
		*mode, *acc, *scbou, *maxfun, *maxit, *iprint);
	if (*iprint != 2)
		goto L_10;
	fprintf(stdout, "\n     OUTPUT IN THE FOLLOWING ORDER:\n        IT    - ITERATION NUMBER\n        F     - OBJECTIVE FUNCTION VALUE\n        SCV   - SUM OF CONSTRAINT VIOLATION");
        fprintf(stdout, "\n        NA    - NUMBER OF ACTIVE CONSTRAINTS\n        I     - NUMBER OF LINE SEARCH ITERATIONS\n        ALPHA - STEPLENGTH PARAMETER\n        DELTA - ADDITIONAL VARIABLE TO PREVENT INCONSISTENCY");
        fprintf(stdout, "\n        DLAN  - MAXIMUM NORM OF LAGRANGIAN GRADIENT\n        KT    - KUHN-TUCKER OPTIMALITY CRITERION\n");
	fprintf(stdout, "\n\n  IT         F          SCV     NA  I    ALPHA    DELTA      DLAN      KT    \n  ---------------------------------------------------------------------------\n");
L_10:
	;
	n10nf.mode1 = *mode;
	n10nf.maxit1 = *maxit;
	maxit2 = *maxit / 2;
	rest = IMSLFALSE;
	if (*mode < 10)
		goto L_20;
	rest = IMSLTRUE;
	n10nf.mode1 = *mode - 10;
L_20:
	;
	/* INITIAL ADDRESSES IN WA */
	n11nf.n1 = *n + 1;
	n11nf.lact = 2 ** mmax + 6;
	if (*ifail < 0)
		n10nf.lmer = active[n11nf.lact - 1];
	m1 = n11nf.lact + 1;
	n11nf.no1 = *n;
	if (*llise)
		n11nf.no1 = 1;
	mo1 = *mmax;
	if (*llise)
		mo1 = 1;
	n11nf.mnn = *m + *n + *n;
	n11nf.nmnn = *n + n11nf.mnn;
	iscg = 1;
	iscf = iscg + *mmax;
	i1 = iscf + 1;
	i2 = i1 + 1;
	i3 = i2 + 1;
	i4 = i3 + 1;
	i5 = i4 + 1;
	i6 = i5 + 1;
	ivmu = i6 + 1;
	idel = ivmu + n11nf.mnn;
	idla = idel + n11nf.n1;
	iclf = idla + *n;
	ibdel = iclf + n11nf.n1;
	ieta = ibdel + *n;
	ixold = ieta + *n;
	idlaol = ixold + *n;
	iv = idlaol + *n;
	iw = iv + n11nf.n1;
	ivmuol = iw + n11nf.n1;
	idphi = ivmuol + n11nf.mnn;
	ipen = idphi + n11nf.mnn + *n;
	ifb = ipen + n11nf.mnn;
	idfb = ifb + 1;
	igb = idfb + n11nf.no1;
	idgb = igb + *mmax;
	iwaql = idgb + mo1 ** n;
	lwaql = *lwa - iwaql;

	/* CHECK INITIAL POINT */
	if (*ifail >= 0) {
		if (n10nf.mode1 == 2 || n10nf.mode1 == 3)
			goto L_40;
		if (n10nf.mode1 == 7 || n10nf.mode1 == 8)
			goto L_40;
                *f = n10nf.alm;
L_40:
		for (i = 1; i <= *n; i++) {
			if (x[i - 1] < xl[i - 1])
				x[i - 1] = xl[i - 1];
			if (x[i - 1] > xu[i - 1])
				x[i - 1] = xu[i - 1];
		}
	}
	lkwa0 = *lkwa - 7;
	/* CALL OF N4ONF */
	if (*ifail >= 0)
		n10nf.lmer = IMSLTRUE;
L_60:
	;
	n10nf.mm1 = *m;
	n10nf.me1 = *me;
	n10nf.f1 = *f;
	n10nf.acc1 = *acc;
	n10nf.scbou1 = *scbou;
	n10nf.maxfu1 = *maxfun;
	n10nf.iprin1 = *iprint;
	n10nf.ifail1 = *ifail;
	n10nf.llise1 = *llise;
	n10nf.lql1 = *lql;
	l_n4onf(fcns, mmax, n, nmax, x, xscale, g, df, dg, lddg, u, xl, xu, c,
                ldc, d, conwk, &wa[ivmu - 1], &wa[idel - 1], &wa[idla - 1],
	         &wa[iclf - 1], &wa[ibdel - 1], &wa[ieta - 1], &wa[ixold - 1],
		 &wa[idlaol - 1], &wa[iv - 1], &wa[iw - 1], &wa[ivmuol - 1],
                 &wa[idphi - 1], &wa[ipen - 1], &wa[iscg - 1], &wa[ifb - 1],
                 &wa[idfb - 1], &wa[igb - 1], &wa[idgb - 1], &wa[iwaql - 1],
                 &lwaql, mnn2, &mo1, &kwa[0], &kwa[1], &kwa[2], &kwa[3],
                 &kwa[4], &kwa[5], &kwa[6], &kwa[7], &lkwa0, &wa[i1 - 1],
                 &wa[i2 - 1], &wa[i3 - 1], &wa[i4 - 1], &wa[i5 - 1],
                 &wa[iscf - 1], &wa[i6 - 1], active, &active[m1 - 1]);
	*f = n10nf.f1;
	*ifail = n10nf.ifail1;

	if (*ifail == 2) {

/*		imsl_ermes(4, 1, "The algorithm calculated an uphill search direction.");
*/
                imsl_ermes(IMSL_FATAL, IMSL_UPHILL_DIRECTION);
	} else if (*ifail == 4) {
		imsl_e1sti(1, *maxfun);

/*		imsl_ermes(4, 2, "The line search used  more than %(i1) function calls, therefore it has been declared unsuccessful.");
*/
                imsl_ermes(IMSL_WARNING, IMSL_TOO_MANY_LINESEARCH);
	} else if (*ifail == 1) {

/*		imsl_ermes(3, 3, "Maximum number of iterations exceeded.");
*/
                imsl_ermes(IMSL_WARNING, IMSL_TOO_MANY_ITN);
	} else if (*ifail == 7) {

/*		imsl_ermes(4, 4, "The search direction is close to zero.");
*/
                imsl_ermes(IMSL_FATAL, IMSL_NO_PROGRESS_MADE);
	} else if (*ifail >= 10) {

/*		imsl_ermes(4, 5, "The constraints for the QP subproblem are inconsistent.");
*/
                imsl_ermes(IMSL_FATAL, IMSL_QP_INCONSISTENT);
	}
	if (imsl_n1rcd(0) != 0)
		goto L_9000;

	if (!n10nf.lmer || *ifail != 4)
		goto L_70;
	n10nf.lmer = IMSLFALSE;
	if (n10nf.mode1 == 0 || n10nf.mode1 == 5)
		n10nf.mode1 += 1;
	if (n10nf.mode1 == 2 || n10nf.mode1 == 7)
		n10nf.mode1 += 1;
	*ifail = 0;
	goto L_60;
L_70:
	;
	if ((!rest || *ifail < 1) || n10nf.maxit1 < maxit2)
		goto L_80;
	n10nf.maxit1 -= 10;
	if (n10nf.maxit1 <= 0)
		goto L_80;
	n10nf.lmer = IMSLTRUE;
	if (*iprint == 0)
		goto L_60;
	fprintf(stdout, "\n\n\n\n     ------------------------------------------------------------\n     RESTART\n     ------------------------------------------------------------\n");
	goto L_60;
L_80:
	;
	active[n11nf.lact - 1] = n10nf.lmer;

L_9000:
	imsl_e1pop("N3ONF ");
	return;
}				/* end of function */
#undef DG
/*----------------------------------------------------------------------- */

/*  IMSL Name:  N3ONG/DN3ONG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    December 2, 1985

    Purpose:    Parse the work space for the main driver.

    Usage:      CALL N3ONG (FCNS, GRAD, M, ME, MMAX, N, NMAX, MNN2, X, F,
                            G, DF, DG, LDDG, U, XL, XU, C, LDC, D, ACC,
                            SCBOU, MAXFUN, MAXIT, IPRINT, MODE, IFAIL,
                            WA, LWA, KWA, LKWA, ACTIVE, LACTIV, LLISE,
                            LQL)

    Arguments:
       FCNS   - User-supplied SUBROUTINE to evaluate the functions at
                a given poMint.  The usage is
                CALL FCNS (M, ME, MMAX, N, X, ACTIVE, F, G), where
                M      - Total number of constraMints.  (Input)
                ME     - Number of equality constraMints.  (Input)
                MMAX   - Maximum of (1, M).  (Input)
                N      - Number of variables.  (Input)
                X      - The poMint at which the function is evaluated.
                         (Input)
                         X should not be changed by FCNS.
                ACTIVE - Logical vector of length MMAX indicating the
                         active constraMints.  (Input)
                F      - The computed function value at the poMint X.
                         (Output)
                G      - Vector of length MMAX containing the values of
                         constraMints at poMint X.  (Output)
                FCNS must be declared EXTERNAL in the calling program.
       GRAD   - User-supplied SUBROUTINE to evaluate the gradients at
                a given poMint.  The usage is
                CALL GRAD (M, ME, MMAX, N, X, ACTIVE, F, G, DF, DG LDDG),
                where
                M      - Total number of constraMints.  (Input)
                ME     - Number of equality constraMints.  (Input)
                MMAX   - Maximum of (1, M).  (Input)
                N      - Number of variables.  (Input)
                X      - The poMint at which the function is evaluated.
                         (Input)
                         X should not be changed by FCNS.
                ACTIVE - Logical vector of length MMAX indicating the
                         active constraMints.  (Input)
                F      - The computed function value at the poMint X.
                         (Output)
                G      - Vector of length MMAX containing the values of
                         constraMints at poMint X.  (Output)
                DF     - Vector of lenght N containing the value of the
                         gradient of the objective function.  (Output)
                DG     - MMAX by N array containing the values of the
                         gradients for the active constraMints.  (Output)
                LDDG   - Leading dimension of DG exactly as specified in
                         the dimension statement of the calling program.
                         (Input)
                GRAD must be declared EXTERNAL in the calling program.
       M      - Number of constraMints.  (Input)
       ME     - Number of equality constraMints.  (Input)
       MMAX   - Leading dimension of DG.  (Input)
                MMAX must be at least MAX(1,M).
       N      - Number of variables.  (Input)
       NMAX   - Leading dimension of C.  (Input)
                NMAX must be at least MAX(2,N+1).
       MNN2   - Length of U where MNN2 = M + N + N + 2.  (Input)
       X      - Vector of length N containing the initial guesses to the
                solution on input and the solution on output.
                (Input/Output)
       F      - Scalar containing the objective function value.  (Output)
       G      - Vector of length MMAX containing constraMint values.
                (Output)
       DF     - Vector of length NMAX containing the gradient of the
                of the objective function.  (Output)
       DG     - Array of dimension MMAX by MMAX containing the gradient
                of the constraMints.  (Output)
       LDDG   - Leading dimension of DG exactly as specified in the
                dimension statement of the calling program.  (Input)
       U      - Vector of length MNN2 containing the multipliers of the
                nonlinear constraMints and the bounds.  (Output)
       XL     - Vector of length N containing the lower bounds for the
                variables.  (Input)
       XU     - Vector of length N containing the upper bounds for the
                variables.  (Input)
       C      - Array of dimension NMAX by NMAX containing an the final
                approximation to the Hessian.  (Output)
       LDC    - Leading dimension of C exactly as specified in the
                dimension statement of the calling program.  (Input)
       D      - Vector of length NMAX containing the diagonal elements of
                the Hessian.  (Output)
       ACC    - Final accuracy.  (Input)
       SCBOU  - Scalar containing the scaling variable for the problem
                function.  (Input)
       MAXFUN - Scalar containing the maximum allowable function calls
                during the line search.  (Input)
       MAXIT  - Scalar containing the maximum allowable iterations.
                (Input)
       IPRINT - Specification for the desired print level.  (Input)
       MODE   - Desired solving version for the NLPQL1 algorithm.
                (Input)
       IFAIL  - Scalar containing error message information.  (Output)
       WA     - Work vector of length LWA.
       LWA    - Length of WA where LWA = 2*N*(N+16) + 4*MMAX + 5*M + 66.
                (Input)
       KWA    - Work vector of length LKWA.
       LKWA   - Length of KWA where LKWA = 19 + M.  (Input)
       ACTIVE - Logical vector of length LACTIV containing information
                as to which constraMints are active.  (Output)
       LACTIV - Length of ACTIVE where LACTIV must be at least 200.
                (Input)
       LLISE  - Logical scalar determining if the line search algorithm
                is to be used if LLISE = TRUE.  (Input)
                Otherwise a user supplied line search is used.
       LQL    - Logical scalar determing whether the qudratic programming
                subproblem is to be solved with a positive definite
                Quasi-Newton matrix if LQL = TRUE.  (Input)
                Otherwise a Cholesky-decomposition is performed and
                updated.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
static struct t_n10ng {
	Mfloat           f1, acc1, scbou1, dbdfac, zefac, rpeno, rpens,
	                rpenu, zefacu, delta, imsl_beta, amue, alm;
	Mint             mm1, me1, maxfu1, maxit1, iprin1, mode1, ifail1;
	Mint            llise1, lql1, lmer;
}               n10ng;
static struct t_n11ng {
	Mint             n1, lact, no1, mnn, nmnn;
}               n11ng;
/* LACTIV is not used, but leave the calling sequence intact */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_n3ong(void (*fcns)(Mint, Mint, Mint, Mfloat*, Mint*, Mfloat*,
                    Mfloat*), void (*grad)(Mint, Mint, Mint, Mint, Mfloat*,
                    Mint*, Mfloat, Mfloat*, Mfloat*, Mfloat*), Mint *m,
                    Mint *me, Mint *mmax, Mint *n, Mint *nmax, Mint *mnn2, 
                    Mfloat x[], Mfloat *f, Mfloat g[], Mfloat df[], Mfloat dg[],
                    Mint *lddg, Mfloat u[], Mfloat xl[], Mfloat xu[], Mfloat *c,
                    Mint *ldc, Mfloat d[], Mfloat *acc, Mfloat *scbou,
                    Mint *maxfun, Mint *maxit, Mint *iprint, Mint *mode,
                    Mint *ifail, Mfloat wa[], Mint *lwa, Mint kwa[],
                    Mint *lkwa, Mint active[], Mint *lactiv, Mint *llise,
                    Mint *lql)
#else
static void l_n3ong(void (*fcns)(Mint, Mint, Mint, Mfloat[], Mint[], Mfloat*,
                    Mfloat[]), void (*grad)(Mint, Mint, Mint, Mint, Mfloat[],
                    Mint[], Mfloat, Mfloat[], Mfloat[], Mfloat[]), Mint *m,
                    Mint *me, Mint *mmax, Mint *n, Mint *nmax, Mint *mnn2, 
                    Mfloat x[], Mfloat *f, Mfloat g[], Mfloat df[], Mfloat dg[],
                    Mint *lddg, Mfloat u[], Mfloat xl[], Mfloat xu[], Mfloat *c,
                    Mint *ldc, Mfloat d[], Mfloat *acc, Mfloat *scbou,
                    Mint *maxfun, Mint *maxit, Mint *iprint, Mint *mode,
                    Mint *ifail, Mfloat wa[], Mint *lwa, Mint kwa[],
                    Mint *lkwa, Mint active[], Mint *lactiv, Mint *llise,
                    Mint *lql)
#endif
#else
static void l_n3ong(fcns, grad, m, me, mmax, n, nmax, mnn2, x, f, g, df, dg,
                    lddg, u, xl, xu, c, ldc, d, acc, scbou, maxfun, maxit, 
                    iprint, mode, ifail, wa, lwa, kwa, lkwa, active, lactiv,
                    llise, lql)
	void            (*fcns) (), (*grad) ();
	Mint            *m, *me, *mmax, *n, *nmax, *mnn2;
	Mfloat          x[], *f, g[], df[], dg[];
	Mint            *lddg;
	Mfloat          u[], xl[], xu[], *c;
	Mint            *ldc;
	Mfloat          d[], *acc, *scbou;
	Mint            *maxfun, *maxit, *iprint, *mode, *ifail;
	Mfloat          wa[];
	Mint            *lwa, kwa[], *lkwa;
	Mint            active[];
	Mint            *lactiv;
	Mint            *llise, *lql;
#endif
{
#define DG(I_,J_)	(dg+(I_)*(*lddg)+(J_))
#define C(I_,J_)	(c+(I_)*(*ldc)+(J_))
	Mint            rest;
        FILE            *iout;
	Mint             i, i1, i2, i3, i4, i5, i6, ibdel, iclf, idel, idfb,
	                idgb, idla, idlaol, idphi, ieta, ifb, igb,
	                ipen, iscf, iscg, iv, ivmu, ivmuol, iw,
	                iwaql, ixold, lkwa0, lwaql, m1, maxit2, mo1;
        Mint            IMSLFALSE = 0, IMSLTRUE = 1;


	/*
	 * PROVIDE CONSTANT PARAMETERS FOR CALLING N4ONG
	 */
	imsl_e1psh("N3ONG ");
	n10ng.dbdfac = 2.0e-1;
	n10ng.zefac = F_ONE;
	if (!*lql)
		n10ng.zefac = sqrt(n10ng.zefac);
	n10ng.rpeno = 1.0e1;
	n10ng.rpens = 1.0e-2;
	n10ng.zefacu = 1.0e6;
	if (!*lql)
		n10ng.zefacu = sqrt(n10ng.zefacu);
	n10ng.rpenu = 1.0e9;
	n10ng.delta = 9.0e-1;
	n10ng.imsl_beta = 1.0e-1;
	n10ng.amue = 1.0e-4;
	n10ng.alm = F_ONE;
	imsl_umach(2, &iout);

	/*
	 * PRINT PARAMETERS
	 */
	if (*iprint == 0)
		goto L_10;
	if (*ifail < 0)
		goto L_10;
	fprintf(stdout, "\n\n\n\n   ---------------------------------------------------------------------------\n     START OF THE SEQUENTIAL QUADRATIC PROGRAMMING ALGORITHM\n   ---------------------------------------------------------------------------\n");
	fprintf(stdout, "\n     PARAMETERS:\n        MODE =%3ld\n        ACC =%13.4e\n        SCBOU =%13.4e\n        MAXFUN =%3ld\n        MAXIT =%5ld\n        IPRINT =%4ld\n",
		*mode, *acc, *scbou, *maxfun, *maxit, *iprint);
	if (*iprint != 2)
		goto L_10;
	fprintf(stdout, "\n     OUTPUT IN THE FOLLOWING ORDER:\n        IT    - ITERATION NUMBER\n        F     - OBJECTIVE FUNCTION VALUE\n        SCV   - SUM OF CONSTRAINT VIOLATION");
        fprintf(stdout, "\n        NA    - NUMBER OF ACTIVE CONSTRAINTS\n        I     - NUMBER OF LINE SEARCH ITERATIONS\n        ALPHA - STEPLENGTH PARAMETER\n        DELTA - ADDITIONAL VARIABLE TO PREVENT INCONSISTENCY");
        fprintf(stdout, "\n        DLAN  - MAXIMUM NORM OF LAGRANGIAN GRADIENT\n        KT    - KUHN-TUCKER OPTIMALITY CRITERION\n");
	fprintf(stdout, "\n\n  IT         F          SCV     NA  I    ALPHA    DELTA      DLAN      KT    \n  ---------------------------------------------------------------------------\n");
L_10:
	;
	n10ng.mode1 = *mode;
	n10ng.maxit1 = *maxit;
	maxit2 = *maxit / 2;
	rest = IMSLFALSE;
	if (*mode < 10)
		goto L_20;
	rest = IMSLTRUE;
	n10ng.mode1 = *mode - 10;
L_20:
	;

	/*
	 * INITIAL ADDRESSES IN WA
	 */
	n11ng.n1 = *n + 1;
	n11ng.lact = 2 ** mmax + 6;
	if (*ifail < 0)
		n10ng.lmer = active[n11ng.lact - 1];
	m1 = n11ng.lact + 1;
	n11ng.no1 = *n;
	if (*llise)
		n11ng.no1 = 1;
	mo1 = *mmax;
	if (*llise)
		mo1 = 1;
	n11ng.mnn = *m + *n + *n;
	n11ng.nmnn = *n + n11ng.mnn;
	iscg = 1;
	iscf = iscg + *mmax;
	i1 = iscf + 1;
	i2 = i1 + 1;
	i3 = i2 + 1;
	i4 = i3 + 1;
	i5 = i4 + 1;
	i6 = i5 + 1;
	ivmu = i6 + 1;
	idel = ivmu + n11ng.mnn;
	idla = idel + n11ng.n1;
	iclf = idla + *n;
	ibdel = iclf + n11ng.n1;
	ieta = ibdel + *n;
	ixold = ieta + *n;
	idlaol = ixold + *n;
	iv = idlaol + *n;
	iw = iv + n11ng.n1;
	ivmuol = iw + n11ng.n1;
	idphi = ivmuol + n11ng.mnn;
	ipen = idphi + n11ng.mnn + *n;
	ifb = ipen + n11ng.mnn;
	idfb = ifb + 1;
	igb = idfb + n11ng.no1;
	idgb = igb + *mmax;
	iwaql = idgb + mo1 ** n;
	lwaql = *lwa - iwaql;


	/*
	 * CHECK INITIAL POINT
	 */
	if (*ifail >= 0) {
		if (n10ng.mode1 == 2 || n10ng.mode1 == 3)
			goto L_40;
		if (n10ng.mode1 == 7 || n10ng.mode1 == 8)
			goto L_40;
                *f = n10ng.alm;
L_40:
		for (i = 1; i <= *n; i++) {
			if (x[i - 1] < xl[i - 1])
				x[i - 1] = xl[i - 1];
			if (x[i - 1] > xu[i - 1])
				x[i - 1] = xu[i - 1];
		}
	}
	lkwa0 = *lkwa - 7;

	/*
	 * CALL OF N4ONG
	 */
	if (*ifail >= 0)
		n10ng.lmer = IMSLTRUE;
L_60:
	;
	n10ng.mm1 = *m;
	n10ng.me1 = *me;
	n10ng.f1 = *f;
	n10ng.acc1 = *acc;
	n10ng.scbou1 = *scbou;
	n10ng.maxfu1 = *maxfun;
	n10ng.iprin1 = *iprint;
	n10ng.ifail1 = *ifail;
	n10ng.llise1 = *llise;
	n10ng.lql1 = *lql;
	l_n4ong(fcns, grad, mmax, n, nmax, x, g, df, dg, lddg, u, xl, xu, c,
                ldc, d, &wa[ivmu - 1], &wa[idel - 1], &wa[idla - 1],
                &wa[iclf - 1], &wa[ibdel - 1], &wa[ieta - 1], &wa[ixold - 1],
                &wa[idlaol - 1], &wa[iv - 1], &wa[iw - 1], &wa[ivmuol - 1],
                &wa[idphi - 1], &wa[ipen - 1], &wa[iscg - 1], &wa[ifb - 1],
                &wa[idfb - 1], &wa[igb - 1], &wa[idgb - 1], &wa[iwaql - 1],
                &lwaql, mnn2, &mo1, &kwa[0], &kwa[1], &kwa[2], &kwa[3], &kwa[4],
                &kwa[5], &kwa[6], &kwa[7], &lkwa0, &wa[i1 - 1], &wa[i2 - 1],
                &wa[i3 - 1], &wa[i4 - 1], &wa[i5 - 1], &wa[iscf - 1],
		&wa[i6 - 1], active, &active[m1 - 1]);
	*f = n10ng.f1;
	*ifail = n10ng.ifail1;

	if (*ifail == 2) {
/*		imsl_ermes(4, 1, "The algorithm calculated an uphill search direction."); */
                imsl_ermes(IMSL_FATAL, IMSL_UPHILL_DIRECTION);
	} else if (*ifail == 4) {
		imsl_e1sti(1, *maxfun);
/*		imsl_ermes(4, 2, "The line search used more than %(i1) function calls, therefore it has been declared unsuccessful."); */
                imsl_ermes(IMSL_FATAL, IMSL_TOO_MANY_LINESEARCH);

	} else if (*ifail == 1) {
/*		imsl_ermes(3, 3, "Maximum number of iterations exceeded."); */
                imsl_ermes(IMSL_WARNING, IMSL_TOO_MANY_ITN);

	} else if (*ifail == 7) {
/*		imsl_ermes(4, 4, "The search direction is close to zero."); */
          imsl_ermes(IMSL_FATAL, IMSL_NO_PROGRESS_MADE);

	} else if (*ifail >= 10) {
/*		imsl_ermes(4, 5, "The constraints for the QP subproblem are inconsistent."); */
                imsl_ermes(IMSL_FATAL, IMSL_QP_INCONSISTENT);
	}
	if (imsl_n1rcd(0) != 0)
		goto L_9000;

	if (!n10ng.lmer || *ifail != 4)
		goto L_70;
	n10ng.lmer = IMSLFALSE;
	if (n10ng.mode1 == 0 || n10ng.mode1 == 5)
		n10ng.mode1 += 1;
	if (n10ng.mode1 == 2 || n10ng.mode1 == 7)
		n10ng.mode1 += 1;
	*ifail = 0;
	goto L_60;
L_70:
	;
	if ((!rest || *ifail < 1) || n10ng.maxit1 < maxit2)
		goto L_80;
	n10ng.maxit1 -= 10;
	if (n10ng.maxit1 <= 0)
		goto L_80;
	n10ng.lmer = IMSLTRUE;
	if (*iprint == 0)
		goto L_60;
	fprintf(stdout,"\n\n\n\n     ------------------------------------------------------------\n     RESTART\n     ------------------------------------------------------------\n");
	goto L_60;
L_80:
	;
	active[n11ng.lact - 1] = n10ng.lmer;

	/*
	 * END OF NLPQL1
	 */
L_9000:
	imsl_e1pop("N3ONG ");
	return;
}				/* end of function */
#undef DG
#undef C
/*----------------------------------------------------------------------- */

/*  IMSL Name:  N4ONF/DN4ONF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    December 2, 1985

    Purpose:    Main driver for the successive quadratic programming
                algorithm.

    Usage:      CALL N4ONF (FCNS, MMAX, N, NMAX, X, XS, G, DF, DG, LDDG,
                            U, XL, XU, DCL, LDDCL, CD, CWK, VMU, DEL,
                            DLA, DCLF, BDEL, ETA, XOLD, DLAOLD, V, W,
                            VMUOLD, DPHI, RPEN, SCG, FBEST, DFBEST,
                            GBEST, DGBEST, WA, LWA, MNN2,
                            MO1, NFUNC, NGRAD, ITER, NQL, ILINE,
                            IFLISE, NOPT, IW, LIW, PHI, DFDEL, DBD,
                            ALPHAM, ALPHAO, SCF, PRD, ACTIVE, L7)

    Arguments:
       FCNS   - User-supplied SUBROUTINE to evaluate the functions at
                a given poMint.  The usage is
                CALL FCNS (M, ME, N, X, ACTIVE, F, G), where
                M      - Total number of constraMints.  (Input)
                ME     - Number of equality constraMints.  (Input)
                N      - Number of variables.  (Input)
                X      - The poMint at which the function is evaluated.
                         (Input)
                         X should not be changed by FCNS.
                ACTIVE - Logical vector of length MMAX indicating the
                         active constraMints.  (Input)
                F      - The computed function value at the poMint X.
                         (Output)
                G      - Vector of length MMAX containing the values of
                         constraMints at poMint X.  (Output)
                FCNS must be declared EXTERNAL in the calling program.
       MMAX   - Order of the array DG.  (Input)
                MMAX must be at least MAX(1,M).
       N      - Number of variables.  (Input)
       NMAX   - Order of DCL where NMAX must be at least MAX(2,N+1).
                (Input)
       X      - Vector of length N containing the initial guesses to the
                solution on input and the solution on output.
                (Input/Output)
       XS     - Vector of length N containing the diagonal scaling
                matrix.  (Input)
       G      - Vector of length MMAX containing constraMint values.
                (Output)
       DF     - Vector of length N+1 containing the gradient of the
                of the objective function.  (Output)
       DG     - Array of dimension MMAX by MMAX containing the gradient
                of the constraMints.  (Output)
       LDDG   - Leading dimension of DG exactly as specified in the
                dimension statement in the calling program.  (Input)
       U      - Vector of length MNN2 containing the multipliers of the
                nonlinear constraMints and the bounds.  (Output)
       XL     - Vector of length N containing the lower bounds for the
                variables.  (Input)
       XU     - Vector of length N containing the upper bounds for the
                variables.  (Input)
       DCL    - Array of dimension NMAX by NMAX containing an the final
                approximation to the Hessian.  (Output)
       LDDCL  - Leading dimension of DCL exactly as specified in the
                dimension statement in the calling program.  (Input)
       CD     - Vector of length NMAX containing the diagonal elements of
                the Hessian.  (Output)
       CWK    - Work vector of length M used in gradient evaluation.
       VMU    - Work vector of length M + 2*N.
       DEL    - Work vector of length N + 1.
       DLA    - Work vector of length N.
       DCLF   - Work vector of length N + 1.
       BDEL   - Work vector of length N.
       ETA    - Work vector of length N.
       XOLD   - Work vector of length N.
       DLAOLD - Work vector of length N.
       V      - Work vector of length N + 1.
       W      - Work vector of length N + 1.
       VMUOLD - Work vector of length M + 2*N.
       DPHI   - Work vector of length M + 3*N.
       RPEN   - Work vector of length M + 2*N.
       SCG    - Work vector of length MMAX.
       FBEST  - Work scalar.
       DFBEST - Work vector of length NMAX.
       GBEST  - Work vector of length MMAX.
       DGBEST - Work array of dimension MO1 by N.
       WA     - Work vector of length LWA.
       LWA    - Length of WA where LWA = N*(2*N+13) + M + MMAX + 12.
                (Input)
       N1     - Scalar containing the value N + 1.  (Input)
       MNN    - Scalar containing the value M + 2*N.  (Input)
       MNN2   - Scalar containing the value M + 2*N + 2.  (Input)
       NMNN   - Scalar containing the value M + 3*N.  (Input)
       NO1    - Scalar containing the value 1 when LLISE is true or N
                when LLISE is false.  (Input)
       MO1    - Scalar containing the value 1 when LLISE is true or MMAX
                when LLISE is false.  (Input)
       NFUNC  - Number of function evaluations.  (Output)
       NGRAD  - Number of gradient evaluations.  (Output)
       ITER   - Number of iterations.  (Output)
       NQL    - Number of QL algorithm evaluations.  (Output)
       ILINE  - Number of line search evaluations.  (Output)
       IFLISE - Error parameter for line search algorithm.  (Output)
       NOPT   - Number of optimality iterations.  (Output)
       IW     - Work vector of length LIW.
       LIW    - Length of IW where LIW = 12.  (Input)
       PHI    - Scalar variable.
       DFDEL  - Scalar variable.
       DBD    - Scalar variable.
       ALPHAM - Scalar variable.
       ALPHAO - Scalar variable.
       SCF    - Scalar variable.
       PRD    - Scalar variable.
       ACTIVE - Logical vector of length LACTIV indicating which
                constraMints are active.  (Output)
       LACTIV - Length of ACTIVE where LACTIV must be at least 200.
                (Input)
       L7     - Logical vector of length 7.

    Remark:
       The NLPQL algorithm was designed by K. Schittkowski.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_n4onf(void (*fcns)(Mint, Mint, Mint, Mfloat*, Mint*, Mfloat*,
                    Mfloat*), Mint *mmax, Mint *n, Mint *nmax, Mfloat x[],
                    Mfloat xs[], Mfloat g[], Mfloat df[], Mfloat dg[],
                    Mint *lddg, Mfloat u[], Mfloat xl[], Mfloat xu[],
                    Mfloat *dcl, Mint *lddcl, Mfloat cd[], Mfloat cwk[],
                    Mfloat vmu[], Mfloat del[], Mfloat dla[], Mfloat dclf[],
                    Mfloat bdel[], Mfloat eta[], Mfloat xold[], Mfloat dlaold[],
                    Mfloat v[], Mfloat w[], Mfloat vmuold[], Mfloat dphi[],
                    Mfloat rpen[], Mfloat scg[], Mfloat *fbest, Mfloat dfbest[],
                    Mfloat gbest[], Mfloat *dgbest, Mfloat wa[], Mint *lwa, 
                    Mint *mnn2, Mint *mo1, Mint *nfunc, Mint *ngrad, Mint *iter,
                    Mint *nql, Mint *iline, Mint *iflise, Mint *nopt, Mint iw[],
                    Mint *liw, Mfloat *phi, Mfloat *dfdel, Mfloat *dbd, 
                    Mfloat *alpham, Mfloat *alphao, Mfloat *scf, Mfloat *prd,
                    Mint active[], Mint l7[])
#else
static void l_n4onf(void (*fcns)(Mint, Mint, Mint, Mfloat[], Mint[], Mfloat*,
                    Mfloat[]), Mint *mmax, Mint *n, Mint *nmax, Mfloat x[],
                    Mfloat xs[], Mfloat g[], Mfloat df[], Mfloat dg[],
                    Mint *lddg, Mfloat u[], Mfloat xl[], Mfloat xu[],
                    Mfloat *dcl, Mint *lddcl, Mfloat cd[], Mfloat cwk[],
                    Mfloat vmu[], Mfloat del[], Mfloat dla[], Mfloat dclf[],
                    Mfloat bdel[], Mfloat eta[], Mfloat xold[], Mfloat dlaold[],
                    Mfloat v[], Mfloat w[], Mfloat vmuold[], Mfloat dphi[],
                    Mfloat rpen[], Mfloat scg[], Mfloat *fbest, Mfloat dfbest[],
                    Mfloat gbest[], Mfloat *dgbest, Mfloat wa[], Mint *lwa, 
                    Mint *mnn2, Mint *mo1, Mint *nfunc, Mint *ngrad, Mint *iter,
                    Mint *nql, Mint *iline, Mint *iflise, Mint *nopt, Mint iw[],
                    Mint *liw, Mfloat *phi, Mfloat *dfdel, Mfloat *dbd, 
                    Mfloat *alpham, Mfloat *alphao, Mfloat *scf, Mfloat *prd,
                    Mint active[], Mint l7[])
#endif
#else
static void l_n4onf(fcns, mmax, n, nmax, x, xs, g, df, dg, lddg, u, xl, xu, dcl,
                    lddcl, cd, cwk, vmu, del, dla, dclf, bdel, eta, xold,
                    dlaold, v, w, vmuold, dphi, rpen, scg, fbest, dfbest, gbest,
                    dgbest, wa, lwa, mnn2, mo1, nfunc, ngrad, iter, nql, iline,
                    iflise, nopt, iw, liw, phi, dfdel, dbd, alpham, alphao, scf,
                    prd, active, l7)
	void            (*fcns) ();
	Mint            *mmax, *n, *nmax;
	Mfloat           x[], xs[], g[], df[], dg[];
	Mint            *lddg;
	Mfloat           u[], xl[], xu[], *dcl;
	Mint            *lddcl;
	Mfloat           cd[], cwk[], vmu[], del[], dla[], dclf[], bdel[],
	                eta[], xold[], dlaold[], v[], w[], vmuold[], dphi[],
	                rpen[], scg[], *fbest, dfbest[], gbest[], *dgbest,
	                wa[];
	Mint            *lwa, *mnn2, *mo1, *nfunc, *ngrad, *iter, *nql,
	               *iline, *iflise, *nopt, iw[], *liw;
	Mfloat          *phi, *dfdel, *dbd, *alpham, *alphao, *scf, *prd;
	Mint            active[], l7[];
#endif
{
#define DG(I_,J_)	(dg+(I_)*(*lddg)+(J_))
#define DCL(I_,J_)	(dcl+(I_)*(*lddcl)+(J_))
#define DGBEST(I_,J_)	(dgbest+(I_)*(*mo1)+(J_))
        FILE            *iout;
	Mint             _l0, _l1, _l2, i, ifail1, ilwls, imerit,
	                ipr, irpmax, j, liwql, lwaql, 
	                mmax2, mn, mn1, n2, nact;
	Mfloat           dbd1, dbdi, dcl11, delnm, dlan, edel,
	                edeli, eps0, fact, ff, of, on, rpmax, sdcl11,
	                sqacc, sqd, sres, sum, theta, theta1, tw, uad,
	                uf, xnm, ze;
        Mint            IMSLFALSE = 0, IMSLTRUE = 1;
	static struct {
		Mfloat           f, acc, scbou, dbdfac, zefac, rpeno, rpens,
		                rpenu, zefacu, delta, imsl_beta, amue,
		                alm;
		Mint             m, me, maxfun, maxit, iprint, mode, ifail;
		Mint            llise, lql, lmerit;
	}              *_n10nf = (void *) &n10nf;
	static struct {
		Mint             n1, lact, no1, mnn, nmnn;
	}              *_n11nf = (void *) &n11nf;

	/* CONSTANT DATA */
	ze = F_ZERO;
	on = F_ONE;
	tw = F_TWO;
	eps0 = 100.0e0 * imsl_amach(4);
	uf = eps0 * eps0;
	of = on / uf;
	imsl_umach(2, &iout);
	/* INITIAL DEFINITIONS */
	mn = _n10nf->m + *n;
	n2 = *n + *n;
	lwaql = *lwa - *mmax - 40;
	liwql = *liw - 10;
	ilwls = 2 ** mmax + 1;
	imerit = 0;
	if (!_n10nf->lmerit)
		imerit = 4;
	l7[5] = IMSLFALSE;
	l7[3] = IMSLFALSE;
	l7[4] = IMSLFALSE;
	sqacc = sqrt(_n10nf->acc);
	if (((_n10nf->mode == 2 || _n10nf->mode == 7) || _n10nf->mode ==
	     3) || _n10nf->mode == 8) {
		l7[5] = IMSLTRUE;
		if (_n10nf->ifail == -1)
			goto L_610;
		if (_n10nf->ifail == -2)
			goto L_650;
	}
	*iline = 0;
	*alphao = ze;
	*nfunc = 0;
	*ngrad = 0;
	*iter = 0;
	*nql = 0;
	*nopt = 0;
	if (_n10nf->m != 0) {
		mmax2 = *mmax + *mmax;
		for (j = 1; j <= mmax2; j++) {
			active[j - 1] = IMSLTRUE;
		}
	}
	if (!l7[5]) {
		imsl_e1usr("ON");
		(*fcns) (_n10nf->m, _n10nf->me, *n, x, &active[*mmax],
                         &_n10nf->f, g);
		imsl_e1usr("OFF");
		l_n5onf(fcns, &_n10nf->m, &_n10nf->me, mmax, n, x, xs,
                        active, &_n10nf->f, g, df, dg, cwk);
	}
	l7[0] = IMSLFALSE;
	l7[1] = IMSLFALSE;
	if (fabs(_n10nf->f) >= _n10nf->scbou) {
		l7[0] = IMSLTRUE;
		if (_n10nf->scbou > ze)
			*scf = F_ONE / sqrt(fabs(_n10nf->f));
		_n10nf->f *= *scf;
		sscal(*n, *scf, df, 1);
	}
	if (_n10nf->m != 0) {
		for (j = 1; j <= _n10nf->m; j++) {
			if (fabs(g[j - 1]) >= _n10nf->scbou)
				l7[1] = IMSLTRUE;
		}
	}
	if (l7[1]) {
		for (j = 1; j <= _n10nf->m; j++) {
			if (_n10nf->scbou > ze)
				scg[j-1] = F_ONE / imsl_f_max(F_ONE,
                                                            sqrt(fabs(g[j-1])));
			g[j - 1] *= scg[j - 1];
			sscal(*n, scg[j - 1], DG(0, j - 1), *lddg);
		}
	}
	if (_n10nf->iprint >= 1) {
		if (l7[0] && !l7[1]) {
			fprintf(stdout, "\n     OBJECTIVE FUNCTION WILL BE SCALED\n");
		}
		if (l7[0] && l7[1]) {
			fprintf(stdout, "\n     OBJECTIVE AND CONSTRAINT FUNCTIONS WILL BE SCALED\n");
		}
		if (!l7[0] && l7[1]) {
			fprintf(stdout, "\n     CONSTRAINT FUNCTIONS WILL BE SCALED\n");
		}
	}
	*nfunc += 1;
	*ngrad += 1;
	dclf[n11nf.n1 - 1] = F_ZERO;
	scopy(*n, df, 1, del, 1);
	sscal(*n, -F_ONE, del, 1);
	sset(*n, F_ZERO, DCL(0, n11nf.n1 - 1), *lddcl);
	sset(*n, F_ZERO, DCL(n11nf.n1 - 1, 0), 1);
	*DCL(n11nf.n1 - 1, n11nf.n1 - 1) = _n10nf->zefac;
	if (((_n10nf->mode == 1 || _n10nf->mode == 6) || _n10nf->mode ==
	      3) || _n10nf->mode == 8) {
		if (_n10nf->lql)
			goto L_50;
		goto L_750;
	}
	sset(*n, F_ONE, cd, 1);
	for (i = 1; i <= *n; i++) {
		sset(*n, F_ZERO, DCL(i - 1, 0), 1);
	}
	sset(*n, F_ONE, DCL(0, 0), *lddcl + 1);
L_50:
	sset(n11nf.mnn, _n10nf->rpens, rpen, 1);
	sset(n11nf.mnn, F_ZERO, vmu, 1);
	if (((_n10nf->mode == 1 || _n10nf->mode == 6) || _n10nf->mode ==
	      3) || _n10nf->mode == 8)
		scopy(n11nf.mnn, u, 1, vmu, 1);
        _l0 = imerit + 3;
        _l1 = 4;
	l_n5ong(&_l0, &_n10nf->m, &_n10nf->me, n, &n11nf.mnn, &n11nf.nmnn,
                &_n10nf->acc, rpen, &_n10nf->f, df, g, dg, lddg, vmu, u, x,
                xl, xu, phi, dphi, active, wa, &_l1);
	/*
	 * START MAIN LOOP, PRINT INTERMEDIATE ITERATES
	 */
L_60:
	;
	l7[2] = IMSLFALSE;
	if (_n10nf->iprint < 3)
		goto L_90;
	if (l7[0])
		_n10nf->f /= *scf;
	fprintf(stdout, "\n\n     ITERATION%3ld\n\n        FUNCTION VALUE:  F(X) =%16.8e",
		*iter, _n10nf->f);
	fprintf(stdout, "\n        VARIABLE:  X =\n         ");
	for (i = 1; i <= *n; i++) {
		fprintf(stdout, "%16.8e", x[i - 1]);
	}
	fprintf(stdout, "\n");
	if (l7[0])
		_n10nf->f *= *scf;
	if (_n10nf->m != 0 && (l7[0] || l7[1])) {
		if (l7[0])
			sscal(_n10nf->m, F_ONE / *scf, vmu, 1);
		if (l7[1]) {
			for (j = 1; j <= _n10nf->m; j++) {
				vmu[j - 1] *= scg[j - 1];
				g[j - 1] /= scg[j - 1];
			}
		}
	}
	fprintf(stdout, "        MULTIPLIERS:  U =\n         ");
	for (j = 1; j <= n11nf.mnn; j++) {
		fprintf(stdout, "%16.8e", vmu[j - 1]);
	}
	fprintf(stdout, "\n");
	if (_n10nf->m != 0) {
		fprintf(stdout, "        CONSTRAINTS: G(X) =\n         ");
		for (j = 1; j <= _n10nf->m; j++) {
			fprintf(stdout, "%16.8e", g[j - 1]);
		}
		fprintf(stdout, "\n");
		if (l7[0] || l7[1]) {
			if (l7[0])
				sscal(_n10nf->m, *scf, vmu, 1);
			if (l7[1]) {
				for (j = 1; j <= _n10nf->m; j++) {
					vmu[j - 1] /= scg[j - 1];
					g[j - 1] *= scg[j - 1];
				}
			}
		}
	}
L_90:
	*iter += 1;
	if (*iter < _n10nf->maxit)
		goto L_100;
	_n10nf->ifail = 1;
	if (_n10nf->iprint == 0)
		goto L_350;
	fprintf(stdout, "        **MORE THAN MAXIT ITERATIONS\n");
	goto L_350;
L_100:
	;
	/* SEARCH DIRECTION */
	scopy(*n, df, 1, dclf, 1);
	for (i = 1; i <= *n; i++) {
		v[i - 1] = xl[i - 1] - x[i - 1];
	}
	for (i = 1; i <= *n; i++) {
		w[i - 1] = xu[i - 1] - x[i - 1];
	}
	ipr = 0;
	if (_n10nf->iprint > 10 && _n10nf->iprint < 1000)
		ipr = _n10nf->iprint - 10;
	if (_n10nf->mode >= 5)
		goto L_130;
	ifail1 = *iter;
	if (l7[3] || l7[4])
		ifail1 = 1;
	iw[10] = 0;
	if (_n10nf->lql)
		iw[10] = 1;
	iw[11] = 0;
	l_n6ong(&_n10nf->m, &_n10nf->me, mmax, n, nmax, &n11nf.mnn, dcl, lddcl,
                dclf, dg, lddg, g, v, w, del, u, &ifail1, &ipr, &wa[*mmax + 40],
		&lwaql, &iw[10], &liwql);
	del[n11nf.n1 - 1] = ze;
	*nql += 1;
	l7[3] = IMSLFALSE;
	if (ifail1 == 0)
		goto L_220;
L_130:
	;
	if (*iter == 1)
		goto L_140;
	fact = tw * fabs(*dbd ** dfdel) / (sqrt(*dbd) * (on - del[n11nf.n1 - 1]));
	if (_n10nf->lql)
		fact *= fact;
	dcl11 = imsl_f_max(_n10nf->zefac, fact);
	*DCL(n11nf.n1 - 1, n11nf.n1 - 1) = imsl_f_min(_n10nf->zefacu, dcl11);
L_140:
	;
	sset(*n, F_ZERO, del, 1);
	del[n11nf.n1 - 1] = F_ONE;

	if (_n10nf->m != 0) {
		scopy(_n10nf->m, g, 1, DG(n11nf.n1 - 1, 0), 1);
		sscal(_n10nf->m, -F_ONE, DG(n11nf.n1 - 1, 0), 1);
		for (j = 1; j <= _n10nf->m; j++) {
			if (!active[j - 1])
				*DG(n11nf.n1 - 1, j - 1) = F_ZERO;
		}
	}
	v[n11nf.n1 - 1] = F_ZERO;
	w[n11nf.n1 - 1] = F_ONE;
	ifail1 = -*iter;
	if (!l7[3] || l7[4])
		ifail1 = -1;
	iw[10] = 0;
	if (_n10nf->lql)
		iw[10] = 1;
	iw[11] = 1;
	l_n6ong(&_n10nf->m, &_n10nf->me, mmax, &n11nf.n1, nmax, mnn2, dcl,
                lddcl, dclf, dg, lddg, g, v, w, del, u, &ifail1, &ipr,
                &wa[*mmax + 40], &lwaql, &iw[10], &liwql);
	*nql += 1;
	mn1 = _n10nf->m + n11nf.n1 + 1;
	l7[3] = IMSLTRUE;
	if (ifail1 == 0)
		goto L_170;
L_160:
	_n10nf->ifail = 10 + ifail1;
	if (_n10nf->iprint == 0)
		goto L_350;
	fprintf(stdout, "        **ERROR IN QL. IFAIL(QL) =%3ld\n",
		ifail1);
	goto L_350;
L_170:
	;
	scopy(*n + 1, &u[mn1 - 1], 1, &u[mn1 - 2], 1);
	if (_n10nf->iprint < 3)
		goto L_180;
	fprintf(stdout, "        ADDITIONAL VARIABLE TO PREVENT INCONSISTENCY:  DELTA =%13.4e\n",
		del[n11nf.n1 - 1]);
	sdcl11 = *DCL(n11nf.n1 - 1, n11nf.n1 - 1);
	if (!_n10nf->lql)
		sdcl11 = sqrt(sdcl11);
	fprintf(stdout, "        PENALTY PARAMETER FOR DELTA:  RHO =%13.4e\n",
		sdcl11);
L_180:
	;
	dcl11 = *DCL(n11nf.n1 - 1, n11nf.n1 - 1);
	if (del[n11nf.n1 - 1] < _n10nf->delta)
		goto L_220;
	*DCL(n11nf.n1 - 1, n11nf.n1 - 1) = dcl11 * _n10nf->rpeno;
	if (_n10nf->lql)
		*DCL(n11nf.n1 - 1, n11nf.n1 - 1) *= _n10nf->rpeno;
	if (dcl11 < _n10nf->zefacu)
		goto L_140;
	/*
	 * AUGMENTED LAGRANGIAN TYPE SEARCH DIRECTION
	 */
L_190:
	l7[4] = IMSLTRUE;
	if (_n10nf->iprint < 3)
		goto L_200;
	fprintf(stdout, "        **WARNING: AUGMENTED LAGRANGIAN SEARCH DIRECTION\n");
L_200:
        _l0 = 4;
	l_n5ong(&_l0, &_n10nf->m, &_n10nf->me, n, &n11nf.mnn, &n11nf.nmnn,
	        &_n10nf->acc, rpen, &_n10nf->f, df, g, dg, lddg, vmu, u, x, xl,
		xu, phi, dphi, active, wa, &_l0);
	scopy(*n, dphi, 1, &wa[40], 1);
	scopy(*n, dphi, 1, dclf, 1);
	ifail1 = 1;
	iw[10] = 0;
	if (_n10nf->lql)
		iw[10] = 1;
	iw[11] = 0;
        _l0 = 0;
	l_n6ong(&_l0, &_l0, mmax, n, nmax, mnn2, dcl, lddcl, dclf, dg, lddg,
                g, v, w, del, u, &ifail1, &ipr, &wa[*mmax + 40], &lwaql,
                &iw[10], &liwql);
	if (ifail1 > 0)
		goto L_160;
	if (_n10nf->m == 0)
		goto L_230;
	scopy(n2, u, -1, &u[_n10nf->m], -1);

	for (j = 1; j <= _n10nf->m; j++) {
		u[j - 1] = vmu[j - 1] - dphi[*n + j - 1];
	}
	goto L_230;
L_220:
	l7[4] = IMSLFALSE;
	/*
	 * PROJECTION OF DEL, MAXIMAL STEPLENGTH, AND NORM OF X,DEL
	 */
L_230:
	*alpham = of;
	xnm = F_ZERO;
	delnm = F_ZERO;
	for (i = 1; i <= *n; i++) {
		if (w[i - 1] < del[i - 1])
			del[i - 1] = w[i - 1];
		if (v[i - 1] > del[i - 1])
			del[i - 1] = v[i - 1];
		uad = fabs(del[i - 1]);
		if (del[i - 1] > uf)
			*alpham = imsl_f_min(*alpham, w[i - 1] / del[i - 1]);
		if (del[i - 1] < -uf)
			*alpham = imsl_f_min(*alpham, v[i - 1] / del[i - 1]);
		xnm = imsl_f_max(fabs(x[i - 1]), xnm);
		delnm = imsl_f_max(uad, delnm);
	}
	*alpham = imsl_f_max(on, *alpham);
	*alpham = imsl_f_min(*alpham, _n10nf->alm);
	/* GRADIENT OF LAGRANGIAN */
L_250:
	for (i = 1; i <= *n; i++) {
		uad = df[i - 1];
		if (l7[4])
			uad = dphi[i - 1];
		dla[i - 1] = uad - u[_n10nf->m + i - 1] - u[mn + i - 1];
	}

	if (_n10nf->m != 0 && !l7[4]) {
		for (j = 1; j <= _n10nf->m; j++) {
			if (u[j - 1] != ze)
				saxpy(*n, -u[j - 1], DG(0, j - 1), *lddg, dla, 1);
		}
	}
	if (l7[2])
		goto L_680;
	/* STORE SOME DATA */
	scopy(*n, dla, 1, dlaold, 1);
	scopy(*n, x, 1, xold, 1);
	*dfdel = imsl_sdot(*n, df, 1, del, 1);
	irpmax = imsl_isamax(*n, dla, 1);
	dlan = fabs(dla[irpmax - 1]);
	scopy(n11nf.mnn, vmu, 1, vmuold, 1);
	/* DETERMINE B*D AND D(T)*B*D */
	if (!_n10nf->lql) {
		for (i = 1; i <= *n; i++) {
			eta[i - 1] = imsl_sdot(*n - i + 1, DCL(i - 1, i - 1),
                                               *lddcl, &del[i-1], 1);
		}
		*dbd = imsl_sdot(*n, eta, 1, eta, 1);
		for (i = 1; i <= *n; i++) {
			bdel[i - 1] = imsl_sdot(i, DCL(i - 1, 0), 1, eta, 1);
		}
	} else {
		for (i = 1; i <= *n; i++) {
			bdel[i - 1] = imsl_sdot(*n, DCL(0, i - 1), *lddcl,
                                                del, 1);
		}
		*dbd = imsl_sdot(*n, bdel, 1, del, 1);
	}
	/* TEST FOR OPTIMALITY AND FINAL OUTPUT */
	sres = ze;
	sum = fabs(*dfdel);
	if (l7[0])
		sum /= *scf;
	nact = 0;
	if (_n10nf->m != 0) {
		for (j = 1; j <= _n10nf->m; j++) {
			if (active[j - 1])
				nact += 1;
			uad = fabs(g[j - 1]);
			if (l7[1])
				uad /= scg[j - 1];
			if (j <= _n10nf->me || g[j - 1] < ze)
				sres += uad;
		}
                _l0 = 1;
		sum += l_a1ot(_n10nf->m, u, _l0, g, _l0);
		if (_n10nf->iprint == 3) {
			fprintf(stdout, "        SUM OF CONSTRAINT VIOLATIONS:                    SCV =%13.4e\n",
				sres);
			fprintf(stdout, "        NUMBER OF ACTIVE CONSTRAINTS:                    NAC =%4ld\n",
				nact);
		}
	}
	for (i = 1; i <= *n; i++) {
		sum += fabs(u[_n10nf->m + i - 1] * v[i - 1]) + fabs(u[mn + i - 1] *
								  w[i - 1]);
	}
	if (_n10nf->iprint == 2) {
		ff = _n10nf->f;
		if (l7[0])
			ff = _n10nf->f / *scf;
		fprintf(stdout, " %3ld%16.8e%10.2e%4ld%3ld%10.2e%10.2e%10.2e%10.2e\n",
		*iter, ff, sres, nact, *iline, *alphao, del[n11nf.n1 - 1],
			dlan, sum);
	}
	if (_n10nf->iprint == 3) {
		fprintf(stdout, "        KUHN-TUCKER OPTIMALITY CONDITION:                KTO =%13e\n",
			sum);
		fprintf(stdout, "        NORM OF LAGRANGIAN GRADIENT:                     NLG =%13e\n",
			dlan);
	}
	if (*dbd >= uf)
		goto L_330;
	if (sres < sqacc)
		goto L_340;
	if (*dbd > ze)
		goto L_390;
	if (!l7[4])
		goto L_190;
	_n10nf->ifail = 7;
	if (_n10nf->iprint == 0)
		goto L_350;
	fprintf(stdout, "        **UNDERFLOW IN D(T)*B*D AND INFEASIBLE ITERATE X\n");
	goto L_350;
L_330:
	;
	if (sum >= _n10nf->acc || sres > sqacc)
		goto L_390;
	if (dlan <= sqrt(sqacc) || *dbd <= _n10nf->acc)
		goto L_340;
	*nopt += 1;
	if (*nopt < 3)
		goto L_390;
L_340:
	_n10nf->ifail = 0;
L_350:
	;
	if (l7[0])
		_n10nf->f /= *scf;
	if (_n10nf->m == 0 || (!l7[0] && !l7[1]))
		goto L_370;
	if (l7[0])
		sscal(*n, F_ONE / *scf, u, 1);
	if (l7[1]) {
		for (j = 1; j <= _n10nf->m; j++) {
			u[j - 1] *= scg[j - 1];
			g[j - 1] /= scg[j - 1];
		}
	}
L_370:
	;

	if (_n10nf->iprint == 0)
		goto L_9000;
	fprintf(stdout, "\n\n     * FINAL CONVERGENCE ANALYSIS\n\n");
	fprintf(stdout, "        OBJECTIVE FUNCTION VALUE:  F(X) =%16.8e\n",
		_n10nf->f);
	fprintf(stdout, "        APPROXIMATION OF SOLUTION:  X =\n         ");
	for (i = 1; i <= *n; i++) {
		fprintf(stdout, "%16.8e", x[i - 1]);
	}
	fprintf(stdout, "\n");
	fprintf(stdout, "        APPROXIMATION OF MULTIPLIERS:  U =\n         ");
	for (j = 1; j <= n11nf.mnn; j++) {
		fprintf(stdout, "%16.8e", u[j - 1]);
	}
	fprintf(stdout, "\n");
	if (_n10nf->m == 0)
		goto L_380;
	fprintf(stdout, "        CONSTRAINT VALUES:  G(X) =\n         ");
	for (j = 1; j <= _n10nf->m; j++) {
		fprintf(stdout, "%16.8e", g[j - 1]);
	}
	fprintf(stdout, "\n");
L_380:
	fprintf(stdout, "        DISTANCE FROM LOWER BOUND:  XL-X =\n         ");
	for (i = 1; i <= *n; i++) {
		fprintf(stdout, "%16.8e", v[i - 1]);
	}
	fprintf(stdout, "\n");
	fprintf(stdout, "        DISTANCE FROM UPPER BOUND:  XU-X =\n         ");
	for (i = 1; i <= *n; i++) {
		fprintf(stdout, "%16.8e", w[i - 1]);
	}
	fprintf(stdout, "\n");
	if (!_n10nf->llise) {
		fprintf(stdout, "        NUMBER OF ITERATIONS:  ITER =%4ld\n",
			*iter);
	}
	fprintf(stdout, "        NUMBER OF FUNC-CALLS:  NFUNC =%4ld\n",
		*nfunc);
	fprintf(stdout, "        NUMBER OF GRAD-CALLS:  NGRAD =%4ld\n",
		*ngrad);
	fprintf(stdout, "        NUMBER OF QL-CALLS:    NQL   =%4ld\n\n\n\n",
		*nql);
	goto L_9000;
L_390:
	;
	/* CORRECT PENALTY PARAMETER */
	if (l7[4])
		goto L_400;
	wa[0] = *dbd;
	wa[1] = del[n11nf.n1 - 1];
	wa[2] = _n10nf->rpenu;
	wa[3] = (Mfloat) (*iter);
        _l0 = imerit + 2;
        _l1 = 4;
	l_n5ong(&_l0, &_n10nf->m, &_n10nf->me, n, &n11nf.mnn, &n11nf.nmnn,
                &_n10nf->acc, rpen, &_n10nf->f, df, g, dg, lddg, vmu, u, x,
                xl, xu, phi, dphi, active, wa, &_l1);
	goto L_430;
L_400:
	sum = ze;
	for (i = 1; i <= *n; i++) {
		sum += dphi[i - 1] * del[i - 1] + fabs(u[_n10nf->m + i - 1] *
				 v[i - 1]) + fabs(u[mn + i - 1] * w[i - 1]);
	}
	if (sum > sqrt(sqacc))
		goto L_430;
	for (j = 1; j <= n11nf.mnn; j++) {
		rpen[j - 1] = imsl_f_min(_n10nf->zefacu, rpen[j - 1] * _n10nf->rpeno);
	}
        _l0 = imerit + 4;
        _l1 = 4;
	l_n5ong(&_l0, &_n10nf->m, &_n10nf->me, n, &n11nf.mnn, &n11nf.nmnn,
                &_n10nf->acc, rpen, &_n10nf->f, df, g, dg, lddg, vmu, u, x,
                xl, xu, phi, dphi, active, wa, &_l1);
L_430:
	if (_n10nf->iprint < 3)
		goto L_440;
	fprintf(stdout, "        PRODUCT OF SEARCH DIRECTION WITH BFGS-MATRIX:    DBD =%13.4e\n",
		*dbd);
	fprintf(stdout, "        PENALTY PARAMETER:  R =\n         ");
	for (j = 1; j <= n11nf.mnn; j++) {
		fprintf(stdout, "%16.8e", rpen[j - 1]);
	}
	fprintf(stdout, "\n");
L_440:
	;
	/*
	 * EVALUATION OF MERIT FUNCTION
	 */
L_450:
        _l0 = imerit + 3;
        _l1 = 4;
	l_n5ong(&_l0, &_n10nf->m, &_n10nf->me, n, &n11nf.mnn, &n11nf.nmnn,
                &_n10nf->acc, rpen, &_n10nf->f, df, g, dg, lddg, vmu, u, x,
                xl, xu, phi, dphi, active, wa, &_l1);
	if (!l7[4])
                _l0 = imerit + 4;
                _l1 = 4;
		l_n5ong(&_l0, &_n10nf->m, &_n10nf->me, n, &n11nf.mnn,
		        &n11nf.nmnn, &_n10nf->acc, rpen, &_n10nf->f, df, g,
                        dg, lddg, vmu, u, x, xl, xu, phi, dphi, active, wa, 
                        &_l1);
	*prd = imsl_sdot(*n, dphi, 1, del, 1);
	for (j = 1; j <= n11nf.mnn; j++) {
		*prd += dphi[j + *n - 1] * (u[j - 1] - vmu[j - 1]);
	}
	if (*prd < ze)
		goto L_480;
	sscal(n11nf.mnn, _n10nf->rpeno, rpen, 1);
        _l0 = 1;
	irpmax = l_ismax(n11nf.mnn, rpen, _l0);
	rpmax = imsl_f_max(rpen[irpmax - 1], F_ZERO);
	if (rpmax < _n10nf->rpenu)
		goto L_450;
	if (l7[4])
		goto L_470;
	if (!l7[3] || *dbd < _n10nf->acc)
		goto L_190;
	dcl11 = *DCL(n11nf.n1 - 1, n11nf.n1 - 1);
	if (dcl11 >= _n10nf->zefacu)
		goto L_190;
	dcl11 *= _n10nf->rpeno;
	if (_n10nf->lql)
		dcl11 *= _n10nf->rpeno;
	*DCL(n11nf.n1 - 1, n11nf.n1 - 1) = dcl11;
	goto L_140;
L_470:
	;
	_n10nf->ifail = 2;
	if (_n10nf->iprint == 0)
		goto L_350;
	fprintf(stdout, "        **SEARCH DIRECTION NOT PROFITABLE:  DPHI*P =%13.4e\n",
		*prd);
	goto L_350;
L_480:
	;
	if (_n10nf->iprint < 3)
		goto L_490;
	fprintf(stdout, "        PRODUCT LAGRANGIAN GRADIENT WITH SEARCH DIRECTION:  DLP =%13.4e\n",
		*prd);
L_490:
	;
	/* LINE SEARCH */
	wa[5] = xnm;
	wa[6] = delnm;
	l7[6] = IMSLFALSE;
	*iflise = 0;
L_500:
	ipr = 0;
	if (_n10nf->iprint >= 1000)
		ipr = _n10nf->iprint - 1000;
        _l0 = 35;
        _l1 = 10;
        _l2 = 5;
	l_n8ong(alphao, alpham, phi, prd, &_n10nf->amue, &_n10nf->imsl_beta,
	        iline, &_n10nf->maxfun, iflise, &ipr, &wa[5], &_l0, iw, &_l1,
                &active[ilwls - 1], &_l2);
	if (*iflise > -2)
		goto L_520;
	l7[6] = IMSLTRUE;
	*fbest = _n10nf->f;
	scopy(_n10nf->m, g, 1, gbest, 1);
	if (_n10nf->llise)
		goto L_500;
	scopy(*n, df, 1, dfbest, 1);
	for (i = 1; i <= *n; i++) {
		scopy(_n10nf->m, DG(i - 1, 0), 1, DGBEST(i - 1, 0), 1);
	}
	goto L_500;
L_520:
	;
	for (i = 1; i <= *n; i++) {
		x[i - 1] = xold[i - 1] + *alphao * del[i - 1];
	}
	for (j = 1; j <= n11nf.mnn; j++) {
		vmu[j - 1] = vmuold[j - 1] + *alphao * (u[j - 1] - vmuold[j - 1]);
	}
	if (*iflise == 0)
		goto L_570;
	if (*iflise == 1)
		goto L_560;
	if (*iflise > 1)
		goto L_550;
	goto L_600;
L_550:
	_n10nf->ifail = 1000 + *iflise;
	if (_n10nf->iprint == 0)
		goto L_350;
	fprintf(stdout, "        **ERROR IN LINE SEARCH. IFLISE =%4ld\n",
		*iflise);
	goto L_350;
L_560:
	_n10nf->ifail = 4;
	if (_n10nf->iprint == 0)
		goto L_350;
	fprintf(stdout, "        **MORE THAN MAXFUN FUNC-CALLS IN LINE SEARCH\n");
	goto L_350;
L_570:
	l7[2] = IMSLTRUE;
	if (_n10nf->iprint < 3)
		goto L_580;
	if (*iline == 1) {
		fprintf(stdout, "        LINE SEARCH SUCCESSFUL AFTER ONE STEP:  ALPHA = 1.\n");
	}
	if (*iline > 1) {
		fprintf(stdout, "        LINE SEARCH SUCCESSFUL AFTER%3ld STEPS:  ALPHA =%13.4e\n",
			*iline, *alphao);
	}
L_580:
	;
	if (!l7[6] && _n10nf->llise)
		goto L_630;
	if (!l7[6] && !_n10nf->llise)
		goto L_250;
	_n10nf->f = *fbest;
	scopy(_n10nf->m, gbest, 1, g, 1);
	if (_n10nf->llise)
		goto L_630;
	scopy(*n, dfbest, 1, df, 1);
	for (i = 1; i <= *n; i++) {
		scopy(_n10nf->m, DGBEST(i - 1, 0), 1, DG(i - 1, 0), 1);
	}
        _l0 = imerit + 1;
        _l1 = 4;
	l_n5ong(&_l0, &_n10nf->m, &_n10nf->me, n, &n11nf.mnn, &n11nf.nmnn,
                &_n10nf->acc, rpen, &_n10nf->f, df, g, dg, lddg, vmu, u, x,
                xl, xu, phi, dphi, active, wa, &_l1);
	goto L_250;
L_600:
	;
	/* NEW FUNCTION AND GRADIENT VALUES */
	if (l7[5]) {
		_n10nf->ifail = -1;
		goto L_9000;
	}
	imsl_e1usr("ON");
	(*fcns) (_n10nf->m, _n10nf->me, *n, x, &active[*mmax], &_n10nf->f,
		 g);
	imsl_e1usr("OFF");
L_610:
	;
	if (l7[0])
		_n10nf->f *= *scf;

	if (_n10nf->m != 0 && l7[1]) {
		for (j = 1; j <= _n10nf->m; j++) {
			g[j - 1] *= scg[j - 1];
		}
	}
	*nfunc += 1;
        _l0 = imerit + 3;
        _l1 = 4;
	l_n5ong(&_l0, &_n10nf->m, &_n10nf->me, n, &n11nf.mnn, &n11nf.nmnn,
                &_n10nf->acc, rpen, &_n10nf->f, df, g, dg, lddg, vmu, u, x,
                xl, xu, phi, dphi, active, wa, &_l1);
	if (_n10nf->llise && !l7[2])
		goto L_500;
L_630:
	;
        _l0 = imerit + 1;
        _l1 = 4;
	l_n5ong(&_l0, &_n10nf->m, &_n10nf->me, n, &n11nf.mnn, &n11nf.nmnn,
                &_n10nf->acc, rpen, &_n10nf->f, df, g, dg, lddg, vmu, u, x,
                xl, xu, phi, dphi, active, wa, &_l1);
	if (l7[0])
		_n10nf->f /= *scf;
	if (_n10nf->m != 0 && l7[1]) {
		for (j = 1; j <= _n10nf->m; j++) {
			g[j - 1] /= scg[j - 1];
		}
	}
	if (l7[5]) {
		_n10nf->ifail = -2;
		goto L_9000;
	}
	l_n5onf(fcns, &_n10nf->m, &_n10nf->me, mmax, n, x, xs, active,
		   &_n10nf->f, g, df, dg, cwk);
L_650:
	;
	*ngrad += 1;
	if (l7[0]) {
		_n10nf->f *= *scf;
		sscal(*n, *scf, df, 1);
	}
	if (_n10nf->m != 0 && l7[1]) {
		for (j = 1; j <= _n10nf->m; j++) {
			g[j - 1] *= scg[j - 1];
			if (active[j - 1])
				sscal(*n, scg[j - 1], DG(0, j - 1), *lddg);
		}
	}
	if (l7[2])
		goto L_250;
        _l0 = imerit + 4;
        _l1 = 4;
	l_n5ong(&_l0, &_n10nf->m, &_n10nf->me, n, &_n11nf->mnn,
	     &n11nf.nmnn, &_n10nf->acc, rpen, &_n10nf->f, df, g, dg, lddg,
		   vmu, u, x, xl, xu, phi, dphi, active, wa, &_l1);
	*prd = imsl_sdot(*n, dphi, 1, del, 1);
	for (j = 1; j <= n11nf.mnn; j++) {
		*prd += dphi[*n + j - 1] * (u[j - 1] - vmuold[j - 1]);
	}
	goto L_500;
	/* UPDATE HESSIAN OF LAGRANGIAN */
L_680:
	*dbd *= *alphao ** alphao;
	sscal(*n, *alphao, bdel, 1);
	for (i = 1; i <= *n; i++) {
		eta[i - 1] = dla[i - 1] - dlaold[i - 1];
	}
	edel = *alphao * imsl_sdot(*n, del, 1, eta, 1);
	dbd1 = _n10nf->dbdfac ** dbd;
	if (edel >= dbd1)
		goto L_720;
	theta = (*dbd - dbd1) / (*dbd - edel);
	theta1 = on - theta;
	for (i = 1; i <= *n; i++) {
		eta[i - 1] = theta * eta[i - 1] + theta1 * bdel[i - 1];
	}

	edel = dbd1;
L_720:
	;
	dbdi = sqrt(on / *dbd);
	edeli = sqrt(on / edel);
	/* UPDATE FACTORIZATION */
	sscal(*n, dbdi, bdel, 1);
	sscal(*n, edeli, eta, 1);
	if (_n10nf->lql) {
		for (i = 1; i <= *n; i++) {
			for (j = 1; j <= i; j++) {
				*DCL(i - 1, j - 1) += eta[i - 1] * eta[j - 1] - bdel[i - 1] *
					bdel[j - 1];
			}
		}
		imsl_csfrg(n, dcl, lddcl);
		goto L_60;
	}
	l_n7ong(n, dcl, lddcl, cd, eta, bdel);
	/* CORRECT DATA FOR QL-SOLUTION */
L_750:
	for (i = 1; i <= *n; i++) {
		sqd = sqrt(cd[i - 1]);
		if (sqd > uf)
			goto L_760;
		_n10nf->ifail = 3;
		if (_n10nf->iprint == 0)
			goto L_350;
		fprintf(stdout, "        **UNDERFLOW IN BFGS-UPDATE\n");
		goto L_350;
L_760:
		;
		if (i < *n)
                        _l0 = *n - i;
                        _l1 = 1;
			imsl_svcal(_l0, sqd, DCL(i - 1, i), _l1, 
                                   DCL(i, i - 1), *lddcl);
		*DCL(i - 1, i - 1) = sqd;
	}
	if (*iter == 0)
		goto L_50;
	/* PERFORM NEXT ITERATION */
	goto L_60;
L_9000:
	return;
}				/* end of function */
#undef DG
#undef DCL
#undef DGBEST
/*----------------------------------------------------------------------- */

/*  IMSL Name:  N4ONG/DN4ONG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    December 2, 1985

    Purpose:    Main driver for the successive quadratic programming
                algorithm.

    Usage:      CALL N4ONG (FCNS, GRAD, MMAX, N, NMAX, X, G, DF, DG,
                            LDDG, U, XL, XU, DCL, LDDCL, CD, VMU,
                            DEL, DLA, DCLF, BDEL, ETA, XOLD, DLAOLD, V,
                            W, VMUOLD, DPHI, RPEN, SCG, FBEST, DFBEST,
                            GBEST, DGBEST, WA, LWA, MNN2, MO1, NFUNC,
                            NGRAD, ITER, NQL, ILINE, IFLISE, NOPT, IW,
                            LIW, PHI, DFDEL, DBD, ALPHAM, ALPHAO, SCF,
                            PRD, ACTIVE, L7)

    Arguments:
       FCNS   - User-supplied SUBROUTINE to evaluate the functions at
                a given poMint.  The usage is
                CALL FCNS (M, ME, N, X, ACTIVE, F, G), where
                M      - Total number of constraMints.  (Input)
                ME     - Number of equality constraMints.  (Input)
                N      - Number of variables.  (Input)
                X      - The poMint at which the function is evaluated.
                         (Input)
                         X should not be changed by FCNS.
                ACTIVE - Logical vector of length MMAX indicating the
                         active constraMints.  (Input)
                F      - The computed function value at the poMint X.
                         (Output)
                G      - Vector of length MMAX containing the values of
                         constraMints at poMint X.  (Output)
                FCNS must be declared EXTERNAL in the calling program.
       GRAD   - User-supplied SUBROUTINE to evaluate the gradients at
                a given poMint.  The usage is
                CALL GRAD (M, ME, MMAX, N, X, ACTIVE, F, G, DF, DG),
                where
                M      - Total number of constraMints.  (Input)
                ME     - Number of equality constraMints.  (Input)
                MMAX   - Maximum of (1, M).  (Input)
                N      - Number of variables.  (Input)
                X      - The poMint at which the function is evaluated.
                         (Input)
                         X should not be changed by FCNS.
                ACTIVE - Logical vector of length MMAX indicating the
                         active constraMints.  (Input)
                F      - The computed function value at the poMint X.
                         (Output)
                G      - Vector of length MMAX containing the values of
                         constraMints at poMint X.  (Output)
                DF     - Vector of lenght N containing the value of the
                         gradient of the objective function.  (Output)
                DG     - MMAX by N array containing the values of the
                         gradients for the active constraMints.  (Output)
                GRAD must be declared EXTERNAL in the calling program.
       MMAX   - Order of the array DG.  (Input)
                MMAX must be at least MAX(1,M).
       N      - Number of variables.  (Input)
       NMAX   - Order of DCL where NMAX must be at least MAX(2,N+1).
                (Input)
       X      - Vector of length N containing the initial guesses to the
                solution on input and the solution on output.
                (Input/Output)
       G      - Vector of length MMAX containing constraMint values.
                (Output)
       DF     - Vector of length N+1 containing the gradient of the
                of the objective function.  (Output)
       DG     - Array of dimension MMAX by MMAX containing the gradient
                of the constraMints.  (Output)
       LDDG   - Leading dimension of DG exactly as specified in the
                dimension statement in the calling program.  (Input)
       U      - Vector of length MNN2 containing the multipliers of the
                nonlinear constraMints and the bounds.  (Output)
       XL     - Vector of length N containing the lower bounds for the
                variables.  (Input)
       XU     - Vector of length N containing the upper bounds for the
                variables.  (Input)
       DCL    - Array of dimension NMAX by NMAX containing an the final
                approximation to the Hessian.  (Output)
       LDDCL  - Leading dimension of DCL exactly as specified in the
                dimension statement in the calling program.  (Input)
       CD     - Vector of length NMAX containing the diagonal elements of
                the Hessian.  (Output)
       VMU    - Work vector of length M + 2*N.
       DEL    - Work vector of length N + 1.
       DLA    - Work vector of length N.
       DCLF   - Work vector of length N + 1.
       BDEL   - Work vector of length N.
       ETA    - Work vector of length N.
       XOLD   - Work vector of length N.
       DLAOLD - Work vector of length N.
       V      - Work vector of length N + 1.
       W      - Work vector of length N + 1.
       VMUOLD - Work vector of length M + 2*N.
       DPHI   - Work vector of length M + 3*N.
       RPEN   - Work vector of length M + 2*N.
       SCG    - Work vector of length MMAX.
       FBEST  - Work scalar.
       DFBEST - Work vector of length NMAX.
       GBEST  - Work vector of length MMAX.
       DGBEST - Work array of dimension MMAX by N.
       WA     - Work vector of length LWA.
       LWA    - Length of WA where LWA = N*(2*N+13) + M + MMAX + 12.
                (Input)
       MNN2   - Scalar containing the value M + 2*N + 2.  (Input)
       MO1    - Scalar containing the value 1 when LLISE is true or MMAX
                when LLISE is false.  (Input)
       NFUNC  - Number of function evaluations.  (Output)
       NGRAD  - Number of gradient evaluations.  (Output)
       ITER   - Number of iterations.  (Output)
       NQL    - Number of QL algorithm evaluations.  (Output)
       ILINE  - Number of line search evaluations.  (Output)
       IFLISE - Error parameter for line search algorithm.  (Output)
       NOPT   - Number of optimality iterations.  (Output)
       IW     - Work vector of length LIW.
       LIW    - Length of IW where LIW = 12.  (Input)
       PHI    - Scalar variable.
       DFDEL  - Scalar variable.
       DBD    - Scalar variable.
       ALPHAM - Scalar variable.
       ALPHAO - Scalar variable.
       SCF    - Scalar variable.
       PRD    - Scalar variable.
       ACTIVE - Logical vector of length LACTIV indicating which
                constraMints are active.  (Output)
       L7     - Logical vector of length 7.

    Remark:
       The NLPQL algorithm was designed by K. Schittkowski.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_n4ong(void (*fcns)(Mint, Mint, Mint, Mfloat*, Mint*, Mfloat*,
                    Mfloat*), void (*grad)(Mint, Mint, Mint, Mint, Mfloat*,
                    Mint*, Mfloat, Mfloat*, Mfloat*, Mfloat*), Mint *mmax,
                    Mint *n, Mint *nmax, Mfloat x[], Mfloat g[], Mfloat df[],
                    Mfloat dg[], Mint *lddg, Mfloat u[], Mfloat xl[],
                    Mfloat xu[], Mfloat *dcl, Mint *lddcl, Mfloat cd[],
                    Mfloat vmu[], Mfloat del[], Mfloat dla[], Mfloat dclf[],
                    Mfloat bdel[], Mfloat eta[], Mfloat xold[], Mfloat dlaold[],
                    Mfloat v[], Mfloat w[], Mfloat vmuold[], Mfloat dphi[],
                    Mfloat rpen[], Mfloat scg[], Mfloat *fbest, Mfloat dfbest[],
                    Mfloat gbest[], Mfloat *dgbest, Mfloat wa[], Mint *lwa, 
                    Mint *mnn2, Mint *mo1, Mint *nfunc, Mint *ngrad, Mint *iter,
                    Mint *nql, Mint *iline, Mint *iflise, Mint *nopt, Mint iw[],
                    Mint *liw, Mfloat *phi, Mfloat *dfdel, Mfloat *dbd, 
                    Mfloat *alpham, Mfloat *alphao, Mfloat *scf, Mfloat *prd,
                    Mint active[], Mint l7[])
#else
static void l_n4ong(void (*fcns)(Mint, Mint, Mint, Mfloat[], Mint[], Mfloat*,
                    Mfloat[]), void (*grad)(Mint, Mint, Mint, Mint, Mfloat[],
                    Mint[], Mfloat, Mfloat[], Mfloat[], Mfloat[]), Mint *mmax,
                    Mint *n, Mint *nmax, Mfloat x[], Mfloat g[], Mfloat df[],
                    Mfloat dg[], Mint *lddg, Mfloat u[], Mfloat xl[],
                    Mfloat xu[], Mfloat *dcl, Mint *lddcl, Mfloat cd[],
                    Mfloat vmu[], Mfloat del[], Mfloat dla[], Mfloat dclf[],
                    Mfloat bdel[], Mfloat eta[], Mfloat xold[], Mfloat dlaold[],
                    Mfloat v[], Mfloat w[], Mfloat vmuold[], Mfloat dphi[],
                    Mfloat rpen[], Mfloat scg[], Mfloat *fbest, Mfloat dfbest[],
                    Mfloat gbest[], Mfloat *dgbest, Mfloat wa[], Mint *lwa, 
                    Mint *mnn2, Mint *mo1, Mint *nfunc, Mint *ngrad, Mint *iter,
                    Mint *nql, Mint *iline, Mint *iflise, Mint *nopt, Mint iw[],
                    Mint *liw, Mfloat *phi, Mfloat *dfdel, Mfloat *dbd, 
                    Mfloat *alpham, Mfloat *alphao, Mfloat *scf, Mfloat *prd,
                    Mint active[], Mint l7[])
#endif
#else
static void l_n4ong(fcns, grad, mmax, n, nmax, x, g, df, dg, lddg, u, xl, xu,
                    dcl, lddcl, cd, vmu, del, dla, dclf, bdel, eta, xold,
                    dlaold, v, w, vmuold, dphi, rpen, scg, fbest, dfbest, gbest,
                    dgbest, wa, lwa, mnn2, mo1, nfunc, ngrad, iter, nql, iline,
                    iflise, nopt, iw, liw, phi, dfdel, dbd, alpham, alphao, scf,
                    prd, active, l7)
	void            (*fcns) (), (*grad) ();
	Mint            *mmax, *n, *nmax;
	Mfloat           x[], g[], df[], dg[];
	Mint            *lddg;
	Mfloat           u[], xl[], xu[], *dcl;
	Mint            *lddcl;
	Mfloat           cd[], vmu[], del[], dla[], dclf[], bdel[], eta[],
	                xold[], dlaold[], v[], w[], vmuold[], dphi[], rpen[],
	                scg[], *fbest, dfbest[], gbest[], *dgbest, wa[];
	Mint            *lwa, *mnn2, *mo1, *nfunc, *ngrad, *iter, *nql,
	               *iline, *iflise, *nopt, iw[], *liw;
	Mfloat          *phi, *dfdel, *dbd, *alpham, *alphao, *scf, *prd;
	Mint            active[], l7[];
#endif
{
#define DG(I_,J_)	(dg+(I_)*(*lddg)+(J_))
#define DCL(I_,J_)	(dcl+(I_)*(*lddcl)+(J_))
#define DGBEST(I_,J_)	(dgbest+(I_)*(*mo1)+(J_))
        FILE            *iout;
	Mint            _l0, _l1, _l2, i, ifail1, ilwls, imerit,
	                ipr, irpmax, j, k, liwql, lwaql, 
	                mmax2, mn, mn1, n2, nact;
	Mfloat          dbd1, dbdi, dcl11, delnm, dlan, edel,
	                edeli, eps0, fact, ff, of, on, rpmax, sdcl11,
	                sqacc, sqd, sres, sum, theta, theta1, uad, uf,
	                xnm, ze;
        Mfloat          *dgt = NULL;
        Mint            IMSLFALSE = 0, IMSLTRUE = 1;
	static struct {
		Mfloat           f, acc, scbou, dbdfac, zefac, rpeno, rpens,
		                rpenu, zefacu, delta, imsl_beta, amue,
		                alm;
		Mint             m, me, maxfun, maxit, iprint, mode, ifail;
		Mint            llise, lql, lmerit;
	}              *_n10ng = (void *) &n10ng;
	static struct {
		Mint             n1, lact, no1, mnn, nmnn;
	}              *_n11ng = (void *) &n11ng;

        dgt = (Mfloat *) imsl_malloc (*mmax * *n * sizeof(*dgt));

	/* CONSTANT DATA */
	ze = F_ZERO;
	on = F_ONE;
	eps0 = imsl_amach(4);
	uf = eps0 * eps0;
	of = on / uf;
	imsl_umach(2, &iout);
	/* INITIAL DEFINITIONS */
	mn = _n10ng->m + *n;
	n2 = *n + *n;
	lwaql = *lwa - *mmax - 40;
	liwql = *liw - 10;
	ilwls = 2 ** mmax + 1;
	imerit = 0;
	if (!_n10ng->lmerit)
		imerit = 4;
	l7[5] = IMSLFALSE;
	l7[3] = IMSLFALSE;
	l7[4] = IMSLFALSE;
	sqacc = sqrt(_n10ng->acc);
	if (((_n10ng->mode == 2 || _n10ng->mode == 7) || _n10ng->mode ==
	     3) || _n10ng->mode == 8) {
		l7[5] = IMSLTRUE;
		if (_n10ng->ifail == -1)
			goto L_630;
		if (_n10ng->ifail == -2)
			goto L_660;
	}
	*iline = 0;
	*alphao = ze;
	*nfunc = 0;
	*ngrad = 0;
	*iter = 0;
	*nql = 0;
	*nopt = 0;
	if (_n10ng->m != 0) {
		mmax2 = *mmax + *mmax;
		for (j = 1; j <= mmax2; j++) {
			active[j - 1] = IMSLTRUE;
		}
	}
	if (!l7[5]) {
		imsl_e1usr("ON");
		(*fcns) (_n10ng->m, _n10ng->me, *n, x, &active[*mmax],
                         &_n10ng->f, g);
		(*grad) (_n10ng->m, _n10ng->me, *mmax, *n, x, active,
                         _n10ng->f, g, df, dgt);
		imsl_e1usr("OFF");
	}

        for (i=0; i<*n; i++) {
             k = i * *mmax;
             scopy (*mmax, &dgt[i], *n, &dg[k], 1);
        }

	l7[0] = IMSLFALSE;
	l7[1] = IMSLFALSE;
	if (fabs(_n10ng->f) >= _n10ng->scbou) {
		l7[0] = IMSLTRUE;
		if (_n10ng->scbou > ze)
			*scf = F_ONE / sqrt(fabs(_n10ng->f));
		_n10ng->f *= *scf;
		sscal(*n, *scf, df, 1);
	}
	if (_n10ng->m != 0) {
		for (j = 1; j <= _n10ng->m; j++) {
			if (fabs(g[j - 1]) >= _n10ng->scbou)
				l7[1] = IMSLTRUE;
		}
	}
	if (l7[1]) {
		for (j = 1; j <= _n10ng->m; j++) {
			if (_n10ng->scbou > ze)
				scg[j - 1] = F_ONE / imsl_f_max(F_ONE, sqrt(fabs(g[j - 1])));
			g[j - 1] *= scg[j - 1];
			sscal(*n, scg[j - 1], DG(0, j - 1), *lddg);
		}
	}
	if (_n10ng->iprint >= 1) {
		if (l7[0] && !l7[1]) {
			fprintf(stdout, "\n     OBJECTIVE FUNCTION WILL BE SCALED\n");
		}
		if (l7[0] && l7[1]) {
			fprintf(stdout, "\n     OBJECTIVE AND CONSTRAINT FUNCTIONS WILL BE SCALED\n");
		}
		if (!l7[0] && l7[1]) {
			fprintf(stdout, "\n     CONSTRAINT FUNCTIONS WILL BE SCALED\n");
		}
	}
	*nfunc += 1;
	*ngrad += 1;
	dclf[_n11ng->n1 - 1] = F_ZERO;
	scopy(*n, df, 1, del, 1);
	sscal(*n, -F_ONE, del, 1);
	sset(*n, F_ZERO, DCL(0, _n11ng->n1 - 1), *lddcl);
	sset(*n, F_ZERO, DCL(_n11ng->n1 - 1, 0), 1);
	*DCL(_n11ng->n1 - 1, _n11ng->n1 - 1) = _n10ng->zefac;
	if (((_n10ng->mode == 1 || _n10ng->mode == 6) || _n10ng->mode ==
	     3) || _n10ng->mode == 8) {
		if (_n10ng->lql)
			goto L_50;
		goto L_760;
	}
	sset(*n, F_ONE, cd, 1);
	for (i = 1; i <= *n; i++) {
		sset(*n, F_ZERO, DCL(i - 1, 0), 1);
	}
	sset(*n, F_ONE, DCL(0, 0), *lddcl + 1);
L_50:
	sset(_n11ng->mnn, _n10ng->rpens, rpen, 1);
	sset(_n11ng->mnn, F_ZERO, vmu, 1);
	if (((_n10ng->mode == 1 || _n10ng->mode == 6) || _n10ng->mode ==
	     3) || _n10ng->mode == 8)
		scopy(_n11ng->mnn, u, 1, vmu, 1);
        _l0 = imerit + 3;
        _l1 = 4;
	l_n5ong(&_l0, &_n10ng->m, &_n10ng->me, n, &_n11ng->mnn, &_n11ng->nmnn,
                &_n10ng->acc, rpen, &_n10ng->f, df, g, dg, lddg, vmu, u, x,
                xl, xu, phi, dphi, active, wa, &_l1);
	/*
	 * START MAIN LOOP, PRINT INTERMEDIATE ITERATES
	 */
L_60:
	;
	l7[2] = IMSLFALSE;
	if (_n10ng->iprint < 3)
		goto L_90;
	if (l7[0])
		_n10ng->f /= *scf;
	fprintf(stdout, "\n\n     ITERATION%3ld\n\n        FUNCTION VALUE:  F(X) =%16.8e",
		*iter, _n10ng->f);
	fprintf(stdout, "\n        VARIABLE:  X =\n         ");
	for (i = 1; i <= *n; i++) {
		fprintf(stdout, "%16.8e", x[i - 1]);
	}
	fprintf(stdout, "\n");
	if (l7[0])
		_n10ng->f *= *scf;
	if (_n10ng->m != 0 && (l7[0] || l7[1])) {
		if (l7[0])
			sscal(_n10ng->m, F_ONE / *scf, vmu, 1);
		if (l7[1]) {
			for (j = 1; j <= _n10ng->m; j++) {
				vmu[j - 1] *= scg[j - 1];
				g[j - 1] /= scg[j - 1];
			}
		}
	}
	fprintf(stdout, "        MULTIPLIERS:  U =\n         ");
	for (j = 1; j <= _n11ng->mnn; j++) {
		fprintf(stdout, "%16.8e", vmu[j - 1]);
	}
	fprintf(stdout, "\n");
	if (_n10ng->m != 0) {
		fprintf(stdout, "        CONSTRAINTS: G(X) =\n         ");
		for (j = 1; j <= _n10ng->m; j++) {
			fprintf(stdout, "%16.8e", g[j - 1]);
		}
		fprintf(stdout, "\n");
		if (l7[0] || l7[1]) {
			if (l7[0])
				sscal(_n10ng->m, *scf, vmu, 1);
			if (l7[1]) {
				for (j = 1; j <= _n10ng->m; j++) {
					vmu[j - 1] /= scg[j - 1];
					g[j - 1] *= scg[j - 1];
				}
			}
		}
	}
L_90:
	*iter += 1;
	if (*iter < _n10ng->maxit)
		goto L_100;
	_n10ng->ifail = 1;
	if (_n10ng->iprint == 0)
		goto L_370;
	fprintf(stdout, "        **MORE THAN MAXIT ITERATIONS\n");
	goto L_370;
L_100:
	;
	/* SEARCH DIRECTION */
	scopy(*n, df, 1, dclf, 1);
	for (i = 1; i <= *n; i++) {
		v[i - 1] = xl[i - 1] - x[i - 1];
		w[i - 1] = xu[i - 1] - x[i - 1];
	}
	ipr = 0;
	if (_n10ng->iprint > 10 && _n10ng->iprint < 1000)
		ipr = _n10ng->iprint - 10;
	if (_n10ng->mode >= 5)
		goto L_120;
	ifail1 = *iter;
	if (l7[3] || l7[4])
		ifail1 = 1;
	iw[10] = 0;
	if (_n10ng->lql)
		iw[10] = 1;
	iw[11] = 0;
	l_n6ong(&_n10ng->m, &_n10ng->me, mmax, n, nmax, &_n11ng->mnn, dcl,
		   lddcl, dclf, dg, lddg, g, v, w, del, u, &ifail1, &ipr, &wa[*mmax + 40],
		   &lwaql, &iw[10], &liwql);
	del[_n11ng->n1 - 1] = ze;
	*nql += 1;
	l7[3] = IMSLFALSE;
	if (ifail1 == 0)
		goto L_210;
L_120:
	;
	if (*iter == 1)
		goto L_130;
	fact = F_TWO * fabs(*dbd ** dfdel) / (sqrt(*dbd) * (on - del[_n11ng->n1 - 1]));
	if (_n10ng->lql)
		fact *= fact;
	dcl11 = imsl_f_max(_n10ng->zefac, fact);
	*DCL(_n11ng->n1 - 1, _n11ng->n1 - 1) = imsl_f_min(_n10ng->zefacu, dcl11);
L_130:
	;
	sset(*n, F_ZERO, del, 1);
	del[_n11ng->n1 - 1] = F_ONE;

	if (_n10ng->m != 0) {
		scopy(_n10ng->m, g, 1, DG(_n11ng->n1 - 1, 0), 1);
		sscal(_n10ng->m, -F_ONE, DG(_n11ng->n1 - 1, 0), 1);
		for (j = 1; j <= _n10ng->m; j++) {
			if (!active[j - 1])
				*DG(_n11ng->n1 - 1, j - 1) = F_ZERO;
		}
	}
	v[_n11ng->n1 - 1] = F_ZERO;
	w[_n11ng->n1 - 1] = F_ONE;
	ifail1 = -*iter;
	if (!l7[3] || l7[4])
		ifail1 = -1;
	iw[10] = 0;
	if (_n10ng->lql)
		iw[10] = 1;
	iw[11] = 1;
	l_n6ong(&_n10ng->m, &_n10ng->me, mmax, &_n11ng->n1, nmax, mnn2,
		   dcl, lddcl, dclf, dg, lddg, g, v, w, del, u, &ifail1, &ipr, &wa[*mmax + 40],
		   &lwaql, &iw[10], &liwql);
	*nql += 1;
	mn1 = _n10ng->m + _n11ng->n1 + 1;
	l7[3] = IMSLTRUE;
	if (ifail1 == 0)
		goto L_160;
L_150:
	_n10ng->ifail = 10 + ifail1;
	if (_n10ng->iprint == 0)
		goto L_370;
	fprintf(stdout, "        **ERROR IN QL. IFAIL(QL) =%3ld\n",
		ifail1);
	goto L_370;
L_160:
	;
	scopy(*n + 1, &u[mn1 - 1], 1, &u[mn1 - 2], 1);
	if (_n10ng->iprint < 3)
		goto L_170;
	fprintf(stdout, "        ADDITIONAL VARIABLE TO PREVENT INCONSISTENCY:  DELTA =%13.4e\n",
		del[_n11ng->n1 - 1]);
	sdcl11 = *DCL(_n11ng->n1 - 1, _n11ng->n1 - 1);
	if (!_n10ng->lql)
		sdcl11 = sqrt(sdcl11);
	fprintf(stdout, "        PENALTY PARAMETER FOR DELTA:  RHO =%13.4e\n",
		sdcl11);
L_170:
	;
	dcl11 = *DCL(_n11ng->n1 - 1, _n11ng->n1 - 1);
	if (del[_n11ng->n1 - 1] < _n10ng->delta)
		goto L_210;
	*DCL(_n11ng->n1 - 1, _n11ng->n1 - 1) = dcl11 * _n10ng->rpeno;
	if (_n10ng->lql)
		*DCL(_n11ng->n1 - 1, _n11ng->n1 - 1) *= _n10ng->rpeno;
	if (dcl11 < _n10ng->zefacu)
		goto L_130;
	/*
	 * AUGMENTED LAGRANGIAN TYPE SEARCH DIRECTION
	 */
L_180:
	l7[4] = IMSLTRUE;
	if (_n10ng->iprint < 3)
		goto L_190;
	fprintf(stdout, "        **WARNING: AUGMENTED LAGRANGIAN SEARCH DIRECTION\n");
L_190:
        _l0 = 4;
	l_n5ong(&_l0, &_n10ng->m, &_n10ng->me, n, &_n11ng->mnn, &_n11ng->nmnn,
	        &_n10ng->acc, rpen, &_n10ng->f, df, g, dg, lddg, vmu, u, x, xl,
		xu, phi, dphi, active, wa, &_l0);
	scopy(*n, dphi, 1, &wa[40], 1);
	scopy(*n, dphi, 1, dclf, 1);
	ifail1 = 1;
	iw[10] = 0;
	if (_n10ng->lql)
		iw[10] = 1;
	iw[11] = 0;
        _l0 = 0;
	l_n6ong(&_l0, &_l0, mmax, n, nmax, mnn2, dcl, lddcl, dclf, dg, lddg,
                g, v, w, del, u, &ifail1, &ipr, &wa[*mmax + 40], &lwaql,
                &iw[10], &liwql);
	if (ifail1 > 0)
		goto L_150;
	if (_n10ng->m == 0)
		goto L_220;
	scopy(n2, u, -1, &u[_n10ng->m], -1);

	for (j = 1; j <= _n10ng->m; j++) {
		u[j - 1] = vmu[j - 1] - dphi[*n + j - 1];
	}
	goto L_220;
L_210:
	l7[4] = IMSLFALSE;
	/*
	 * PROJECTION OF DEL, MAXIMAL STEPLENGTH, AND NORM OF X,DEL
	 */
L_220:
	*alpham = of;
	xnm = F_ZERO;
	delnm = F_ZERO;
	for (i = 1; i <= *n; i++) {
		if (w[i - 1] < del[i - 1])
			del[i - 1] = w[i - 1];
		if (v[i - 1] > del[i - 1])
			del[i - 1] = v[i - 1];
		uad = fabs(del[i - 1]);
		if (del[i - 1] > uf)
			*alpham = imsl_f_min(*alpham, w[i - 1] / del[i - 1]);
		if (del[i - 1] < -uf)
			*alpham = imsl_f_min(*alpham, v[i - 1] / del[i - 1]);
		xnm = imsl_f_max(fabs(x[i - 1]), xnm);
		delnm = imsl_f_max(uad, delnm);
	}
	*alpham = imsl_f_max(on, *alpham);
	*alpham = imsl_f_min(*alpham, _n10ng->alm);
	/* GRADIENT OF LAGRANGIAN */
L_240:
	for (i = 1; i <= *n; i++) {
		uad = df[i - 1];
		if (l7[4])
			uad = dphi[i - 1];
		dla[i - 1] = uad - u[_n10ng->m + i - 1] - u[mn + i - 1];
	}

	if (_n10ng->m != 0 && !l7[4]) {
		for (j = 1; j <= _n10ng->m; j++) {
			if (u[j - 1] != ze)
				saxpy(*n, -u[j - 1], DG(0, j - 1), *lddg, dla, 1);
		}
	}
	if (l7[2])
		goto L_710;
	/* STORE SOME DATA */
	scopy(*n, dla, 1, dlaold, 1);
	scopy(*n, x, 1, xold, 1);
	*dfdel = imsl_sdot(*n, df, 1, del, 1);
	irpmax = imsl_isamax(*n, dla, 1);
	dlan = fabs(dla[irpmax - 1]);
	scopy(_n11ng->mnn, vmu, 1, vmuold, 1);
	/* DETERMINE B*D AND D(T)*B*D */
	if (!_n10ng->lql) {
		for (i = 1; i <= *n; i++) {
			eta[i - 1] = imsl_sdot(*n - i + 1, DCL(i - 1, i - 1), *lddcl,
					       &del[i - 1], 1);
		}
		*dbd = imsl_sdot(*n, eta, 1, eta, 1);
		for (i = 1; i <= *n; i++) {
			bdel[i - 1] = imsl_sdot(i, DCL(i - 1, 0), 1, eta, 1);
		}
	} else {
		for (i = 1; i <= *n; i++) {
			bdel[i - 1] = imsl_sdot(*n, DCL(0, i - 1), *lddcl, del, 1);
		}
		*dbd = imsl_sdot(*n, bdel, 1, del, 1);
	}
	/* TEST FOR OPTIMALITY AND FINAL OUTPUT */
	sres = ze;
	sum = fabs(*dfdel);
	if (l7[0])
		sum /= *scf;
	nact = 0;
	if (_n10ng->m == 0)
		goto L_310;
	for (j = 1; j <= _n10ng->m; j++) {
		if (active[j - 1])
			nact += 1;
		uad = fabs(g[j - 1]);
		if (l7[1])
			uad /= scg[j - 1];
		if (j <= _n10ng->me || g[j - 1] < ze)
			sres += uad;
	}
        _l0 = 1;
	sum += l_a1ot(_n10ng->m, u, _l0, g, _l0);
	if (_n10ng->iprint != 3)
		goto L_310;
	fprintf(stdout, "        SUM OF CONSTRAINT VIOLATIONS:                    SCV =%13.4e\n",
		sres);
	fprintf(stdout, "        NUMBER OF ACTIVE CONSTRAINTS:                    NAC =%4ld\n",
		nact);
L_310:
	;
	for (i = 1; i <= *n; i++) {
		sum += fabs(u[_n10ng->m + i - 1] * v[i - 1]) + fabs(u[mn + i - 1] *
								  w[i - 1]);
	}
	if (_n10ng->iprint != 2)
		goto L_330;
	ff = _n10ng->f;
	if (l7[0])
		ff = _n10ng->f / *scf;
	fprintf(stdout, " %3ld%16.8e%10.2e%4ld%3ld%10.2e%10.2e%10.2e%10.2e\n",
		*iter, ff, sres, nact, *iline, *alphao, del[_n11ng->n1 - 1],
		dlan, sum);
L_330:
	;
	if (_n10ng->iprint != 3)
		goto L_340;
	fprintf(stdout, "        KUHN-TUCKER OPTIMALITY CONDITION:                KTO =%13e\n",
		sum);
	fprintf(stdout, "        NORM OF LAGRANGIAN GRADIENT:                     NLG =%13e\n",
		dlan);
L_340:
	;
	if (*dbd >= uf)
		goto L_350;
	if (sres < sqacc)
		goto L_360;
	if (*dbd > ze)
		goto L_410;
	if (!l7[4])
		goto L_180;
	_n10ng->ifail = 7;
	if (_n10ng->iprint == 0)
		goto L_370;
	fprintf(stdout, "        **UNDERFLOW IN D(T)*B*D AND INFEASIBLE ITERATE X\n");
	goto L_370;
L_350:
	;
	if (sum >= _n10ng->acc || sres > sqacc)
		goto L_410;
	if (dlan <= sqrt(sqacc) || *dbd <= _n10ng->acc)
		goto L_360;
	*nopt += 1;
	if (*nopt < 3)
		goto L_410;
L_360:
	_n10ng->ifail = 0;
L_370:
	;
	if (l7[0])
		_n10ng->f /= *scf;
	if (_n10ng->m == 0 || (!l7[0] && !l7[1]))
		goto L_390;
	if (l7[0])
		sscal(_n10ng->m, F_ONE / *scf, u, 1);
	if (l7[1]) {
		for (j = 1; j <= _n10ng->m; j++) {
			u[j - 1] *= scg[j - 1];
			g[j - 1] /= scg[j - 1];
		}
	}
L_390:
	;

	if (_n10ng->iprint == 0)
		goto L_9000;
	fprintf(stdout, "\n\n     * FINAL CONVERGENCE ANALYSIS\n\n");
	fprintf(stdout, "        OBJECTIVE FUNCTION VALUE:  F(X) =%16.8e\n",
		_n10ng->f);
	fprintf(stdout, "        APPROXIMATION OF SOLUTION:  X =\n         ");
	for (i = 1; i <= *n; i++) {
		fprintf(stdout, "%16.8e", x[i - 1]);
	}
	fprintf(stdout, "\n");
	fprintf(stdout, "        APPROXIMATION OF MULTIPLIERS:  U =\n         ");
	for (j = 1; j <= _n11ng->mnn; j++) {
		fprintf(stdout, "%16.8e", u[j - 1]);
	}
	fprintf(stdout, "\n");
	if (_n10ng->m == 0)
		goto L_400;
	fprintf(stdout, "        CONSTRAINT VALUES:  G(X) =\n         ");
	for (j = 1; j <= _n10ng->m; j++) {
		fprintf(stdout, "%16.8e", g[j - 1]);
	}
	fprintf(stdout, "\n");
L_400:
	fprintf(stdout, "        DISTANCE FROM LOWER BOUND:  XL-X =\n         ");
	for (i = 1; i <= *n; i++) {
		fprintf(stdout, "%16.8e", v[i - 1]);
	}
	fprintf(stdout, "\n");
	fprintf(stdout, "        DISTANCE FROM UPPER BOUND:  XU-X =\n         ");
	for (i = 1; i <= *n; i++) {
		fprintf(stdout, "%16.8e", w[i - 1]);
	}
	fprintf(stdout, "\n");
	if (!_n10ng->llise) {
		fprintf(stdout, "        NUMBER OF ITERATIONS:  ITER =%4ld\n",
			*iter);
	}
	fprintf(stdout, "        NUMBER OF FUNC-CALLS:  NFUNC =%4ld\n",
		*nfunc);
	fprintf(stdout, "        NUMBER OF GRAD-CALLS:  NGRAD =%4ld\n",
		*ngrad);
	fprintf(stdout, "        NUMBER OF QL-CALLS:    NQL   =%4ld\n\n\n\n",
		*nql);
	goto L_9000;
L_410:
	;
	/* CORRECT PENALTY PARAMETER */
	if (l7[4])
		goto L_420;
	wa[0] = *dbd;
	wa[1] = del[_n11ng->n1 - 1];
	wa[2] = _n10ng->rpenu;
	wa[3] = (Mfloat) (*iter);
        _l0 = imerit + 2;
        _l1 = 4;
	l_n5ong(&_l0, &_n10ng->m, &_n10ng->me, n, &_n11ng->mnn, &_n11ng->nmnn,
                &_n10ng->acc, rpen, &_n10ng->f, df, g, dg, lddg, vmu, u, x,
                xl, xu, phi, dphi, active, wa, &_l1);
	goto L_450;
L_420:
	sum = ze;
	for (i = 1; i <= *n; i++) {
		sum += dphi[i - 1] * del[i - 1] + fabs(u[_n10ng->m + i - 1] *
				 v[i - 1]) + fabs(u[mn + i - 1] * w[i - 1]);
	}
	if (sum > sqrt(sqacc))
		goto L_450;
	for (j = 1; j <= _n11ng->mnn; j++) {
		rpen[j - 1] = imsl_f_min(_n10ng->zefacu, rpen[j - 1] * _n10ng->rpeno);
	}
        _l0 = imerit + 4; 
        _l1 = 4;
	l_n5ong(&_l0, &_n10ng->m, &_n10ng->me, n, &_n11ng->mnn, &_n11ng->nmnn,
                &_n10ng->acc, rpen, &_n10ng->f, df, g, dg, lddg, vmu, u, x,
                xl, xu, phi, dphi, active, wa, &_l1);
L_450:
	if (_n10ng->iprint < 3)
		goto L_460;
	fprintf(stdout, "        PRODUCT OF SEARCH DIRECTION WITH BFGS-MATRIX:    DBD =%13.4e\n",
		*dbd);
	fprintf(stdout, "        PENALTY PARAMETER:  R =\n         ");
	for (j = 1; j <= _n11ng->mnn; j++) {
		fprintf(stdout, "%16.8e", rpen[j - 1]);
	}
	fprintf(stdout, "\n");
L_460:
	;
	/*
	 * EVALUATION OF MERIT FUNCTION
	 */
L_470:
        _l0 = imerit + 3; 
        _l1 = 4;
	l_n5ong(&_l0, &_n10ng->m, &_n10ng->me, n, &_n11ng->mnn, &_n11ng->nmnn,
                &_n10ng->acc, rpen, &_n10ng->f, df, g, dg, lddg, vmu, u, x,
                xl, xu, phi, dphi, active, wa, &_l1);
	if (!l7[4])
                _l0 = imerit + 4; 
                _l1 = 4;
		l_n5ong(&_l0, &_n10ng->m, &_n10ng->me, n, &_n11ng->mnn,
		        &_n11ng->nmnn, &_n10ng->acc, rpen, &_n10ng->f, df, g,
                        dg, lddg, vmu, u, x, xl, xu, phi, dphi, active, wa,
                        &_l1);
	*prd = imsl_sdot(*n, dphi, 1, del, 1);
	for (j = 1; j <= _n11ng->mnn; j++) {
		*prd += dphi[j + *n - 1] * (u[j - 1] - vmu[j - 1]);
	}
	if (*prd < ze)
		goto L_500;
	sscal(_n11ng->mnn, _n10ng->rpeno, rpen, 1);
        _l0 = 1;
	irpmax = l_ismax(_n11ng->mnn, rpen, _l0);
	rpmax = imsl_f_max(rpen[irpmax - 1], F_ZERO);
	if (rpmax < _n10ng->rpenu)
		goto L_470;
	if (l7[4])
		goto L_490;
	if (!l7[3] || *dbd < _n10ng->acc)
		goto L_180;
	dcl11 = *DCL(_n11ng->n1 - 1, _n11ng->n1 - 1);
	if (dcl11 >= _n10ng->zefacu)
		goto L_180;
	dcl11 *= _n10ng->rpeno;
	if (_n10ng->lql)
		dcl11 *= _n10ng->rpeno;
	*DCL(_n11ng->n1 - 1, _n11ng->n1 - 1) = dcl11;
	goto L_130;
L_490:
	;
	_n10ng->ifail = 2;
	if (_n10ng->iprint == 0)
		goto L_370;
	fprintf(stdout, "        **SEARCH DIRECTION NOT PROFITABLE:  DPHI*P =%13.4e\n",
		*prd);
	goto L_370;
L_500:
	;
	if (_n10ng->iprint < 3)
		goto L_510;
	fprintf(stdout, "        PRODUCT LAGRANGIAN GRADIENT WITH SEARCH DIRECTION:  DLP =%13.4e\n",
		*prd);
L_510:
	;
	/* LINE SEARCH */
	wa[5] = xnm;
	wa[6] = delnm;
	l7[6] = IMSLFALSE;
	*iflise = 0;
L_520:
	ipr = 0;
	if (_n10ng->iprint >= 1000)
		ipr = _n10ng->iprint - 1000;
        _l0 = 35;
        _l1 = 10;
        _l2 = 5;
	l_n8ong(alphao, alpham, phi, prd, &_n10ng->amue, &_n10ng->imsl_beta,
	        iline, &_n10ng->maxfun, iflise, &ipr, &wa[5], &_l0, iw, &_l1,
                &active[ilwls - 1], &_l2);
	if (*iflise > -2)
		goto L_540;
	l7[6] = IMSLTRUE;
	*fbest = _n10ng->f;
	scopy(_n10ng->m, g, 1, gbest, 1);
	if (_n10ng->llise)
		goto L_520;
	scopy(*n, df, 1, dfbest, 1);
	for (i = 1; i <= *n; i++) {
		scopy(_n10ng->m, DG(i - 1, 0), 1, DGBEST(i - 1, 0), 1);
	}
	goto L_520;
L_540:
	;
	for (i = 1; i <= *n; i++) {
		x[i - 1] = xold[i - 1] + *alphao * del[i - 1];
	}
	for (j = 1; j <= _n11ng->mnn; j++) {
		vmu[j - 1] = vmuold[j - 1] + *alphao * (u[j - 1] - vmuold[j - 1]);
	}
	if (*iflise == 0)
		goto L_590;
	if (*iflise == 1)
		goto L_580;
	if (*iflise > 1)
		goto L_570;
	goto L_620;
L_570:
	_n10ng->ifail = 1000 + *iflise;
	if (_n10ng->iprint == 0)
		goto L_370;
	fprintf(stdout, "        **ERROR IN LINE SEARCH. IFLISE =%4ld\n",
		*iflise);
	goto L_370;
L_580:
	_n10ng->ifail = 4;
	if (_n10ng->iprint == 0)
		goto L_370;
	fprintf(stdout, "        **MORE THAN MAXFUN FUNC-CALLS IN LINE SEARCH\n");
	goto L_370;
L_590:
	l7[2] = IMSLTRUE;
	if (_n10ng->iprint < 3)
		goto L_600;
	if (*iline == 1) {
		fprintf(stdout, "        LINE SEARCH SUCCESSFUL AFTER ONE STEP:  ALPHA = 1.\n");
	}
	if (*iline > 1) {
		fprintf(stdout, "        LINE SEARCH SUCCESSFUL AFTER%3ld STEPS:  ALPHA =%13.4e\n",
			*iline, *alphao);
	}
L_600:
	;
	if (!l7[6] && _n10ng->llise)
		goto L_640;
	if (!l7[6] && !_n10ng->llise)
		goto L_240;
	_n10ng->f = *fbest;
	scopy(_n10ng->m, gbest, 1, g, 1);
	if (_n10ng->llise)
		goto L_640;
	scopy(*n, dfbest, 1, df, 1);
	for (i = 1; i <= *n; i++) {
		scopy(_n10ng->m, DGBEST(i - 1, 0), 1, DG(i - 1, 0), 1);
	}
        _l0 = imerit + 1;
        _l1 = 4;
	l_n5ong(&_l0, &_n10ng->m, &_n10ng->me, n, &_n11ng->mnn, &_n11ng->nmnn,
                &_n10ng->acc, rpen, &_n10ng->f, df, g, dg, lddg, vmu, u, x,
                xl, xu, phi, dphi, active, wa, &_l1);
	goto L_240;
L_620:
	;
	/* NEW FUNCTION AND GRADIENT VALUES */
	if (l7[5]) {
		_n10ng->ifail = -1;
		goto L_9000;
	}
	imsl_e1usr("ON");
	(*fcns) (_n10ng->m, _n10ng->me, *n, x, &active[*mmax], &_n10ng->f,
		 g);
	imsl_e1usr("OFF");
L_630:
	;
	if (l7[0])
		_n10ng->f *= *scf;
	if (_n10ng->m != 0 && l7[1]) {
                _l0 = 1;
		l_shprod(&_n10ng->m, g, &_l0, scg, &_l1, g, &_l2);
	}
	*nfunc += 1;
        _l0 = imerit + 3; 
        _l1 = 4; 
	l_n5ong(&_l0, &_n10ng->m, &_n10ng->me, n, &_n11ng->mnn, &_n11ng->nmnn,
                &_n10ng->acc, rpen, &_n10ng->f, df, g, dg, lddg, vmu, u, x,
                xl, xu, phi, dphi, active, wa, &_l1);
	if (_n10ng->llise && !l7[2])
		goto L_520;
L_640:
	;
        _l0 = imerit + 1; 
        _l1 = 4; 
	l_n5ong(&_l0, &_n10ng->m, &_n10ng->me, n, &_n11ng->mnn, &_n11ng->nmnn,
                &_n10ng->acc, rpen, &_n10ng->f, df, g, dg, lddg, vmu, u, x,
                xl, xu, phi, dphi, active, wa, &_l1);
	if (l7[0])
		_n10ng->f /= *scf;
	if (_n10ng->m != 0 && l7[1]) {
		for (j = 1; j <= _n10ng->m; j++) {
			g[j - 1] /= scg[j - 1];
		}
	}
	if (l7[5]) {
		_n10ng->ifail = -2;
		goto L_9000;
	}
	imsl_e1usr("ON");
	(*grad) (_n10ng->m, _n10ng->me, *mmax, *n, x, active, _n10ng->f, g,
                 df, dgt);
	imsl_e1usr("OFF");
        
        for (i=0; i<*n; i++) {
             k = i * *mmax;
             scopy (*mmax, &dgt[i], *n, &dg[k], 1);
        }

L_660:
	;
	*ngrad += 1;
	if (!l7[0])
		goto L_670;
	_n10ng->f *= *scf;
	sscal(*n, *scf, df, 1);
L_670:
	if (_n10ng->m == 0 || !l7[1])
		goto L_690;
	for (j = 1; j <= _n10ng->m; j++) {
		g[j - 1] *= scg[j - 1];
		if (active[j - 1])
			sscal(*n, scg[j - 1], DG(0, j - 1), *lddg);
	}
L_690:
	;
	if (l7[2])
		goto L_240;
        _l0 = imerit + 4;
        _l1 = 4;
	l_n5ong(&_l0, &_n10ng->m, &_n10ng->me, n, &_n11ng->mnn, &_n11ng->nmnn,
                &_n10ng->acc, rpen, &_n10ng->f, df, g, dg, lddg, vmu, u, x, 
                xl, xu, phi, dphi, active, wa, &_l1);
	*prd = imsl_sdot(*n, dphi, 1, del, 1);
	for (j = 1; j <= _n11ng->mnn; j++) {
		*prd += dphi[*n + j - 1] * (u[j - 1] - vmuold[j - 1]);
	}
	goto L_520;
	/* UPDATE HESSIAN OF LAGRANGIAN */
L_710:
	*dbd *= *alphao ** alphao;
	sscal(*n, *alphao, bdel, 1);
	for (i = 1; i <= *n; i++) {
		eta[i - 1] = dla[i - 1] - dlaold[i - 1];
	}
	edel = *alphao * imsl_sdot(*n, del, 1, eta, 1);
	dbd1 = _n10ng->dbdfac ** dbd;

	if (edel < dbd1) {
		theta = (*dbd - dbd1) / (*dbd - edel);
		theta1 = on - theta;
		for (i = 1; i <= *n; i++) {
			eta[i - 1] = theta * eta[i - 1] + theta1 * bdel[i - 1];
		}
		edel = dbd1;
	}
	dbdi = sqrt(on / *dbd);
	edeli = sqrt(on / edel);
	/* UPDATE FACTORIZATION */
	sscal(*n, dbdi, bdel, 1);
	sscal(*n, edeli, eta, 1);
	if (!_n10ng->lql)
		goto L_750;
	for (i = 1; i <= *n; i++) {
		for (j = 1; j <= i; j++) {
			*DCL(i - 1, j - 1) += eta[i - 1] * eta[j - 1] - bdel[i - 1] *
				bdel[j - 1];
		}
	}
	imsl_csfrg(n, dcl, lddcl);
	goto L_60;
L_750:
	;
	l_n7ong(n, dcl, lddcl, cd, eta, bdel);
	/* CORRECT DATA FOR QL-SOLUTION */
L_760:
	for (i = 1; i <= *n; i++) {
		sqd = sqrt(cd[i - 1]);
		if (sqd > uf)
			goto L_770;
		_n10ng->ifail = 3;
		if (_n10ng->iprint == 0)
			goto L_370;
		fprintf(stdout, "        **UNDERFLOW IN BFGS-UPDATE\n");
		goto L_370;
L_770:
		;
		if (i <= *n - 1)
                        _l0 = *n - i;
                        _l1 = 1;
			imsl_svcal(_l0, sqd, DCL(i - 1, i), _l1,
				   DCL(i, i - 1), *lddcl);
		sset(*n - i, sqd, DCL(0, 0), *lddcl + 1);
	}
	if (*iter == 0)
		goto L_50;
	/* PERFORM NEXT ITERATION */
	goto L_60;
L_9000:
        if (dgt != NULL) imsl_free (dgt);

	return;
}				/* end of function */
#undef DG
#undef DCL
#undef DGBEST
/*----------------------------------------------------------------------- */

/*  IMSL Name:  N5ONF/DN5ONF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 1, 1985

    Purpose:    Forward difference approximation for the gradient of the
                objective function and the constraMints.

    Usage:      CALL N5ONF (FCN, M, ME, MMAX, N, XC, XSCALE, ACTIVE, FC,
                            GC, DF, DG, WORK)

    Arguments:
       FCN    - User-supplied SUBROUTINE to evaluate the functions at
                a given poMint.  The usage is
                CALL FCN (M, ME, N, X, ACTIVE, F, G), where
                M      - Total number of constraMints.  (Input)
                ME     - Number of equality constraMints.  (Input)
                N      - Number of variables.  (Input)
                X      - The poMint at which the function is evaluated.
                         (Input)
                         X should not be changed by FCN.
                ACTIVE - Logical vector of length MMAX indicating the
                         active constraMints.  (Input)
                F      - The computed function value at the poMint X.
                         (Output)
                G      - Vector of length MMAX containing the values of
                         constraMints at poMint X.  (Output)
                FCN must be declared EXTERNAL in the calling program.
       M      - Number of constraMints.  (Input)
       ME     - Number of equality constraMints.  (Input)
       MMAX   - Leading dimension of DG.  (Input)
                MMAX = MAX(1,M)
       N      - Number of variables.  (Input)
       XC     - Vector of length N containing the current poMint.  (Input)
       XSCALE - Vector of length N containing the diagonal scaling
                matrix.  (Input)
       ACTIVE - Logical vector of length M determining which constraMints
                are active.  (Input)
       FC     - Current objective function value.  (Input)
       GC     - Vector of length M containing the current constraMint
                values.  (Input)
       DF     - Vector of length M containing the estimate of the
                gradient of the objective function.  (Output)
       DG     - Matrix of dimension MMAX by N containing the estimates of
                the gradients of the constraMints.  (Output)
       WORK   - Work vector of length M.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
#if defined (COMPUTER_HP97C)
static void l_n5onf(void (*fcn)(Mint, Mint, Mint, Mfloat*, Mint*, Mfloat*,
                    Mfloat*), Mint *m, Mint *me, Mint *mmax, Mint *n,
                    Mfloat xc[], Mfloat xscale[], Mint active[], Mfloat *fc,
                    Mfloat gc[], Mfloat df[], Mfloat dg[], Mfloat work[])
#else
static void l_n5onf(void (*fcn)(Mint, Mint, Mint, Mfloat[], Mint[], Mfloat*,
                    Mfloat[]), Mint *m, Mint *me, Mint *mmax, Mint *n,
                    Mfloat xc[], Mfloat xscale[], Mint active[], Mfloat *fc,
                    Mfloat gc[], Mfloat df[], Mfloat dg[], Mfloat work[])
#endif
#else
static void l_n5onf(fcn, m, me, mmax, n, xc, xscale, active, fc, gc, df, dg,
                    work)
	void            (*fcn) ();
	Mint            *m, *me, *mmax, *n;
	Mfloat           xc[], xscale[];
	Mint            active[];
	Mfloat          *fc, gc[], df[], dg[], work[];
#endif
{
#define DG(I_,J_)	(dg+(I_)*(*mmax)+(J_))
	Mint             i, j;
	Mfloat           eps, fnew, stepsz, xtempj;


	eps = sqrt(imsl_amach(4));
	for (j = 1; j <= *n; j++) {
		stepsz = eps * imsl_f_max(fabs(xc[j-1]), F_ONE / xscale[j-1]);
		if (xc[j-1] < F_ZERO)
			stepsz = -stepsz;
		xtempj = xc[j-1];
		xc[j-1] = xtempj + stepsz;
		imsl_e1usr("ON");
		(*fcn) (*m, *me, *n, xc, active, &fnew, work);
		imsl_e1usr("OFF");
		xc[j-1] = xtempj;
		df[j-1] = (fnew - *fc) / stepsz;
		for (i = 1; i <= *m; i++) {
			if (active[i-1])
				*DG(j-1, i-1) = (work[i-1] - gc[i-1]) / stepsz;
		}
	}

	return;
}				/* end of function */
#undef DG
/*----------------------------------------------------------------------- */

/*  IMSL Name:  N5ONG/DN5ONG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    December 2, 1985

    Purpose:    Evaluate the function and gradient of the augmented
                Lagrangian merit function, compute penalty parameters or
                set active constraMints.

    Usage:      CALL N5ONG (MODE, M, ME, MMAX, N, MNN, NMNN, ACC, R, F,
                            DF, G, DG, LDDG, V, U, X, XL, XU, PHI,
                            ACTIVE, WA, LWA)

    Arguments:
       MODE   - Parameter determining type of computation.  (Input)
       M      - Number of constraMints.  (Input)
       ME     - Number of equality constraMints.  (Input)
       MMAX   - Leading dimension of DG.  (Input)
                MMAX = MAX(1,M)
       N      - Number of variables.  (Input)
       MNN    - Scalar variable such that MNN = M + 2*N.  (Input)
       NMNN   - Scalar variable such that NMNN = M + 3*N.  (Input)
       ACC    - Final accuracy.  (Input)
       R      - Vector of length MNN containing the penalty function
                parameters.  (Input/Output)
       F      - Current function value.  (Input)
       DF     - Vector of length MAX(2,N+1) containing the current
                gradient of the objective function.  (Input)
       G      - Vector of length MMAX containing current constraMint
                values.  (Input)
       DG     - Array of dimension MMAX by MMAX containing the current
                gradient of the constraMints.  (Input)
       LDDG   - Leading dimension of DG exactly as specified in the
                dimension statement of the calling program.  (Input)
       V      - Vector of length MNN.  (Input)
       U      - Vector of length MNN containing the multipliers of the
                nonlinear constraMints and the bounds.  (Input)
       X      - Vector of length N containing the current poMint being
                evaluated.  (Input)
       XL     - Vector of length N containing the lower bounds for the
                variables.  (Input)
       XU     - Vector of length N containing the upper bounds for the
                variables.  (Input)
       PHI    - Scalar variable.  (Output)
       DPHI   - Vector of length NMNN.
       ACTIVE - Logical vector of length MMAX containing information
                as to which constraMints are active.  (Output)
       WA     - Work vector of length LWA.
       LWA    - Length of WA where LWA = 4.  (Input)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_n5ong(Mint *mode, Mint *m, Mint *me, Mint *n, Mint *mnn,
                    Mint *nmnn, Mfloat *acc, Mfloat r[], Mfloat *f, Mfloat df[],
                    Mfloat g[], Mfloat dg[], Mint *lddg, Mfloat v[], Mfloat u[],
                    Mfloat x[], Mfloat xl[], Mfloat xu[], Mfloat *phi, 
                    Mfloat dphi[], Mint active[], Mfloat wa[], Mint *lwa)
#else
static void l_n5ong(mode, m, me, n, mnn, nmnn, acc, r, f, df, g, dg, lddg,
                    v, u, x, xl, xu, phi, dphi, active, wa, lwa)
	Mint            *mode, *m, *me, *n, *mnn, *nmnn;
	Mfloat          *acc, r[], *f, df[], g[], dg[];
	Mint            *lddg;
	Mfloat           v[], u[], x[], xl[], xu[], *phi, dphi[];
	Mint            active[];
	Mfloat           wa[];
	Mint            *lwa;
#endif
{
#define DG(I_,J_)	(dg+(I_)*(*lddg)+(J_))
	Mint            index[1000];
        Mint            IMSLFALSE = 0, IMSLTRUE = 1;
	Mint             i, j, me1, mn;
	Mfloat           eps, fact, phi1, rmi, rmj, rmni, rmnj, rnew, sum,
	                uad, uf, vmi, vmni, xmxl, xumx;


	/*
	 * CONSTANT DATA
	 */
	eps = imsl_amach(4);
	uf = eps * eps;
	if (*mode == 2)
		goto L_30;
	if (*mode == 3 || *mode == 4)
		goto L_50;
	if (*mode == 5)
		goto L_130;
	if (*mode == 6)
		goto L_150;
	if (*mode == 7)
		goto L_170;
	if (*mode == 8)
		goto L_190;

	/*
	 * EVALUATION OF ACTIVE SET
	 */

	if (*m == 0 || *me == *m)
		goto L_9000;
	me1 = *me + 1;
	for (j = me1; j <= *m; j++) {
		active[j - 1] = IMSLTRUE;
		if (v[j - 1] == F_ZERO && g[j - 1] > *acc)
			active[j - 1] = IMSLFALSE;
	}
	goto L_9000;

	/*
	 * EVALUATION OF PENALTY PARAMETER
	 */
L_30:
	fact = wa[2];
	if (wa[0] > uf)
		fact = (Mfloat) (2 ** m) / (wa[0] + (F_ONE - wa[1]));
	for (j = 1; j <= *mnn; j++) {
		rnew = imsl_f_min(wa[2], fact * imsl_fi_power(u[j - 1] - v[j - 1], 2));
		uad = imsl_f_min(r[j - 1], wa[3] * sqrt(r[j - 1]));
		r[j - 1] = imsl_f_max(rnew, uad);
	}
	goto L_9000;

	/*
	 * PREPARE EVALUATION OF AUGMENTED LAGRANGIAN
	 */
L_50:
	for (j = 1; j <= *m; j++) {
		index[j - 1] = IMSLTRUE;
		if (j > *me && g[j - 1] > v[j - 1] / r[j - 1])
			index[j - 1] = IMSLFALSE;
	}
	if (*mode != 4) {

		/*
		 * EVALUATION OF AUGMENTED LAGRANGIAN VALUE
		 */
		*phi = *f;
		for (j = 1; j <= *m; j++) {
			if (index[j - 1])
				phi1 = (v[j - 1] - 0.5e0 * r[j - 1] * g[j - 1]) * g[j - 1];
			if (!index[j - 1])
				phi1 = 0.5e0 * v[j - 1] * v[j - 1] / r[j - 1];
			*phi -= phi1;
		}
		mn = *m + *n;
		for (i = 1; i <= *n; i++) {
			rmj = r[*m + i - 1];
			rmi = F_ONE / rmj;
			rmnj = r[mn + i - 1];
			rmni = F_ONE / rmnj;
			vmi = v[*m + i - 1];
			vmni = v[mn + i - 1];
			xmxl = x[i - 1] - xl[i - 1];
			xumx = xu[i - 1] - x[i - 1];
			if (xmxl <= vmi * rmi)
				*phi += -(vmi - 0.5e0 * rmj * xmxl) * xmxl;
			if (xmxl > vmi * rmi)
				*phi += -0.5e0 * rmi * vmi * vmi;
			if (xumx <= vmni * rmni)
				*phi += -(vmni - 0.5e0 * rmnj * xumx) * xumx;
			if (xumx > vmni * rmni)
				*phi += -0.5e0 * rmni * vmni * vmni;
		}
		goto L_9000;
	}
	/*
	 * GRADIENT EVALUATION OF AUGMENTED LAGRANGIAN
	 */
	for (i = 1; i <= *n; i++) {
		dphi[i - 1] = df[i - 1];
		sum = F_ZERO;
		for (j = 1; j <= *m; j++) {
			if (index[j - 1])
				sum += *DG(i - 1, j - 1) * (v[j - 1] - r[j - 1] * g[j - 1]);
		}
		dphi[i - 1] -= sum;
	}

	for (j = 1; j <= *m; j++) {
		if (index[j - 1])
			dphi[*n + j - 1] = -g[j - 1];
		if (!index[j - 1])
			dphi[*n + j - 1] = -v[j - 1] / r[j - 1];
	}
	mn = *m + *n;
	for (i = 1; i <= *n; i++) {
		vmi = v[*m + i - 1];
		vmni = v[mn + i - 1];
		xmxl = x[i - 1] - xl[i - 1];
		xumx = xu[i - 1] - x[i - 1];
		rmj = r[*m + i - 1];
		rmi = F_ONE / rmj;
		rmnj = r[mn + i - 1];
		rmni = F_ONE / rmnj;
		if (xmxl <= vmi * rmi)
			dphi[i - 1] += -vmi + rmj * xmxl;
		if (xumx <= vmni * rmni)
			dphi[i - 1] += vmni - rmnj * xumx;
		if (xmxl <= vmi * rmi)
			dphi[mn + i - 1] = -xmxl;
		if (xmxl > vmi * rmi)
			dphi[mn + i - 1] = -vmi * rmi;
		if (xumx <= vmni * rmni)
			dphi[*mnn + i - 1] = -xumx;
		if (xumx > vmni * rmni)
			dphi[*mnn + i - 1] = -vmni * rmni;
	}
	goto L_9000;

	/*
	 * EVALUATION OF L1-PENALTY FUNCTION
	 * 
	 * ACTIVE SET
	 */
L_130:
	for (j = 1; j <= *m; j++) {
		active[j - 1] = IMSLTRUE;
	}
	goto L_9000;

	/*
	 * PENALTY PARAMETER
	 */
L_150:
	for (j = 1; j <= *mnn; j++) {
		phi1 = fabs(u[j - 1]);
		r[j - 1] = imsl_f_max(phi1, 0.5e0 * (r[j - 1] + phi1));
	}
	goto L_9000;

	/*
	 * MERIT FUNCTION VALUE
	 */
L_170:
	*phi = *f;
	for (j = 1; j <= *m; j++) {
		phi1 = fabs(g[j - 1]);
		if (j > *me && g[j - 1] > F_ZERO)
			phi1 = F_ZERO;
		*phi += r[j - 1] * phi1;
	}

	/*
	 * GRADIENT OF MERIT FUNCTION
	 */
L_190:
	sset(*nmnn - *n + 1, F_ZERO, &dphi[*n - 1], 1);
	for (i = 1; i <= *n; i++) {
		dphi[i - 1] = df[i - 1];
		if (*m != 0) {
			phi1 = F_ZERO;
			for (j = 1; j <= *m; j++) {
				if (g[j - 1] > F_ZERO && j <= *me)
					phi1 += r[j - 1] ** DG(i - 1, j - 1);
				if (g[j - 1] <= F_ZERO)
					phi1 += -r[j - 1] ** DG(i - 1, j - 1);
			}
			dphi[i - 1] += phi1;
		}
	}

L_9000:
	return;
}				/* end of function */
#undef DG
/*----------------------------------------------------------------------- */

/*  IMSL Name:  N6ONG/DN6ONG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 24, 1987

    Purpose:    Solve the quadratic programming problem.

    Usage:      CALL N6ONG (M, ME, MMAX, N, NMAX, C, LDC, D, A, LDA, B,
                            XL, XU, X, U, IFAIL, IPRINT, WAR, LWAR,
                            IWAR, LIWAR)

    Arguments:
       M      - Number of constraMints.  (Input)
       ME     - Number of equality constraMints.  (Input)
       MMAX   - Leading dimension of A.  (Input)
                MMAX must be at least MAX(1,M).
       N      - Number of variables.  (Input)
       NMAX   - Leading dimension of C.  (Input)
                NMAX must be at least MAX(2,N+1).
       C      - Array of dimension NMAX by NMAX containing the objective
                function matrix.  (Output)
       LDC    - Leading dimension of C exactly as specified in the
                dimension statement of the calling program.  (Input)
       D      - Vector of length NMAX containing the objective function
                gradient.  (Input)
       A      - Array of dimension MMAX by MMAX containing the linear
                constraMints.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       B      - Vector of length MMAX containing the constant data of the
                linear constraMints.  (Input)
       XL     - Vector of length N containing the lower bounds for the
                variables.  (Input)
       XU     - Vector of length N containing the upper bounds for the
                variables.  (Input)
       X      - Vector of length N containing the current poMint being
                evaluated.  (Input)
       U      - Vector of length M + 2*(N+1) containing the multipliers
                of the nonlinear constraMints and the bounds.  (Output)
       IFAIL  - Scalar containing error message information.  (Output)
       IPRINT - Specification for the desired print level.  (Input)
       WAR    - Work vector of length LWA.
       LWAR   - Length of WA where LWA = N*(2*N+13) + MMAX + M + 12.
                (Input)
       IWAR   - Work vector of length LKWA.
       LIWAR  - Length of KWA where LKWA = MAX(N,M).  (Input)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* LIWAR is not used here, but leave the calling sequence */
#ifdef ANSI
static void l_n6ong(Mint *m, Mint *me, Mint *mmax, Mint *n, Mint *nmax,
                    Mint *mnn, Mfloat *c, Mint *ldc, Mfloat d[], Mfloat *a,
                    Mint *lda, Mfloat b[], Mfloat xl[], Mfloat xu[], Mfloat x[],
                    Mfloat u[], Mint *ifail, Mint *iprint, Mfloat war[],
                    Mint *lwar, Mint iwar[], Mint *liwar)
#else
static void l_n6ong(m, me, mmax, n, nmax, mnn, c, ldc, d, a, lda, b, xl, xu,
                    x, u, ifail, iprint, war, lwar, iwar, liwar)
	Mint            *m, *me, *mmax, *n, *nmax, *mnn;
	Mfloat          *c;
	Mint            *ldc;
	Mfloat           d[], *a;
	Mint            *lda;
	Mfloat           b[], xl[], xu[], x[], u[];
	Mint            *ifail, *iprint;
	Mfloat           war[];
	Mint            *lwar, iwar[], *liwar;
#endif
{
#define C(I_,J_)	(c+(I_)*(*ldc)+(J_))
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
	Mint            lql;
        Mint            IMSLFALSE = 0, IMSLTRUE = 1;
        FILE            *iout;
	Mint             _l0, i, idiag, in, info, inw1, inw2, 
	                lw, mn, nact;
	Mfloat           _f0, diag;

	/*
	 * CONSTANT DATA
	 */
	info = 0;
	lql = IMSLFALSE;
	if (iwar[0] == 1)
		lql = IMSLTRUE;
	imsl_umach(2, &iout);
	inw1 = 1;
	inw2 = inw1 + *mmax;

	/*
	 * PREPARE PROBLEM DATA FOR EXECUTION
	 */
        _l0 = 1;
        _f0 = -F_ONE;
	imsl_svcal(*m, _f0, b, _l0, &war[inw1 - 1], _l0);
	lw = 3 ** nmax ** nmax / 2 + 10 ** nmax + *m;
	if ((inw2 + lw) <= *lwar) {
		mn = *m + *n;

		/*
		 * CALL OF N9ONG
		 */
		l_n9ong(n, m, me, mmax, &mn, mnn, nmax, &lql, a, lda, &war[inw1 - 1],
			   d, c, ldc, xl, xu, x, &nact, &iwar[0], &info, &diag, &war[inw2 - 1],
			   &lw);

		/*
		 * TEST OF MATRIX CORRECTIONS
		 */
		*ifail = 0;
		if (lql) {
			idiag = diag;
			if ((*iprint > 0) && (diag > F_ZERO)) {
				fprintf(stdout, "\n        ***QL: MATRIX G WAS ENLARGED%3ld-TIMES BY UNITMATRIX\n",
					idiag);
			}
		}
		/*
		 * REORDER MULTIPLIER
		 */
		if (info >= 0) {
			sset(*mnn, F_ZERO, u, 1);
			in = inw2 - 1;
			for (i = 1; i <= nact; i++) {
				u[iwar[i - 1] - 1] = war[in + i - 1];
			}
			goto L_9000;
		}
	}
	/*
	 * ERROR MESSAGES
	 */
	*ifail = -info + 10;
	if ((*iprint > 0) && (nact > 0)) {
		fprintf(stdout, "\n        ***QL: CONDITION %2ld", info);
		fprintf(stdout, " NOT CONSISTENT TO \n          ");
		for (i = 1; i <= nact; i++) {
			fprintf(stdout, "%12e", war[inw2 + i - 2]);
		}
		fprintf(stdout, "\n");
	}
L_9000:
	return;
}				/* end of function */
#undef A
#undef C
/*----------------------------------------------------------------------- */

/*  IMSL Name:  N7ONG/DN7ONG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    December 2, 1985

    Purpose:    Update the Cholesky factorization, L*D*TRANS(L) =
                L*D*TRANS(L) + U*TRANS(U) - V*TRANS(V)

    Usage:      CALL N7ONG (N, CL, LDCL, D, U, V)

    Arguments:
       N      - Number of variables.  (Input)
       CL     - Matrix of dimension MAX(2,N+1) by N containing the
                approximate Hessian.  (Input/Output)
       LDCL   - Leading dimension of CL exactly as specified in the
                dimension statement of the calling program.  (Input)
       D      - Vector of length N containing the diagonal elements of
                the Hessian.  (Input/Output)
       U      - Vector of length N.  (Input/Output)
       V      - Vector of length N.  (Input/Output)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_n7ong(Mint *n, Mfloat *cl, Mint *ldcl, Mfloat d[], Mfloat u[],
                    Mfloat v[])
#else
static void l_n7ong(n, cl, ldcl, d, u, v)
	Mint            *n;
	Mfloat          *cl;
	Mint            *ldcl;
	Mfloat           d[], u[], v[];
#endif
{
#define CL(I_,J_)	(cl+(I_)*(*ldcl)+(J_))
	Mint             i, j, j1, k, n1;
	Mfloat           imsl_beta, clkj, dj, djnew, eps, pdp, pj, t, tnew,
	                tt, tti, ui;


	eps = imsl_amach(4);
	/* UPDATE L*D*L(T) + U*U(T) */
	t = F_ONE;
	for (j = 1; j <= *n; j++) {
		pj = u[j - 1];
		dj = d[j - 1];
		tnew = t + pj * pj / dj;
		tt = tnew / t;
		tti = F_ONE / tt;
		djnew = dj * tt;
		imsl_beta = pj / (dj * tnew);
		j1 = j + 1;
		if (j1 <= *n) {
			if (djnew > F_FOUR * dj) {
				for (k = j1; k <= *n; k++) {
					clkj = *CL(j - 1, k - 1);
					*CL(j - 1, k - 1) = tti * clkj + imsl_beta * u[k - 1];
					u[k - 1] += -pj * clkj;
				}
			} else {
				saxpy(*n - j, -pj, CL(j - 1, j1 - 1), 1, &u[j1 - 1],
					   1);
				saxpy(*n - j, imsl_beta, &u[j1 - 1], 1, CL(j - 1, j1 - 1),
					   1);
			}
		}
		t = tnew;
		d[j - 1] = djnew;
	}

	/*
	 * UPDATE L*D*L(T) - V*V(T)
	 */
	u[0] = v[0];
	pdp = u[0] * u[0] / d[0];
	for (i = 2; i <= *n; i++) {
		ui = v[i - 1] - imsl_sdot(i - 1, CL(0, i - 1), *ldcl, u, 1);
		pdp += ui * ui / d[i - 1];
		u[i - 1] = ui;
	}

	t = F_ONE - pdp;
	if (t <= eps)
		t = eps;

	n1 = *n + 1;
	for (i = 1; i <= *n; i++) {
		j = n1 - i;
		pj = u[j - 1];
		dj = d[j - 1];
		tnew = t + pj * pj / dj;
		tt = t / tnew;
		djnew = dj * tt;
		imsl_beta = -pj / (dj * t);
		v[j - 1] = pj;
		t = tnew;
		d[j - 1] = djnew;
		j1 = j + 1;
		for (k = j1; k <= *n; k++) {
			clkj = *CL(j - 1, k - 1);
			*CL(j - 1, k - 1) = clkj + imsl_beta * v[k - 1];
			v[k - 1] += pj * clkj;
		}
	}

	return;
}				/* end of function */
#undef CL
/*----------------------------------------------------------------------- */

/*  IMSL Name:  N8ONG/DN8ONG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    December 2, 1985

    Purpose:    Line search based on quadratic Minterpolation and
                Armijo-type stopping criterion.

    Usage:      CALL N8ONG (ALPHA, ALPHAM, PHI, DPHI, AMUE, BETA, ILINE,
                            MAXFUN, IFAIL, IPRINT, WA, LWA, KWA, LKWA,
                            LOWA, LLOWA)

    Arguments:
       ALPHA  - Scalar variable.  (Output)
       ALPHAM - Scalar variable.  (Input)
       PHI    - Scalar variable.  (Input)
       DPHI   - Scalar variable.  (Input/Output)
       AMUE   - Scalar variable.  (Input)
       BETA   - Scalar variable.  (Input)
       ILINE  - Number of line search iterations.  (Output)
       MAXFUN - Maximum number of function evaluations allowed.  (Input)
       IFAIL  - Error indicator.  (Input/Output)
       IPRINT - PrMinting parameter.  (Input)
       WA     - Work vector of length LWA.
       LWA    - Length of vector WA where LWA = 35.  (Input)
       KWA    - Work vector of length LKWA.
       LKWA   - Length of vector KWA where LKWA = 10.  (Input)
       LOWA   - Logical vector of length LLOWA.
       LLOWA  - Length of vector LOWA where LLOWA = 5.  (Input)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* IPRINT, LWA, KWA, LKWA, LOWA, LLOWA are not used here, but leave 
    the calling sequence intact. */
#ifdef ANSI
static void l_n8ong(Mfloat *alpha, Mfloat *alpham, Mfloat *phi, Mfloat *dphi,
                    Mfloat *amue, Mfloat *imsl_beta, Mint *iline, Mint *maxfun,
                    Mint *ifail, Mint *iprint, Mfloat wa[], Mint *lwa,
                    Mint kwa[], Mint *lkwa, Mint lowa[], Mint *llowa)
#else
static void l_n8ong(alpha, alpham, phi, dphi, amue, imsl_beta, iline, maxfun,
                    ifail, iprint, wa, lwa, kwa, lkwa, lowa, llowa)
	Mfloat          *alpha, *alpham, *phi, *dphi, *amue, *imsl_beta;
	Mint            *iline, *maxfun, *ifail, *iprint;
	Mfloat           wa[];
	Mint            *lwa, kwa[], *lkwa;
	Mint            lowa[];
	Mint            *llowa;
#endif
{
	Mfloat           imsl_diff;

	/* START LINE SEARCH */
	if (*ifail == -1)
		goto L_20;
	wa[0] = *phi;
	*iline = 0;
	*alpha = *alpham;
	wa[1] = *alpha;
L_10:
	*iline += 1;
	*dphi *= wa[1];
	*ifail = -1;
	goto L_9000;
L_20:
	if (*iline >= *maxfun) {
		*ifail = 1;
		goto L_9000;
	}
	imsl_diff = *phi - wa[0];
	if (imsl_diff <= *amue ** dphi) {
		*ifail = 0;
		goto L_9000;
	}
	wa[1] = imsl_f_max(*imsl_beta, 0.5e0 ** dphi / (*dphi - imsl_diff));
	*alpha *= wa[1];
	goto L_10;

L_9000:
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  N9ONG/DN9ONG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 24, 1987

    Purpose:    Compute minimum of the unconstrained problem.

    Usage:      CALL N9ONG (N, M, MEQ, MMAX, MN, MNN, NMAX, LQL, A, B,
                            GRAD, G, XL, XU, X, NACT, IACT, INFO, DIAG,
                            W, LW)

    Arguments:
       N      - Number of variables.  (Input)
       M      - Number of constraMints.  (Input)
       MEQ    - Number of equality constraMints.  (Input)
       MMAX   - Leading dimension of A.  (Input)
                MMAX must be at least MAX(1,M).
       MN     - Scalar variable suxh that MN = M + N.  (Input)
       MNN    - Scalar variable suxh that MNN = M + 2*N.  (Input)
       NMAX   - Leading dimension of G.  (Input)
                NMAX must be at least MAX(2,N).
       LQL    - Logical scalar determining the initial decomposition.
                (Input)
                If LQL is true, the initial Cholesky-factorization of G
                is performed.  If LQL is false, the upper triangle of G
                contains the Cholesky-factor of a suitable decomposition.
       A      - Array of dimension MMAX by NMAX containing the constraMint
                normals in the columns.  (Output)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       B      - Vector of length MMAX containing the right-hand-sides of
                the constraMints.  (Input)
       GRAD   - Vector of length N containing the objective function
                gradient.  (Input)
       G      - Array of dimension NMAX by N containing symmetric
                objective function matrix.  (Input)
       XL     - Vector of length N containing the lower bounds for the
                variables.  (Input)
       XU     - Vector of length N containing the upper bounds for the
                variables.  (Input)
       X      - Vector of length N containing the current poMint being
                evaluated.  (Input)
       NACT   - Number of active constraMints.  (Output)
       IACT   - Vector of length NACT indicating the final active
                constraMints.   (Output)
       INFO   - Scalar containing exiting information.  (Output)
       DIAG   - Scalar containing multiple of the unit matrix that was
                added to G to achieve positive definiteness.  (Output)
       W      - Work vector of length LW.
       LW     - Length of W where LW = NMAX*(2*NMAX+10) + M.
                (Input)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* MMAX, MNN, and LW are not used, but leave the calling sequence intact.*/
#ifdef ANSI
static void l_n9ong(Mint *n, Mint *m, Mint *meq, Mint *mmax, Mint *mn,
                    Mint *mnn, Mint *nmax, Mint *lql, Mfloat *a, Mint *lda,
                    Mfloat b[], Mfloat grad[], Mfloat *g, Mint *ldg, 
                    Mfloat xl[], Mfloat xu[], Mfloat x[], Mint *nact,
                    Mint iact[], Mint *info, Mfloat *diag, Mfloat w[],
                    Mint *lw)
#else
static void l_n9ong(n, m, meq, mmax, mn, mnn, nmax, lql, a, lda, b, grad, g,
                    ldg, xl, xu, x, nact, iact, info, diag, w, lw)
	Mint            *n, *m, *meq, *mmax, *mn, *mnn, *nmax;
	Mint           *lql;
	Mfloat          *a;
	Mint            *lda;
	Mfloat           b[], grad[], *g;
	Mint            *ldg;
	Mfloat           xl[], xu[], x[];
	Mint            *nact, iact[], *info;
	Mfloat          *diag, w[];
	Mint            *lw;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define G(I_,J_)	(g+(I_)*(*ldg)+(J_))
	Mint            lower;
        Mint            IMSLFALSE = 0, IMSLTRUE = 1;
	Mint             _l0, i, ia, id, ifinc, iflag, ii, il, ip, ipp,
	                ir, ira, irb, is, iterc, itref, iu, iw, iwa, iwd,
	                iwr, iws, iww, iwwn, iwx, iwy, iwz, ix, iy, iz,
	                iza, j, jfinc, jflag, jl, k, k1, kdrop, kfinc, kflag,
	                kk, knext, lflag, mflag, nflag, nm, nu;
	Mfloat          big, cvmax, diagr, fdiff, fdiffa,
	                ga, gb, parinc, parnew, ratio, res, small,
	                step, sum, suma, sumb, sumc, sumx, sumy, temp,
	                tempa, vfact, vsmall, xmag, xmagr, zero;

	/* INITIAL ADDRESSES */
	vsmall = imsl_amach(4);
	small = imsl_amach(1);
	big = imsl_amach(2);
	if (small * big < F_ONE)
		small = F_ONE / big;

	zero = F_ZERO;
	iwz = *nmax;
	iwr = iwz + *nmax ** nmax;
	iww = iwr + (*nmax * (*nmax + 3)) / 2;
	iwd = iww + *nmax;
	iwx = iwd + *nmax;
	iwa = iwx + *nmax;
	/* SET SOME CONSTANTS. */
	vfact = F_ONE;
	/*
	 * SET SOME PARAMETERS. NUMBER LESS THAN VSMALL ARE ASSUMED TO BE
	 * NEGLIGIBLE. THE MULTIPLE OF I THAT IS ADDED TO G IS AT MOST DIAGR
	 * TIMES THE LEAST MULTIPLE OF I THAT GIVES POSITIVE DEFINITENESS. X
	 * IS RE-INITIALISED IF ITS MAGNITUDE IS REDUCED BY THE FACTOR XMAGR.
	 * A CHECK IS MADE FOR AN INCREASE IN F EVERY IFINC ITERATIONS, AFTER
	 * KFINC ITERATIONS ARE COMPLETED.
	 */
	diagr = F_TWO;
	xmagr = 1.0e-2;
	ifinc = 3;
	kfinc = imsl_i_max(10, *n);
	/*
	 * FIND THE RECIPROCALS OF THE LENGTHS OF THE CONSTRAINT NORMALS.
	 * RETURN IF A CONSTRAINT IS INFEASIBLE DUE TO A F_ZERO NORMAL.
	 */
	*nact = 0;
	for (k = 1; k <= *m; k++) {
		sum = imsl_sdot(*n, A(0, k - 1), *lda, A(0, k - 1), *lda);
		if (sum > F_ZERO)
			goto L_10;
		if (b[k - 1] == F_ZERO)
			goto L_20;
		*info = -k;
		if (k <= *meq)
			goto L_1020;
		if (b[k - 1] > 0)
			goto L_1020;
		goto L_20;
L_10:
		sum = F_ONE / sqrt(sum);
L_20:
		ia = iwa + k;
		w[ia - 1] = sum;
	}
	sset(*n, F_ONE, &w[iwa + *m], 1);
	/*
	 * IF NECESSARY INCREASE THE DIAGONAL ELEMENTS OF G.
	 */
	if (!*lql)
		goto L_150;
	*diag = F_ZERO;
	for (i = 1; i <= *n; i++) {
		id = iwd + i;
		w[id - 1] = *G(i - 1, i - 1);
		*diag = imsl_f_max(*diag, vsmall - w[id - 1]);
		if (i == *n)
			goto L_50;
		ii = i + 1;
		for (j = ii; j <= *n; j++) {
			ga = -imsl_f_min(w[id - 1], *G(j - 1, j - 1));
			gb = fabs(w[id - 1] - *G(j - 1, j - 1)) + fabs(*G(j - 1, i - 1));
			if (gb > small)
				ga += *G(j - 1, i - 1) ** G(j - 1, i - 1) / gb;
			*diag = imsl_f_max(*diag, ga);
		}
L_50:
		;
	}
	if (*diag <= F_ZERO)
		goto L_80;
L_60:
	*diag *= diagr;
	for (i = 1; i <= *n; i++) {
		id = iwd + i;
		*G(i - 1, i - 1) = *diag + w[id - 1];
	}
	/*
	 * FORM THE CHOLESKY FACTORISATION OF G. THE TRANSPOSE OF THE FACTOR
	 * WILL BE PLACED IN THE R-PARTITION OF W.
	 */
L_80:
	ir = iwr;
	for (j = 1; j <= *n; j++) {
		ira = iwr;
		irb = ir + 1;
		for (i = 1; i <= j; i++) {
			temp = *G(j - 1, i - 1);
			if (i != 1) {
				temp -= imsl_sdot(ir - irb + 1, &w[irb - 1], 1, &w[ira],
						  1);
				ira += ir - irb + 1;
			}
	
			ir += 1;
			ira += 1;
			if (i < j)
				w[ir - 1] = temp / w[ira - 1];
		}
		if (temp < vsmall)
			goto L_120;
		w[ir - 1] = sqrt(temp);
	}
	goto L_170;
	/*
	 * INCREASE FURTHER THE DIAGONAL ELEMENT OF G.
	 */
L_120:
	w[j - 1] = F_ONE;
	sumx = F_ONE;
	k = j;
L_130:
	sum = F_ZERO;
	ira = ir - 1;
	for (i = k; i <= j; i++) {
		sum += -w[ira - 1] * w[i - 1];
		ira += i;
	}
	ir -= k;
	k -= 1;
	w[k - 1] = sum / w[ir - 1];
	sumx += w[k - 1] * w[k - 1];
	if (k >= 2)
		goto L_130;
	*diag += vsmall - temp / sumx;
	goto L_60;
	/*
	 * STORE THE CHOLESKY FACTORISATION IN THE R-PARTITION OF W.
	 */
L_150:
	ir = iwr;
	for (i = 1; i <= *n; i++) {
		scopy(i, G(i - 1, 0), 1, &w[ir + i * (i - 1) / 2], 1);
	}
	/*
	 * SET Z THE INVERSE OF THE MATRIX IN R.
	 */
L_170:
	nm = *n - 1;
	for (i = 1; i <= *n; i++) {
		iz = iwz + i;
		sset(i - 1, F_ZERO, &w[iz - 1], *n);
		iz += *n * (i - 1);
		ir = iwr + (i + i * i) / 2;
		w[iz - 1] = F_ONE / w[ir - 1];
		if (i == *n)
			goto L_190;
		iza = iz;
		for (j = i; j <= nm; j++) {
			ir += i;
			sum = imsl_sdot((iz - iza) / *n + 1, &w[iza - 1], *n, &w[ir - 1],
					1);
			ir += (iz - iza) / *n + 1;
			iz += *n;
			w[iz - 1] = -sum / w[ir - 1];
		}
L_190:
		;
	}
	/*
	 * SET THE INITIAL VALUES OF SOME VARIABLES. ITERC COUNTS THE NUMBER
	 * OF ITERATIONS. ITREF IS SET TO F_ONE WHEN ITERATIVE REFINEMENT IS
	 * REQUIRED. JFINC INDICATES WHEN TO TEST FOR AN INCREASE IN F.
	 */
	iterc = 1;
	itref = 0;
	jfinc = -kfinc;
	/*
	 * SET X TO F_ZERO AND SET THE CORRESPONDING RESIDUALS OF THE
	 * KUHN-TUCKER CONDITIONS.
	 */
L_200:
	iflag = 1;
	iws = iww - *n;
	sset(*n, F_ZERO, x, 1);
	for (i = 1; i <= *n; i++) {
		iw = iww + i;
		w[iw - 1] = grad[i - 1];
		if (i > *nact)
			goto L_230;
		w[i - 1] = F_ZERO;
		is = iws + i;
		k = iact[i - 1];
		if (k <= *m)
			goto L_220;
		if (k > *mn)
			goto L_210;
		k1 = k - *m;
		w[is - 1] = xl[k1 - 1];
		goto L_230;
L_210:
		k1 = k - *mn;
		w[is - 1] = -xu[k1 - 1];
		goto L_230;
L_220:
		w[is - 1] = b[k - 1];
L_230:
		;
	}
	xmag = F_ZERO;
	vfact = F_ONE;
	if (*nact > 0)
		goto L_340;
	goto L_390;
	/*
	 * SET THE RESIDUALS OF THE KUHN-TUCKER CONDITIONS FOR GENERAL X.
	 */
L_240:
	iflag = 2;
	iws = iww - *n;
	for (i = 1; i <= *n; i++) {
		iw = iww + i;
		w[iw - 1] = grad[i - 1];
		if (*lql)
			goto L_270;
		id = iwd + i;
		w[id - 1] = F_ZERO;
		for (j = i; j <= *n; j++) {
			w[id - 1] += *G(j - 1, i - 1) * x[j - 1];
		}
		for (j = 1; j <= i; j++) {
			id = iwd + j;
			w[iw - 1] += *G(i - 1, j - 1) * w[id - 1];
		}
		goto L_290;
L_270:
		for (j = 1; j <= *n; j++) {
			w[iw - 1] += *G(j - 1, i - 1) * x[j - 1];
		}
L_290:
		;
	}
	if (*nact == 0)
		goto L_390;
	for (k = 1; k <= *nact; k++) {
		kk = iact[k - 1];
		is = iws + k;
		if (kk > *m)
			goto L_310;
		w[is - 1] = b[kk - 1];
		for (i = 1; i <= *n; i++) {
			iw = iww + i;
			w[iw - 1] += -w[k - 1] ** A(i - 1, kk - 1);
			w[is - 1] += -x[i - 1] ** A(i - 1, kk - 1);
		}
		goto L_330;
L_310:
		if (kk > *mn)
			goto L_320;
		k1 = kk - *m;
		iw = iww + k1;
		w[iw - 1] -= w[k - 1];
		w[is - 1] = xl[k1 - 1] - x[k1 - 1];
		goto L_330;
L_320:
		k1 = kk - *mn;
		iw = iww + k1;
		w[iw - 1] += w[k - 1];
		w[is - 1] = -xu[k1 - 1] + x[k1 - 1];
L_330:
		;
	}
	/*
	 * PRE-MULTIPLY THE VECTOR IN THE S-PARTITION OF W BY THE INVERS OF R
	 * TRANSPOSE.
	 */
L_340:
	ir = iwr;
	ip = iww + 1;
	ipp = iww + *n;
	il = iws + 1;
	iu = iws + *nact;
	for (i = il; i <= iu; i++) {
		sum = imsl_sdot(i - il, &w[ir], 1, &w[il - 1], 1);
		ir += i - il + 1;
		w[i - 1] = (w[i - 1] - sum) / w[ir - 1];
	}
	/*
	 * SHIFT X TO SATISFY THE ACTIVE CONSTRAINTS AND MAKE THE
	 * CORRESPONDING CHANGE TO THE GRADIENT RESIDUALS.
	 */
	for (i = 1; i <= *n; i++) {
		iz = iwz + i;
		sum = imsl_sdot(iu - il + 1, &w[il - 1], 1, &w[iz - 1], *n);
		iz += (iu - il + 1) ** n;
		x[i - 1] += sum;
		if (!*lql) {
			id = iwd + i;
                        _l0 = *n - i + 1;
			w[id - 1] = sum * l_ssum(&_l0, G(i - 1, i - 1), ldg);
			iw = iww + i;
			for (j = 1; j <= i; j++) {
				id = iwd + j;
				w[iw - 1] += *G(i - 1, j - 1) * w[id - 1];
			}
		} else {
			for (j = 1; j <= *n; j++) {
				iw = iww + j;
				w[iw - 1] += sum ** G(j - 1, i - 1);
			}
		}
	}
	/*
	 * FORM THE SCALAR PRODUCT OF THE CURRENT GRADIENT RESIDUALS WITH
	 * EACH COLUMN OF Z.
	 */
L_390:
	kflag = 1;
	goto L_1260;
L_400:
	if (*nact != *n) {
		/*
		 * SHIFT X SO THAT IT SATISFIES THE REMAINING KUHN-TUCKER
		 * CONDITIONS.
		 */
		il = iws + *nact + 1;
		iza = iwz + *nact ** n;
		for (i = 1; i <= *n; i++) {
			iz = iza + i;
			sum = imsl_sdot(iww - il + 1, &w[iz - 1], *n, &w[il - 1],
					1);
			iz += (iww - il + 1) ** n;
			x[i - 1] -= sum;
		}
		*info = iterc;
		if (*nact == 0)
			goto L_440;
	}
	/* UPDATE THE LAGRANGE MULTIPLIERS. */
	lflag = 3;
	goto L_1030;
L_420:
	for (k = 1; k <= *nact; k++) {
		iw = iww + k;
		w[k - 1] += w[iw - 1];
	}
	/*
	 * REVISE THE VALUES OF XMAG. BRANCH IF ITERATIVE REFINEMENT IS
	 * REQUIRED.
	 */
L_440:
	jflag = 1;
	goto L_1230;
L_450:
	if (iflag == itref)
		goto L_240;
	/*
	 * DELETE A CONSTRAINT IF A LAGRANGE MULTIPLIER OF AN INEQUALITY
	 * CONSTRAINT IS NEGATIVE.
	 */
	kdrop = 0;
	goto L_470;
L_460:
	kdrop += 1;
	if (w[kdrop - 1] >= F_ZERO)
		goto L_470;
	if (iact[kdrop - 1] <= *meq)
		goto L_470;
	nu = *nact;
	mflag = 1;
	goto L_1120;
L_470:
	if (kdrop < *nact)
		goto L_460;
	/*
	 * SEEK THE GREATEAST NORMALISED CONSTRAINT VIOLATION, DISREGARDING
	 * ANY THAT MAY BE DUE TO COMPUTER ROUNDING ERRORS.
	 */
L_480:
	cvmax = F_ZERO;

	for (k = 1; k <= *m; k++) {
		ia = iwa + k;
		if (w[ia - 1] > F_ZERO) {
			sum = imsl_sdot(*n, x, 1, A(0, k - 1), *lda) - b[k - 1];
			sumx = -sum * w[ia - 1];
			if (k <= *meq)
				sumx = fabs(sumx);
			if (sumx > cvmax) {
                                _l0 = 1;
				temp = fabs(b[k - 1]) + l_a1ot(*n, x, _l0,
					    A(0, k - 1), *lda);
				tempa = temp + fabs(sum);
				if (tempa > temp) {
					temp += 1.5e0 * fabs(sum);
					if (temp > tempa) {
						cvmax = sumx;
						res = sum;
						knext = k;
					}
				}
			}
		}
	}

	for (k = 1; k <= *n; k++) {
		lower = IMSLTRUE;
		ia = iwa + *m + k;
		if (w[ia - 1] <= F_ZERO)
			goto L_520;
		sum = xl[k - 1] - x[k - 1];
		if (sum < 0)
			goto L_500;
		if (sum > 0)
			goto L_510;
		goto L_520;
L_500:
		sum = x[k - 1] - xu[k - 1];
		lower = IMSLFALSE;
L_510:
		if (sum <= cvmax)
			goto L_520;
		cvmax = sum;
		res = -sum;
		knext = k + *m;
		if (lower)
			goto L_520;
		knext = k + *mn;
L_520:
		;
	}
	/* TEST FOR CONVERGENCE */
	*info = iterc;
	if (cvmax <= vsmall)
		goto L_990;
	/*
	 * RETURN IF, DUE TO ROUNDING ERRORS, THE ACTUAL CHANGE IN X MAY NOT
	 * INCREASE THE OBJECTIVE FUNCTION
	 */
	jfinc += 1;
	if (jfinc == 0)
		goto L_590;
	if (jfinc != ifinc)
		goto L_610;
	fdiff = F_ZERO;
	fdiffa = F_ZERO;
	for (i = 1; i <= *n; i++) {
		sum = F_TWO * grad[i - 1];
		sumx = fabs(sum);
		if (*lql)
			goto L_550;
		id = iwd + i;
		w[id - 1] = F_ZERO;
		for (j = i; j <= *n; j++) {
			ix = iwx + j;
			w[id - 1] += *G(j - 1, i - 1) * (w[ix - 1] + x[j - 1]);
		}

		for (j = 1; j <= i; j++) {
			id = iwd + j;
			temp = *G(i - 1, j - 1) * w[id - 1];
			sum += temp;
			sumx += fabs(temp);
		}
		goto L_570;
L_550:
		for (j = 1; j <= *n; j++) {
			ix = iwx + j;
			temp = *G(j - 1, i - 1) * (w[ix - 1] + x[j - 1]);
			sum += temp;
			sumx += fabs(temp);
		}
L_570:
		ix = iwx + i;
		fdiff += sum * (x[i - 1] - w[ix - 1]);
		fdiffa += sumx * fabs(x[i - 1] - w[ix - 1]);
	}
	*info = 0;
	sum = fdiffa + fdiff;
	if (sum <= fdiffa)
		goto L_990;
	temp = fdiffa + 1.5e0 * fdiff;
	if (temp <= sum)
		goto L_990;
	jfinc = 0;
L_590:
	for (i = 1; i <= *n; i++) {
		ix = iwx + i;
		w[ix - 1] = x[i - 1];
	}
	/*
	 * FORM THE SCALAR PRODUCT OF THE NEW CONSTRAINT NORMAL WITH EACH
	 * COLUMN OF Z. PARNEW WILL BECOME THE LAGRANGE MULTIPLIER OF THE NEW
	 * CONSTRAINT.
	 */
L_610:
	iterc += 1;
	iws = iwr + (*nact + *nact ** nact) / 2;
	if (knext > *m)
		goto L_630;
	for (i = 1; i <= *n; i++) {
		iw = iww + i;
		w[iw - 1] = *A(i - 1, knext - 1);
	}
	goto L_680;
L_630:
	for (i = 1; i <= *n; i++) {
		iw = iww + i;
		w[iw - 1] = F_ZERO;
	}
	k1 = knext - *m;
	if (k1 > *n)
		goto L_660;
	iw = iww + k1;
	w[iw - 1] = F_ONE;
	iz = iwz + k1;
	for (i = 1; i <= *n; i++) {
		is = iws + i;
		w[is - 1] = w[iz - 1];
		iz += *n;
	}
	goto L_690;
L_660:
	k1 = knext - *mn;
	iw = iww + k1;
	w[iw - 1] = -F_ONE;
	iz = iwz + k1;
	for (i = 1; i <= *n; i++) {
		is = iws + i;
		w[is - 1] = -w[iz - 1];
		iz += *n;
	}
	goto L_690;
L_680:
	kflag = 2;
	goto L_1260;
L_690:
	parnew = F_ZERO;
	/*
	 * APPLY GIVENS ROTATIONS TO MAKE THE LAST (N-NACT-2) SCALAR PRODUCTS
	 * EQUAL TO F_ZERO.
	 */
	if (*nact == *n)
		goto L_740;
	nu = *n;
	nflag = 1;
	goto L_1180;
	/*
	 * BRANCH IF THERE IS NO NEED TO DELETE A CONSTRAINT.
	 */
L_700:
	is = iws + *nact;
	if (*nact == 0)
		goto L_930;
	suma = F_ZERO;
	sumb = F_ZERO;
	iz = iwz + *nact ** n;
	sumc = imsl_sdot(*n, &w[iz], 1, &w[iz], 1);
	for (i = 1; i <= *n; i++) {
		iz += 1;
		iw = iww + i;
		suma += w[iw - 1] * w[iz - 1];
		sumb += fabs(w[iw - 1] * w[iz - 1]);
	}
	temp = sumb + .1e0 * fabs(suma);
	tempa = sumb + .2e0 * fabs(suma);
	if (temp <= sumb)
		goto L_740;
	if (tempa <= temp)
		goto L_740;
	if (sumb > vsmall)
		goto L_720;
	goto L_740;
L_720:
	sumc = sqrt(sumc);
	ia = iwa + knext;
	if (knext <= *m)
		sumc /= w[ia - 1];
	temp = sumc + .1e0 * fabs(suma);
	tempa = sumc + .2e0 * fabs(suma);
	if (temp <= sumc)
		goto L_730;
	if (tempa <= temp)
		goto L_730;
	goto L_930;
	/*
	 * CALCULATE THE MULTIPLIERS FOR THE NEW CONSTRAINT NORMAL EXPRESSED
	 * IN TERMS OF THE ACTIVE CONSTRAINT NORMALS. THEN WORK OUT WHICH
	 * CONTRAINT TO DROP.
	 */
L_730:
	lflag = 4;
	goto L_1030;
L_740:
	lflag = 1;
	goto L_1030;
	/*
	 * COMPLETE THE TEST FOR LINEARLY DEPENDENT CONSTRAINTS.
	 */
L_750:
	if (knext > *m)
		goto L_790;
	for (i = 1; i <= *n; i++) {
		suma = *A(i - 1, knext - 1);
		sumb = fabs(suma);
		if (*nact == 0)
			goto L_770;
		for (k = 1; k <= *nact; k++) {
			kk = iact[k - 1];
			if (kk <= *m)
				goto L_758;
			kk -= *m;
			temp = zero;
			if (kk == i)
				temp = w[iww + kk - 1];
			kk -= *n;
			if (kk == i)
				temp = -w[iww + kk - 1];
			goto L_759;
	L_758:
			;
			iw = iww + k;
			temp = w[iw - 1] ** A(i - 1, kk - 1);
	L_759:
			;
			suma -= temp;
			sumb += fabs(temp);
		}
L_770:
		if (suma <= vsmall)
			goto L_780;
		temp = sumb + .1e0 * fabs(suma);
		tempa = sumb + .2e0 * fabs(suma);
		if (temp <= sumb)
			goto L_780;
		if (tempa <= temp)
			goto L_780;
		goto L_920;
L_780:
		;
	}
	lflag = 1;
	goto L_1080;
L_790:
	k1 = knext - *m;
	if (k1 > *n)
		k1 -= *n;
	for (i = 1; i <= *n; i++) {
		suma = F_ZERO;
		if (i != k1)
			goto L_800;
		suma = F_ONE;
		if (knext > *mn)
			suma = -F_ONE;
L_800:
		sumb = fabs(suma);
		if (*nact == 0)
			goto L_840;
		for (k = 1; k <= *nact; k++) {
			kk = iact[k - 1];
			if (kk <= *m)
				goto L_810;
			kk -= *m;
			temp = F_ZERO;
			if (kk == i)
				temp = w[iww + kk - 1];
			kk -= *n;
			if (kk == i)
				temp = -w[iww + kk - 1];
			goto L_820;
	L_810:
			iw = iww + k;
			temp = w[iw - 1] ** A(i - 1, kk - 1);
	L_820:
			suma -= temp;
			sumb += fabs(temp);
		}
L_840:
		temp = sumb + .1e0 * fabs(suma);
		tempa = sumb + .2e0 * fabs(suma);
		if (temp <= sumb)
			goto L_850;
		if (tempa <= temp)
			goto L_850;
		goto L_920;
L_850:
		;
	}
	lflag = 1;
	goto L_1080;
	/*
	 * BRANCH IF THE CONTRAINTS ARE INCONSISTENT.
	 */
L_860:
	*info = -knext;
	if (kdrop == 0)
		goto L_990;
	parinc = ratio;
	parnew = parinc;
	/*
	 * REVISE THE LAGRANGE MULTIPLIERS OF THE ACTIVE CONSTRAINTS.
	 */
L_870:
	if (*nact == 0)
		goto L_890;
	for (k = 1; k <= *nact; k++) {
		iw = iww + k;
		w[k - 1] += -parinc * w[iw - 1];
		if (iact[k - 1] > *meq)
			w[k - 1] = imsl_f_max(F_ZERO, w[k - 1]);
	}
L_890:
	if (kdrop == 0)
		goto L_970;
	/*
	 * DELETE THE CONSTRAINT TO BE DROPPED. SHIFT THE VECTOR OF SCALAR
	 * PRODUCTS. THEN, IF APPROPRIATE, MAKE ONE MORE SCALAR PRODUCT
	 * F_ZERO.
	 */
	nu = *nact + 1;
	mflag = 2;
	goto L_1120;
L_900:
	iws += -*nact - 1;
	nu = imsl_i_min(*n, nu);
	for (i = 1; i <= nu; i++) {
		is = iws + i;
		j = is + *nact;
		w[is - 1] = w[j];
	}
	nflag = 2;
	goto L_1180;
	/*
	 * CALCULATE THE STEP TO THE VIOLATED CONSTRAINT.
	 */
L_920:
	is = iws + *nact;
L_930:
	sumy = w[is];
	step = -res / sumy;
	parinc = step / sumy;
	if (*nact == 0)
		goto L_950;
	/*
	 * CALCULATE THE CHANGES TO THE LAGRANGE MULTIPLIERS, AND REDUCE THE
	 * STEP ALONG THE NEW SEARCH DIRECTION IF NECESSARY.
	 */
	lflag = 2;
	goto L_1030;
L_940:
	if (kdrop == 0)
		goto L_950;
	temp = F_ONE - ratio / parinc;
	if (temp <= F_ZERO)
		kdrop = 0;
	if (kdrop == 0)
		goto L_950;
	step = ratio * sumy;
	parinc = ratio;
	res *= temp;
	/*
	 * UPDATE X AND THE LAGRANGE MULTIPIERS. DROP A CONSTRAINT IF THE
	 * FULL STEP IS NOT TAKEN.
	 */
L_950:
	iwy = iwz + *nact ** n;
	for (i = 1; i <= *n; i++) {
		iy = iwy + i;
		x[i - 1] += step * w[iy - 1];
	}
	parnew += parinc;
	if (*nact >= 1)
		goto L_870;
	/*
	 * ADD THE NEW CONSTRAINT TO THE ACTIVE SET.
	 */
L_970:
	*nact += 1;
	w[*nact - 1] = parnew;
	iact[*nact - 1] = knext;
	ia = iwa + knext;
	if (knext > *mn)
		ia -= *n;
	w[ia - 1] = -w[ia - 1];
	/*
	 * ESTIMATE THE MAGNITUDE OF X. THEN BEGIN A NEW ITERATION,
	 * RE-INITILISING X IF THIS MAGNITUDE IS SMALL.
	 */
	jflag = 2;
	goto L_1230;
L_980:
	if (sum < (xmagr * xmag))
		goto L_200;
	if (itref > 0)
		goto L_240;
	goto L_480;
	/*
	 * INITIATE ITERATIVE REFINEMENT IF IT HAS NOT YET BEEN USED, OR
	 * RETURN AFTER RESTORING THE DIAGONAL ELEMENTS OF G.
	 */
L_990:
	if (iterc == 0)
		goto L_1000;
	itref += 1;
	jfinc = -1;
	if (itref == 1)
		goto L_240;
L_1000:
	if (!*lql)
		return;
	for (i = 1; i <= *n; i++) {
		id = iwd + i;
		*G(i - 1, i - 1) = w[id - 1];
	}
L_1020:
	return;
	/*
	 * THE REMAINIG INSTRUCTIONS ARE USED AS SUBROUTINES. CALCULATE THE
	 * LAGRANGE MULTIPLIERS BY PRE-MULTIPLYING THE VECTOR IN THE
	 * S-PARTITION OF W BY THE INVERSE OF R.
	 */
L_1030:
	ir = iwr + (*nact + *nact ** nact) / 2;
	i = *nact;
	sum = F_ZERO;
	goto L_1070;
L_1040:
	ira = ir - 1;
	sum = F_ZERO;
	if (*nact == 0)
		goto L_1060;
	for (j = i; j <= *nact; j++) {
		iw = iww + j;
		sum += w[ira - 1] * w[iw - 1];
		ira += j;
	}
L_1060:
	ir -= i;
	i -= 1;
L_1070:
	iw = iww + i;
	is = iws + i;
	w[iw - 1] = (w[is - 1] - sum) / w[ir - 1];
	if (i > 1)
		goto L_1040;
	if (lflag == 3)
		goto L_420;
	if (lflag == 4)
		goto L_750;
	/*
	 * CALCULATE THE NEXT CONSTRAINT TO DROP.
	 */
L_1080:
	ip = iww + 1;
	ipp = iww + *nact;
	kdrop = 0;
	if (*nact == 0)
		goto L_1110;
	for (k = 1; k <= *nact; k++) {
		if (iact[k - 1] <= *meq)
			goto L_1100;
		iw = iww + k;
		if ((res * w[iw - 1]) >= F_ZERO)
			goto L_1100;
		temp = w[k - 1] / w[iw - 1];
		if (kdrop == 0)
			goto L_1090;
		if (fabs(temp) >= fabs(ratio))
			goto L_1100;
L_1090:
		kdrop = k;
		ratio = temp;
L_1100:
		;
	}
L_1110:
	if (lflag == 1)
		goto L_860;
	if (lflag == 2)
		goto L_940;
	/*
	 * DROP THE CONSTRAINT IN POSITION KDROP IN THE ACTIVE SET.
	 */
L_1120:
	ia = iwa + iact[kdrop - 1];
	if (iact[kdrop - 1] > *mn)
		ia -= *n;
	w[ia - 1] = -w[ia - 1];
	if (kdrop == *nact)
		goto L_1170;
	/*
	 * SET SOME INDICES AND CALCULATE THE ELEMENTS OF THE NEXT GIVENS
	 * ROTATION.
	 */
	iz = iwz + kdrop ** n;
	ir = iwr + (kdrop + kdrop * kdrop) / 2;
L_1130:
	ira = ir;
	ir += kdrop + 1;
	temp = imsl_f_max(fabs(w[ir - 2]), fabs(w[ir - 1]));
	sum = temp * sqrt(imsl_fi_power(w[ir - 2] / temp, 2) + imsl_fi_power(w[ir - 1] / temp, 2));
	ga = w[ir - 2] / sum;
	gb = w[ir - 1] / sum;
	/* EXCHANGE THE COLUMNS OF R. */
	for (i = 1; i <= kdrop; i++) {
		ira += 1;
		j = ira - kdrop;
		temp = w[ira - 1];
		w[ira - 1] = w[j - 1];
		w[j - 1] = temp;
	}
	w[ir - 1] = F_ZERO;
	/* APPLY THE ROTATION TO THE ROWS OF R. */
	w[j - 1] = sum;
	kdrop += 1;
	for (i = kdrop; i <= nu; i++) {
		temp = ga * w[ira - 1] + gb * w[ira];
		w[ira] = ga * w[ira] - gb * w[ira - 1];
		w[ira - 1] = temp;
		ira += i;
	}
	/*
	 * APPLY THE ROTATION TO THE COLUMNS OF Z.
	 */
	for (i = 1; i <= *n; i++) {
		iz += 1;
		j = iz - *n;
		temp = ga * w[j - 1] + gb * w[iz - 1];
		w[iz - 1] = ga * w[iz - 1] - gb * w[j - 1];
		w[j - 1] = temp;
	}
	/*
	 * REVISE IACT AND THE LAGRANGE MULTIPLIERS.
	 */
	iact[kdrop - 2] = iact[kdrop - 1];
	w[kdrop - 2] = w[kdrop - 1];
	if (kdrop < *nact)
		goto L_1130;
L_1170:
	*nact -= 1;
	if (mflag == 1)
		goto L_240;
	if (mflag == 2)
		goto L_900;
	/*
	 * APPLY GIVENS ROTATION TO REDUCE SOME OF THE SCALAR PRODUCTS IN THE
	 * S-PARTITION OF W TO F_ZERO.
	 */
L_1180:
	iz = iwz + nu ** n;
L_1190:
	iz -= *n;
L_1200:
	is = iws + nu;
	nu -= 1;
	if (nu == *nact)
		goto L_1220;
	if (w[is - 1] == F_ZERO)
		goto L_1190;
	temp = imsl_f_max(fabs(w[is - 2]), fabs(w[is - 1]));
	sum = temp * sqrt(imsl_fi_power(w[is - 2] / temp, 2) + imsl_fi_power(w[is - 1] / temp, 2));
	ga = w[is - 2] / sum;
	gb = w[is - 1] / sum;
	w[is - 2] = sum;
	for (i = 1; i <= *n; i++) {
		k = iz + *n;
		temp = ga * w[iz - 1] + gb * w[k - 1];
		w[k - 1] = ga * w[k - 1] - gb * w[iz - 1];
		w[iz - 1] = temp;
		iz -= 1;
	}
	goto L_1200;
L_1220:
	if (nflag == 1)
		goto L_700;
	if (nflag == 2)
		goto L_920;
	/*
	 * CALCULATE THE MAGNITUDE OF X AN REVISE XMAG.
	 */
L_1230:
	sum = F_ZERO;
	for (i = 1; i <= *n; i++) {
		sum += fabs(x[i - 1]) * vfact * (fabs(grad[i - 1]) + fabs(*G(i - 1, i - 1) *
								 x[i - 1]));
		if (*lql)
			goto L_1240;
		if (sum < 1.0e-30)
			goto L_1240;
		vfact *= 1.0e-10;
		sum *= 1.0e-10;
		xmag *= 1.0e-10;
L_1240:
		;
	}

	xmag = imsl_f_max(xmag, sum);
	if (jflag == 1)
		goto L_450;
	if (jflag == 2)
		goto L_980;
	/*
	 * PRE-MULTIPLY THE VECTOR IN THE W-PARTITION OF W BY Z TRANSPOSE.
	 */
L_1260:
	jl = iww + 1;
	iz = iwz;
	for (i = 1; i <= *n; i++) {
		is = iws + i;
		w[is - 1] = F_ZERO;
		iwwn = iww + *n;
		w[is - 1] = imsl_sdot(iwwn - jl + 1, &w[jl - 1], 1, &w[iz], 1);
		iz += iwwn - jl + 1;
	}
	if (kflag == 1)
		goto L_400;
	if (kflag == 2)
		goto L_690;
	return;
}				/* end of function */
#undef A
#undef G
/*
  -----------------------------------------------------------------------
    IMSL Name:  A1OT/DA1OT (Single/Double precision version)
 
    Computer:   FORC/SINGLE
 
    Revised:    August 9, 1986
 
    Purpose:    Compute sum of absolute values of products.
 
    Usage:      A1OT(N, SX, INCX, SY, INCY)
 
    Arguments:
       N      - Length of vectors X and Y.  (Input)
       SX     - Real vector of length MAX(N*IABS(INCX),1).  (Input)
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be.. SX(1+(I-1)*INCX) if INCX .GE. 0
                or SX(1+(I-N)*INCX) if INCX .LT. 0.
       SY     - Real vector of length MAX(N*IABS(INCY),1).  (Input)
       INCY   - Displacement between elements of SY.  (Input)
                Y(I) is defined to be.. SY(1+(I-1)*INCY) if INCY .GE. 0
                or SY(1+(I-N)*INCY) if INCY .LT. 0.
       A1OT   - Sum from I=1 to N of ABS(X(I)*Y(I)).  (Output)
                X(I) and Y(I) refer to specific elements of SX and SY,
                respectively.  See INCX and INCY argument descriptions.
 
    GAMS:       D1a4
 
    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mfloat l_a1ot(Mint n, Mfloat sx[], Mint incx, Mfloat sy[], Mint incy)
#else
static Mfloat l_a1ot(n, sx, incx, sy, incy)
    Mint            n;
    Mfloat          sx[];
    Mint            incx;
    Mfloat          sy[];
    Mint            incy;
#endif
{
        Mint            i, ix, iy;
        Mfloat          a1ot_v;


        a1ot_v = F_ZERO;
        if (n > 0) {
                if (incx != 1 || incy != 1) {
                        /* CODE FOR UNEQUAL INCREMENTS */
                        ix = 1;
                        iy = 1;
                        if (incx < 0)
                                ix = (-n + 1) * incx + 1;
                        if (incy < 0)
                                iy = (-n + 1) * incy + 1;
                        for (i = 1; i <= n; i++) {
                                a1ot_v += fabs(sx[ix - 1] * sy[iy - 1]);
                                ix += incx;
                                iy += incy;
                        }
                } else {
                        for (i = 1; i <= n; i++) {
                                a1ot_v += fabs(sx[i - 1] * sy[i - 1]);
                        }
                }
        }
        return (a1ot_v);
}                               /* end of function */
/* -----------------------------------------------------------------------
    IMSL Name:  SHPROD (Single precision version)
 
    Computer:   FORC/SINGLE
 
    Revised:    August 9, 1986
 
    Purpose:    Compute the Hadamard product of two single precision
                vectors.
 
    Usage:      CALL SHPROD (N, SX, INCX, SY, INCY, SZ, INCZ)
 
    Arguments:
       N      - Length of vectors X, Y and Z.  (Input)
       SX     - Real vector of length MAX(N*IABS(INCX),1).  (Input)
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be.. SX(1+(I-1)*INCX) if INCX .GE. 0
                or SX(1+(I-N)*INCX) if INCX .LT. 0.
       SY     - Real vector of length MAX(N*IABS(INCY),1).  (Input)
       INCY   - Displacement between elements of SY.  (Input)
                Y(I) is defined to be.. SY(1+(I-1)*INCY) if INCY .GE. 0
                or SY(1+(I-N)*INCY) if INCY .LT. 0.
       SZ     - Real vector of length MAX(N*IABS(INCZ),1).  (Output)
                SZ returns the Hadamard product of SX and SY,
                     Z(I) = X(I)*Y(I) for I=1,...,N.
       INCZ   - Displacement between elements of SZ.  (Input)
                Z(I) is defined to be.. SZ(1+(I-1)*INCZ) if INCZ .GE. 0
                or SZ(1+(I-N)*INCZ) if INCZ .LT. 0.
 
    Keyword:    Level 1 BLAS
 
    GAMS:       D1a
 
    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support
 
    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
 
    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_shprod(Mint *n, Mfloat sx[], Mint *incx, Mfloat sy[], Mint *incy,
                    Mfloat sz[], Mint *incz)
#else
static void l_shprod(n, sx, incx, sy, incy, sz, incz)
        Mint            *n;
        Mfloat           sx[];
        Mint            *incx;
        Mfloat           sy[];
        Mint            *incy;
        Mfloat           sz[];
        Mint            *incz;
#endif
{
        Mint             i, ix, iy, iz;


        if (*n <= 0)
                goto L_9000;
        if ((*incx != 1 || *incy != 1) || *incz != 1) {

                /*
                 * CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS NOT EQUAL
                 * TO 1
                 */
                ix = 1;
                iy = 1;
                iz = 1;
                if (*incx < 0)
                        ix = (-*n + 1) ** incx + 1;
                if (*incy < 0)
                        iy = (-*n + 1) ** incy + 1;
                if (*incz < 0)
                        iz = (-*n + 1) ** incz + 1;
                for (i = 1; i <= *n; i++) {
                        sz[iz - 1] = sx[ix - 1] * sy[iy - 1];
                        ix += *incx;
                        iy += *incy;
                        iz += *incz;
                }
        } else {
                for (i = 1; i <= *n; i++) {
                        sz[i - 1] = sx[i - 1] * sy[i - 1];
                }
        }
L_9000:
        return;
}                               /* end of function */

/*
  -----------------------------------------------------------------------
    IMSL Name:  ISMAX (Single precision version)
 
    Computer:   FORC/SINGLE
 
    Revised:    August 9, 1986
 
    Purpose:    Find the smallest index of the component of a
                single-precision vector having maximum value.
 
    Usage:      ISMAX(N, SX, INCX)
 
    Arguments:
       N      - Length of vector X.  (Input)
       SX     - Real vector of length N*INCX.  (Input)
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be SX(1+(I-1)*INCX). INCX must be
                greater than zero.
       ISMAX  - The smallest index I such that X(I)  is the maximum of
                X(J) for J=1 to N.  (Output)
                X(I) refers to a specific element of SX. See INCX
                argument description.
 
    Keyword:    Level 1 BLAS
 
    GAMS:       D1a2
 
    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support
 
    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
 
    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.
 
  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mint l_ismax(Mint n, Mfloat sx[], Mint incx)
#else
static Mint l_ismax(n, sx, incx)
        Mint           n;
        Mfloat           sx[];
        Mint           incx;
#endif
{
        Mint            i, ismax_v, ix;
        Mfloat           smax;


        ismax_v = 0;
        if (n >= 1) {
                ismax_v = 1;
                if (n != 1) {
                        if (incx != 1) {
                                /* CODE FOR INCREMENT NOT EQUAL TO 1 */
                                ix = 1;
                                smax = sx[0];
                                ix += incx;
                                for (i = 2; i <= n; i++) {
                                        if (sx[ix - 1] > smax) {
                                                ismax_v = i;
                                                smax = sx[ix - 1];
                                        }
                                        ix += incx;
                                }
                        } else {
                                /* CODE FOR INCREMENT EQUAL TO 1 */
                                smax = sx[0];
                                for (i = 2; i <= n; i++) {
                                        if (sx[i - 1] > smax) {
                                                ismax_v = i;
                                                smax = sx[i - 1];
                                        }
                                }
                        }
                }
        }
        return (ismax_v);
}                               /* end of function */

/*  -----------------------------------------------------------------------
    IMSL Name:  SSUM (Single precision version)
 
    Computer:   FORC/SINGLE
 
    Revised:    August 9, 1986
 
    Purpose:    Sum the values of a single precision vector.
 
    Usage:      SSUM(N, SX, INCX)
 
    Arguments:
       N      - Length of vectors X.  (Input)
       SX     - Real vector of length N*INCX.  (Input)
       INCX   - Displacement between elements of SX.  (Input)
                X(I) is defined to be SX(1+(I-1)*INCX). INCX must be
                greater than 0.
       SSUM   - Single precision sum from I=1 to N of X(I).  (Output)
                X(I) refers to a specific element of SX.
 
    Keyword:    Level 1 BLAS
 
    GAMS:       D1a
 
    Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
                STAT/LIBRARY Mathematical Support
 
    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
 
    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.
 
  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mfloat l_ssum(Mint *n, Mfloat sx[], Mint *incx)
#else
static Mfloat l_ssum(n, sx, incx)
        Mint            *n;
        Mfloat          *sx;
        Mint            *incx;
#endif
{
        Mint            _d_l, _d_m, _do0, _do1, i, nincx;
        Mfloat          ssum_v;


        ssum_v = F_ZERO;
        if (*n > 0) {
                if (*incx != 1) {
                        /* CODE FOR INCREMENT NOT EQUAL TO 1 */
                        nincx = *n * *incx;
                        for (i = 1, _do0 = DOCNT(1, nincx, _do1 = *incx); _do0 > 0; i += _do1, _do0--) {
                                ssum_v += sx[i - 1];
                        }
                } else {
                        for (i = 1; i <= *n; i++) {
                                ssum_v += sx[i - 1];
                        }
                }
        }
        return (ssum_v);
}  
