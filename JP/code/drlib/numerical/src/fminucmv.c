#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK   l_min_uncon_multivar(Mfloat (*fcn)(Mint, Mfloat[]),
                                      Mint n, va_list argptr);
static void      l_u2inf(Mfloat (*fcn) (Mint, Mfloat[]), Mint *n, 
                         Mfloat xguess[], Mfloat xscale[],
                         Mfloat *fscale, Mfloat *grad_tol,
                         Mfloat *step_tol, Mfloat *fcn_tol,
                         Mfloat *max_step, Mint *ndigit,
                         Mint *maxitn, Mint *maxfcn, Mint *maxgrad,
                         Mint *ihess, Mfloat x[], Mfloat *fvalue,
                         Mfloat wk[]); 
static void      l_u2ing(Mfloat (*fcn) (Mint, Mfloat[]),
                         void (*grad) (Mint, Mfloat[], Mfloat[]),
                         Mint *n, Mfloat xguess[], Mfloat xscale[],
                         Mfloat *fscale, Mfloat *grad_tol,
                         Mfloat *step_tol, Mfloat *fcn_tol,
                         Mfloat *max_step, Mint *ndigit,
                         Mint *maxitn, Mint *maxfcn, Mint *maxgrad,
                         Mint *ihess, Mfloat x[], Mfloat *fvalue,
                         Mfloat wk[]);
static void      l_u3inf(Mfloat (*fcn)(Mint, Mfloat[]), Mint *n,
                         Mfloat xc[], Mfloat xscale[],
                         Mfloat *fscale, Mint iparam[],
                         Mfloat rparam[], Mfloat *fvalue,
                         Mfloat xp[], Mfloat sc[], Mfloat snwtn[],
                         Mfloat gc[], Mfloat gp[], Mfloat *h,
                         Mint *ldh, Mfloat wk1[], Mfloat wk2[],
                         Mfloat wk3[]); 
static void      l_u3ing(Mfloat (*fcn)(Mint, Mfloat[]),
                         void (*grad) (Mint, Mfloat[], Mfloat[]),
                         Mint *n, Mfloat xc[], Mfloat xscale[],
                         Mfloat *fscale, Mint iparam[],
                         Mfloat rparam[], Mfloat *fvalue,
                         Mfloat xp[], Mfloat sc[], Mfloat snwtn[],
                         Mfloat gc[], Mfloat gp[], Mfloat *h,
                         Mint *ldh, Mfloat wk1[], Mfloat wk2[],
                         Mfloat wk3[]);
static void      l_u5inf(Mint *n, Mfloat xc[], Mfloat xscale[],
                         Mfloat *fscale, Mint *usrhes, Mint iparam[],
                         Mfloat rparam[]); 
static void      l_u6inf(Mint *n, Mfloat xp[], Mfloat sc[],
                         Mfloat *fp, Mfloat gp[], Mfloat xscale[],
                         Mfloat *fscale, Mint *icode, Mint *iter,
                         Mint *nfcn, Mint *ngrad, Mint *nhess,
                         Mint *usrhes, Mint *mxtake);  
static void      l_u8inf(Mint *n, Mfloat sc[], Mfloat gc[],
                         Mfloat gp[], Mfloat *epsfcn, Mint *usrder,
                         Mfloat *a, Mint *lda, Mfloat y[],
                         Mfloat t[], Mfloat u[]);  
static void      l_u9inf(Mint *n, Mfloat *fx0, Mfloat *fscale,
                         Mfloat xscale[], Mint *ihess, Mfloat *h,
                         Mint *ldh);
static void      l_u10nf(Mint *n, Mfloat *h, Mint *ldh, Mfloat gc[],
                         Mfloat snwtn[]); 
static void      l_u11nf(Mint *n, Mfloat y[], Mint *k, Mfloat z[],
                         Mfloat x[]); 
static void      l_u13nf(Mint *n, Mfloat *h, Mint *ldh, Mfloat gc[],
                         Mfloat snwtn[]); 
static void      l_u14nf(Mint *n, Mfloat *h, Mint *ldh, Mfloat y[],
                         Mfloat snwtn[]);  
static void      l_u15nf(Mint *method, Mint *n, Mfloat u[],
                         Mfloat v[], Mfloat *q, Mint *ldq,
                         Mfloat *a, Mint *lda);  
static void      l_u17nf(Mfloat (*fcn)(Mint, Mfloat[]), Mint *fdiff,
                         Mint *n, Mfloat xc[], Mfloat *fc,
                         Mfloat gc[], Mfloat sn[], Mfloat xscale[],
                         Mfloat *stepmx, Mfloat *steptl, Mint *icode,
                         Mfloat xp[], Mfloat *fp, Mfloat gp[],
                         Mfloat sc[], Mint *mxtake, Mfloat *rnwtnl,
                         Mfloat *epsfcn, Mint *nfcn, Mint *ngrad); 
static void      l_u17ng(Mfloat (*fcn)(Mint, Mfloat[]),
                         void (*grad) (Mint, Mfloat[], Mfloat[]),
                         Mint *fdiff, Mint *n, Mfloat xc[],
                         Mfloat *fc, Mfloat gc[], Mfloat sn[],
                         Mfloat xscale[], Mfloat *stepmx,
                         Mfloat *steptl, Mint *icode, Mfloat xp[],
                         Mfloat *fp, Mfloat gp[], Mfloat sc[],
                         Mint *mxtake, Mfloat *rnwtnl,
                         Mfloat *epsfcn, Mint *nfcn, Mint *ngrad);
static void      l_u18nf(Mint *icode);    
static void      l_u19nf(Mint *icode);
static void      l_cdgrd(Mfloat (*fcn)(Mint, Mfloat[]), Mint *n,
                         Mfloat xc[], Mfloat xscale[],
                         Mfloat *epsfcn, Mfloat gc[]);
static void      l_fdgrd(Mfloat (*fcn)(Mint, Mfloat[]), Mint *n,
                         Mfloat xc[], Mfloat xscale[], Mfloat *fc,
                         Mfloat *epsfcn, Mfloat gc[]);
void             imsl_srot(Mint n, Mfloat sx[], Mint incx,
                           Mfloat sy[], Mint incy, Mfloat c,
                           Mfloat s);
void             imsl_srotg(Mfloat *sa, Mfloat *sb, Mfloat *sc,
                            Mfloat *ss);
#else
static VA_LIST_HACK   l_min_uncon_multivar();
static void      l_u2inf();
static void      l_u2ing();
static void      l_u3inf(); 
static void      l_u3ing();
static void      l_u5inf(); 
static void      l_u6inf();  
static void      l_u8inf();  
static void      l_u9inf();
static void      l_u10nf(); 
static void      l_u11nf(); 
static void      l_u13nf(); 
static void      l_u14nf();  
static void      l_u15nf();  
static void      l_u17nf(); 
static void      l_u17ng();
static void      l_u18nf();    
static void      l_u19nf();
static void      l_cdgrd();
static void      l_fdgrd();
void             imsl_srot();
void             imsl_srotg();
#endif

static Mfloat       *lv_value;

#ifdef ANSI
#if defined(COMPUTER_HP97C)
Mfloat *imsl_f_min_uncon_multivar(Mfloat (*fcn) (Mint, Mfloat*),
                                  Mint n, ...)
#else
Mfloat *imsl_f_min_uncon_multivar(Mfloat (*fcn) (Mint, Mfloat[]),
                                  Mint n, ...)
#endif
#else
Mfloat *imsl_f_min_uncon_multivar(fcn, n, va_alist)
    Mfloat (*fcn) ();
    Mint n;
    va_dcl
#endif
{
    va_list argptr;

    VA_START(argptr, n);
    E1PSH("imsl_f_min_uncon_multivar", "imsl_d_min_uncon_multivar");
    lv_value = NULL;
    IMSL_CALL(l_min_uncon_multivar(fcn, n, argptr));
    va_end(argptr);
    E1POP("imsl_f_min_uncon_multivar", "imsl_d_min_uncon_multivar"); 
    return lv_value;
}


#ifdef ANSI
#if defined(COMPUTER_HP97C)
static VA_LIST_HACK l_min_uncon_multivar(Mfloat (*fcn) (Mint, Mfloat*),
                                    Mint n, va_list argptr)
#else
static VA_LIST_HACK l_min_uncon_multivar(Mfloat (*fcn) (Mint, Mfloat[]),
                                    Mint n, va_list argptr)
#endif
#else
static VA_LIST_HACK l_min_uncon_multivar(fcn, n, argptr)
    Mfloat        (*fcn) ();
    Mint          n;
    va_list       argptr;
#endif
{
    Mint          i, code;
    Mint          arg_number    = 2;
    Mint          maxitn        = 100;
    Mint          maxfcn        = 400;
    Mint          maxgrad       = 400;
    Mint          ihess         = 0;
    Mint          ndigit;
    Mfloat        *xguess       = NULL;
    Mfloat        *xguess_float = NULL;
    Mfloat        *xscale       = NULL;
    Mfloat        *xscale_float = NULL;
    Mfloat        *work         = NULL;
    Mfloat        fscale        = F_ONE;
    Mfloat        grad_tol;
    Mfloat        step_tol;
    Mfloat        fcn_tol;
    Mfloat        max_step = -9999.0e0;
    Mfloat        eps;
    Mfloat        eps_onet;
    Mfloat        eps_twot;
    Mint          user_xguess   = 0;
    Mint          user_xscale   = 0;
    Mfloat        *fvalue       = NULL;
    Mfloat        f;
#ifdef ANSI
    void          (*grad)(Mint, Mfloat[], Mfloat[]);
#else
    void          (*grad)();
#endif
    Mint          user_gradient = 0;
    Mint          fvalue_user   = 0;
    Mint          user_return   = 0;

    eps  = imsl_amach(4);
    eps_onet = pow(eps, F_ONE/F_THREE);
    eps_twot = pow(eps, F_TWO/F_THREE);

#ifdef DOUBLE
    grad_tol = eps_onet; 
    fcn_tol  = imsl_f_max(1.0e-20, eps_twot);
#else
    grad_tol = sqrt(eps);
    fcn_tol  = imsl_f_max(1.0e-10, eps_twot);
#endif
    ndigit   = (Mint)(-log10(eps) + 0.1e0);
    step_tol = eps_twot;

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
            case IMSL_GRAD:
                arg_number++;
                grad = va_arg(argptr, void*);
                user_gradient = 1;
                break;
            case IMSL_XSCALE:
                arg_number++;
                xscale = va_arg(argptr, Mfloat*);
                user_xscale = 1;
                break;
            case IMSL_FSCALE:
                arg_number++;
                fscale = (Mfloat) va_arg(argptr, Mdouble);
                break;
            case IMSL_GRAD_TOL:
                arg_number++; 
                grad_tol = (Mfloat) va_arg(argptr, Mdouble);
                break;
            case IMSL_STEP_TOL:
                arg_number++;
                step_tol = (Mfloat) va_arg(argptr, Mdouble);
                break; 
            case IMSL_REL_FCN_TOL: 
                arg_number++; 
                fcn_tol = (Mfloat) va_arg(argptr, Mdouble); 
                break;  
            case IMSL_MAX_STEP:
                arg_number++;
                max_step = (Mfloat) va_arg(argptr, Mdouble);
                break;
            case IMSL_FSCALE_ADR:
                arg_number++;
                fscale = *(va_arg(argptr, Mfloat *));
                break;
            case IMSL_GRAD_TOL_ADR:
                arg_number++; 
                grad_tol = *(va_arg(argptr, Mfloat *));
                break;
            case IMSL_STEP_TOL_ADR:
                arg_number++;
                step_tol = *(va_arg(argptr, Mfloat *));
                break; 
            case IMSL_REL_FCN_TOL_ADR: 
                arg_number++; 
                fcn_tol = *(va_arg(argptr, Mfloat *)); 
                break;  
            case IMSL_MAX_STEP_ADR:
                arg_number++;
                max_step = *(va_arg(argptr, Mfloat *));
                break;
            case IMSL_GOOD_DIGIT:
                arg_number++; 
                ndigit = va_arg(argptr, Mint);
                break;
            case IMSL_MAX_ITN:
                arg_number++; 
                maxitn = va_arg(argptr, Mint);
                break;
            case IMSL_MAX_FCN:
                arg_number++; 
                maxfcn = va_arg(argptr, Mint);
                break;
            case IMSL_MAX_GRAD:
                arg_number++; 
                maxgrad = va_arg(argptr, Mint);
                break;
            case IMSL_INIT_HESSIAN:
                arg_number++; 
                ihess = va_arg(argptr, Mint);
                break;
            case IMSL_RETURN_USER:
                arg_number++;
                lv_value = va_arg(argptr, Mfloat*);
                user_return = 1;
                break;
            case IMSL_FVALUE:
                arg_number++;
                fvalue  = va_arg(argptr, Mfloat*);
                fvalue_user = 1;
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

    if (n <= 0) {
        imsl_e1sti(1, n);
        imsl_ermes(IMSL_TERMINAL, IMSL_N_MUST_BE_POSITIVE);
    }

    if (imsl_n1rty(0)) goto RETURN;

    work   = (Mfloat *) imsl_malloc ((n+8)*n*sizeof(*work));
    xguess_float = (Mfloat *) imsl_malloc (n*sizeof(*xguess_float)); 
    xscale_float = (Mfloat *) imsl_malloc (n*sizeof(*xscale_float));

    if (work==NULL || xguess_float==NULL || xscale_float==NULL){
        imsl_e1stl(1, "n");
        imsl_e1sti(1, n);
        imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
        goto FREE_SPACE;
    }

    if (user_xguess) 
       for (i=0; i<n; i++)  xguess_float[i] = (Mfloat) xguess[i];
    else
       for (i=0; i<n; i++)  xguess_float[i] = F_ZERO;

    if (user_xscale) 
       for (i=0; i<n; i++)
            xscale_float[i] = (Mfloat) xscale[i];
    else
       for (i=0; i<n; i++)  xscale_float[i] = F_ONE;

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
       l_u2ing(fcn, grad, &n, xguess_float, xscale_float, &fscale,
               &grad_tol, &step_tol, &fcn_tol, &max_step, &ndigit,
               &maxitn, &maxfcn, &maxgrad, &ihess, lv_value, &f, 
               work);
    else
       l_u2inf(fcn, &n, xguess_float, xscale_float, &fscale,
               &grad_tol, &step_tol, &fcn_tol, &max_step, &ndigit,
               &maxitn, &maxfcn, &maxgrad, &ihess, lv_value, &f,
               work);

    if (fvalue_user) *fvalue = f;

FREE_SPACE:
    if (work != NULL)          imsl_free (work);
    if (xscale_float != NULL)  imsl_free (xscale_float);
    if (xguess_float != NULL)  imsl_free (xguess_float);

RETURN:
    if (imsl_n1rty(0) == 5) {
       if (!user_return && lv_value != NULL) imsl_free(lv_value);
       lv_value = NULL;
    }
    return (argptr);
}



static struct t_u16nf {
              Mfloat      gradtl, steptl, rftol, aftol, falstl;
              Mint        mxiter, maxfcn, maxgrd, maxhes;
}       lv_u16nf;



/* -----------------------------------------------------------------------
    IMSL Name:  U10NF/DU10NF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Solve (L*TRANS(L))*s = -g for s.

    Usage:      CALL U10NF (N, H, LDH, GC, SNWTN)

    Arguments:
       N      - Length of the vectors GC, SNWTN. (Input)
       H      - N by N matrix containing the Cholesky factor of the
                Hessian in the lower triangle and diagonal.  (Input)
       LDH    - Leading dimension of H exactly as specified in the
                dimension statement of the calling program.  (Input)
       GC     - Vector of length N containing the current gradient.
                (Input)
       SNWTN  - Vector of length N containing Newton's step.  (Output)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_u10nf(Mint *n, Mfloat *h, Mint *ldh, Mfloat gc[], Mfloat snwtn[])
#else
static void l_u10nf(n, h, ldh, gc, snwtn)
	Mint            *n;
	Mfloat          *h;
	Mint            *ldh;
	Mfloat          gc[], snwtn[];
#endif
{
#define H(I_,J_)	(h+(I_)*(*ldh)+(J_))


	l_u13nf(n, h, ldh, gc, snwtn);
	/* Solve L**Ts=y */
	l_u14nf(n, h, ldh, snwtn, snwtn);
	sscal(*n, -F_ONE, snwtn, 1);

	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  U11NF/DU11NF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Compute X = Y**K * Z.

    Usage:      CALL U11NF (N, Y, K, Z, X)

    Arguments:
       N      - Length of the vectors X, Y and Z.  (Input)
       Y      - Vector of length N.  (Input)
       K      - Integer specifying the exponent.  (Input)
       Z      - Vector of length N.  (Input)
       X      - Vector of length N.  (Output)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_u11nf(Mint *n, Mfloat y[], Mint *k, Mfloat z[], Mfloat x[])
#else
static void l_u11nf(n, y, k, z, x)
	Mint            *n;
	Mfloat           y[];
	Mint            *k;
	Mfloat           z[], x[];
#endif
{
	Mint             i;


	if (*k < 0) {
		if (*k == -1) {
			for (i = 1; i <= *n; i++) {
				x[i - 1] = z[i - 1] / y[i - 1];
			}
		} else {
			for (i = 1; i <= *n; i++) {
				x[i - 1] = z[i - 1] / (imsl_fi_power(y[i - 1], -*k));
			}
		}
	} else {
		if (*k == 1) {
			for (i = 1; i <= *n; i++) {
				x[i - 1] = z[i - 1] * y[i - 1];
			}
		} else {
			for (i = 1; i <= *n; i++) {
				x[i - 1] = imsl_fi_power(y[i - 1], *k) * z[i - 1];
			}
		}
	}

	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  U13NF/DU13NF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Solve L*s = g for s.

    Usage:      CALL U13NF (N, H, LDH, GC, SNWTN)

    Arguments:
       N      - Length of the vectors GC, SNWTN.  (Input)
       H      - N by N matrix containing the Cholesky factor of the
                Hessian in the lower triangle and diagonal.  (Input)
       LDH    - Leading dimension of H exactly as specified in the
                dimension statement of the calling program.  (Input)
       GC     - Vector of length N containing the current gradient.
                (Input)
       SNWTN  - Vector of length N containing the solution.  (Output)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_u13nf(Mint *n, Mfloat *h, Mint *ldh, Mfloat gc[], Mfloat snwtn[])
#else
static void l_u13nf(n, h, ldh, gc, snwtn)
	Mint            *n;
	Mfloat          *h;
	Mint            *ldh;
	Mfloat           gc[], snwtn[];
#endif
{
#define H(I_,J_)	(h+(I_)*(*ldh)+(J_))
	Mint             i;
	Mfloat           sum;


	snwtn[0] = gc[0] / *H(0, 0);
	for (i = 2; i <= *n; i++) {
		sum = imsl_sdot(i - 1, H(0, i - 1), *ldh, snwtn, 1);
		snwtn[i - 1] = (gc[i - 1] - sum) / *H(i - 1, i - 1);
	}

	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  U14NF/DU14NF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Solve TRANS(L)*s = y for s.

    Usage:      CALL U14NF (N, H, LDH, Y, SNWTN)

    Arguments:
       N      - Length of the vector SNWTN.  (Input)
       H      - N by N matrix containing the Cholesky factor of the
                Hessian in the lower triangle and diagonal.  (Input)
       LDH    - Leading dimension of H exactly as specified in the
                dimension statement of the calling program.  (Input)
       Y      - Vector of length N containing the right-hand-side.
                (Input)
       SNWTN  - Vector of length N containing the solution.  (Output)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_u14nf(Mint *n, Mfloat *h, Mint *ldh, Mfloat y[], Mfloat snwtn[])
#else
static void l_u14nf(n, h, ldh, y, snwtn)
	Mint            *n;
	Mfloat          *h;
	Mint            *ldh;
	Mfloat           y[], snwtn[];
#endif
{
#define H(I_,J_)	(h+(I_)*(*ldh)+(J_))
	Mint             i;
	Mfloat           sum;


	snwtn[*n - 1] = y[*n - 1] / *H(*n - 1, *n - 1);
	for (i = *n - 1; i >= 1; i--) {
		sum = imsl_sdot(*n - i, H(i - 1, i), 1, &snwtn[i], 1);
		snwtn[i - 1] = (y[i - 1] - sum) / *H(i - 1, i - 1);
	}

	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  U15NF/DU15NF  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    QR factorization update after a rank-1 matrix is added.

    Usage:      CALL U15NF (METHOD, N, U, V, Q, LDQ, A, LDA)

    Arguments:
       METHOD - Integer flag.  (Input)
                METHOD = 1 if the orthogonal matrix Q is to be updated,
                METHOD = 2 if only R is to be updated.
       N      - Dimension of the problem.  (Input)
       U      - Vector of length N determining the rank-1 matrix to be
                added.  (Input)
       V      - Vector of length N determining the rank-1 matrix to be
                added.  (Input)
       Q      - N by N matrix used only when METHOD = 1.  (Input/Output)
                On input, Q contains the orthogonal matrix Q.
                On output, Q contains the updated orthogonal matrix Q(+).
       LDQ    - Leading dimension of Q exactly as specified in the
                dimension statement of the calling program.  (Input)
       A      - N by N matrix.  (Input/Output)
                On input, the upper triangle of A contains the upper
                   triangular matrix R.
                On output, the upper triangle of A contains the updated
                   upper triangular matrix R(+).
       LDA    - Row dimension of A exactly as specified in the dimension
                statement of the calling program.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_u15nf(Mint *method, Mint *n, Mfloat u[], Mfloat v[], Mfloat *q,
                    Mint *ldq, Mfloat *a, Mint *lda)
#else
static void l_u15nf(method, n, u, v, q, ldq, a, lda)
	Mint            *method, *n;
	Mfloat           u[], v[], *q;
	Mint            *ldq;
	Mfloat          *a;
	Mint            *lda;
#endif
{
#define Q(I_,J_)	(q+(I_)*(*ldq)+(J_))
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
	Mint             _l0, _l1, i, k;
	Mfloat           _f0, sc, ss, temp1, temp2;


	k = *n;
	/* DO WHILE (U(K)=0.0 and K>1) */
L_10000:
	if (!(u[k - 1] == F_ZERO && k > 1))
		goto L_10002;
	k -= 1;
	/* ENDWHILE */

	goto L_10000;
L_10002:
	;
	/*
	 * TRANSFORM R + U*V**T TO UPPER HESSENBERG
	 */
	for (i = k; i >= 2; i--) {
		imsl_srotg(&u[i - 2], &u[i - 1], &sc, &ss);
                _l0 = k - i + 2;
		imsl_srot(_l0, A(i - 2, i - 2), *lda, A(i - 2, i - 1),
			  *lda, sc, ss);
		if (*method == 1) {
                        _l0 = 1;
                        _l1 = 1;
                        _f0 = -ss;
			imsl_srot(*n, Q(i - 2, 0), _l0, Q(i - 1, 0), _l1,
				  sc, _f0);
		}
	}
	/*
	 * ADD THE REMAINING CONTRIBUTION OF U*(V**T) TO THE FIRST ROW OF R
	 */
	saxpy(*n, u[0], v, 1, A(0, 0), *lda);
	/*
	 * TRANSFORM UPPER HESSENBERG MATRIX TO UPPER TRIANGULAR
	 */
	for (i = 1; i <= (k - 1); i++) {
		temp1 = *A(i - 1, i - 1);
		temp2 = *A(i - 1, i);
		imsl_srotg(&temp1, &temp2, &sc, &ss);
                _l0 = *n - i + 1;
		imsl_srot(_l0, A(i - 1, i - 1), *lda, A(i - 1, i), *lda, 
                          sc, ss);
		if (*method == 1) {
                        _l0 = 1;
                        _l1 = 1;
                        _f0 = -ss;
			imsl_srot(*n, Q(i - 1, 0), _l0, Q(i, 0), _l1, sc, _f0);
		}
	}

	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  U17NF/DU17NF  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    November 1, 1984

    Purpose:    Find a next iterate by a line search for the minimization
                problem.

    Usage:      CALL U17NF (FCN, N, XC, FC, GC, SN, XSCALE, STEPMX,
                            STEPTL, ICODE, XP, FP, SC, MXTAKE, NFCN,
                            NGRAD)

    Arguments
       FCN    - A user supplied subroutine to evaluate the function at a
                point X.  FCN must be declared EXTERNAL in the calling
                program and must have the form,
                   SUBROUTINE FCN (N, X, F)
                   INTEGER    N
                   REAL       X(*), F
                      .
                      .
                   RETURN
                   END
                X should not be changed by FCN.
       N      - Dimension of the problem.  (Input)
       XC     - Real vector of length N containing the current iterate.
                   (Input)
       FC     - Real scalar containing the function value at XC.  (Input)
       GC     - Real vector of length N containing the gradient at XC.
                   (Input)
       SN     - Real vector of length N containing a descent direction
                which is the Newton's direction for most cases.  (Input)
       XSCALE - Real vector of length N containing the diagonal scaling
                matrix for the variables.  (Input)
       STEPMX - Real scalar containing the MAXimum allowable step size.
                   (Input)
       STEPTL - Real scalar containing the relative step size at which
                successive iterates are considered close enough to stop
                the algorithm.  (Input)
       ICODE  - Integer return code.  (Output)
                ICODE = 0 means a satisfactory new iterate is found.
                ICODE = 1 means the routine failed to locate a new point
                          sufficiently different from the current point.
       XP     - Real vector of length N containing the new iterate.
                   (Output)
       FP     - Real scalar containing the function value at XP. (Output)
       SC     - Real vector of length N containing the step taken.
                   (Output)
       MXTAKE - Logical variable.  (Output)
                MXTAKE = .TRUE. if a step of MAXimum length was taken.
                MXTAKE = .FALSE. otherwise.
       NFCN   - Number of function evaluations.  (Input/Output)
       NGRAD  - Number of gradient evaluations.  (Input/Output)

    Remarks:
       This is based on Algorithm A6.3.1, imsl_page 325, Dennis-Schnabel book.
       The algorithm performs a bactracking linesearch with an
       alpha-condition.

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied to
                this code.  No other warranty, expressed or implied, is
                applicable.

  -----------------------------------------------------------------------
 */
#define	ALPHA	1.0e-4

#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_u17nf(Mfloat (*fcn)(Mint, Mfloat*), Mint *fdiff, Mint *n,
                    Mfloat xc[], Mfloat *fc, Mfloat gc[], Mfloat sn[],
                    Mfloat xscale[], Mfloat *stepmx, Mfloat *steptl,
                    Mint *icode, Mfloat xp[], Mfloat *fp, Mfloat gp[],
                    Mfloat sc[], Mint *mxtake, Mfloat *rnwtnl, Mfloat *epsfcn,
                    Mint *nfcn, Mint *ngrad)
#else
static void l_u17nf(Mfloat (*fcn)(Mint, Mfloat[]), Mint *fdiff, Mint *n,
                    Mfloat xc[], Mfloat *fc, Mfloat gc[], Mfloat sn[],
                    Mfloat xscale[], Mfloat *stepmx, Mfloat *steptl,
                    Mint *icode, Mfloat xp[], Mfloat *fp, Mfloat gp[],
                    Mfloat sc[], Mint *mxtake, Mfloat *rnwtnl, Mfloat *epsfcn,
                    Mint *nfcn, Mint *ngrad)
#endif
#else
static void l_u17nf(fcn, fdiff, n, xc, fc, gc, sn, xscale, stepmx, steptl, 
                    icode, xp, fp, gp, sc, mxtake, rnwtnl, epsfcn, nfcn, ngrad)
	Mfloat         (*fcn) ();
	Mint           *fdiff;
	Mint           *n;
	Mfloat         xc[], *fc, gc[], sn[], xscale[], *stepmx, *steptl;
	Mint           *icode;
	Mfloat         xp[], *fp, gp[], sc[];
	Mint           *mxtake;
	Mfloat         *rnwtnl, *epsfcn;
	Mint           *nfcn, *ngrad;
#endif
{
	Mint            i;
	Mfloat          a, adiffa, aincra, alamda, alomda, amaxla, aminla,
	                b, betb, disc, factor, fhi, flo, pfp, plamda, rellen,
	                slope, slopn, t1, t2, t3, t4, t5, temp, temp1,
	                temp2, templ;
        Mint            IMSLFALSE = 0, IMSLTRUE = 1;
        Mint            a_changed = 0, b_changed = 0;


	imsl_e1psh("U17NF ");

	*mxtake = IMSLFALSE;
	*icode = 2;
	/* INITIAL STEP TO TRY IS NEWTON STEP */
	scopy(*n, sn, 1, sc, 1);
	/* NEWTON STEP LONGER THAN MAXSTEP */
	if (*rnwtnl > *stepmx) {
		factor = *stepmx / *rnwtnl;
		sscal(*n, factor, sc, 1);
		*rnwtnl = *stepmx;
	}
	/*
	 * RELLEN = RELATIVE LENGTH OF SC AS CALCULATED IN STOPPING RULE
	 */
	slope = imsl_sdot(*n, gc, 1, sc, 1);
	rellen = F_ZERO;
	for (i = 1; i <= *n; i++) {
		temp = fabs(sc[i - 1]) / (imsl_f_max(fabs(xc[i - 1]), F_ONE / xscale[i - 1]));
		rellen = imsl_f_max(rellen, temp);
	}
	aminla = imsl_f_max(*steptl, *epsfcn) / rellen;
	alamda = F_ONE;
	/* MAIN LOOP */
L_20:
	;
	scopy(*n, xc, 1, xp, 1);
	saxpy(*n, alamda, sc, 1, xp, 1);
	imsl_e1usr("ON");
	*fp = (*fcn) (*n, xp);
	imsl_e1usr("OFF");
	*nfcn += 1;
	if (*fp <= (*fc + ALPHA * alamda * slope)) {
		/* XP SATISFIES ALPHA CONDITION */
		if (*fdiff) {
			l_fdgrd(fcn, n, xp, xscale, fp, epsfcn, gp);
		} else {
			l_cdgrd(fcn, n, xp, xscale, epsfcn, gp);
		}
		*ngrad += 1;
		betb = 0.9e0;
		slopn = imsl_sdot(*n, gp, 1, sc, 1);
		if (slopn < betb * slope) {
			if (alamda == F_ONE && *rnwtnl < *stepmx) {
				amaxla = *stepmx / *rnwtnl;
		L_30:
				plamda = alamda;
				pfp = *fp;
				alamda = imsl_f_min(F_TWO * alamda, amaxla);
				scopy(*n, xc, 1, xp, 1);
				saxpy(*n, alamda, sc, 1, xp, 1);
				imsl_e1usr("ON");
				*fp = (*fcn) (*n, xp);
				imsl_e1usr("OFF");
				*nfcn += 1;
				if (*fp <= *fc + ALPHA * alamda * slope) {
				     if (*fdiff) {
					 l_fdgrd(fcn, n, xp, xscale, fp, epsfcn,
                                                 gp);
				      } else {
					  l_cdgrd(fcn, n, xp, xscale, epsfcn,
                                                  gp);
				      }
				      *ngrad += 1;
				      slopn = imsl_sdot(*n, gp, 1, sc, 1);
				}
				if ((*fp <= *fc + alamda * ALPHA * slope
                                     && alamda < amaxla) &&
                                     slopn < betb * slope)
					goto L_30;
			}
			if (alamda < F_ONE || (alamda > F_ONE && *fp > *fc +
                            ALPHA * alamda * slope)) {
				alomda = imsl_f_min(alamda, plamda);
				adiffa = fabs(plamda - alamda);
				if (alamda < plamda) {
					flo = *fp;
					fhi = pfp;
				} else {
					flo = pfp;
					fhi = *fp;
				}
		L_40:
				temp1 = slopn * adiffa * adiffa;
				temp2 = flo + slopn * adiffa;
				aincra = -temp1 / (2 * (fhi - temp2));
				if (aincra < 0.2 * adiffa)
					aincra = 0.2 * adiffa;
				alamda = alomda + aincra;
				scopy(*n, xc, 1, xp, 1);
				saxpy(*n, alamda, sc, 1, xp, 1);
				imsl_e1usr("ON");
				*fp = (*fcn) (*n, xp);
				imsl_e1usr("OFF");
				*nfcn += 1;
				if (*fp > *fc + ALPHA * alamda * slope) {
					adiffa = aincra;
					fhi = *fp;
				} else {
					if (*fdiff) {
						l_fdgrd(fcn, n, xp, xscale, fp,
                                                        epsfcn, gp);
					} else {
						l_cdgrd(fcn, n, xp, xscale,
                                                        epsfcn, gp);
					}
					*ngrad += 1;
					slopn = imsl_sdot(*n, gp, 1, sc, 1);
					if (slopn < betb * slope) {
						alomda = alamda;
						adiffa -= aincra;
						flo = *fp;
					}
				}
				if (slopn < betb * slope && adiffa >= aminla)
					goto L_40;
				if (slopn < betb * slope) {
					*fp = flo;
					scopy(*n, xc, 1, xp, 1);
					saxpy(*n, alomda, sc, 1, xp, 1);
				}
			}
		}
		*icode = 0;
		if (alamda ** rnwtnl > 0.99 ** stepmx)
			*mxtake = IMSLTRUE;
	} else if (alamda < aminla) {
		/*
		 * NO SATISFACTORY XP CAN BE FOUND SUFFICIENTLY DISTINCT FROM
		 * XC
		 */
		*icode = 1;
		scopy(*n, xc, 1, xp, 1);
	} else {
		/* REDUCE THE STEP SIZE LAMBDA */
		if (alamda == F_ONE) {
			templ = -slope / (F_TWO * (*fp - *fc - slope));
		} else {
			t1 = *fp - *fc - alamda * slope;
			t2 = pfp - *fc - plamda * slope;
			t3 = F_ONE / (alamda - plamda);
			t4 = alamda * alamda;
			t5 = plamda * plamda;
#if defined(COMPUTER_VAX) || defined(COMPUTER_ALFAC_IEEE)
                        if (fabs(t3) > fabs(imsl_amach(2)/(t1 / t4 - t2 / t5)))
                                t3 = imsl_amach(2)/(t1 / t4 - t2 / t5);
#endif
			a = t3 * (t1 / t4 - t2 / t5);
			b = t3 * (-t1 * plamda / t4 + t2 * alamda / t5);
			/*
			 * DISC WILL BE NONNEGATIVE AS LONG AS ALPHA IS .LT.
			 * 0.25
			 */
#if defined(COMPUTER_VAX) || defined(COMPUTER_ALFAC_IEEE)
                        {int asign = a > 0 ? 1 : -1;
                        int bsign = b > 0 ? 1 : -1;
#if defined(COMPUTER_ALFAC_IEEE) && defined(DOUBLE)
			if (fabs(b) > 1.340781e+154) {
				b =bsign*(1.340781e+154);
#else
                        if (fabs(b) > sqrt(imsl_amach(2))) {
				b =bsign*sqrt(imsl_amach(2));
#endif
				b_changed = 1;
			}
                        if (fabs(a) > imsl_amach(2)/(F_THREE*slope)) {
                                a = .5*asign*imsl_amach(2)/(F_THREE*slope);
                                a_changed = 1;
                        }
                        }
#endif
#if defined(COMPUTER_VAX) || defined(COMPUTER_ALFAC_IEEE)
                        if (a_changed && b_changed && (a*slope < 0))
                                disc = .5*b*b - 1.5*a*slope;
                        else
                                disc = b * b - F_THREE * a * slope;
#else
			disc = b * b - F_THREE * a * slope;
#endif
			if (a == F_ZERO) {
				templ = -slope / (F_TWO * b);
			} else {
				templ = (-b + sqrt(disc)) / (F_THREE * a);
			}
			if (templ > F_HALF * alamda)
				templ = F_HALF * alamda;
		}
		plamda = alamda;
		pfp = *fp;
		if (templ <= 0.1e0 * alamda) {
			alamda *= 0.1e0;
		} else {
			alamda = templ;
		}
	}
	if (*icode >= 2)
		goto L_20;
	/* OUTPUT THE STEP TAKEN */
	for (i = 1; i <= *n; i++) {
		sc[i - 1] = xp[i - 1] - xc[i - 1];
	}

	imsl_e1pop("U17NF ");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  U17NG/DU17NG  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    November 1, 1984

    Purpose:    Find a next iterate by a line search for the minimization
                problem.

    Usage:      CALL U17NG (FCN, N, XC, FC, GC, SN, XSCALE, STEPMX,
                            STEPTL, ICODE, XP, FP, SC, MXTAKE, NFCN,
                            NGRAD)

    Arguments
       FCN    - A user supplied subroutine to evaluate the function at a
                point X.  FCN must be declared EXTERNAL in the calling
                program and must have the form,
                   SUBROUTINE FCN (N, X, F)
                   INTEGER    N
                   REAL       X(*), F
                      .
                      .
                   RETURN
                   END
                X should not be changed by FCN.
       GRAD   - A user supplied subroutine to evaluate the gradient at a
                point X.  GRAD must be declared EXTERNAL in the calling
                program and must have the form,
                   SUBROUTINE GRAD (N, X, G)
                   INTEGER    N
                   REAL       X(*), G(*)
                      .
                      .
                   RETURN
                   END
                X should not be changed by GRAD.
       N      - Dimension of the problem.  (Input)
       XC     - Real vector of length N containing the current iterate.
                   (Input)
       FC     - Real scalar containing the function value at XC.  (Input)
       GC     - Real vector of length N containing the gradient at XC.
                   (Input)
       SN     - Real vector of length N containing a descent direction
                which is the Newton's direction for most cases.  (Input)
       XSCALE - Real vector of length N containing the diagonal scaling
                matrix for the variables.  (Input)
       STEPMX - Real scalar containing the MAXimum allowable step size.
                   (Input)
       STEPTL - Real scalar containing the relative step size at which
                successive iterates are considered close enough to stop
                the algorithm.  (Input)
       ICODE  - Integer return code.  (Output)
                ICODE = 0 means a satisfactory new iterate is found.
                ICODE = 1 means the routine failed to locate a new point
                          sufficiently different from the current point.
       XP     - Real vector of length N containing the new iterate.
                   (Output)
       FP     - Real scalar containing the function value at XP. (Output)
       SC     - Real vector of length N containing the step taken.
                   (Output)
       MXTAKE - Logical variable.  (Output)
                MXTAKE = .TRUE. if a step of MAXimum length was taken.
                MXTAKE = .FALSE. otherwise.
       NFCN   - Number of function evaluations.  (Input/Output)
       NGRAD  - Number of gradient evaluations.  (Input/Output)

    Remarks:
       This is based on Algorithm A6.3.1, imsl_page 325, Dennis-Schnabel book.
       The algorithm performs a bactracking linesearch with an
       alpha-condition.

    Copyright:  1984 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied to
                this code.  No other warranty, expressed or implied, is
                applicable.

  -----------------------------------------------------------------------
 */
/* FDIFF and EPSFCN are not used here, but leave the calling sequence intact.*/
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_u17ng(Mfloat (*fcn)(Mint, Mfloat*), void (*grad) (Mint, Mfloat*,
                    Mfloat*), Mint *fdiff, Mint *n, Mfloat xc[], Mfloat *fc,
                    Mfloat gc[], Mfloat sn[], Mfloat xscale[], Mfloat *stepmx,
                    Mfloat *steptl, Mint *icode, Mfloat xp[], Mfloat *fp,
                    Mfloat gp[], Mfloat sc[], Mint *mxtake, Mfloat *rnwtnl,
                    Mfloat *epsfcn, Mint *nfcn, Mint *ngrad)
#else
static void l_u17ng(Mfloat (*fcn)(Mint, Mfloat[]), void (*grad) (Mint, Mfloat[],
                    Mfloat[]), Mint *fdiff, Mint *n, Mfloat xc[], Mfloat *fc,
                    Mfloat gc[], Mfloat sn[], Mfloat xscale[], Mfloat *stepmx,
                    Mfloat *steptl, Mint *icode, Mfloat xp[], Mfloat *fp,
                    Mfloat gp[], Mfloat sc[], Mint *mxtake, Mfloat *rnwtnl,
                    Mfloat *epsfcn, Mint *nfcn, Mint *ngrad)
#endif
#else
static void l_u17ng(fcn, grad, fdiff, n, xc, fc, gc, sn, xscale,
	   stepmx, steptl, icode, xp, fp, gp, sc, mxtake, rnwtnl, epsfcn,
	   nfcn, ngrad)
	Mfloat          (*fcn) ();
	Mvoid		(*grad) ();
	Mint            *fdiff;
	Mint            *n;
	Mfloat           xc[], *fc, gc[], sn[], xscale[], *stepmx, *steptl;
	Mint            *icode;
	Mfloat           xp[], *fp, gp[], sc[];
	Mint            *mxtake;
	Mfloat          *rnwtnl, *epsfcn;
	Mint            *nfcn, *ngrad;
#endif
{
	Mint            i;
	Mfloat          a, adiffa, aincra, alamda, alomda, amaxla, aminla,
	                b, betb, disc, factor, fhi, flo, pfp, plamda, rellen,
	                slope, slopn, t1, t2, t3, t4, t5, temp, temp1,
	                temp2, templ;
        Mint            IMSLFALSE = 0, IMSLTRUE = 1;
        Mint            a_changed = 0, b_changed = 0;


	imsl_e1psh("U17NG ");

	*mxtake = IMSLFALSE;
	*icode = 2;
	/* INITIAL STEP TO TRY IS NEWTON STEP */
	scopy(*n, sn, 1, sc, 1);
	/* NEWTON STEP LONGER THAN MAXSTEP */
	if (*rnwtnl > *stepmx) {
		factor = *stepmx / *rnwtnl;
		sscal(*n, factor, sc, 1);
		*rnwtnl = *stepmx;
	}
	/*
	 * RELLEN = RELATIVE LENGTH OF SC AS CALCULATED IN STOPPING RULE
	 */
	slope = imsl_sdot(*n, gc, 1, sc, 1);
	rellen = F_ZERO;
	for (i = 1; i <= *n; i++) {
		temp = fabs(sc[i - 1]) / (imsl_f_max(fabs(xc[i - 1]), F_ONE / xscale[i - 1]));
		rellen = imsl_f_max(rellen, temp);
	}
	aminla = *steptl / rellen;
	alamda = F_ONE;
	/* MAIN LOOP */
L_20:
	;
	scopy(*n, xc, 1, xp, 1);
	saxpy(*n, alamda, sc, 1, xp, 1);
	imsl_e1usr("ON");
	*fp = (*fcn) (*n, xp);
	imsl_e1usr("OFF");
	*nfcn += 1;
	if (*fp <= (*fc + ALPHA * alamda * slope)) {
		/* XP SATISFIES ALPHA CONDITION */
		imsl_e1usr("ON");
		(*grad) (*n, xp, gp);
		imsl_e1usr("OFF");
		*ngrad += 1;
		betb = 0.9e0;
		slopn = imsl_sdot(*n, gp, 1, sc, 1);
		if (slopn < betb * slope) {
			if (alamda == F_ONE && *rnwtnl < *stepmx) {
				amaxla = *stepmx / *rnwtnl;
		L_30:
				plamda = alamda;
				pfp = *fp;
				alamda = imsl_f_min(F_TWO * alamda, amaxla);
				scopy(*n, xc, 1, xp, 1);
				saxpy(*n, alamda, sc, 1, xp, 1);
				imsl_e1usr("ON");
				*fp = (*fcn) (*n, xp);
				imsl_e1usr("OFF");
				*nfcn += 1;
				if (*fp <= *fc + ALPHA * alamda * slope) {
					imsl_e1usr("ON");
					(*grad) (*n, xp, gp);
					imsl_e1usr("OFF");
					*ngrad += 1;
					slopn = imsl_sdot(*n, gp, 1, sc, 1);
				}
				if ((*fp <= *fc + alamda * ALPHA * slope && alamda < amaxla) &&
				    slopn < betb * slope)
					goto L_30;
			}
			if (alamda < F_ONE || (alamda > F_ONE && *fp > *fc + ALPHA *
					       alamda * slope)) {
				alomda = imsl_f_min(alamda, plamda);
				adiffa = fabs(plamda - alamda);
				if (alamda < plamda) {
					flo = *fp;
					fhi = pfp;
				} else {
					flo = pfp;
					fhi = *fp;
				}
		L_40:
				temp1 = slopn * adiffa * adiffa;
				temp2 = flo + slopn * adiffa;
				aincra = -temp1 / (2 * (fhi - temp2));
				if (aincra < 0.2 * adiffa)
					aincra = 0.2 * adiffa;
				alamda = alomda + aincra;
				scopy(*n, xc, 1, xp, 1);
				saxpy(*n, alamda, sc, 1, xp, 1);
				imsl_e1usr("ON");
				*fp = (*fcn) (*n, xp);
				imsl_e1usr("OFF");
				*nfcn += 1;
				if (*fp > *fc + ALPHA * alamda * slope) {
					adiffa = aincra;
					fhi = *fp;
				} else {
					imsl_e1usr("ON");
					(*grad) (*n, xp, gp);
					imsl_e1usr("OFF");
					*ngrad += 1;
					slopn = imsl_sdot(*n, gp, 1, sc, 1);
					if (slopn < betb * slope) {
						alomda = alamda;
						adiffa -= aincra;
						flo = *fp;
					}
				}
				if (slopn < betb * slope && adiffa >= aminla)
					goto L_40;
				if (slopn < betb * slope) {
					*fp = flo;
					scopy(*n, xc, 1, xp, 1);
					saxpy(*n, alomda, sc, 1, xp, 1);
				}
			}
		}
		*icode = 0;
		if (alamda ** rnwtnl > 0.99 ** stepmx)
			*mxtake = IMSLTRUE;
	} else if (alamda < aminla) {
		/*
		 * NO SATISFACTORY XP CAN BE FOUND SUFFICIENTLY DISTINCT FROM
		 * XC
		 */
		*icode = 1;
		scopy(*n, xc, 1, xp, 1);
	} else {
		/* REDUCE THE STEP SIZE LAMBDA */
		if (alamda == F_ONE) {
			templ = -slope / (F_TWO * (*fp - *fc - slope));
		} else {
			t1 = *fp - *fc - alamda * slope;
			t2 = pfp - *fc - plamda * slope;
			t3 = F_ONE / (alamda - plamda);
			t4 = alamda * alamda;
			t5 = plamda * plamda;
#if defined(COMPUTER_VAX) || defined(COMPUTER_ALFAC_IEEE)
                        if (fabs(t3) > fabs(imsl_amach(2)/(t1 / t4 - t2 / t5)))
                                t3 = imsl_amach(2)/(t1 / t4 - t2 / t5);
#endif
			a = t3 * (t1 / t4 - t2 / t5);
			b = t3 * (-t1 * plamda / t4 + t2 * alamda / t5);
			/*
			 * DISC WILL BE NONNEGATIVE AS LONG AS ALPHA IS .LT.
			 * 0.25
			 */
#if defined(COMPUTER_VAX) || defined(COMPUTER_ALFAC_IEEE)
                        {int asign = a > 0 ? 1 : -1;
                        int bsign = b > 0 ? 1 : -1;
#if defined(COMPUTER_ALFAC_IEEE) && defined(DOUBLE)
			if (fabs(b) > 1.340781e+154) {
				b =bsign*(1.340781e+154);
#else
                        if (fabs(b) > sqrt(imsl_amach(2))) {
                                b =bsign*sqrt(imsl_amach(2));
#endif
                                b_changed = 1;
                        }
                        if (fabs(a) > fabs(imsl_amach(2)/(F_THREE*slope))) {
                                a = .5*asign*imsl_amach(2)/(F_THREE*slope);
                                a_changed = 1;
                        }
                        }
#endif
#if defined(COMPUTER_VAX) || defined(COMPUTER_ALFAC_IEEE)
                        if (a_changed && b_changed && (a*slope < 0))
                                disc = .5*b*b - 1.5*a*slope;
                        else
                                disc = b * b - F_THREE * a * slope;
#else
			disc = b * b - F_THREE * a * slope;
#endif
			if (a == F_ZERO) {
				templ = -slope / (F_TWO * b);
			} else {
				templ = (-b + sqrt(disc)) / (F_THREE * a);
			}
			if (templ > F_HALF * alamda)
				templ = F_HALF * alamda;
		}
		plamda = alamda;
		pfp = *fp;
		if (templ <= 0.1e0 * alamda) {
			alamda *= 0.1e0;
		} else {
			alamda = templ;
		}
	}
	if (*icode >= 2)
		goto L_20;
	/* OUTPUT THE STEP TAKEN */
	for (i = 1; i <= *n; i++) {
		sc[i - 1] = xp[i - 1] - xc[i - 1];
	}

	imsl_e1pop("U17NG ");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  U18NF/DU18NF  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Stopping conditions for unconstrained minimization.

    Usage:      CALL U18NF (ICODE)

    Arguments:
       ICODE  - Integer flag containing an error code.  (Input)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_u18nf(Mint *icode)
#else
static void l_u18nf(icode)
	Mint            *icode;
#endif
{


	if (*icode == 3) {
                imsl_ermes(IMSL_WARNING, IMSL_TOO_MANY_ITN);

	} else if (*icode == 4) {
                imsl_ermes(IMSL_WARNING, IMSL_TOO_MANY_FCN_EVAL);

	} else if (*icode == 5) {
                imsl_ermes(IMSL_WARNING, IMSL_TOO_MANY_GRAD_EVAL);

	} else if (*icode == 7) {
                imsl_ermes(IMSL_WARNING, IMSL_TOO_MANY_HESSIAN_EVAL);

	} else if (*icode == 6) {

/*		(3, 6, " Five consecutive steps of length MAX_STEP have been taken; either the function is unbounded below, or has a finite asymptote in some direction or the maximum allowable step size MAX_STEP is too small.");
*/
                imsl_ermes(IMSL_WARNING, IMSL_UNBOUNDED);
	}
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  U19NF/DU19NF  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Check validity of input to unconstrained minimization.

    Usage:      CALL U19NF (ICODE)

    Arguments:
       ICODE  - Integer flag containing an error code.  (Input)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_u19nf(Mint *icode)
#else
static void l_u19nf(icode)
	Mint            *icode;
#endif
{


	if (*icode == 0) {
                imsl_ermes(IMSL_TERMINAL, IMSL_N_MUST_BE_POSITIVE);

	} else if (*icode == 1) {
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_INEFFICIENT_PROB_SIZE);

	} else if (*icode == 2) {

/*		(6, 2, "The diagonal scaling matrix for the variables must be positive while some of the entries are less than or equal to zero.  The algorithm will use the identity scaling matrix for XSCALE.");
*/
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_SCALING_WITH_IDENT_MATRIX);
	} else if (*icode == 4) {

/*		(6, 4, "The estimate of the number of good digits in the function must be positive while NDIGIT = %(i1) is given.  The algorithm will assume that the function is accurate to the precision of the arithmetic.");
*/
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEED_POSITIVE_NDIGIT);
	} else if (*icode == 5) {

/*		(6, 5, "The maximum number of iterations must be positive while MAXITN = %(i1) is given.  The algorithm will use MAXITN = 100.");
*/
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEED_POSITIVE_MXITER);
	} else if (*icode == 6) {

/*		(6, 6, "The maximum number of function evaluations must be positive while MAXFCN = %(i1) is given.  The algorithm will use MAXFCN = 400.");
*/
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEED_POSITIVE_MAXFCN);
	} else if (*icode == 7) {

/*		(6, 7, "The maximum number of gradient evaluations must be positive while MAXGRAD = %(i1) is given.  The algorithm will use MAXGRAD = 400.");
*/
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_MAXGRAD_VALUE_TOO_SMALL);
	} else if (*icode == 8) {

/*		(6, 8, "The maximum number of Hessian evaluations must be positive while MAXHES = %(i1) is given.  The algorithm will use MAXHES = 100.");
*/
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_MAXHES_VALUE_TOO_SMALL);
	}
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  U2INF/DU2INF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Minimize a function of N variables using a quasi-Newton
                method and a finite difference gradient.

    Usage:      CALL U2INF (FCN, N, XGUESS, XSCALE, FSCALE, IPARAM,
                            RPARAM, X, FVALUE, WK)

    Arguments:  See UMINF/DUMINF.

    Remarks:    See UMINF/DUMINF.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_u2inf(Mfloat (*fcn) (Mint, Mfloat*), Mint *n, Mfloat xguess[],
                    Mfloat xscale[], Mfloat *fscale, Mfloat *grad_tol,
                    Mfloat *step_tol, Mfloat *fcn_tol, Mfloat *max_step,
                    Mint *ndigit, Mint *maxitn, Mint *maxfcn, Mint *maxgrad,
                    Mint *ihess, Mfloat x[], Mfloat *fvalue, Mfloat wk[])
#else
static void l_u2inf(Mfloat (*fcn) (Mint, Mfloat[]), Mint *n, Mfloat xguess[],
                    Mfloat xscale[], Mfloat *fscale, Mfloat *grad_tol,
                    Mfloat *step_tol, Mfloat *fcn_tol, Mfloat *max_step,
                    Mint *ndigit, Mint *maxitn, Mint *maxfcn, Mint *maxgrad,
                    Mint *ihess, Mfloat x[], Mfloat *fvalue, Mfloat wk[])
#endif
#else
static void l_u2inf(fcn, n, xguess, xscale, fscale,grad_tol, step_tol,
                    fcn_tol, max_step, ndigit, maxitn, maxfcn, maxgrad,
                    ihess, x, fvalue, wk)
	Mfloat          (*fcn) ();
	Mint            *n;
	Mfloat           xguess[], xscale[], *fscale, *grad_tol, *step_tol,
                        *fcn_tol, *max_step;
        Mint            *ndigit, *maxitn, *maxfcn, *maxgrad, *ihess;
	Mfloat           x[], *fvalue, wk[];
#endif
{
        Mint             iparam[7];
        Mfloat           rparam[7], tv3, eps;


	imsl_e1psh("U2INF ");

        iparam[0] = 1;
        iparam[1] = *ndigit;
        iparam[2] = *maxitn;
        iparam[3] = *maxfcn;
        iparam[4] = *maxgrad;
        iparam[5] = *ihess;
        iparam[6] = 100;

        eps = imsl_amach(4);
        tv3 = F_TWO/F_THREE;
        rparam[0] = *grad_tol;
        rparam[1] = *step_tol;
        rparam[2] = *fcn_tol;
        rparam[3] = imsl_f_max(1.0e-10, pow(eps, tv3));
        rparam[4] = 100*eps;
        rparam[5] = *max_step;
        rparam[6] = -9999.0e0;

	scopy(*n, xguess, 1, x, 1);
	/*
	* Call unconstrained minimization solver using F.D. gradient
	*/
	l_u3inf(fcn, n, x, xscale, fscale, iparam, rparam, fvalue, &wk[0],
                &wk[*n], &wk[*n * 2], &wk[*n * 3], &wk[*n * 4], &wk[*n * 8], n,
                &wk[*n * 5], &wk[*n * 6], &wk[*n * 7]);

	imsl_e1pop("U2INF ");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  U2ING/DU2ING (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Minimize a function of N variables using a quasi-Newton
                method and a user-supplied gradient.

    Usage:      CALL U2ING (FCN, GRAD, N, XGUESS, XSCALE, FSCALE, IPARAM,
                            RPARAM, X, FVALUE, WK)

    Arguments:  See UMING/DUMING.

    Remarks:    See UMING/DUMING.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_u2ing(Mfloat (*fcn) (Mint, Mfloat*), void (*grad) (Mint,
                    Mfloat*, Mfloat*), Mint *n, Mfloat xguess[],
                    Mfloat xscale[], Mfloat *fscale, Mfloat *grad_tol,
                    Mfloat *step_tol, Mfloat *fcn_tol, Mfloat *max_step,
                    Mint *ndigit, Mint *maxitn, Mint *maxfcn, Mint *maxgrad,
                    Mint *ihess, Mfloat x[], Mfloat *fvalue, Mfloat wk[])
#else
static void l_u2ing(Mfloat (*fcn) (Mint, Mfloat[]), void (*grad) (Mint,
                    Mfloat[], Mfloat[]), Mint *n, Mfloat xguess[],
                    Mfloat xscale[], Mfloat *fscale, Mfloat *grad_tol,
                    Mfloat *step_tol, Mfloat *fcn_tol, Mfloat *max_step,
                    Mint *ndigit, Mint *maxitn, Mint *maxfcn, Mint *maxgrad,
                    Mint *ihess, Mfloat x[], Mfloat *fvalue, Mfloat wk[])
#endif
#else
static void l_u2ing(fcn, grad, n, xguess, xscale, fscale, grad_tol, step_tol,
                    fcn_tol, max_step, ndigit, maxitn, maxfcn, maxgrad, ihess,
                    x, fvalue, wk)
	Mfloat          (*fcn) ();
	Mvoid		(*grad) ();
	Mint            *n;
	Mfloat           xguess[], xscale[], *fscale, *grad_tol, *step_tol,
                        *fcn_tol, *max_step;
        Mint            *ndigit, *maxitn, *maxfcn, *maxgrad, *ihess;
	Mfloat           x[], *fvalue, wk[];
#endif
{
        Mint             iparam[7];
        Mfloat           rparam[7], tv3, eps;


	imsl_e1psh("U2ING ");

        iparam[0] = 1;
        iparam[1] = *ndigit;
        iparam[2] = *maxitn;
        iparam[3] = *maxfcn;
        iparam[4] = *maxgrad;
        iparam[5] = *ihess;
        iparam[6] = 100;

        eps = imsl_amach(4);
        tv3 = F_TWO/F_THREE;
        rparam[0] = *grad_tol;
        rparam[1] = *step_tol;
        rparam[2] = *fcn_tol;
        rparam[3] = imsl_f_max(1.0e-10, pow(eps, tv3));
        rparam[4] = 100*eps;
        rparam[5] = *max_step;
        rparam[6] = -9999.0e0;
	/*
	* CALL UNCONSTRAINED MINIMIZATION SOLVER USING ANALYTIC GRADIENT
	*/
	scopy(*n, xguess, 1, x, 1);
	l_u3ing(fcn, grad, n, x, xscale, fscale, iparam, rparam, fvalue, &wk[0],
                &wk[*n], &wk[*n * 2], &wk[*n * 3], &wk[*n * 4], &wk[*n * 8], n,
                &wk[*n * 5], &wk[*n * 6], &wk[*n * 7]);

	imsl_e1pop("U2ING ");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  U3INF/DU3INF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Driver for unconstrained minimization solver using finite
                difference or analytic gradient.

    Usage:      CALL U3INF (FCN, N, XC, XSCALE, FSCALE, IPARAM, RPARAM,
                            FVALUE, XP, SC, SNWTN, GC, GP, H, LDH, WK1,
                            WK2, WK3)

    Arguments:
       FCN    - A user-supplied subroutine to evaluate the function at a
                point X.  FCN must be declared EXTERNAL in the calling
                program and must have the form,
                   SUBROUTINE FCN (N, X, F)
                   INTEGER    N
                   REAL       X(*), F
                      .
                      .
                   RETURN
                   END
                X should not be changed by FCN.
       N      - Dimension of the problem.  (Input)
       XC     - Real vector of length N containing initial guess on input
                and approximate solution on output.  (Input / Output)
       XSCALE - Real vector of length N containing the diagonal scaling
                matrix for the variables.  (Input)
       FSCALE - Real scalar containing the function scaling.  (Input)
       IPARAM - Integer parameters vector of length 7.  (Input)
                IPARAM(1) = initialization flag
                IPARAM(2) = number of good digits in the function
                IPARAM(3) = maximum number of iterations.
                IPARAM(4) = maximum number of function evaluations
                IPARAM(5) = maximum number of gradient evaluations
                IPARAM(6) = initialization of Hessian parameter
                IPARAM(7) = maximum number of Hessian evaluations
       RPARAM - Real parameters vector of length 7.  (Input)
                RPARAM(1) = scaled gradient tolerance
                RPARAM(2) = scaled step tolerance
                RPARAM(3) = relative function tolerance
                RPARAM(4) = absolute function tolerance
                RPARAM(5) = false convergence tolerance
                RPARAM(6) = maximum allowable step size
                RPARAM(7) = size of initial trust region radius
       FVALUE - Real scalar containing the value of the function at the
                solution.  (Output)
       XP     - Real vector of length N containing the updated point.
                   (Output)
       SC     - Real vector of length N containing the last step taken.
                   (Output)
       SNWTN  - Real vector of length N containing the last Newton step.
                   (Output)
       GC     - Real vector of length N containing an estimate of the
                gradient at the current point.  (Output)
       GP     - Real vector of length N containing an estimate of the
                gradient at the updated point.  (Output)
       H      - Real N by N matrix containing an estimate of the Hessian
                at the approximate solution.  (Output)
       LDH    - Leading dimension of H exactly as specified in the
                dimension statement of the calling program.  (Input)
       WK1    - Real work vector of length N.  (Output)
       WK2    - Real work vector of length N.  (Output)
       WK3    - Real work vector of length N.  (Output)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */

#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_u3inf(Mfloat (*fcn)(Mint, Mfloat*), Mint *n, Mfloat xc[],
                    Mfloat xscale[], Mfloat *fscale, Mint iparam[],
                    Mfloat rparam[], Mfloat *fvalue, Mfloat xp[], Mfloat sc[],
                    Mfloat snwtn[], Mfloat gc[], Mfloat gp[], Mfloat *h, 
                    Mint *ldh, Mfloat wk1[], Mfloat wk2[], Mfloat wk3[])
#else
static void l_u3inf(Mfloat (*fcn)(Mint, Mfloat[]), Mint *n, Mfloat xc[],
                    Mfloat xscale[], Mfloat *fscale, Mint iparam[],
                    Mfloat rparam[], Mfloat *fvalue, Mfloat xp[], Mfloat sc[],
                    Mfloat snwtn[], Mfloat gc[], Mfloat gp[], Mfloat *h, 
                    Mint *ldh, Mfloat wk1[], Mfloat wk2[], Mfloat wk3[])
#endif
#else
static void l_u3inf(fcn, n, xc, xscale, fscale, iparam, rparam, fvalue, xp,
                    sc, snwtn, gc, gp, h, ldh, wk1, wk2, wk3)
	Mfloat          (*fcn) ();
	Mint            *n;
	Mfloat          xc[], xscale[], *fscale;
	Mint            iparam[];
	Mfloat          rparam[], *fvalue, xp[], sc[], snwtn[], gc[], gp[], *h;
	Mint            *ldh;
	Mfloat          wk1[], wk2[], wk3[];
#endif
{
#define H(I_,J_)	(h+(I_)*(*ldh)+(J_))
	Mint            fdiff, mxtake;
	Mint            _l0, i, icode, ihess, iter, nfcn, ngrad, nhess;
	Mfloat          delta, eps, epsfcn, fc, fdigit, fp, rnwtnl, scgrad,
	                stepmx, valmax;
        Mint            IMSLFALSE = 0, IMSLTRUE = 1;


	imsl_e1psh("U3INF ");
	/*
	 * CHECK THE VALIDITY OF THE USER SPECIFIED PARAMETERS
	 */
        _l0 = IMSLFALSE;
	l_u5inf(n, xc, xscale, fscale, (Mint *)&_l0, iparam, rparam);
	if (imsl_n1rty(1) == 5)
		goto L_9000;

	fdigit = iparam[1];
	lv_u16nf.mxiter = iparam[2];
	lv_u16nf.maxfcn = iparam[3];
	lv_u16nf.maxgrd = iparam[4];
	ihess = iparam[5];
	lv_u16nf.maxhes = iparam[6];

	lv_u16nf.gradtl = rparam[0];
	lv_u16nf.steptl = rparam[1];
	lv_u16nf.rftol = rparam[2];
	lv_u16nf.aftol = rparam[3];
	lv_u16nf.falstl = rparam[4];
	stepmx = rparam[5];
	delta = rparam[6];
	/*
	 * SET EPSFCN TO ESTIMATE OF RELATIVE NOISE IN FUNCTION
	 */
	eps = imsl_amach(4);
	epsfcn = imsl_f_max(eps, pow(F_TEN, -fdigit));
	/*
	 * INITIALIZE ITERATION, FUNCTION & GRADIENT EVALUATION COUNTER
	 */
	iter = 0;
	nfcn = 0;
	ngrad = 0;
	nhess = 0;
	icode = 0;
	fdiff = IMSLFALSE;
	/*
	 * EVALUATE THE FUNCTION & GRADIENT AT THE INITIAL POINT; ALSO SET
	 * FVALUE = FUNCTION VALUE AT INITIAL GUESS FOR THE CASE THAT THE
	 * INITIAL GUESS IS THE SOLUTION.
	 */
        imsl_e1usr("ON");
	fc = (*fcn) (*n, xc);
        imsl_e1usr("OFF");
	nfcn += 1;
	*fvalue = fc;
	l_fdgrd(fcn, n, xc, xscale, &fc, &epsfcn, gc);
	fdiff = IMSLTRUE;
	ngrad += 1;
	/*
	 * CHECK STOPPING CRITERIA AT THE INITIAL POINT
	 */
        _l0 = IMSLFALSE;
	l_u6inf(n, xc, sc, &fc, gc, xscale, fscale, &icode, &iter, &nfcn,
		   &ngrad, &nhess, (Mint *)&_l0, &mxtake);
	if (imsl_n1rcd(1) != 0 || icode == -999)
		goto L_60;
	/*
	 * GET THE (APPROXIMATE) HESSIAN AT THE INITAL POINT
	 */
	l_u9inf(n, &fc, fscale, xscale, &ihess, h, ldh);
	/* MAIN ITERATION LOOP */
L_10:
	;
	iter += 1;
	/*
	 * COMPUTE NEWTON STEP AND LENGTH OF SCALED NEWTON STEP
	 */
L_20:
	l_u10nf(n, h, ldh, gc, snwtn);
        _l0 = 1;
	l_u11nf(n, xscale, &_l0, snwtn, wk1);
	rnwtnl = imsl_snrm2(*n, wk1, 1);
	/* COMPUTE THE SEARCH DIRECTION */
	l_u17nf(fcn, &fdiff, n, xc, &fc, gc, snwtn, xscale, &stepmx, &lv_u16nf.steptl,
	 &icode, xp, &fp, gp, sc, &mxtake, &rnwtnl, &epsfcn, &nfcn, &ngrad);
	/*
	 * USE CENTRAL DIFFERENCE IF NECESSARY
	 */
	if ((icode == 1) && (fdiff)) {
		l_cdgrd(fcn, n, xc, xscale, &epsfcn, gc);
		fdiff = IMSLFALSE;
		valmax = F_ZERO;
		for (i = 1; i <= *n; i++) {
			scgrad = fabs(gc[i - 1]) * imsl_f_max(fabs(xc[i - 1]), F_ONE /
			     xscale[i - 1]) / imsl_f_max(fabs(fc), *fscale);
			valmax = imsl_f_max(scgrad, valmax);
		}
		if (valmax <= lv_u16nf.gradtl) {
			*fvalue = fc;
			goto L_60;
		}
		goto L_20;
	}
	/*
	 * EVALUATE THE GRADIENT AT NEW POINT
	 */
	if (icode == 1)
		goto L_50;
	if (fdiff) {
		valmax = F_ZERO;
		for (i = 1; i <= *n; i++) {
			scgrad = fabs(gp[i - 1]) * imsl_f_max(fabs(xp[i - 1]), F_ONE /
			     xscale[i - 1]) / imsl_f_max(fabs(fp), *fscale);
			valmax = imsl_f_max(scgrad, valmax);
		}

		if (valmax <= 0.1e0) {
			l_cdgrd(fcn, n, xp, xscale, &epsfcn, gp);
			ngrad += 1;
			fdiff = IMSLFALSE;
		}
	}
	/* CHECK STOPPING CRITERIA AT NEW POINT */
L_50:
        _l0 = IMSLFALSE;
	l_u6inf(n, xp, sc, &fp, gp, xscale, fscale, &icode, &iter, &nfcn,
		   &ngrad, &nhess, (Mint *)&_l0, &mxtake);
	if (imsl_n1rcd(1) == 0 && icode != -999) {
		/*
		 * UPDATE THE HESSIAN APPROXIMATION, XC, GC, FC; NEXT
		 * ITERATION
		 */
                _l0 = IMSLFALSE; 
		l_u8inf(n, sc, gc, gp, &epsfcn, (Mint *)&_l0, h, ldh, wk1,
			   wk2, wk3);
		scopy(*n, xp, 1, xc, 1);
		scopy(*n, gp, 1, gc, 1);
		fc = fp;
		goto L_10;
	}
	/*
	 * OTHERWISE THE STOPPING CRITERIA IS SATISFIED; RETURN
	 */
	scopy(*n, xp, 1, xc, 1);
	scopy(*n, gp, 1, gc, 1);
	*fvalue = fp;
L_60:
	iparam[2] = iter;
	iparam[3] = nfcn;
	iparam[4] = ngrad;

L_9000:
	imsl_e1pop("U3INF ");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  U3ING/DU3ING (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Driver for unconstrained minimization solver using finite
                difference or analytic gradient.

    Usage:      CALL U3ING (FCN, GRAD, N, XC, XSCALE, FSCALE, IPARAM,
                            RPARAM, FVALUE, XP, SC, SNWTN, GC, GP, H,
                            LDH, WK1, WK2, WK3)

    Arguments:
       FCN    - User-supplied SUBROUTINE to evaluate the function to be
                minimized.  The usage is
                CALL FCN (N, X, F), where
                N      - Length of X.  (Input)
                X      - The point at which the function is evaluated.
                         (Input)
                         X should not be changed by FCN.
                F      - The computed function value at the point X.
                         (Output)
                FCN must be declared EXTERNAL in the calling program.
       GRAD   - User-supplied SUBROUTINE to compute the gradient at the
                point X.  The usage is
                CALL GRAD (N, X, G), where
                N      - Length of X and G.  (Input)
                X      - The point at which the gradient is evaluated.
                         (Input)
                         X should not be changed by GRAD.
                G      - The gradient evaluated at the point X.
                         (Output)
                GRAD must be declared EXTERNAL in the calling program.
       N      - Dimension of the problem.  (Input)
       XC     - Real vector of length N containing initial guess on input
                and approximate solution on output.  (Input / Output)
       XSCALE - Real vector of length N containing the diagonal scaling
                matrix for the variables.  (Input)
       FSCALE - Real scalar containing the function scaling.  (Input)
       IPARAM - Integer parameters vector of length 7.  (Input)
                IPARAM(1) = initialization flag
                IPARAM(2) = number of good digits in the function
                IPARAM(3) = maximum number of iterations.
                IPARAM(4) = maximum number of function evaluations
                IPARAM(5) = maximum number of gradient evaluations
                IPARAM(6) = initialization of Hessian parameter
                IPARAM(7) = maximum number of Hessian evaluations
       RPARAM - Real parameters vector of length 7.  (Input)
                RPARAM(1) = scaled gradient tolerance
                RPARAM(2) = scaled step tolerance
                RPARAM(3) = relative function tolerance
                RPARAM(4) = absolute function tolerance
                RPARAM(5) = false convergence tolerance
                RPARAM(6) = maximum allowable step size
                RPARAM(7) = size of initial trust region radius
       FVALUE - Real scalar containing the value of the function at the
                solution.  (Output)
       XP     - Real vector of length N containing the updated point.
                   (Output)
       SC     - Real vector of length N containing the last step taken.
                   (Output)
       SNWTN  - Real vector of length N containing the last Newton step.
                   (Output)
       GC     - Real vector of length N containing an estimate of the
                gradient at the current point.  (Output)
       GP     - Real vector of length N containing an estimate of the
                gradient at the updated point.  (Output)
       H      - Real N by N matrix containing an estimate of the Hessian
                at the approximate solution.  (Output)
       LDH    - Leading dimension of H exactly as specified in the
                dimension statement of the calling program.  (Input)
       WK1    - Real work vector of length N.  (Output)
       WK2    - Real work vector of length N.  (Output)
       WK3    - Real work vector of length N.  (Output)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_u3ing(Mfloat (*fcn)(Mint, Mfloat*), void (*grad)(Mint, Mfloat*,
                    Mfloat*), Mint *n, Mfloat xc[], Mfloat xscale[],
                    Mfloat *fscale, Mint iparam[], Mfloat rparam[],
                    Mfloat *fvalue, Mfloat xp[], Mfloat sc[], Mfloat snwtn[],
                    Mfloat gc[], Mfloat gp[], Mfloat *h, Mint *ldh,
                    Mfloat wk1[], Mfloat wk2[], Mfloat wk3[])
#else
static void l_u3ing(Mfloat (*fcn)(Mint, Mfloat[]), void (*grad)(Mint, Mfloat[],
                    Mfloat[]), Mint *n, Mfloat xc[], Mfloat xscale[],
                    Mfloat *fscale, Mint iparam[], Mfloat rparam[],
                    Mfloat *fvalue, Mfloat xp[], Mfloat sc[], Mfloat snwtn[],
                    Mfloat gc[], Mfloat gp[], Mfloat *h, Mint *ldh,
                    Mfloat wk1[], Mfloat wk2[], Mfloat wk3[])
#endif
#else
static void l_u3ing(fcn, grad, n, xc, xscale, fscale, iparam, rparam, fvalue,
                    xp, sc, snwtn, gc, gp, h, ldh, wk1, wk2, wk3)
	Mfloat          (*fcn) ();
	Mvoid		(*grad) ();
	Mint            *n;
	Mfloat          xc[], xscale[], *fscale;
	Mint            iparam[];
	Mfloat          rparam[], *fvalue, xp[], sc[], snwtn[], gc[], gp[], *h;
	Mint            *ldh;
	Mfloat          wk1[], wk2[], wk3[];
#endif
{
#define H(I_,J_)	(h+(I_)*(*ldh)+(J_))
	Mint            fdiff, mxtake;
	Mint            _l0, icode, ihess, iter, nfcn, ngrad, nhess;
	Mfloat          eps, epsfcn, fc, fdigit, fp, rnwtnl, stepmx;
        Mint            IMSLFALSE = 0, IMSLTRUE = 1;


	imsl_e1psh("U3ING ");
	/*
	 * CHECK THE VALIDITY OF THE USER SPECIFIED PARAMETERS
	 */
        _l0 = IMSLFALSE;
	l_u5inf(n, xc, xscale, fscale, (Mint *)&_l0, iparam, rparam);
	if (imsl_n1rty(1) == 5)
		goto L_9000;

	fdigit = iparam[1];
	lv_u16nf.mxiter = iparam[2];
	lv_u16nf.maxfcn = iparam[3];
	lv_u16nf.maxgrd = iparam[4];
	ihess = iparam[5];
	lv_u16nf.maxhes = iparam[6];

	lv_u16nf.gradtl = rparam[0];
	lv_u16nf.steptl = rparam[1];
	lv_u16nf.rftol = rparam[2];
	lv_u16nf.aftol = rparam[3];
	lv_u16nf.falstl = rparam[4];
	stepmx = rparam[5];
	/*
	 * SET EPSFCN TO ESTIMATE OF RELATIVE NOISE IN FUNCTION
	 */
	eps = imsl_amach(4);
	epsfcn = imsl_f_max(eps, pow(F_TEN, -fdigit));
	/*
	 * INITIALIZE ITERATION, FUNCTION & GRADIENT EVALUATION COUNTER
	 */
	iter = 0;
	nfcn = 0;
	ngrad = 0;
	nhess = 0;
	/*
	 * EVALUATE THE FUNCTION & GRADIENT AT THE INITIAL POINT; ALSO SET
	 * FVALUE = FUNCTION VALUE AT INITIAL GUESS FOR THE CASE THAT THE
	 * INITIAL GUESS IS THE SOLUTION.
	 */
	imsl_e1usr("ON");
	fc = (*fcn) (*n, xc);
	imsl_e1usr("OFF");
	nfcn += 1;
	*fvalue = fc;
	imsl_e1usr("ON");
	(*grad) (*n, xc, gc);
	imsl_e1usr("OFF");
	ngrad += 1;
	icode = 0;
	fdiff = IMSLFALSE;
	/*
	 * CHECK STOPPING CRITERIA AT THE INITIAL POINT
	 */
        _l0 = IMSLFALSE;
	l_u6inf(n, xc, sc, &fc, gc, xscale, fscale, &icode, &iter, &nfcn,
		   &ngrad, &nhess, (Mint *)&_l0, &mxtake);
	if (imsl_n1rcd(1) != 0 || icode == -999)
		goto L_40;
	/*
	 * GET THE (APPROXIMATE) HESSIAN AT THE INITAL POINT
	 */
	l_u9inf(n, &fc, fscale, xscale, &ihess, h, ldh);
	/* MAIN ITERATION LOOP */
L_10:
	;
	iter += 1;
	/*
	 * COMPUTE NEWTON STEP AND LENGTH OF SCALED NEWTON STEP
	 */

	l_u10nf(n, h, ldh, gc, snwtn);
        _l0 = 1;
	l_u11nf(n, xscale, &_l0, snwtn, wk1);
	rnwtnl = imsl_snrm2(*n, wk1, 1);
	/* LINE SEARCH */
	l_u17ng(fcn, grad, &fdiff, n, xc, &fc, gc, snwtn, xscale, &stepmx,
	  &lv_u16nf.steptl, &icode, xp, &fp, gp, sc, &mxtake, &rnwtnl, &epsfcn,
		   &nfcn, &ngrad);
	/* CHECK STOPPING CRITERIA AT NEW POINT */

        _l0 = IMSLFALSE;
	l_u6inf(n, xp, sc, &fp, gp, xscale, fscale, &icode, &iter, &nfcn,
		   &ngrad, &nhess, (Mint *)&_l0, &mxtake);
	if (imsl_n1rcd(1) == 0 && icode != -999) {
		/*
		 * UPDATE THE HESSIAN APPROXIMATION, XC, GC, FC; NEXT
		 * ITERATION
		 */
                _l0 = IMSLTRUE;
		l_u8inf(n, sc, gc, gp, &epsfcn, (Mint *)&_l0, h, ldh, wk1,
			   wk2, wk3);
		scopy(*n, xp, 1, xc, 1);
		scopy(*n, gp, 1, gc, 1);
		fc = fp;
		goto L_10;
	}
	/*
	 * OTHERWISE THE STOPPING CRITERIA IS SATISFIED; RETURN
	 */
	scopy(*n, xp, 1, xc, 1);
	scopy(*n, gp, 1, gc, 1);
	*fvalue = fp;
L_40:
	iparam[2] = iter;
	iparam[3] = nfcn;
	iparam[4] = ngrad;

L_9000:
	imsl_e1pop("U3ING ");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  U5INF/DU5INF  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Check validity of input to unconstrained minimization.

    Usage:      CALL U5INF (N, XC, XSCALE, FSCALE, USRHES, IPARAM,
                            RPARAM)

    Arguments:
       N      - Size of the problem.  (Input)
       XC     - Vector of length N containing the initial point.  (Input)
       XSCALE - Vector of length N containing the diagonal scaling matrix
                for the variables.  (Input/Output)
       FSCALE - Estimate of the scale of the objective function.
                (Input/Output)
       USRHES - Logical variable.  (Input)
                USRHES = .TRUE. if analytic Hessian or finite difference
                         Hessian is used.
                USRHES = .FALSE. otherwise.
       IPARAM - Integer parameters vector of length 6.  (Input/Output)
                See UMINF or UMIDH for details.
       RPARAM - Real parameters vector of length 7.  (Input/Output)
                See UMINF or UMIDH for details.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_u5inf(Mint *n, Mfloat xc[], Mfloat xscale[], Mfloat *fscale,
                    Mint *usrhes, Mint iparam[], Mfloat rparam[])
#else
static void l_u5inf(n, xc, xscale, fscale, usrhes, iparam, rparam)
	Mint            *n;
	Mfloat          xc[], xscale[], *fscale;
	Mint            *usrhes;
	Mint            iparam[];
	Mfloat          rparam[];
#endif
{
	Mint            _l0, i;
	Mfloat          machep, temp, temp1, temp2;


	imsl_e1psh("U5INF ");
	machep = imsl_amach(4);
	/*
	 * FDIGIT  good digits in F 
         * MXITER  Min number of iterations
         * MAXFCN  Max fcn. evaluations 
         * MAXGRD  Max grad. evaluations
         * IHESS   Hessian initial. parameter 
         * MAXHES  Max Hess. evaluations 
         * GRADTL  Scaled gradient tolerance 
         * STEPTL  Scaled step tolerance
         * RFTOL   Relative function tolerance 
         * AFTOL   Absolute function tolerance
         * FALSTL  False convergence tolerance 
         * STEPMX  Maximum allowable step size
	 * DELTA   Size of initial trust region
	 */
	if (*n <= 0) {
		imsl_e1sti(1, *n);
                _l0 = 0;
		l_u19nf(&_l0);
                /* Print error message:  (5, 1, 'The size of the problem must
                 * be positive while N = %(I1) is given.') */
		goto L_9000;

	} else if (*n == 1) {
                _l0 = 1;
		l_u19nf(&_l0);
                /* Print error message:  (6, 1, 'This routine may be inefficient                 *                        for a problem of size N = 1.') */
	}

	/* CHECK VARIABLE SCALING MATRIX */
	for (i = 1; i <= *n; i++) {
		if (xscale[i - 1] <= F_ZERO)
			goto L_20;
	}
	goto L_30;
L_20:
        _l0 = 2;
	l_u19nf(&_l0);
        /* Print error message:  (6, 2, 'The diagonal scaling matrix for the
         *                        variables must be positive while some of the
         *                        entries are less than or equal to zero.  The
         *                        algorithm will use the identity scaling
         *                        matrix for XSCALE.') */
	sset(*n, F_ONE, xscale, 1);

	/* CHECK FUNCTION SCALING */
L_30:
	if (*fscale <= F_ZERO) {
		imsl_e1str(1, *fscale);
/*		(6, 3, "The estimate of the scale of the objective function must be positive while FSCALE = %(r1) is given.  The algorithm will use FSCALE = 1.0.");
*/
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_FSCALE_VALUE_TOO_SMALL);
		*fscale = F_ONE;
	}
	/* CHECK ACCURACY OF THE PROBLEM */
	if (iparam[1] <= 0) {
		imsl_e1sti(1, iparam[1]);
                _l0 = 4;
		l_u19nf(&_l0);
		/* Print error message:  (6, 4, 'The estimate of the number of
		 *                        good digits in the function must be 
                 *                        positive while NDIGIT = %(I1) is 
                 *                        given.  The algorithm will assume that
                 *                        the function is accurate to the 
                 *                        precision of the arithmetic.')
		 */
		iparam[1] = imsl_i_machine(7);
	}
	/* CHECK MAXIMUM NUMBER OF ITERATIONS */
	if (iparam[2] <= 0) {
		imsl_e1sti(1, iparam[2]);
                _l0 = 5;
		l_u19nf(&_l0);
		/* Print error message:  (6, 5, 'The maximum number of 
                 *                        iterations must be positive while 
                 *                        MXITER = %(I1) is given.  The 
                 *                        algorithm will use MXITER = 100.')
		 */
		iparam[2] = 100;
	}
	/* CHECK MAXIMUM FUNCTION EVALUATIONS */
	if (iparam[3] <= 0) {
		imsl_e1sti(1, iparam[3]);
                _l0 = 6;
		l_u19nf(&_l0);
		/* Print error message:  (6, 6, 'The maximum number of function
                 *                        evaluations must be positive while 
                 *                        MAXFCN = %(I1) is given.  The 
                 *                        algorithm will use MAXFCN = 400.')
		 */
		iparam[3] = 400;
	}
	/* CHECK MAXIMUM GRADIENT EVALUATIONS */
	if (iparam[4] <= 0) {
		imsl_e1sti(1, iparam[4]);
                _l0 = 7;
		l_u19nf(&_l0);
		/* Print error message:  (6, 7, 'The maximum number of gradient
                 *                        evaluations must be positive while 
                 *                        MAXGRD = %(I1) is given.  The 
                 *                        algorithm will use MAXGRD = 400.')
		 */
		iparam[4] = 400;
	}
	/*
	 * CHECK MAXIMUM HESSIAN EVALUATIONS IF A NEWTON METHOD IS USED
	 */
	if (*usrhes) {
		if (iparam[6] <= 0) {
			imsl_e1sti(1, iparam[6]);
                        _l0 = 8;
			l_u19nf(&_l0);
			/* Print error message:  (6, 8, 'The maximum number of 
                         *                        Hessian evaluations must be
                         *                        positive while MAXHES = %(I1)
                         *                        is given.  The algorithm will
                         *                        use MAXHES = 100.')
			 */
			iparam[6] = 100;
		}
	}
	temp1 = pow(machep, F_TWO / F_THREE);
	temp2 = pow(machep, F_TWO / F_THREE);
	/* CHECK THE GRADIENT TOLERANCE */
	if (rparam[0] < F_ZERO) {
		imsl_e1str(1, rparam[0]);
		imsl_e1str(2, temp1);
/*		(6, 9, "The gradient tolerance must be nonnegative while   */
/*                      GRAD_TOL = %(r1) is given.  The algorithm will use */
/*                      GRAD_TOL = %(r2)."); */
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEED_POSITIVE_GRADTL);
		rparam[0] = temp1;
	}
	/* CHECK THE STEP TOLERANCE */
	if (rparam[1] < F_ZERO) {
		imsl_e1str(1, rparam[1]);
		imsl_e1str(2, temp2);
/*		(6, 10, "The step tolerance must be nonnegative while */
/*                       STEP_TOL = %(r1) is given.  The algorithm will use */
/*                       STEP_TOL = %(r2)."); */
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEGATIVE_STEP_TOL);
		rparam[1] = temp2;
	}
	/* CHECK RELATIVE FUNCTION TOLERANCE */
	if (rparam[2] < F_ZERO) {
		temp = imsl_f_max(1.0e-10, temp2);
		imsl_e1str(1, rparam[2]);
		imsl_e1str(2, temp);
/*		(6, 11, "The relative function tolerance must be nonnegative */
/*                       while RFCN_TOL = %(r1) is given.  The algorithm will */
/*                       use RFCN_TOL = %(r2)."); */
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEGATIVE_REL_FCN_TOL);
		rparam[2] = temp;
	}
	/* CHECK ABSOLUTE FUNCTION TOLERANCE */
	if (rparam[3] < F_ZERO) {
		temp = imsl_f_max(1.0e-20, machep * machep);
		imsl_e1str(1, rparam[3]);
		imsl_e1str(2, temp);
/*		(6, 12, "The absolute function tolerance must be nonnegative */
/*                       while AFCN_TOL = %(r1) is given.  The algorithm will */
/*                       use AFCN_TOL = %(r2).");                            */
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEGATIVE_ABS_FCN_TOL);
		rparam[3] = temp;
	}
	/* CHECK FALSE CONVERGENCE TOLERANCE */
	if (rparam[4] < F_ZERO) {
		temp = 1.0e2 * machep;
		imsl_e1str(1, rparam[4]);
		imsl_e1str(2, temp);
/*		(6, 13, "The false convergence tolerance must be nonnegative */
/*              while FALSTOL = %(r1) is given.  The algorithm will use */
/*              FALSTOL = %(r2)."); */
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEGATIVE_FALSE_CONV_TOL);
		rparam[4] = temp;
	}
	/* CHECK MAXIMUM ALLOWED STEP SIZE */
	if (rparam[5] <= F_ZERO) {
	    temp1 = F_ZERO;
	    for (i = 1; i <= *n; i++) {
		temp1 += imsl_fi_power(xscale[i - 1] * xc[i - 1], 2);
	    }
	    temp1 = sqrt(temp1);
	    temp2 = imsl_snrm2(*n, xscale, 1);
	    temp = 1.0e3 * imsl_f_max(temp1, temp2);
	    if (iparam[0] != 0 && rparam[5] != -9999.0e0) {
		imsl_e1str(1, rparam[5]);
		imsl_e1str(2, temp);
/*		(6, 14, "The maximum allowable scaled step length must be */
/*                       positive while MAX_STEP = %(r1) is given.  The   */
/*                       algorithm will use MAX_STEP =  %(r2).");         */
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEED_NONNEGATIVE_STEPMX);
	    }
	    rparam[5] = temp;
	}
	/* CHECK INITIAL TRUST REGION RADIUS */
	if (rparam[6] <= F_ZERO) {
	    if (iparam[0] != 0 && rparam[6] != -9999.0e0) {
		imsl_e1str(1, rparam[6]);
/*		(6, 15, "The initial trust region radius must be  */
/*                       positive while DELTA = %(r1) is given.   */
/*                       The algorithm will use the length of the */
/*                       initial scaled Cauchy step for DELTA."); */
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEED_NONNEGATIVE_DELTA);
	    }
	    rparam[6] = -F_ONE;
	}
L_9000:
	imsl_e1pop("U5INF ");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  U6INF/DU6INF  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Stopping conditions for unconstrained minimization.

    Usage:      CALL U6INF (N, XP, SC, FP, GP, XSCALE, FSCALE, ICODE,
                            ITER, NFCN, NGRAD, NHESS, USRHES, MXTAKE)

    Arguments:
       N      - Dimension of the problem.  (Input)
       XP     - Vector of length N containing the new iterate.
                (Input)
       SC     - Vector of length N containing step taken.  (Input)
       FP     - Scalar containing the function value at XP.  (Input)
       GP     - Vector of length N containing the gradient at XP.
                (Input)
       XSCALE - Vector of length N containing the diagonal scaling
                matrix for the variables.  (Input)
       FSCALE - Scalar containing the function scaling.  (Input)
       ICODE  - Return code from the global strategy algorithm.  (Input)
       ITER   - Number of iterations.  (Input)
       NFCN   - Number of function evaluations.  (Input)
       NGRAD  - Number of gradient evaluations.  (Input)
       NHESS  - Number of Hessian evaluations.   (Input)
       USRHES - Logical variable.  (Input)
                USRHES = .TRUE. if Newton's method is used.
                USRHES = .FALSE. otherwise.
       MXTAKE - Logical variable indicating a step of maximum length was
                taken.  (Input)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_u6inf(Mint *n, Mfloat xp[], Mfloat sc[], Mfloat *fp,
                    Mfloat gp[], Mfloat xscale[], Mfloat *fscale,
                    Mint *icode, Mint *iter, Mint *nfcn, Mint *ngrad,
                    Mint *nhess, Mint *usrhes, Mint *mxtake)
#else
static void l_u6inf(n, xp, sc, fp, gp, xscale, fscale, icode, iter, nfcn,
                    ngrad, nhess, usrhes, mxtake)
	Mint            *n;
	Mfloat           xp[], sc[], *fp, gp[], xscale[], *fscale;
	Mint            *icode, *iter, *nfcn, *ngrad, *nhess;
	Mint           *usrhes, *mxtake;
#endif
{
	Mint            _l0, i;
	static Mint     nmaxs;
	Mfloat          scgrad, scstep, valmax;

	imsl_e1psh("U6INF ");

	/* TEST OF NORM OF SCALED GRADIENT */
	valmax = F_ZERO;
	for (i = 1; i <= *n; i++) {
		scgrad = fabs(gp[i - 1]) * imsl_f_max(fabs(xp[i - 1]), F_ONE /
			    xscale[i - 1]) / imsl_f_max(fabs(*fp), *fscale);
		valmax = imsl_f_max(scgrad, valmax);
	}
	if (valmax <= lv_u16nf.gradtl) {
		*icode = -999;
		goto L_9000;
	}
	/*
	 * IF FIRST ITER., INITIALIZE COUNTER FOR MX. STEP TAKEN AND RETURN
	 */
	if (*iter == 0) {
		nmaxs = 0;
		goto L_9000;
	}
	/* CHECK LAST GLOBAL STEP */
	if (*icode == 1) {
		imsl_e1str(1, lv_u16nf.steptl);

/*		(3, 8, "The last global step failed to locate a lower point than the current X value.  The current X may be an approximate local minimizer and no more accuracy is possible or the step tolerance may be too large where STEP_TOL = %(r1) is given.");
*/
               imsl_ermes(IMSL_WARNING, IMSL_NO_FURTHER_PROGRESS);
		goto L_9000;
	}
	/* TEST NORM OF SCALED STEP */
	valmax = F_ZERO;
	for (i = 1; i <= *n; i++) {
		scstep = fabs(sc[i - 1]) / imsl_f_max(fabs(xp[i - 1]), F_ONE
                          / xscale[i - 1]);
		valmax = imsl_f_max(scstep, valmax);
	}
	if (valmax <= lv_u16nf.steptl) {
		*icode = -999;
                imsl_ermes(IMSL_NOTE, IMSL_STEP_TOLERANCE);
		goto L_9000;
	}
	/*
	 * CHECK RELATIVE FUNCTION CONVERGENCE TOLERANCE
	 */
	if (*icode == 2) {
		imsl_e1str(1, lv_u16nf.rftol);

/*		(3, 1, "RELATIVE FUNCTION CONVERGENCE - Both the actual and predicted relative reductions in the function are less than or equal to the relative function convergence tolerance FCN_TOL = %(r1).");
*/
                imsl_ermes(IMSL_WARNING, IMSL_REL_FCN_TOLERANCE);
		goto L_9000;
	}
	/* CHECK FALSE CONVERGENCE TOLERANCE */
	if (*icode == 3) {

/*		(4, 2, "FALSE CONVERGENCE - The iterates appear to be converging to a noncritical point.  Possibly incorrect gradient information is used, or the function is discontinuous, or the other stopping tolerances are too tight.");
*/
                imsl_ermes(IMSL_FATAL, IMSL_FALSE_CONVERGENCE);
		goto L_9000;
	}
	/*
	 * CHECK ITERATION, FUNCTION, GRADIENT & HESSIAN EVALUATIONS LIMIT
	 */
	if (*iter >= lv_u16nf.mxiter) {
		/* Print error message */
                _l0 = 3;
		l_u18nf(&_l0);
		/*
		 * C         CALL E1MES (4, 3, ' Maximum number of iterations
		 * '// C     &               'exceeded.')
		 */
	} else if (*nfcn >= lv_u16nf.maxfcn) {
		/* Print error message */
                _l0 = 4;
		l_u18nf(&_l0);
		/*
		 * C         CALL E1MES (4, 4, ' Maximum number of function
		 * evaluations'// C     &               ' exceeded.')
		 */
	} else if (*ngrad >= lv_u16nf.maxgrd) {
		/* Print error message */
                _l0 = 5;
		l_u18nf(&_l0);
		/*
		 * C         CALL E1MES (4, 5, ' Maximum number of gradient
		 * evaluations'// C     &               ' exceeded.')
		 */
	} else if (*usrhes && (*nhess >= lv_u16nf.maxhes)) {
		/* Print error message */
                _l0 = 7;
		l_u18nf(&_l0);
		/*
		 * C         CALL E1MES (4, 7, 'Maximum number of Hessian
		 * evaluations '// C     &               'exceeded.')
		 */
	} else if (*mxtake) {
		nmaxs += 1;
		if (nmaxs == 5) {
			/* Print error message */
                        _l0 = 6;
			l_u18nf(&_l0);
			/*
			 * C            CALL E1MES (4, 6, ' Five consecutive
			 * steps of '// C     &                  'length
			 * STEPMX have been taken; either the '// C     &
			 * 'function is unbounded below, or has a '// C     &
			 * 'finite asymptote in some direction or the '// C
			 * &                  'maximum allowable step size
			 * STEPMX is too '// C     &
			 * 'small.')
			 */
		}
	}
L_9000:
	imsl_e1pop("U6INF ");
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  U8INF/DU8INF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Compute the Dennis-Schnabel BFGS factor update.

    Usage:      CALL U8INF (N, SC, GC, GP, EPSFCN, USRDER, A, LDA, Y,
                            T, U)

    Arguments:
       N      - The size of the problem. (Input)
       SC     - Real vector of length N containing the current step.
                   (Input)
       GC     - Real vector of length N containing the gradient at the
                current point.  (Input)
       GP     - Real vector of length N containing the gradient at the
                new point.  (Input)
       EPSFCN - An estimated bound on the relative noise in the function
                value. (Input)
       USRDER - Logical variable. (Input)
                USRDER = .TRUE.  means that analytic gradient is used
                USRDER = .FALSE. means that finite differences gradient
                                 is used.
       A      - Real N by N matrix.  (Input/Output)
                On input the Cholesky decomposition of the current appro-
                   ximate Hessian is in the lower part and diagonal of A.
                On output the updated Cholesky decomposition is in lower
                   part and diagonal of A.
       LDA    - Row dimension of A exactly as specified in the dimension
                statement of the calling program.  (Input)
       Y      - Real work vector of length N.
       T      - Real work vector of length N.
       U      - Real work vector of length N.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_u8inf(Mint *n, Mfloat sc[], Mfloat gc[], Mfloat gp[],
                    Mfloat *epsfcn, Mint *usrder, Mfloat *a, Mint *lda,
                    Mfloat y[], Mfloat t[], Mfloat u[])
#else
static void l_u8inf(n, sc, gc, gp, epsfcn, usrder, a, lda, y, t, u)
	Mint            *n;
	Mfloat           sc[], gc[], gp[], *epsfcn;
	Mint           *usrder;
	Mfloat          *a;
	Mint            *lda;
	Mfloat           y[], t[], u[];
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
	Mint            skpupd;
	Mint             _l0, _l1, i;
	Mfloat           alpha, gmax, machep, *qdum, snorm, temp1,
	                temp3, tnorm, tol, ynorm;
        Mint            IMSLFALSE = 0, IMSLTRUE = 1;


        qdum = (Mfloat *) imsl_malloc (2* *n * sizeof(*qdum));
	machep = imsl_amach(4);
	/* COMPUTE SECANT CONDITION */
	for (i = 1; i <= *n; i++) {
		y[i - 1] = gp[i - 1] - gc[i - 1];
	}
	temp1 = imsl_sdot(*n, y, 1, sc, 1);
	snorm = imsl_snrm2(*n, sc, 1);
	ynorm = imsl_snrm2(*n, y, 1);
	/* DETERMINE TO SKIP UPDATE OR NOT */
	if (temp1 >= (sqrt(machep) * snorm * ynorm)) {
		/* COMPUTE (L**T) * SC AND STORE IT IN T */
		for (i = 1; i <= *n; i++) {
			t[i - 1] = imsl_sdot(*n - i + 1, A(i - 1, i - 1), 1, &sc[i - 1],
					     1);
		}
		tnorm = imsl_snrm2(*n, t, 1);
		alpha = sqrt(temp1) / tnorm;
		if (*usrder) {
			tol = *epsfcn;
		} else {
			tol = sqrt(*epsfcn);
		}
		/* DETERMINE TO SKIP UPDATE OR NOT */
		skpupd = IMSLTRUE;
		for (i = 1; i <= *n; i++) {
			temp3 = imsl_sdot(i, A(0, i - 1), *lda, t, 1);
			gmax = imsl_f_max(fabs(gc[i - 1]), fabs(gp[i - 1]));
			if ((fabs(y[i - 1] - temp3)) >= (tol * gmax))
				skpupd = IMSLFALSE;
			u[i - 1] = y[i - 1] - alpha * temp3;
		}
		/* PERFORM UPDATE STEP */
		if (!skpupd) {
			temp3 = F_ONE / (sqrt(temp1) * tnorm);
			sscal(*n, temp3, t, 1);
			/*
			 * COPY L(transpose) INTO UPPER TRIANGLE
			 */
			for (i = 2; i <= *n; i++) {
				scopy(i - 1, A(0, i - 1), *lda, A(i - 1, 0), 1);
				sset(i - 1, F_ZERO, A(0, i - 1), *lda);
			}
			/* QR FACTORIZATION OF L(t) + U*V(t) */
                        _l0 = 2;
                        _l1 = 1;
			l_u15nf(&_l0, n, t, u, qdum, &_l1, a, lda);
			/*
			 * COPY TRANSPOSE OF UPPER TRIANGLE OF A INTO L
			 */
			for (i = 2; i <= *n; i++) {
				scopy(i - 1, A(i - 1, 0), 1, A(0, i - 1), *lda);
			}
		}
	}
        imsl_free(qdum);
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  U9INF/DU9INF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Compute the Hessian at the initial point.

    Usage:      CALL U9INF (N, FX0, FSCALE, XSCALE, IHESS, H, LDH)

    Arguments:
       N      - Order of H.  (Input)
       FX0    - Function value at the initial guess.  (Input)
       FSCALE - Scaling for the function.  (Input)
       XSCALE - Real vector of length N containing the diagonal scaling
                matrix for the variables.  (Input)
       IHESS  - Hessian initialization parameter.  (Input)
                If IHESS = 0 the Hessian is initialized to the identity
                matrix, otherwise it is initialized to a diagonal matrix
                containing the Cholesky factor on the diagonal.
       H      - N by N matrix containing the initialized Hessian.
                (Output)
       LDH    - Row dimension of H exactly as specified in the dimension
                statement of the calling program.  (Input)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_u9inf(Mint *n, Mfloat *fx0, Mfloat *fscale, Mfloat xscale[],
                    Mint *ihess, Mfloat *h, Mint *ldh)
#else
static void l_u9inf(n, fx0, fscale, xscale, ihess, h, ldh)
	Mint            *n;
	Mfloat          *fx0, *fscale, xscale[];
	Mint            *ihess;
	Mfloat          *h;
	Mint            *ldh;
#endif
{
#define H(I_,J_)	(h+(I_)*(*ldh)+(J_))
	Mint             j;
	Mfloat           temp;


	temp = sqrt(imsl_f_max(fabs(*fx0), *fscale));
	for (j = 1; j <= *n; j++) {
		sset(*n, F_ZERO, H(j - 1, 0), 1);
		if (*ihess == 0) {
			*H(j - 1, j - 1) = F_ONE;
		} else {
			*H(j - 1, j - 1) = temp * xscale[j - 1];
		}
	}

	return;
}				/* end of function */
/* -----------------------------------------------------------------------
    IMSL Name:  FDGRD/DFDGRD  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Approximate the gradient using forward differences.

    Usage:      CALL FDGRD (FCN, N, XC, XSCALE, FC, EPSFCN, GC)

    Arguments:
       FCN    - User-supplied SUBROUTINE to evaluate the function to be
                minimized.  The usage is
                CALL FCN (N, X, F), where
                N      - Length of X.  (Input)
                X      - The point at which the function is evaluated.
                         (Input)
                         X should not be changed by FCN.
                F      - The computed function value at the point X.
                         (Output)
                FCN must be declared EXTERNAL in the calling program.
       N      - Dimension of the problem.  (Input)
       XC     - Vector of length N containing the point at which the
                gradient is to be estimated.  (Input)
       XSCALE - Vector of length N containing the diagonal scaling matrix
                for the variables.  (Input)
                In the absence of other information, set all entries
                to 1.0.
       FC     - Scalar containing the value of the function at XC.
                (Input)
       EPSFCN - Estimate of the relative noise in the function.  (Input)
                EPSFCN must be less than or equal to 0.1.  In the absence
                of other information, set EPSFCN to 0.0.
       GC     - Vector of length N containing the estimated gradient
                at XC.  (Output)

    Remark:
       This is Algorithm A5.6.3, Dennis and Schnabel, 1983, page 322.

    Keywords:   Forward difference; Gradient; Service routine

    GAMS:       G4f
 
    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_fdgrd(Mfloat (*fcn)(Mint, Mfloat*), Mint *n, Mfloat xc[],
                    Mfloat xscale[], Mfloat *fc, Mfloat *epsfcn, Mfloat gc[])
#else
static void l_fdgrd(Mfloat (*fcn)(Mint, Mfloat[]), Mint *n, Mfloat xc[],
                    Mfloat xscale[], Mfloat *fc, Mfloat *epsfcn, Mfloat gc[])
#endif
#else
static void l_fdgrd(fcn, n, xc, xscale, fc, epsfcn, gc)
        Mfloat          (*fcn) ();
        Mint            *n;
        Mfloat          xc[], xscale[], *fc, *epsfcn, gc[];
#endif
{
        Mint             j;
        Mfloat           eps, fnew, stepsz, xtempj;


        imsl_e1psh("FDGRD ");

        if (*epsfcn > 0.1e0 || *epsfcn < F_ZERO) {
             imsl_e1str(1, *epsfcn);
/*             (5, 2, "The estimate for the relative noise in the function must be between 0.0 and 0.1 while EPSFCN = %(r1) is given.");
*/
               imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_EPSFCN_VALUE);
        }

        if (imsl_n1rcd(0) == 0) {
             eps = sqrt(imsl_f_max(*epsfcn, imsl_amach(4)));
             for (j = 0; j < *n; j++) {
                   stepsz = (eps) * imsl_f_max(fabs(xc[j]), F_ONE / xscale[j]);
                   if (xc[j] < F_ZERO)
                       stepsz = -stepsz;
                   xtempj = xc[j];
                   xc[j] = xtempj + stepsz;
                   imsl_e1usr("ON");
                   fnew = (*fcn) (*n, xc);
                   imsl_e1usr("OFF");
                   xc[j] = xtempj;
                   gc[j] = (fnew - *fc) / stepsz;
             }
        }
        imsl_e1pop("FDGRD ");
        return;
}                               /* end of function */
/* -----------------------------------------------------------------------
    IMSL Name:  CDGRD/DCDGRD (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Approximate the gradient using central differences.

    Usage:      CALL CDGRD (FCN, N, XC, XSCALE, EPSFCN, GC)

    Arguments:
       FCN    - User-supplied SUBROUTINE to evaluate the function to be
                minimized.  The usage is
                CALL FCN (N, X, F), where
                N      - Length of X.  (Input)
                X      - The point at which the function is evaluated.
                         (Input)
                         X should not be changed by FCN.
                F      - The computed function value at the point X.
                         (Output)
                FCN must be declared EXTERNAL in the calling program.
       N      - Dimension of the problem.  (Input)
       XC     - Vector of length N containing the point at which the
                gradient is to be estimated.  (Input)
       XSCALE - Vector of length N containing the diagonal scaling
                matrix for the variables.  (Input)
                In the absence of other information, set all entries
                to 1.0.
       EPSFCN - Estimate for the relative noise in the function.  (Input)
                EPSFCN must be less than or equal to 0.1.  In the absence
                of other information, set EPSFCN to 0.0.
       GC     - Vector of length N containing the estimated gradient
                at XC.  (Output)

    Remark:
       This is Algorithm A5.6.4, Dennis and Schnabel, 1983, imsl_page 323.

    Keywords:   Central difference; Gradient; Service routine

    GAMS:       G4f

    Chapter:    MATH/LIBRARY Optimization
 
    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_cdgrd(Mfloat (*fcn)(Mint, Mfloat*), Mint *n, Mfloat xc[],
                    Mfloat xscale[], Mfloat *epsfcn, Mfloat gc[])
#else
static void l_cdgrd(Mfloat (*fcn)(Mint, Mfloat[]), Mint *n, Mfloat xc[],
                    Mfloat xscale[], Mfloat *epsfcn, Mfloat gc[])
#endif
#else
static void l_cdgrd(fcn, n, xc, xscale, epsfcn, gc)
        Mfloat          (*fcn) ();
        Mint            *n;
        Mfloat          xc[], xscale[], *epsfcn, gc[];
#endif
{
        Mint             j;
        Mfloat           eps, fnew1, fnew2, stepsz, xtempj;


        imsl_e1psh("CDGRD ");

        if (*epsfcn > 0.1e0 || *epsfcn < F_ZERO) {
             imsl_e1str(1, *epsfcn);
/*           (5, 2, "The estimate for the relative noise in the function must be between 0.0 and 0.1 while EPSFCN = %(r1) is given.");
*/
             imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_EPSFCN_VALUE);
        }

        if (imsl_n1rcd(0) == 0) {
             eps = sqrt(imsl_f_max(*epsfcn, imsl_amach(4)));
             for (j = 0; j < *n; j++) {
                   stepsz = (eps) * imsl_f_max(fabs(xc[j]), F_ONE / xscale[j]);
                   if (xc[j] < F_ZERO)
                       stepsz = -stepsz;
                   xtempj = xc[j];
                   xc[j] = xtempj + stepsz;
                   imsl_e1usr("ON");
                   fnew1 = (*fcn) (*n, xc);
                   imsl_e1usr("OFF");
                   xc[j] = xtempj - stepsz;
                   imsl_e1usr("ON");
                   fnew2 = (*fcn) (*n, xc);
                   imsl_e1usr("OFF");
                   xc[j] = xtempj;
                   gc[j] = (fnew1 - fnew2) / (F_TWO * stepsz);
             }
        }
        imsl_e1pop("CDGRD ");
        return;
}                               /* end of function */
