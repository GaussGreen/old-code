#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK   l_nonlin_least_squares(void (*fcn)(Mint, Mint, 
                                        Mfloat[], Mfloat[]), Mint m,
                                        Mint n, va_list argptr);
static void      l_u2lsf(void (*fcn)(Mint, Mint, Mfloat[], Mfloat[]), 
                         Mint m, Mint n, Mfloat xguess[], Mfloat 
                         xscale[], Mfloat fscale[], Mint iparam[],
                         Mfloat rparam[], Mfloat x[], Mfloat fvec[],
                         Mfloat fjac[], Mint ldfjac, Mfloat wk[],
                         Mint iwk[]); 
static void      l_u2lsj(void (*fcn)(Mint, Mint, Mfloat[], Mfloat[]), 
                         void (*jac)(Mint, Mint, Mfloat[], Mfloat[],
                         Mint), Mint m, Mint n, Mfloat xguess[],
                         Mfloat xscale[], Mfloat fscale[],
                         Mint iparam[], Mfloat rparam[], Mfloat x[],
                         Mfloat fvec[], Mfloat fjac[], Mint ldfjac,
                         Mfloat wk[], Mint iwk[], Mfloat[], Mint);
static void      l_u3lsf(void (*fcn)(Mint, Mint, Mfloat[], Mfloat[]),
                         Mint m, Mint n, Mfloat xguess[],
                         Mfloat xscale[], Mfloat fscale[],
                         Mint iparam[], Mfloat rparam[], Mfloat xc[],
                         Mfloat fc[], Mfloat fp[], Mfloat xp[],
                         Mfloat sc[], Mfloat gnstep[], Mfloat gc[],
                         Mfloat *fjc, Mint ldfjc, Mint ipvt[],
                         Mfloat wk1[], Mfloat wk2[], Mfloat wk3[],
                         Mfloat wk4[], Mfloat wk5[], Mfloat wk6[]);
static void      l_u3lsj(void (*fcn)(Mint, Mint, Mfloat[], Mfloat[]),
                         void (*jac)(Mint, Mint, Mfloat[], Mfloat[],
                         Mint), Mint m, Mint n, Mfloat xguess[],
                         Mfloat xscale[], Mfloat fscale[],
                         Mint iparam[], Mfloat rparam[], Mfloat xc[],
                         Mfloat fc[], Mfloat fp[], Mfloat xp[],
                         Mfloat sc[], Mfloat gnstep[], Mfloat gc[],
                         Mfloat *fjc, Mint ldfjc, Mint ipvt[],
                         Mfloat wk1[], Mfloat wk2[], Mfloat wk3[],
                         Mfloat wk4[], Mfloat wk5[], Mfloat wk6[],
                         Mfloat[], Mint);
static void      l_u5lsf(Mint m, Mint n, Mfloat xc[], Mfloat xscale[], 
                         Mfloat fscale[], Mint usrjac, Mint iparam[],
	                 Mfloat rparam[]);
static void      l_u6lsf(Mint m, Mint n, Mfloat xp[], Mfloat sc[],
                         Mfloat fp[], Mfloat fpnorm, Mfloat gp[],
                         Mfloat xscale[], Mfloat fscale[],
                         Mint *icode, Mint iter, Mint nfcn,
                         Mint njac, Mint usrjac, Mint mxtake);
static void      l_u7lsf(Mint n, Mfloat gc[], Mfloat *a, Mint lda,
                         Mint ipvt[], Mfloat xscale[], Mfloat qtf[],
                         Mfloat stepmx, Mfloat *delta, Mfloat *amu,
                         Mint *first, Mfloat sc[], Mfloat gnstep[],
                         Mint *gauss, Mfloat diag[], Mfloat wk1[],
                         Mfloat wk2[]);
static void      l_u8lsf(Mint iter, Mint n, Mfloat cjnorm[],
                         Mfloat xscale[]);
static void      l_u9lsf(void (*fcn) (Mint, Mint, Mfloat[], Mfloat[]), 
                         Mint m, Mint n, Mfloat xc[], Mfloat fcnorm, 
                         Mfloat gc[], Mfloat *a, Mint lda, Mint ipvt[], 
                         Mfloat sc[], Mfloat xscale[], Mint gauss,
                         Mfloat stepmx, Mfloat *delta, Mint *icode,
                         Mfloat xpprev[], Mfloat fpprev[], Mfloat xp[],
                         Mfloat fc[], Mfloat fp[], Mfloat *fpnorm,
                         Mint *mxtake, Mint *nfcn);
static void      l_u10sf(Mint m, Mint n, Mfloat *fjac, Mint ldfjac, 
                         Mfloat qraux[], Mfloat f[], Mfloat qtf[]);
static void      l_u11nf(Mint n, Mfloat y[], Mint k, Mfloat z[],
                         Mfloat x[]); 
static void      l_u11sf(Mint n, Mfloat *r, Mint ldr, Mint ipvt[], 
                         Mfloat diag[], Mfloat qtb[], Mfloat x[],
                         Mfloat sdiag[], Mfloat wa[]);
static void      l_u12sf(Mint n, Mfloat *h, Mint ldh, Mfloat y[], 
                         Mfloat snwtn[]);
static void      l_u13sf(Mint icode);
static void      l_f2jac(void (*fcn) (Mint, Mint, Mfloat[], Mfloat[]),
                         Mint m, Mint n, Mfloat xc[], Mfloat xscale[],
                         Mfloat fc[], Mfloat epsfcn, Mfloat *fjac,
                         Mint ldfjac, Mfloat work[]);
#else
static VA_LIST_HACK   l_nonlin_least_squares();
static void      l_u2lsf();
static void      l_u2lsj();
static void      l_u3lsf(); 
static void      l_u3lsj();
static void      l_u5lsf(); 
static void      l_u6lsf();  
static void      l_u7lsf();  
static void      l_u8lsf();  
static void      l_u9lsf();
static void      l_u10sf(); 
static void      l_u11sf(); 
static void      l_u11nf(); 
static void      l_u12sf(); 
static void      l_u13sf();  
static void      l_f2jac();
#endif

static Mfloat       *lv_value;
static Mint	     loop_counter_maximum = 100;

#ifdef ANSI
#if defined(COMPUTER_HP97C)
Mfloat *imsl_f_nonlin_least_squares(void (*fcn)(Mint, Mint, Mfloat*, 
                                    Mfloat*), Mint m, Mint n, ...)
#else
Mfloat *imsl_f_nonlin_least_squares(void (*fcn)(Mint, Mint, Mfloat[], 
                                    Mfloat[]), Mint m, Mint n, ...)
#endif
#else
Mfloat *imsl_f_nonlin_least_squares(fcn, m, n, va_alist)
    void (*fcn) ();
    Mint m;
    Mint n;
    va_dcl
#endif
{
    va_list argptr;

    VA_START(argptr, n);
    E1PSH ("imsl_f_nonlin_least_squares", "imsl_d_nonlin_least_squares");
    lv_value = NULL;
    IMSL_CALL (l_nonlin_least_squares(fcn, m, n, argptr));
    va_end(argptr);
    E1POP ("imsl_f_nonlin_least_squares", "imsl_d_nonlin_least_squares");
    return lv_value;
}


#ifdef ANSI
#if defined(COMPUTER_HP97C)
static VA_LIST_HACK l_nonlin_least_squares(void (*fcn) (Mint, Mint, Mfloat*, 
                                      Mfloat*), Mint m, Mint n, 
                                      va_list argptr)
#else
static VA_LIST_HACK l_nonlin_least_squares(void (*fcn) (Mint, Mint, Mfloat[], 
                                      Mfloat[]), Mint m, Mint n, 
                                      va_list argptr)
#endif
#else
static VA_LIST_HACK l_nonlin_least_squares(fcn, m, n, argptr)
    void          (*fcn) ();
    Mint          m;
    Mint          n;
    va_list       argptr;
#endif
{
    Mint          i, code, k, j, l, ll, il, ll1, il1, irank;
    Mint          arg_number      = 3;
    Mint          maxitn          = 100;
    Mint          maxfcn          = 400;
    Mint          maxjacobian     = 400;
    Mint          iscale          = 0;
    Mint          ndigit;
    Mint          fjac_col_dim    = n;
    Mint          jtj_inv_col_dim = n;
    Mfloat        *xguess         = NULL;
    Mfloat        *fvec           = NULL;
    Mfloat        **pfvec         = NULL;
    Mfloat        *xscale         = NULL;
    Mfloat        *jac            = NULL;
    Mfloat        *fjac           = NULL;
    Mfloat        **pfjac         = NULL;
    Mfloat        *jtj_inv        = NULL;
    Mfloat        **pjtj_inv      = NULL;
    Mfloat        *work_pjtj      = NULL;
    Mfloat        *work           = NULL;
    Mfloat        *fscale         = NULL;
    Mint          *iwork          = NULL;
    Mint          *rank           = NULL;
    Mfloat        max_step        = -999.0e0;
    Mfloat        trust_region    = -999.0e0;
    Mfloat        *tmp            = NULL;
    Mfloat        rparam[7];
    Mint          iparam[7];
    Mint          user_xguess     = 0;
    Mint          user_xscale     = 0;
    Mint          user_fscale     = 0;
    Mint          user_tol        = 0;
    Mint          user_pfjac      = 0;
    Mint          user_fjac       = 0;
    Mint          return_rank     = 0;
    Mint          return_user     = 0;
    Mint          user_pjtj       = 0;
    Mint          user_jtj        = 0;
    Mint          user_pfvec      = 0;
    Mint          user_fvec       = 0;
#ifdef ANSI
    void	  (*jacobian)(Mint, Mint, Mfloat[], Mfloat[], Mint) = NULL;
#else
    void          (*jacobian)()   =NULL;
#endif
    Mint          user_jacobian   = 0;
    Mfloat        *jtwork         = NULL;
    Mfloat        grad_tol, step_tol, rfcn_tol, afcn_tol;
    Mfloat        eps, sqrt_eps, sqr_eps, eps_onet, eps_twot;
    Mfloat        onet, twot, s2, tol, tolsq, temp, c, s;

    onet     = 1.0e0/3.0e0;
    twot     = 2.0e0/3.0e0;
    eps      = imsl_amach(4);
    sqrt_eps = sqrt(eps);
    sqr_eps  = pow(eps, 2.0);
    eps_onet = pow(eps, onet);
    eps_twot = pow(eps, twot);

#ifdef DOUBLE
    grad_tol = eps_onet; 
    rfcn_tol = imsl_f_max(1.0e-20, eps_twot);
    afcn_tol = imsl_f_max(1.0e-40, sqr_eps);
#else
    grad_tol = sqrt_eps;
    rfcn_tol = imsl_f_max(1.0e-10, eps_twot);
    afcn_tol = imsl_f_max(1.0e-20, sqr_eps);
#endif
    ndigit   = (Mint)(-log10(eps) + 0.1e0);
    step_tol = eps_twot;
    afcn_tol = imsl_f_max(1.0e-10, pow(eps, twot));
    tol = sqrt_eps;

    code = 1;
    while (code > 0) {
        code = va_arg(argptr, Mint);
        arg_number++;
        switch(code) {
            case IMSL_XGUESS:
                xguess = va_arg(argptr, Mfloat*);
                user_xguess = 1;
                arg_number++;
                break;
            case IMSL_JACOBIAN:
#ifdef ANSI
                jacobian = (void (*)(Mint, Mint, Mfloat[], Mfloat[], Mint))
				va_arg(argptr, void*);
#else
                jacobian = (void (*)()) va_arg(argptr, void*);
#endif
                user_jacobian = 1;
                if (!user_tol) tol = 100.0 * eps;
                arg_number++;
                break;
            case IMSL_XSCALE:
                xscale = va_arg(argptr, Mfloat*);
                user_xscale = 1;
                arg_number++;
                break;
            case IMSL_FSCALE:
                user_fscale = 1;
                fscale = va_arg(argptr, Mfloat*);
                arg_number++;
                break;

            case IMSL_GRAD_TOL:
                grad_tol = (Mfloat) va_arg(argptr, Mdouble);
                arg_number++; 
                break;
            case IMSL_STEP_TOL:
                step_tol = (Mfloat) va_arg(argptr, Mdouble);
                arg_number++;
                break; 
            case IMSL_REL_FCN_TOL: 
                rfcn_tol = (Mfloat) va_arg(argptr, Mdouble); 
                arg_number++; 
                break;  
            case IMSL_ABS_FCN_TOL: 
                afcn_tol = (Mfloat) va_arg(argptr, Mdouble); 
                arg_number++; 
                break;  
            case IMSL_MAX_STEP:
                max_step = (Mfloat) va_arg(argptr, Mdouble);
                arg_number++;
                break;

            case IMSL_GRAD_TOL_ADR:
                grad_tol = *(va_arg(argptr, Mfloat *));
                arg_number++; 
                break;
            case IMSL_STEP_TOL_ADR:
                step_tol = *(va_arg(argptr, Mfloat *));
                arg_number++;
                break; 
            case IMSL_REL_FCN_TOL_ADR: 
                rfcn_tol = *(va_arg(argptr, Mfloat *)); 
                arg_number++; 
                break;  
            case IMSL_ABS_FCN_TOL_ADR: 
                afcn_tol = *(va_arg(argptr, Mfloat *)); 
                arg_number++; 
                break;  
            case IMSL_MAX_STEP_ADR:
                max_step = *(va_arg(argptr, Mfloat *));
                arg_number++;
                break;
            case IMSL_GOOD_DIGIT:
                ndigit = va_arg(argptr, Mint);
                arg_number++; 
                break;
            case IMSL_MAX_ITN:
                maxitn = va_arg(argptr, Mint);
                arg_number++; 
                break;
            case IMSL_MAX_FCN:
                maxfcn = va_arg(argptr, Mint);
                arg_number++; 
                break;
            case IMSL_MAX_JACOBIAN:
                maxjacobian = va_arg(argptr, Mint);
                arg_number++; 
                break;
            case IMSL_INIT_TRUST_REGION:
                trust_region = (Mfloat) va_arg(argptr, Mdouble);
                arg_number++;
                break;
            case IMSL_INIT_TRUST_REGION_ADR:
                trust_region = *(va_arg(argptr, Mfloat *));
                arg_number++;
                break;
            case IMSL_INTERN_SCALE:
                iscale = 1;
                arg_number++;
                break;
            case IMSL_TOLERANCE:
                user_tol = 1;
                tol = (Mfloat) va_arg(argptr, Mdouble);
                arg_number++; 
                break;
            case IMSL_TOLERANCE_ADR:
                user_tol = 1;
                tol = *(va_arg(argptr, Mfloat *));
                arg_number++; 
                break;
            case IMSL_RETURN_USER:
                lv_value = va_arg(argptr, Mfloat*);
                arg_number++;
                return_user = 1;
                break;
            case IMSL_FVEC:
                user_pfvec = 1;
                arg_number++;
                pfvec = va_arg(argptr, Mfloat**);
                break;
            case IMSL_FVEC_USER:
                user_fvec = 1;
                fvec = va_arg(argptr, Mfloat*);
                arg_number++;
                break;
            case IMSL_FJAC:
                user_pfjac = 1;
                pfjac = va_arg(argptr, Mfloat**);
                arg_number++;
                break;
            case IMSL_FJAC_USER:
                user_fjac = 1;
                fjac = va_arg(argptr, Mfloat*);
                arg_number++;
                break;
            case IMSL_FJAC_COL_DIM:
                fjac_col_dim = va_arg(argptr, Mint);
                arg_number++; 
                break;
            case IMSL_RANK:
	        rank = va_arg (argptr, Mint *);
                return_rank = 1;
                arg_number++;
                break;
            case IMSL_JTJ_INVERSE:
                user_pjtj = 1;
                pjtj_inv = va_arg(argptr, Mfloat**);
                arg_number++;
                break;
            case IMSL_JTJ_INVERSE_USER:
                user_jtj = 1;
                jtj_inv = va_arg(argptr, Mfloat*);
                arg_number++;
                break;
            case IMSL_JTJ_INV_COL_DIM:
                jtj_inv_col_dim = va_arg(argptr, Mint);
                arg_number++; 
                break;
	    case IMSL_CHANGE_LOOP_MAXIMUM:
		loop_counter_maximum = va_arg (argptr, Mint);
		arg_number++;
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

    iparam[0] = 1;
    iparam[1] = ndigit;
    iparam[2] = maxitn;
    iparam[3] = maxfcn;
    iparam[4] = maxjacobian;
    iparam[5] = iscale;
    iparam[6] = 100;

    rparam[0] = grad_tol;
    rparam[1] = step_tol;
    rparam[2] = rfcn_tol;
    rparam[3] = afcn_tol;
    rparam[4] = 100*eps;
    rparam[5] = max_step;
    rparam[6] = trust_region;

    if (m <= 0) {
        imsl_e1sti(1, m);
        imsl_ermes(IMSL_TERMINAL, IMSL_NEED_POSITIVE_NUM_FCNS);
    } 
    if (n <= 0) {
        imsl_e1sti(1, n);
        imsl_ermes(IMSL_TERMINAL, IMSL_N_MUST_BE_POSITIVE);
    } 
    if (m < n) {
        imsl_e1sti(1, n);
        imsl_e1sti(2, m);
        /* (5, 4, "The number of variables must be less than or equal */
        /*         to the number of functions while N = %(i1) and     */
        /*         M = %(i2) are given.");                            */
        imsl_ermes(IMSL_TERMINAL, IMSL_TOO_MANY_VARIABLES);
    }
    if (imsl_n1rty(0)) goto RETURN;

    work   = (Mfloat *) imsl_malloc ((n*9+m*3-1)*sizeof(*work));
    iwork  = (Mint *)   imsl_malloc (n*sizeof(*iwork));
    jac    = (Mfloat *) imsl_malloc (m*n*sizeof(*jac));
    jtwork = (Mfloat *) imsl_malloc (m*fjac_col_dim*sizeof(*jtwork));

    if (!user_fvec)
       fvec = (Mfloat *) imsl_malloc (m*sizeof(*fvec));

    if (!user_xguess) {
       xguess = (Mfloat *) imsl_malloc (n*sizeof(*xguess));
       for (i=0; i<n; i++)  xguess[i] = 0.0;
    }

    if (!user_xscale) { 
       xscale = (Mfloat *) imsl_malloc (n*sizeof(*xscale));
       for (i=0; i<n; i++)  xscale[i] = 1.0;
    }

    if (!user_fscale) {
       fscale = (Mfloat *) imsl_malloc (m*sizeof(*fscale));
       for (i=0; i<m; i++)  fscale[i] = 1.0;
    }

    if (lv_value == NULL) {
       lv_value = (Mfloat *) imsl_malloc (n*sizeof(*lv_value));
       if (lv_value == NULL) {
          imsl_e1sti(1, n);
          imsl_e1stl(1, "n");
          imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
          goto FREE_SPACE;
       }
    }

    if (user_jacobian) {
        l_u2lsj(fcn, jacobian, m, n, xguess, xscale, fscale,
                iparam, rparam, lv_value, fvec, jac, m, work,
                iwork, jtwork, fjac_col_dim);
    }
    else {
        l_u2lsf(fcn, m, n, xguess, xscale, fscale, iparam,
                rparam, lv_value, fvec, jac, m, work, iwork);
    }

    if (user_pfvec) *pfvec = fvec;

FREE_SPACE:
    if (work != NULL){
	imsl_free (work);
	work = NULL;
    }
    if (jtwork != NULL)                    imsl_free (jtwork);
    if (fjac != NULL && !user_fjac)        imsl_free (fjac);
    if (!user_xguess && xguess != NULL)    imsl_free (xguess);
    if (!user_xscale && xscale != NULL)    imsl_free (xscale);
    if (!user_fscale && fscale != NULL)    imsl_free (fscale);

RETURN:
    if (user_fjac) {
       for (i = 0; i < m; i++) {
           k = i * fjac_col_dim;
           scopy (n, &jac[i], m, &fjac[k], 1);
       }  
    }
    if (user_pfjac) {
       work   = (Mfloat *) imsl_malloc ((m*fjac_col_dim)*sizeof(*work));
       for (i = 0; i < m; i++) {
           k = i * fjac_col_dim;
           scopy (n, &jac[i], m, &work[k], 1);
       }  
       *pfjac = work;
    }
    if (user_pjtj)
        work_pjtj = (Mfloat *) imsl_malloc ((n*jtj_inv_col_dim)*sizeof(*work_pjtj));

    if (user_pjtj || user_jtj || return_rank) {
       tmp  = (Mfloat *) imsl_malloc ((3*n)*sizeof(*tmp));
       k = 0;
       imsl_l2rrr(&m, &n, jac, &m, &k, iwork, jac, &m, tmp, &tmp[n], &tmp[n*2]);
       for (i = 0; i < n; i++) {
	   k = i * (m+1);
           if (jac[k] < 0.0) sscal(n-i, -1.0, &jac[k], m);
       }
       tolsq = imsl_fi_power(tol, 2);
       irank = n;
       for (i = 0; i < n; i++) {
	   k = i * m;
           j = k + i;
           temp = imsl_fi_power (imsl_snrm2(i+1, &jac[k], 1), 2);
           if (imsl_fi_power(jac[j],2) < tolsq * temp) {
               irank = irank - 1;
	       for (l = i+1; l < n; l++) {
		   ll = l * (m+1);
                   il = i + l * m;
                   imsl_srotg (&jac[ll], &jac[il], &c, &s);
                   if (l != n-1) {
		       ll1 = ll + m;
                       il1 = il + m;
                       imsl_srot (n-l-1, &jac[ll1], m, &jac[il1], m, c, s);
		   }
	       }
               sset(n-i, 0.0, &jac[j], m);
	   }
       }
       if (return_rank) *rank = irank;
       if (user_pjtj || user_jtj) {
          s2 = 1.0;
          imsl_rcovb (n, jac, m, s2, jac, m);
          if (user_pjtj) {
             for (i = 0; i < n; i++) {
                 k = i * jtj_inv_col_dim;
                 scopy (n, &jac[i], m, &work_pjtj[k], 1);
	     }
             *pjtj_inv = work_pjtj;
          }
          else {
             for (i = 0; i < n; i++) {
                 k = i * jtj_inv_col_dim;
                 scopy (n, &jac[i], m, &jtj_inv[k], 1);
             }
          }
       }
    }
    if (tmp != NULL)                               imsl_free (tmp);
    if (iwork != NULL)                             imsl_free (iwork);
    if (jac != NULL)                               imsl_free (jac);
    if (work != NULL && !user_pfjac)               imsl_free (work);
    if (fvec != NULL && !user_fvec && !user_pfvec) imsl_free (fvec);

    if (imsl_n1rty(0) == 5)  {
        if (!return_user && lv_value != NULL)  imsl_free (lv_value);
        lv_value = NULL;
    }
    return (argptr);
}

/* -----------------------------------------------------------------------
    IMSL Name:  U2LSF/DU2LSF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Solve a nonlinear least squares problem using a modified
                Levenberg-Marquardt algorithm and a finite difference
                Jacobian.

    Usage:      CALL U2LSF (FCN, M, N, XGUESS, XSCALE, FSCALE, IPARAM,
                            RPARAM, X, FVEC, FJAC, LDFJAC, WK, IWK)

    Arguments:  See UNLSF/DUNLSF.

    Remarks:    See UNLSF/DUNLSF.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  ------------------------------------------------------------------- */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_u2lsf(void (*fcn) (Mint, Mint, Mfloat*, Mfloat*), Mint 
                    m, Mint n, Mfloat xguess[], Mfloat xscale[], 
                    Mfloat fscale[], Mint iparam[], Mfloat rparam[], 
                    Mfloat x[], Mfloat fvec[], Mfloat fjac[], Mint 
                    ldfjac, Mfloat wk[], Mint iwk[])
#else
static void l_u2lsf(void (*fcn) (Mint, Mint, Mfloat[], Mfloat[]), Mint 
                    m, Mint n, Mfloat xguess[], Mfloat xscale[], 
                    Mfloat fscale[], Mint iparam[], Mfloat rparam[], 
                    Mfloat x[], Mfloat fvec[], Mfloat fjac[], Mint 
                    ldfjac, Mfloat wk[], Mint iwk[])
#endif
#else
static void l_u2lsf(fcn, m, n, xguess, xscale, fscale, iparam,
                    rparam, x, fvec, fjac, ldfjac, wk, iwk)
	void            (*fcn) ();
	Mint            m, n;
	Mfloat           xguess[], xscale[], fscale[];
	Mint             iparam[];
	Mfloat           rparam[], x[], fvec[], fjac[];
	Mint            ldfjac;
	Mfloat           wk[];
	Mint             iwk[];
#endif
{
	imsl_e1psh("U2LSF ");

	/* CALL THE UNCONSTRAINED NONLINEAR LEAST SQUARES SOLVER */

	l_u3lsf(fcn, m, n, xguess, xscale, fscale, iparam, rparam,
		x, fvec, &wk[n * 9 + m * 2 - 1], &wk[0], &wk[n], &wk[n * 2],
		&wk[n * 3], fjac, ldfjac, iwk, &wk[n * 4], &wk[n * 5],
                &wk[n * 6], &wk[n * 8 - 1], &wk[n * 9 - 1], &wk[n * 9 + m - 1]);

	imsl_e1pop("U2LSF ");
	return;
}

/* -----------------------------------------------------------------------
    IMSL Name:  U2LSJ/DU2LSJ (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Solve a nonlinear least squares problem using a modified
                Levenberg-Marquardt algorithm and a user-supplied
                Jacobian.

    Usage:      CALL U2LSJ (FCN, JAC, M, N, XGUESS, XSCALE, FSCALE,
                            IPARAM, RPARAM, X, FVEC, FJAC, LDFJAC, WK,
                            IWK)

    Arguments:  See UNLSJ/DUNLSJ.

    Remarks:    See UNLSJ/DUNLSJ.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_u2lsj(void (*fcn) (Mint, Mint, Mfloat*, Mfloat*), void
                    (*jac) (Mint, Mint, Mfloat*, Mfloat*, Mint), Mint  
                    m, Mint n, Mfloat xguess[], Mfloat xscale[], 
                    Mfloat fscale[], Mint iparam[], Mfloat rparam[], 
                    Mfloat x[], Mfloat fvec[], Mfloat fjac[], Mint 
                    ldfjac, Mfloat wk[], Mint iwk[], Mfloat jtwork[],
                    Mint fjaccd)
#else
static void l_u2lsj(void (*fcn) (Mint, Mint, Mfloat[], Mfloat[]), void
                    (*jac) (Mint, Mint, Mfloat[], Mfloat[], Mint), Mint  
                    m, Mint n, Mfloat xguess[], Mfloat xscale[], 
                    Mfloat fscale[], Mint iparam[], Mfloat rparam[], 
                    Mfloat x[], Mfloat fvec[], Mfloat fjac[], Mint 
                    ldfjac, Mfloat wk[], Mint iwk[], Mfloat jtwork[],
                    Mint fjaccd)
#endif
#else
static void l_u2lsj(fcn, jac, m, n, xguess, xscale, fscale, iparam,
	            rparam, x, fvec, fjac, ldfjac, wk, iwk, jtwork, fjaccd)
	void            (*fcn) (), (*jac) ();
	Mint            m, n;
	Mfloat           xguess[], xscale[], fscale[];
	Mint             iparam[];
	Mfloat           rparam[], x[], fvec[], *fjac;
	Mint            ldfjac;
	Mfloat           wk[];
	Mint             iwk[];
        Mfloat           jtwork[];
        Mint             fjaccd;
#endif
{
#define FJAC(I_,J_)	(fjac+(I_)*(ldfjac)+(J_))


	imsl_e1psh("U2LSJ ");

	/* CALL THE UNCONSTRAINED NONLINEAR LEAST SQUARES SOLVER */

	l_u3lsj(fcn, jac, m, n, xguess, xscale, fscale, iparam, rparam,
		x, fvec, &wk[n * 9 + m * 2 - 1], &wk[0], &wk[n], &wk[n * 2],
		&wk[n * 3], fjac, ldfjac, iwk, &wk[n * 4], &wk[n * 5],
                &wk[n * 6], &wk[n * 8 - 1], &wk[n * 9 - 1], 
                &wk[n * 9 + m - 1], jtwork, fjaccd);

	imsl_e1pop("U2LSJ ");
	return;
}

/* -----------------------------------------------------------------------
    IMSL Name:  U3LSF/DU3LSF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 26, 1988

    Purpose:    Unconstrained nonlinear least squares solver using finite
                difference.

    Usage:      CALL U3LSF (FCN, M, N, XGUESS, XSCALE, FSCALE, IPARAM,
                            RPARAM, XC, FC, FP, XP, SC, GNSTEP, GC, FJC,
                            LDFJC, IPVT, WK1, WK2, WK3, WK4, WK5, WK6)

    Arguments:
       FCN    - User-supplied SUBROUTINE to evaluate the function to be
                minimized.  The usage is
                CALL FCN (M, N, X, F), where
                M      - Length of F.  (Input)
                N      - Length of X.  (Input)
                X      - The point at which the function is evaluated.
                         (Input)
                         X should not be changed by FCN.
                F      - The computed function value at the point X.
                         (Output)
                FCN must be declared EXTERNAL in the calling program.
       M      - Number of functions.  (Input)
       N      - Number of variables.  (Input)
       XGUESS - Vector of length N containing initial guess.  (Input)
       XSCALE - Vector of length N containing the diagonal scaling
                matrix for the variables.  (Input)
       FSCALE - Vector of length N containing the diagonal scaling
                matrix for the functions.  (Input)
       IPARAM - Parameter vector of length 6.  (Input/Output)
                See UNLSF for details.
       RPARAM - Parameter vector of length 7.  (Input/Output)
                See UNLSF for details.
       XC     - Vector of length N containing the approximate solution.
                (Output)
       FC     - Vector of length M containing the residuals at the
                solution.  (Output)
       FP     - Vector of length M containing the updated residual.
                (Output)
       XP     - Vector of length N containing the updated point.
                (Output)
       SC     - Vector of length N containing the last step taken.
                (Output)
       GNSTEP - Vector of length N containing the last Gauss-Newton
                step.  (Output)
       GC     - Vector of length N containing an estimate of the
                scaled gradient at the approximate solution.  (Output)
       FJC    - M by N matrix containing an estimate of the Jacobian
                at the approximate solution.  (Output)
       LDFJC  - Leading dimension of FJC exactly as specified in the
                dimension statement of the calling program.  (Input)
       IPVT   - Vector of length N containing the permutation matrix used
                in the QR factorization of the Jacobian at the
                approximate solution.
       WK1    - Real work vector of length N.
       WK2    - Real work vector of length N.
       WK3    - Real work vector of length 2*N - 1.
       WK4    - Real work vector of length N.
       WK5    - Real work vector of length N.
       WK6    - Real work vector of length M.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
static struct t_u14sf {
	Mfloat           fjactl, steptl, rftol, aftol, falstl;
	Mint             mxiter, maxfcn, maxjac;
}               u14sf;
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_u3lsf(void (*fcn) (Mint, Mint, Mfloat*, Mfloat*),Mint m,
                    Mint n, Mfloat xguess[], Mfloat xscale[], Mfloat 
                    fscale[],Mint iparam[], Mfloat rparam[], Mfloat xc[],
                    Mfloat fc[], Mfloat fp[], Mfloat xp[], Mfloat sc[], 
                    Mfloat gnstep[], Mfloat gc[], Mfloat *fjc, Mint 
                    ldfjc, Mint ipvt[], Mfloat wk1[], Mfloat wk2[], 
                    Mfloat wk3[], Mfloat wk4[], Mfloat wk5[], Mfloat wk6[])
#else
static void l_u3lsf(void (*fcn) (Mint, Mint, Mfloat[], Mfloat[]),Mint m,
                    Mint n, Mfloat xguess[], Mfloat xscale[], Mfloat 
                    fscale[],Mint iparam[], Mfloat rparam[], Mfloat xc[],
                    Mfloat fc[], Mfloat fp[], Mfloat xp[], Mfloat sc[], 
                    Mfloat gnstep[], Mfloat gc[], Mfloat *fjc, Mint 
                    ldfjc, Mint ipvt[], Mfloat wk1[], Mfloat wk2[], 
                    Mfloat wk3[], Mfloat wk4[], Mfloat wk5[], Mfloat wk6[])
#endif
#else
static void l_u3lsf(fcn, m, n, xguess, xscale, fscale, iparam, rparam,
	            xc, fc, fp, xp, sc, gnstep, gc, fjc, ldfjc, ipvt, 
                    wk1, wk2, wk3, wk4, wk5, wk6)
	void            (*fcn) ();
	Mint            m, n;
	Mfloat           xguess[], xscale[], fscale[];
	Mint             iparam[];
	Mfloat           rparam[], xc[], fc[], fp[], xp[], sc[], gnstep[],
	                gc[], *fjc;
	Mint            ldfjc, ipvt[];
	Mfloat           wk1[], wk2[], wk3[], wk4[], wk5[], wk6[];
#endif
{
#define FJC(I_,J_)	(fjc+(I_)*(ldfjc)+(J_))
	Mint            first, gauss, mxtake;
	Mint             _l0, i, icode, iter, j, mode, nfcn, njac;
	Mfloat           amu, delta, eps, epsfcn, fcnorm, fdigit, fpnorm,
	                stepmx;

	imsl_e1psh("U3LSF ");
	/*
	 * CHECK THE VALIDITY OF THE USER SPECIFIED PARAMETERS
	 */
	scopy(n, xguess, 1, xc, 1);
	if (iparam[5] == 1)
		sset(n, 1.0e0, xscale, 1);
	l_u5lsf(m, n, xc, xscale, fscale, 0, iparam, rparam);
	if (imsl_n1rty(1) == 5)
		goto L_9000;

	fdigit = iparam[1];
	u14sf.mxiter = iparam[2];
	u14sf.maxfcn = iparam[3];
	u14sf.maxjac = iparam[4];
	mode = iparam[5];

	u14sf.fjactl = rparam[0];
	u14sf.steptl = rparam[1];
	u14sf.rftol = rparam[2];
	u14sf.aftol = rparam[3];
	u14sf.falstl = rparam[4];
	stepmx = rparam[5];
	delta = rparam[6];
	/*
	 * SET EPSFCN TO ESTIMATE OF RELATIVE NOISE IN FUNCTION
	 */
	eps = imsl_amach(4);
	epsfcn = imsl_f_max(eps, pow(10.0e0, -fdigit));
	/*
	 * INITIALIZE ITERATION, FUNCTION & JACOBIAN EVALUATION COUNTER
	 */
	iter = 0;
	nfcn = 0;
	njac = 0;
	/*
	 * EVALUATE THE FUNCTION & JACOBIAN AT THE INITIAL POINT
	 */
	icode = 0;
	imsl_e1usr("ON");
	(*fcn) (m, n, xc, fc);
	imsl_e1usr("OFF");
	fcnorm = 0.5e0 * imsl_fi_power(imsl_snrm2(m, fc, 1), 2);
	nfcn += 1;
	l_f2jac(fcn, m, n, xc, xscale, fc, epsfcn, fjc, ldfjc, wk6);
	njac += 1;
	mxtake = 0;
	/*
	 * CALCULATE THE SCALED GRADIENT = TRANS(FJC) * FSCALE**2 * FC
	 */
	sset(n, 0.0e0, gc, 1);
	for (j = 1; j <= n; j++) {
		for (i = 1; i <= m; i++) {
			gc[j - 1] += *FJC(j - 1, i - 1) * fscale[i - 1] * fscale[i - 1] *
				fc[i - 1];
		}
	}
	/*
	 * CHECK STOPPING CRITERIA AT THE INITIAL POINT
	 */
	l_u6lsf(m, n, xc, sc, fc, fcnorm, gc, xscale, fscale, &icode,
		   iter, nfcn, njac, 0, mxtake);
	if (imsl_n1rcd(1) != 0 || icode == -999)
		goto L_80;

	/* MAIN ITERATION LOOP */
L_30:
	;
	iter += 1;
	/*
	 * SCALE THE JACOBIAN OF THE RESIDUAL = FSCALE * FJC
	 */
	for (i = 1; i <= m; i++) {
		sscal(n, fscale[i - 1], FJC(0, i - 1), ldfjc);
	}
	/*
	 * COMPUTE THE QR FACTORIZATION OF THE SCALED JACOBIAN AND UPDATE
	 * XSCALE IF VARIABLE INTERNAL SCALING IS REQUESTED.
	 */
	iset(n, 0, ipvt, 1);
        _l0 = 1;
	imsl_l2rrr(&m, &n, fjc, &ldfjc, &_l0, ipvt, fjc, &ldfjc, wk1,
		   wk2, wk3);
	if (mode == 1)
		l_u8lsf(iter, n, wk2, xscale);
	/*
	 * FORM (Q**T)*FSCALE*F AND STORE THE FIRST N COMPONENTS IN WK6
	 */
	l_u11nf(m, fscale, 1, fc, wk6);
	l_u10sf(m, n, fjc, ldfjc, wk1, wk6, wk6);
	/* COMPUTE THE LEVENBERG-MARQUARDT STEP */
	icode = 6;
	first = 1;
L_50:
	l_u7lsf(n, gc, fjc, ldfjc, ipvt, xscale, wk6, stepmx, &delta,
		   &amu, &first, sc, gnstep, &gauss, wk1, wk2, wk3);

	if (imsl_n1rty(1) >= 4)
		goto L_9000;
	/*
	 * CHECK NEW POINT AND UPDATE THE TRUST REGION
	 */
	l_u9lsf(fcn, m, n, xc, fcnorm, gc, fjc, ldfjc, ipvt, sc, xscale,
	     gauss, stepmx, &delta, &icode, wk4, wk5, xp, fc, fp, &fpnorm,
		   &mxtake, &nfcn);
	if (icode >= 4)
		goto L_50;
	/*
	 * EVALUATE THE JACOBIAN AT THE NEW POINT
	 */
	l_f2jac(fcn, m, n, xp, xscale, fp, epsfcn, fjc, ldfjc, wk6);
	njac += 1;
	/*
	 * CALCULATE THE NEW SCALED GRADIENT = TRANS(FJC) * FSCALE**2 * FP
	 */
	sset(n, 0.0e0, gc, 1);
	for (j = 1; j <= n; j++) {
		for (i = 1; i <= m; i++) {
			gc[j - 1] += *FJC(j - 1, i - 1) * fscale[i - 1] * fscale[i - 1] *
				fp[i - 1];
		}
	}
	/* CHECK STOPPING CRITERIA AT NEW POINT */
	l_u6lsf(m, n, xp, sc, fp, fpnorm, gc, xscale, fscale, &icode,
		   iter, nfcn, njac, 0, mxtake);
	if (imsl_n1rcd(1) == 0 && icode != -999) {
		/* UPDATE XC, FC THEN NEXT ITERATION */
		scopy(n, xp, 1, xc, 1);
		scopy(m, fp, 1, fc, 1);
		fcnorm = fpnorm;
		goto L_30;
	}
	/* STOPPING CRITERIA SATISFIED */
	scopy(n, xp, 1, xc, 1);
	scopy(m, fp, 1, fc, 1);
	fcnorm = fpnorm;
L_80:
	iparam[2] = iter;
	iparam[3] = nfcn;
	iparam[4] = njac;

L_9000:
	imsl_e1pop("U3LSF ");
	return;
}				/* end of function */
/* Structured by FOR_STRUCT, v0.2, on 10/02/90 at 10:51:11
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  U3LSJ/DU3LSJ (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 26, 1988

    Purpose:    Unconstrained nonlinear least squares solver using an
                analytic Jacobian.

    Usage:      CALL U3LSJ (FCN, JAC, M, N, XGUESS, XSCALE, FSCALE,
                            IPARAM, RPARAM, XC, FC, FP, XP, SC, GNSTEP,
                            GC, FJC, LDFJC, IPVT, WK1, WK2, WK3, WK4,
                            WK5, WK6)

    Arguments:
       FCN    - User-supplied SUBROUTINE to evaluate the function to be
                minimized.  The usage is
                CALL FCN (M, N, X, F), where
                M      - Length of F.  (Input)
                N      - Length of X.  (Input)
                X      - The point at which the function is evaluated.
                         (Input)
                         X should not be changed by FCN.
                F      - The computed function value at the point X.
                         (Output)
                FCN must be declared EXTERNAL in the calling program.
       JAC    - User-supplied SUBROUTINE to evaluate the Jacobian at a
                point X.  The usage is
                CALL JAC (M, N, X, FJAC, LDFJAC), where
                M      - Length of F.  (Input)
                N      - Length of X.  (Input)
                X      - The point at which the function is evaluated.
                         (Input)
                         X should not be changed by FCN.
                FJAC   - The computed M by N Jacobian at the point X.
                         (Output)
                LDFJAC - Leading dimension of FJAC.  (Input)
                JAC must be declared EXTERNAL in the calling program.
       M      - Number of functions.  (Input)
       N      - Number of variables.  (Input)
       XGUESS - Vector of length N containing initial guess.  (Input)
       XSCALE - Vector of length N containing the diagonal scaling
                matrix for the variables.  (Input)
       FSCALE - Vector of length N containing the diagonal scaling
                matrix for the functions.  (Input)
       IPARAM - Parameter vector of length 6.  (Input/Output)
                See UNLSF for details.
       RPARAM - Parameter vector of length 7.  (Input/Output)
                See UNLSF for details.
       XC     - Vector of length N containing the approximate solution.
                (Output)
       FC     - Vector of length M containing the residuals at the
                solution.  (Output)
       FP     - Vector of length M containing the updated residual.
                (Output)
       XP     - Vector of length N containing the updated point.
                (Output)
       SC     - Vector of length N containing the last step taken.
                (Output)
       GNSTEP - Vector of length N containing the last Gauss-Newton
                step.  (Output)
       GC     - Vector of length N containing an estimate of the
                scaled gradient at the approximate solution.  (Output)
       FJC    - M by N matrix containing an estimate of the Jacobian
                at the approximate solution.  (Output)
       LDFJC  - Leading dimension of FJC exactly as specified in the
                dimension statement of the calling program.  (Input)
       IPVT   - Vector of length N containing the permutation matrix used
                in the QR factorization of the Jacobian at the
                approximate solution.
       WK1    - Real work vector of length N.
       WK2    - Real work vector of length N.
       WK3    - Real work vector of length 2*N - 1.
       WK4    - Real work vector of length N.
       WK5    - Real work vector of length N.
       WK6    - Real work vector of length M.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_u3lsj(void (*fcn) (Mint, Mint, Mfloat*, Mfloat*), void
                    (*jac) (Mint, Mint, Mfloat*, Mfloat*, Mint), Mint m,
                    Mint n, Mfloat xguess[], Mfloat xscale[], Mfloat 
                    fscale[],Mint iparam[], Mfloat rparam[], Mfloat xc[],
                    Mfloat fc[], Mfloat fp[], Mfloat xp[], Mfloat sc[], 
                    Mfloat gnstep[], Mfloat gc[], Mfloat *fjc, Mint 
                    ldfjc, Mint ipvt[], Mfloat wk1[], Mfloat wk2[], 
                    Mfloat wk3[], Mfloat wk4[], Mfloat wk5[], Mfloat wk6[],
                    Mfloat fjct[], Mint fjct_col_dim)
#else
static void l_u3lsj(void (*fcn) (Mint, Mint, Mfloat[], Mfloat[]), void
                    (*jac) (Mint, Mint, Mfloat[], Mfloat[], Mint), Mint m,
                    Mint n, Mfloat xguess[], Mfloat xscale[], Mfloat 
                    fscale[],Mint iparam[], Mfloat rparam[], Mfloat xc[],
                    Mfloat fc[], Mfloat fp[], Mfloat xp[], Mfloat sc[], 
                    Mfloat gnstep[], Mfloat gc[], Mfloat *fjc, Mint 
                    ldfjc, Mint ipvt[], Mfloat wk1[], Mfloat wk2[], 
                    Mfloat wk3[], Mfloat wk4[], Mfloat wk5[], Mfloat wk6[],
                    Mfloat fjct[], Mint fjct_col_dim)
#endif
#else
static void l_u3lsj(fcn, jac, m, n, xguess, xscale, fscale, iparam,
	            rparam, xc, fc, fp, xp, sc, gnstep, gc, fjc, ldfjc,
                    ipvt, wk1, wk2, wk3, wk4, wk5, wk6, fjct, fjct_col_dim)
	void            (*fcn) (), (*jac) ();
	Mint            m, n;
	Mfloat          xguess[], xscale[], fscale[];
	Mint            iparam[];
	Mfloat          rparam[], xc[], fc[], fp[], xp[], sc[], gnstep[],
	                gc[], *fjc;
	Mint            ldfjc, ipvt[];
	Mfloat          wk1[], wk2[], wk3[], wk4[], wk5[], wk6[], fjct[];
        Mint            fjct_col_dim;
#endif
{
#define FJC(I_,J_)	(fjc+(I_)*(ldfjc)+(J_))
	Mint             first, gauss, mxtake;
	Mint             _l0, i, icode, iter, j, mode, nfcn, njac;
	Mfloat           amu, delta, eps, fcnorm, fdigit, fpnorm,
	                 stepmx;

	imsl_e1psh("U3LSJ ");

	/* CHECK THE VALIDITY OF THE USER SPECIFIED PARAMETERS */

	scopy(n, xguess, 1, xc, 1);
	if (iparam[5] == 1)
		sset(n, 1.0e0, xscale, 1);
	l_u5lsf(m, n, xc, xscale, fscale, 1, iparam, rparam);
	if (imsl_n1rty(1) == 5)
		goto L_9000;

	fdigit = iparam[1];
	u14sf.mxiter = iparam[2];
	u14sf.maxfcn = iparam[3];
	u14sf.maxjac = iparam[4];
	mode = iparam[5];

	u14sf.fjactl = rparam[0];
	u14sf.steptl = rparam[1];
	u14sf.rftol = rparam[2];
	u14sf.aftol = rparam[3];
	u14sf.falstl = rparam[4];
	stepmx = rparam[5];
	delta = rparam[6];
	/*
	 * SET EPSFCN TO ESTIMATE OF RELATIVE NOISE IN FUNCTION
	 */
	eps = imsl_amach(4);
	/*
	 * INITIALIZE ITERATION, FUNCTION & JACOBIAN EVALUATION COUNTER
	 */
	iter = 0;
	nfcn = 0;
	njac = 0;
	/*
	 * EVALUATE THE FUNCTION & JACOBIAN AT THE INITIAL POINT
	 */
	icode = 0;
	imsl_e1usr("ON");
	(*fcn) (m, n, xc, fc);
	(*jac) (m, n, xc, fjct, fjct_col_dim);
	imsl_e1usr("OFF");
	fcnorm = 0.5e0 * imsl_fi_power(imsl_snrm2(m, fc, 1), 2);
	nfcn += 1;
	njac += 1;
	mxtake = 0;

        /* fjc = transpose(fjct) */

#if 0
        for (i=0; i<m; i++) {
             j = i * fjct_col_dim;
             scopy(m, &fjct[i], n, &fjc[j], 1);
        }
#endif
	imsl_f_m1ran (m, n, fjct, fjc);

	/*
	 * CALCULATE THE SCALED GRADIENT = TRANS(FJC) * FSCALE^2 * FC
	 */
	sset(n, 0.0e0, gc, 1);
	for (j = 0; j < n; j++) {
	     for (i = 0; i < m; i++) {
	          gc[j] += *FJC(j, i) * fscale[i] * fscale[i] * fc[i];
		}
	}
	/*
	 * CHECK STOPPING CRITERIA AT THE INITIAL POINT
	 */
	l_u6lsf(m, n, xc, sc, fc, fcnorm, gc, xscale, fscale, &icode,
		   iter, nfcn, njac, 1, mxtake);
	if (imsl_n1rcd(1) != 0 || icode == -999)
		goto L_80;

	/* MAIN ITERATION LOOP */
L_30:
	;
	iter += 1;
	/*
	 * SCALE THE JACOBIAN OF THE RESIDUAL = FSCALE * FJC
	 */
	for (i = 1; i <= m; i++) {
		sscal(n, fscale[i - 1], FJC(0, i - 1), ldfjc);
	}
	/*
	 * COMPUTE THE QR FACTORIZATION OF THE SCALED JACOBIAN AND UPDATE
	 * XSCALE IF VARIABLE INTERNAL SCALING IS REQUESTED.
	 */
	iset(n, 0, ipvt, 1);
        _l0 = 1;
	imsl_l2rrr(&m, &n, fjc, &ldfjc, &_l0, ipvt, fjc, &ldfjc, wk1,
		   wk2, wk3);
	if (mode == 1)
		l_u8lsf(iter, n, wk2, xscale);
	/*
	 * FORM (Q**T)*FSCALE*F AND STORE THE FIRST N COMPONENTS IN WK6
	 */
	l_u11nf(m, fscale, 1, fc, wk6);
	l_u10sf(m, n, fjc, ldfjc, wk1, wk6, wk6);
	/* COMPUTE THE LEVENBERG-MARQUARDT STEP */
	icode = 6;
	first = 1;
L_50:
	l_u7lsf(n, gc, fjc, ldfjc, ipvt, xscale, wk6, stepmx, &delta,
		   &amu, &first, sc, gnstep, &gauss, wk1, wk2, wk3);

	if (imsl_n1rty(1) >= 4)
		goto L_9000;

	/*
	 * CHECK NEW POINT AND UPDATE THE TRUST REGION
	 */
	l_u9lsf(fcn, m, n, xc, fcnorm, gc, fjc, ldfjc, ipvt, sc, xscale,
	     gauss, stepmx, &delta, &icode, wk4, wk5, xp, fc, fp, &fpnorm,
		   &mxtake, &nfcn);
	if (icode >= 4)
		goto L_50;
	/*
	 * EVALUATE THE JACOBIAN AT THE NEW POINT
	 */
	imsl_e1usr("ON");
	(*jac) (m, n, xp, fjct, fjct_col_dim);
	imsl_e1usr("OFF");
	njac += 1;

#if 0
        /* fjc = transpose(fjct) */
        for (i=0; i<m; i++) {
             j = i * fjct_col_dim;
             scopy(m, &fjct[i], n, &fjc[j], 1);
        }
#endif
	imsl_f_m1ran (m, n, fjct, fjc);

	/*
	 * CALCULATE THE NEW SCALED GRADIENT = TRANS(FJC) * FSCALE**2 * FP
	 */
	sset(n, 0.0e0, gc, 1);
	for (j = 1; j <= n; j++) {
		for (i = 1; i <= m; i++) {
			gc[j - 1] += *FJC(j - 1, i - 1) * imsl_fi_power(fscale[i - 1], 2) *
				fp[i - 1];
		}
	}
	/* CHECK STOPPING CRITERIA AT NEW POINT */
	l_u6lsf(m, n, xp, sc, fp, fpnorm, gc, xscale, fscale, &icode,
		   iter, nfcn, njac, 1, mxtake);
	if (imsl_n1rcd(1) == 0 && icode != -999) {
		/* UPDATE XC, FC THEN NEXT ITERATION */
		scopy(n, xp, 1, xc, 1);
		scopy(m, fp, 1, fc, 1);
		fcnorm = fpnorm;
		goto L_30;
	}
	/* STOPPING CRITERIA SATISFIED */
	scopy(n, xp, 1, xc, 1);
	scopy(m, fp, 1, fc, 1);
	fcnorm = fpnorm;
L_80:
	iparam[2] = iter;
	iparam[3] = nfcn;
	iparam[4] = njac;

L_9000:
	imsl_e1pop("U3LSJ ");
	return;
}				/* end of function */
/* Structured by FOR_STRUCT, v0.2, on 09/28/90 at 16:07:29
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  U5LSF/DU5LSF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Check validity of input to nonlinear least squares solver

    Usage:      CALL U5LSF (M, N, XC, XSCALE, FSCALE, USRJAC, IPARAM,
                            RPARAM)

    Arguments:
       M      - Number of functions.  (Input)
       N      - Number of variables.  (Input)
       XC     - Real vector of length N containing the current point.
                 (Input)
       XSCALE - Real vector of length N containing the diagonal scaling
                matrix for the variables.  (Input/Output)
       FSCALE - Real vector of length M containing the diagonal scaling
                matrix for the functions.  (Input/Output)
       USRJAC - Logical variable.  (Input)
                USRJAC = .TRUE. if analytic Jacobian is used.
                USRJAC = .FALSE. otherwise.
       IPARAM - Integer parameters vector of length 6.  (Input / Output)
                See UNLSF or UNLSJ for details.
       RPARAM - Real parameters vector of length 7.  (Input / Output)
                See UNLSF or UNLSJ for details.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_u5lsf(Mint m, Mint n, Mfloat xc[], Mfloat xscale[], 
                    Mfloat fscale[], Mint usrjac, Mint iparam[],
	            Mfloat rparam[])
#else
static void l_u5lsf(m, n, xc, xscale, fscale, usrjac, iparam, rparam)
	Mint            m, n;
	Mfloat           xc[], xscale[], fscale[];
	Mint             usrjac;
	Mint             iparam[];
	Mfloat           rparam[];
#endif
{
	Mint             i; 
	Mfloat           arftol, eps, temp, temp1, temp2;


	imsl_e1psh("U5LSF ");

	eps = imsl_amach(4);
	/*
	 * FDIGIT Number of good digits in F MXITER Maximum number of
	 * iterations MAXFCN Max number of fcn. evaluations MAXJAC Max number
	 * of Jac. evaluations MODE   Flag for internal scaling GRADTL Scaled
	 * gradient tolerance STEPTL Scaled step tolerance RFTOL  Relative
	 * function tolerance AFTOL  Absolute function tolerance FALSTL False
	 * convergence tolerance STEPMX Maximum allowable step size DELTA
	 * Size of initial trust region
	 * 
	 * CHECK VARIABLE SCALING MATRIX ONLY IF EXTERNAL SCALING REQUESTED
	 */
	if (iparam[5] != 1) {
		for (i = 1; i <= n; i++) {
			if (xscale[i - 1] <= 0.0e0) {
				/* Print error message */
				l_u13sf(2);
				/*
				 * C               CALL E1MES (6, 2, 'The
				 * values for the scaling matrix, ' C     &
				 * 'XSCALE, must be positive while at
				 * least'// C     &                    ' one
				 * entry is less than or equal to zero'// C
				 * &                     '.  The algorithm
				 * will use the identity '// C     &
				 * 'scaling matrix for XSCALE.')
				 */
				sset(n, 1.0e0, xscale, 1);
				goto L_20;
			}
		}
	}
	/* CHECK FUNCTION SCALING */
L_20:
	for (i = 1; i <= m; i++) {
		if (fscale[i - 1] <= 0.0e0) {
			/* Print error message */
			l_u13sf(3);
			/*
			 * C            CALL E1MES (6, 3, 'The values for the
			 * diagonal matrix, '// C     &
			 * 'FSCALE, must be positive while at least '// C
			 * &                  'one entry is less than or
			 * equal to zero. '// C     &                'The
			 * algorithm will use the identity scaling'// C     &
			 * ' matrix for FSCALE.')
			 */
			sset(m, 1.0e0, fscale, 1);
			goto L_40;
		}
	}
	/* CHECK ACCURACY OF THE PROBLEM */
L_40:
	if (iparam[1] <= 0) {
		imsl_e1sti(1, iparam[1]);
		/* Print error message */
		l_u13sf(4);
		/*
		 * C         CALL E1MES (6, 4, 'The estimate of the number of
		 * '// C     &               'good digits in the functions
		 * must be positive'// C     &              ' while FDIGIT =
		 * %(I1) is given.  The algorithm'// C     &               '
		 * will assume that the function is accurate to'// C     &
		 * ' the precision of the arithmetic.')
		 */
		iparam[1] = imsl_i_machine(7);
	}
	/* CHECK MAXIMUM NUMBER OF ITERATIONS */
	if (iparam[2] <= 0) {
		imsl_e1sti(1, iparam[2]);
		/* Print error message */
		l_u13sf(5);
		/*
		 * C         CALL E1MES (6, 5, 'The maximum number of
		 * iterations '// C     &               'must be positive
		 * while MXITER = %(I1) is '// C     &               'given.
		 * The algorithm will use MXITER = 100.')
		 */
		iparam[2] = 100;
	}
	/* CHECK MAXIMUM FUNCTION EVALUATIONS */
	if (iparam[3] <= 0) {
		imsl_e1sti(1, iparam[3]);
		/* Print error message */
		l_u13sf(6);
		/*
		 * C         CALL E1MES (6, 6, 'The maximum number of
		 * function '// C     &               'evaluations must be
		 * positive while MAXFCN = '// C     &               '%(I1)
		 * is given.  The algorithm will use '// C     &
		 * 'MAXFCN = 400.')
		 */
		iparam[3] = 400;
	}
	/* CHECK MAXIMUM JACOBIAN EVALUATIONS */
	if (usrjac) {
		if (iparam[4] <= 0) {
			imsl_e1sti(1, iparam[4]);
			/* Print error message */
			l_u13sf(7);
			/*
			 * C            CALL E1MES (6, 7, 'The maximum number
			 * of Jacobian '// C     &
			 * 'evaluations must be positive while MAXJAC ='// C
			 * &                  ' %(I1) is given.  The
			 * algorithm will use '// C     &
			 * 'MAXJAC = 100.')
			 */
			iparam[4] = 400;
		}
	}
	temp1 = pow(eps, 1.0e0 / 3.0e0);
	temp2 = pow(eps, 2.0e0 / 3.0e0);
	/* CHECK THE GRADIENT TOLERANCE */
	if (rparam[0] < 0.0e0) {
		imsl_e1str(1, rparam[0]);
		imsl_e1str(2, temp2);
		/* Print error message */
		l_u13sf(8);
		/*
		 * C         CALL E1MES (6, 8, 'The gradient tolerance must
		 * be '// C     &               'nonnegative while GRADTL =
		 * %(R1) is given.  '// C     &               'The algorithm
		 * will use GRADTL = %(R2).')
		 */
		rparam[0] = temp2;
	}
	/* CHECK THE STEP TOLERANCE */
	if (rparam[1] < 0.0e0) {
		imsl_e1str(1, rparam[1]);
		imsl_e1str(2, temp2);
/*		(6, 9, "The step tolerance must be nonnegative while STEPTOL = %(r1) is given.  The algorithm will use STEPTOL = %(r2).");
*/
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEGATIVE_STEP_TOL);
		rparam[1] = temp2;
	}
	/* CHECK RELATIVE FUNCTION TOLERANCE */
	arftol = 1.0e-10;
	if (rparam[2] < 0.0e0) {
		temp = imsl_f_max(arftol, temp2);
		imsl_e1str(1, rparam[2]);
		imsl_e1str(2, temp);

/*		(6, 10, "The relative function tolerance must be nonnegative while RFCFTOL = %(r1) is given.  The algorithm will use RFCNTOL = %(r2).");
*/
                imsl_ermes(IMSL_WARNING_IMMEDIATE,  IMSL_NEGATIVE_REL_FCN_TOL);
		rparam[2] = temp;
	}
	/* CHECK ABSOLUTE FUNCTION TOLERANCE */
	if (rparam[3] < 0.0e0) {
		temp = imsl_f_max(arftol * arftol, eps * eps);
		imsl_e1str(1, rparam[3]);
		imsl_e1str(2, temp);

/*		(6, 11, "The absolute function tolerance must be nonnegative while AFCNTOL = %(r1) is given.  The algorithm will use AFCNTOL = %(r2).");
*/
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEGATIVE_ABS_FCN_TOL);
		rparam[3] = temp;
	}
	/* CHECK FALSE CONVERGENCE TOLERANCE */
	if (rparam[4] < 0.0e0) {
		temp = 1.0e2 * eps;
		rparam[4] = temp;
	}
	/* CHECK MAXIMUM ALLOWED STEP SIZE */
	if (rparam[5] <= 0.0e0) {
		temp1 = 0.0e0;
		for (i = 1; i <= n; i++) {
			temp1 += imsl_fi_power(xscale[i - 1] * xc[i - 1], 2);
		}
		temp1 = sqrt(temp1);
		temp2 = imsl_snrm2(n, xscale, 1);
		temp = 1.0e3 * imsl_f_max(temp1, temp2);
		if (iparam[0] != 0 && rparam[5] != -999.0e0) {
			imsl_e1str(1, rparam[5]);
			imsl_e1str(2, temp);

/*			(6, 13, "The maximum allowable scaled step length must be positive while STEPMAX = %(r1) is given.  The algorithm will use STEPMAX = %(r2).");
*/
                        imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEED_NONNEGATIVE_STEPMX);
		}
		rparam[5] = temp;
	}
	/* CHECK INITIAL TRUST REGION RADIUS */
	if (rparam[6] <= 0.0e0) {
		if (iparam[0] != 0 && rparam[6] != -999.0e0) {
			imsl_e1str(1, rparam[6]);

/*			(6, 14, "The initial trust region radius must be positive while DELTA = %(r1) is given.  The algorithm will use the length of the initial scaled Cauchy step for DELTA.");
*/
                        imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEED_NONNEGATIVE_DELTA);
		}
		rparam[6] = -999.0e0;
	}
	imsl_e1pop("U5LSF ");
	return;
}				/* end of function */
/* Structured by FOR_STRUCT, v0.2, on 09/28/90 at 16:08:39
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  U6LSF/DU6LSF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Stopping conditions for the nonlinear least squares
                problem.

    Usage:      CALL U6LSF (M, N, XP, SC, FP, FPNORM, GP, XSCALE, FSCALE,
                            ICODE, ITER, NFCN, NJAC, USRJAC, MXTAKE)

    Arguments:
       M      - Number of functions.  (Input)
       N      - Number of variables.  (Input)
       XP     - Real vector of length N containing the new iterate.
                   (Input)
       SC     - Real vector of length N containing step taken.  (Input)
       FP     - Real vector of length M containing F(XP).  (Input)
       FPNORM - Real scalar containing the 2-norm of F(XP).  (Input)
       GP     - Real vector of length N containing the gradient at XP.
                   (Input)
       XSCALE - Real vector of length N containing the diagonal scaling
                matrix for the variables.  (Input)
       FSCALE - Real vector of length M containing the diagonal scaling
                matrix for the functions.  (Input)
       ICODE  - Return code from the global strategy algorithm.  (Input)
       ITER   - Number of iterations.  (Input)
       NFCN   - Number of function evaluations.  (Input)
       NJAC   - Number of Jacobian evaluations.  (Input)
       USRJAC - Logical variable.  (Input)
                USRJAC = .TRUE. if analytic Jacobian is used.
                USRJAC = .FALSE. otherwise.
       MXTAKE - Logical variable indicating a step of maximum length was
                taken.  (Input)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
/* M, FP, and FSCALE are not used, but leave the calling sequence intact.*/
#ifdef ANSI
static void l_u6lsf(Mint m, Mint n, Mfloat xp[], Mfloat sc[], Mfloat 
                    fp[], Mfloat fpnorm, Mfloat gp[], Mfloat xscale[], 
                    Mfloat fscale[], Mint *icode, Mint iter, Mint nfcn,
                    Mint njac, Mint usrjac, Mint mxtake)
#else
static void l_u6lsf(m, n, xp, sc, fp, fpnorm, gp, xscale, fscale,
	            icode, iter, nfcn, njac, usrjac, mxtake)
	Mint            m, n;
	Mfloat           xp[], sc[], fp[], fpnorm, gp[], xscale[], fscale[];
	Mint            *icode, iter, nfcn, njac;
        Mint            usrjac, mxtake;
#endif
{
	Mint             i;
	static Mint      nmaxs;
	Mfloat           big, scgrad, scstep, small, valmax;


	imsl_e1psh("U6LSF ");
	/*
	 * CHECK ABSOLUTE FUNCTION CONVERGENCE TOLERANCE
	 */
	if (fpnorm <= u14sf.aftol) {
		*icode = -999;
		goto L_9000;
	}
	/* TEST OF NORM OF SCALED GRADIENT */
	valmax = 0.0e0;
	small = imsl_amach(1);
	big = imsl_amach(2);
	if (big * small < 1.0e0)
		small = 1.0e0 / big;
	for (i = 1; i <= n; i++) {
		if (fpnorm <= small) {
			scgrad = fabs(gp[i - 1]) * imsl_f_max(fabs(xp[i - 1]), 1.0e0 /
							      xscale[i - 1]);
		} else {
			scgrad = fabs(gp[i - 1]) * imsl_f_max(fabs(xp[i - 1]), 1.0e0 /
						   xscale[i - 1]) / fpnorm;
		}
		valmax = imsl_f_max(scgrad, valmax);
	}
	if (valmax <= u14sf.fjactl) {
		*icode = -999;
		goto L_9000;
	}
	/*
	 * IF FIRST ITER., INITIALIZE COUNTER FOR AMAX1. STEP TAKEN AND
	 * RETURN
	 */
	if (iter == 0) {
		nmaxs = 0;
		goto L_9000;
	}
	/* TEST OF NORM OF SCALED STEP */
	valmax = 0.0e0;
	for (i = 1; i <= n; i++) {
		scstep = fabs(sc[i - 1]) / imsl_f_max(fabs(xp[i - 1]), 1.0e0 /
						      xscale[i - 1]);
		valmax = imsl_f_max(scstep, valmax);
	}
	if (valmax <= u14sf.steptl) {
		*icode = -999;
                imsl_ermes(IMSL_NOTE, IMSL_STEP_TOLERANCE);
		goto L_9000;
	}
	/*
	 * CHECK RELATIVE FUNCTION CONVERGENCE TOLERANCE
	 */
	if (*icode == 2) {
		imsl_e1str(1, u14sf.rftol);

/*		(3, 1, "RELATIVE FUNCTION CONVERGENCE - Both the scaled actual and predicted reductions in the function are less than or equal to the relative function convergence tolerance RFTOL = %(r1)."); */
                imsl_ermes(IMSL_WARNING, IMSL_LITTLE_FCN_CHANGE);
		goto L_9000;
	}
	/* CHECK FALSE CONVERGENCE TOLERANCE */
	if (*icode == 3) {

/*		(3, 2, "FALSE CONVERGENCE - The iterates appear to be converging to a noncritical point.  Incorrect gradient information, a discontinuous function, or stopping tolerances being too tight may be the cause.");
*/
                imsl_ermes(IMSL_FATAL, IMSL_FALSE_CONVERGENCE);
		goto L_9000;
	}
	/*
	 * CHECK ITERATION, FUNCTION, GRADIENT & HESSIAN EVALUATIONS LIMIT
	 */
	if (iter >= u14sf.mxiter) {
                imsl_ermes(IMSL_WARNING, IMSL_TOO_MANY_ITN);

	} else if (nfcn >= u14sf.maxfcn) {
                imsl_ermes(IMSL_WARNING, IMSL_TOO_MANY_FCN_EVAL);

	} else if (usrjac && (njac >= u14sf.maxjac)) {
                imsl_ermes(IMSL_WARNING, IMSL_TOO_MANY_JACOBIAN_EVAL);

	} else if (mxtake) {
		nmaxs += 1;
		if (nmaxs == 5) {

/*			(3, 6, " Five consecutive steps of length STEPMAX have been taken; either the function is unbounded below, or has a finite asymptote in some direction, or STEPMAX is too small."); */
                        imsl_ermes(IMSL_WARNING, IMSL_UNBOUNDED);
		}
	}
L_9000:
	imsl_e1pop("U6LSF ");
	return;
}				/* end of function */
/* Structured by FOR_STRUCT, v0.2, on 09/28/90 at 16:09:30
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  U7LSF/DU7LSF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Compute the Levenberg-Marquardt step.

    Usage:      CALL U7LSF (N, GC, A, LDA, IPVT, XSCALE, QTF, STEPMX,
                            DELTA, AMU, FIRST, SC, GNSTEP, GAUSS, DIAG,
                            WK1, WK2)

    Arguments:
       N      - Dimension of the problem.  (Input)
       GC     - Real vector of length N containing the scaled gradient
                of the residual vector.  (Input)
       A      - Real N by N array.  (Input/Output)
                On input the full upper triangle must contain the full
                   upper triangle of the matrix R resulting from the QR
                   factorization of the scaled Jacobian.
                On output the full upper triangle is unaltered, and the
                   strict lower triangle contains the strict lower
                   triangle of the lower triangular matrix L which is
                   essentially equal to the Cholesky factor of the
                   matrix (J**T)*J + AMU*XSCALE.
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       IPVT   - Integer array of length N containing the pivoting info-
                mations from the QR factorization routine.  (Input)
       XSCALE - Real vector of length N containing the diagonal scaling
                matrix for the variables.  (Input)
       QTF    - Real vector of length N containing the first N elements
                of Q(transpose) * (the scaled residual vector).  (Input)
       STEPMX - Real scalar containing the maximum allowable step size.
                   (Input)
       DELTA  - Real scalar containing the trust region radius with value
                retained between calls.  (Input/Output)
       AMU    - Scalar containing an initial estimate of the Levenberg
                Marquardt parameter on.  (Input/Output)
                On output, AMU contains the final estimate of the L-M
                parameter
       FIRST  - Logical variable.  (Input/Output)
                FIRST = .TRUE. if this is the first call to this routine
                               in this iteration.
                FIRST = .FALSE. otherwise.
       SC     - Real vector of length N containing the LM step.  (Output)
       GNSTEP - Real vector of length N containing the Gauss-Newton step.
                   (Output)
       GAUSS  - Logical variable.  (Output)
                GAUSS = .TRUE. if Gauss-Newton step was taken.
                GAUSS = .FALSE. if Levenberg-Marquardt step was taken.
       DIAG   - Vector of length N containing the diagonal elements of
                the Cholesky factor of (J**T)*J + AMU*XSCALE.  (Output)
       WK1    - Real work vector of length N.
       WK2    - Real work vector of length N.

    Remark:
       This is based on the combined idea of the MINPACK routine LMPAR
       together with the internal doubling strategy in Dennis-Schnabel.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_u7lsf(Mint n, Mfloat gc[], Mfloat *a, Mint lda, Mint 
                    ipvt[], Mfloat xscale[], Mfloat qtf[], Mfloat 
                    stepmx, Mfloat *delta, Mfloat *amu, Mint *first, 
                    Mfloat sc[], Mfloat gnstep[], Mint *gauss, Mfloat 
                    diag[], Mfloat wk1[], Mfloat wk2[])
#else
static void l_u7lsf(n, gc, a, lda, ipvt, xscale, qtf, stepmx,
	   delta, amu, first, sc, gnstep, gauss, diag, wk1, wk2)
	Mint            n;
	Mfloat           gc[], *a;
	Mint            lda, ipvt[];
	Mfloat           xscale[], qtf[], stepmx, *delta, *amu;
	Mint             *first;
	Mfloat           sc[], gnstep[];
	Mint             *gauss;
	Mfloat           diag[], wk1[], wk2[];
#endif
{
#define A(I_,J_)	(a+(I_)*(lda)+(J_))
	Mint            done;
	Mint             i, j, l;
	static Mint      nsing;
	Mfloat           alow, alpha, amulow, amuup, imsl_beta, big, eps,
	                high, small, stplen, temp;
	static Mfloat    deltap, gnleng, phi, phip, phipi, sgnorm;
	Mint 		loop_counter;

	imsl_e1psh ("l_u7lsf");

	small = imsl_amach(1);
	big = imsl_amach(2);
	if (big * small < 1.0e0)
		small = 1.0e0 / big;
	eps = imsl_amach(4);
	/*
	 * IF INITIAL TRUST REGION IS NOT PROVIDED BY THE USER, COMPUTE AND
	 * USE LENGTH OF CAUCHY STEP. BETA = NORM2(R*TRANS(P)*D**(-2)*G)**2
	 */
	if (*delta == -999.0e0) {
		*amu = 0.0e0;
		l_u11nf(n, xscale, -1, gc, wk1);
		alpha = imsl_snrm2(n, wk1, 1);
		imsl_beta = 0.0e0;

		for (i = 1; i <= n; i++) {
			temp = 0.0e0;
			for (j = i; j <= n; j++) {
				l = ipvt[j - 1];
				temp += (*A(j - 1, i - 1) * gc[l - 1]) / (xscale[l - 1] *
							     xscale[l - 1]);
			}
                        wk1[i - 1] = temp;
		}
		imsl_beta = imsl_snrm2(n, wk1, 1);

		if (imsl_beta <= small) {
			*delta = alpha * sqrt(alpha);
		} else {
			*delta = alpha / imsl_beta;
			*delta = *delta * alpha;
			*delta = *delta / imsl_beta;
			*delta = *delta * alpha;
		}
		*delta = imsl_f_min(*delta, stepmx);
	}
	/*
	 * THE FOLLOWING IS DONE ONLY ON THE FIRST TIME THROUGH THIS ITER.:
	 * 1. COMPUTE THE GAUSS-NEWTON STEP. IF JACOBIAN IS RANK-DEFICIENT,
	 * OBTAIN A LEAST SQUARES SOLUTION. 2. COMPUTE THE LENGTH OF THE
	 * SCALED GAUSS-NEWTON STEP. 3. COMPUTE THE 2-NORM OF THE SCALED
	 * GRADIENT USED IN COMPUTING AN UPPER BOUND FOR AMU.
	 */
	if (*first) {
		nsing = n;
		for (j = 1; j <= n; j++) {
			if (*A(j - 1, j - 1) == 0.0e0 && nsing == n)
				nsing = j - 1;
			if (nsing < n)
				wk1[j - 1] = 0.0e0;
		}
		l_u12sf(nsing, a, lda, qtf, wk1);
		for (j = 1; j <= n; j++) {
			gnstep[ipvt[j - 1] - 1] = -wk1[j - 1];
		}
		/* LENGTH OF SCALED GAUSS-NEWTON STEP */
		l_u11nf(n, xscale, 1, gnstep, wk1);
		gnleng = imsl_snrm2(n, wk1, 1);
		/* 2-NORM OF THE SCALED GRADIENT */
		l_u11nf(n, xscale, -1, gc, wk1);
		sgnorm = imsl_snrm2(n, wk1, 1);
	}
	/* BOUNDS ON THE COMPUTED STEP */
	high = 1.5e0;
	alow = 0.75e0;

	if (gnleng <= high ** delta) {
		/* THE STEP IS A GAUSS-NEWTON STEP */
		*gauss = 1;
		scopy(n, gnstep, 1, sc, 1);
		*amu = 0.0e0;
		*delta = imsl_f_min(*delta, gnleng);
	} else {
		/*
		 * FIND A NONTRIVIAL STEP; COMPUTE STARTING VALUE OF AMU IF
		 * PREVIOUS STEP WASN'T A GAUSS-NEWTON STEP
		 */
		*gauss = 0;
		if (*amu > 0.0e0)
			*amu += -((phi + deltap) / *delta) * (((deltap - *delta) +
							       phi) / phip);
		phi = gnleng - *delta;
		/*
		 * IF THE JACOBIAN IS NOT RANK DEFICIENT, THE NEWTON STEP
		 * PROVIDES A LOWER BOUND FOR AMU. OTHERWISE SET THIS BOUND
		 * TO ZERO.
		 */
		if (nsing == n) {
			if (*first) {
				*first = 0;
				for (j = 1; j <= n; j++) {
					l = ipvt[j - 1];
					wk1[j - 1] = xscale[l - 1] * xscale[l - 1] * gnstep[l - 1];
				}
				/*
				 * OBTAIN TRANS(R**-1)*(TRANS(P)*S) BY
				 * SOLVING TRANS(R)*WK1 = WK1
				 */
				l_u12sf(n, a, lda, wk1, wk1);
				phipi = -imsl_fi_power(imsl_snrm2(n, wk1, 1), 2) / gnleng;
			}
			amulow = -phi / phipi;
		} else {
			*first = 0;
			amulow = 0.0e0;
		}
		amuup = sgnorm / *delta;
		/*
		 * LOOP TO CHECK AMU AND GENERATE NEXT AMU IF NECESSARY
		 */
		done = 0;
		loop_counter = 0;
		/* DOUNTIL */
L_10022:
		;

		if (*amu < amulow || *amu > amuup)
			*amu = imsl_f_max(sqrt(amulow * amuup), 1.0e-3 * amuup);
		temp = sqrt(*amu);
		imsl_svcal(n, temp, xscale, 1, wk1, 1);
		/*
		 * SOLVE THE DAMPED LEAST SQUARES SYSTEM FOR THE L-M STEP
		 * USING MORE*S TECHNIQUE IN MINPACK
		 */
		l_u11sf(n, a, lda, ipvt, wk1, qtf, sc, diag, wk2);
		sscal(n, -1.0e0, sc, 1);
		l_u11nf(n, xscale, 1, sc, wk2);
		stplen = imsl_snrm2(n, wk2, 1);
		phi = stplen - *delta;

		for (j = 1; j <= n; j++) {
			l = ipvt[j - 1];
			wk1[j - 1] = xscale[l - 1] * wk2[l - 1];
		}
		for (j = 1; j <= n; j++) {
			if (fabs(diag[j - 1]) >= small)
				wk1[j - 1] /= diag[j - 1];
			if (j < n)
				saxpy(n - j, -wk1[j - 1], A(j - 1, j), 1, &wk1[j],
					   1);
		}
		phip = -imsl_fi_power(imsl_snrm2(n, wk1, 1), 2) / stplen;

		if ((stplen >= alow ** delta && stplen <= high ** delta) || (amuup -
							   amulow) <= eps) {
			done = 1;
		} else {
			/*
			 * SC IS NOT ACCEPTABLE; UPDATE AMU, AMULOW AND AMUUP
			 */
			amulow = imsl_f_max(amulow, *amu - (phi / phip));
			if (phi < 0.0e0)
				amuup = *amu;
			*amu += -(stplen / *delta) * (phi / phip);
		}
		/* ENDUNTIL (DONE) */

		if (loop_counter++ > loop_counter_maximum) {
			imsl_e1sti (1, loop_counter_maximum);
			imsl_ermes (IMSL_FATAL, 
				IMSL_EXCEEDED_MAX_LOOP_COUNTER);
			done = 1;
		}

		if (!(done))
			goto L_10022;
	}
	deltap = *delta;

	
	imsl_e1pop ("l_u7lsf");
	return;
}				/* end of function */
/* Structured by FOR_STRUCT, v0.2, on 10/02/90 at 10:55:24
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  U8LSF/DU8LSF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Unconstrained nonlinear least squares solver using finite
                difference Jacobian.

    Usage:      CALL U8LSF (ITER, N, CJNORM, XSCALE)

    Arguments:
       ITER   - Current iteration of the QR factorization.  (Input)
       N      - Length of XSCALE and CJNORM.  (Input)
       CJNORM - Vector of length N containing the norms of the columns of
                the QR matrix.  (Input)
       XSCALE - Vector of length N containing the diagonal scaling matrix
                for the variables.  (Output)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_u8lsf(Mint iter, Mint n, Mfloat cjnorm[], Mfloat xscale[])
#else
static void l_u8lsf(iter, n, cjnorm, xscale)
	Mint             iter, n;
	Mfloat           cjnorm[], xscale[];
#endif
{
	Mint             i;


	if (iter == 1) {
		scopy(n, cjnorm, 1, xscale, 1);
	} else {
		for (i = 1; i <= n; i++) {
			xscale[i - 1] = imsl_f_max(xscale[i - 1], cjnorm[i - 1]);
		}
	}

	for (i = 1; i <= n; i++) {
		if (xscale[i - 1] <= 1.0e-6)
			xscale[i - 1] = 1.0e0;
	}

	return;
}				/* end of function */
/* Structured by FOR_STRUCT, v0.2, on 10/02/90 at 10:58:12
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  U9LSF/DU9LSF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Updating the model trust region for the nonlinear least
                squares problem.

    Usage:      CALL U9LSF (FCN, M, N, XC, FCNORM, GC, A, LDA, IPVT, SC,
                            XSCALE, GAUSS, STEPMX, DELTA, ICODE, XPPREV,
                            FPPREV, XP, FC, FP, FPNORM, MXTAKE, NFCN)

    Arguments:
       FCN    - User-supplied SUBROUTINE to evaluate the function to be
                minimized.  The usage is
                CALL FCN (M, N, X, F), where
                M      - Length of F.  (Input)
                N      - Length of X.  (Input)
                X      - The point at which the function is evaluated.
                         (Input)
                         X should not be changed by FCN.
                F      - The computed function at the point X.
                         (Output)
                FCN must be declared EXTERNAL in the calling program.
       M      - Number of functions.  (Input)
       N      - Number of variables.  (Input)
       XC     - Real vector of length N containing the current iterate.
                   (Input)
       FCNORM - Real scalar containing the 2-norm of F(XC).  (Input)
       GC     - Real vector of length N containing (JC**T)*FSCALE**2*FC.
                   (Input)
       A      - Real M by N matrix containing the upper triangular matrix
                R from the QR factorization of the current Jacobian in
                the upper triangular part and diagonal.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       IPVT   - Integer vector of length N containing the permutation
                matrix from QR factorization of the Jacobian.  (Output)
       SC     - Real vector of length N containing current step.  (Input)
       XSCALE - Real vector of length N containing the diagonal scaling
                matrix for X.  (Input)
       GAUSS  - Logical variable equal to .TRUE. when the Gauss-Newton
                step is taken.  (Input)
       STEPMX - Maximum allowable step size.  (Input)
       DELTA  - Trust region radius with value retained between calls.
                   (Input/Output)
       ICODE  - Return code.  (Output)
                ICODE = 0 means XP accepted as next iterate, DELTA is
                          trust region for next iteration.
                ICODE = 1 means the algorithm was unable to find a
                          satisfactory XP sufficiently distinct from XC.
                ICODE = 2 means both the scaled actual and predicted
                          function reductions are smaller than RFTOL.
                ICODE = 3 means that False Convergence is detected.
                ICODE = 4 means FPNORM is too large, current iteration is
                          continued with a new, reduced trust region.
                ICODE = 5 means FPNORM is sufficiently small, but the
                          chance of taking a longer successful step seems
                          good that the current iteration is to be
                          continued with a new, doubled trust region.
       XPPREV - Real vector of length N containing the value of XP at the
                previous call within this iteration.  (Input/Output)
       FPPREV - Real vector of length M containing F(XPPREV).
                   (Input/Output)
       XP     - Real vector of length N containing the new iterate.
                   (Output)
       FP     - Real vector of length M containing the functions at XP.
                   (Output)
       FPNORM - Real scalar containing the 2-norm of F(XP).  (Output)
       MXTAKE - Logical variable.  (Output)
                MXTAKE = .TRUE. if a step of maximum length was taken.
                MXTAKE = .FALSE. otherwise.
       NFCN   - Number of function evaluations used.  (Input/Output)

    Remark:
       This is based on the idea of NL2SOL and Dennis-Schnabel.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#define	ALPHA	1.0e-4

#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_u9lsf(void (*fcn) (Mint, Mint, Mfloat*, Mfloat*), 
                    Mint m, Mint n, Mfloat xc[], Mfloat fcnorm, 
                    Mfloat gc[], Mfloat *a, Mint lda, Mint ipvt[], 
                    Mfloat sc[], Mfloat xscale[], Mint gauss, Mfloat 
                    stepmx, Mfloat *delta, Mint *icode, Mfloat 
                    xpprev[], Mfloat fpprev[], Mfloat xp[], Mfloat fc[],
	            Mfloat fp[], Mfloat *fpnorm, Mint *mxtake, Mint *nfcn)
#else
static void l_u9lsf(void (*fcn) (Mint, Mint, Mfloat[], Mfloat[]), 
                    Mint m, Mint n, Mfloat xc[], Mfloat fcnorm, 
                    Mfloat gc[], Mfloat *a, Mint lda, Mint ipvt[], 
                    Mfloat sc[], Mfloat xscale[], Mint gauss, Mfloat 
                    stepmx, Mfloat *delta, Mint *icode, Mfloat 
                    xpprev[], Mfloat fpprev[], Mfloat xp[], Mfloat fc[],
	            Mfloat fp[], Mfloat *fpnorm, Mint *mxtake, Mint *nfcn)
#endif
#else
static void l_u9lsf(fcn, m, n, xc, fcnorm, gc, a, lda, ipvt, sc, xscale,
                    gauss, stepmx, delta, icode, xpprev, fpprev, xp, fc,
	            fp, fpnorm, mxtake, nfcn)
	void            (*fcn) ();
	Mint            m, n;
	Mfloat           xc[], fcnorm, gc[], *a;
	Mint            lda, ipvt[];
	Mfloat           sc[], xscale[];
	Mint            gauss;
	Mfloat          stepmx, *delta;
	Mint            *icode;
	Mfloat           xpprev[], fpprev[], xp[], fc[], fp[], *fpnorm;
	Mint            *mxtake;
	Mint            *nfcn;
#endif
{
#define A(I_,J_)	(a+(I_)*(lda)+(J_))
	Mint            ltemp;
	Mint             i, l;
	Mfloat           actred, big, prered, rellen, slope, small, stplen,
	                temp;
	static Mfloat    fpnrmp;


	*mxtake = 0;
	l_u11nf(n, xscale, 1, sc, xp);
	stplen = imsl_snrm2(n, xp, 1);
	/*
	 * COMPUTE NEW TRIAL POINT AND NEW FUNCTION VALUES
	 */
	for (i = 1; i <= n; i++) {
		xp[i - 1] = xc[i - 1] + sc[i - 1];
	}
	imsl_e1usr("ON");
	(*fcn) (m, n, xp, fp);
	imsl_e1usr("OFF");
	*nfcn += 1;
	*fpnorm = 0.5e0 * imsl_fi_power(imsl_snrm2(m, fp, 1), 2);
	actred = *fpnorm - fcnorm;
	slope = imsl_sdot(n, gc, 1, sc, 1);

	if (*icode != 5)
		fpnrmp = 0.0e0;
	if (*icode == 5 && ((*fpnorm >= fpnrmp) || (actred > ALPHA * slope))) {
		/*
		 * INTERNAL DOUBLING NO GOOD; RESET XP TO XPPREV AND
		 * TERMINATE
		 */
		*icode = 0;
		scopy(n, xpprev, 1, xp, 1);
		scopy(m, fpprev, 1, fp, 1);
		*fpnorm = fpnrmp;
		*delta *= 0.5e0;
	} else if (actred >= ALPHA * slope) {
		/*
		 * FPNORM IS TOO LARGE; THE STEP IS UNACCEPTABLE
		 */
		rellen = 0.0e0;
		for (i = 1; i <= n; i++) {
			temp = fabs(sc[i - 1]) / imsl_f_max(fabs(xp[i - 1]), 1.0e0 /
							    xscale[i - 1]);
			rellen = imsl_f_max(rellen, temp);
		}
		if (rellen <= u14sf.steptl) {
			/*
			 * XP - XC IS TOO SMALL, TERMINATE THE GLOBAL STEP
			 */
			*icode = 1;
			scopy(n, xc, 1, xp, 1);
			scopy(m, fc, 1, fp, 1);
		} else {
			/*
			 * QUADRATIC INTERPOLATION STEP; REDUCE DELTA AND
			 * CONTINUE
			 */
			*icode = 4;
			small = imsl_amach(1);
			big = imsl_amach(2);
			if (big * small < 1.0e0)
				small = 1.0e0 / big;

			if (fabs(actred - slope) > small) {
				temp = (-slope * stplen) / (2.0e0 * (actred - slope));
			} else {
				temp = (-slope * stplen) / 2.0e0;
			}
			if (temp < 0.1e0 ** delta) {
				*delta *= 0.1e0;
			} else if (temp > 0.5e0 ** delta) {
				*delta *= 0.5e0;
			} else {
				*delta = temp;
			}
		}
	} else {
		/*
		 * FPNORM IS SUFFICIENTLY SMALL; THE STEP IS ACCEPTABLE
		 * COMPUTE PREDICTED REDUCTION PRERED = G(t)*S +
		 * (1/2)*S(t)*H*S WITH H = P R**T * R * P**T.
		 */
		prered = slope;
		for (i = 1; i <= n; i++) {
			l = ipvt[i - 1];
			temp = imsl_sdot(n - i + 1, A(i - 1, i - 1), lda, &sc[l - 1],
					 0);
			prered += 0.5e0 * temp * temp;
		}

		ltemp = fabs(prered - actred) <= 0.1e0 * fabs(actred);
		if (((*icode != 4 && (ltemp || (actred <= slope))) && !gauss
		     ) && (*delta <= 0.99e0 * stepmx)) {

			/*
			 * IF ACTRED AND PRERED AGREE TO WITHIN RELATIVE
			 * ERROR 0.1 OR NEGATIVE CURVATURE IS INDICATED, AND
			 * A LONGER STEP IS POSSIBLE AND DELTA HAS NOT BEEN
			 * DECREASED THIS ITERATION, THEN DOUBLE TRUST REGION
			 * AND CONTINUE GLOBAL STEP
			 */
			*icode = 5;
			scopy(n, xp, 1, xpprev, 1);
			scopy(m, fp, 1, fpprev, 1);
			fpnrmp = *fpnorm;
			*delta = imsl_f_min(2.0e0 ** delta, stepmx);
		} else {
			/*
			 * ACCEPT XP AND CHOOSE NEW TRUST REGION FOR NEXT
			 * ITERATION
			 */
			*icode = 0;
			if (stplen > 0.99e0 * stepmx)
				*mxtake = 1;
			if (actred >= 0.1e0 * prered) {
				*delta *= 0.5e0;
			} else if (actred <= 0.75e0 * prered) {
				*delta = imsl_f_min(2.0e0 ** delta, stepmx);
			} else {
				/* DELTA = DELTA !* KEEP SAME DELTA */
			}
		}
		/*
		 * CHECK RELATIVE FUNCTION CONVERGENCE AND FALSE CONVERGENCE
		 */
		if (actred <= 2.0e0 * prered) {
			if (fabs(actred) <= u14sf.rftol * fabs(fcnorm) && fabs(prered) <=
			    u14sf.rftol * fabs(fcnorm)) {
				*icode = 2;
			}
		} else {
			rellen = 0.0e0;
			for (i = 1; i <= n; i++) {
				temp = fabs(sc[i - 1]) / imsl_f_max(fabs(xp[i - 1]),
						     1.0e0 / xscale[i - 1]);
				rellen = imsl_f_max(rellen, temp);
			}
			if (rellen < u14sf.falstl)
				*icode = 3;
		}
	}

	return;
}				/* end of function */
/* Structured by FOR_STRUCT, v0.2, on 10/02/90 at 10:57:20
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  U10SF/DU10SF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Compute QTF = TRANS(Q)*F

    Usage:      CALL U10SF (M, N, FJAC, LDFJAC, QRAUX, F, QTF)

    Arguments:
       M      - Number of rows in FJAC.  (Input)
       N      - Number of columns in FJAC.  (Input)
       FJAC   - M by N matrix containing information about the
                QR factorization as output from routine LQRRR.  (Input)
       LDFJAC - Leading dimension of FJAC exactly as specified in the
                dimension statement of the calling program.  (Input)
       QRAUX  - Vector of length N containing information about the
                orthogonal part of the decomposition as output from
                routine LQRRR.  (Input)
       F      - Vector of length N to be multiplied by Q**T.  (Input)
       QTF    - Vector of length N containing the result.  (Output)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_u10sf(Mint m, Mint n, Mfloat *fjac, Mint ldfjac, 
                    Mfloat qraux[], Mfloat f[], Mfloat qtf[])
#else
static void l_u10sf(m, n, fjac, ldfjac, qraux, f, qtf)
	Mint            m, n;
	Mfloat          *fjac;
	Mint            ldfjac;
	Mfloat           qraux[], f[], qtf[];
#endif
{
#define FJAC(I_,J_)	(fjac+(I_)*(ldfjac)+(J_))
	Mint             j;
	Mfloat           sum, temp;


	for (j = 1; j <= n; j++) {
		if (qraux[j - 1] != 0.0e0) {
			sum = qraux[j - 1] * f[j - 1];
			if (j < m)
				sum += imsl_sdot(m - j, FJAC(j - 1, j), 1, &f[j], 1);
			temp = -sum / qraux[j - 1];
			f[j - 1] += qraux[j - 1] * temp;
			if (j < m)
				saxpy(m - j, temp, FJAC(j - 1, j), 1, &f[j], 1);
		}
		qtf[j - 1] = f[j - 1];
	}

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
static void l_u11nf(Mint n, Mfloat y[], Mint k, Mfloat z[], Mfloat x[])
#else
static void l_u11nf(n, y, k, z, x)
	Mint            n;
	Mfloat           y[];
	Mint            k;
	Mfloat           z[], x[];
#endif
{
	Mint             i;


	if (k < 0) {
		if (k == -1) {
			for (i = 1; i <= n; i++) {
				x[i - 1] = z[i - 1] / y[i - 1];
			}
		} else {
			for (i = 1; i <= n; i++) {
				x[i - 1] = z[i - 1] / (imsl_fi_power(y[i - 1], -k));
			}
		}
	} else {
		if (k == 1) {
			for (i = 1; i <= n; i++) {
				x[i - 1] = z[i - 1] * y[i - 1];
			}
		} else {
			for (i = 1; i <= n; i++) {
				x[i - 1] = imsl_fi_power(y[i - 1], k) * z[i - 1];
			}
		}
	}

	return;
}				/* end of function */
/* Structured by FOR_STRUCT, v0.2, on 10/02/90 at 11:17:40
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  U11SF/DU11SF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Solve the systems A*X = B, and D*X = 0 in the least
                squares sense.

    Usage:      CALL U11SF (N, R, LDR, IPVT, DIAG, QTB, X, SDIAG, WA)

    Arguments:
       N      - Order of R.  (Input)
       R      - N by N array containing the upper triangular matrix R.
                (Input/Output)
                On output the full triangle is unaltered, and the strict
                lower triangle contains the transpose of the strict upper
                triangluar matrix S.
       LDR    - Leading dimension of R exactly as specified in the
                dimension statement of the calling program.  (Input)
       IPVT   - Vector of length N which defines the permutation matrix P
                such that A*P = Q*R.  (Input)
                Column J of P is column IPVT(J) of the identity matrix.
       DIAG   - Vector of length N containing the diagonal elements of
                the matrix D.  (Input)
       QTB    - Vector of length N which must contain the first N
                elements of the vectorS TRANS(Q)*B, where TRANS(Q)
                represents the transpose of Q.  (Input)
       X      - Vector of length N containing the least squares solution
                of the systems A*X = B, D*X = 0.  (Output)
       SDIAG  - Vector of length N containing the diagonal elements of
                the upper triangular matrix S.  (Output)
       WA     - Work vector of length N.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_u11sf(Mint n, Mfloat *r, Mint ldr, Mint ipvt[], 
                    Mfloat diag[], Mfloat qtb[], Mfloat x[], Mfloat 
                    sdiag[], Mfloat wa[])
#else
static void l_u11sf(n, r, ldr, ipvt, diag, qtb, x, sdiag, wa)
	Mint            n;
	Mfloat          *r;
	Mint            ldr, ipvt[];
	Mfloat           diag[], qtb[], x[], sdiag[], wa[];
#endif
{
#define R(I_,J_)	(r+(I_)*(ldr)+(J_))
	Mint             i, j, jp1, k, kp1, l, nsing;
	Mfloat           cos1, cot1, qtbpj, sin1, sum, tan1, temp;

	/*
	 * COPY R AND (Q TRANSPOSE)*B TO PRESERVE INPUT AND INITIALIZE S. IN
	 * PARTICULAR, SAVE THE DIAGONAL ELEMENTS OF R IN X.
	 */
	imsl_csfrg(&n, r, &ldr);
	scopy(n, R(0, 0), ldr + 1, x, 1);
	scopy(n, qtb, 1, wa, 1);
	/*
	 * ELIMINATE THE DIAGONAL MATRIX D USING A GIVENS ROTATION.
	 */
	for (j = 1; j <= n; j++) {
		/*
		 * PREPARE THE ROW OF D TO BE ELIMINATED, LOCATING THE
		 * DIAGONAL ELEMENT USING P FROM THE QR FACTORIZATION.
		 */
		l = ipvt[j - 1];
		if (diag[l - 1] != 0.0e0) {
			sset(n - j + 1, 0.0e0, &sdiag[j - 1], 1);
			sdiag[j - 1] = diag[l - 1];
			/*
			 * THE TRANSFORMATIONS TO ELIMINATE THE ROW OF D
			 * MODIFY ONLY A SINGLE ELEMENT OF (Q TRANSPOSE)*B
			 * BEYOND THE FIRST N, WHICH IS INITIALLY ZERO.
			 */
			qtbpj = 0.0e0;
			for (k = j; k <= n; k++) {
				/*
				 * DETERMINE A GIVENS ROTATION WHICH
				 * ELIMINATES THE APPROPRIATE ELEMENT IN THE
				 * CURRENT ROW OF D.
				 */
				if (sdiag[k - 1] != 0.0e0) {
					if (fabs(*R(k - 1, k - 1)) < fabs(sdiag[k - 1])) {
						cot1 = *R(k - 1, k - 1) / sdiag[k - 1];
						sin1 = 0.5e0 / sqrt(0.25e0 + 0.25e0 * cot1 * cot1);
						cos1 = sin1 * cot1;
					} else {
						tan1 = sdiag[k - 1] / *R(k - 1, k - 1);
						cos1 = 0.5e0 / sqrt(0.25e0 + 0.25e0 * tan1 * tan1);
						sin1 = cos1 * tan1;
					}
					/*
					 * COMPUTE THE MODIFIED DIAGONAL
					 * ELEMENT OF R AND THE MODIFIED
					 * ELEMENT OF ((Q TRANSPOSE)*B,0).
					 */
					*R(k - 1, k - 1) = cos1 ** R(k - 1, k - 1) + sin1 *
						sdiag[k - 1];
					temp = cos1 * wa[k - 1] + sin1 * qtbpj;
					qtbpj = -sin1 * wa[k - 1] + cos1 * qtbpj;
					wa[k - 1] = temp;
					/*
					 * ACCUMULATE THE TRANFORMATION IN
					 * THE ROW OF S.
					 */
					kp1 = k + 1;
					if (n >= kp1) {
						for (i = kp1; i <= n; i++) {
							temp = cos1 ** R(k - 1, i - 1) + sin1 * sdiag[i - 1];
							sdiag[i - 1] = -sin1 ** R(k - 1, i - 1) +
								cos1 * sdiag[i - 1];
							*R(k - 1, i - 1) = temp;
						}
					}
				}
			}
		}
		/*
		 * STORE THE DIAGONAL ELEMENT OF S AND RESTORE THE
		 * CORRESPONDING DIAGONAL ELEMENT OF R.
		 */
		sdiag[j - 1] = *R(j - 1, j - 1);
		*R(j - 1, j - 1) = x[j - 1];
	}
	/*
	 * SOLVE THE TRIANGULAR SYSTEM FOR Z. IF THE SYSTEM IS SINGULAR, THEN
	 * OBTAIN A LEAST SQUARES SOLUTION.
	 */
	nsing = n;
	for (j = 1; j <= n; j++) {
		if (sdiag[j - 1] == 0.0e0 && nsing == n)
			nsing = j - 1;
		if (nsing < n)
			wa[j - 1] = 0.0e0;
	}
	if (nsing >= 1) {
		for (k = 1; k <= nsing; k++) {
			j = nsing - k + 1;
			sum = 0.0e0;
			jp1 = j + 1;
			if (nsing >= jp1)
				sum = imsl_sdot(nsing - j, R(j - 1, jp1 - 1), 1, &wa[jp1 - 1],
						1);
			wa[j - 1] = (wa[j - 1] - sum) / sdiag[j - 1];
		}
	}
	/*
	 * PERMUTE THE COMPONENTS OF Z BACK TO COMPONENTS OF X.
	 */
	for (j = 1; j <= n; j++) {
		l = ipvt[j - 1];
		x[l - 1] = wa[j - 1];
	}

	return;
}				/* end of function */
/* Structured by FOR_STRUCT, v0.2, on 10/02/90 at 11:18:30
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  U12SF/DU12SF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Solve TRANS(R)*s=y for s.

    Usage:      CALL U12SF (N, H, LDH, Y, SNWTN)

    Arguments:
       N      - Length of the vector SNWTN.  (Input)
       H      - N by N matrix containing the Cholesky factor of the
                Hessian in the upper triangle and diagonal.  (Input)
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
static void l_u12sf(Mint n, Mfloat *h, Mint ldh, Mfloat y[], 
                    Mfloat snwtn[])
#else
static void l_u12sf(n, h, ldh, y, snwtn)
	Mint            n;
	Mfloat          *h;
	Mint            ldh;
	Mfloat           y[], snwtn[];
#endif
{
#define H(I_,J_)	(h+(I_)*(ldh)+(J_))
	Mint             i;
	Mfloat           sum;


	snwtn[n - 1] = y[n - 1] / *H(n - 1, n - 1);

	for (i = n - 1; i >= 1; i--) {
		sum = imsl_sdot(n - i, H(i, i - 1), ldh, &snwtn[i], 1);
		snwtn[i - 1] = (y[i - 1] - sum) / *H(i - 1, i - 1);
	}

	return;
}				/* end of function */
/* Structured by FOR_STRUCT, v0.2, on 10/02/90 at 11:19:21
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  U13SF/DU13SF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    September 16, 1985

    Purpose:    Check validity of input to nonlinear least squares solver

    Usage:      CALL U13SF (ICODE)

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
static void l_u13sf(Mint icode)
#else
static void l_u13sf(icode)
	Mint            icode;
#endif
{


	if (icode == 2) {

/*		(6, 2, "The values for the scaling matrix, XSCALE, must be positive while at least one entry is less than or equal to zero. The algorithm will use the identity scaling matrix for XSCALE."); */
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEED_POSITIVE_XSCALE_ELEM);
	} else if (icode == 3) {

/*		(6, 3, "The values for the diagonal matrix, FSCALE, must be positive while at least one entry is less than or equal to zero. The algorithm will use the identity scaling matrix for FSCALE."); */
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEED_POSITIVE_FSCALE_ELEM);
	} else if (icode == 4) {

/*		(6, 4, "The estimate of the number of good digits in the functions must be positive while NDIGIT = %(i1) is given.  The algorithm will assume that the function is accurate to the precision of the arithmetic."); */
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEED_POSITIVE_NDIGIT);
	} else if (icode == 5) {

/*		(6, 5, "The maximum number of iterations must be positive while MAXITN = %(i1) is given.  The algorithm will use MAXITN = 100."); */
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEED_POSITIVE_MXITER);
	} else if (icode == 6) {

/*		(6, 6, "The maximum number of function evaluations must be positive while MAXFCN = %(i1) is given.  The algorithm will use MAXFCN = 400.");
*/
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEED_POSITIVE_MAXFCN);
	} else if (icode == 7) {

/*		(6, 7, "The maximum number of Jacobian evaluations must be positive while MAXJAC = %(i1) is given.  The algorithm will use MAXJAC = 100.");
*/
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEED_POSITIVE_MAXJAC);
	} else if (icode == 8) {

/*		(6, 8, "The gradient tolerance must be nonnegative while GRADTOL = %(r1) is given.  The algorithm will use GRADTOL = %(r2).");
*/
                imsl_ermes(IMSL_WARNING_IMMEDIATE, IMSL_NEED_POSITIVE_GRADTL);
	}
	return;
}				/* end of function */
/* Structured by FOR_STRUCT, v0.2, on 10/10/90 at 09:19:58
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  F2JAC/DF2JAC  (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 1, 1985

    Purpose:    Approximate the Jacobian of M functions in N unknowns
                using forward differences.

    Usage:      CALL F2JAC (FCN, M, N, XC, XSCALE, FC, EPSFCN, FJAC,
                            LDFJAC, WK)

    Arguments:  See FDJAC/DFDJAC.

    Remarks:    See FDJAC/DFDJAC.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_f2jac(void (*fcn) (Mint, Mint, Mfloat*, Mfloat*), Mint m, 
                    Mint n, Mfloat xc[], Mfloat xscale[], Mfloat fc[], 
                    Mfloat epsfcn, Mfloat *fjac, Mint ldfjac, Mfloat work[])
#else
static void l_f2jac(void (*fcn) (Mint, Mint, Mfloat[], Mfloat[]), Mint m, 
                    Mint n, Mfloat xc[], Mfloat xscale[], Mfloat fc[], 
                    Mfloat epsfcn, Mfloat *fjac, Mint ldfjac, Mfloat work[])
#endif
#else
static void l_f2jac(fcn, m, n, xc, xscale, fc, epsfcn, fjac, ldfjac,
	            work)
	void            (*fcn) ();
	Mint            m, n;
	Mfloat           xc[], xscale[], fc[], epsfcn, *fjac;
	Mint            ldfjac;
	Mfloat           work[];
#endif
{
#define FJAC(I_,J_)	(fjac+(I_)*(ldfjac)+(J_))
	Mint             i, j;
	Mfloat           eps, stepsz, xtempj;


	imsl_e1psh("F2JAC ");

	if (m <= 0) {
		imsl_e1sti(1, m);
/*		(5, 1, "The number of functions must be positive while M = %(i1) is given."); */
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_POSITIVE_NUM_FCNS);

	} else if (n <= 0) {
		imsl_e1sti(1, n);
/*		(5, 2, "The number of variables must be positive while N = %(i1) is given."); */
                imsl_ermes(IMSL_TERMINAL, IMSL_N_MUST_BE_POSITIVE);

	} else if (epsfcn >= 0.1e0 || epsfcn < 0.0e0) {
		imsl_e1str(1, epsfcn);
/*		(5, 4, "The estimate for the relative noise in the function must be between 0.0 and 0.1 while EPSFCN = %(r1) is given."); */
                imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_EPSFCN_VALUE);

	} else {
		for (i = 0; i < n; i++) {
			if (xscale[i] <= 0.0e0) {
				imsl_e1sti(1, i);
				imsl_e1str(1, xscale[i]);
/*				(5, 4, "The values for the diagonal scaling matrix must be positive while XSCALE(%(i1)) = %(r1) is given."); */
                                imsl_ermes(IMSL_TERMINAL, IMSL_POS_XSCALE_ELMNTS_NEEDED);
				goto L_9000;
			}
		}
	}

	if (imsl_n1rcd(0) == 0) {
		eps = sqrt(imsl_f_max(epsfcn, imsl_amach(4)));
		for (j = 1; j <= n; j++) {
			stepsz = (eps) * imsl_f_max(fabs(xc[j - 1]), 1.0e0 / xscale[j - 1]);
			if (xc[j - 1] < 0.0)
				stepsz = -stepsz;
			xtempj = xc[j - 1];
			xc[j - 1] = xtempj + stepsz;
                        imsl_e1usr("ON"); 
			(*fcn) (m, n, xc, work);
                        imsl_e1usr("OFF"); 
			xc[j - 1] = xtempj;
			for (i = 1; i <= m; i++) {
				*FJAC(j - 1, i - 1) = (work[i - 1] - fc[i - 1]) / stepsz;
			}
		}
	}
L_9000:
	imsl_e1pop("F2JAC ");
	return;
}				/* end of function */
