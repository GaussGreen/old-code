#define ADR(t,x)    ( t = x, &t )

#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK  l_ode_adams_gear_mgr(Mint ido, Mchar **state,
                        va_list argptr);
static void 	l_init_state(void);
static void 	l_i3prk(Mint neq, Mfloat v[], Mfloat y[], Mfloat ymax[],
                        Mfloat *enorm);
static void 	l_i4pag(void (*fcn) (Mint, Mfloat, Mfloat*, Mfloat*),
                        void (*fcnj) (Mint, Mfloat, Mfloat[], Mfloat*, Mfloat*),
 			Mint *n, Mfloat *a,
                        Mint *lda, Mfloat *y, Mfloat ymax[], Mfloat error[],
                        Mfloat save1[], Mfloat save2[], Mfloat pw[],
                        Mint ipvt[], Mint *kflag, Mint *ido,
                        void (*vnorm) ());
static void 	l_i5pag(void (*fcn) (Mint, Mfloat, Mfloat*, Mfloat*),
			Mint *n, Mfloat *pw, Mfloat save1[],
                        Mfloat save2[], Mfloat *y, Mfloat ymax[]);
static void 	l_i6pag(void (*fcn) (Mint, Mfloat, Mfloat*, Mfloat*),
			Mint *n, Mfloat *pw, Mint *ldpw,
                        Mfloat save1[], Mfloat save2[], Mfloat *y,
                        Mfloat ymax[], Mfloat wk[]);
static void 	l_i7pag(void (*fcn) (Mint, Mfloat, Mfloat*, Mfloat*),
			Mint *n, Mfloat pw[], Mfloat save1[],
                        Mfloat save2[], Mfloat *y, Mfloat ymax[],
                        Mfloat *hl0, Mint *ier);
static void 	l_i8pag(void (*fcn) (Mint, Mfloat, Mfloat*, Mfloat*), 
			Mint *n, Mfloat pw[], Mfloat save1[],
                        Mfloat *y, Mfloat *hl0, Mint *ier);
static void 	l_i9pag(Mint *n, Mfloat *tout, Mfloat *y, Mfloat y0[]);
static void	l_i10ag(Mint *meth, Mint *nq, Mfloat el[], Mfloat tq[],
                        Mint *maxder);
static void 	l_i12ag(Mint *n, Mfloat *a, Mint *lda, Mint *nlc, Mint *nuc,
                        Mfloat b[], Mfloat x[], Mfloat imsl_fac[],
                        Mint ipvt[], Mfloat wk[], Mint *mtype, Mint *istop);
static void	l_i13ag(Mint *icode, Mint *i1, Mfloat *r1, Mfloat *r2,
                        Mfloat *r3);
static void 	l_ssbmv(Mchar *uplo, unsigned uplo_s, Mint *n, Mint *ncoda,
                        Mfloat *alpha, Mfloat *a, Mint *lda, Mfloat x[],
                        Mint *incx, Mfloat *imsl_beta, Mfloat y[],
                        Mint *incy);
static void 	l_sgbmv(Mchar *trans, unsigned trans_s, Mint *m, Mint *n,
                        Mint *nlca, Mint *nuca, Mfloat *alpha,
                        Mfloat *a, Mint *lda, Mfloat x[], Mint *incx,
                        Mfloat *imsl_beta, Mfloat y[], Mint *incy);
#else
static VA_LIST_HACK 	l_ode_adams_gear_mgr();
static void	l_init_state();
static void 	l_i3prk();
static void 	l_i4pag();
static void	l_i5pag();
static void 	l_i6pag();
static void	l_i7pag();
static void	l_i8pag();
static void	l_i9pag();
static void	l_i10ag();
static void	l_i12ag();
static void	l_i13ag();
static void	l_ssbmv();
static void	l_sgbmv();
#endif


typedef struct {
#ifdef ANSI
#if defined(COMPUTER_HP98C) || defined(COMPUTER_HP97C) || defined(COMPUTER_PMXUX)
#ifdef DOUBLE
	void (*fcnj) (int, double, double[], double*, double*);
#else
	void (*fcnj) (int, float, float[], float*, float*);
#endif /*DOUBLE*/
#else
        void (*fcnj) (Mint, Mfloat, Mfloat[], Mfloat*, Mfloat*);
#endif /* COMPUTER_HP9?C */
#else
	void (*fcnj) ();
#endif /* ANSI */
	Mint method;
	Mint maxord;
	Mint miter;
	Mfloat tol;
	Mfloat hinit;
	Mfloat hmin;
	Mfloat hmax;
	Mint max_steps;
	Mint max_fcn_evals;
	Mfloat scale;
	Mint norm;
	Mfloat floor;
	Mint *nstep;
	Mint *nfcn;
	Mint *nfcnj;
	Mint ido;
	Mfloat *wk;
	Mint *iwk;
	Mfloat *iymax;
	Mfloat *ierror;
	Mfloat *isave1;
	Mfloat *isave2;
	Mfloat *iytemp;
	Mfloat *ipw;
	Mint istat;
} Env;

static Env *lv_state;

static struct t_i11ag {
	Mfloat           t, h, hmin, hmax, eps, uround, hused, el[13], oldl0,
	                epsj;
	Mint             intrp1, intrp2, iatype, mf, jstart, nfcn, nje,
	                meth, miter, mtype, nlc, nuc, nstep, mxfcn, maxord,
	                ido1, istop;
}               i11ag;
static struct t_i5prk {
	Mfloat           floor;
	Mint             inorm;
}               i5prk;


#ifdef ANSI
void imsl_f_ode_adams_gear_mgr(Mint ido, Mchar **state, ...)
#else
void imsl_f_ode_adams_gear_mgr(ido, state, va_alist)
    Mchar        **state;
    Mint        ido;
    va_dcl
#endif
{
    va_list     argptr;

    VA_START(argptr, state);

    E1PSH("imsl_f_ode_adams_gear_mgr", "imsl_d_ode_adams_gear_mgr");
    IMSL_CALL(l_ode_adams_gear_mgr(ido, state, argptr));
    va_end(argptr);
    E1POP("imsl_f_ode_adams_gear_mgr", "imsl_d_ode_adams_gear_mgr");
}


#ifdef ANSI
static VA_LIST_HACK l_ode_adams_gear_mgr(Mint ido, Mchar **state, va_list argptr)
#else
static VA_LIST_HACK l_ode_adams_gear_mgr(ido, state, argptr)
    Mchar        **state;
    Mint        ido;
    va_list     argptr;
#endif
{
    Mint            code;
    Mint            arg_number  = 2;

    if (ido == 1) {
        *state = (Mchar*)imsl_malloc(sizeof(Env));
        lv_state = (Env *)*state;
        l_init_state();
        lv_state->ido = 1;
	if (lv_state->istat == 0) 
		lv_state->istat = 1;
	else {
		Mint _l0 = 19, _l1 = 0; 
		Mfloat _f0 = F_ZERO;
                l_i13ag(&_l0, &_l1, &_f0, &_f0,
                          &_f0);
		goto RETURN;
	}
    }
    else if (ido == 3) {
        lv_state = (Env *)*state;
	if (lv_state!=NULL) lv_state->istat = 0;
        if (lv_state!=NULL && lv_state->wk!=NULL) imsl_free(lv_state->wk);
        if (lv_state!=NULL && lv_state->iwk!=NULL) imsl_free(lv_state->iwk);
	if (lv_state!=NULL)imsl_free(lv_state);
	goto RETURN;
#if 0
	if (lv_state != NULL) {
		if (lv_state->wk != NULL) imsl_free(lv_state->wk);
		if (lv_state->iwk != NULL) imsl_free(lv_state->iwk);
		if (lv_state->iymax != NULL) imsl_free(lv_state->iymax);
		if (lv_state->ierror != NULL) imsl_free(lv_state->ierror);
		if (lv_state->isave1 != NULL) imsl_free(lv_state->isave1);
		if (lv_state->isave2 != NULL) imsl_free(lv_state->isave2);
		if (lv_state->iytemp != NULL) imsl_free(lv_state->iytemp);
		if (lv_state->ipw != NULL) imsl_free(lv_state->ipw);
		imsl_free(lv_state);
	}
#endif
    }

    code = 1;
    while (code > 0) {
        code = va_arg(argptr, Mint);
        arg_number++;
        switch (code) {
	    case IMSL_JACOBIAN:
#ifdef ANSI
		lv_state->fcnj = (void (*) (Mint, Mfloat, Mfloat[], 
			Mfloat*, Mfloat*)) va_arg(argptr, void *);
#else
		lv_state->fcnj = (void (*) ()) va_arg(argptr, void *);
#endif
		arg_number++;
		break;
	    case IMSL_METHOD:
		lv_state->method = va_arg(argptr, Mint);
		arg_number++;
		break;
	    case IMSL_MAXORD:
		lv_state->maxord = va_arg(argptr, Mint);
		arg_number++;
		break;
	    case IMSL_MITER:
		lv_state->miter = va_arg(argptr, Mint);
		arg_number++;
		break;
	    case IMSL_TOL:
		lv_state->tol = (Mfloat) va_arg(argptr, Mdouble);
		arg_number++;
		break;
	    case IMSL_HINIT:
		lv_state->hinit = (Mfloat) va_arg(argptr, Mdouble);
		arg_number++;
		break;
	    case IMSL_HMIN:
		lv_state->hmin = (Mfloat) va_arg(argptr, Mdouble);
		arg_number++;
		break;
	    case IMSL_HMAX:
		lv_state->hmax = (Mfloat) va_arg(argptr, Mdouble);
		arg_number++;
		break;
	    case IMSL_TOL_ADR:
		lv_state->tol = *(va_arg(argptr, Mfloat *));
		arg_number++;
		break;
	    case IMSL_HINIT_ADR:
		lv_state->hinit = *(va_arg(argptr, Mfloat *));
		arg_number++;
		break;
	    case IMSL_HMIN_ADR:
		lv_state->hmin = *(va_arg(argptr, Mfloat *));
		arg_number++;
		break;
	    case IMSL_HMAX_ADR:
		lv_state->hmax = *(va_arg(argptr, Mfloat *));
		arg_number++;
		break;
	    case IMSL_MAX_NUMBER_STEPS:
		lv_state->max_steps = va_arg(argptr, Mint);
		arg_number++;
		break;
	    case IMSL_MAX_NUMBER_FCN_EVALS:
		lv_state->max_fcn_evals = va_arg(argptr, Mint);
		arg_number++;
		break;
	    case IMSL_SCALE:
		lv_state->scale = (Mfloat) va_arg(argptr, Mdouble);
		arg_number++;
		break;
	    case IMSL_SCALE_ADR:
		lv_state->scale = *(va_arg(argptr, Mfloat *));
		arg_number++;
		break;
	    case IMSL_NORM:
		lv_state->norm = va_arg(argptr, Mint);
		arg_number++;
		break;
	    case IMSL_FLOOR:
		lv_state->floor = (Mfloat) va_arg(argptr, Mdouble);
		arg_number++;
		break;
	    case IMSL_FLOOR_ADR:
		lv_state->floor = *(va_arg(argptr, Mfloat *));
		arg_number++;
		break;
	    case IMSL_NSTEP:
		lv_state->nstep = va_arg(argptr, Mint*);
		arg_number++;
		break;
            case IMSL_NFCN:
                lv_state->nfcn = va_arg(argptr, Mint*);
                arg_number++;
                break;
	    case IMSL_NFCNJ:
		lv_state->nfcnj = va_arg(argptr, Mint*);
		arg_number++;
		break;
	    case 0:
		break;
	    default:
		imsl_e1sti(1, code);
		imsl_e1sti(2, arg_number);
		imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
		break;
	}
    }

    if (lv_state->method == 1) lv_state->maxord = 12;
    if (lv_state->method == 2) lv_state->maxord = 5;

        /* Error return if HMIN > HMAX */
    if (lv_state->hmin > lv_state->hmax) {
            /* HMIN = %(r1) is greater than HMAX = %(r2). */
        imsl_e1str(1, lv_state->hmin);
        imsl_e1str(2, lv_state->hmax);
        imsl_ermes(IMSL_TERMINAL, IMSL_HMIN_GT_HMAX);
    }
RETURN:
    return (argptr);
}

#ifdef ANSI
static void l_init_state(void)
#else
static void l_init_state()
#endif
{
	lv_state->method = 2;
	lv_state->miter = 3;
	lv_state->tol = .001;
	lv_state->hmax = imsl_amach(2);
	lv_state->hmin = F_ZERO;
	lv_state->max_steps = 500;
	lv_state->max_fcn_evals = imsl_i_machine(5);
	lv_state->scale = 1;
	lv_state->norm = 0;
	lv_state->floor = F_ONE;
	lv_state->hinit = 0.0;
#if 0
        if (lv_state->istat != 1)
#endif
            lv_state->istat = 0;

	lv_state->nstep = NULL;
	lv_state->nfcn = NULL;
	lv_state->nfcnj = NULL;
	lv_state->wk = NULL;
	lv_state->iwk = NULL;
#ifdef ANSI
	lv_state->fcnj = (void (*) (Mint, Mfloat, Mfloat[], 
		Mfloat*, Mfloat*)) 0;
#else
	lv_state->fcnj = (void (*) ()) 0;
#endif
} 






/* Structured by FOR_STRUCT, v0.2, on 10/08/90 at 17:36:19
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  I10AG/DI10AG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Solve an initial value problem for ordinary differential
                equations using an Adams-Moulton or Gear method.

    Usage:      CALL I10AG (METH, NQ, EL, TQ, MAXDER)

    Arguments:
       METH   - Method indicator.  (Input)
                METH=1 selects Adams method
                METH=2 selects Gear's backward
                difference method.
       NQ     - Current order of the solver.  (Input)
       EL     - Vector of length MAXDER containing coefficients used in
                error estimates.  (Output)
       TQ     - Vector of length 4 containing coefficients used in
                error estimates.  (Output)
       MAXDER - Maximum order of the method.  (Output)

    Chapter:    MATH/LIBRARY Differential Equations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_i10ag(Mint *meth, Mint *nq, Mfloat el[], Mfloat tq[],
			Mint *maxder)
#else
static void l_i10ag(meth, nq, el, tq, maxder)
	Mint            *meth, *nq;
	Mfloat           el[], tq[];
	Mint            *maxder;
#endif
{
	static Mint             _i, _r;
	static Mfloat    elco[12][13][2], tesco[12][3][2];
	static Mint      maxnq[2] = {12, 5};
	static Mint      _aini = 1;
	if (_aini) {
		elco[0][0][0] = F_ONE;
		elco[1][0][0] = F_HALF;
		elco[1][2][0] = F_HALF;
		elco[2][0][0] = 4.16666666666666666666666666667e-01;
		elco[2][2][0] = 7.5e-01;
		elco[2][3][0] = 1.66666666666666666666666666667e-01;
		elco[3][0][0] = 3.75e-01;
		elco[3][2][0] = 9.16666666666666666666666666667e-01;
		elco[3][3][0] = 3.33333333333333333333333333333e-01;
		elco[3][4][0] = 4.16666666666666666666666666667e-02;
		elco[4][0][0] = 3.48611111111111111111111111111e-01;
		elco[4][2][0] = 1.04166666666666666666666666667e00;
		elco[4][3][0] = 4.86111111111111111111111111111e-01;
		elco[4][4][0] = 1.04166666666666666666666666667e-01;
		elco[4][5][0] = 8.33333333333333333333333333333e-03;
		elco[5][0][0] = 3.29861111111111111111111111111e-01;
		elco[5][2][0] = 1.14166666666666666666666666667e00;
		elco[5][3][0] = 6.25e-01;
		elco[5][4][0] = 1.77083333333333333333333333333e-01;
		elco[5][5][0] = 2.5e-02;
		elco[5][6][0] = 1.38888888888888888888888888889e-03;
		elco[6][0][0] = 3.15591931216931216931216931217e-01;
		elco[6][2][0] = 1.225e00;
		elco[6][3][0] = 7.51851851851851851851851851852e-01;
		elco[6][4][0] = 2.55208333333333333333333333333e-01;
		elco[6][5][0] = 4.86111111111111111111111111111e-02;
		elco[6][6][0] = 4.86111111111111111111111111111e-03;
		elco[6][7][0] = 1.98412698412698412698412698413e-04;
		elco[7][0][0] = 3.04224537037037037037037037037e-01;
		elco[7][2][0] = 1.29642857142857142857142857143e00;
		elco[7][3][0] = 8.68518518518518518518518518519e-01;
		elco[7][4][0] = 3.35763888888888888888888888889e-01;
		elco[7][5][0] = 7.77777777777777777777777777778e-02;
		elco[7][6][0] = 1.06481481481481481481481481482e-02;
		elco[7][7][0] = 7.93650793650793650793650793651e-04;
		elco[7][8][0] = 2.48015873015873015873015873016e-05;
		elco[8][0][0] = 2.94868000440917107583774250441e-01;
		elco[8][2][0] = 1.35892857142857142857142857143e00;
		elco[8][3][0] = 9.76554232804232804232804232804e-01;
		elco[8][4][0] = 4.171875e-01;
		elco[8][5][0] = 1.11354166666666666666666666667e-01;
		elco[8][6][0] = 1.875e-02;
		elco[8][7][0] = 1.93452380952380952380952380952e-03;
		elco[8][8][0] = 1.11607142857142857142857142857e-04;
		elco[8][9][0] = 2.7557319223985890652557319224e-06;
		elco[9][0][0] = 2.86975446428571428571428571429e-01;
		elco[9][2][0] = 1.41448412698412698412698412698e00;
		elco[9][3][0] = 1.07721560846560846560846560847e00;
		elco[9][4][0] = 4.985670194003527336860670194e-01;
		elco[9][5][0] = 1.484375e-01;
		elco[9][6][0] = 2.90605709876543209876543209877e-02;
		elco[9][7][0] = 3.72023809523809523809523809524e-03;
		elco[9][8][0] = 2.99685846560846560846560846561e-04;
		elco[9][9][0] = 1.3778659611992945326278659612e-05;
		elco[9][10][0] = 2.7557319223985890652557319224e-07;
		elco[10][0][0] = 2.80189596443936721714499492277e-01;
		elco[10][2][0] = 1.46448412698412698412698412698e00;
		elco[10][3][0] = 1.17151455026455026455026455027e00;
		elco[10][4][0] = 5.79358190035273368606701940035e-01;
		elco[10][5][0] = 1.88322861552028218694885361552e-01;
		elco[10][6][0] = 4.14303626543209876543209876543e-02;
		elco[10][7][0] = 6.21114417989417989417989417989e-03;
		elco[10][8][0] = 6.25206679894179894179894179894e-04;
		elco[10][9][0] = 4.04174015285126396237507348619e-05;
		elco[10][10][0] = 1.51565255731922398589065255732e-06;
		elco[10][11][0] = 2.50521083854417187750521083854e-08;
		elco[11][0][0] = 2.74265540031599059376837154615e-01;
		elco[11][2][0] = 1.50993867243867243867243867244e00;
		elco[11][3][0] = 1.26027116402116402116402116402e00;
		elco[11][4][0] = 6.59234182098765432098765432099e-01;
		elco[11][5][0] = 2.30458002645502645502645502646e-01;
		elco[11][6][0] = 5.56972461052322163433274544386e-02;
		elco[11][7][0] = 9.43948412698412698412698412698e-03;
		elco[11][8][0] = 1.11927496693121693121693121693e-03;
		elco[11][9][0] = 9.09391534391534391534391534392e-05;
		elco[11][10][0] = 4.8225308641975308641975308642e-06;
		elco[11][11][0] = 1.50312650312650312650312650313e-07;
		elco[11][12][0] = 2.08767569878680989792100903212e-09;
		tesco[0][0][0] = F_ZERO;
		tesco[0][1][0] = F_TWO;
		tesco[0][2][0] = 1.2e01;
		tesco[1][0][0] = F_ONE;
		tesco[1][1][0] = 1.2e01;
		tesco[1][2][0] = 2.4e01;
		tesco[2][0][0] = F_TWO;
		tesco[2][1][0] = 2.4e01;
		tesco[2][2][0] = 3.78947368421052631578947368421e01;
		tesco[3][0][0] = F_ONE;
		tesco[3][1][0] = 3.78947368421052631578947368421e01;
		tesco[3][2][0] = 5.33333333333333333333333333333e01;
		tesco[4][0][0] = 3.15789473684210526315789473684e-01;
		tesco[4][1][0] = 5.33333333333333333333333333333e01;
		tesco[4][2][0] = 7.00811123986095017381228273465e01;
		tesco[5][0][0] = 7.40740740740740740740740740741e-02;
		tesco[5][1][0] = 7.00811123986095017381228273465e01;
		tesco[5][2][0] = 8.79709090909090909090909090909e01;
		tesco[6][0][0] = 1.39049826187717265353418308227e-02;
		tesco[6][1][0] = 8.79709090909090909090909090909e01;
		tesco[6][2][0] = 1.06877153712484905604806644479e02;
		tesco[7][0][0] = 2.18181818181818181818181818182e-03;
		tesco[7][1][0] = 1.06877153712484905604806644479e02;
		tesco[7][2][0] = 1.267016986435292679946229989e02;
		tesco[8][0][0] = 2.9452478426059552911377492416e-04;
		tesco[8][1][0] = 1.267016986435292679946229989e02;
		tesco[8][2][0] = 1.47365474076838378148388230122e02;
		tesco[9][0][0] = 3.4915591557409961418271329062e-05;
		tesco[9][1][0] = 1.47365474076838378148388230122e02;
		tesco[9][2][0] = 1.68803254121173196317704988225e02;
		tesco[10][0][0] = 3.69181582884495696419523183527e-06;
		tesco[10][1][0] = 1.68803254121173196317704988225e02;
		tesco[10][2][0] = 1.90960201551200647713309316882e02;
		tesco[11][0][0] = 3.52406451504907700345270220862e-07;
		tesco[11][1][0] = 1.90960201551200647713309316882e02;
		tesco[11][2][0] = F_ZERO;
		elco[0][0][1] = F_ONE;
		elco[1][0][1] = 6.66666666666666666666666666667e-01;
		elco[1][2][1] = 3.33333333333333333333333333333e-01;
		elco[2][0][1] = 5.45454545454545454545454545455e-01;
		elco[2][2][1] = 5.45454545454545454545454545455e-01;
		elco[2][3][1] = 9.09090909090909090909090909091e-02;
		elco[3][0][1] = 4.8e-01;
		elco[3][2][1] = 7.0e-01;
		elco[3][3][1] = 2.0e-01;
		elco[3][4][1] = 2.0e-02;
		elco[4][0][1] = 4.37956204379562043795620437956e-01;
		elco[4][2][1] = 8.21167883211678832116788321168e-01;
		elco[4][3][1] = 3.10218978102189781021897810219e-01;
		elco[4][4][1] = 5.47445255474452554744525547445e-02;
		elco[4][5][1] = 3.64963503649635036496350364964e-03;
		tesco[0][0][1] = F_ONE;
		tesco[0][1][1] = F_TWO;
		tesco[0][2][1] = F_THREE;
		tesco[1][0][1] = F_ONE;
		tesco[1][1][1] = 4.5e00;
		tesco[1][2][1] = F_SIX;
		tesco[2][0][1] = F_HALF;
		tesco[2][1][1] = 7.33333333333333333333333333333e00;
		tesco[2][2][1] = 9.16666666666666666666666666667e00;
		tesco[3][0][1] = 1.66666666666666666666666666667e-01;
		tesco[3][1][1] = 1.04166666666666666666666666667e01;
		tesco[3][2][1] = 1.25e01;
		tesco[4][0][1] = 4.16666666666666666666666666667e-02;
		tesco[4][1][1] = 1.37e01;
		tesco[4][2][1] = 1.59833333333333333333333333333e01;
		_aini = 0;
	}
	*maxder = maxnq[*meth - 1];
	el[0] = elco[*nq - 1][0][*meth - 1];
	scopy(*nq - 1, &elco[*nq - 1][2][*meth - 1], 2, &el[2], 1);
	scopy(3, &tesco[*nq - 1][0][*meth - 1], 2, tq, 1);
	tq[3] = F_HALF * tq[1] / (*nq + 2);
	return;
}				/* end of function */


/* Structured by FOR_STRUCT, v0.2, on 10/09/90 at 09:11:15
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  I12AG/DI12AG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Solve the linear system A*X = B.

    Usage:      CALL I12AG (N, A, LDA, NLC, NUC, B, X, FAC, IPVT, WK,
                            MTYPE, ISTOP)

    Arguments:
       N      - Number of differential equations.  (Input)
       A      - Matrix used when ODE system is implicit.  (Input)
                A is referenced only if PARAM(19) = IATYPE is nonzero.
                Its data structure is determined by PARAM(14) = MTYPE.
                A must be non-singular.  Use only if MITER is 1 or 2.
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       NLC    - Number of lower codiagonals of A.  (Input)
       NUC    - Number of upper codiagonals of A.  (Input)
       B      - Vector of length N containing the right-hand side of the
                linear system.  (Input)
       X      - Vector of length N containing the solution to the linear
                system.  (Output)
       FAC    - Work vector of length N**2 if A is in full form, of
                length (2*NLC+NUC+1)*N if A is in band form.
       IPVT   - Work vector of length N.
       WK     - Work vector of length N.
       MTYPE  - Matrix type for A.  (Input)
                MTYPE=0 selects full matrices.
                MTYPE=1 selects band matrices.
       ISTOP  - Error parameter.  (Output)

    Chapter:    MATH/LIBRARY Differential Equations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_i12ag(Mint *n, Mfloat *a, Mint *lda, Mint *nlc, Mint *nuc,
			Mfloat b[], Mfloat x[], Mfloat imsl_fac[],
			Mint ipvt[], Mfloat wk[], Mint *mtype, Mint *istop)
#else
static void l_i12ag(n, a, lda, nlc, nuc, b, x, imsl_fac, ipvt, wk,
	   mtype, istop)
	Mint            *n;
	Mfloat          *a;
	Mint            *lda, *nlc, *nuc;
	Mfloat           b[], x[], imsl_fac[];
	Mint             ipvt[];
	Mfloat           wk[];
	Mint            *mtype, *istop;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
	static Mint             _l0, ldfac, nra;


	*istop = 0;
	if (*mtype == 0) {
		/* FACTOR */
		imsl_l2trg(*n, a, *lda, imsl_fac, *n, ipvt, wk);
		if (imsl_n1rty(1) != 4) {
			/* SOLUTION STEP */
			_l0 = 1;
			imsl_lfsrg(*n, imsl_fac, *n, ipvt, b, &_l0, x);
		}
	} else if (*mtype == 1) {
		/* FACTOR */

#if 0
		nra = *nlc + *nuc + 1;
		ldfac = 2 ** nlc + *nuc + 1;
		imsl_l2trb(n, a, lda, nlc, nuc, imsl_fac, &ldfac, ipvt, wk);
		/* SOLUTION STEP */
		if (imsl_n1rty(1) != 4) {
			imsl_lfsrb(n, imsl_fac, &ldfac, nlc, nuc, ipvt, b, ADR(_l0, 1),
				   x);
		}
#endif
	} else if (*mtype == 2) {
		/* FACTOR */
#if 0
		imsl_lftds(n, a, lda, imsl_fac, n);
		if (imsl_n1rty(1) != 4) {
			/* SOLUTION STEP */
			imsl_lfsds(n, imsl_fac, n, b, x);
		}
#endif
	} else if (*mtype == 3) {
		nra = *nuc + 1;
		/* FACTOR */
#if 0
		imsl_lftqs(n, a, lda, nuc, imsl_fac, &nra);
		if (imsl_n1rty(1) != 4) {
			/* SOLUTION STEP */
			imsl_lfsqs(n, imsl_fac, &nra, nuc, b, x);
		}
#endif
	}
	if (imsl_n1rty(1) == 4) {

		imsl_ermes(IMSL_FATAL, IMSL_INPUT_MATRIX_A_IS_SINGULAR);
		*istop = 1;
	}
	return;
}				/* end of function */

#undef A

/* Structured by FOR_STRUCT, v0.2, on 10/09/90 at 09:13:32
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  I13AG/DI13AG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Solve an initial value problem for ordinary differential
                equations using an Adams-Moulton or Gear method.

    Usage:      CALL I13AG (ICODE,I1,R1,R2,R3)

    Arguments:  (See IVPAG)

    Chapter:    MATH/LIBRARY Differential Equations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_i13ag(Mint *icode, Mint *i1, Mfloat *r1, Mfloat *r2,
			Mfloat *r3)
#else
static void l_i13ag(icode, i1, r1, r2, r3)
	Mint            *icode, *i1;
	Mfloat          *r1, *r2, *r3;
#endif
{


	if (*icode == 12) {
		imsl_e1sti(1, *i1);
		imsl_e1str(1, *r1);

		imsl_ermes(IMSL_TERMINAL, IMSL_INCORRECT_LDA_VALUE_GIVEN);
	} else if (*icode == 16) {
		imsl_e1sti(1, *i1);
		imsl_e1str(1, *r1);

		imsl_ermes(IMSL_TERMINAL, IMSL_INCORRECT_LDA_GIVEN);
	} else if (*icode == 14) {
		imsl_e1str(1, *r1);

		imsl_ermes(IMSL_TERMINAL, IMSL_ARGUMENT_X_CHANGED_VALUE);
	} else if (*icode == 15) {
		imsl_e1str(1, *r1);

		imsl_ermes(IMSL_TERMINAL, IMSL_ARGUMENT_XEND_IS_UNCHANGED);
	} else if (*icode == 43) {
		imsl_e1sti(1, *i1);

		imsl_ermes(IMSL_FATAL, IMSL_ODE_TOO_MANY_STEPS);
	} else if (*icode == 42) {
		imsl_e1sti(1, (*i1 + 1));
		imsl_e1sti(2, *i1);

		imsl_ermes(IMSL_FATAL, IMSL_ODE_TOO_MANY_EVALS);
	} else if (*icode == 41) {

		imsl_ermes(IMSL_FATAL, IMSL_REPEATED_ERR_TEST_FAILURE);
	} else if (*icode == 45) {
		imsl_e1str(1, *r1);

		imsl_ermes(IMSL_FATAL, IMSL_ODE_FAIL);
	} else if (*icode == 46) {
		imsl_e1str(1, *r1);
		imsl_e1str(2, *r2);

		imsl_ermes(IMSL_FATAL, IMSL_INTEGRATION_HALTED_1);
	} else if (*icode == 47) {
		imsl_e1str(1, *r1);
		imsl_e1str(2, *r2);

		imsl_ermes(IMSL_FATAL, IMSL_INTEGRATION_HALTED_2);
	} else if (*icode == 44) {
		imsl_e1str(1, *r1);
		imsl_e1str(2, *r2);
		imsl_e1str(3, *r3);

		imsl_ermes(IMSL_FATAL, IMSL_TOL_TOO_SMALL_OR_STIFF);
	} else if (*icode == 57) {
		imsl_e1sti(1, *i1);

                imsl_ermes(IMSL_TERMINAL, IMSL_IDO_OUT_OF_RANGE);
	} else if (*icode == 19) {

                imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_IDO_VALUE);
	}
	return;
}				/* end of function */

/* Structured by FOR_STRUCT, v0.2, on 10/08/90 at 17:20:27
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  I2PAG/DI2PAG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Solve an initial value problem for ordinary differential
                equations using an Adams-Moulton or Gear method.

    Usage:      CALL I2PAG (IDO, NEQ, FCN, FCNJ, A, X, XEND, TOL,
                            PARAM, Y, YTEMP, YMAX, ERROR, SAVE1,
                            SAVE2, PW, IPVT, VNORM)

    Arguments:  (See IVPAG)

    Chapter:    MATH/LIBRARY Differential Equations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */

#ifdef ANSI
#if defined (COMPUTER_HP97C)
static void l_i2pag(Mint *ido, Mint *neq, 
			void (*fcn) (Mint, Mfloat, Mfloat*, Mfloat*),
			void (*fcnj) (Mint, Mfloat, Mfloat*, Mfloat*, Mfloat*),
			Mfloat a[], Mfloat *x, Mfloat *xend, Mfloat *tol,
			Mfloat y[], Mfloat ytemp[], Mfloat ymax[],
			Mfloat error[], Mfloat save1[], Mfloat save2[],
			Mfloat pw[], Mint ipvt[], void (*vnorm) ())
#else
static void l_i2pag(Mint *ido, Mint *neq, 
			void (*fcn) (Mint, Mfloat, Mfloat*, Mfloat*),
			void (*fcnj) (Mint, Mfloat, Mfloat[], Mfloat*, Mfloat*),
			Mfloat a[], Mfloat *x, Mfloat *xend, Mfloat *tol,
			Mfloat y[], Mfloat ytemp[], Mfloat ymax[],
			Mfloat error[], Mfloat save1[], Mfloat save2[],
			Mfloat pw[], Mint ipvt[], void (*vnorm) ())
#endif
#else
static void l_i2pag(ido, neq, fcn, fcnj, a, x, xend, tol,
	   y, ytemp, ymax, error, save1, save2, pw, ipvt, vnorm)
	Mint            *ido, *neq;
	void            (*fcn) (), (*fcnj) ();
	Mfloat           a[], *xend, *tol, y[], ytemp[], ymax[],
	                error[], save1[], save2[], pw[];
	Mfloat *x;
	Mint             ipvt[];
	void            (*vnorm) ();
#endif
{
	static Mint             _l0, _l1, i, ixend, jer, k, ker, kflag, lda, mxstep,
	                n, nhcut, nra;
	static Mfloat           _f0, _f1, _f2, ayi, d, dd, harg, hmaxu, hminu,
	                toutp, xendpv, xsave;
	static Mint      istat = 0;



	imsl_e1psh("l_i2pag");

	if (*ido < 1 || *ido > 7) {
		_l0 = 57; _f0 = F_ZERO; _f1 = F_ZERO; _f2 = F_ZERO;
		l_i13ag(&_l0, ido, &_f0, &_f1, &_f2);
		goto L_9001;
	}
/* checked in adams_gear_mgr() */
#if 0
	if (*ido == 1) {
		if (istat == 0) {
			istat = 1;
		} else {
			_l0 = 19; _l1 = 0; _f0 = F_ZERO;
			l_i13ag(&_l0, &_l1, &_f0, &_f0,
				   &_f0);
			goto L_9001;
		}
	} else if (*ido == 2 || *ido == 3) {
		if (istat == 1) {
			if (*ido == 3) {
				istat = 0;
				goto L_9000;
			}
		} else {
			imsl_e1sti(1, *ido);

                        imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_IDO_VALLUE_2);
			goto L_9001;
		}
	}
#endif
	if (*ido == 2)
		goto L_50;
	if (*ido == 3)
		goto L_9000;
	if (*ido == 4 || *ido == 7)
		goto L_70;
	if (*ido == 5)
		goto L_80;
	if (*ido == 6)
		goto L_130;
	/* IDO = 1 (Initial call) */
L_10:
	;
	/* Initialize */
	i11ag.istop = 0;
	i11ag.uround = imsl_amach(4);
	n = *neq;
	hminu = F_ZERO;
	hmaxu = imsl_amach(2);
	mxstep = 500;
	i11ag.mxfcn = 0;
	i11ag.maxord = 12;
	i11ag.intrp1 = 0;
	i11ag.intrp2 = 0;
	i5prk.inorm = 0;
	i5prk.floor = F_ONE;
	i11ag.meth = 1;
	i11ag.miter = 0;
	i11ag.mtype = 0;
	i11ag.epsj = sqrt(i11ag.uround);

#if 0

	/* Check PARAM */
	for (k = 1; k <= 20; k++) {
		if (nint(param[k - 1]) < 0) {
			imsl_e1sti(1, k);
			imsl_e1str(1, param[k - 1]);

			imsl_ermes(IMSL_TERMINAL,
			IMSL_PARAM_CANNOT_BE_NEGATIVE);
		}
	}
	/* Check INORM */
	if (nint(param[9]) > 3) {
		imsl_e1str(1, param[9]);

		imsl_ermes(IMSL_TERMINAL, IMSL_NEED_INORM_LESS_THAN_FOUR);
	}
	if (nint(param[11]) > 2) {
		imsl_e1str(1, param[11]);

		imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_METHOD_INDICATOR);
	}
	for (k = 13; k <= 14; k++) {
		if (nint(param[k - 1]) > 3) {
/*			imsl_e1sti(1, k);
			imsl_e1str(1, param[k - 1]);

			imsl_ermes(5, 4, "The argument PARAM(%(i1)) = %(r1).  Both MITER and MTYPE must be less than four.");
*/
		}
	}

	/* Check IATYPE */
	if (nint(param[18]) >= 3) {
/*		imsl_e1str(1, param[18]);

		imsl_ermes(5, 6, "The type of the matrix A must be 0, 1, or 2, while PARAM(19) = IATYPE = %(r1) is given.");
*/
	}
	if (imsl_n1rty(0) != 0)
		goto L_9001;

#endif

	/* Set HARG */
	harg = sign(lv_state->hinit, *xend - *x);
	hminu = lv_state->hmin;
	hmaxu = lv_state->hmax;
	mxstep = lv_state->max_steps;
	i11ag.mxfcn = lv_state->max_fcn_evals;
	i11ag.maxord = lv_state->maxord;
	i11ag.intrp1 = 0;
	i11ag.intrp2 = 0;
	i5prk.inorm = lv_state->norm;
	i5prk.floor = lv_state->floor;
	i11ag.meth = lv_state->method;
	i11ag.miter = (lv_state->miter)-1;
	i11ag.mtype = 0;
	i11ag.nlc = 0;
	i11ag.nuc = 0;
	i11ag.epsj = sqrt(imsl_amach(4));
	i11ag.iatype = 0;
	lda = 1;
	if (i11ag.miter == 0 || i11ag.miter == 3)
		i11ag.mtype = -1;

	ker = 0;
	i11ag.mf = 10 * i11ag.meth + i11ag.miter;
	/* Check TOL */
	if (*tol <= F_ZERO) {
		imsl_e1str(1, *tol);

		imsl_ermes(IMSL_TERMINAL, IMSL_ODE_NEG_TOL);
	}
	/* Check N */
	if (n < 1) {
		imsl_e1sti(1, n);

		imsl_ermes(IMSL_TERMINAL, IMSL_ODE_NEG_NEQ);
	}

		nra = i11ag.nlc + i11ag.nuc + 1;

	if (imsl_n1rty(0) > 0)
		goto L_9001;

	for (i = 1; i <= n; i++) {
		ymax[i - 1] = fabs(y[i - 1]);
		if (ymax[i - 1] == F_ZERO)
			ymax[i - 1] = F_ONE;
	}
	scopy(n, y, 1, ytemp, 1);
	i11ag.t = *x;
	i11ag.h = harg;
	if ((i11ag.t + i11ag.h) == i11ag.t)
		ker = 33;
	i11ag.hmin = imsl_f_max(fabs(harg), hminu);
	i11ag.hmax = imsl_f_min(fabs(*x - *xend) * F_TEN, hmaxu);
	i11ag.eps = *tol;
	i11ag.jstart = 0;
	nhcut = 0;
	i11ag.el[1] = F_ONE;
	i11ag.oldl0 = F_ONE;
	/*
	 * Set previous XEND initially to initial value of X
	 */
	xendpv = *x;
	ixend = 0;
	goto L_60;
	/* IDO = 2  (Middle calls) */
L_50:
	;
	/*
	 * Case 2 - Normal re-entry (IDO.EQ.2) abort if XEND reached, and
	 * either X changed or XEND not changed
	 */
#if defined(COMPUTER_PMXUX) || defined(COMPUTER_RTXLXS)
        if (ixend != 0 && (fabs(*x - xendpv) > imsl_amach(1))) {
#else
	if (ixend != 0 && *x != xendpv) {
#endif

		/* Print error message */
		l_i13ag(ADR(_l0, 14), ADR(_l1, 0), x, ADR(_f0, F_ZERO), ADR(_f1, F_ZERO));
		goto L_9001;
	} else if (ixend != 0 && *xend == xendpv) {
		/* Print error message */
		l_i13ag(ADR(_l0, 15), ADR(_l1, 0), xend, ADR(_f0, F_ZERO), ADR(_f1, F_ZERO));
		goto L_9001;
	}
	/* Re-initialize flag IXEND */
	ixend = 0;
	/*
	 * TOUTP is the previous value of XEND for use in HMAX.
	 */
	i11ag.hmax = imsl_f_min(fabs(*xend - toutp) * F_TEN, hmaxu);
	goto L_100;

L_60:
	;
	/* Break if INTERUPT 1 */
	if (i11ag.intrp1 > 0) {
		*ido = 4;
		goto L_9000;
	}
	/* Re-enter at 70 */
L_70:
	;
	if (*ido == 7)
		*x = xsave;

	l_i4pag(fcn, fcnj, &n, a, &lda, ytemp, ymax, error, save1, save2,
		   pw, ipvt, &kflag, ido, vnorm);
	if (i11ag.istop == 1)
		goto L_9000;
	if (*ido == 7) {
		xsave = *x;
		*x = i11ag.t;
		goto L_9000;
	}
	/* KFLAG = 0, -1, -2, -3 */
	if (i11ag.nstep >= mxstep) {
		/* Print error message */
		_l0 = 43; _f0 = F_ZERO; _f1 = F_ZERO; _f2 = F_ZERO;
		l_i13ag(&_l0, &mxstep, &_f0, &_f1, &_f2);
		goto L_9000;
	}
	if (i11ag.mxfcn > 0 && i11ag.nfcn >= i11ag.mxfcn) {
		/* Print error message */
		_l0 = 42; _f0 = F_ZERO; _f1 = F_ZERO; _f2 = F_ZERO;
		l_i13ag(&_l0, &i11ag.mxfcn, &_f0, &_f1, &_f2);
		goto L_9000;
	}
	if (kflag == -2) {
		/* Print error message */
		_l0 = 41; _l1 = 0; _f0 = F_ZERO; _f1 = F_ZERO; _f2 = F_ZERO;
		l_i13ag(&_l0, &_l1, &_f0, &_f1, &_f2);
		goto L_140;
	}
	if (kflag != 0)
		goto L_110;

	/*
	 * Normal return from integrator.  The weights YMAX(I) are updated. A
	 * test is made for TOL being too small for the machine precision. If
	 * INDEX = 3, Y is set to the current solution on return. If INDEX =
	 * 2, H is controlled to hit XEND (within roundoff error), and then
	 * the current solution is put in Y on return.  For any other value
	 * of INDEX, control returns to the integrator unless XEND has been
	 * reached.  Then interpolated values of the solution are computed
	 * and stored in Y on return.  If interpolation is not desired, the
	 * call to I9PAG should be removed and control transferred to
	 * statement 95 instead of 105.
	 * 
	 * Break if INTERUPT 2
	 */
	if (i11ag.intrp2 > 0) {
		*ido = 5;
		goto L_9000;
	}
	/* Re-enter at 38 */
L_80:
	;

	d = F_ZERO;
	for (i = 1; i <= n; i++) {
		ayi = fabs(ytemp[i - 1]);
		ymax[i - 1] = imsl_f_max(ymax[i - 1], ayi);
	}
	(*vnorm) (n, ytemp, ytemp, ymax, &dd);
#ifdef COMPUTER_APLC
        if (dd >= sqrt(imsl_f_machine(2)))
                dd = sqrt(imsl_f_machine(2));
#endif
	d = dd * dd;
	d *= imsl_fi_power(i11ag.uround / *tol, 2);
	if (d > (Mfloat) (n)) {
		/* Print error message */
		_l0 = 45; _l1 = 0; _f0 = F_ZERO; _f1 = F_ZERO;
		l_i13ag(&_l0, &_l1, tol, &_f0, &_f1);
		kflag = -2;
		goto L_140;
	}
L_100:
	if ((i11ag.t - *xend) * i11ag.h < F_ZERO) {
		if ((i11ag.t + i11ag.h) == i11ag.t)
			ker = 33;
		goto L_60;
	} else {
		l_i9pag(&n, xend, ytemp, y);
		*x = *xend;
		/*
		 * Set IDO=2, successful.  Save XEND, set flag
		 */
		*ido = 2;
		lv_state->ido = 2;
		xendpv = *xend;
		ixend = 1;
		goto L_150;
	}

	/*
	 * On an error return from integrator, an immediate return occurs if
	 * KFLAG = -2, and recovery attempts are made otherwise. To recover,
	 * H and HMIN are reduced by a factor of .1 up to 10 times before
	 * giving up.
	 */
L_110:
	;
	if (kflag != -1 && kflag != -3)
		goto L_140;
	if (nhcut != 10)
		goto L_120;
	if (kflag == -1) {
		/* Print error message */
		_l0 = 46; _l1 = 0; _f0 = F_ZERO;
		l_i13ag(&_l0, &_l1, &i11ag.h, tol, &_f0);
		goto L_9000;
	} else {
		/* Print error message */
		_l0 = 47; _l1 = 0; _f0 = F_ZERO;
		l_i13ag(&_l0, &_l1, &i11ag.h, tol, &_f0);
		jer = 133;
		goto L_9000;
	}
L_120:
	;
	/*
	 * Break with IDO = 6 (unsuccessful attempt) if INTERUPT 2
	 */
	if (i11ag.intrp2 > 0) {
		*ido = 6;
		goto L_9000;
	}
	/* Re-enter at 91 */
L_130:
	;

	nhcut += 1;
	i11ag.hmin = imsl_f_max(i11ag.hmin * .1, hminu);
	i11ag.h *= .1;
	i11ag.jstart = -1;
	if ((i11ag.t + i11ag.h) == i11ag.t)
		ker = 33;
	goto L_60;

L_140:
	*x = i11ag.t;
	scopy(n, ytemp, 1, y, 1);
L_150:
	;

	toutp = *x;
L_9000:
	;
	if (*ido >= 4) {
		scopy(n, ytemp, 1, y, 1);
		*x = i11ag.t;
	}
	harg = i11ag.hused;
	if (kflag != 0)
		harg = i11ag.h;
/*
	param[30] = harg;
	param[31] = i11ag.hmin;
	param[32] = i11ag.hmax;
	param[33] = i11ag.nstep;
	param[34] = i11ag.nfcn;
	param[35] = i11ag.nje;
*/
	if (lv_state->nstep != NULL) *lv_state->nstep = i11ag.nstep;
	if (lv_state->nfcn != NULL)  *lv_state->nfcn = i11ag.nfcn;
	if (lv_state->nfcnj != NULL) *lv_state->nfcnj = i11ag.nje;

	if (imsl_n1rty(0) < 3 && ker == 33) {
		/* Print error message */
		_l0 = 44; _l1 = 0;
		l_i13ag(&_l0, &_l1, &i11ag.t, &harg, tol);
	}
L_9001:
	imsl_e1pop("l_i2pag");

	return;
}				/* end of function */


/* Structured by FOR_STRUCT, v0.2, on 10/08/90 at 17:23:03
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  I4PAG/DI4PAG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Solve an initial value problem for ordinary differential
                equations using an Adams-Moulton or Gear method.

    Usage:      CALL I4PAG (FCN, FCNJ, N, A, LDA, Y, YMAX, ERROR, SAVE1,
                            SAVE2, PW, IPVT, KFLAG, IDO, VNORM)

    Arguments:
       FCN    - User-supplied SUBROUTINE to evaluate functions.
                The usage is
                CALL FCN (NEQ, X, Y, YPRIME), where
                NEQ    - Number of equations.  (Input)
                X      - Independent variable.  (Input)
                Y      - Array of length NEQ containing the dependent
                         variable values.  (Input)
                YPRIME - Array of length NEQ containing the values of
                         the right hand side of the system of
                         equations.  (Output)
                         See remark 3.
                FCN must be declared EXTERNAL in the calling program.
       FCNJ   - User-supplied SUBROUTINE to compute the Jacobian.
                The usage is
                CALL FCNJ (NEQ, X, Y, DYPDY), where
                NEQ    - Number of equations.  (Input)
                X      - Independent variable.  (Input)
                Y      - Array of length NEQ containing the dependent
                         variable values.  (Input)
                DYPDY  - Matrix (whose type is determined by MTYPE,
                         containing dYPRIME(i)/dY(j).  (Output)
                FCNJ must be declared EXTERNAL in the calling program.
                If IATYPE is nonzero, then FCNJ should
                compute the Jacobian of the right hand side of the
                equation Ay'=f(x,y).
                FCNJ is used only if MITER = 1.
       N      - Number of differential equations.  (Input)
       A      - Matrix used when ODE system is implicit.  (Input)
                A is referenced only if IATYPE is nonzero.
                Its data structure is determined by MTYPE.
                A must be non-singular.  Use only if MITER is 1 or 2.
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
       Y      - Array of size N by IYCOL, whose first column holds the
                present solution, and whose other columns are used as
                workspace.  If METH is 1, then IYCOL = 13.  If METH is 2,
                then IYCOL = 6.
       YMAX   - Vector of length N containing the maximum Y values
                computed so far.  (Output)
       ERROR  - Vector of length N containing error estimates at each
                element of Y.  (Output)
       SAVE1  - Work vector of length N.
       SAVE2  - Work vector of length N.
       PW     - Work vector of length NPW.  PW is used both to store the
                Jacobian and as workspace.
                See IVPAG for the value of NPW.
       IPVT   - Work vector of length N.
       KFLAG  - Scalar to indicate status of the result.  (Output)
       IDO    - Flag indicating the state of the computation.
                (Input/Output)
                1        Initial entry
                2        Normal reentry
                3        Final call
                4        Return because of INTERRUPT 1
                5        Return because of INTERRUPT 2 with step accepted
                6        Return because of INTERRUPT 2 with step rejected
                7        Return for new value of A.  The matrix A at X is
                         must be recomputed and IVPAG/DIVPAG called
                         again.  No other argument (including IDO) should
                         be changed.  This value of IDO is returned only
                         if IATYPE=2.
       VNORM  - User-supplied SUBROUTINE to compute the norm of the
                error.
                The routine may be provided by the user or the IMSL
                routine I3PRK/DI3PRK may be used.
                The usage is
                CALL VNORM (NEQ, V, Y, YMAX, ENORM), where
                NEQ    - Number of equations.  (Input)
                V      - Vector of length NEQ containing the vector whose
                         norm is to be computed.  (Input)
                Y      - Vector of length NEQ containing the values of
                         the dependent variable.  (Input)
                YMAX   - Vector of length NEQ containing the maximum Y
                         values computed so far.  (Input)
                ENORM  - Norm of the vector V.  (Output)
                VNORM must be declared EXTERNAL in the calling program.

    Chapter:    MATH/LIBRARY Differential Equations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */

#ifdef ANSI
#if defined(COMPUTER_HP97C)
static void l_i4pag(void (*fcn) (Mint, Mfloat, Mfloat*, Mfloat*),
                    void (*fcnj) (Mint, Mfloat, Mfloat*, Mfloat*, Mfloat*),
			Mint *n, Mfloat *a,
			Mint *lda, Mfloat *y, Mfloat ymax[], Mfloat error[],
			Mfloat save1[], Mfloat save2[], Mfloat pw[],
			Mint ipvt[], Mint *kflag, Mint *ido,
			void (*vnorm) ())
#else
static void l_i4pag(void (*fcn) (Mint, Mfloat, Mfloat*, Mfloat*),
                    void (*fcnj) (Mint, Mfloat, Mfloat[], Mfloat*, Mfloat*),
			Mint *n, Mfloat *a,
			Mint *lda, Mfloat *y, Mfloat ymax[], Mfloat error[],
			Mfloat save1[], Mfloat save2[], Mfloat pw[],
			Mint ipvt[], Mint *kflag, Mint *ido,
			void (*vnorm) ())
#endif
#else
static void l_i4pag(fcn, fcnj, n, a, lda, y, ymax, error, save1,
	   save2, pw, ipvt, kflag, ido, vnorm)
	void            (*fcn) (), (*fcnj) ();
	Mint            *n;
	Mfloat          *a;
	Mint            *lda;
	Mfloat          *y, ymax[], error[], save1[], save2[], pw[];
	Mint             ipvt[], *kflag, *ido;
	void            (*vnorm) ();
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
#define Y(I_,J_)	(y+(I_)*(*n)+(J_))
	static Mint      _l0, _l1, i, idoub, ier, iredo, iret, iweval, j,
	                j1, j2, k, l, ldfac, ldpw, lmax_, m, maxder, meo,
	                mfold, mio, newq, nold, nq, nstepj, nwk, nwka;
	static Mfloat    _f0, _f1, bnd, crate, d, d1, dd, dd1, e, edn, epsold, eup,
	                hl0, hold, pr1, pr2, pr3, r1, rc, rh, rmax, told,
	                tq[4];
	static Mfloat 		*yprime_in_fcnj;

	/*
	 * THIS ROUTINE PERFORMS ONE STEP OF THE INTEGRATION OF AN INITIAL
	 * VALUE PROBLEM FOR A SYSTEM OF ORDINARY DIFFERENTIAL EQUATIONS.
	 */
	if (*ido == 7)
		goto L_110;

	*kflag = 0;
	told = i11ag.t;

	if (i11ag.jstart > 0) {
		goto L_60;
	} else if (i11ag.jstart == 0) {
		/*
		 * On the first call, the ORDER is set to 1 and the initial
		 * YDOT is calculated.  RMAX is the maximum ratio by which H
		 * can be increased in a single step.  It is initially 1.E4
		 * to compensate for the small initial H, but then is
		 * normally equal to 10.  If a failure occurs (in corrector
		 * convergence or error test), RMAX is set at 2 for the next
		 * increase. If MITER is 0 or 3, then MTYPE = -1
		 */
		if (i11ag.mtype == 0) {
			nwk = *n ** n + 1;
		} else if (i11ag.mtype == 1) {
			ldpw = i11ag.nuc + i11ag.nlc + 1;
			ldfac = 2 * i11ag.nlc + i11ag.nuc + 1;
			nwk = ldfac ** n + 1;
		} else if (i11ag.mtype == 2) {
			nwk = *n ** n + 1;
		} else if (i11ag.mtype == 3) {
			ldpw = i11ag.nuc + 1;
			nwk = *n * (i11ag.nuc + 1) + 1;
		} else if (i11ag.miter == 3) {
			nwk = *n + 1;
		} else {
			nwk = 1;
		}
		nwka = nwk + 2 ** n;

		(*fcn) (*n, i11ag.t, y, save1);
		for (i = 1; i <= *n; i++) {
			*Y(1, i - 1) = i11ag.h * save1[i - 1];
		}
		/* Multiply by inv(A) */
		if (i11ag.iatype == 1 || i11ag.iatype == 2) {
			l_i12ag(n, a, lda, &i11ag.nlc, &i11ag.nuc, Y(1, 0), save2,
				   &pw[nwka - 1], ipvt, error, &i11ag.mtype, &i11ag.istop);
			if (i11ag.istop == 1)
				goto L_9000;
			scopy(*n, save2, 1, Y(1, 0), 1);
		}
		i11ag.meth = i11ag.mf / 10;
		i11ag.miter = i11ag.mf - 10 * i11ag.meth;
		nq = 1;
		l = 2;
		idoub = 3;
		rmax = 1.0e4;
		rc = F_ZERO;
		crate = F_ONE;
		hold = i11ag.h;
		mfold = i11ag.mf;
		i11ag.nstep = 0;
		nstepj = 0;
		i11ag.nfcn = 1;
		i11ag.nje = 0;
		iret = 3;
	} else {
		/*
		 * If the caller has changed METH, I10AG is called to set the
		 * coefficients of the method.  If the caller has changed N,
		 * EPS, or METH, the constants E, EDN, EUP, and BND must be
		 * reset.  E is a comparison for errors of the current order
		 * NQ.  EUP is to test for increasing the order, EDN for
		 * decreasing the order.  BND is used to test for convergence
		 * of the corrector iterates.  If the caller has changed H, Y
		 * must be rescaled.  If H or METH has been changed, IDOUB is
		 * reset to L + 1 to prevent further changes in H for that
		 * many steps.
		 */
		if (i11ag.mf == mfold) {
			if ((i11ag.eps == epsold) && (*n == nold)) {
				if (i11ag.h == hold)
					goto L_60;
				rh = i11ag.h / hold;
				i11ag.h = hold;
				iredo = 3;
				goto L_40;
			}
			if (*n == nold)
				iweval = i11ag.miter;
			iret = 1;
			goto L_30;
		}
		meo = i11ag.meth;
		mio = i11ag.miter;
		i11ag.meth = i11ag.mf / 10;
		i11ag.miter = mod(i11ag.mf, 10);
		mfold = i11ag.mf;
		if (i11ag.miter != mio)
			iweval = i11ag.miter;
		if (i11ag.meth == meo) {
			if ((i11ag.eps == epsold) && (*n == nold)) {
				if (i11ag.h == hold)
					goto L_60;
				rh = i11ag.h / hold;
				i11ag.h = hold;
				iredo = 3;
				goto L_40;
			}
			if (*n == nold)
				iweval = i11ag.miter;
			iret = 1;
			goto L_30;
		}
		idoub = l + 1;
		iret = 1;
	}
L_20:
	;
	l_i10ag(&i11ag.meth, &nq, i11ag.el, tq, &maxder);
	lmax_ = imsl_i_min(maxder, i11ag.maxord) + 1;
	rc = rc * i11ag.el[0] / i11ag.oldl0;
	i11ag.oldl0 = i11ag.el[0];
L_30:
	;
	edn = *n * imsl_fi_power(tq[0] * i11ag.eps, 2);
	e = *n * imsl_fi_power(tq[1] * i11ag.eps, 2);
	eup = *n * imsl_fi_power(tq[2] * i11ag.eps, 2);
	bnd = *n * imsl_fi_power(tq[3] * i11ag.eps, 2);
	epsold = i11ag.eps;
	nold = *n;
	if (iret == 1) {
		if (i11ag.h == hold)
			goto L_60;
		rh = i11ag.h / hold;
		i11ag.h = hold;
		iredo = 3;
	} else if (iret == 2) {
		rh = imsl_f_max(rh, i11ag.hmin / fabs(i11ag.h));
	} else if (iret == 3) {
		goto L_60;
	}
L_40:
	;
	rh = imsl_f_vmin(3, rh, i11ag.hmax / fabs(i11ag.h), rmax);
	r1 = F_ONE;
	for (j = 2; j <= l; j++) {
		r1 *= rh;
		sscal(*n, r1, Y(j - 1, 0), 1);
	}
	i11ag.h = imsl_f_min(i11ag.h * rh, i11ag.hmax);
	rc *= rh;
	idoub = l + 1;
	if (iredo == 0) {
		rmax = F_TEN;
		hold = i11ag.h;
		i11ag.jstart = nq;
		goto L_9000;
	}
	/*
	 * RC is the ratio of new to old values of the coefficient H*EL(1).
	 * Force the partials to be updated when RC differs from 1. by more
	 * than 30 percent or after 20 steps
	 */
L_60:
	;
	if (fabs(rc - F_ONE) > 0.3 || i11ag.nstep >= nstepj + 20) {
		iweval = i11ag.miter;
	}
	i11ag.t += i11ag.h;
	/*
	 * Compute the predicted values by effectively multiplying the Y
	 * array by the Pascal triangle matrix.
	 */
	for (j1 = 1; j1 <= nq; j1++) {
		for (j2 = j1; j2 <= nq; j2++) {
			j = (nq + j1) - j2;
			for (i = 1; i <= *n; i++) {
				*Y(j - 1, i - 1) += *Y(j, i - 1);
			}
		}
	}

	/*
	 * Up to 3 corrector iterations are taken.  A convergence test is
	 * made on the r.m.s. norm of each correction, using BND, which is
	 * dependent on EPS. The sum of the corrections is accumulated in the
	 * vector ERROR(I).  The Y array is not altered in the corrector
	 * loop.  The updated Y vector is stored temporarily in SAVE1.
	 */
L_100:
	;
	m = 0;
	if (i11ag.iatype == 0) {
		(*fcn) (*n, i11ag.t, y, save2);
		scopy(*n, save2, 1, &pw[nwk - 1], 1);
	} else {
		(*fcn) (*n, i11ag.t, y, &pw[nwk - 1]);
	}
	/*
	 * If (IATYPE = 2) return for A with IDO = 7
	 */
	if (i11ag.iatype == 2) {
		i11ag.ido1 = *ido;
		*ido = 7;
		goto L_9000;
	}
	/* Return here on re-entry */
L_110:
	;
	if (i11ag.iatype == 2)
		*ido = i11ag.ido1;
	if (i11ag.iatype == 1 || i11ag.iatype == 2) {
		scopy(*n, &pw[nwk - 1], 1, save2, 1);
	}
	sset(*n, F_ZERO, error, 1);
	i11ag.nfcn += 1;
	if (i11ag.mxfcn > 0 && i11ag.nfcn >= i11ag.mxfcn)
		goto L_9000;
	/*
	 * If indicated, the matrix P = A - H*EL(1)*J is reevaluated before
	 * starting the corrector iteration.  IWEVAL is set to 0 as an
	 * indicator that this has been done.
	 */
	if (iweval > 0) {
		iweval = 0;
		rc = F_ONE;
		i11ag.nje += 1;
		nstepj = i11ag.nstep;
		/*
		 * Evaluate the Jacobian Zero the Jacobian
		 */
		if (i11ag.mtype == 0) {
			sset(imsl_ii_power(*n, 2), F_ZERO, pw, 1);
		} else if (i11ag.mtype == 1) {
			sset(*n * (i11ag.nlc + i11ag.nuc + 1), F_ZERO, pw, 1);
		} else if (i11ag.mtype == 2) {
			sset(imsl_ii_power(*n, 2), F_ZERO, pw, 1);
		} else if (i11ag.mtype == 3) {
			sset(*n * (i11ag.nuc + 1), F_ZERO, pw, 1);
		}
		if (i11ag.miter == 1) {
			/* Evaluate user supplied Jacobian */
			yprime_in_fcnj = (Mfloat *) 
				imsl_malloc(*n * sizeof(Mfloat));
			(*fcn) (*n, i11ag.t, y, yprime_in_fcnj);
			(*fcnj) (*n, i11ag.t, y, yprime_in_fcnj, pw);
			if (yprime_in_fcnj != NULL)
				imsl_free(yprime_in_fcnj);
		} else if (i11ag.miter == 2) {
			/* Evaluate finite difference Jacobian */
			if (i11ag.mtype == 0 || i11ag.mtype == 2) {
				/* Full matrix Jacobian */
				l_i5pag(fcn, n, pw, save1, &pw[nwk - 1], y, ymax);
			} else if (i11ag.mtype == 1 || i11ag.mtype == 3) {
				if (i11ag.mtype == 3) {
					i11ag.nlc = 0;
				}
				/* Band matrix Jacobian */
				l_i6pag(fcn, n, pw, ADR(_l0, i11ag.nuc + i11ag.nlc +
							   1), save1, &pw[nwk - 1], y, ymax, &pw[nwk + *n - 1]);
				if (i11ag.mtype == 3)
					i11ag.nlc = i11ag.nuc;
			}
			/* Diagonal Jacobian */
		} else if (i11ag.miter == 3) {
			l_i7pag(fcn, n, pw, save1, save2, y, ymax, &hl0, &ier);
			if (ier == 2)
				goto L_260;
			goto L_200;
		}
		/* Form PW = A - h*el(1)*J */
		if (i11ag.mtype == 0 || i11ag.mtype == 2) {
			sscal(imsl_ii_power(*n, 2), -i11ag.h * i11ag.el[0], pw, 1);
		} else if (i11ag.mtype == 1) {
			sscal((i11ag.nuc + i11ag.nlc + 1) ** n, -i11ag.h * i11ag.el[0],
				   pw, 1);
		} else if (i11ag.mtype == 3) {
			sscal((i11ag.nuc + 1) ** n, -i11ag.h * i11ag.el[0], pw, 1);
		}
		/*
		 * Add matrix A
		 */
		if (i11ag.iatype == 0) {
			/* A is the identity matrix */
			if (i11ag.mtype == 0 || i11ag.mtype == 2) {
				sadd(*n, F_ONE, pw, *n + 1);
			} else if (i11ag.mtype == 1) {
				sadd(*n, F_ONE, &pw[i11ag.nuc], i11ag.nuc + i11ag.nlc +
					  1);
			} else if (i11ag.mtype == 3) {
				sadd(*n, F_ONE, &pw[i11ag.nuc], i11ag.nuc + 1);
			}
		} else {
			/* A is user supplied */
			if (i11ag.mtype == 0 || i11ag.mtype == 2) {
				k = 1;
				for (j = 1; j <= *n; j++) {
					for (i = 1; i <= *n; i++) {
						pw[k - 1] += *A(j - 1, i - 1);
						k += 1;
					}
				}
			} else if (i11ag.mtype == 1 || i11ag.mtype == 3) {
				for (i = 1; i <= ldpw; i++) {
					for (j = imsl_i_max(1, i11ag.nuc + 2 - i); j <= imsl_i_min(*n,
					     *n + i11ag.nuc + 1 - i); j++) {
						pw[i + (j - 1) * ldpw - 1] += *A(j - 1, i - 1);
					}
				}
			}
		}
		/*
		 * Factor matrix Full matrix
		 */
		if (i11ag.mtype == 0) {
			imsl_l2trg(*n, pw, *n, pw, *n, ipvt, &pw[nwk - 1]);
			/* Band matrix */
		} else if (i11ag.mtype == 1) {
			/*
			 * Set up PW so that L2TRB can be used with PW and
			 * FAC sharing the same storage space
			 */
			for (i = *n - 1; i >= 0; i--) {
				scopy(i11ag.nlc + i11ag.nuc + 1, &pw[(i) * (i11ag.nlc + i11ag.nuc + 1)],
					   -1, &pw[(i) * (2 * i11ag.nlc + i11ag.nuc + 1)], -1);
			}
			for (i = *n - 1; i >= 0; i--) {
				sset(i11ag.nlc, F_ZERO, &pw[i11ag.nlc + i11ag.nuc + (i) * (2 * i11ag.nlc + i11ag.nuc + 1) + 1],
					  1);
			}

			imsl_l2trb(n, pw, &ldfac, &i11ag.nlc, &i11ag.nuc, pw, &ldfac,
				   ipvt, &pw[nwk - 1]);
			/* Positive definite matrix */
		} else if (i11ag.mtype == 2) {
/*
			imsl_lftds(n, pw, n, pw, n);
*/
			/* Positive definite band matrix */
		} else if (i11ag.mtype == 3) {
/*
			imsl_lftqs(n, pw, &ldpw, &i11ag.nuc, pw, &ldpw);
*/
		}
		/*
		 * If the matrix PW is singular and an error condition
		 * arises, a new step size will be chosen without issuing
		 * this error message.
		 */
		if (imsl_n1rty(1) != 0) {
			imsl_e1mes(0, 0, " ");
			goto L_260;
		}
	}
	/* Update, compute inv(A-hbJ)*(hf-hAy') */
L_180:
	;
	/* Set SAVE1 to an estimate of -hy' */
	for (i = 1; i <= *n; i++) {
		save1[i - 1] = -*Y(1, i - 1) - error[i - 1];
	}

	if (i11ag.iatype > 0) {
		/*
		 * Multiply SAVE1 by A if A is user-provided
		 */
		scopy(*n, save1, 1, &pw[nwk - 1], 1);
		if (i11ag.mtype == 0) {
			imsl_sgemv("N", sizeof("N"), n, n, ADR(_f0, F_ONE), a, lda, &pw[nwk - 1],
			    ADR(_l0, 1), ADR(_f1, F_ZERO), save1, ADR(_l1, 1));
		} else if (i11ag.mtype == 2) {
			imsl_ssymv("U", sizeof("U"), n, ADR(_f0, F_ONE), a, lda, &pw[nwk - 1],
			    ADR(_l0, 1), ADR(_f1, F_ZERO), save1, ADR(_l1, 1));
		} else if (i11ag.mtype == 1) {
			l_sgbmv("N", sizeof("N"), n, n, &i11ag.nlc, &i11ag.nuc,
				   ADR(_f0, F_ONE), a, lda, &pw[nwk - 1], ADR(_l0, 1), ADR(_f1, F_ZERO),
				   save1, ADR(_l1, 1));
		} else if (i11ag.mtype == 3) {
			l_ssbmv("U", sizeof("U"), n, &i11ag.nuc, ADR(_f0, F_ONE), a,
				   lda, &pw[nwk - 1], ADR(_l0, 1), ADR(_f1, F_ZERO), save1, ADR(_l1, 1));
		}
	}
	/* Add h*SAVE2 to SAVE1 */
	saxpy(*n, i11ag.h, save2, 1, save1, 1);
	/* Compute (A-hbJ)**(-1)*SAVE1 */
	if (i11ag.miter == 1 || i11ag.miter == 2) {
		/* Full matrix */
		if (i11ag.mtype == 0) {
			_l0 = 1;
			imsl_lfsrg(*n, pw, *n, ipvt, save1, &_l0, save1);
		} else if (i11ag.mtype == 2) {
/*
			imsl_lfsds(n, pw, n, save1, save1);
*/
			/* Band matrix */
		} else if (i11ag.mtype == 1) {
/*
			imsl_lfsrb(n, pw, &ldfac, &i11ag.nlc, &i11ag.nuc, ipvt, save1,
				   ADR(_l0, 1), save1);
*/
		} else if (i11ag.mtype == 3) {
/*
			imsl_lfsqs(n, pw, &ldpw, &i11ag.nuc, save1, save1);
*/
			/*
			 * It is not possible for an error condition to arise
			 * from the above linear solvers
			 */
		}
	} else if (i11ag.miter == 3) {
		l_i8pag(fcn, n, pw, save1, y, &hl0, &ier);
		if (ier == 3) {
			iweval = i11ag.miter;
			goto L_100;
		}
	}
L_200:
	;
	(*vnorm) (*n, save1, y, ymax, &dd);
#ifdef COMPUTER_APLC
	if (dd >= sqrt(imsl_f_machine(2)))
		dd = sqrt(imsl_f_machine(2));
#endif
	d = dd * dd;
	for (i = 1; i <= *n; i++) {
		error[i - 1] += save1[i - 1];
		save1[i - 1] = *Y(0, i - 1) + i11ag.el[0] * error[i - 1];
	}
	/*
	 * Test for convergence.  If M.GT.0, the square of the convergence
	 * rate constant is estimated as CRATE, and this is used in the test.
	 */
	if (m != 0)
		crate = imsl_f_max(.9 * crate, d / d1);
	if ((d * imsl_f_min(F_ONE, F_TWO * crate)) <= bnd) {
		/*
		 * The corrector has converged.  IWEVAL is set to -1 if
		 * partial derivatives were used, to signal that they may
		 * need updating on subsequent steps.
		 */
		if (i11ag.miter != 0)
			iweval = -1;
		i11ag.nfcn += m;
		if (i11ag.mxfcn > 0 && i11ag.nfcn >= i11ag.mxfcn)
			goto L_9000;
		(*vnorm) (*n, error, y, ymax, &dd);
		d = dd * dd;
		if (d <= e) {
			/*
			 * A successful step Increase on the next step.  If A
			 * change in H is considered, an increase or decrease
			 * in order by one is considered also.  a change in H
			 * is made only if it is by a factor of at least 1.1.
			 * If not, IDOUB is set to 10 to prevent testing for
			 * that many steps.
			 */
			*kflag = 0;
			iredo = 0;
			i11ag.nstep += 1;
			i11ag.hused = i11ag.h;
			/* Update the Y array */
			imsl_sger(*n, l, F_ONE, error, 1, i11ag.el, 1, y, *n);
			if (idoub > 1) {
				/* Consider changing H */
				idoub -= 1;
				if (idoub <= 1 && l < lmax_) {
					/*
					 * The error is saved for a possible
					 * order increase on the next step
					 */
					scopy(*n, error, 1, Y(lmax_ - 1, 0), 1);
				}
				hold = i11ag.h;
				i11ag.jstart = nq;
				goto L_9000;
			}
			if (l < lmax_) {
				for (i = 1; i <= *n; i++) {
					pw[nwk + i - 2] = error[i - 1] - *Y(lmax_ - 1, i - 1);
				}
				(*vnorm) (*n, &pw[nwk - 1], y, ymax, &dd1);
				d1 = dd1 * dd1;
				pr3 = (pow(d1 / eup, F_HALF / (l + 1))) * 1.4 + 1.4e-6;
			} else {
				/* To avoid an order increase */
				pr3 = 1.0e20;
			}
			/*
			 * The error test failed.  KFLAG keeps track of
			 * multiple failures.
			 */
		} else {
			*kflag -= 1;
			/* Restore T to its previous value */
			i11ag.t = told;
			/* Restore the Y array */
			for (j1 = 1; j1 <= nq; j1++) {
				for (j2 = j1; j2 <= nq; j2++) {
					j = (nq + j1) - j2;
					for (i = 1; i <= *n; i++) {
						*Y(j - 1, i - 1) -= *Y(j, i - 1);
					}
				}
			}
			/*
			 * Prepare to try step again compute the optimum step
			 * size for this or one lower order
			 */
			rmax = F_TWO;
			if (fabs(i11ag.h) <= i11ag.hmin * 1.00001) {
				*kflag = -1;
				hold = i11ag.h;
				i11ag.jstart = nq;
				goto L_9000;
			}
			if (*kflag <= -3) {
				/* 3 more failures have occured. */
				if (*kflag == -7) {
					/*
					 * After a total of 7 failures, an
					 * exit is taken with KFLAG = -2
					 */
					*kflag = -2;
					hold = i11ag.h;
					i11ag.jstart = nq;
					goto L_9000;
				}
				/* H is reduced by a factor of 10 */
				rh = imsl_f_max(i11ag.hmin / fabs(i11ag.h), 0.1);
				i11ag.h = imsl_f_min(i11ag.h * rh, i11ag.hmax);
				/*
				 * It is assumed that the derivatives that
				 * have accumulated in the Y array have
				 * errors of the wrong order.  Hence the
				 * first derivative is recomputed, and the
				 * order is set to 1.
				 */
				(*fcn) (*n, i11ag.t, y, save1);
				i11ag.nfcn += 1;
				if (i11ag.mxfcn > 0 && i11ag.nfcn >= i11ag.mxfcn)
					goto L_9000;
				for (i = 1; i <= *n; i++) {
					*Y(1, i - 1) = i11ag.h * save1[i - 1];
				}
				/* The step is retried */
				iweval = i11ag.miter;
				idoub = 10;
				if (nq == 1) {
					goto L_60;
				} else {
					nq = 1;
					l = 2;
					iret = 3;
					goto L_20;
				}
			}
			iredo = 2;
			pr3 = 1.0e20;
		}
		/*
		 * Regardless of the success or failure of the step, factors
		 * PR1, PR2, and PR3 are computed, by which H could be
		 * divided at order NQ - 1, order NQ, or order NQ + 1,
		 * respectively. The smallest of these is determined and the
		 * new order chosen accordingly.  If the order is to be
		 * increased, we compute one additional scaled derivative.
		 */
		pr2 = (pow(d / e, F_HALF / l)) * 1.2 + 1.2e-6;
		pr1 = 1.0e20;
		if (nq != 1) {
			(*vnorm) (*n, Y(l - 1, 0), y, ymax, &dd);
			d = dd * dd;
			pr1 = (pow(d / edn, F_HALF / nq)) * 1.3 + 1.3e-6;
		}
		if (pr3 > pr2) {
			if (pr2 > pr1) {
				newq = nq - 1;
				rh = F_ONE / pr1;
				if (*kflag != 0)
					rh = imsl_f_min(rh, F_ONE);
			} else {
				newq = nq;
				rh = F_ONE / pr2;
			}
		} else {
			if (pr3 < pr1) {
				newq = l;
				rh = F_ONE / pr3;
				if (rh >= 1.1) {
					for (i = 1; i <= *n; i++) {
						*Y(newq, i - 1) = error[i - 1] * i11ag.el[l - 1] /
							l;
					}
					nq = newq;
					l = nq + 1;
					iret = 2;
					goto L_20;
				} else {
					idoub = 10;
					hold = i11ag.h;
					i11ag.jstart = nq;
					goto L_9000;
				}
			}
		}
		if ((*kflag == 0) && (rh < 1.1)) {
			idoub = 10;
			hold = i11ag.h;
			i11ag.jstart = nq;
			goto L_9000;
		}
		/*
		 * If there is a change of order, reset NQ, L, and the
		 * coefficients. in any case H is reset according to RH and
		 * the Y array is rescaled. then exit if the step was OK, or
		 * redo the step otherwise.
		 */
		if (newq == nq) {
			rh = imsl_f_max(rh, i11ag.hmin / fabs(i11ag.h));
			goto L_40;
		} else {
			nq = newq;
			l = nq + 1;
			iret = 2;
			goto L_20;
		}
	}
	d1 = d;
	m += 1;
	if (m < 3) {
		(*fcn) (*n, i11ag.t, save1, save2);
		goto L_180;
		/*
		 * The corrector iteration failed to converge in 3 tries
		 */
	} else {
		i11ag.nfcn += 2;
		if (i11ag.mxfcn > 0 && i11ag.nfcn >= i11ag.mxfcn)
			goto L_9000;
		if (iweval == -1) {
			/*
			 * Partials are involved but are not up to date, they
			 * are reevaluated for the next try
			 */
			iweval = i11ag.miter;
			goto L_100;
		}
	}
L_260:
	;
	/* H is reduced */
	rmax = F_TWO;
	/*
	 * The Y array is retracted to its values before prediction
	 */
	i11ag.t = told;
	for (j1 = 1; j1 <= nq; j1++) {
		for (j2 = j1; j2 <= nq; j2++) {
			j = (nq + j1) - j2;
			for (i = 1; i <= *n; i++) {
				*Y(j - 1, i - 1) -= *Y(j, i - 1);
			}
		}
	}
	if (fabs(i11ag.h) <= i11ag.hmin * 1.00001) {
		/* No-convergence exit is taken */
		*kflag = -3;
		hold = i11ag.h;
		i11ag.jstart = nq;
		goto L_9000;
	}
	rh = .25;
	iredo = 1;
	rh = imsl_f_max(rh, i11ag.hmin / fabs(i11ag.h));
	goto L_40;

L_9000:
	;
	return;
}				/* end of function */

#undef A
#undef Y

/* Structured by FOR_STRUCT, v0.2, on 10/08/90 at 17:24:07
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  I5PAG/DI5PAG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Compute the Jacobian, in full matrix mode, using
                finite differences.

    Usage:      CALL I5PAG (FCN, N, PW, SAVE1, SAVE2, Y, YMAX)

    Arguments:
       FCN    - User-supplied SUBROUTINE to evaluate functions.
                The usage is
                CALL FCN (NEQ, X, Y, YPRIME), where
                NEQ    - Number of equations.  (Input)
                X      - Independent variable.  (Input)
                Y      - Array of length NEQ containing the dependent
                         variable values.  (Input)
                YPRIME - Array of length NEQ containing the values of
                         the right hand side of the system of
                         equations.  (Output)
                         See remark 3.
                FCN must be declared EXTERNAL in the calling program.
       N      - Number of differential equations.  (Input)
       PW     - Work vector of length NPW.  PW is used both to store the
                Jacobian and as workspace.
                See IVPAG for the value of NPW.
       SAVE1  - Work vector of length N.
       SAVE2  - Work vector of length N.
       Y      - Work array of size N by IYCOL, whose first column holds
                the present solution, and whose other columns are used as
                workspace.  If METH is 1, then IYCOL = 13.  If METH is 2,
                then IYCOL = 6.
       YMAX   - Vector of length N containing the maximum Y values
                computed so far.  (Output)


    Chapter:    MATH/LIBRARY Differential Equations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */

#ifdef ANSI
static void l_i5pag(void (*fcn) (Mint, Mfloat, Mfloat*, Mfloat*), Mint *n, 
		    Mfloat *pw, Mfloat save1[],	Mfloat save2[], Mfloat *y, 
		    Mfloat ymax[])
#else
static void l_i5pag(fcn, n, pw, save1, save2, y, ymax)
	void            (*fcn) ();
	Mint            *n;
	Mfloat          *pw, save1[], save2[], *y, ymax[];
#endif
{
#define PW(I_,J_)	(pw+(I_)*(*n)+(J_))
#define Y(I_,J_)	(y+(I_)*(*n)+(J_))
	static Mint             i, j;
	static Mfloat           r, r0, yi;

	/*
	 * Finite difference full matrix Jacobian
	 */
	i11ag.nfcn += *n;
	r0 = fabs(i11ag.h) * imsl_snrm2(*n, save2, 1) * 1000.0 * i11ag.uround;
	for (j = 1; j <= *n; j++) {
		yi = *Y(0, j - 1);
		r = imsl_f_max(i11ag.epsj * imsl_f_max(ymax[j - 1], fabs(*Y(0, j - 1))),
			       r0);
		*Y(0, j - 1) += r;
		(*fcn) (*n, i11ag.t, y, save1);
		for (i = 1; i <= *n; i++) {
			*PW(j - 1, i - 1) = (save1[i - 1] - save2[i - 1]) / r;
		}
		*Y(0, j - 1) = yi;
	}
L_9000:
	;
	return;
}				/* end of function */

#undef PW
#undef Y

/* Structured by FOR_STRUCT, v0.2, on 10/08/90 at 17:25:07
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  I6PAG/DI6PAG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Compute the Jacobian, in band storage mode, using
                finite differences.

    Usage:      CALL I6PAG (FCN, N, PW, LDPW, SAVE1, SAVE2, Y, YMAX, WK)

    Arguments:
       FCN    - User-supplied SUBROUTINE to evaluate functions.
                The usage is
                CALL FCN (NEQ, X, Y, YPRIME), where
                NEQ    - Number of equations.  (Input)
                X      - Independent variable.  (Input)
                Y      - Array of length NEQ containing the dependent
                         variable values.  (Input)
                YPRIME - Array of length NEQ containing the values of
                         the right hand side of the system of
                         equations.  (Output)
                         See remark 3.
                FCN must be declared EXTERNAL in the calling program.
       N      - Number of differential equations.  (Input)
       PW     - Work vector of length NPW.  PW is used both to store the
                Jacobian and as workspace.
                See IVPAG for the value of NPW.
       LDPW   - NUCPW + NLCPW + 1, the leading dimension of the band
                matrix formed from PW.  (Input)
       SAVE1  - Work vector of length N.
       SAVE2  - Work vector of length N.
       Y      - Work array of size N by IYCOL, whose first column holds
                the present solution, and whose other columns are used as
                workspace.  If METH is 1, then IYCOL = 13.  If METH is 2,
                then IYCOL = 6.
       YMAX   - Vector of length N containing the maximum Y values
                computed so far.  (Output)
       WK     - Work vector of length N.

    Chapter:    MATH/LIBRARY Differential Equations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */

#ifdef ANSI
static void l_i6pag(void (*fcn) (Mint, Mfloat, Mfloat*, Mfloat*), Mint *n, 
		    Mfloat *pw, Mint *ldpw, Mfloat save1[], Mfloat save2[], 
		    Mfloat *y, Mfloat ymax[], Mfloat wk[])
#else
static void l_i6pag(fcn, n, pw, ldpw, save1, save2, y, ymax, wk)
	void            (*fcn) ();
	Mint            *n;
	Mfloat          *pw;
	Mint            *ldpw;
	Mfloat           save1[], save2[], *y, ymax[], wk[];
#endif
{
#define PW(I_,J_)	(pw+(I_)*(*ldpw)+(J_))
#define Y(I_,J_)	(y+(I_)*(*n)+(J_))
	static Mint             _d_l, _d_m, _do0, _do1, _do2, _do3, i, j, kband,
	                nb;
	static Mfloat           r, r0, yi;

	/*
	 * Finite difference band matrix Jacobian
	 */
	r0 = fabs(i11ag.h) * imsl_snrm2(*n, save2, 1) * 1000.0 * i11ag.uround;
	if (i11ag.mtype == 3) {
		nb = 2 * i11ag.nuc + 1;
	} else {
		nb = *ldpw;
	}
	if (nb > *n) {
		/* Dense approximate Jacobian */
		i11ag.nfcn += *n;
		for (j = 1; j <= *n; j++) {
			yi = *Y(0, j - 1);
			r = imsl_f_max(i11ag.epsj * imsl_f_max(ymax[j - 1], fabs(*Y(0, j - 1))),
				       r0);
			*Y(0, j - 1) = yi + r;
			(*fcn) (*n, i11ag.t, y, save1);
			for (i = imsl_i_max(1, j - i11ag.nuc); i <= imsl_i_min(*n, j + i11ag.nlc); i++) {
				*PW(j - 1, i - j + i11ag.nuc) = (save1[i - 1] - save2[i - 1]) /
					r;
			}
			*Y(0, j - 1) = yi;
		}
	} else {
		/* Banded approximate Jacobian */
		i11ag.nfcn += nb;
		for (kband = 1; kband <= nb; kband++) {
			for (j = kband, _do0 = DOCNT(kband, *n, _do1 = nb); _do0 > 0; j += _do1, _do0--) {
				wk[j - 1] = *Y(0, j - 1);
				r = imsl_f_max(i11ag.epsj * imsl_f_max(ymax[j - 1], fabs(*Y(0, j - 1))),
					       r0);
				*Y(0, j - 1) = wk[j - 1] + r;
			}
			(*fcn) (*n, i11ag.t, y, save1);
			for (j = kband, _do2 = DOCNT(kband, *n, _do3 = nb); _do2 > 0; j += _do3, _do2--) {
				r = *Y(0, j - 1) - wk[j - 1];
				for (i = imsl_i_max(j - i11ag.nuc, 1); i <= imsl_i_min(j + i11ag.nlc,
								 *n); i++) {
					*PW(j - 1, i - j + i11ag.nuc) = (save1[i - 1] -
							  save2[i - 1]) / r;
				}
				*Y(0, j - 1) = wk[j - 1];
			}
		}
	}

L_9000:
	;
	return;
}				/* end of function */

#undef PW
#undef Y

/* Structured by FOR_STRUCT, v0.2, on 10/08/90 at 17:26:37
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  I7PAG/DI7PAG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Solve an initial value problem for ordinary differential
                equations using an Adams-Moulton or Gear method.

    Usage:      CALL I7PAG (FCN, N, PW, SAVE1, SAVE2, Y, YMAX, HL0, IER)

    Arguments:
       FCN    - User-supplied SUBROUTINE to evaluate functions.
                The usage is
                CALL FCN (NEQ, X, Y, YPRIME), where
                NEQ    - Number of equations.  (Input)
                X      - Independent variable.  (Input)
                Y      - Array of length NEQ containing the dependent
                         variable values.  (Input)
                YPRIME - Array of length NEQ containing the values of
                         the right hand side of the system of
                         equations.  (Output)
                         See remark 3.
                FCN must be declared EXTERNAL in the calling program.
       N      - Number of differential equations.  (Input)
       PW     - Work vector of length NPW.  PW is used both to store the
                Jacobian and as workspace.
                See IVPAG for the value of NPW.
       SAVE1  - Work vector of length N.
       SAVE2  - Work vector of length N.
       Y      - Work array of size N by IYCOL, whose first column holds
                the present solution, and whose other columns are used as
                workspace.  If METH is 1, then IYCOL = 13.  If METH is 2,
                then IYCOL = 6.
       YMAX   - Vector of length N containing the maximum Y values
                computed so far.  (Output)
       HL0    - Real scalar.  (Output)
       IER    - Error parameter.  (Output)

    Chapter:    MATH/LIBRARY Differential Equations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */

#ifdef ANSI
static void l_i7pag(void (*fcn) (Mint, Mfloat, Mfloat*, Mfloat*), Mint *n, 
		    Mfloat pw[], Mfloat save1[],
		    Mfloat save2[], Mfloat *y, Mfloat ymax[],
		    Mfloat *hl0, Mint *ier)
#else
static void l_i7pag(fcn, n, pw, save1, save2, y, ymax, hl0, ier)
	void            (*fcn) ();
	Mint            *n;
	Mfloat           pw[], save1[], save2[], *y, ymax[], *hl0;
	Mint            *ier;
#endif
{
#define Y(I_,J_)	(y+(I_)*(*n)+(J_))
	static Mint             i;
	static Mfloat           d, r0;

	/*
	 * Matrix used PW = I - H*EL(1)*D, where D is a diagonal matrix
	 */
	*ier = 0;
	for (i = 1; i <= *n; i++) {
		pw[i - 1] = *Y(0, i - 1) + i11ag.el[0] * .1 * (i11ag.h * save2[i - 1] -
							       *Y(1, i - 1));
	}
	(*fcn) (*n, i11ag.t, pw, save1);
	i11ag.nfcn += 1;
	*hl0 = i11ag.h * i11ag.el[0];
	for (i = 1; i <= *n; i++) {
		r0 = i11ag.h * save2[i - 1] - *Y(1, i - 1);
		d = .1 * r0 - i11ag.h * (save1[i - 1] - save2[i - 1]);
		if (fabs(r0) < i11ag.uround * ymax[i - 1]) {
			pw[i - 1] = F_ONE;
			save1[i - 1] = F_ZERO;
		} else {
			if (fabs(d) == F_ZERO) {
				*ier = 2;
				goto L_9000;
			}
			pw[i - 1] = .1 * r0 / d;
			save1[i - 1] = pw[i - 1] * r0;
		}
	}

L_9000:
	;
	return;
}				/* end of function */

#undef Y

/* Structured by FOR_STRUCT, v0.2, on 10/08/90 at 17:29:37
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  I8PAG/DI8PAG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Solve an initial value problem for ordinary differential
                equations using an Adams-Moulton or Gear method.

    Usage:      CALL I8PAG (FCN, N, PW, SAVE1, Y, HL0, IER)

    Arguments:
       FCN    - User-supplied SUBROUTINE to evaluate functions.
                The usage is
                CALL FCN (NEQ, X, Y, YPRIME), where
                NEQ    - Number of equations.  (Input)
                X      - Independent variable.  (Input)
                Y      - Array of length NEQ containing the dependent
                         variable values.  (Input)
                YPRIME - Array of length NEQ containing the values of
                         the right hand side of the system of
                         equations.  (Output)
                         See remark 3.
                FCN must be declared EXTERNAL in the calling program.
       N      - Number of differential equations.  (Input)
       PW     - Work vector of length NPW.  PW is used both to store the
                Jacobian and as workspace.
                See IVPAG for the value of NPW.
       SAVE1  - Work vector of length N.
       Y      - Work array of size N by IYCOL, whose first column holds
                the present solution, and whose other columns are used as
                workspace.  If METH is 1, then IYCOL = 13.  If METH is 2,
                then IYCOL = 6.
       HL0    - Real scalar.  (Output)
       IER    - Error parameter.  (Output)

    Chapter:    MATH/LIBRARY Differential Equations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */

#ifdef ANSI
static void l_i8pag(void (*fcn) (Mint, Mfloat, Mfloat*, Mfloat*), Mint *n, 
		    Mfloat pw[], Mfloat save1[], Mfloat *y, Mfloat *hl0, Mint *ier)
#else
static void l_i8pag(fcn, n, pw, save1, y, hl0, ier)
	void            (*fcn) ();
	Mint            *n;
	Mfloat           pw[], save1[], *y, *hl0;
	Mint            *ier;
#endif
{
#define Y(I_,J_)	(y+(I_)*(*n)+(J_))
	static Mint             i;
	static Mfloat           d, phl0, r;

	/*
	 * The coefficient H*EL(1) in P is updated.
	 */
	phl0 = *hl0;
	*hl0 = i11ag.h * i11ag.el[0];
	if (*hl0 != phl0) {
		r = *hl0 / phl0;
		for (i = 1; i <= *n; i++) {
			d = F_ONE - r * (F_ONE - F_ONE / pw[i - 1]);
			if (fabs(d) == F_ZERO) {
				*ier = 3;
				goto L_9000;
			}
			pw[i - 1] = F_ONE / d;
		}
	}
	for (i = 1; i <= *n; i++) {
		save1[i - 1] *= pw[i - 1];
	}

L_9000:
	;
	return;
}				/* end of function */

#undef Y

/* Structured by FOR_STRUCT, v0.2, on 10/08/90 at 17:30:54
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  I9PAG/DI9PAG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Solve an initial value problem for ordinary differential
                equations using an Adams-Moulton or Gear method.

    Usage:      CALL I9PAG (N, TOUT, Y, Y0)

    Arguments:
       N      - Number of differential equations.  (Input)
       TOUT   - Real scalar.  (Output)
       Y      - Work array of size N by IYCOL, whose first column holds
                the present solution, and whose other columns are used as
                workspace.  If METH is 1, then IYCOL = 13.  If METH is 2,
                then IYCOL = 6.
       Y0     - Vector of length N, holding the interpolated values of
                Y.  (Output)

    Chapter:    MATH/LIBRARY Differential Equations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */

#ifdef ANSI
static void l_i9pag(Mint *n, Mfloat *tout, Mfloat *y, Mfloat y0[])
#else
static void l_i9pag(n, tout, y, y0)
	Mint            *n;
	Mfloat          *tout, *y, y0[];
#endif
{
#define Y(I_,J_)	(y+(I_)*(*n)+(J_))
	static Mint             j, l;
	static Mfloat           s, s1;


	scopy(*n, y, 1, y0, 1);
	/*
	 * This subroutine computes interpolated values of the dependent
	 * variable Y and stores them in Y0. The interpolation is to the
	 * point T = TOUT, and uses the Nordsieck history array Y, as
	 * follows. NQ Y0(I)  =  SUM  Y(I,J+1)*S**J , J=0 where S =
	 * -(T-TOUT)/H.
	 */
	l = i11ag.jstart + 1;
	s = (*tout - i11ag.t) / i11ag.h;
	s1 = F_ONE;
	for (j = 2; j <= l; j++) {
		s1 *= s;
		saxpy(*n, s1, Y(j - 1, 0), 1, y0, 1);
	}
	return;
}				/* end of function */

#undef Y

/* Structured by FOR_STRUCT, v0.2, on 10/22/90 at 15:18:05
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  IVPAG/DIVPAG (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 14, 1986

    Purpose:    Solve an initial-value problem for ordinary differential
                equations using an Adams-Moulton or Gear method.

    Usage:      CALL IVPAG (IDO, NEQ, FCN, FCNJ, A, X, XEND, TOL,
                            PARAM, Y)

    Arguments:
       IDO    - Flag indicating the state of the computation.
                (Input/Output)
                1        Initial entry
                2        Normal reentry
                3        Final call to release workspace
                4        Return because of interrupt 1
                5        Return because of interrupt 2 with step accepted
                6        Return because of interrupt 2 with step rejected
                7        Return for new value of A.  The matrix A at X
                         must be recomputed and IVPAG/DIVPAG called
                         again.  No other argument (including IDO) should
                         be changed.  This value of IDO is returned only
                         if PARAM(19)=2.
                Normally, the initial call is made with IDO=1.  The
                routine then sets IDO=2 and this value is then used for
                all but the last call which is made with IDO=3.  This
                final call is only used to release workspace, which was
                automatically allocated by the initial call with IDO=1.
                See Remark 5 for a description of the interrupts.
       NEQ    - Number of differential equations.  (Input)
       FCN    - User-supplied SUBROUTINE to evaluate functions.
                The usage is
                CALL FCN (NEQ, X, Y, YPRIME), where
                NEQ    - Number of equations.  (Input)
                X      - Independent variable.  (Input)
                Y      - Array of length NEQ containing the dependent
                         variable values.  (Input)
                YPRIME - Array of length NEQ containing the values of
                         the right-hand side of the system of
                         equations.  (Output)
                         See Remark 3.
                FCN must be declared EXTERNAL in the calling program.
       FCNJ   - User-supplied SUBROUTINE to compute the Jacobian.
                The usage is
                CALL FCNJ (NEQ, X, Y, DYPDY), where
                NEQ    - Number of equations.  (Input)
                X      - Independent variable.  (Input)
                Y      - Array of length NEQ containing the dependent
                         variable values.  (Input)
                DYPDY  - Matrix (whose type is determined by PARAM(14) =
                         MTYPE), containing the derivatives
                         dYPRIME(i)/dY(j).  (Output)
                         YPRIME is the right-hand side of the ODE.
                FCNJ must be declared EXTERNAL in the calling program.
                If PARAM(19) = IATYPE is nonzero, then FCNJ should
                compute the Jacobian of the right-hand side of the
                equation Ay'=f(x,y).
                FCNJ is used only if PARAM(13) = MITER = 1.
       A      - Matrix used when ODE system is implicit.  (Input)
                A is referenced only if PARAM(19) = IATYPE is nonzero.
                Its data structure is determined by PARAM(14) = MTYPE.
                A must be nonsingular and MITER must be 1 or 2.
                See Remark 3.
       X      - Independent variable.  (Input/Output)
                On input, X supplies the initial value.
                On output, X is replaced by XEND unless error conditions
                arise.  See IDO for details.
       XEND   - Value of X at which the solution is desired.  (Input)
                XEND may be less than the initial value of X.
       TOL    - Tolerance for error control.  (Input)
                An attempt is made to control the norm of the local error
                such that the global error is proportional to TOL.
                More than one run, with different values of TOL, can be
                used to estimate the global error.
                Generally, it should not be greater than 0.001.
       PARAM  - Vector of length 50 containing optional parameters.
                (Input/Output)
                If a parameter is zero then the default value is used.
                The following parameters must be set by the user.
                   PARAM              Meaning
                     1   HINIT  - Initial value of the step size H.
                                  Always nonnegative.
                                  Default: ABS(.001*(XEND-X)).
                     2   HMIN   - Minimum value of the step size H.
                                  Default: 0.0.
                     3   HMAX   - Maximum value of the step size H.
                                  Default: No limit is imposed on the
                                  step size.
                     4   MXSTEP - Maximum number of steps allowed.
                                  Default: 500.
                     5   MXFCN  - Maximum number of function evaluations
                                  allowed.
                                  Default: No limit.
                     6   MAXORD - Maximum order of the method.
                                  Default: If Adams' method, 12.
                                  If Gear's method, 5.
                                  The defaults are also the maximum
                                  values allowed.
                     7   INTRP1 - If nonzero then the subroutine will
                                  return, with IDO=4, before every step.
                                  See Remark 5.
                                  Default: 0.
                     8   INTRP2 - If nonzero then the subroutine will
                                  return, with IDO=5, after every
                                  successful step and with IDO=6 after
                                  every unsuccessful step.
                                  See Remark 5.
                                  Default: 0
                     9   SCALE  - A measure of the scale of the problem,
                                  such as an approximation to the average
                                  value of a norm of the Jacobian along
                                  the trajectory.
                                  Default: 1.0
                    10   INORM  - Switch determining error norm.
                                  In the following Ei is the absolute
                                  value of an estimate of the error in
                                  Y(i), called Yi here.
                                  0 - imsl_i_min(absolute error, relative error)
                                      = imsl_i_max(Ei/Wi), i=1,2,...,NEQ, where
                                      Wi = imsl_i_max(abs(Yi,1.0),
                                  1 - absolute error = imsl_i_max(Ei), i=1,2,...
                                  2 - imsl_i_max(Ei/Wi), i=1,2,..., where
                                      Wi = imsl_i_max(abs(Yi),FLOOR),
                                      and FLOOR is PARAM(11).
                                  3 - Euclidean norm scaled by YMAX
                                      =sqrt(sum(Ei**2/Wi**2)), where
                                      Wi = imsl_i_max(abs(Yi),1.0); for YMAX,
                                      see Remark 1.
                    11   FLOOR  - Used in the norm computation.
                                  Default: 1.0
                    12   METH   - Method indicator.
                                  METH=1 selects Adams' method
                                  METH=2 selects Gear's backward
                                  difference method.  Default: 1
                    13   MITER  - Iteration method indicator.
                                  MITER=0 selects functional iteration.
                                     IATYPE must be set to zero with this
                                     option.
                                  MITER=1 selects the chord method with
                                     an analytic, user-supplied Jacobian.
                                  MITER=2 selects the chord method with
                                     a finite-difference Jacobian.
                                  MITER=3 selects the chord method with
                                     the Jacobian replaced by a diagonal
                                     approximation based on a directional
                                     derivative.  IATYPE must be set to
                                     zero with this option.
                                  Default: 0
                    14   MTYPE  - Matrix type for A (if used) and the
                                  Jacobian (if MITER is 1 or 2).
                                  These must have the same type.
                                  MTYPE=0 selects full matrices.
                                  MTYPE=1 selects band matrices.
                                  MTYPE=2 selects symmetric positive
                                     definite matrices
                                  MTYPE=3 selects band symmetric
                                     positive definite matrices.
                                  Default: 0
                    15   NLC    - Number of lower codiagonals, used
                                  if MTYPE=1.  Default: 0
                    16   NUC    - Number of upper codiagonals, used
                                  if MTYPE is 1 or 3.  Default: 0.
                    17          - Not used.
                    18   EPSJ   - Epsilon used in computing finite-
                                  difference Jacobians.
                                  Default: SQRT(eps), where eps is
                                     the machine precision.
                    19   IATYPE - Type of the matrix A.
                                  IATYPE=0 implies A is not used (the
                                     ODE system is explicit).
                                  IATYPE=1 implies A is constant.
                                  IATYPE=2 implies A depends on X.
                                  Default: 0
                    20   LDA    - Leading dimension of A exactly as
                                  specified in the dimension statement
                                  in the calling program.  Used if
                                  IATYPE is not zero.
                                  Default: N,         if MTYPE=0 or 2
                                           NUC+NLC+1, if MTYPE=1
                                           NUC+1,     if MTYPE=3
                    21-30       - Not used.
                The following entries in PARAM are set by the program.
                    31   HTRIAL - Current trial step size.
                    33   HMINC  - Computed minimum step size.
                    33   HMAXC  - Computed maximum step size.
                    34   NSTEP  - Number of steps taken.
                    35   NFCN   - Number of function evaluations used.
                    36   NJE    - Number of Jacobian evaluations.
                    37-50       - Not used.
       Y      - Vector of length NEQ of dependent variables.
                (Input/Output)
                On input, Y contains the initial values.  On output,
                Y contains the approximate solution.

    Remarks:
    1. Automatic workspace usage is
                IVPAG    4*NEQ + NMETH + NPW + NIPVT, or
                DIVPAG   8*NEQ + 2*NMETH + 2*NPW + NIPVT units.
       Here
       NMETH = 13*NEQ if METH is 1,
       NMETH =  6*NEQ if METH is 2.
       NPW = 2*NEQ + NPWM + NPWA
       where
       NPWM = 0 if MITER is 0 or 3,
       NPWM = NEQ**2 if MITER is 1 or 2, and if MTYPE is 0 or 2.
       NPWM = NEQ*(2*NLC+NUC+1) if MITER is 1 or 2 and MTYPE=1.
       NPWM = NEQ*(NLC+1) if MITER is 1 or 2 and if MTYPE=3.
       NPWA = 0 if IATYPE is 0.
       NPWA = NEQ**2 if IATYPE is nonzero and MTYPE=0,
       NPWA = NEQ*(2*NLC+NUC+1) if IATYPE is nonzero and MTYPE=1
       NIPVT = NEQ if MITER is 1 or 2 and MTYPE is 0 or 1,
       NIPVT = 1, otherwise.
       Workspace may be explicitly provided, if desired, by use of
       I2PAG/DI2PAG.  The reference is
                CALL I2PAG (IDO, NEQ, FCN, FCNJ, A, X, XEND, TOL, PARAM,
                            Y, YTEMP, YMAX, ERROR, SAVE1, SAVE2, PW,
                            IPVT, VNORM)
       None of the additional array arguments should be changed from the
       first call with IDO=1 until after the final call with IDO=3.
       The additional arguments are as follows:
       YTEMP  - Vector of length NMETH.  (Workspace)
       YMAX   - Vector of length NEQ containing the maximum Y-values
                computed so far.  (Output)
       ERROR  - Vector of length NEQ containing error estimates for each
                component of Y.  (Output)
       SAVE1  - Vector of length NEQ.  (Workspace)
       SAVE2  - Vector of length NEQ.  (Workspace)
       PW     - Vector of length NPW.  PW is used both to store the
                Jacobian and as workspace.  (Workspace)
       IPVT   - Vector of length NEQ.  (Workspace)
       VNORM  - User-supplied SUBROUTINE to compute the norm of the
                error.  (Input)
                The routine may be provided by the user, or the IMSL
                routine I3PRK/DI3PRK may be used.
                The usage is
                CALL VNORM (NEQ, V, Y, YMAX, ENORM), where
                NEQ    - Number of equations.  (Input)
                V      - Vector of length NEQ containing the vector whose
                         norm is to be computed.  (Input)
                Y      - Vector of length NEQ containing the values of
                         the dependent variable.  (Input)
                YMAX   - Vector of length NEQ containing the maximum Y
                         values computed so far.  (Input)
                ENORM  - Norm of the vector V.  (Output)
                VNORM must be declared EXTERNAL in the calling program.

    2. Informational errors
  *     Type*//* Code 4   1  After some initial success, the integration was
  * halted by repeated error-test failures. 4   2  The maximum number of
  * function evaluations have been used. 4   3  The maximum number of steps
  * allowed have been used. The problem may be stiff. 4   4  On the next step
  * X+H will equal X.  Either TOL is too small or the problem is stiff. 4   5
  * After some initial success, the integration was halted by a test on TOL.
  * 4   6  Integration was halted after failing to pass the error test, even
  * after reducing the step size by a factor of 1.0E+10.  TOL may be too
  * small. 4   7  Integration was halted after failing to achieve corrector
  * convergence, even after reducing the step size by a factor of 1.0E+10.
  * TOL may be too small. 4   8  IATYPE is nonzero and the input matrix A is
  * singular.
  * 
  * 3. Both explicit ODE systems, of the form y'=f(x,y), and implicit ODE
  * systems, of the form Ay'=f(x,y) can be solved.  If the system is explicit
  * then PARAM(19)=0 and the matrix A is not referenced.  If the system is
  * implicit then PARAM(14) determines the data structure of the matrix A.
  * If PARAM(19)=1 then A is assumed to be a constant matrix (not depending
  * on x or y).  The value of A used on the first call (with IDO=1) is used
  * until after a call with IDO=3.  The value of A must not be changed
  * between these calls. If PARAM(19)=2 then the matrix is assumed to be a
  * function of x.
  * 
  * 4. If MTYPE is greater than zero, then MITER must equal 1 or 2.
  * 
  * 5. If PARAM(7) is nonzero, the subroutine returns with IDO = 4, and will
  * resume calculation at the point of interruption if reentered with IDO =
  * 4.  If PARAM(8) is nonzero, the subroutine will interrupt the
  * calculations immediately after it decides whether or not to accept the
  * result of the most recent trial step.  IDO = 5 if the routine plans to
  * accept, or IDO = 6 if it plans to reject.  IDO may be changed by the user
  * in order to force acceptance of a step (by changing IDO from 6 to 5) that
  * would otherwise be rejected, or vice versa. Relevant parameters to
  * observe after return from an interrupt are IDO, HTRIAL, NSTEP, NFCN, NJE,
  * and Y.  Y is the newly computed trial value, accepted or not.
  * 
  * Keywords:   ODE; Ordinary differential equation; Stiff; First order;
  * Predictor-corrector; Multi-step; Backward differentiation
  * 
  * GAMS:       I1a1b
  * 
  * Chapter:    MATH/LIBRARY Differential Equations
  * 
  * Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
  * 
  * Warranty:   IMSL warrants only that IMSL testing has been applied to this
  * code.  No other warranty, expressed or implied, is applicable.
  * 
  * ----------------------------------------------------------------------- */
#ifdef ANSI
void imsl_f_ode_adams_gear(Mint neq, Mfloat *x, Mfloat xend,
			Mfloat *y, Mchar *state, 
			void (*fcn) (Mint, Mfloat, Mfloat*, Mfloat*))
#else
void imsl_f_ode_adams_gear(neq, x, xend, y, state, fcn)
	Mint             neq;
	void            (*fcn) ();
	Mfloat           xend, y[];
	Mfloat 	*x;
	Mchar		*state;
#endif
{
	static Mfloat		xend_float;
	static Mint      _l0, _l1, i, iatype,
	                 k, meth, miter,
	                mtype, nipvt, nlc, nmeth, npw, nuc;
	static Mfloat *iymax, *ierror, *isave1, *isave2, *iytemp, *ipw, *iw, *a;
	static Mint *iipvt;


#ifdef DOUBLE
	imsl_e1psh("imsl_d_ode_adams_gear");
#else
	imsl_e1psh("imsl_f_ode_adams_gear");
#endif

	lv_state = (Env *) state;

	/* Check IDO */
	if (lv_state->ido == 1) {
		/* Check NEQ */
		if (neq < 1) {
			imsl_e1sti(1, neq);

          		imsl_ermes(IMSL_TERMINAL, IMSL_ODE_NEG_NEQ);
			goto L_9000;
		}

		iatype = 0;
		meth = lv_state->method;
		miter = (lv_state->miter) - 1;
		mtype = 0;
		if (meth == 1) {
			nmeth = 13 * neq;
		} else if (meth == 2) {
			nmeth = 6 * neq;
		}
		if (miter == 1 || miter == 2) {
			if (mtype == 0) {
				npw = neq * (neq + 2);
				if (iatype != 0)
					npw += neq * neq;
				nipvt = neq;
			}
		} else if (miter == 0) {
			npw = neq;
			nipvt = 0;
		} else if (miter == 3) {
			npw = 2 * neq;
			nipvt = 0;
		}
		/* Allocate workspace */
/* 
		iw = imsl_i1kgt(ADR(_l0, 4 ** neq + nmeth + npw), ADR(_l1, 3));
		iipvt = imsl_i1kgt(ADR(_l0, imsl_i_max(nipvt, 1)), ADR(_l1, 2));
*/

		if (lv_state->ido == 1) {
		lv_state->wk = (Mfloat *) imsl_malloc((4*neq+nmeth+npw)
		*sizeof(*lv_state->wk));
		lv_state->iwk = (Mint *) imsl_malloc(imsl_i_max(nipvt,1)
		*sizeof(*lv_state->iwk));
		}

/*
		if (imsl_n1rty(1) != 0) {
			imsl_e1mes(5, 5, " ");
			imsl_e1sti(1, neq);
			imsl_e1sti(2, nmeth);
			imsl_e1sti(3, npw);

			imsl_ermes(IMSL_TERMINAL,
                        IMSL_WORKSPACE_REQUIREMENT);
			goto L_9000;
		}
*/
		if (lv_state->ido == 1) {
			lv_state->iymax = lv_state->wk;
			lv_state->ierror = lv_state->iymax + neq;
			lv_state->isave1 = lv_state->ierror + neq;
			lv_state->isave2 = lv_state->isave1 + neq;
			lv_state->iytemp = lv_state->isave2 + neq;
			lv_state->ipw = lv_state->iytemp + nmeth;
		}
	}

	if (lv_state->hinit == 0.0) 
		lv_state->hinit = .001*fabs(xend - *x);
	xend_float = xend;

	l_i2pag(&(lv_state->ido), &neq, fcn, lv_state->fcnj, a, x, 
			&xend_float, 
			&(lv_state->tol), y, lv_state->iytemp,
	       		lv_state->iymax, lv_state->ierror, lv_state->isave1,
		   	lv_state->isave2, lv_state->ipw, lv_state->iwk, 
			l_i3prk);

L_9000:
	;
#ifdef DOUBLE
	imsl_e1pop("imsl_d_ode_adams_gear");
#else
	imsl_e1pop("imsl_f_ode_adams_gear");
#endif
	return;
}				/* end of function */
/* Structured by FOR_STRUCT, v0.2, on 10/24/90 at 14:52:06
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  I3PRK/DI3PRK (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 29, 1985

    Purpose:    Solve an initial value problem for ordinary differential
                equations using the Runge-Kutta-Verner fifth and sixth
                order method.

    Usage:      CALL I3PRK (NEQ, V, Y, YMAX, ENORM)

    Arguments:
       NEQ    - Number of equations.  (Input)
       V      - Vector of length NEQ containing the vector whose
                norm is to be computed.  (Input)
       Y      - Vector of length NEQ containing the values of the
                dependent variable.  (Input)
       YMAX   - Vector of length NEQ containing the maximum Y values
                computed so far.  (Input)
       ENORM  - Norm of the vector V.  (Output)

    Chapter:    MATH/LIBRARY Integration and Differentiation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_i3prk(Mint neq, Mfloat v[], Mfloat y[], Mfloat ymax[],
			Mfloat *enorm)
#else
static void l_i3prk(neq, v, y, ymax, enorm)
	Mint             neq;
	Mfloat           v[], y[], ymax[], *enorm;
#endif
{
	Mint             k;
	Mfloat           weight;


	if (i5prk.inorm == 0) {
		/* imsl_i_max (absolute, relative) */
		*enorm = F_ZERO;
		for (k = 1; k <= neq; k++) {
			weight = imsl_f_max(F_ONE, fabs(y[k - 1]));
			*enorm = imsl_f_max(*enorm, fabs(v[k - 1]) / weight);
		}
		/* Absolute error control */
	} else if (i5prk.inorm == 1) {
		*enorm = fabs(v[imsl_isamax(neq, v, 1) - 1]);
		/* Relative error control */
	} else if (i5prk.inorm == 2) {
		*enorm = F_ZERO;
		for (k = 1; k <= neq; k++) {
			weight = imsl_f_max(fabs(y[k - 1]), i5prk.floor);
			*enorm = imsl_f_max(*enorm, fabs(v[k - 1]) / weight);
		}
	} else if (i5prk.inorm == 3) {
		/* Same as DGEAR's error control */
		*enorm = F_ZERO;
		for (k = 1; k <= neq; k++) {
			weight = imsl_f_max(F_ONE, fabs(ymax[k - 1]));
			*enorm += imsl_fi_power(v[k - 1] / weight, 2);
		}
		*enorm = sqrt(*enorm);
	}
L_9000:
	;
	return;
}				/* end of function */
/* Structured by FOR_STRUCT, v0.2, on 10/24/90 at 15:31:05
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  SSBMV (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    September 26, 1989

    Purpose:    Perform the matrix-vector operation y = alpha*A*x +
                imsl_beta*y, where A is a symmetric matrix in band symmetric
                storage mode.

    Usage:      CALL SSBMV (UPLO, N, NCODA, ALPHA, A, LDA, X, INCX,
                            BETA, Y, INCY )

    Arguments:
       UPLO   - Character specifing the storage structure.
                (Input)
                   UPLO              Structure
                'U' or 'u'      Symmetric matrix A is referenced using
                                its upper triangular part.
                'L' or 'l'      Symmetric matrix A is referenced using
                                its lower triangular part.
       N      - Order of the matrix A.  (Input)
       NCODA  - Number of codiagonals of A.  (Input)
       ALPHA  - Scalar multiplier for the matrix-vector product.
                (Input)
       A      - LDA by N array containing the matrix of order N stored
                in symmetric band form.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling routine.  (Input)
                LDA must be greater than or equal to NCODA + 1.
       X      - Vector of length (N-1)*IABS(INCX)+1.  (Input)
       INCX   - Displacement between elements of X.  (Input)
       BETA   - Scalar multiplier for Y.  (Input)
                When BETA is 0.0E0, Y is not referenced.  In that case,
                BETA*Y is defined as the zero vector.
       Y      - Vector of length (N-1)*IABS(INCY)+1.
                (Input/Output)
                On output, Y is replaced by the updated Y.
       INCY   - Displacement between elements of Y.  (Input)

    GAMS:       D1b

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_ssbmv(Mchar *uplo, unsigned uplo_s, Mint *n, Mint *ncoda,
			Mfloat *alpha, Mfloat *a, Mint *lda, Mfloat x[],
			Mint *incx, Mfloat *imsl_beta, Mfloat y[],
			Mint *incy)
#else
static void l_ssbmv(uplo, uplo_s, n, ncoda, alpha, a, lda, x,
	   incx, imsl_beta, y, incy)
	Mchar           *uplo;
	unsigned        uplo_s;
	Mint            *n, *ncoda;
	Mfloat          *alpha, *a;
	Mint            *lda;
	Mfloat           x[];
	Mint            *incx;
	Mfloat          *imsl_beta, y[];
	Mint            *incy;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
	Mint             lower, upper;
	Mint             ix, iy, j, jx, jy, m, ncodax, ncoday;
	Mfloat           temp1, temp2;


	upper = imsl_l1ame(uplo, uplo_s, "U", sizeof("U"));
	lower = imsl_l1ame(uplo, uplo_s, "L", sizeof("L"));
	/*
	 * Test the input parameters.
	 */
	if (*n < 0) {
		imsl_e1psh("SSBMV ");
		imsl_e1sti(1, *n);

		imsl_ermes(IMSL_TERMINAL, IMSL_NEED_N_GE_ZERO);
		imsl_e1pop("SSBMV ");
		goto L_9000;
	} else if (((*ncoda < 0) || (*ncoda > *n - 1)) && (*n > 0)) {
		imsl_e1psh("SSBMV ");
		imsl_e1sti(1, *ncoda);
		imsl_e1sti(2, *n);

		imsl_ermes(IMSL_TERMINAL, IMSL_BAD_NCODA_VALUE_GIVEN);
		imsl_e1pop("SSBMV ");
		goto L_9000;
	} else if (*lda < *ncoda + 1) {
		imsl_e1psh("SSBMV ");
		imsl_e1sti(1, *lda);
		imsl_e1sti(2, *ncoda);

		imsl_ermes(IMSL_TERMINAL, IMSL_NEED_LDA_GT_NCODA);
		imsl_e1pop("SSBMV ");
		goto L_9000;
	} else if (*incx == 0) {
		imsl_e1psh("SSBMV ");
		imsl_e1sti(1, *incx);

		imsl_ermes(IMSL_TERMINAL, IMSL_INCX_EQUALS_ZERO);
		imsl_e1pop("SSBMV ");
		goto L_9000;
	} else if (*incy == 0) {
		imsl_e1psh("SSBMV ");
		imsl_e1sti(1, *incy);

		imsl_ermes(IMSL_TERMINAL, IMSL_INCY_EQUALS_ZERO);
		imsl_e1pop("SSBMV ");
		goto L_9000;
	} else if ((!upper) && (!lower)) {
		imsl_e1psh("SSBMV ");
		imsl_e1stl(1, uplo);

		imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_UPLO_VALUE);
		imsl_e1pop("SSBMV ");
		goto L_9000;
	}
	/* Quick return if possible. */
	if (*n == 0 || (*alpha == F_ZERO && *imsl_beta == F_ONE))
		goto L_9000;
	/* First form Y = BETA*Y */
	if (*imsl_beta == F_ZERO) {
		sset(*n, F_ZERO, y, abs(*incy));
	} else if (*imsl_beta != F_ONE) {
		sscal(*n, *imsl_beta, y, abs(*incy));
	}
	/* Set up the start points in X and Y */
	ncodax = 1;
	ncoday = 1;
	if (*incx < 0)
		ncodax = 1 - (*n - 1) ** incx;
	if (*incy < 0)
		ncoday = 1 - (*n - 1) ** incy;

	if (*alpha == F_ZERO)
		goto L_9000;

	if (upper) {
		/*
		 * Form Y when upper triangle of A is stored.
		 */
		ix = ncodax - *incx;
		for (j = 1; j <= *n; j++) {
			temp1 = *alpha * x[ix + *incx - 1];
			m = imsl_i_max(*ncoda + 1 - j, 0);
			ix = ncodax + (*ncoda - m) ** incx;
			iy = ncoday + (*ncoda - m) ** incy;
			jx = ncodax + (*ncoda - m - 1) * imsl_i_min(*incx, 0);
			jy = ncoday + (*ncoda - m - 1) * imsl_i_min(*incy, 0);
			saxpy(*ncoda - m, temp1, A(j - 1, m), 1, &y[jy - 1], *incy);
			temp2 = imsl_sdot(*ncoda - m, A(j - 1, m), 1, &x[jx - 1], *incx);
			y[iy - 1] += temp1 ** A(j - 1, *ncoda) + *alpha * temp2;
			if (j > *ncoda) {
				ncodax += *incx;
				ncoday += *incy;
			}
		}
	} else {
		/*
		 * Form Y when lower triangle of A is stored.
		 */
		jx = ncodax;
		jy = ncoday;
		for (j = 1; j <= *n; j++) {
			temp1 = *alpha * x[jx - 1];
			m = imsl_i_min(*ncoda, *n - j);
			ix = jx + (m - 1) * imsl_i_min(*incx, 0) + *incx;
			iy = jy + (m - 1) * imsl_i_min(*incy, 0) + *incy;
			y[jy - 1] += temp1 ** A(j - 1, 0);
			saxpy(m, temp1, A(j - 1, 1), 1, &y[iy - 1], *incy);
			temp2 = imsl_sdot(m, A(j - 1, 1), 1, &x[ix - 1], *incx);
			y[jy - 1] += *alpha * temp2;
			jx += *incx;
			jy += *incy;
		}
	}
L_9000:
	return;
}				/* end of function */

#undef A

/* Structured by FOR_STRUCT, v0.2, on 10/24/90 at 15:31:59
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  SGBMV (Single precision version)

    Computer:   FORC/SINGLE

    Revised:    September 26, 1989

    Purpose:    Perform one of the matrix-vector operations:
                    y = alpha*A*x + imsl_beta*y
                 or
                    y = alpha*trans(A)*x + imsl_beta*y
                where A is a matrix stored in band storage mode, and
                trans(A) is the transpose of the matrix.

    Usage:      CALL SGBMV (TRANS, M, N, NLCA, NUCA, ALPHA, A, LDA,
                            X, INCX, BETA, Y, INCY )

    Arguments:
       TRANS  - Character specifying if the transpose product is to be
                computed.  (Input)
                   TRANS              Meaning
                'N' or 'n'      Compute y = alpha*A*x + imsl_beta*y
                'T' or 't'      Compute x = alpha*trans(A)*x + imsl_beta*y
                'C' or 'c'      Compute x = alpha*trans(A)*x + imsl_beta*y
       M      - Number of rows in A.  (Input)
       N      - Number of columns in A.  (Input)
       NLCA   - Number of lower diagonals of A.  (Input)
       NUCA   - Number of upper diagonals of A.  (Input)
       ALPHA  - Scalar multiplier for a matrix-vector product.  (Input)
       A      - M by N band matrix stored in band storage mode.
                (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling routine.  (Input)
       X      - Vector of length (N-1)*IABS(INCX)+1.  (Input)
       INCX   - Displacement between elements of X.  (Input)
       BETA   - Scalar multiplier for Y.  (Input)
                When BETA is 0.0E0, Y is not referenced.  In that case,
                BETA*Y is defined as the zero vector.
       Y      - Vector of length (N-1)*IABS(INCY)+1.
                (Input/Output)
                On output, Y is replaced by the updated Y.
       INCY   - Displacement between elements of Y.  (Input)

    GAMS:       D1b

    Chapter:    MATH/LIBRARY Basic Matrix/Vector Operations

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_sgbmv(Mchar *trans, unsigned trans_s, Mint *m, Mint *n,
			Mint *nlca, Mint *nuca, Mfloat *alpha,
			Mfloat *a, Mint *lda, Mfloat x[], Mint *incx,
			Mfloat *imsl_beta, Mfloat y[], Mint *incy)
#else
static void l_sgbmv(trans, trans_s, m, n, nlca, nuca, alpha, a,
	   lda, x, incx, imsl_beta, y, incy)
	Mchar           *trans;
	unsigned        trans_s;
	Mint            *m, *n, *nlca, *nuca;
	Mfloat          *alpha, *a;
	Mint            *lda;
	Mfloat           x[];
	Mint            *incx;
	Mfloat          *imsl_beta, y[];
	Mint            *incy;
#endif
{
#define A(I_,J_)	(a+(I_)*(*lda)+(J_))
	Mint            imsl_ctran, ntran, tran;
	Mint            i, i1, i2, ix, iy, j, kx, lenx, leny, ndiags, nlow,
	                nup, nx;



	ntran = imsl_l1ame(trans, trans_s, "N", sizeof("N"));
	tran = imsl_l1ame(trans, trans_s, "T", sizeof("T"));
	imsl_ctran = imsl_l1ame(trans, trans_s, "C", sizeof("C"));
	/*
	 * Test the input parameters.
	 */
	if (*m < 0) {
		imsl_e1psh("SGBMV ");
		imsl_e1sti(1, *m);

		imsl_ermes(IMSL_TERMINAL, IMSL_NEED_M_GE_ZERO);
		imsl_e1pop("SGBMV ");
		goto L_9000;
	} else if (*n < 0) {
		imsl_e1psh("SGBMV ");
		imsl_e1sti(1, *n);

		imsl_ermes(IMSL_TERMINAL, IMSL_NEED_N_GE_ZERO);
		imsl_e1pop("SGBMV ");
		goto L_9000;
	} else if (((*nlca < 0) || (*nlca > *m - 1)) && (*m > 0)) {
		imsl_e1psh("SGBMV ");
		imsl_e1sti(1, *nlca);
		imsl_e1sti(2, *m);

		imsl_ermes(IMSL_TERMINAL, IMSL_NEED_ZERO_LE_NLCA_LT_M);
		imsl_e1pop("SGBMV ");
		goto L_9000;
	} else if (((*nuca < 0) || (*nuca > *n - 1)) && (*n > 0)) {
		imsl_e1psh("SGBMV ");
		imsl_e1sti(1, *nuca);
		imsl_e1sti(2, *n);

		imsl_ermes(IMSL_TERMINAL, IMSL_WANT_ZERO_LE_NUCA_LT_N);
		imsl_e1pop("SGBMV ");
		goto L_9000;
	} else if (*lda < *nlca + *nuca + 1) {
		imsl_e1psh("SGBMV ");
		imsl_e1sti(1, *lda);
		imsl_e1sti(2, *nlca);
		imsl_e1sti(3, *nuca);

		imsl_ermes(IMSL_TERMINAL, IMSL_INCREASE_LDA_VALUE);
		imsl_e1pop("SGBMV ");
		goto L_9000;
	} else if (*incx == 0) {
		imsl_e1psh("SGBMV ");
		imsl_e1sti(1, *incx);

		imsl_ermes(IMSL_TERMINAL, IMSL_INCX_EQUALS_ZERO);
		imsl_e1pop("SGBMV ");
		goto L_9000;
	} else if (*incy == 0) {
		imsl_e1psh("SGBMV ");
		imsl_e1sti(1, *incy);

		imsl_ermes(IMSL_TERMINAL, IMSL_INCY_EQUALS_ZERO);
		imsl_e1pop("SGBMV ");
		goto L_9000;
	} else if (((!ntran) && (!tran)) && (!imsl_ctran)) {
		imsl_e1psh("SGBMV ");
		imsl_e1stl(1, trans);

		imsl_ermes(IMSL_TERMINAL, IMSL_TRANS_MUST_EQUAL_N_T_OR_C);
		imsl_e1pop("SGBMV ");
		goto L_9000;
	}
	/* Quick return if possible. */
	if ((*m == 0 || *n == 0) || (*alpha == F_ZERO && *imsl_beta == F_ONE))
		goto L_9000;
	/*
	 * Set LENX and LENY, the lengths of the vectors x and y.
	 */
	if (ntran) {
		lenx = *n;
		leny = *m;
	} else {
		lenx = *m;
		leny = *n;
	}

	ix = 1;
	iy = 1;
	if (*incx < 0)
		ix = (-lenx + 1) ** incx + 1;
	if (*incy < 0)
		iy = (-leny + 1) ** incy + 1;
	ndiags = *nuca + *nlca + 1;

	if (*imsl_beta == F_ZERO) {
		sset(leny, F_ZERO, y, abs(*incy));
	} else if (*imsl_beta != F_ONE) {
		sscal(leny, *imsl_beta, y, abs(*incy));
	}
	if (*alpha == F_ZERO)
		goto L_9000;
	/* Y = A*X */
	if (ntran) {
		/* Work down first column */
		for (i = 1; i <= (*nlca + 1); i++) {
			nup = imsl_i_min(*nuca, *n - i) + 1;
			nlow = imsl_i_min(*nlca, i - 1);
			kx = ix + (nup + nlow - 1) * imsl_i_min(*incx, 0);
			y[iy - 1] += *alpha * imsl_sdot(nup + nlow, A(0, *nuca + i - 1),
					       *lda - 1, &x[kx - 1], *incx);
			iy += *incy;
		}
		/* Work across bottom row */
		i = *nlca + 2;
		for (j = 2; j <= (*m - *nlca); j++) {
			nup = imsl_i_min(*nuca, *n - i) + 1;
			nlow = imsl_i_min(*nlca, i - 1);
			kx = ix + (j - 1) ** incx + (nup + nlow - 1) * imsl_i_min(*incx,
									 0);
			y[iy - 1] += *alpha * imsl_sdot(nup + nlow, A(j - 1, ndiags - 1),
					       *lda - 1, &x[kx - 1], *incx);
			i += 1;
			iy += *incy;
		}
		/* Y = trans(A)*X */
	} else {
		for (j = 1; j <= *n; j++) {
			i1 = imsl_i_max(1, j - *nuca);
			i2 = imsl_i_min(*m, j + *nlca);
			nx = i2 - i1 + 1;
			kx = ix + (i1 - 1) ** incx + (nx - 1) * imsl_i_min(*incx, 0);
			y[iy - 1] += *alpha * imsl_sdot(nx, A(j - 1, i1 - j + *nuca),
						      1, &x[kx - 1], *incx);
			iy += *incy;
		}
	}

L_9000:
	return;
}				/* end of function */

#undef A
