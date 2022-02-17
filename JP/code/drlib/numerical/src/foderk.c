#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

typedef void    PROTO((*Vnorm),(Mint,Mfloat*,Mfloat*,Mfloat*,Mfloat*));
static void	PROTO(l_i3prk,(Mint neq, Mfloat *v, Mfloat *y, Mfloat *ymax,
			Mfloat *enorm));
static void	l_init_state();
static VA_LIST_HACK	PROTO(l_ode_runge_kutta_mgr,(Mint ido, Mchar **state,
                        va_list argptr));

static Mfloat	lv_tiny = 0.0;
static Mfloat	lv_eps;

typedef struct {
    Mint	ido;
    Vnorm	vnorm;
    Mfloat	*wk;
    Mfloat	tol;
    Mint	inorm;
    Mfloat	hmin;
    Mfloat	hinit;
    Mfloat	scale;
    Mfloat	hmag;
    Mfloat	hmax;
    Mint	mxstep;
    Mint	nstep;
    Mint	mxfcn;
    Mint	nfcn;
    Mint	intrp1;
    Mint	intrp2;
    Mfloat	floor;
    Mint	*p_nstep;
    Mfloat	*p_htrial;
    Mfloat	*p_hmax;
    Mint	*p_nfcn;
    Mint	ixend;
    Mint	newfcn;
    Mint	nsucfl;
    Mint	nsucst;
    Mfloat      est;
    Mfloat	htrial;
    Mfloat	xendpv;
    Mfloat	xtrial;
} Env;

static Env  *lv_state;

static Mfloat	lv_rk[43] = {
	1.0 / 6.0,
	4.0 / 75.0,
	16.0 / 75.0,
	5.0 / 6.0,
	-8.0 / 3.0,
	5.0 / 2.0,
	-165.0 / 64.0,
	55.0 / 6.0,
	-425.0 / 64.0,
	85.0 / 96.0,
	12.0 / 5.0,
	-8.0,
	4015.0 / 612.0,
	-11.0 / 36.0,
	88.0 / 255.0,
	-8263.0 / 15000.0,
	124.0 / 75.0,
	-4501.0 / 4760.0,
	-81.0 / 250.0,
	2484.0 / 10625.0,
	3501.0 / 1720.0,
	-300.0 / 43.0,
	297275.0 / 52632.0,
	-319.0 / 2322.0,
	24068.0 / 84065.0,
	0.0,
	3850.0 / 26703.0,
	3.0 / 40.0,
	0.0,
	875.0 / 2244.0,
	23.0 / 72.0,
	264.0 / 1955.0,
	0.0,
	125.0 / 11592.0,
	43.0 / 616.0,
	1.0 / 160.0,
	0.0,
	125.0 / 17952.0,
	-1.0 / 144.0,
	84.0 / 13685.0,
	3.0 / 44.0,
	-125.0 / 11592.0,
	-43.0 / 616.0
};

/* Structured by FOR_STRUCT, v0.2, on 01/03/90 at 15:51:27
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  IVPRK/DIVPRK (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 14, 1986

    Purpose:    Solve an initial-value problem for ordinary differential
                equations using the Runge-Kutta-Verner fifth-order and
                sixth-order method.

    Usage:      CALL IVPRK (IDO, NEQ, FCN, X, xend, TOL, PARAM, Y)

    Arguments:
       IDO    - Flag indicating the state of the computation.
                (Input/Output)
                1        Initial entry
                2        Normal reentry
                3        Final call to release workspace
                4        Return because of interrupt 1
                5        Return because of interrupt 2 with step accepted
                6        Return because of interrupt 2 with step rejected
                Normally, the initial call is made with IDO=1.  The
                routine then sets IDO=2 and this value is then used for
                all but the last call which is made with IDO=3.  This
                final call is only used to release workspace, which was
                automatically allocated by the initial call with IDO=1.
                See Remark 3 for a description of the interrupts.
       NEQ    - Number of differential equations.  (Input)
       FCN    - User-supplied SUBROUTINE to evaluate functions.
                The usage is
                CALL FCN (NEQ, X, Y, YPRIME), where
                NEQ    - Number of equations.  (Input)
                X      - Independent variable.  (Input)
                Y      - Array of length NEQ containing the dependent
                         variable values.  (Input)
                YPRIME - Array of length NEQ containing the values of
                         dY/dX at (X,Y).  (Output)
                FCN must be declared EXTERNAL in the calling program.
       X      - Independent variable.  (Input/Output)
                On input, X supplies the initial value.
                On output, X is replaced by xend unless error conditions
                arise.  See IDO for details.
       xend   - Value of X at which the solution is desired.  (Input)
                xend may be less than the initial value of X.
       TOL    - Tolerance for error control.  (Input)
                An attempt is made to control the norm of the local error
                such that the global error is proportional to TOL.
                More than one run, with different values of TOL, can be
                used to estimate the global error.
                Generally, it should not be greater than 0.001.
       PARAM  - Vector of length 50 containing optional parameters.
                (Input/Output)
                If a parameter is zero then a default value is used.
                The following parameters must be set by the user.
                   PARAM              Meaning
                     1   HINIT  - Initial value of the step size H.
                                  Default: 10.0*MAX(AMACH(1),
                                  AMACH(4)*MAX(ABS(xend),ABS(X)))
                     2   HMIN   - Minimum value of the step size H.
                                  Default: 0.0
                     3   HMAX   - Maximum value of the step size H.
                                  Default: 2.0
                     4   MXSTEP - Maximum number of steps allowed.
                                  Default: 500
                     5   MXFCN  - Maximum number of function evaluations
                                  allowed.
                                  Default: No limit
                     6          - Not used.
                     7   INTRP1 - If nonzero then return with IDO=4,
                                  before each step.
                                  See Remark 3.
                                  Default: 0.
                     8   INTRP2 - If nonzero then return with IDO=5,
                                  after every successful step and with
                                  IDO=6 after every unsuccessful step.
                                  See Remark 3.
                                  Default: 0.
                     9   SCALE  - A measure of the scale of the problem,
                                  such as an approximation to the average
                                  value of a norm of the Jacobian along
                                  the trajectory.
                                  Default: 1.0
                    10   INORM  - Switch determining error norm.
                                  In the following Ei is the absolute
                                  value of an estimate of the error in
                                  Y(i), called Yi here.
                                  0 - min(absolute error, relative error)
                                      = max(Ei/Wi), i=1,2,...,NEQ, where
                                      Wi = max(abs(Yi,1.0),
                                  1 - absolute error = max(Ei), i=1,2,...
                                  2 - max(Ei/Wi), i=1,2,..., where
                                      Wi = max(abs(Yi),FLOOR),
                                      and FLOOR is PARAM(11).
                                  3 - Euclidean norm scaled by YMAX
                                      = sqrt(sum(Ei**2/Wi**2)), where
                                      Wi = max(abs(Yi),1.0); for YMAX,
                                      see Remark 1.
                    11   FLOOR  - Used in the norm computation.
                                  Default: 1.0
                    12-30       - Not used.
                The following entries in PARAM are set by the program.
                    31   HTRIAL - Current trial step size.
                    32   HMINC  - Computed minimum step size allowed.
                    33   HMAXC  - Computed maximum step size allowed.
                    34   NSTEP  - Number of steps taken.
                    35   NFCN   - Number of function evaluations used.
                    36-50       - Not used.
       Y      - Vector of length NEQ of dependent variables.
                (Input/Output)
                On input, Y contains the initial values.  On output,
                Y contains the approximate solution.

    Remarks:
    1. Automatic workspace usage is
                IVPRK    10*NEQ units, or
                DIVPRK   20*NEQ units.
       Workspace may be explicitly provided, if desired, by use of
       I2PRK/DI2PRK.  The reference is
                CALL I2PRK (IDO, NEQ, FCN, X, xend, TOL, PARAM, Y,
                            VNORM, WK)
       The additional arguments are as follows:
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
       WK     - Work array of length 10*NEQ.  WK must not be changed
                from the first call with IDO=1 until after the final call
                with IDO=3.

    2. Informational errors
       Type Code
         4   1  Cannot satisfy error condition.  TOL may be too small.
         4   2  Too many function evaluations needed.
         4   3  Too many steps needed.  The problem may be stiff.

    3. If PARAM(7) is nonzero, the subroutine returns with
       IDO = 4, and will resume calculation at the point of interruption
       if reentered with IDO = 4.  If PARAM(8) is nonzero, the
       subroutine will interrupt the calculations immediately after it
       decides whether or not to accept the result of the most
       recent trial step.  IDO = 5 if the routine plans to accept,
       or IDO = 6 if it plans to reject.  IDO may be changed by the user
       in order to force acceptance of a step (by changing IDO from 6
       to 5) that would otherwise be rejected, or vice versa.
       Relevant parameters to observe after return from an interrupt
       are IDO, HTRIAL, NSTEP, NFCN, and Y.  Y is the newly computed
       trial value, accepted or not.

    Keywords:   ODE; Ordinary differential equation; One-step;
                First order

    GAMS:       I1a1a

    Chapter:    MATH/LIBRARY Differential Equations

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */


#ifdef ANSI
void imsl_f_ode_runge_kutta(Mint neq,
			    Mfloat *x,
			    Mfloat xend, Mfloat *y, Mchar *state,
			    void (*fcn)(Mint,Mfloat,Mfloat*,Mfloat*))
#else
void imsl_f_ode_runge_kutta(neq, x, xend, y, state, fcn)
    Mint            neq;
    void            (*fcn) ();
    Mfloat          *x, xend, y[];
    Mchar	    *state;
#endif
{
#define WK(I_,J_)	(lv_state->wk+(I_)*(neq)+(J_))
	Mint		 i, k;
	Mfloat		 temp, temp1;

	E1PSH("imsl_f_ode_runge_kutta","imsl_d_ode_runge_kutta");

	    /* Check NEQ */
        if (neq < 1) {
	    imsl_e1sti(1, neq);
	    imsl_ermes(IMSL_TERMINAL, IMSL_ODE_NEG_NEQ);
            goto L_9000;
	}
	lv_state = (Env*)state;
	/*
	 * Cases - Initial entry, normal re-entry, interrupt re-entries
	 */
	if (lv_state->ido == 2)
		goto L_40;
	if (lv_state->ido == 3)
		goto L_9000;
	if (lv_state->ido == 4)
		goto L_70;
	if (lv_state->ido == 5 || lv_state->ido == 6)
		goto L_160;
	/* Case 1 - initial entry */
	;
	if (lv_state->wk == NULL) {
	    lv_state->wk = (Mfloat*)imsl_malloc(10*neq*sizeof(Mfloat));
	    /********** NEED WORKSPACE CHECK HERE  *******/
            if (lv_state->wk == NULL) {
                imsl_e1sti(1, neq);
                imsl_e1stl(1, "neq");
                imsl_ermes(IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
                goto L_9000;
                }
	}

	/* Initial WK(*,10) = YMAX(*) */
	for (i = 1; i <= neq; i++) {
		*WK(9, i - 1) = fabs(y[i - 1]);
	}
	/* Initialize EPS and HMAG */
	lv_state->est = F_ZERO;
	if (lv_state->hinit == F_ZERO) {
		lv_state->hmag = F_HALF * lv_state->hmax * imsl_ff_power(lv_state->tol, F_ONE / F_SIX);
	} else {
		lv_state->hmag = F_HALF * lv_state->hinit;
	}
	/*
	 * Set previous xend initially to initial value of X
	 */
	lv_state->xendpv = *x;
	lv_state->ixend = 0;
	lv_state->nsucst = 0;
	lv_state->nsucfl = 0;
	goto L_50;
	/*
	 * Case 2 - Normal re-entry (IDO.EQ.2) abort if xend reached, and
	 * either X changed or xend not changed
	 */
L_40:
	;
	if (lv_state->ixend != 0 && *x != lv_state->xendpv) {
		    /* Independent variable X = %(r1) is changed from */
		    /* the previous call. */
		imsl_e1str(1, *x);
		imsl_ermes(IMSL_TERMINAL, IMSL_ODE_T_CHANGED);
		goto L_9000;
	} else if (lv_state->ixend != 0 && xend == lv_state->xendpv) {
		    /* Final point xend = %(r1) is not changed from */
		    /* the previous call. */
		imsl_e1str(1, xend);
		imsl_ermes(IMSL_TERMINAL, IMSL_ODE_TEND_UNCHANGED);
		goto L_9000;
	}
	/* Re-initialize flag Ixend */
	lv_state->ixend = 0;
	/*
	 * Case 3 - Re-entry following an interrupt (IDO .EQ. 4 to 6)
	 */
L_50:
	;
	/*
	 * Loop through the following four stages, once for each trial step
	 * until the occurrence of one of the following (A) Normal return
	 * (with IDO .EQ. 2) on reaching xend in stage 4, or (B) An error
	 * return in stage 1 or 4, or (C) An interrupt return (with IDO .EQ.
	 * 4, 5 or 6) in stage 1 or 4.
	 */
L_60:
	;
	/*
	 * Stage 1 -- Prepare Do calculations of HMIN, HMAX, etc., and some
	 * parameter checking. End up with suitable values of HMAG, XTRIAL
	 * and HTRIAL for use in the integration step. Check MXFCN
	 */
	lv_state->newfcn = 7;
	if (lv_state->ido != 6)
		lv_state->newfcn++;
	if (lv_state->mxfcn > 0 && lv_state->nfcn + lv_state->newfcn > lv_state->mxfcn) {
		    /* Completion of the next step would make the */
		    /* number of function evaluations %(i1), but  */
		    /* only %(i2) evaluations are allowed.	  */
		imsl_e1sti(1, lv_state->nfcn + lv_state->newfcn);
		imsl_e1sti(2, lv_state->mxfcn);
		imsl_ermes(IMSL_FATAL, IMSL_ODE_TOO_MANY_EVALS);
		goto L_9000;
	}
	/* Check MXSTEP */
	if (lv_state->mxstep > 0 && lv_state->nstep >= lv_state->mxstep) {
		    /* Maximum number of steps allowed,		*/
		    /* %(i1), used.  The problem may be stiff.	*/
		imsl_e1sti(1, lv_state->mxstep);
		imsl_ermes(IMSL_FATAL, IMSL_ODE_TOO_MANY_STEPS);
		goto L_9000;
	}
	/* Calculate slope */
	if (lv_state->ido != 6) {
		imsl_e1usr("ON");
		(*fcn) (neq, *x, y, WK(0, 0));
		imsl_e1usr("OFF");
		lv_state->nfcn++;
		if (imsl_n1rty(0) > 3)
			goto L_9000;
	}
	/*
	 * Calculate HMIN - Use default unless value prescribed
	 */
	if (lv_state->hmin == F_ZERO) {
		/* Calculate default value of HMIN */
		lv_state->hmin = F_TEN * imsl_f_max(lv_tiny,
		    lv_eps * imsl_f_max(fabs(xend), fabs(*x)));
	}
	/* Calculate preliminary HMAG - */
	if (lv_state->nsucfl <= 1) {
		/*
		 * After a successful step, or at most one failure
		 */
		temp = F_TWO * lv_state->hmag;
		/*
		 * Avoid possible overflow
		 */
		if (lv_state->tol < imsl_fi_power(F_TWO / .9, 6) * lv_state->est)
			temp = 0.9 * imsl_ff_power(lv_state->tol / lv_state->est, F_ONE / F_SIX) * lv_state->hmag;
		/* Avoid reduction by more than half */
		lv_state->hmag = imsl_f_max(temp, F_HALF * lv_state->hmag);
	} else {
		/*
		 * After two or more successive failures
		 */
		lv_state->hmag *= F_HALF;
	}
	/* Check against HMAX and HMIN */
	lv_state->hmag = imsl_f_max(imsl_f_min(lv_state->hmag, lv_state->hmax), lv_state->hmin);
	/*
	 * Interrupt 1 (with IDO=4) if requested
	 */
	if (lv_state->intrp1 != 0) {
		lv_state->ido = 4;
		goto L_9000;
	}
	/* Resume here on re-entry */
L_70:
	;
	/*
	 * Calculate HMAG, XTRIAL - depending on preliminary HMAG, xend
	 */
	if (lv_state->hmag < fabs(xend - *x)) {
		/*
		 * Do not step more than half way to xend
		 */
		lv_state->hmag = imsl_f_min(lv_state->hmag, F_HALF * fabs(xend - *x));
		lv_state->xtrial = *x + sign(lv_state->hmag, xend - *x);
	} else {
		/* Hit xend exactly */
		lv_state->hmag = fabs(xend - *x);
		lv_state->xtrial = xend;
	}
	/* Calculate HTRIAL */
	lv_state->htrial = lv_state->xtrial - *x;
	/*
	 * Stage 2 -- Calculate YTRIAL WK(*,2), ..., WK(*,8) hold
	 * intermediate results needed in stage 3. WK(*,9) is temporary
	 * storage until finally it holds YTRIAL.
	 */
	scopy(neq, y, 1, WK(8, 0), 1);
	saxpy(neq, lv_state->htrial * lv_rk[0], WK(0, 0), 1, WK(8, 0), 1);
	imsl_e1usr("ON");
	(*fcn) (neq, *x + lv_state->htrial / F_SIX, WK(8, 0), WK(1, 0));
	imsl_e1usr("OFF");
	lv_state->nfcn++;
	if (imsl_n1rty(0) > 3)
		goto L_9000;

	for (k = 1; k <= neq; k++) {
		temp1 = imsl_sdot(2, WK(0, k - 1), neq, &lv_rk[1], 1);
		*WK(8, k - 1) = y[k - 1] + lv_state->htrial * temp1;
	}
	imsl_e1usr("ON");
	(*fcn) (neq, *x + F_FOUR * lv_state->htrial / 15.0, WK(8, 0), WK(2, 0));
	imsl_e1usr("OFF");
	lv_state->nfcn++;
	if (imsl_n1rty(0) > 3)
		goto L_9000;

	for (k = 1; k <= neq; k++) {
		temp1 = imsl_sdot(3, WK(0, k - 1), neq, &lv_rk[3], 1);
		*WK(8, k - 1) = y[k - 1] + lv_state->htrial * temp1;
	}
	imsl_e1usr("ON");
	(*fcn) (neq, *x + F_TWO * lv_state->htrial / F_THREE, WK(8, 0), WK(3, 0));
	imsl_e1usr("OFF");
	lv_state->nfcn++;
	if (imsl_n1rty(0) > 3)
		goto L_9000;

	for (k = 1; k <= neq; k++) {
		temp1 = imsl_sdot(4, WK(0, k - 1), neq, &lv_rk[6], 1);
		*WK(8, k - 1) = y[k - 1] + lv_state->htrial * temp1;
	}
	imsl_e1usr("ON");
	(*fcn) (neq, *x + F_FIVE * lv_state->htrial / F_SIX, WK(8, 0), WK(4, 0));
	imsl_e1usr("OFF");
	lv_state->nfcn++;
	if (imsl_n1rty(0) > 3)
		goto L_9000;

	for (k = 1; k <= neq; k++) {
		temp1 = imsl_sdot(5, WK(0, k - 1), neq, &lv_rk[10], 1);
		*WK(8, k - 1) = y[k - 1] + lv_state->htrial * temp1;
	}
	imsl_e1usr("ON");
	(*fcn) (neq, *x + lv_state->htrial, WK(8, 0), WK(5, 0));
	imsl_e1usr("OFF");
	lv_state->nfcn++;
	if (imsl_n1rty(0) > 3)
		goto L_9000;

	for (k = 1; k <= neq; k++) {
		temp1 = imsl_sdot(5, WK(0, k - 1), neq, &lv_rk[15], 1);
		*WK(8, k - 1) = y[k - 1] + lv_state->htrial * temp1;
	}
	imsl_e1usr("ON");
	(*fcn) (neq, *x + lv_state->htrial / 15.0, WK(8, 0), WK(6, 0));
	imsl_e1usr("OFF");
	lv_state->nfcn++;
	if (imsl_n1rty(0) > 3)
		goto L_9000;

	for (k = 1; k <= neq; k++) {
		temp1 = imsl_sdot(7, WK(0, k - 1), neq, &lv_rk[20], 1);
		*WK(8, k - 1) = y[k - 1] + lv_state->htrial * temp1;
	}
	imsl_e1usr("ON");
	(*fcn) (neq, *x + lv_state->htrial, WK(8, 0), WK(7, 0));
	imsl_e1usr("OFF");
	lv_state->nfcn++;
	if (imsl_n1rty(0) > 3)
		goto L_9000;
	/*
	 * Calculate YTRIAL, the extrapolated approximation and store in
	 * WK(*,9)
	 */
	for (k = 1; k <= neq; k++) {
		temp1 = imsl_sdot(8, WK(0, k - 1), neq, &lv_rk[27], 1);
		*WK(8, k - 1) = y[k - 1] + lv_state->htrial * temp1;
	}
	/*
	 * Stage 3 -- Calculate the error estimate EST Calculate the
	 * unweighted absolute error estimate vector
	 */
	for (k = 1; k <= neq; k++) {
		*WK(1, k - 1) = imsl_sdot(8, WK(0, k - 1), neq, &lv_rk[35], 1);
	}
	/*
	 * Calculate the weighted max norm of WK(*,2) as specified by the
	 * error control indicator INORM
	 */
	imsl_e1usr("ON");
	(*lv_state->vnorm) (neq, WK(1, 0), y, WK(9, 0), &temp);
	imsl_e1usr("OFF");
	if (imsl_n1rty(0) > 3)
		goto L_9000;
	/*
	 * Calculate EST - (The weighted max norm of WK(*,2))*HMAG*SCALE -
	 * EST is inxended to be a measure of the error per unit step in
	 * YTRIAL
	 */
	lv_state->est = temp * lv_state->hmag * lv_state->scale;
	/*
	 * Stage 4 -- Make decisions Set IDO=5 if step acceptable, else set
	 * IDO=6
	 */
	if (lv_state->est <= lv_state->tol) {
		lv_state->ido = 5;
	} else {
		lv_state->ido = 6;
	}
	/* Interrupt 2 if requested */
	if (lv_state->intrp2 != 0) {
		sswap(neq, y, 1, WK(8, 0), 1);
		goto L_9000;
	}
	/* Resume here on re-entry */
L_160:
	;
	if (lv_state->intrp2 != 0)
		sswap(neq, y, 1, WK(8, 0), 1);

	if (lv_state->ido == 5) {
		/*
		 * Step accepted, so update X, Y from XTRIAL, YTRIAL
		 */
		*x = lv_state->xtrial;
		scopy(neq, WK(8, 0), 1, y, 1);
		/* Update YMAX values */
		for (i = 1; i <= neq; i++) {
			*WK(9, i - 1) = imsl_f_max(*WK(9, i - 1), fabs(y[i - 1]));
		}
		lv_state->nsucst++;
		lv_state->nsucfl = 0;
		lv_state->nstep++;
		/*
		 * Return with IDO=2, xend saved, flag set
		 */
		if (*x == xend) {
			lv_state->ido = 2;
			lv_state->xendpv = xend;
			lv_state->ixend = 1;
			goto L_9000;
		}
	} else if (lv_state->ido == 6) {
		/*
		 * Step not accepted - Add 1 to number successive failures
		 */
		lv_state->nsucfl++;
		/* Error return if HMAG .LE. HMIN */
		if (lv_state->hmag <= lv_state->hmin) {
			    /* Unable to satisfy the error requirement. */
			    /* TOL = %(r1) may be too small.		*/
			imsl_e1str(1, lv_state->tol);
			imsl_ermes(IMSL_FATAL, IMSL_ODE_FAIL);
			goto L_9000;
		}
	}
	/* End stage 4 */
	goto L_60;
	/* End loop */
L_9000:
	;
	/* Set values for user */
	if (lv_state->p_nstep  != NULL) *lv_state->p_nstep  = lv_state->nstep;
	if (lv_state->p_htrial != NULL) *lv_state->p_htrial = lv_state->htrial;
	if (lv_state->p_nfcn   != NULL) *lv_state->p_nfcn   = lv_state->nfcn;

	E1POP("imsl_f_ode_runge_kutta","imsl_d_ode_runge_kutta");
	return;
}				/* end of function */
#undef WK


/* Structured by FOR_STRUCT, v0.2, on 01/08/90 at 17:44:57
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
static void l_i3prk(Mint neq, Mfloat *v, Mfloat *y, Mfloat *ymax, Mfloat *enorm)
#else
static void l_i3prk(neq, v, y, ymax, enorm)
	Mint             neq;
	Mfloat           v[], y[], ymax[], *enorm;
#endif
{
	Mint             k;
	Mfloat           weight;


	if (lv_state->inorm == 0) {
		/* max (absolute, relative) */
		*enorm = F_ZERO;
		for (k = 1; k <= neq; k++) {
			weight = imsl_f_max(F_ONE, fabs(y[k - 1]));
			*enorm = imsl_f_max(*enorm, fabs(v[k - 1]) / weight);
		}
		/* Absolute error control */
	} else if (lv_state->inorm == 1) {
		*enorm = fabs(v[imsl_isamax(neq, v, 1) - 1]);
		/* Relative error control */
	} else if (lv_state->inorm == 2) {
		*enorm = F_ZERO;
		for (k = 1; k <= neq; k++) {
			weight = imsl_f_max(fabs(y[k - 1]), lv_state->floor);
			*enorm = imsl_f_max(*enorm, fabs(v[k - 1]) / weight);
		}
	} else if (lv_state->inorm == 3) {
		/* Same as DGEAR's error control */
		*enorm = F_ZERO;
		for (k = 1; k <= neq; k++) {
			weight = imsl_f_max(F_ONE, fabs(ymax[k - 1]));
			*enorm += imsl_fi_power(v[k - 1] / weight, 2);
		}
		*enorm = sqrt(*enorm);
	}
	return;
}


    /*	set state to its initial values */
static void l_init_state()
{
	if (lv_tiny == 0.0) {
	    lv_tiny = imsl_amach(1);
	    lv_eps = imsl_amach(4);
        }
	lv_state->wk = NULL;
	lv_state->vnorm = l_i3prk;
	lv_state->tol = 100*lv_eps;
	lv_state->inorm = 0;
	lv_state->hmin = F_ZERO;
	lv_state->hinit = F_ZERO;
	lv_state->scale = F_ONE;
	lv_state->hmag = F_ZERO;
	lv_state->hmax = F_TWO;
	lv_state->mxstep = 500;
	lv_state->nstep = 0;
	lv_state->mxfcn = 0;
	lv_state->nfcn = 0;
	lv_state->intrp1 = 0;
	lv_state->intrp2 = 0;
	lv_state->inorm = 0;
	lv_state->floor = F_ONE;

	lv_state->p_nstep  = NULL;
	lv_state->p_htrial = NULL;
	lv_state->p_nfcn   = NULL;
}


#ifdef ANSI
void imsl_f_ode_runge_kutta_mgr(Mint ido, Mchar **state, ...)
#else
void imsl_f_ode_runge_kutta_mgr(ido, state, va_alist)
    Mchar	**state;
    Mint        ido;
    va_dcl
#endif
{
    va_list	argptr;

    VA_START(argptr, state);
    E1PSH("imsl_f_ode_runge_kutta_mgr", "imsl_d_ode_runge_kutta_mgr");
    IMSL_CALL(l_ode_runge_kutta_mgr(ido, state, argptr));
    va_end(argptr);
    E1POP("imsl_f_ode_runge_kutta_mgr", "imsl_d_ode_runge_kutta_mgr");
}


#ifdef ANSI
static VA_LIST_HACK l_ode_runge_kutta_mgr(Mint ido, Mchar **state, va_list argptr)
#else
static VA_LIST_HACK l_ode_runge_kutta_mgr(ido, state, argptr)
    Mchar	**state;
    Mint        ido;
    va_list	argptr;
#endif
{
    Mint	    code;
    Mint	    arg_number  = 2;
    Mint            user_floor = 0;

    if (ido == 1) {
        *state = (Mchar*)imsl_malloc(sizeof(Env));
        lv_state = (Env *)*state;
	l_init_state();
	lv_state->ido = 1;


#if 0
	if (lv_state->istat == 0) {
	    lv_state->istat = 1;
	} else {
		/* IDO = 1.  IDO can only be set to 1 in the */
		/* initial call to the routine, or if the    */
		/* previous call was made with IDO = 3.	     */
	    imsl_ermes(IMSL_TERMINAL, IMSL_ODE_IDO_1);
	}
#endif


    } else if (ido == 3) {
	lv_state = (Env *)*state;
	if (lv_state!=NULL && lv_state->wk!=NULL) imsl_free(lv_state->wk);
	imsl_free(lv_state);
        goto RETURN;


#if 0
	lv_state->ido = 3;
	if (lv_state->istat == 1) {
	    if (lv_state->ido == 3) {
		lv_state->istat = 0;
		goto RETURN;
	    }
	} else {
		/* IDO has been set to %(i1), but the routine */
		/* has not been initialized in a call with    */
		/* IDO = 1.				      */
	    imsl_e1sti(1, lv_state->ido);
	    imsl_ermes(IMSL_TERMINAL, IMSL_ODE_IDO_NOT_1);
	    goto RETURN;
	}
#endif
    }

    code = 1;
    while (code > 0) {
	code = va_arg(argptr, Mint);
	arg_number++;
	switch (code) {
	    case IMSL_TOL:
		lv_state->tol = (Mfloat) va_arg(argptr, Mdouble);
		arg_number++;
		if (lv_state->tol <= F_ZERO) {
			/* The tolerance TOL = %(r1).  It must be */
			/* greater than 0.0.			  */
		    imsl_e1str(1, lv_state->tol);
		    imsl_ermes(IMSL_TERMINAL, IMSL_ODE_NEG_TOL);
		}
		break;
	    case IMSL_TOL_ADR:
		lv_state->tol = *(va_arg(argptr, Mfloat *));
		arg_number++;
		if (lv_state->tol <= F_ZERO) {
			/* The tolerance TOL = %(r1).  It must be */
			/* greater than 0.0.			  */
		    imsl_e1str(1, lv_state->tol);
		    imsl_ermes(IMSL_TERMINAL, IMSL_ODE_NEG_TOL);
		}
		break;
	    case IMSL_NORM:
		arg_number++;
		lv_state->inorm = va_arg(argptr, Mint);
		break;
	    case IMSL_HINIT:
		arg_number++;
		lv_state->hinit = (Mfloat) va_arg(argptr, Mdouble);
		break;
	    case IMSL_HMIN:
		arg_number++;
		lv_state->hmin = (Mfloat) va_arg(argptr, Mdouble);
		break;
	    case IMSL_HMAX:
		arg_number++;
		lv_state->hmax = (Mfloat) va_arg(argptr, Mdouble);
		break;
	    case IMSL_SCALE:
		arg_number++;
		lv_state->scale = (Mfloat) va_arg(argptr, Mdouble);
		lv_state->hmax  = imsl_f_min(lv_state->hmax, F_TWO / lv_state->scale);
		break;
	    case IMSL_FLOOR:
		arg_number++;
		lv_state->floor = (Mfloat) va_arg(argptr, Mdouble);
                user_floor = 1;
		break;


	    case IMSL_HINIT_ADR:
		arg_number++;
		lv_state->hinit = *(va_arg(argptr, Mfloat *));
		break;
	    case IMSL_HMIN_ADR:
		arg_number++;
		lv_state->hmin = *(va_arg(argptr, Mfloat *));
		break;
	    case IMSL_HMAX_ADR:
		arg_number++;
		lv_state->hmax = *(va_arg(argptr, Mfloat *));
		break;
	    case IMSL_SCALE_ADR:
		arg_number++;
		lv_state->scale = *(va_arg(argptr, Mfloat *));
		lv_state->hmax  = imsl_f_min(lv_state->hmax, F_TWO / lv_state->scale);
		break;
	    case IMSL_FLOOR_ADR:
		arg_number++;
		lv_state->floor = *(va_arg(argptr, Mfloat *));
                user_floor = 1;
		break;
	    case IMSL_MAX_NUMBER_STEPS:
		arg_number++;
		lv_state->mxstep = va_arg(argptr, Mint);
		break;
	    case IMSL_MAX_NUMBER_FCN_EVALS:
		arg_number++;
		lv_state->mxfcn = va_arg(argptr, Mint);
		break;
	    case IMSL_VNORM:
		arg_number++;
		lv_state->vnorm = va_arg(argptr, Vnorm);
		break;
	    case IMSL_NSTEP:
		arg_number++;
		lv_state->p_nstep = va_arg(argptr, Mint*);
		break;
	    case IMSL_NFCN:
		arg_number++;
		lv_state->p_nfcn = va_arg(argptr, Mint*);
		break;
	    case IMSL_HTRIAL:
		arg_number++;
		lv_state->p_htrial = va_arg(argptr, Mfloat*);
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

	/* Error return if HMIN > HMAX */
    if (lv_state->hmin > lv_state->hmax) {
	    /* HMIN = %(r1) is greater than HMAX = %(r2). */
	imsl_e1str(1, lv_state->hmin);
	imsl_e1str(2, lv_state->hmax);
	imsl_ermes(IMSL_TERMINAL, IMSL_HMIN_GT_HMAX);
    }
    if (lv_state->hmin < 0) {
	imsl_e1str(1, lv_state->hmin);
/*	(5, 1, "hmin = %(r1). The minimum stepsize, hmin, must be greater than or equal to zero.");
*/
        imsl_ermes(IMSL_TERMINAL, IMSL_MIN_STEPSIZE_TOO_SMALL);
    }
    if (lv_state->hmax < 0) {
	imsl_e1str(1, lv_state->hmax);
/*	(5, 2, "hmax = %(r1). The maximum stepsize, hmax, must be greater than or equal to zero.");
*/
        imsl_ermes(IMSL_TERMINAL, IMSL_MAX_STEPSIZE_TOO_SMALL);
    }
    if (user_floor && (lv_state->inorm != 2)) {
/*      (5,3,"The optional argument IMSL_FLOOR may only be specified  when norm = 2.");
*/
        imsl_ermes(IMSL_TERMINAL, IMSL_INVALID_IMSL_FLOOR_USAGE);
    }
    if ((lv_state->inorm <0) || (lv_state->inorm >3)) {
        imsl_e1sti(1,lv_state->inorm); 
/*      (5,4,"inorm = %(i1). Valid values for inorm include 0, 1, 2, and 3.");
*/
        imsl_ermes(IMSL_TERMINAL, IMSL_INORM_OUTSIDE_OF_RANGE);
    }
RETURN:
    return (argptr);
}
