#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static VA_LIST_HACK  PROTO(l_lin_prog,(Mint, Mint, Mfloat*, Mfloat*, 
                      Mfloat*, va_list argptr));
static void	PROTO(l_d3prs,(Mint, Mint, Mfloat*, Mint, Mfloat*,
                      Mfloat*, Mfloat*, Mfloat*, Mfloat*, Mfloat*,
                      Mfloat*, Mfloat*, Mfloat*, Mfloat*, Mint*, Mint,
                      Mint*, Mint*, Mint*, Mint*, Mfloat*, Mfloat*,
                      Mint*, Mint));
static void	PROTO(l_d4prs,(Mint, Mint, Mfloat*, Mint, Mfloat*,
                      Mfloat*, Mint, Mfloat*, Mfloat*, Mint*, Mint*,
                      Mint*));
static void	PROTO(l_d5prs,(Mint, Mint, Mfloat*, Mint, Mfloat*,
                      Mfloat*, Mfloat*, Mfloat*, Mint*, Mint*, 
                      Mfloat*));
static void	PROTO(l_d6prs,(Mint, Mfloat*, Mint, Mint*, Mint*,
                      Mfloat*, Mint, Mfloat*, Mint, Mfloat*, Mfloat*));
static void	PROTO(l_d7prs,(Mint, Mint, Mfloat*, Mfloat*, Mfloat*,
                      Mfloat*, Mint*, Mfloat*, Mfloat*, Mint*, 
                      Mfloat, Mfloat*, Mint*));
static void	PROTO(l_d8prs,(Mint, Mint, Mfloat*, Mint,  Mfloat*,
                      Mfloat*, Mint, Mint*, Mfloat*));
static void	PROTO(l_d9prs,(Mint, Mint, Mint*, Mint*, Mint*, Mint*,
                      Mfloat*, Mfloat*, Mfloat*, Mfloat*, Mfloat, 
                      Mfloat, Mfloat, Mfloat*, Mint*));
static void	PROTO(l_d10rs,(Mint, Mint, Mint, Mint*, Mint*, Mint*,
                      Mfloat*, Mfloat, Mfloat, Mfloat*, Mfloat*,
                      Mfloat*, Mfloat*, Mfloat*, Mfloat*, Mint*));
static void	PROTO(l_d11rs,(Mint, Mint, Mint, Mint, Mint*, Mint*,
                      Mint*, Mfloat , Mfloat*, Mfloat*, Mfloat*,
                      Mfloat*, Mfloat*, Mfloat*, Mfloat,Mfloat*, 
                      Mint*));
static void	PROTO(l_d12rs,(Mint, Mint, Mfloat*, Mint, Mint*, 
                      Mfloat, Mfloat, Mfloat*, Mfloat*, Mfloat*,
                      Mint*));
static void	PROTO(l_d13rs,(Mint, Mfloat*, Mint, Mint*, Mfloat*,
                      Mint, Mint*, Mfloat*, Mint*, Mfloat*));
static void	PROTO(l_d14rs,(Mint, Mfloat*, Mint, Mint*, Mfloat*,
                      Mint, Mfloat*));

static Mfloat	*lv_x;

#ifdef ANSI
Mfloat *imsl_f_lin_prog(Mint m, Mint n, Mfloat *a, Mfloat b[],
                        Mfloat c[], ...)
#else
Mfloat *imsl_f_lin_prog(m, n, a, b, c, va_alist)
    Mint	m;
    Mint	n;
    Mfloat	*a;
    Mfloat	*b;
    Mfloat	*c;
    va_dcl
#endif
{
    va_list	argptr;

    VA_START(argptr,c);
    E1PSH("imsl_f_lin_prog", "imsl_d_lin_prog");
    lv_x = NULL;
    IMSL_CALL(l_lin_prog(m, n, a, b, c, argptr));
    va_end(argptr);
    E1POP("imsl_f_lin_prog", "imsl_d_lin_prog"); 
    return lv_x;
}


#ifdef ANSI
static VA_LIST_HACK l_lin_prog(Mint m, Mint n, Mfloat *a, Mfloat *b,
                          Mfloat *c, va_list argptr)
#else
static VA_LIST_HACK l_lin_prog(m, n, a, b, c, argptr)
    Mint	m;
    Mint	n;
    Mfloat	*a;
    Mfloat	*b;
    Mfloat	*c;
    va_list	argptr;
#endif
{
    Mint            code;
    Mint            arg_number  = 5;
    Mint            a_col_dim   = n;
    Mint            max_itn     = 10000;
    Mfloat          *bu         = NULL;
    Mint            user_bu     = 0;
    Mfloat          *xlb        = NULL;
    Mint            user_xlb    = 0;
    Mfloat          *xub        = NULL;
    Mint            user_xub    = 0;
    Mfloat          *d          = NULL;
    Mint            user_d      = 0;
    Mfloat          *obj        = NULL;
    Mint            user_obj    = 0;
    Mint            *irtype     = NULL;
    Mint            user_basis  = 0;
    Mint            inout       = 0;
    Mint            *ibasis     = NULL;
    Mint            *ibb        = NULL;
    Mint            user_pd     = 0;
    Mint            maximiz     = 0;
    Mfloat          **pd        = NULL;
    Mint            user_irtype = 0;
    Mfloat          *work	= NULL;
    Mfloat          *at  	= NULL;
    Mint            *iwork      = NULL;
    Mint            i, k;
    Mint            user_x      = 0;
    Mfloat          f;

    code = 1;
    while (code > 0) {
        code = va_arg(argptr, Mint);
        arg_number++;
        switch (code) {
            case IMSL_RETURN_USER:
                lv_x = va_arg(argptr, Mfloat*);
                arg_number++;
                user_x = 1;
                break;
            case IMSL_DUAL_USER:
                user_d = 1;
                d = va_arg(argptr, Mfloat*);
                arg_number++;
                break;
            case IMSL_DUAL:
                user_pd = 1;
                pd = va_arg(argptr, Mfloat**);
                arg_number++;
                break;
            case IMSL_A_COL_DIM:
                a_col_dim = va_arg(argptr, Mint);
                arg_number++;
                break;
            case IMSL_UPPER_LIMIT:
                user_bu = 1;
                bu = va_arg(argptr, Mfloat*);
                arg_number++;
                break;
	    case IMSL_CONSTR_TYPE:
		user_irtype = 1;
		irtype = va_arg(argptr, Mint*);
		arg_number ++;
		break;
	    case IMSL_LOWER_BOUND:
		user_xlb = 1;
		xlb = va_arg(argptr, Mfloat*);
		arg_number++;
		break;
	    case IMSL_UPPER_BOUND:
		user_xub = 1;
		xub = va_arg(argptr, Mfloat*);
		arg_number++;
		break;
	    case IMSL_BASIS:
		user_basis = 1;
                inout  = va_arg(argptr, Mint);
		arg_number ++;
		ibasis = va_arg(argptr, Mint*);
		arg_number ++;
		ibb = va_arg(argptr, Mint*);
		arg_number ++;
		break;
	    case IMSL_OBJ:
		user_obj = 1;
		obj = va_arg(argptr, Mfloat*);
		arg_number++;
		break;
            case IMSL_MAX_ITN:
                max_itn = va_arg(argptr, Mint);
                arg_number++;
                break;
            case IMSL_MAXIMIZATION:
                maximiz = 1;
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
        imsl_ermes(IMSL_TERMINAL, IMSL_NEGATIVE_ORDER);

    } else if (n > a_col_dim) {
        imsl_e1sti(1, m);
        imsl_e1sti(2, a_col_dim);
        imsl_e1stl(1, "a");
        imsl_ermes(IMSL_TERMINAL, IMSL_COL_DIM_LESS_ORDER);
    }

    if (m <= 0) {
	imsl_e1sti(1, m);
	imsl_ermes(IMSL_TERMINAL, IMSL_NEG_CONSTRAINT_VALUE);
    }

    if (imsl_n1rty(0)) goto RETURN;

    if (!user_irtype) {
        irtype  = (Mint *) imsl_malloc(m*sizeof(*irtype));
        iset(m, 0, irtype, 1);
    }

    if (!user_bu)  bu = b; 

    if (!user_xlb) {
	xlb  = (Mfloat *) imsl_malloc(n*sizeof(*xlb));
        sset(n, F_ZERO, xlb, 1);
    }

    if (!user_xub) {
	xub  = (Mfloat *) imsl_malloc(n*sizeof(*xub));
        sset(n, -1.0e30, xub, 1);
    }

    if (!user_d)  d = (Mfloat *) imsl_malloc(m*sizeof(*d));
    if (!user_basis) {
        ibasis = (Mint *) imsl_malloc((m+n)*sizeof(*ibasis));
        ibb = (Mint *) imsl_malloc((m+n)*sizeof(*ibb));
    }
    at    = (Mfloat *) imsl_malloc(m*n*sizeof(*at));
    work  = (Mfloat *) imsl_malloc(m*(m+28)*sizeof(*work));
    iwork = (Mint *) imsl_malloc((m*27+n)*sizeof(*iwork));

    if (bu==NULL || irtype==NULL || work==NULL || xlb==NULL || xub==NULL
        || iwork==NULL || at==NULL || ibasis==NULL || ibb==NULL) {
	imsl_e1stl(1, "n");
	imsl_e1sti(1, n);
	imsl_e1stl(2, "m");
	imsl_e1sti(2, m);
	imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
	goto FREE_SPACE;
    }

    if (lv_x == NULL) {
	lv_x = (Mfloat *)imsl_malloc (n*sizeof(*lv_x));
	if (lv_x == NULL) {
	    imsl_e1stl(1, "n");
	    imsl_e1sti(1, n);
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    user_bu = 0;
	    goto FREE_SPACE;
	}
    }
    if (maximiz) {
        sscal(n, -1.0e0, c, 1);
    }

    /* Transpose A */
    for (i = 0; i < n; i++) {
         k = i * m;
         scopy (m, &a[i], a_col_dim, &at[k], 1);
    }

    imsl_d2prs(m, n, at, m, b, bu, c, irtype, xlb, xub, inout, ibasis, 
            ibb, &f, lv_x, d, work, iwork, max_itn);

    if (maximiz) {
        sscal(n, -1.0e0, c, 1);
        sscal(m, -1.0e0, d, 1);
        f = -f;
    }

    if (user_pd)   *pd = d;
    if (user_obj)  *obj = f;

FREE_SPACE:
    if (work != NULL)                      imsl_free(work);
    if (iwork != NULL)                     imsl_free(iwork);
    if (at != NULL)                        imsl_free (at);
    if (!user_xlb && xlb != NULL)          imsl_free(xlb); 
    if (!user_xub && xub != NULL)          imsl_free(xub);
    if (!user_irtype && irtype != NULL)    imsl_free(irtype);
    if (!user_basis && ibasis != NULL) {
        imsl_free(ibasis);
        imsl_free(ibb);
    }
    if (!user_d && !user_pd && d != NULL)  imsl_free(d);

RETURN:
    if (imsl_n1rty(0) > 4) {
        if (!user_x && lv_x != NULL) imsl_free(lv_x);
        lv_x = NULL;
    }
    return (argptr);
}

/* -----------------------------------------------------------------------
    IMSL Name:  D2PRS/DD2PRS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    April 21, 1989

    Purpose:    Solve a linear programming problem via the revised
                simplex algorithm.

    Usage:      CALL D2PRS (M, NVAR, A, LDA, BL, BU, C, IRTYPE, XLB, XUB, 
                            OBJ, XSOL, DSOL, AWK, LDAWK, WK, IWK, MAXITN)

    Arguments:  See DLPRS.

    Remark:     See DLPRS.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1989 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------  */
#ifdef ANSI
void imsl_d2prs(Mint m, Mint nvar, Mfloat *a, Mint lda, Mfloat *bl, 
                    Mfloat *bu, Mfloat *c, Mint *irtype, Mfloat *xlb,
                    Mfloat *xub, Mint inout, Mint *ibasis, Mint *ibb, 
                    Mfloat *obj, Mfloat *xsol, Mfloat *dsol, Mfloat *wk,
                    Mint *iwk, Mint maxitn)
#else
void imsl_d2prs(m, nvar, a, lda, bl, bu, c, irtype, xlb, xub, inout, 
                    ibasis, ibb, obj, xsol, dsol, wk, iwk, maxitn)
	Mint            m, nvar;
	Mfloat          *a;
	Mint            lda;
	Mfloat          bl[], bu[], c[];
	Mint            irtype[], ibasis[], ibb[], inout;
	Mfloat          xlb[], xub[], *obj, xsol[], dsol[], wk[];
	Mint            iwk[], maxitn;
#endif
{
#define A(I_,J_)	(a+(I_)*(lda)+(J_))
	Mint             i, info, itmp, j;
	Mfloat           neg_flag, pos_flag, sqteps;

	imsl_e1psh("D2PRS ");
	sqteps = sqrt(imsl_amach(4));

	if (imsl_n1rcd(0) == 0) {

	    for (j = 0; j < m; j++) {
		 if (bu[j] < bl[j] && irtype[j] == 3) {
		     imsl_e1sti(1, j);
		     imsl_e1str(2, bl[j]);
		     imsl_e1str(3, bu[j]);
                     imsl_e1stl(1, "b");
                     imsl_e1stl(2, "bl");
                     imsl_e1stl(3, "bu");
/*		     imsl_ermes(4, 5, "The bounds given for B(%(i1)) are     */
/*                                     inconsistent:  BL(%(i1)) = %(r2) must */
/*                                     not exceed BU(%(i1)) = %(r3).");      */
                     imsl_ermes(IMSL_FATAL, IMSL_BOUNDS_INCONSISTENT);
		     goto L_9000;
		 }
		 itmp = irtype[j];
		 if (itmp == 0) {
		     iwk[nvar + j] = 3;
		 } else if (itmp == 1) {
		     iwk[nvar + j] = 2;
		 } else if (itmp == 2) {
		     iwk[nvar + j] = 1;
		 } else if (itmp == 3) {
		     iwk[nvar + j] = 3;
		 } else {
		     imsl_e1sti(1, j);
		     imsl_e1sti(2, itmp);
/*		     imsl_ermes(5, 5, "The values for the type of constraint */
/*                                     must be either 0, 1 or 2, yet         */
/*                                     IRTYPE(%(i1)) = %(i2) is given.");    */
                     imsl_ermes(IMSL_TERMINAL, IMSL_WRONG_CONSTRAINT_TYPE);
		     goto L_9000;
		 }
	    }

            neg_flag = -1.0e30 + sqteps;
            pos_flag = 1.0e30 - sqteps;
	    for (i = 0; i < nvar; i++) {
		 if (xlb[i] >= pos_flag && xub[i] <= neg_flag) {
		     iwk[i] = 4;
		 } else {
		     if (xlb[i] >= pos_flag) {
			 iwk[i] = 2;
		     } else {
			 if (xub[i] <= neg_flag) {
			     iwk[i] = 1;
			 } else {
			     iwk[i] = 3;
			     if (xub[i] < xlb[i]) {
				 imsl_e1sti(1, i);
				 imsl_e1str(2, xlb[i]);
				 imsl_e1str(3, xub[i]);
                                 imsl_e1stl(1, "x");
                                 imsl_e1stl(2, "xlb");
                                 imsl_e1stl(3, "xub");
                              /* imsl_ermes(4, 5,"The bounds given for X(i1) */
                              /*                   are inconsistent: XLB(i1) */
                              /*                   = %(r2) must not exceed   */
                              /*                   XUB(i1) = (r3).");        */
                               imsl_ermes(IMSL_FATAL, IMSL_BOUNDS_INCONSISTENT);
				 goto L_9000;
		             }
			 }
		     }
		 }
	    }

	    /* Call D3PRS to solve LP problem. */

	    l_d3prs(m, nvar, a, lda, bl, bu, c, xlb, xub, xsol, dsol, wk,
		    &wk[m * m], &wk[m * (m + 1)], iwk, inout, ibasis,
                    ibb, &iwk[(nvar + m)], &iwk[(nvar + m) + m], 
                    &wk[m * (m + 2)], &wk[m * (m + 3)], &info, maxitn);
	    if (info == 0) {
		*obj = imsl_sdot(nvar, xsol, 1, c, 1);
	    } else {
		if (info == 2)
		    *obj = imsl_sdot(nvar, xsol, 1, c, 1);
	    }
	}
L_9000:
	imsl_e1pop("D2PRS ");
	return;
}				/* end of function */

/* -----------------------------------------------------------------------
    IMSL Name:  D3PRS/DD3PRS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 5, 1989

    Purpose:    Minimize TRANS(C)*X subject to A*X = B.

    Usage:      CALL D3PRS (M, N, A, LDA, BL, BU, C, XLB, XUB, D, PI, B,
                            XB, CB, ISTAT, IBASIS, IBB, IPVT, IETA, Y,
                            ETA, INFO, MAXITN)

    Arguments:
       M      - Number of constraints.  (Input)
       N      - Number of variables.  (Input)
       A      - Matrix of dimension M by N containing the coefficients
                of the M constraints.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
                LDA must be at least M.
       BL     - Vector of length M containing the lower limit of the
                general constraints; if there is no lower limit on the
                I-th constraint, then BL(I) would not be referred.
                (Input)
       BU     - Vector of length M containing the upper limit of the
                general constraints; if there is no upper limit on
                the I-th constraint, then BU(I) would not be referred;
                if there is no range constraint, BL and BU can share
                the same storage locations.  (Input)
       C      - Vector of length N containing the coefficients of
                the objective function.  (Input)
       XLB    - Vector of length N containing the lower bound on the
                variables.  (Input)
                If there is no lower bound on a variable, then 1.0E30
                should be set as the lower bound.
       XUB    - Vector of length N containing the upper bound on the
                variables.  (Input)
       D      - Vector of length N containing reduced costs.
       PI     - Vector of length M containing the simplex multipliers.
       B      - Matrix of dimension M by M which contains the initial
                basis.
       XB     - Vector of length M containing values of the basic
                variables.
       CB     - Vector of length M containing objective coefficients.
       ISTAT  - Vector of length M+N indicating the bound types of the
                variables.
       IBASIS - Vector of length M+N; elements 1-M point to basic
                columns
                others to nonbasic columns.
       IBB    - Vector of length M+N which signifies at which bound,
                if any, a nonbasic variable is set, or that the variable
                is basic.
       IPVT   - Vector of length N containing the pivoting information
                for the LU factorization.
       IETA   - Vector of length 25*m containing pivoting information
                regarding the eta vectors.
       Y      - Vector of length m used for intermediate solutions
                to linear systems.
       ETA    - Vector of length 25*m representing the sequence of
                elementary matrices used in the product form of the
                inverse of B.
       INFO   - A flag which contains the informational error code
                returned from the revised simplex algorithm.  (Output)
       MAXITN - The maximum number of iterations.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1988 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_d3prs(Mint m, Mint n, Mfloat *a, Mint lda, Mfloat *bl,
                    Mfloat *bu, Mfloat *c, Mfloat *xlb, Mfloat *xub,
                    Mfloat *d, Mfloat *pi, Mfloat *b, Mfloat *xb,
                    Mfloat *cb, Mint *istat, Mint inout, Mint *ibasis, 
                    Mint *ibb, Mint *ipvt, Mint *ieta, Mfloat *y, 
                    Mfloat *eta, Mint *info, Mint maxitn)
#else
static void l_d3prs(m, n, a, lda, bl, bu, c, xlb, xub, d, pi, b, xb,
	            cb, istat, inout, ibasis, ibb, ipvt, ieta, y, eta, 
                    info, maxitn)
	Mint            m, n;
	Mfloat          *a;
	Mint            lda;
	Mfloat          bl[], bu[], c[], xlb[], xub[], d[], pi[], *b, xb[],
	                cb[];
	Mint            istat[], ibasis[], ibb[], ipvt[], ieta[];
	Mfloat           y[], eta[];
	Mint            *info, inout, maxitn;
#endif
{
#define A(I_,J_)	(a+(I_)*(lda)+(J_))
#define B(I_,J_)	(b+(I_)*(m)+(J_))
	Mint            artif, bigerr, change, feas, finite, found;
	Mint            i, ichck, icount, idcomp, ienter, ileave, k,
	                imax, iphase, iprint, iredfq, itmp, itn, j, jq;
	Mfloat          dirnrm, eps, obj, pinrm, theta, tinf,
	                toldj, toln, tolx;
        FILE            *nout;

	imsl_e1psh("l_d3prs ");
	eps = imsl_amach(4);
 	tolx = sqrt(eps);
/*	tolx = 100.0*eps;*/
	toldj = tolx;
	toln = tolx;
	artif = 1;
	bigerr = 0;
	iredfq = 40;
	ichck = 20;
	imax = imsl_i_min(25 * m, 1000);
	*info = 0;
	change = 0;
	iphase = 1;
	itn = 0;
	iprint = 0;
	idcomp = 0;
        imsl_umach(2,&nout);

	/* Set up initial B, ibasis, and istat.  */

        if (!inout) {
	   l_d4prs(m, n, a, lda, xlb, xub, artif, b, xb, ibasis, ibb, istat);
          goto L_30;

	} else {
	    artif = 0;
	    for (i = 1; i <= n; i++) {
		k = m + i;
		j = ibasis[k - 1];
		if (istat[j - 1] == 1) {
		    ibb[j - 1] = 1;
		} else if (istat[j - 1] == 2) {
		    ibb[j - 1] = 2;
		}
	    }
	}

	/* Reset B for reinversion.  */
L_10:
	;
	for (i = 1; i <= m; i++) {
	     j = ibasis[i - 1];
	     if (j <= n) {
		 scopy(m, A(j - 1, 0), 1, B(i - 1, 0), 1);
	     } else {
		 sset(m, F_ZERO, B(i - 1, 0), 1);
		 *B(i - 1, j - n - 1) = -F_ONE;
	     }
	}
	/* Decomposition, b = lu. */
L_30:
	imsl_l2trg(m, b, m, b, m, ipvt, y);
	idcomp += 1;
	icount = 1;
	/*
	 * Compute basic variable values. First, compute rhs = -N*xn,
	 * then, solve B xb = rhs.
	 */
	if (artif && itn == 0)
	    goto L_40;
	if (!bigerr) {
	    l_d5prs(m, n, a, lda, xlb, xub, bl, bu, istat, ibb, cb);
	} else {
	    bigerr = 0;
	}
	l_d6prs(m, b, m, ipvt, ieta, eta, icount, cb, 1, xb, y);

	/* Main loop. */
L_40:
	itn += 1;
	/* Check feasibility. */
	l_d7prs(m, n, xlb, xub, bl, bu, istat, cb, xb, ibasis, tolx, &tinf,
		   &feas);

	if (feas && iphase == 1) {
	    iphase = 2;
	} else if (!feas && iphase == 2) {
	    iphase = 1;
	    *info = 4;
	}
	/*
	 * Determine cb, if phase II. Otherwise, cb = -1 :  xb(i) > bu(i) 1 :
	 * xb(i) < bl(i) 0 :  bl(i) <= xb(i) <= bu(i).
	 */
	if (iphase == 2) {
	    for (i = 1; i <= m; i++) {
		 j = ibasis[i - 1];
		 if (j > n) {
		     cb[i - 1] = F_ZERO;
		 } else {
		     cb[i - 1] = c[j - 1];
		 }
	    }
	}

	/* 1.  Compute the simplex multiplier, pi by solving T B pi = cb.  */

	l_d6prs(m, b, m, ipvt, ieta, eta, icount, cb, 2, pi, y);

	pinrm = imsl_sasum(m, pi, 1) / m;

	/* 2.  Price out the nonbasic columns. */

	l_d8prs(m, n, a, lda, c, pi, iphase, ibasis, d);
	/*
	 * a.  If none selected, leave (we've found the optimal feasible
	 * solution-- PHASE II or this is an infeasible problem--PHASE I). b.
	 * Otherwise, say selected variable corresponds to j=q. 3.  Find
	 * variable to enter basis.
	 */
	l_d9prs(m, n, &ienter, ibasis, istat, ibb, xlb, xub, bl, bu, eps,
		toldj, pinrm, d, &found);
	if (!found)
	    goto L_70;

	/* 4.  Update the entering column. Solve B alphaq = a(q).  */

	jq = ibasis[ienter - 1];
	if (jq <= n) {
	    l_d6prs(m, b, m, ipvt, ieta, eta, icount, A(jq - 1, 0), 1, cb, y);
	} else {
	    sset(m, F_ZERO, cb, 1);
	    cb[jq - n - 1] = -F_ONE;
	    l_d6prs(m, b, m, ipvt, ieta, eta, icount, cb, 1, cb, y);
	}

	if (ibb[jq - 1] == 2 || (istat[jq - 1] == 4 && d[ienter - m - 1] >
                                                                       F_ZERO))
	    sscal(m, -F_ONE, cb, 1);
	dirnrm = imsl_sasum(m, cb, 1) / m; 
   /*     dirnrm = 1.0;*/
	/*
	 * 5.  Find the leaving basic variable. a.  If none selected, leave
	 * (we have an unbounded solution for the entire problem). b.
	 * Otherwise, say selected variable corresponds to i=p.
	 */
	l_d10rs(m, n, ienter, &ileave, ibasis, istat, &theta, dirnrm,
		   toln, cb, xlb, xub, bl, bu, xb, &finite);
	if (!finite)
		goto L_130;

	/* 6.  Update xb, ibasis, ibb. */

	l_d11rs(m, n, ienter, ileave, ibasis, istat, ibb, theta, xb,
		   cb, xlb, xub, bl, bu, tolx, d, &change);

	/* Return to Step 0. */

	if (iprint >= 2) {
	    fprintf(nout, "    itn.    obj.    theta    ienter   ileave  \n");
	    fprintf(nout, "   %5ld %8.2e %8.2e   %4ld     %4ld\n", itn, tinf,
                           theta, ienter, ileave);
	}
	if (iprint >= 3) {
	    fprintf(nout, " ibasis = \n");
	    fprintf(nout, "  ");
	    for (i = 1; i <= (m + n); i++) {
		 fprintf(nout, "%4ld", ibasis[i - 1]);
	    }
	    fprintf(nout, "  \n");
	    fprintf(nout, " ibb = \n");
	    fprintf(nout, "  ");
	    for (i = 1; i <= (m + n); i++) {
		 fprintf(nout, "%4ld", ibb[i - 1]);
	    }
	    fprintf(nout, "  \n");
	}
	if (itn > maxitn)
	    goto L_100;

	/* Reinversion every iredfq iterations.  */

	if (mod(itn, iredfq) == 0) {
	    goto L_10;
	}
	/* Check for errors every ichck iterations.  */

	if (mod(itn, ichck) == 0) {
	    l_d5prs(m, n, a, lda, xlb, xub, bl, bu, istat, ibb, cb);
	    l_d12rs(m, n, a, lda, ibasis, tolx, eps, pi, xb, cb, &bigerr);
	    if (bigerr) goto L_10;
	}

	/* Update LU.  */

	if (change) {
	    itmp = icount + m - ileave;
	    if (itmp <= imax) {
		l_d13rs(m, b, m, ipvt, y, ileave, ieta, eta, &icount, pi);
		goto L_40;
	    } else {
		goto L_10;
	    }
	}
	goto L_40;

L_70:
	if (iphase == 2) {
	    for (i = 1; i <= m; i++) {
		 j = ibasis[i - 1];
		 if (j <= n) d[j - 1] = xb[i - 1];
	    }
	    for (i = 1; i <= n; i++) {
		 j = ibasis[m + i - 1];
		 if (j <= n) {
		     if (ibb[j - 1] == 1) {
			 d[j - 1] = xlb[j - 1];
		     } else {
			 d[j - 1] = xub[j - 1];
		     }
		 }
	    }
	} else {
	    *info = 3;
	}
	goto L_9000;
	/* Maximum iteration exceeded. */
L_100:
	*info = 2;
	if (iphase == 2) {
	    for (i = 1; i <= m; i++) {
		 j = ibasis[i - 1];
		 if (j <= n) d[j - 1] = xb[i - 1];
	    }
	    for (i = 1; i <= n; i++) {
		 j = ibasis[m + i - 1];
		 if (j <= n) {
		     if (ibb[j - 1] == 1) {
			 d[j - 1] = xlb[j - 1];
		     } else {
			 d[j - 1] = xub[j - 1];
		     }
		 }
	    }
	} else {
	    *info = -2;
	}
	goto L_9000;

	/* Unbounded solution or numerical error.  */

L_130:
	if (iphase == 2) {
	    *info = 1;
	} else {
	    *info = 4;
	}

L_9000:
	if (iprint >= 1) {
	    obj = imsl_sdot(n, d, 1, c, 1);
	    fprintf(nout, " No. of itn. = %d \n", itn - 1);
	    fprintf(nout, " No. of decomp. = %d \n", idcomp);
	    fprintf(nout, " obj = %g \n", obj);
	}
	/* Error messages */
	if (*info == 1) {
            imsl_ermes(IMSL_WARNING, IMSL_PROB_UNBOUNDED);

	} else if (*info == 2) {
            imsl_ermes(IMSL_WARNING, IMSL_TOO_MANY_ITN);

	} else if (*info == 3) {
            imsl_ermes(IMSL_WARNING, IMSL_PROB_INFEASIBLE);

	} else if (*info == 4) {
/*	    imsl_ermes(4, 4, "Numerical difficulty occurred; using double    */
/*                            precision may help.");                         */
            imsl_ermes(IMSL_FATAL, IMSL_NUMERIC_DIFFICULTY);

	} else if (*info == -2) {
            imsl_ermes(IMSL_WARNING, IMSL_TOO_MANY_ITN);
	}
	imsl_e1pop("l_d3prs ");
	return;
}				/* end of function */

/* -----------------------------------------------------------------------
    IMSL Name:  D4PRS/DD4PRS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 5, 1989

    Purpose:    Set up a linear programming problem to be solved via the
                revised simplex algorith.

    Usage:      CALL D4PRS (M, NVAR, A, LDA, XLB, XUB, ARTIF, B,
                            XB, IBASIS, IBB, ISTAT)

    Arguments:
       M      - Number of constraints.  (Input)
       NVAR   - Number of variables.  (Input)
       A      - Matrix of dimension M by NVAR containing the
                coefficients of the M constraints.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
                LDA must be at least M.
       XLB    - Vector of length NVAR containing the lower bounds on
                the variables.  (Input)
                If there is no lower bound on a variable, then 1.0E30
                should be set as the lower bound.
       XUB    - Vector of length NVAR containing the upper bounds on
                variables.  (Input)
                If there is no upper bound on a variable, then -1.0E30
                should be set as the upper bound.
       ARTIF  - A flag indicating whether an artificial start to the
                problem is to be employed.  (Input)
       B      - Matrix of dimension M by M which contains the initial
                basis.  (Output)
       XB     - Vector of length M containing values of the basic
                variables.  (Output)
       IBASIS - Vector of length M+NVAR; elements 1-M point to basic
                columns
                others to nonbasic columns.  (Output)
       IBB    - Vector of length M+NVAR which signifies at which bound,
                if any, a nonbasic variable is set, or that the variable
                is basic.  (Output)
       ISTAT  - Vector of length M+N indicating the bound types of the
                variables.  (Input)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1988 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_d4prs(Mint m, Mint nvar, Mfloat *a, Mint lda, 
                    Mfloat *xlb, Mfloat *xub, Mint artif, Mfloat *b,
                    Mfloat *xb, Mint *ibasis, Mint *ibb, Mint *istat)
#else
static void l_d4prs(m, nvar, a, lda, xlb, xub, artif, b, xb, ibasis,
	            ibb, istat)
	Mint            m, nvar;
	Mfloat          *a;
	Mint            lda;
	Mfloat           xlb[], xub[];
	Mint            artif;
	Mfloat          *b, xb[];
	Mint             ibasis[], ibb[], istat[];
#endif
{
#define A(I_,J_)	(a+(I_)*(lda)+(J_))
#define B(I_,J_)	(b+(I_)*(m)+(J_))
	Mint             i, j, k;
	Mfloat           sum;


	imsl_e1psh("l_d4prs ");
	if (artif) {
		/*
		 * Set up initial B, representing a full artificial basis, B
		 * = -I.
		 */
		sset(m * m, F_ZERO, b, 1);
		sset(m, -F_ONE, B(0, 0), m + 1);
		/*
		 * Set up the status vector ISTAT, whose elements signify
		 * whether a variable is: 0 - basic 1 - nonbasic, at lower
		 * bound 2 - nonbasic, at upper bound 3 - nonbasic, no bound
		 * (arbitrarily set to zero) Initially the basic variables
		 * are x(n+1),...,x(n+m) and the nonbasic x(1),...,x(n).
		 */
		for (i = 1; i <= m; i++) {
			k = nvar + i;
			ibasis[i - 1] = k;
			ibb[k - 1] = -1;
		}
		for (j = 1; j <= nvar; j++) {
			k = m + j;
			ibasis[k - 1] = j;
			if (istat[j - 1] == 1) {
				ibb[j - 1] = 1;
			} else if (istat[j - 1] == 2) {
				ibb[j - 1] = 2;
			} else if (istat[j - 1] == 3) {
				if (fabs(xlb[j - 1]) <= fabs(xub[j - 1])) {
					ibb[j - 1] = 1;
				} else {
					ibb[j - 1] = 2;
				}
			} else {
				ibb[j - 1] = 3;
			}
		}
		/*
		 * 0.  Compute basic variable values. Initially, with a full
		 * artificial basis, we have B = -I, N = A, whence B xb = -N
		 * xn implies xb = A xn.
		 */
		for (i = 1; i <= m; i++) {
			sum = F_ZERO;
			for (j = 1; j <= nvar; j++) {
				if (ibb[j - 1] == 1) {
					sum += *A(j - 1, i - 1) * xlb[j - 1];
				} else if (ibb[j - 1] == 2) {
					sum += *A(j - 1, i - 1) * xub[j - 1];
				}
			}
			xb[i - 1] = sum;
		}

	} else {
		/* Future alternate start scheme(s). */
	}
	imsl_e1pop("l_d4prs ");
	return;
}				/* end of function */
/*Translated by FOR_C++, v0.1, on 08/09/90 at 10:47:11 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/09/90 at 10:47:09
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  D5PRS/DD5PRS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 5, 1989

    Purpose:    Compute the right hand side for the current basis.

    Usage:      CALL D5PRS (M, N, A, LDA, XLB, XUB, BL, BU, ISTAT,
                            IBB, RHS)

    Arguments:
       M      - Number of constraints.  (Input)
       N      - Number of variables.  (Input)
       A      - Matrix of dimension M by N containing the
                coefficients of the M constraints.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
                LDA must be at least M.
       XLB    - Vector of length N containing the lower bounds on
                the variables; if there is no lower bound on a variable,
                then 1.0E30 should be set as the lower bound.  (Input)
       XUB    - Vector of length N containing the upper bounds on
                variables; if there is no upper bound on a variable,
                then -1.0E30 should be set as the upper bound.  (Input)
       BL     - Vector of length M containing the lower limit of the
                general constraints; if there is no lower limit on the
                I-th constraint, then BL(I) would not be referred.
                (Input)
       BU     - Vector of length M containing the upper limit of the
                general constraints; if there is no upper limit on the
                I-th constraint, then BU(I) would not be referred;
                if there is no range constraint, BL and BU can share
                the same storage locations.  (Input)
       ISTAT  - Vector of length M+N indicating the bound types of the
                variables.
       IBB    - Vector of length M+N which signifies at which bound,
                if any, a nonbasic variable is set, or that the variable
                is basic.  (Input)
       RHS    - Right-hand side; vector of length m.  (Output)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1988 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_d5prs(Mint m, Mint n, Mfloat *a, Mint lda, Mfloat *xlb,
                    Mfloat *xub, Mfloat *bl, Mfloat *bu, Mint *istat,
                    Mint *ibb, Mfloat *rhs)
#else
static void l_d5prs(m, n, a, lda, xlb, xub, bl, bu, istat, ibb, rhs)
	Mint            m, n;
	Mfloat          *a;
	Mint            lda;
	Mfloat           xlb[], xub[], bl[], bu[];
	Mint             istat[], ibb[];
	Mfloat           rhs[];
#endif
{
#define A(I_,J_)	(a+(I_)*(lda)+(J_))
	Mint             j, j1;
	Mfloat           scalr;


	imsl_e1psh("l_d5prs ");

	/*
	 * Initially set right-hand side vector to zero.
	 */
	sset(m, F_ZERO, rhs, 1);

	/*
	 * Translate RHS according to classification of variables.
	 */
	for (j = 1; j <= n; j++) {
		if (istat[j - 1] == 4)
			goto L_10;
		if (ibb[j - 1] < 0)
			goto L_10;

		if (ibb[j - 1] == 1) {
			scalr = xlb[j - 1];
		} else if (ibb[j - 1] == 2) {
			scalr = xub[j - 1];
		}
		if (scalr != F_ZERO) {
			saxpy(m, -scalr, A(j - 1, 0), 1, rhs, 1);
		}
L_10:
		;
	}

	for (j = n + 1; j <= (n + m); j++) {
		if (istat[j - 1] == 4)
			goto L_20;
		if (ibb[j - 1] < 0)
			goto L_20;

		j1 = j - n;
		if (ibb[j - 1] == 1) {
			scalr = bl[j1 - 1];
		} else if (ibb[j - 1] == 2) {
			scalr = bu[j1 - 1];
		}
		if (scalr != F_ZERO) {
			rhs[j1 - 1] += scalr;
		}
L_20:
		;
	}

	imsl_e1pop("l_d5prs ");
	return;
}				/* end of function */
/*Translated by FOR_C++, v0.1, on 08/09/90 at 10:47:57 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/09/90 at 10:47:55
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  D6PRS/DD6PRS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 5, 1989

    Purpose:    Solve a real general system of linear equations given the
                LU factors and ETA vectors of the coefficient matrix.

    Usage:      CALL D6PRS (N, FAC, LDFAC, IPVT, IETA, ETA, ICOUNT, B,
                            IPATH, X, Y)

    Arguments:
       N      - Number of variables.  (Input)
       FAC    - M by M matrix containing the LU factorization of the
                basis matrix as output from subroutine
                LFCRG/DLFCRG or LFTRG/DLFTRG.  (Input)
       LDFAC  - Leading dimension of FAC exactly as specified in the
                dimension statement of the calling program.  (Input)
       IPVT   - Vector of length M containing the pivoting information
                for the LU factorization of B as output from subroutine
                LFCRG/DLFCRG or LFTRG/DLFTRG.  (Input)
       IETA   - Vector of length 25*m containing pivoting information
                regarding the eta vectors.
       ETA    - Vector of length 25*m representing the sequence of
                elementary matrices used in the product form of the
                inverse of B.
       ICOUNT - Counter used to determine when to reinvert the basis
                matrix B.
       B      - The right-hand side of the linear system.
                (Input)
       IPATH  - Indicates whether to use a coefficient matrix or its
                transpose when solving a linear system.
       X      - Vector of length N containing the solution to the linear
                system.  (Output)  If B is not needed, B and X can share
                the same storage locations.
       Y      - Vector of length m used for intermediate solutions
                to linear systems.

    Chapters:   MATH/LIBRARY Linear Systems
                STAT/LIBRARY Mathematical Support

    Copyright:  1988 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_d6prs(Mint n, Mfloat *imsl_fac, Mint ldfac, Mint *ipvt,
                    Mint *ieta, Mfloat *eta, Mint icount, Mfloat *b,
	            Mint ipath, Mfloat *x, Mfloat *y)
#else
static void l_d6prs(n, imsl_fac, ldfac, ipvt, ieta, eta, icount, b,
	   ipath, x, y)
	Mint            n;
	Mfloat          *imsl_fac;
	Mint            ldfac, ipvt[], ieta[];
	Mfloat           eta[];
	Mint             icount;
	Mfloat           b[];
	Mint            ipath;
	Mfloat           x[], y[];
#endif
{
#define FAC(I_,J_)	(imsl_fac+(I_)*(ldfac)+(J_))
	Mint             i, i1, i2;
	Mfloat           pivt, tmp;

	imsl_e1psh("l_d6prs ");

	if (ipath == 1) {

		/*
		 * Solve An*x = bn; first solve L*y = p*bn.
		 */
		l_d14rs(n, imsl_fac, ldfac, ipvt, b, ipath, y);
		/* Compute y = ln*pn...y. */
		for (i = 1; i <= (icount - 1); i++) {
			i1 = ieta[i - 1];
			if (i1 > 0) {
				tmp = y[i1 - 1];
				y[i1 - 1] = y[i1];
				y[i1] = tmp;
			}
			pivt = eta[i - 1];
			if (pivt != F_ZERO) {
				i2 = abs(i1);
				y[i2] += y[i2 - 1] * pivt;
			}
		}

		/*
		 * Solve Un*x = y.
		 */
	/*	imsl_lslrt(n, imsl_fac, ldfac, y, ADR(_l0, 2), x);*/
                scopy(n, y, 1, x, 1);
                imsl_strsv("u", "n", "n", n, imsl_fac, ldfac, x, 1);

	} else {

		/*
		 * Solve tran(An)*x = b; first solve tran(Un)*x = b.
		 */
	      /* imsl_lslrt(n, imsl_fac, ldfac, b, ADR(_l0, 4), x);*/
                scopy(n, b, 1, x, 1);
                imsl_strsv("u", "t", "n", n, imsl_fac, ldfac, x, 1);
		/*
		 * Compute tran(x)=tran(x)*ln*pn.
		 */
		for (i = icount - 1; i >= 1; i--) {
			pivt = eta[i - 1];
			i1 = ieta[i - 1];
			if (pivt != F_ZERO) {
				i2 = abs(i1);
				x[i2 - 1] += x[i2] * pivt;
			}
			if (i1 > 0) {
				tmp = x[i1 - 1];
				x[i1 - 1] = x[i1];
				x[i1] = tmp;
			}
		}

		/*
		 * Solve tran(L)*x=x.
		 */
		l_d14rs(n, imsl_fac, ldfac, ipvt, x, 2, x);
	}

	imsl_e1pop("l_d6prs ");
	return;
}				/* end of function */
/*Translated by FOR_C++, v0.1, on 08/09/90 at 10:48:35 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/09/90 at 10:48:32
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  D7PRS/DD7PRS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 5, 1989

    Purpose:    Check the feasibility of the current point.

    Usage:      CALL D7PRS (M, N, XLB, XUB, BL, BU, ISTAT, CB, XB,
                            IBASIS, TOLX, TINF, FEAS)

    Arguments:
       M      - Number of constraints.  (Input)
       N      - Number of variables.  (Input)
       XLB    - Vector of length N containing the lower bounds on
                the variables; if there is no lower bound on a variable,
                then 1.0E30 should be set as the lower bound.  (Input)
       XUB    - Vector of length N containing the upper bounds on
                variables; if there is no upper bound on a variable,
                then -1.0E30 should be set as the upper bound.  (Input)
       BL     - Vector of length M containing the lower limit of the
                general constraints; if there is no lower limit on the
                I-th constraint, then BL(I) would not be referred.
                (Input)
       BU     - Vector of length M containing the upper limit of the
                general constraints; if there is no upper limit on the
                I-th constraint, then BU(I) would not be referred;
                if there is no range constraint, BL and BU can share
                the same storage locations.  (Input)
       ISTAT  - Vector of length M+N indicating the bound types of the
                variables.
       CB     - Vector of length M containing objective coefficients.
                (Output)
       XB     - Vector of length M containing values of the basic
                variables.  (Input)
       IBASIS - Vector of length M+N; elements 1-M point to basic
                columns
                others to nonbasic columns.  (Input)
       TOLX   - Tolerance for primal variables.  (Input)
       TINF   - Endpoint tolerance.
       FEAS   - Feasibility flag.  (Output)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1988 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_d7prs(Mint m, Mint n, Mfloat *xlb, Mfloat *xub,
                    Mfloat *bl, Mfloat *bu, Mint *istat, Mfloat *cb,
                    Mfloat *xb, Mint *ibasis, Mfloat tolx, 
                    Mfloat *tinf, Mint *feas)
#else
static void l_d7prs(m, n, xlb, xub, bl, bu, istat, cb, xb, ibasis,
	            tolx, tinf, feas)
	Mint            m, n;
	Mfloat           xlb[], xub[], bl[], bu[];
	Mint             istat[];
	Mfloat           cb[], xb[];
	Mint             ibasis[];
	Mfloat           tolx, *tinf;
	Mint           *feas;
#endif
{
	Mint             i, j, j1;
	Mfloat           t, xbi;


	imsl_e1psh("l_d7prs ");

	/*
	 * Define the classification of the basic variables:  -1 violates
	 * lower bound; 0 is feasible; +1 violates upper bound.  This
	 * information is stored in cb(1)-cb(m).
	 */
	*feas = 1;
	*tinf = F_ZERO;
	sset(m, F_ZERO, cb, 1);
	for (i = 1; i <= m; i++) {
		j = ibasis[i - 1];
		if (istat[j - 1] == 4)
			goto L_20;
		xbi = xb[i - 1];
		if (j > n)
			goto L_10;
		if (istat[j - 1] == 1) {
			t = xlb[j - 1] - xbi;
			if (t > tolx) {
				cb[i - 1] = -F_ONE;
				*feas = 0;
				*tinf += t;
			}
		} else if (istat[j - 1] == 2) {
			t = xbi - xub[j - 1];
			if (t > tolx) {
				cb[i - 1] = F_ONE;
				*feas = 0;
				*tinf += t;
			}
		} else {
			t = xlb[j - 1] - xbi;
			if (t > tolx) {
				cb[i - 1] = -F_ONE;
				*feas = 0;
				*tinf += t;
			} else {
				t = xbi - xub[j - 1];
				if (t > tolx) {
					cb[i - 1] = F_ONE;
					*feas = 0;
					*tinf += t;
				}
			}
		}
		goto L_20;

L_10:
		j1 = j - n;
		if (istat[j - 1] == 1) {
			t = bl[j1 - 1] - xbi;
			if (t > tolx) {
				cb[i - 1] = -F_ONE;
				*feas = 0;
				*tinf += t;
			}
		} else if (istat[j - 1] == 2) {
			t = xbi - bu[j1 - 1];
			if (t > tolx) {
				cb[i - 1] = F_ONE;
				*feas = 0;
				*tinf += t;
			}
		} else {
			t = bl[j1 - 1] - xbi;
			if (t > tolx) {
				cb[i - 1] = -F_ONE;
				*feas = 0;
				*tinf += t;
			} else {
				t = xbi - bu[j1 - 1];
				if (t > tolx) {
					cb[i - 1] = F_ONE;
					*feas = 0;
					*tinf += t;
				}
			}
		}
L_20:
		;
	}
	imsl_e1pop("l_d7prs ");
	return;
}				/* end of function */
/*Translated by FOR_C++, v0.1, on 08/09/90 at 10:49:07 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/09/90 at 10:49:06
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  D8PRS/DD8PRS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 5, 1989

    Purpose:    Price out the nonbasic columns.

    Usage:      CALL D8PRS (M, NVAR, A, LDA, C, PI, IPHASE, IBASIS, D)

    Arguments:
       M      - Number of constraints.  (Input)
       NVAR   - Number of variables.  (Input)
       A      - Matrix of dimension M by NVAR containing the
                coefficients of the M constraints.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
                LDA must be at least M.
       C      - Objective coefficients.  (Input)
       PI     - Vector of length M containing the simplex multipliers.
                (Input)
       IPHASE - Indicates Phase 1 or Phase 2 of the revised simplex
                algorithm.  (Input)
       IBASIS - Vector of length M+NVAR; elements 1-M point to basic
                columns, others to nonbasic columns.  (Input)
       D      - Vector of length NVAR containing reduced costs.
                (Output)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1988 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_d8prs(Mint m, Mint nvar, Mfloat *a, Mint lda, 
                    Mfloat *c, Mfloat *pi, Mint iphase, Mint *ibasis,
                    Mfloat *d)
#else
static void l_d8prs(m, nvar, a, lda, c, pi, iphase, ibasis, d)
	Mint            m, nvar;
	Mfloat          *a;
	Mint            lda;
	Mfloat           c[], pi[];
	Mint             iphase, ibasis[];
	Mfloat           d[];
#endif
{
#define A(I_,J_)	(a+(I_)*(lda)+(J_))
	Mint             j, k;


	imsl_e1psh("l_d8prs ");

	for (j = 1; j <= nvar; j++) {
		k = ibasis[m + j - 1];
		if (k > nvar) {
			d[j - 1] = pi[k - nvar - 1];
		} else {
			if (iphase == 2) {
				d[j - 1] = c[k - 1] - imsl_sdot(m, pi, 1, A(k - 1, 0),
								1);
			} else {
				d[j - 1] = -imsl_sdot(m, pi, 1, A(k - 1, 0), 1);
			}
		}
	}
	imsl_e1pop("l_d8prs ");
	return;
}				/* end of function */
/*Translated by FOR_C++, v0.1, on 08/09/90 at 10:49:41 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/09/90 at 10:49:38
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  D9PRS/DD9PRS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 5, 1989

    Purpose:    Find variable to enter basis.

    Usage:      CALL D9PRS (MRELAS, NVARS, IENTER, IBASIS, ISTAT, IBB,
                            XLB, XUB, BL, BU, EPS, TOLDJ, PINRM, D,
                            FOUND)

    Arguments:
       MRELAS - Number of constraints.  (Input)
       NVARS  - Number of variables.  (Input)
       IENTER - Next column to enter basis.  (Output)
       IBASIS - Vector of length M+NVARS; elements 1-M point to basic
                columns, others to nonbasic columns.  (Input)
       ISTAT  - Vector of length M+N indicating the bound types of the
                variables.
       IBB    - Indicator for non-basic variables, polarity of variables,
                and potentially infinite variables.
       XLB    - Vector of length NVARS containing the lower bound on the
                variables.  (Input)
       XUB    - Vector of length NVARS containing the upper bound on the
                variables.  (Input)
       BL     - Vector of length M containing the lower limit of the
                general constraints; if there is no lower limit on the
                I-th constraint, then BL(I) would not be referred.
                (Input)
       BU     - Vector of length M containing the upper limit of the
                general constraints; if there is no upper limit on
                the I-th constraint, then BU(I) would not be referred;
                if there is no range constraint, BL and BU can share
                the same storage locations.  (Input)
       EPS    - Tolerance.
       TOLDJ  - Tolerance for reduced costs.
       PINRM  - Norm of simplex multiplier vector.
       D      - Vector of length NVARS containing reduced costs.  (Input)
       FOUND  - Flag to indicate whether an entering variable has been
                found.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1988 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_d9prs(Mint mrelas, Mint nvars, Mint *ienter, 
                    Mint *ibasis, Mint *istat, Mint *ibb,
	            Mfloat *xlb, Mfloat *xub, Mfloat *bl, Mfloat *bu,
                    Mfloat eps, Mfloat toldj, Mfloat pinrm, Mfloat *d,
                    Mint *found)
#else
static void l_d9prs(mrelas, nvars, ienter, ibasis, istat, ibb,
                    xlb, xub, bl, bu, eps, toldj, pinrm, d, found)
	Mint            mrelas, nvars, *ienter, ibasis[], istat[], ibb[];
	Mfloat           xlb[], xub[], bl[], bu[], eps, toldj, pinrm,
	                d[];
	Mint           *found;
#endif
{
	Mint             i, j;
	Mfloat           blj, buj, ratio, rcost, rmax;


	imsl_e1psh("l_d9prs ");
	rmax = F_ZERO;
	*found = 0;
	for (i = mrelas + 1; i <= (mrelas + nvars); i++) {
		j = ibasis[i - 1];
		/*
		 * If J=IBASIS(I) < 0, then the variable is left at a zero
		 * level and is not considered as a candidate to enter.
		 */
		if (j > 0) {

			/*
			 * Do not consider variables which correspond to
			 * unbounded step lengths.
			 */
			if (ibb[j - 1] == 0)
				goto L_10;
			/*
			 * Do not consider equations as candidates to enter
			 * the basis.
			 */
			if (istat[j - 1] == 3) {
				if (j > nvars) {
					buj = bu[j - nvars - 1];
					blj = bl[j - nvars - 1];
				} else {
					buj = xub[j - 1];
					blj = xlb[j - 1];
				}
				if (buj - blj <= eps * (fabs(buj) + fabs(blj)))
					goto L_10;
			}
			/*
			 * Compute the reduced cost.
			 */
			rcost = d[i - mrelas - 1];
			/*
			 * If variable is at upper bound, it can only
			 * decrease.  This accounts for the possible change
			 * of sign.
			 */
			if (mod(ibb[j - 1], 2) == 0)
				rcost = -rcost;

			/*
			 * If the variable is free, use the negative
			 * magnitude of the reduced cost for that variable.
			 */
			if (istat[j - 1] == 4)
				rcost = -fabs(rcost);

			/*
			 * Test for negativity of reduced costs.
			 */
			if (rcost + pinrm * toldj < F_ZERO) {
				ratio = fabs(rcost);
				if (ratio > rmax) {
					*found = 1;
					rmax = ratio;
					*ienter = i;
				}
			}
		}
L_10:
		;
	}

	imsl_e1pop("l_d9prs ");
	return;
}				/* end of function */
/*Translated by FOR_C++, v0.1, on 08/09/90 at 10:50:26 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/09/90 at 10:50:22
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  D10RS/DD10RS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 5, 1989

    Purpose:    Choose variable to leave the basis.

    Usage:      CALL D10RS (MRELAS, NVARS, IENTER, ILEAVE, IBASIS, ISTAT,
                            THETA, DIRNRM, TOLN, ALPHAQ, XLB, XUB,
                            BL, BU, XB, FINITE)

    Arguments:
       MRELAS - Number of constraints.  (Input)
       NVARS  - Number of variables.  (Input)
       IENTER - Next column to enter basis.  (Input)
       ILEAVE - Next column to leave basis.  (Output)
       IBASIS - Vector of length M+NVARS; elements 1-M point to basic
                columns, others to nonbasic columns.  (Input)
       ISTAT  - Vector of length M+N indicating the bound types of the
                variables.
       THETA  - Scaling variable.  (Output)
       DIRNRM - Norm of the direction vector.
       TOLN   - Tolerance used in connection with the direction
                component.
       ALPHAQ - Work array of length MRELAS.
       XLB    - Vector of length NVARS containing the lower bound on the
                variables.  (Input)
       XUB    - Vector of length NVARS containing the upper bound on the
                variables.  (Input)
       BL     - Vector of length M containing the lower limit of the
                general constraints; if there is no lower limit on the
                I-th constraint, then BL(I) would not be referred.
                (Input)
       BU     - Vector of length M containing the upper limit of the
                general constraints; if there is no upper limit on
                the I-th constraint, then BU(I) would not be referred;
                if there is no range constraint, BL and BU can share
                the same storage locations.  (Input)
       XB     - Vector of length M containing values of the basic
                variables.  (Input)
       FINITE - Logical variable.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1988 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_d10rs(Mint mrelas, Mint nvars, Mint ienter, 
                    Mint *ileave, Mint *ibasis, Mint *istat,
	            Mfloat *theta, Mfloat dirnrm, Mfloat toln, 
                    Mfloat *alphaq, Mfloat *xlb, Mfloat *xub,
                    Mfloat *bl, Mfloat *bu, Mfloat *xb, Mint *finite)
#else
static void l_d10rs(mrelas, nvars, ienter, ileave, ibasis, istat,
                    theta, dirnrm, toln, alphaq, xlb, xub, bl, bu,
                    xb, finite)
	Mint            mrelas, nvars, ienter, *ileave, ibasis[], istat[];
	Mfloat          *theta, dirnrm, toln, alphaq[], xlb[], xub[],
	                bl[], bu[], xb[];
	Mint           *finite;
#endif
{
	Mint             i, ileav1, j, j1;
	Mfloat           d, ratio, tmpp, theta1;


	imsl_e1psh("l_d10rs ");
	*theta = 1.0e36;
	theta1 = F_ZERO;

	/*
	 * See if the entering variable is restricting the step length
	 * because of an upper bound.
	 */
	*finite = 1;
	j = ibasis[ienter - 1];
	if (istat[j - 1] == 3) {
		if (j > nvars) {
			j1 = j - nvars;
			*theta = bu[j1 - 1] - bl[j1 - 1];
		} else {
			*theta = xub[j - 1] - xlb[j - 1];
		}
		*ileave = ienter;
	}
	/*
	 * Now use the basic variables to possibly restrict the step length
	 * even further.
	 */
	for (i = 1; i <= mrelas; i++) {
		j = ibasis[i - 1];

		/*
		 * If this is a free variable, do not use it to restrict the
		 * step length.
		 */
		if (istat[j - 1] == 4) {
			goto L_10;

		}
		/*
		 * If direction component is about zero, ignore it for
		 * computing the step length.
		 */
		ratio = alphaq[i - 1];
		if (fabs(ratio) <= dirnrm * toln) {
			goto L_10;

		}
		if (ratio > F_ZERO) {

			/* If no lower bound go to end of loop. */
			if (istat[j - 1] == 2) {
				if (j > nvars) {
					d = xb[i - 1] - bu[j - nvars - 1];
				} else {
					d = xb[i - 1] - xub[j - 1];
				}
				if (d <= toln)
					goto L_10;
				ratio = d / ratio;
				if (theta1 < ratio) {
					theta1 = ratio;
					ileav1 = i;
				}
				goto L_10;
			}
			if (j > nvars) {
				d = xb[i - 1] - bl[j - nvars - 1];
			} else {
				d = xb[i - 1] - xlb[j - 1];
			}
			if (d <= -toln)
				goto L_10;
		} else {
			if (istat[j - 1] == 1) {
				if (j > nvars) {
					d = xb[i - 1] - bl[j - nvars - 1];
				} else {
					d = xb[i - 1] - xlb[j - 1];
				}
				if (d >= -toln)
					goto L_10;
				ratio = d / ratio;
				if (theta1 < ratio) {
					theta1 = ratio;
					ileav1 = i;
				}
				goto L_10;
			}
			if (j > nvars) {
				d = xb[i - 1] - bu[j - nvars - 1];
			} else {
				d = xb[i - 1] - xub[j - 1];
			}
			if (d >= toln)
				goto L_10;
		}
		ratio = d / ratio;
		if (ratio < F_ZERO)
			ratio = F_ZERO;
		if (*theta < ratio) goto L_10;
                if (*theta==ratio &&
		fabs(alphaq[*ileave-1])>=fabs(alphaq[i-1])) goto L_10;
		*theta = ratio;
		*ileave = i;
L_10:
		;
	}
        tmpp = 1.0e36 - imsl_amach(4);
	if (*theta >= tmpp) {
		if (theta1 == F_ZERO) {
			*finite = 0;
		} else {
			*ileave = ileav1;
			*theta = theta1;
		}
	}
	imsl_e1pop("l_d10rs ");
	return;
}				/* end of function */
/*Translated by FOR_C++, v0.1, on 08/09/90 at 10:51:00 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/09/90 at 10:50:57
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  D11RS/DD11RS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 5, 1989

    Purpose:    Updates the primal solution, basis indicators, and
                variable status indicators.

    Usage:      CALL D11RS (M, N, IENTER, ILEAVE, IBASIS, ISTAT, IBB,
                            THETA, XB, ALPHAQ, XLB, XUB, BL, BU, TOLX,
                            CHANGE)

    Arguments:
       M      - Number of constraints.  (Input)
       N      - Number of variables.  (Input)
       IENTER - Next column to enter the basis.  (Input)
       ILEAVE - Next column to leave the basis.  (Input)
       IBASIS - Vector of length M+N; elements 1-M point to basic
                columns
                others to nonbasic columns.  (Input)
       ISTAT  - Vector of length M+N indicating the bound types of the
                variables.
       IBB    - Indicator for non-basic variables, polarity of variables,
                and potentially infinite variables.
       THETA  - Scaling variable.  (Input)
       XB     - Internal primal solution vector of length M.
       ALPHAQ - Work array of length M.
       XLB    - Vector of length NVAR containing the lower bound on the
                variables; if there is no lower bound on a variable,
                then 1.0E30 should be set as the lower bound.  (Input)
       XUB    - Vector of length NVAR containing the upper bound on the
                variables; if there is no upper bound on a variable,
                then -1.0E30 should be set as the upper bound.  (Input)
       BL     - Vector of length M containing the lower limit of the
                general constraints; if there is no lower limit on the
                I-th constraint, then BL(I) would not be referred.
                (Input)
       BU     - Vector of length M containing the upper limit of the
                general constraints; if there is no upper limit on
                the I-th constraint, then BU(I) would not be referred;
                if there is no range constraint, BL and BU can share
                the same storage locations.  (Input)
       TOLX   - Tolerance for primal variables.  (Input)
       CHANGE - Logical variable indicating whether a basis change is
                needed.  (Output)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1988 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_d11rs(Mint m, Mint n, Mint ienter, Mint ileave, 
                    Mint *ibasis, Mint *istat, Mint *ibb, 
                    Mfloat theta, Mfloat *xb, Mfloat *alphaq,
                    Mfloat *xlb, Mfloat *xub, Mfloat *bl, Mfloat *bu,
                    Mfloat tolx,Mfloat * d, Mint *change)
#else
static void l_d11rs(m, n, ienter, ileave, ibasis, istat, ibb,
                    theta, xb, alphaq, xlb, xub, bl, bu, tolx, d,
                    change)
	Mint            m, n, ienter, ileave, ibasis[], istat[], ibb[];
	Mfloat          theta, xb[], alphaq[], xlb[], xub[], bl[], bu[],
	               tolx, d[];
	Mint           *change;
#endif
{
	Mint             j, j1;
	Mfloat           t;


	imsl_e1psh("l_d11rs ");
	*change = 1;

	/*
	 * Update the primal solution with a multiple of the search
	 * direction.
	 */
	saxpy(m, -theta, alphaq, 1, xb, 1);

	/*
	 * Check feasibility.  If IENTER equals ILEAVE, then move the
	 * entering column to the other bound; no basis change needed.
	 */
	if (ienter == ileave) {
		*change = 0;
		j = ibasis[ienter - 1];
		if (ibb[j - 1] == 1) {
			ibb[j - 1] = 2;
		} else {
			ibb[j - 1] = 1;
		}
		goto L_10;
	}
	/*
	 * Set the leaving variable to be the nonbasic variable.
	 */
	j = ibasis[ileave - 1];
	if (istat[j - 1] == 1) {
		ibb[j - 1] = 1;
	} else if (istat[j - 1] == 2) {
		ibb[j - 1] = 2;
	} else {
		if (j > n) {
			t = xb[ileave - 1] - bl[j - n - 1];
		} else {
			t = xb[ileave - 1] - xlb[j - 1];
		}
		if (fabs(t) < tolx) {
			ibb[j - 1] = 1;
		} else {
			ibb[j - 1] = 2;
		}
	}

	/* Update xb. */
	scopy(m - ileave, &xb[ileave], 1, &xb[ileave - 1], 1);

	/*
	 * Put the value in for the new basic variable IENTER.
	 */
	j = ibasis[ienter - 1];
	if (ibb[j - 1] == 3) {
		if (d[ienter - m - 1] > F_ZERO) {
			xb[m - 1] = -theta;
		} else {
			xb[m - 1] = theta;
		}
		goto L_15;
	}
	if (j > n) {
		j1 = j - n;
		if (ibb[j - 1] == 1) {
			xb[m - 1] = bl[j1 - 1] + theta;
		} else {
			xb[m - 1] = bu[j1 - 1] - theta;
		}
	} else {
		if (ibb[j - 1] == 1) {
			xb[m - 1] = xlb[j - 1] + theta;
		} else {
			xb[m - 1] = xub[j - 1] - theta;
		}
	}
L_15:
	ibb[j - 1] = -1;

	/*
	 * Update column pointers to note exchange of columns.
	 */
	ibasis[ienter - 1] = ibasis[abs(ileave) - 1];
	icopy(m - ileave, &ibasis[ileave], 1, &ibasis[ileave - 1],
		   1);
	ibasis[m - 1] = j;

	/*
	 * If variable was exchanged at a zero level, mark it so that it
	 * can't be brought back in.  This is to help prevent cycling.
	 * Update the matrix decomposition; column ABS(ILEAVE) is leaving.
	 */
L_10:
	;
	imsl_e1pop("l_d11rs ");
	return;
}				/* end of function */
/*Translated by FOR_C++, v0.1, on 08/09/90 at 10:51:34 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/09/90 at 10:51:31
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  D12RS/DD12RS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 5, 1989

    Purpose:    Checks errors to see if reinversion is needed.

    Usage:      CALL D12RS (M, NVARS, A, LDA, IBASIS, TOLX, EPS, WK,
                            XB, RHS, BIGERR)

    Arguments:
       M      - Number of constraints.  (Input)
       NVARS  - Number of variables.  (Input)
       A      - Matrix of dimension M by N containing the
                coefficients of the M constraints.  (Input)
       LDA    - Leading dimension of A exactly as specified in the
                dimension statement of the calling program.  (Input)
                LDA must be at least M.
       IBASIS - Vector of length M+NVARS; elements 1-M point to basic
                columns, others to nonbasic columns.  (Input)
       TOLX   - Tolerance for primal variables.  (Input)
       EPS    - Tolerance.  (Input)
       WK     - Work vector.  (Output)
       XB     - Vector of length M containing values of the basic
                variables.  (Input)
       RHS    - Right-hand side; vector of length m.  (Input)
       BIGERR - Flag indicating singular system.  (Output)

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1988 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_d12rs(Mint m, Mint nvars, Mfloat *a, Mint lda, 
                    Mint *ibasis, Mfloat tolx, Mfloat eps, Mfloat *wk,
                    Mfloat *xb, Mfloat *rhs, Mint *bigerr)
#else
static void l_d12rs(m, nvars, a, lda, ibasis, tolx, eps, wk, xb,
	            rhs, bigerr)
	Mint            m, nvars;
	Mfloat          *a;
	Mint            lda, ibasis[];
	Mfloat          tolx, eps, wk[], xb[], rhs[];
	Mint           *bigerr;
#endif
{
#define A(I_,J_)	(a+(I_)*(lda)+(J_))
	Mint             i, j, l;
	Mfloat           erpmax, xbl;


	imsl_e1psh("l_d12rs ");

	scopy(m, rhs, 1, wk, 1);
	for (l = 1; l <= m; l++) {
		xbl = xb[l - 1];
		if (xbl <= eps)
			goto L_10;
		j = ibasis[l - 1];
		if (j > nvars) {
			i = j - nvars;
			wk[i - 1] += xbl;

		} else {
			saxpy(m, -xbl, A(j - 1, 0), 1, wk, 1);
		}

L_10:
		;
	}
	erpmax = F_ZERO;
	for (i = 1; i <= m; i++) {
		erpmax = imsl_f_max(fabs(wk[i - 1] / (fabs(rhs[i - 1]) + F_ONE)),
				    erpmax);
	}

	/*
	 * System becomes singular when the accuracy of the solution is
	 * greater than TOLX, which may need to be changed on some machines.
	 */
	*bigerr = erpmax > tolx;

	imsl_e1pop("l_d12rs ");
	return;
}				/* end of function */
/*Translated by FOR_C++, v0.1, on 08/09/90 at 10:52:08 */
/*FOR_C++ Options SET: cio no=dkrp op=an pf=/home/usr2/imsl1/clib/newclib/include/imsl_int.h c - prototypes */
/* Structured by FOR_STRUCT, v0.2, on 08/09/90 at 10:52:07
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  D13RS/DD13RS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 5, 1989

    Purpose:    Update U and ETA vectors.

    Usage:      CALL D13RS (N, FAC, LDFAC, IPVT, UQ, ILEAVE, IETA,
                            ETA, ICOUNT, DIAG)

    Arguments:
       N      - Number of variables.  (Input)
       FAC    - Matrix of dimension M by M which contains the LU
                factorization of the basis matrix.  (Input)
       LDFAC  - Leading dimension of FAC exactly as specified in the
                dimension statement of the calling program.  (Input)
       IPVT   - Vector of length M containing the pivoting information
                for the LU factorization of B as output from subroutine
                LFCRG/DLFCRG or LFTRG/DLFTRG.  (Input)
       UQ     - Entering column of length N which will form the new U.
                (Input)
       ILEAVE - Next column to leave basis.  (Output)
       IETA   - Vector of length 25*m containing pivoting information
                regarding the eta vectors.
       ETA    - Vector of length 25*m representing the sequence of
                elementary matrices used in the product form of the
                inverse of B.
       ICOUNT - Counter used to determine when to reinvert the basis
                matrix B.
       DIAG   - Vector of diagonal elements.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1988 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_d13rs(Mint n, Mfloat *imsl_fac, Mint ldfac, Mint *ipvt,
                    Mfloat *uq, Mint ileave, Mint *ieta, Mfloat *eta,
                    Mint *icount, Mfloat *diag)
#else
static void l_d13rs(n, imsl_fac, ldfac, ipvt, uq, ileave, ieta, eta,
                    icount, diag)
	Mint            n;
	Mfloat          *imsl_fac;
	Mint            ldfac, ipvt[];
	Mfloat           uq[];
	Mint             ileave, ieta[];
	Mfloat           eta[];
	Mint            *icount;
	Mfloat           diag[];
#endif
{
#define FAC(I_,J_)	(imsl_fac+(I_)*(ldfac)+(J_))
	Mint             i1, j, k;
	Mfloat           pivt;


	imsl_e1psh("l_d13rs ");
	/* Uj <- U(j+1) */
	i1 = ileave + 1;
	scopy(n - ileave, FAC(i1 - 1, i1 - 1), ldfac + 1, diag, 1);
	for (j = ileave; j <= (n - 1); j++) {
		scopy(j, FAC(j, 0), 1, FAC(j - 1, 0), 1);
	}
	scopy(n, uq, 1, FAC(n - 1, 0), 1);
	k = 1;
	for (j = ileave; j <= (n - 1); j++) {
		if (fabs(*FAC(j - 1, j - 1)) >= fabs(diag[k - 1])) {
			ieta[*icount - 1] = -j;
			pivt = -diag[k - 1] / *FAC(j - 1, j - 1);
		} else {
			ieta[*icount - 1] = j;
			pivt = -*FAC(j - 1, j - 1) / diag[k - 1];
			*FAC(j - 1, j - 1) = diag[k - 1];
			sswap(n - j, FAC(j, j - 1), ldfac, FAC(j, j), ldfac);
		}
		if (pivt != F_ZERO) {
			saxpy(n - j, pivt, FAC(j, j - 1), ldfac, FAC(j, j), ldfac);
		}
		eta[*icount - 1] = pivt;
		*icount += 1;
		k += 1;
	}

	imsl_e1pop("l_d13rs ");
	return;
}				/* end of function */


/* -----------------------------------------------------------------------
    IMSL Name:  D14RS/DD14RS (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    January 5, 1989

    Purpose:    Solve L*y = b or Trans(L)*y = b.

    Usage:      CALL D14RS (N, FAC, LDFAC, IPVT, B, IPATH, X)

    Arguments:
       N      - Number of variables.  (Input)
       FAC    - Matrix of dimension M by M which contains the LU
                factorization of the basis matrix.  (Input)
       LDFAC  - Leading dimension of FAC exactly as specified in the
                dimension statement of the calling program.  (Input)
       IPVT   - Vector of length M containing the pivoting information
                for the LU factorization of B as output from subroutine
                LFCRG/DLFCRG or LFTRG/DLFTRG.  (Input)
       B      - Right-hand side of the linear system; vector of
                length N.  (Input)
       IPATH  - Indicates whether to use a coefficient matrix or its
                transpose when solving a linear system.
       X      - Preserves RHS of linear system.

    Chapter:    MATH/LIBRARY Optimization

    Copyright:  1988 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_d14rs(Mint n, Mfloat *imsl_fac, Mint ldfac, Mint *ipvt,
                    Mfloat *b, Mint ipath, Mfloat *x)
#else
static void l_d14rs(n, imsl_fac, ldfac, ipvt, b, ipath, x)
	Mint            n;
	Mfloat          *imsl_fac;
	Mint            ldfac, ipvt[];
	Mfloat           b[];
	Mint            ipath;
	Mfloat           x[];
#endif
{
#define FAC(I_,J_)	(imsl_fac+(I_)*(ldfac)+(J_))
	Mint             k, l;
	Mfloat           t;


	imsl_e1psh("l_d14rs ");

	/* Copy B into X and use X to preserve input.  */

	scopy(n, b, 1, x, 1);

	if (ipath == 1) {
	    /* IPATH = 1 , solve  L*Y = B. */
	    for (k = 1; k <= (n - 1); k++) {
		 l = ipvt[k - 1];
		 t = x[l - 1];
		 if (l != k) {
		     x[l - 1] = x[k - 1];
		     x[k - 1] = t;
		 }
		 saxpy(n - k, t, FAC(k - 1, k), 1, &x[k], 1);
	    }

	} else if (ipath == 2) {
	    /* IPATH = 2, solve TRANS(L)*X = B. */
	    for (k = n - 1; k >= 1; k--) {
		 x[k - 1] += imsl_sdot(n - k, FAC(k - 1, k), 1, &x[k], 1);
		 l = ipvt[k - 1];
		 if (l != k) {
		     t = x[l - 1];
		     x[l - 1] = x[k - 1];
		     x[k - 1] = t;
		 }
	    }

	} else {
	    imsl_e1sti(1, ipath);
/*	    imsl_ermes(5, 5, "IPATH must be either 1 or 2 while a value of  */
/*                            %(i1) is given.");                            */
            imsl_ermes(IMSL_TERMINAL, IMSL_IPATH_RANGE_3);
	}

	imsl_e1pop("l_d14rs ");
	return;
}				/* end of function */
