#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static void PROTO(l_s2rf,(Mint*, Mfloat*, Mfloat[], Mint*,
                   Mint*, Mfloat[], Mfloat[], Mfloat*,
                   Mint*, Mint[], Mfloat[]));
static Mint PROTO(l_s3rf,(Mfloat*, Mint*, Mint*, Mint*, Mint*));
static void PROTO(l_s4rf,(Mint*, Mfloat*, Mfloat[], Mfloat[]));
static void PROTO(l_s5rf,(Mfloat*, Mfloat [], Mint*, Mint[],
                   Mint*, Mint[], Mfloat[], Mint*,
	           Mfloat*, Mfloat*, Mfloat*));
static void PROTO(l_s6rf,(Mint*, Mfloat*, Mint*, Mint[], Mint*,
                   Mint[], Mint[], Mint[], Mfloat[]));
static void PROTO(l_s7rf,(Mfloat*, Mint*, Mint[], Mint*, Mint[],
                   Mint*, Mint*, Mfloat[], Mfloat[],
                   Mint[], Mint[]));

static VA_LIST_HACK  PROTO(l_scattered_2d_interp,(Mint ndata, Mfloat *xydata,
                                             Mfloat fdata[], Mint nx_out, Mint ny_out,
                                             Mfloat x_out[], Mfloat y_out[],
                                             va_list argptr));
                    
static Mfloat           *lv_result = NULL;
#ifdef ANSI
Mfloat  *imsl_f_scattered_2d_interp(Mint ndata, Mfloat xydata[],Mfloat fdata[], 
                                    Mint nx_out, Mint ny_out, Mfloat x_out[],
                                    Mfloat y_out[], ...)
#else
Mfloat  *imsl_f_scattered_2d_interp( ndata, xydata, fdata, nx_out, ny_out,
                                     x_out, y_out, va_alist)
   Mint          ndata;
   Mfloat        *xydata;
   Mfloat        fdata[];
   Mint          nx_out;
   Mint          ny_out;
   Mfloat        x_out[];
   Mfloat        y_out[];
   va_dcl
#endif
{
    va_list     argptr;
    VA_START(argptr,y_out);
#ifdef DOUBLE
   imsl_e1psh("imsl_d_scattered_2d_interp");
#else
   imsl_e1psh("imsl_f_scattered_2d_interp");
#endif
    IMSL_CALL(l_scattered_2d_interp( ndata, xydata, fdata, nx_out, ny_out, x_out, y_out, argptr));
    va_end(argptr);
#ifdef DOUBLE
   imsl_e1pop("imsl_d_scattered_2d_interp");
#else
   imsl_e1pop("imsl_f_scattered_2d_interp");
#endif
    return lv_result;
}

#ifdef ANSI
static VA_LIST_HACK   l_scattered_2d_interp(Mint ndata, Mfloat *xydata,Mfloat fdata[], Mint nx_out, Mint ny_out, Mfloat x_out[], Mfloat y_out[], va_list argptr)
#else
static VA_LIST_HACK   l_scattered_2d_interp( ndata, xydata, fdata, nx_out, ny_out, x_out, y_out, argptr)
   Mint          ndata;
   Mfloat        *xydata;
   Mfloat        fdata[];
   Mint          nx_out;
   Mint          ny_out;
   Mfloat        x_out[];
   Mfloat        y_out[];
   va_list       argptr;
#endif

{
   Mint          arg_number = 7;
   Mint          code;
   Mint          users_space = 0;
   Mint          free_result = 0;
   Mint          sur_col_dim;
   Mint          sur_col_dim_given = 0;
   Mfloat        *wk = NULL;
   Mint          *iwk = NULL;
   Mfloat        *result = NULL;
    code = 1;
    while (code > 0) {
        code = va_arg(argptr, int);
        arg_number++;
        switch (code) {
            case IMSL_RETURN_USER:
                result =  va_arg(argptr, Mfloat*);
                arg_number++;
                users_space = 1;
                break;
            case IMSL_SUR_COL_DIM:
                sur_col_dim = va_arg(argptr, Mint);
                arg_number++;
                sur_col_dim_given = 1;
                break;
            case 0:
                break;
            default:
                imsl_e1sti (1, code);
                imsl_e1sti (2, arg_number);
                imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
                goto RETURN;
                break;
        }
    } 
  

  if (sur_col_dim_given == 0) sur_col_dim = ny_out;

	if (ndata <= 3) {
		imsl_e1sti(1, ndata);
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_AT_LEAST_4_PTS); 
	}
	/* Check NXOUT. */
	if (nx_out <= 0) {
		imsl_e1sti(1, nx_out);
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_NXOUT_GT_ZERO);
	}
	/* Check NYOUT. */
	if (ny_out <= 0) {
		imsl_e1sti(1, ny_out);
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_NYOUT_GT_ZERO);
	}
  /* Check SUR_COL_DIM */
        if (sur_col_dim < ny_out) {
                imsl_e1sti(1, sur_col_dim);
                imsl_e1sti(2, ny_out);
                imsl_ermes(IMSL_TERMINAL, IMSL_COL_DIM_SUR);
        }
        if (imsl_n1rty(0) > 0) goto FREE_SPACE;
            
               /* GET WORKSPACE NEEDED IN S2RF */
 
     wk = (Mfloat *) imsl_malloc((6*ndata)*sizeof(*wk));
     iwk = (Mint *) imsl_malloc((31*ndata+nx_out*ny_out)*sizeof(*iwk));
        if ((wk == NULL)||(iwk == NULL)) {
            imsl_e1stl(1, "nx_out");
            imsl_e1sti(1, nx_out);
            imsl_e1stl(2, "ny_out");
            imsl_e1sti(2, ny_out);
            imsl_e1stl(3, "ndata");
            imsl_e1sti(3, ndata);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_3);
            goto FREE_SPACE;
        }
               /* GET SPACE NEED FOR RESULT IF NEEDED */
  if (users_space == 0){
      result = (Mfloat *) imsl_malloc(nx_out*sur_col_dim*sizeof(*result));
        if (result == NULL) {
            imsl_e1stl(1, "nx_out");
            imsl_e1sti(1, nx_out);
            imsl_e1stl(2, "sur_col_dim");
            imsl_e1sti(2, sur_col_dim);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
            goto FREE_SPACE;
        }
  }
               /* CALL S2RF */
     l_s2rf(&ndata, xydata, fdata, &nx_out, &ny_out, x_out, y_out,
            result, &nx_out, iwk, wk);
     if (imsl_n1rty(1) > 3) {
	 free_result = 1;
         goto FREE_SPACE;
     }
               /* TRANSPOSE THE RESULT SO IT WILL BE IN C STORAGE */ 
    imsl_f_m1ran(sur_col_dim, nx_out, result, result);
    lv_result = result;
               /* FREE THE SPACE */
FREE_SPACE:
   if ((free_result == 1)&&(result != NULL)&&(users_space == 0)) imsl_free(result);
   if (wk  != NULL)                          imsl_free(wk);
   if (iwk != NULL)                          imsl_free(iwk);
RETURN:
    return argptr;
}




















/* Structured by FOR_STRUCT, v0.2, on 08/23/90 at 14:29:06
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  S2RF/DS2RF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 12, 1986

    Purpose:    Compute a smooth bivariate interpolant to scattered data
                which is locally a quintic polynomial in two variables.

    Usage:      CALL S2RF (NDATA, XYDATA, FDATA, NXOUT, NYOUT,
                           XOUT, YOUT, SUR, LDSUR, IWK, WK)

    Arguments:  (See SURF)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
static struct t_s8rf {
	Mfloat           dmmy[27];
	Mint             itpv;
}               lv_s8rf;

#ifdef ANSI
static void l_s2rf(Mint *ndata, Mfloat *xydata, Mfloat fdata[], Mint *nxout,
                   Mint *nyout, Mfloat xout[], Mfloat yout[], Mfloat *sur,
                   Mint *ldsur, Mint iwk[], Mfloat wk[])
#else
static void l_s2rf(ndata, xydata, fdata, nxout, nyout, xout, yout,
	  sur, ldsur, iwk, wk)
	Mint            *ndata;
	Mfloat          * xydata, fdata[];
	Mint            *nxout, *nyout;
	Mfloat           xout[], yout[], *sur;
	Mint            *ldsur, iwk[];
	Mfloat           wk[];
#endif
{
#define XYDATA(I_,J_)   (xydata+(I_)*(2)+(J_))
#define SUR(I_,J_)	(sur+(I_)*(*ldsur)+(J_))
	Mint             i, il1, il2, iti, ixi, iyi, iz, jig0mn, jig0mx,
	                jig1mn, jig1mx, jigp, jngp, jwigp, jwigp0, jwipl,
	                jwipt, jwiwl, jwiwp, jwngp, jwngp0, ndp0, ngp0,
	                ngp1, nl, nngp, nt, nxi0, nyi0;


	imsl_e1psh("l_s2rf");
	/*
	 * Setting of some input parameters to local variables.
	 */
	ndp0 = *ndata;
	nxi0 = *nxout;
	nyi0 = *nyout;
	/* Check NDATA */
	if (*ndata <= 3) {
		imsl_e1sti(1, *ndata);

/*		imsl_ermes(5, 1, "The number of data points must be at least 4 while NDATA = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_AT_LEAST_4_PTS); 
	}
	/* Check NXOUT. */
	if (*nxout <= 0) {
		imsl_e1sti(1, *nxout);

/*		imsl_ermes(5, 2, "The number of elements in XOUT must be greater than zero while NXOUT = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_NXOUT_GT_ZERO);
	}
	/* Check NYOUT. */
	if (*nyout <= 0) {
		imsl_e1sti(1, *nyout);

/*		imsl_ermes(5, 3, "The number of elements in YOUT must be greater than zero while NYOUT = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_NYOUT_GT_ZERO);
	}
	/* Check LDSUR */
	if (*ldsur < *nxout) {
		imsl_e1sti(1, *nxout);
		imsl_e1sti(2, *ldsur);

/*		imsl_ermes(5, 4, "The leading dimension of SUR must be at least as big as the number of elements in XOUT, while LDSUR = %(i2) and NXOUT = %(i1) are given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_LDSUR_TOO_SMALL);
	}
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* Check XOUT */
	for (i = 1; i <= (*nxout - 1); i++) {
		if (xout[i] <= xout[i - 1]) {
			imsl_e1sti(1, i-1);
			imsl_e1sti(2, i );
			imsl_e1str(1, xout[i - 1]);
			imsl_e1str(2, xout[i]);

/*			imsl_ermes(4, 6, "The elements in XOUT must be increasing, while XOUT(%(i1)) = %(r1) and XOUT(%(i2)) = %(r2) are given.");
*/
                        imsl_ermes(IMSL_FATAL,
			IMSL_XOUT_NOT_STRICTLY_INCRSING);
			goto L_9000;
		}
	}
	/* Check YOUT */
	for (i = 1; i <= (*nyout - 1); i++) {
		if (yout[i] <= yout[i - 1]) {
			imsl_e1sti(1, i-1);
			imsl_e1sti(2, i );
			imsl_e1str(1, yout[i - 1]);
			imsl_e1str(2, yout[i]);

/*			imsl_ermes(4, 7, "The elements in YOUT must be increasing, while YOUT(%(i1)) = %(r1) and YOUT(%(i2)) = %(r2) are given.");
*/
                        imsl_ermes(IMSL_FATAL,
			IMSL_YOUT_NOT_STRICTLY_INCRSING);
			goto L_9000;
		}
	}

	iwk[0] = ndp0;
	iwk[2] = nxi0;
	iwk[3] = nyi0;
	/* Partition storage in IWK. */
	jwipt = 16;
	jwiwl = 6 * ndp0 + 1;
	jwngp0 = jwiwl - 1;
	jwipl = 24 * ndp0 + 1;
	jwiwp = 30 * ndp0 + 1;
	jwigp0 = 31 * ndp0;
	/* Triangulates the X-Y plane. */
	l_s6rf(&ndp0, xydata, &nt, &iwk[jwipt - 1], &nl, &iwk[jwipl - 1],
		  &iwk[jwiwl - 1], &iwk[jwiwp - 1], wk);
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	iwk[4] = nt;
	iwk[5] = nl;
	if (nt != 0) {
		/*
		 * Sorts output grid points in ascending order of the
		 * triangle number and the border line segment number.
		 */
		l_s7rf(xydata, &nt, &iwk[jwipt - 1], &nl, &iwk[jwipl - 1],
		      &nxi0, &nyi0, xout, yout, &iwk[jwngp0], &iwk[jwigp0]);

		/*
		 * Estimates partial derivatives at all data points.
		 */
		l_s4rf(&ndp0, xydata, fdata, wk);

		/* Interpolates the ZI values. */
		lv_s8rf.itpv = 0;
		jig0mx = 0;
		jig1mn = nxi0 * nyi0 + 1;
		nngp = nt + 2 * nl;
		for (jngp = 1; jngp <= nngp; jngp++) {
			iti = jngp;
			if (jngp > nt) {
				il1 = (jngp - nt + 1) / 2;
				il2 = (jngp - nt + 2) / 2;
				if (il2 > nl)
					il2 = 1;
				iti = il1 * (nt + nl) + il2;
			}
			jwngp = jwngp0 + jngp;
			ngp0 = iwk[jwngp - 1];
			if (ngp0 != 0) {
				jig0mn = jig0mx + 1;
				jig0mx += ngp0;
				for (jigp = jig0mn; jigp <= jig0mx; jigp++) {
					jwigp = jwigp0 + jigp;
					iz = iwk[jwigp - 1];
					iyi = (iz - 1) / nxi0 + 1;
					ixi = iz - nxi0 * (iyi - 1);
					l_s5rf(xydata, fdata, &nt, &iwk[jwipt - 1], &nl,
						  &iwk[jwipl - 1], wk, &iti, &xout[ixi - 1], &yout[iyi - 1],
						  SUR(iyi - 1, ixi - 1));
				}
			}
			jwngp = jwngp0 + 2 * nngp + 1 - jngp;
			ngp1 = iwk[jwngp - 1];
			if (ngp1 != 0) {
				jig1mx = jig1mn - 1;
				jig1mn -= ngp1;
				for (jigp = jig1mn; jigp <= jig1mx; jigp++) {
					jwigp = jwigp0 + jigp;
					iz = iwk[jwigp - 1];
					iyi = (iz - 1) / nxi0 + 1;
					ixi = iz - nxi0 * (iyi - 1);
					l_s5rf(xydata, fdata, &nt, &iwk[jwipt - 1], &nl,
						  &iwk[jwipl - 1], wk, &iti, &xout[ixi - 1], &yout[iyi - 1],
						  SUR(iyi - 1, ixi - 1));
				}
			}
		}

	}
L_9000:
	imsl_e1pop("l_s2rf");
	return;
}				/* end of function */
#undef SUR
/*----------------------------------------------------------------------- */

/*  IMSL Name:  S3RF/DS3RF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 12, 1986

    Purpose:    Compute a smooth bivariate interpolant to scattered data
                which is locally a quintic polynomial in two variables.

    Usage:      S3RF(XYDATA, I1, I2, I3, I4)

    Arguments:
       XYDATA - A 2 by NDATA Array containing the coordinates of
                the interpolation points.  (Input)
                These points must be distinct.  The x-coordinate
                of the Ith data point is stored in XYDATA(1,I) and the
                y-coordinate of the Ith data point is stored in
                XYDATA(2,I).
       I1     - A pointer to XYDATA.  (Input)
       I2     - A pointer to XYDATA.  (Input)
       I3     - A pointer to XYDATA.  (Input)
       I4     - A pointer to XYDATA.  (Input)
       S3RF   - Integer function value.  (Output)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static Mint l_s3rf(Mfloat *xydata, Mint *i1, Mint *i2, Mint *i3, Mint *i4)
#else
static Mint l_s3rf(xydata, i1, i2, i3, i4)
        Mfloat          * xydata;
        Mint            *i1, *i2, *i3, *i4;
#endif
{
#define XYDATA(I_,J_)   (xydata+(I_)*(2)+(J_))
        Mint             iiidx, s3rf_v;
        Mfloat          a1sq, a2sq, a3sq, a4sq, c1sq, c3sq, epsln, s1sq,
                        s2sq, s3sq, s4sq, u1, u2, u3, u4, x1, x2, x3, x4,
                        y1, y2, y3, y4;



        epsln = 100.0 * imsl_amach(4);
        /* PRELIMINARY PROCESSING */
        x1 = *XYDATA(*i1 - 1, 0);
        y1 = *XYDATA(*i1 - 1, 1);
        x2 = *XYDATA(*i2 - 1, 0);
        y2 = *XYDATA(*i2 - 1, 1);
        x3 = *XYDATA(*i3 - 1, 0);
        y3 = *XYDATA(*i3 - 1, 1);
        x4 = *XYDATA(*i4 - 1, 0);
        y4 = *XYDATA(*i4 - 1, 1);
        /* CALCULATION */
        iiidx = 0;
        u3 = (y2 - y3) * (x1 - x3) - (x2 - x3) * (y1 - y3);
        u4 = (y1 - y4) * (x2 - x4) - (x1 - x4) * (y2 - y4);
        if (u3 * u4 > F_ZERO) {
                u1 = (y3 - y1) * (x4 - x1) - (x3 - x1) * (y4 - y1);
                u2 = (y4 - y2) * (x3 - x2) - (x4 - x2) * (y3 - y2);
                a1sq = imsl_fi_power(x1 - x3, 2) + imsl_fi_power(y1 - y3, 2);
                a4sq = imsl_fi_power(x4 - x1, 2) + imsl_fi_power(y4 - y1, 2);
                c1sq = imsl_fi_power(x3 - x4, 2) + imsl_fi_power(y3 - y4, 2);
                a2sq = imsl_fi_power(x2 - x4, 2) + imsl_fi_power(y2 - y4, 2);
                a3sq = imsl_fi_power(x3 - x2, 2) + imsl_fi_power(y3 - y2, 2);
                c3sq = imsl_fi_power(x2 - x1, 2) + imsl_fi_power(y2 - y1, 2);
                s1sq = u1 * u1 / (c1sq * imsl_f_max(a1sq, a4sq));
                s2sq = u2 * u2 / (c1sq * imsl_f_max(a2sq, a3sq));
                s3sq = u3 * u3 / (c3sq * imsl_f_max(a3sq, a1sq));
                s4sq = u4 * u4 / (c3sq * imsl_f_max(a4sq, a2sq));
                if ((imsl_f_min(s3sq, s4sq) - imsl_f_min(s1sq, s2sq)) > epsln)
                        iiidx = 1;
        }
        s3rf_v = iiidx;
        return (s3rf_v);
}                               /* end of function */


/*----------------------------------------------------------------------- */

/*  IMSL Name:  S5RF/DS5RF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 12, 1986

    Purpose:    Compute a smooth bivariate interpolant to scattered data
                which is locally a quintic polynomial in two variables.

    Usage:      CALL S5RF (XYDATA, FDATA, NT, IPT, NL, IPL, PDD, ITI,
                           XII, YII, ZII)

    Arguments:
       XYDATA - A 2 by NDATA Array containing the coordinates of
                the interpolation points.  (Input)
                These points must be distinct.  The x-coordinate
                of the Ith data point is stored in XYDATA(1,I) and the
                y-coordinate of the Ith data point is stored in
                XYDATA(2,I).
       FDATA  - Array of length NDATA containing the interpolation
                values.  (Input)
       NT     - A counter.  (Input)
       IPT    - Work array of length 6*NDATA-14.
       NL     - A counter.  (Input)
       IPL    - Work array of length 6*NDATA.
       PDD    - Work array of length 6*NDATA.
       ITI    - A counter.  (Input)
       XII    - Value of x-coordinate where interpolant is to be
                evaluated.  (Input)
       YII    - Value of y-coordinate where interpolant is to be
                evaluated.  (Input)
       ZII    - Value of of interpolant at (XII,YII).  (Output)

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_s5rf(Mfloat *xydata, Mfloat fdata[], Mint *nt, Mint ipt[],
                   Mint *nl, Mint ipl[], Mfloat pdd[], Mint *iti,
	           Mfloat *xii, Mfloat *yii, Mfloat *zii)
#else
static void l_s5rf(xydata, fdata, nt, ipt, nl, ipl, pdd, iti,
        	  xii, yii, zii)
	Mfloat          *xydata, fdata[];
	Mint            *nt, ipt[], *nl, ipl[];
	Mfloat           pdd[];
	Mint            *iti;
	Mfloat          *xii, *yii, *zii;
#endif
{
	Mint             i, idp, il1, il2, it0, jipl, jipt, jpd, jpdd, kpd,
	                ntl;
	Mfloat           a, aa, ab, ac, act2, ad, adbc, b, bb, bc, bdt2,
	                c, cc, cd, csuv, d, dd, dlt, dx, dy, g1, g2, h1,
	                h2, h3, lu, lv, p0, p1, p2, p3, p4, *p5, pd[15],
	                thsv, thus, thuv, thxu, u, v, x[3], y[3], z[3], z0,
	                zu[3], zuu[3], zuv[3], zv[3], zvv[3];
	struct {
		Mfloat           x0, y0, ap, bp, cp, dp, p00, p10, p20,
		                p30, p40, p50, p01, p11, p21, p31, p41,
		                p02, p12, p22, p32, p03, p13, p23, p04,
		                p14, p05;
		Mint             itpv;
	}              *_s8rf = (void *) &lv_s8rf;
	p5 = (Mfloat *) &_s8rf->p50;

	/* Preliminary processing */
	it0 = *iti;
	ntl = *nt + *nl;
	if (it0 > ntl) {
		il1 = it0 / ntl;
		il2 = it0 - il1 * ntl;
		if (il1 == il2)
			goto L_40;
		goto L_80;
		/*
		 * Calculation of ZII by interpolation. Checks if the
		 * necessary coefficients have been calculated.
		 */
	}
	if (it0 != _s8rf->itpv) {
		/*
		 * Loads coordinates and partial derivative values at the
		 * vertexes.
		 */
		jipt = 3 * (it0 - 1);
		jpd = 0;
		for (i = 1; i <= 3; i++) {
			jipt += 1;
			idp = ipt[jipt - 1];
			x[i - 1] = *XYDATA(idp - 1, 0);
			y[i - 1] = *XYDATA(idp - 1, 1);
			z[i - 1] = fdata[idp - 1];
			jpdd = 5 * (idp - 1);
			for (kpd = 1; kpd <= 5; kpd++) {
				jpd += 1;
				jpdd += 1;
				pd[jpd - 1] = pdd[jpdd - 1];
			}
		}
		/*
		 * Determines the coefficients for the coordinate system
		 * transformation from the X-Y system to the U-V system and
		 * vice-versa.
		 */
		_s8rf->x0 = x[0];
		_s8rf->y0 = y[0];
		a = x[1] - _s8rf->x0;
		b = x[2] - _s8rf->x0;
		c = y[1] - _s8rf->y0;
		d = y[2] - _s8rf->y0;
		ad = a * d;
		bc = b * c;
		dlt = ad - bc;
		_s8rf->ap = d / dlt;
		_s8rf->bp = -b / dlt;
		_s8rf->cp = -c / dlt;
		_s8rf->dp = a / dlt;
		/*
		 * Converts the partial derivatives at the vertexes of the
		 * triangle for the U-V coordinate system.
		 */
		aa = a * a;
		act2 = F_TWO * a * c;
		cc = c * c;
		ab = a * b;
		adbc = ad + bc;
		cd = c * d;
		bb = b * b;
		bdt2 = F_TWO * b * d;
		dd = d * d;
		for (i = 1; i <= 3; i++) {
			jpd = 5 * i;
			zu[i - 1] = a * pd[jpd - 5] + c * pd[jpd - 4];
			zv[i - 1] = b * pd[jpd - 5] + d * pd[jpd - 4];
			zuu[i - 1] = aa * pd[jpd - 3] + act2 * pd[jpd - 2] + cc * pd[jpd - 1];
			zuv[i - 1] = ab * pd[jpd - 3] + adbc * pd[jpd - 2] + cd * pd[jpd - 1];
			zvv[i - 1] = bb * pd[jpd - 3] + bdt2 * pd[jpd - 2] + dd * pd[jpd - 1];
		}
		/*
		 * Calculates the coefficients of the polynomial.
		 */
		_s8rf->p00 = z[0];
		_s8rf->p10 = zu[0];
		_s8rf->p01 = zv[0];
		_s8rf->p20 = F_HALF * zuu[0];
		_s8rf->p11 = zuv[0];
		_s8rf->p02 = F_HALF * zvv[0];
		h1 = z[1] - _s8rf->p00 - _s8rf->p10 - _s8rf->p20;
		h2 = zu[1] - _s8rf->p10 - zuu[0];
		h3 = zuu[1] - zuu[0];
		_s8rf->p30 = F_TEN * h1 - F_FOUR * h2 + F_HALF * h3;
		_s8rf->p40 = -15.0 * h1 + F_SEVEN * h2 - h3;
		_s8rf->p50 = F_SIX * h1 - F_THREE * h2 + F_HALF * h3;
		h1 = z[2] - _s8rf->p00 - _s8rf->p01 - _s8rf->p02;
		h2 = zv[2] - _s8rf->p01 - zvv[0];
		h3 = zvv[2] - zvv[0];
		_s8rf->p03 = F_TEN * h1 - F_FOUR * h2 + F_HALF * h3;
		_s8rf->p04 = -15.0 * h1 + F_SEVEN * h2 - h3;
		_s8rf->p05 = F_SIX * h1 - F_THREE * h2 + F_HALF * h3;
		lu = sqrt(aa + cc);
		lv = sqrt(bb + dd);
		thxu = atan2(c, a);
		thuv = atan2(d, b) - thxu;
		csuv = cos(thuv);
		_s8rf->p41 = F_FIVE * lv * csuv / lu * _s8rf->p50;
		_s8rf->p14 = F_FIVE * lu * csuv / lv * _s8rf->p05;
		h1 = zv[1] - _s8rf->p01 - _s8rf->p11 - _s8rf->p41;
		h2 = zuv[1] - _s8rf->p11 - F_FOUR * _s8rf->p41;
		_s8rf->p21 = F_THREE * h1 - h2;
		_s8rf->p31 = -F_TWO * h1 + h2;
		h1 = zu[2] - _s8rf->p10 - _s8rf->p11 - _s8rf->p14;
		h2 = zuv[2] - _s8rf->p11 - F_FOUR * _s8rf->p14;
		_s8rf->p12 = F_THREE * h1 - h2;
		_s8rf->p13 = -F_TWO * h1 + h2;
		thus = atan2(d - c, b - a) - thxu;
		thsv = thuv - thus;
		aa = sin(thsv) / lu;
		bb = -cos(thsv) / lu;
		cc = sin(thus) / lv;
		dd = cos(thus) / lv;
		ac = aa * cc;
		ad = aa * dd;
		bc = bb * cc;
		g1 = aa * ac * (F_THREE * bc + F_TWO * ad);
		g2 = cc * ac * (F_THREE * ad + F_TWO * bc);
		h1 = -aa * aa * aa * (F_FIVE * aa * bb * _s8rf->p50 + (F_FOUR * bc + ad) * _s8rf->p41) -
			cc * cc * cc * (F_FIVE * cc * dd * _s8rf->p05 + (F_FOUR * ad + bc) * _s8rf->p14);
		h2 = F_HALF * zvv[1] - _s8rf->p02 - _s8rf->p12;
		h3 = F_HALF * zuu[2] - _s8rf->p20 - _s8rf->p21;
		_s8rf->p22 = (g1 * h2 + g2 * h3 - h1) / (g1 + g2);
		_s8rf->p32 = h2 - _s8rf->p22;
		_s8rf->p23 = h3 - _s8rf->p22;
		_s8rf->itpv = it0;
		/* Converts XII and YII to U-V system. */
	}
	dx = *xii - _s8rf->x0;
	dy = *yii - _s8rf->y0;
	u = _s8rf->ap * dx + _s8rf->bp * dy;
	v = _s8rf->cp * dx + _s8rf->dp * dy;
	/* Evaluates the polynomial. */
	p0 = _s8rf->p00 + v * (_s8rf->p01 + v * (_s8rf->p02 + v * (_s8rf->p03 +
				       v * (_s8rf->p04 + v * _s8rf->p05))));
	p1 = _s8rf->p10 + v * (_s8rf->p11 + v * (_s8rf->p12 + v * (_s8rf->p13 +
							  v * _s8rf->p14)));
	p2 = _s8rf->p20 + v * (_s8rf->p21 + v * (_s8rf->p22 + v * _s8rf->p23));
	p3 = _s8rf->p30 + v * (_s8rf->p31 + v * _s8rf->p32);
	p4 = _s8rf->p40 + v * _s8rf->p41;
	*zii = p0 + u * (p1 + u * (p2 + u * (p3 + u * (p4 + u ** p5))));
	goto L_9000;
	/*
	 * Calculation of ZII by extrapolation in the rectangle. Checks if
	 * the necessary coefficients have been calculated.
	 */
L_40:
	;
	if (it0 != _s8rf->itpv) {
		/*
		 * Loads coordinate and partial derivative values at the end
		 * points of the border line segment.
		 */
		jipl = 3 * (il1 - 1);
		jpd = 0;
		for (i = 1; i <= 2; i++) {
			jipl += 1;
			idp = ipl[jipl - 1];
			x[i - 1] = *XYDATA(idp - 1, 0);
			y[i - 1] = *XYDATA(idp - 1, 1);
			z[i - 1] = fdata[idp - 1];
			jpdd = 5 * (idp - 1);
			for (kpd = 1; kpd <= 5; kpd++) {
				jpd += 1;
				jpdd += 1;
				pd[jpd - 1] = pdd[jpdd - 1];
			}
		}
		/*
		 * Determines the coefficients for the coordinate system
		 * transformation from the X-Y system to the U-V system and
		 * vice-versa.
		 */
		_s8rf->x0 = x[0];
		_s8rf->y0 = y[0];
		a = y[1] - y[0];
		b = x[1] - x[0];
		c = -b;
		d = a;
		ad = a * d;
		bc = b * c;
		dlt = ad - bc;
		_s8rf->ap = d / dlt;
		_s8rf->bp = -b / dlt;
		_s8rf->cp = -_s8rf->bp;
		_s8rf->dp = _s8rf->ap;
		/*
		 * Converts the partial derivatives at the end points of the
		 * border line segemnt for the U-V coordinate system.
		 */
		aa = a * a;
		act2 = F_TWO * a * c;
		cc = c * c;
		ab = a * b;
		adbc = ad + bc;
		cd = c * d;
		bb = b * b;
		bdt2 = F_TWO * b * d;
		dd = d * d;
		for (i = 1; i <= 2; i++) {
			jpd = 5 * i;
			zu[i - 1] = a * pd[jpd - 5] + c * pd[jpd - 4];
			zv[i - 1] = b * pd[jpd - 5] + d * pd[jpd - 4];
			zuu[i - 1] = aa * pd[jpd - 3] + act2 * pd[jpd - 2] + cc * pd[jpd - 1];
			zuv[i - 1] = ab * pd[jpd - 3] + adbc * pd[jpd - 2] + cd * pd[jpd - 1];
			zvv[i - 1] = bb * pd[jpd - 3] + bdt2 * pd[jpd - 2] + dd * pd[jpd - 1];
		}
		/*
		 * Calcualtes the coefficients of the polynomial.
		 */
		_s8rf->p00 = z[0];
		_s8rf->p10 = zu[0];
		_s8rf->p01 = zv[0];
		_s8rf->p20 = F_HALF * zuu[0];
		_s8rf->p11 = zuv[0];
		_s8rf->p02 = F_HALF * zvv[0];
		h1 = z[1] - _s8rf->p00 - _s8rf->p01 - _s8rf->p02;
		h2 = zv[1] - _s8rf->p01 - zvv[0];
		h3 = zvv[1] - zvv[0];
		_s8rf->p03 = F_TEN * h1 - F_FOUR * h2 + F_HALF * h3;
		_s8rf->p04 = -15.0 * h1 + F_SEVEN * h2 - h3;
		_s8rf->p05 = F_SIX * h1 - F_THREE * h2 + F_HALF * h3;
		h1 = zu[1] - _s8rf->p10 - _s8rf->p11;
		h2 = zuv[1] - _s8rf->p11;
		_s8rf->p12 = F_THREE * h1 - h2;
		_s8rf->p13 = -F_TWO * h1 + h2;
		_s8rf->p21 = F_ZERO;
		_s8rf->p23 = -zuu[1] + zuu[0];
		_s8rf->p22 = -1.5 * _s8rf->p23;
		_s8rf->itpv = it0;
		/* Converts XII and YII to U-V system. */
	}
	dx = *xii - _s8rf->x0;
	dy = *yii - _s8rf->y0;
	u = _s8rf->ap * dx + _s8rf->bp * dy;
	v = _s8rf->cp * dx + _s8rf->dp * dy;
	/* Evaulates the polynomial. */
	p0 = _s8rf->p00 + v * (_s8rf->p01 + v * (_s8rf->p02 + v * (_s8rf->p03 +
				       v * (_s8rf->p04 + v * _s8rf->p05))));
	p1 = _s8rf->p10 + v * (_s8rf->p11 + v * (_s8rf->p12 + v * _s8rf->p13));
	p2 = _s8rf->p20 + v * (_s8rf->p21 + v * (_s8rf->p22 + v * _s8rf->p23));
	*zii = p0 + u * (p1 + u * p2);
	goto L_9000;
	/*
	 * Calculation of ZII by extrapolation in the triangle. Checks if the
	 * necessary coefficients have been calculated.
	 */
L_80:
	;
	if (it0 != _s8rf->itpv) {
		/*
		 * Loads coordinate and partial derivative values at the
		 * vertex of the triangle.
		 */
		jipl = 3 * il2 - 2;
		idp = ipl[jipl - 1];
		_s8rf->x0 = *XYDATA(idp - 1, 0);
		_s8rf->y0 = *XYDATA(idp - 1, 1);
		z0 = fdata[idp - 1];
		jpdd = 5 * (idp - 1);
		for (kpd = 1; kpd <= 5; kpd++) {
			jpdd += 1;
			pd[kpd - 1] = pdd[jpdd - 1];
		}
		/*
		 * Calculates the coefficients of the polynomial.
		 */
		_s8rf->p00 = z0;
		_s8rf->p10 = pd[0];
		_s8rf->p01 = pd[1];
		_s8rf->p20 = F_HALF * pd[2];
		_s8rf->p11 = pd[3];
		_s8rf->p02 = F_HALF * pd[4];
		_s8rf->itpv = it0;
		/* Converts XII and YII to U-V system. */
	}
	u = *xii - _s8rf->x0;
	v = *yii - _s8rf->y0;
	/* Evaluates the polynomial. */
	p0 = _s8rf->p00 + v * (_s8rf->p01 + v * _s8rf->p02);
	p1 = _s8rf->p10 + v * _s8rf->p11;
	*zii = p0 + u * (p1 + u * _s8rf->p20);
L_9000:
	return;
}				/* end of function */
/*----------------------------------------------------------------------- */

/*  IMSL Name:  S6RF/DS6RF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 12, 1986

    Purpose:    Compute a smooth bivariate interpolant to scattered data
                which is locally a quintic polynomial in two variables.

    Usage:      CALL S6RF (NDATA, XYDATA, NT, IPT, NL, IPL, IWL, IWP,
                           WK)

    Arguments:
       NDATA  - Number of data points.  (Input)
                NDATA must be at least 3.
       XYDATA - A 2 by NDATA Array containing the coordinates of
                the interpolation points.  (Input)
                These points must be distinct.  The x-coordinate
                of the Ith data point is stored in XYDATA(1,I) and the
                y-coordinate of the Ith data point is stored in
                XYDATA(2,I).
       NT     - A counter.  (Output)
       IPT    - Work array of length 6*NDATA-14.
       NL     - A counter.  (Output)
       IPL    - Work array of length 6*NDATA.
       IWL    - Work array of length 18*NATA.
       IWP    - Work array of length NDATA.
       WK     - Work array of length 6*NDATA.

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_s6rf(Mint *ndata, Mfloat *xydata, Mint *nt, Mint ipt[], Mint *nl,
                   Mint ipl[], Mint iwl[], Mint iwp[], Mfloat wk[])
#else
static void l_s6rf(ndata, xydata, nt, ipt, nl, ipl, iwl, iwp,
        	  wk)
	Mint            *ndata;
	Mfloat          *xydata;
	Mint            *nt, ipt[], *nl, ipl[], iwl[], iwp[];
	Mfloat           wk[];
#endif
{
#define XYDATA(I_,J_)   (xydata+(I_)*(2)+(J_))
	Mint             i, il, ilf, iliv, ilt3, ilvs, ip, ip1, ip2, ip3,
	                ipl1, ipl2, iplj1, iplj2, ipmn1, ipmn2, ipt1, ipt2,
	                ipt3, ipti, ipti1, ipti2, irep, it, it1t3, it2t3,
	                itf[2], its, itt3, itt3r, ixvs, ixvspv, j, jl1,
	                jl2, jlt3, jmin, jp, jp1, jpmx, jwl, jwl1, jwl1mn,
	                ndp0, nl0, nlf, nlfc, nlft2, nln, nlnt3,
	                nlsh, nlsht3, nlt3, nt0, ntf, ntt3, ntt3p3;
	Mfloat           dsqi, dsqmn, epsln, sp,
	                vp, x1, x2, x3, xdmp, y1, y2, y3, ydmp;
	static Mint      nrep = 100;


	/* STATEMENT FUNCTIONS */
#define DSQF(u1,v1,u2,v2)	(Mfloat)(imsl_fi_power((u2) - (u1), 2) + imsl_fi_power((v2) - \
	 (v1), 2))
#define SPDT(u1,v1,u2,v2,u3,v3)	(Mfloat)(((u2) - (u1))*((u3) - \
	 (u1)) + ((v2) - (v1))*((v3) - (v1)))
#define VPDT(u1,v1,u2,v2,u3,v3)	(Mfloat)(((v3) - (v1))*((u2) - \
	 (u1)) - ((u3) - (u1))*((v2) - (v1)))

	imsl_e1psh("l_s6rf");

	epsln = 100.0e0 * imsl_amach(4);
	/* Preliminary processing */
	ndp0 = *ndata;
	/*
	 * Determines the closest pair of data points and their midpoint.
	 */
	dsqmn = DSQF(*XYDATA(0, 0), *XYDATA(0, 1), *XYDATA(1, 0), *XYDATA(1, 1));
	ipmn1 = 1;
	ipmn2 = 2;
	for (i = 1; i <= (*ndata - 1); i++) {
		x1 = *XYDATA(i - 1, 0);
		y1 = *XYDATA(i - 1, 1);
		for (ip1 = i + 1; ip1 <= ndp0; ip1++) {
			dsqi = DSQF(x1, y1, *XYDATA(ip1 - 1, 0), *XYDATA(ip1 - 1, 1));
			if (dsqi == F_ZERO) {
				imsl_e1sti(1, i-1);
				imsl_e1sti(2, ip1-1);
				imsl_e1str(1, x1);
				imsl_e1str(2, y1);
				imsl_e1str(3, *XYDATA(ip1 - 1, 0));
				imsl_e1str(4, *XYDATA(ip1 - 1, 1));

/*				imsl_ermes(4, 5, "Two data points are equal.  XYDATA(%(i1),I) = (%(r1),%(r2)) and XYDATA(%(i2),I) = (%(r3),%(r4)).  These points must be distinct.");
*/
                                imsl_ermes(IMSL_FATAL,
				IMSL_DUPLICATE_XYDATA_VALUES);
				goto L_9000;
			} else if (dsqi < dsqmn) {
				dsqmn = dsqi;
				ipmn1 = i;
				ipmn2 = ip1;
			}
		}
	}
	xdmp = (*XYDATA(ipmn1 - 1, 0) + *XYDATA(ipmn2 - 1, 0)) / F_TWO;
	ydmp = (*XYDATA(ipmn1 - 1, 1) + *XYDATA(ipmn2 - 1, 1)) / F_TWO;
	/*
	 * Sorts the other (NDP-2) data points in ascending order of distance
	 * from the midpoint and stores the sorted data point numbers in the
	 * IWP array.
	 */
	jp1 = 2;
	for (i = 1; i <= *ndata; i++) {
		if (i != ipmn1 && i != ipmn2) {
			jp1 += 1;
			iwp[jp1 - 1] = i;
			wk[jp1 - 1] = DSQF(xdmp, ydmp, *XYDATA(i - 1, 0), *XYDATA(i - 1, 1));
		}
	}
	for (i = 3; i <= (*ndata - 1); i++) {
		dsqmn = wk[i - 1];
		jmin = i;
		for (j = i; j <= *ndata; j++) {
			if (wk[j - 1] < dsqmn) {
				dsqmn = wk[j - 1];
				jmin = j;
			}
		}
		its = iwp[i - 1];
		iwp[i - 1] = iwp[jmin - 1];
		iwp[jmin - 1] = its;
		wk[jmin - 1] = wk[i - 1];
	}
	/*
	 * If necessary, modifies the ordering in such a way that the first
	 * three data points are not collinear.
	 */
	x1 = *XYDATA(ipmn1 - 1, 0);
	y1 = *XYDATA(ipmn1 - 1, 1);
	x2 = *XYDATA(ipmn2 - 1, 0);
	y2 = *XYDATA(ipmn2 - 1, 1);
	for (j = 3; j <= *ndata; j++) {
		jp = j;
		ip = iwp[jp - 1];
		sp = SPDT(*XYDATA(ip - 1, 0), *XYDATA(ip - 1, 1), x1, y1, x2,
			  y2);
		vp = VPDT(*XYDATA(ip - 1, 0), *XYDATA(ip - 1, 1), x1, y1, x2,
			  y2);
		if (fabs(vp) > (fabs(sp) * epsln))
			goto L_70;
	}


/*	imsl_ermes(5, 8, "All points are collinear."); */
        imsl_ermes(IMSL_TERMINAL, IMSL_ALL_POINTS_COLLINEAR);
	goto L_9000;

L_70:
	if (jp != 3) {
		jpmx = jp;
		for (j = 4; j <= jpmx; j++) {
			jp = jpmx + 4 - j;
			iwp[jp - 1] = iwp[jp - 2];
		}
		iwp[2] = ip;
	}
	/*
	 * Forms the first triangle. Stores point numbers of the vertexes of
	 * the triangle in the IPT array, and stores point numbers of the
	 * border line segments and the triangle number in the IPL array.
	 */
	ip1 = ipmn1;
	ip2 = ipmn2;
	ip3 = iwp[2];
	if (VPDT(*XYDATA(ip1 - 1, 0), *XYDATA(ip1 - 1, 1), *XYDATA(ip2 - 1, 0),
	       *XYDATA(ip2 - 1, 1), *XYDATA(ip3 - 1, 0), *XYDATA(ip3 - 1, 1)) <
	    F_ZERO) {
		ip1 = ipmn2;
		ip2 = ipmn1;
	}
	nt0 = 1;
	ntt3 = 3;
	ipt[0] = ip1;
	ipt[1] = ip2;
	ipt[2] = ip3;
	nl0 = 3;
	nlt3 = 9;
	ipl[0] = ip1;
	ipl[1] = ip2;
	ipl[2] = 1;
	ipl[3] = ip2;
	ipl[4] = ip3;
	ipl[5] = 1;
	ipl[6] = ip3;
	ipl[7] = ip1;
	ipl[8] = 1;
	/*
	 * Adds the remaining (NDP-3) data points, one by one.
	 */
	for (jp1 = 4; jp1 <= ndp0; jp1++) {
		ip1 = iwp[jp1 - 1];
		x1 = *XYDATA(ip1 - 1, 0);
		y1 = *XYDATA(ip1 - 1, 1);
		/*
		 * Determines the first invisible and visible border line
		 * segments, ILIV and ILVS.
		 */
		for (il = 1; il <= nl0; il++) {
			ip2 = ipl[il * 3 - 3];
			ip3 = ipl[il * 3 - 2];
			x2 = *XYDATA(ip2 - 1, 0);
			y2 = *XYDATA(ip2 - 1, 1);
			x3 = *XYDATA(ip3 - 1, 0);
			y3 = *XYDATA(ip3 - 1, 1);
			sp = SPDT(x1, y1, x2, y2, x3, y3);
			vp = VPDT(x1, y1, x2, y2, x3, y3);
			if (il == 1) {
				ixvs = 0;
				if (vp <= (fabs(sp) * (-epsln)))
					ixvs = 1;
				iliv = 1;
				ilvs = 1;
			} else {
				ixvspv = ixvs;
				if (vp <= (fabs(sp) * (-epsln))) {
					ixvs = 1;
					if (ixvspv != 1) {
						ilvs = il;
						if (iliv != 1)
							goto L_100;
					}
				} else {
					ixvs = 0;
					if (ixvspv != 0) {
						iliv = il;
						if (ilvs != 1)
							goto L_100;
					}
				}
			}
		}
		if (iliv == 1 && ilvs == 1)
			ilvs = nl0;
L_100:
		if (ilvs < iliv)
			ilvs += nl0;
		/*
		 * Shifts (Rotates) the IPL array to have the invisible
		 * border line segments contained in the first part of the
		 * IPL array.
		 */
		if (iliv != 1) {
			nlsh = iliv - 1;
			nlsht3 = nlsh * 3;
			for (jl1 = 1; jl1 <= nlsht3; jl1++) {
				jl2 = jl1 + nlt3;
				ipl[jl2 - 1] = ipl[jl1 - 1];
			}
			for (jl1 = 1; jl1 <= nlt3; jl1++) {
				jl2 = jl1 + nlsht3;
				ipl[jl1 - 1] = ipl[jl2 - 1];
			}
			ilvs -= nlsh;
		}
		/*
		 * Adds triangles to the IPT array, updates border line
		 * segments in the IPL array, and sets flags for the border
		 * line segments to be re- examined in the IWL array.
		 */
		jwl = 0;
		for (il = ilvs; il <= nl0; il++) {
			ilt3 = il * 3;
			ipl1 = ipl[ilt3 - 3];
			ipl2 = ipl[ilt3 - 2];
			it = ipl[ilt3 - 1];
			/*
			 * Adds a triangle to the IPT array.
			 */
			nt0 += 1;
			ntt3 += 3;
			ipt[ntt3 - 3] = ipl2;
			ipt[ntt3 - 2] = ipl1;
			ipt[ntt3 - 1] = ip1;
			/*
			 * Updates border line segments in the IPL array.
			 */
			if (il == ilvs) {
				ipl[ilt3 - 2] = ip1;
				ipl[ilt3 - 1] = nt0;
			}
			if (il == nl0) {
				nln = ilvs + 1;
				nlnt3 = nln * 3;
				ipl[nlnt3 - 3] = ip1;
				ipl[nlnt3 - 2] = ipl[0];
				ipl[nlnt3 - 1] = nt0;
			}
			/*
			 * Determines the vertex that does not lie on the
			 * border line segments.
			 */
			itt3 = it * 3;
			ipti = ipt[itt3 - 3];
			if (ipti == ipl1 || ipti == ipl2) {
				ipti = ipt[itt3 - 2];
				if (ipti == ipl1 || ipti == ipl2) {
					ipti = ipt[itt3 - 1];
				}
			}
			/*
			 * Checks if the exchange is necessary.
			 */
			if (l_s3rf(xydata, &ip1, &ipti, &ipl1, &ipl2) != 0) {
				/*
				 * Modifies the IPT array when necessary.
				 */
				ipt[itt3 - 3] = ipti;
				ipt[itt3 - 2] = ipl1;
				ipt[itt3 - 1] = ip1;
				ipt[ntt3 - 2] = ipti;
				if (il == ilvs)
					ipl[ilt3 - 1] = it;
				if (il == nl0 && ipl[2] == it)
					ipl[2] = nt0;

				/* Sets flags in the IWL array. */
				jwl += 4;
				iwl[jwl - 4] = ipl1;
				iwl[jwl - 3] = ipti;
				iwl[jwl - 2] = ipti;
				iwl[jwl - 1] = ipl2;
			}
		}
		nl0 = nln;
		nlt3 = nlnt3;
		nlf = jwl / 2;
		/* No improvement is necessary */
		if (nlf == 0)
			goto L_200;
		/* Improves triangulation. */
		ntt3p3 = ntt3 + 3;
		for (irep = 1; irep <= nrep; irep++) {
			for (ilf = 1; ilf <= nlf; ilf++) {
				ipl1 = iwl[ilf * 2 - 2];
				ipl2 = iwl[ilf * 2 - 1];
				/*
				 * Locates in the IPT array two triangles on
				 * both sides of the flagged line segment.
				 */
				ntf = 0;
				for (itt3r = 3; itt3r <= ntt3; itt3r += 3) {
					itt3 = ntt3p3 - itt3r;
					ipt1 = ipt[itt3 - 3];
					ipt2 = ipt[itt3 - 2];
					ipt3 = ipt[itt3 - 1];
					if ((ipl1 == ipt1 || ipl1 == ipt2) || ipl1 ==
					    ipt3) {
						if ((ipl2 == ipt1 || ipl2 == ipt2) || ipl2 ==
						    ipt3) {
							ntf += 1;
							itf[ntf - 1] = itt3 / 3;
							if (ntf == 2)
								goto L_150;
						}
					}
				}
				if (ntf < 2)
					goto L_170;
				/*
				 * Determines the vertexes of the triangles
				 * that do not lie on the line segment.
				 */
		L_150:
				it1t3 = itf[0] * 3;
				ipti1 = ipt[it1t3 - 3];
				if (ipti1 == ipl1 || ipti1 == ipl2) {
					ipti1 = ipt[it1t3 - 2];
					if (ipti1 == ipl1 || ipti1 == ipl2) {
						ipti1 = ipt[it1t3 - 1];
					}
				}
				it2t3 = itf[1] * 3;
				ipti2 = ipt[it2t3 - 3];
				if (ipti2 == ipl1 || ipti2 == ipl2) {
					ipti2 = ipt[it2t3 - 2];
					if (ipti2 == ipl1 || ipti2 == ipl2) {
						ipti2 = ipt[it2t3 - 1];
					}
				}
				/*
				 * Checks if the exchange is necessary.
				 */
				if (l_s3rf(xydata, &ipti1, &ipti2, &ipl1, &ipl2) !=
				    0) {
					/*
					 * Modifies the IPT array when
					 * necessary.
					 */
					ipt[it1t3 - 3] = ipti1;
					ipt[it1t3 - 2] = ipti2;
					ipt[it1t3 - 1] = ipl1;
					ipt[it2t3 - 3] = ipti2;
					ipt[it2t3 - 2] = ipti1;
					ipt[it2t3 - 1] = ipl2;
					/* Sets new flags. */
					jwl += 8;
					iwl[jwl - 8] = ipl1;
					iwl[jwl - 7] = ipti1;
					iwl[jwl - 6] = ipti1;
					iwl[jwl - 5] = ipl2;
					iwl[jwl - 4] = ipl2;
					iwl[jwl - 3] = ipti2;
					iwl[jwl - 2] = ipti2;
					iwl[jwl - 1] = ipl1;
					for (jlt3 = 3; jlt3 <= nlt3; jlt3 += 3) {
						iplj1 = ipl[jlt3 - 3];
						iplj2 = ipl[jlt3 - 2];
						if ((iplj1 == ipl1 && iplj2 == ipti2) || (iplj2 ==
						    ipl1 && iplj1 == ipti2))
							ipl[jlt3 - 1] = itf[0];
						if ((iplj1 == ipl2 && iplj2 == ipti1) || (iplj2 ==
						    ipl2 && iplj1 == ipti1))
							ipl[jlt3 - 1] = itf[1];
					}
				}
		L_170:
				;
			}
			nlfc = nlf;
			nlf = jwl / 2;
			/* No more improvement necessary */
			if (nlf == nlfc)
				goto L_200;
			/*
			 * Resets the IWL array for the next round.
			 */
			jwl1mn = 2 * nlfc + 1;
			nlft2 = nlf * 2;
			for (jwl1 = jwl1mn; jwl1 <= nlft2; jwl1++) {
				jwl = jwl1 + 1 - jwl1mn;
				iwl[jwl - 1] = iwl[jwl1 - 1];
			}
			nlf = jwl / 2;
		}
		/* Done with improvement */
L_200:
		;
	}
	/*
	 * Rearranges the IPT array so that the vertexes of each triangle are
	 * listed counter-clockwise.
	 */
	for (itt3 = 3; itt3 <= ntt3; itt3 += 3) {
		ip1 = ipt[itt3 - 3];
		ip2 = ipt[itt3 - 2];
		ip3 = ipt[itt3 - 1];
		if (VPDT(*XYDATA(ip1 - 1, 0), *XYDATA(ip1 - 1, 1), *XYDATA(ip2 - 1, 0),
		*XYDATA(ip2 - 1, 1), *XYDATA(ip3 - 1, 0), *XYDATA(ip3 - 1, 1)) <
		    F_ZERO) {
			ipt[itt3 - 3] = ip2;
			ipt[itt3 - 2] = ip1;
		}
	}
	*nt = nt0;
	*nl = nl0;

L_9000:
	imsl_e1pop("l_s6rf");
	return;
}				/* end of function */
#undef	VPDT
#undef	SPDT
#undef	DSQF
/*----------------------------------------------------------------------- */

/*  IMSL Name:  S7RF/DS7RF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 12, 1986

    Purpose:    Compute a smooth bivariate interpolant to scattered data
                which is locally a quintic polynomial in two variables.

    Usage:      CALL  S7RF (XYDATA, NT, IPT, NL, IPL, NXOUT, NYOUT,
                            XOUT, YOUT, NGP, IGP)

    Arguments:
       XYDATA - A 2 by NDATA Array containing the coordinates of
                the interpolation points.  (Input)
                These points must be distinct.  The x-coordinate
                of the Ith data point is stored in XYDATA(1,I) and the
                y-coordinate of the Ith data point is stored in
                XYDATA(2,I).
       NT     - A counter.  (Input)
       IPT    - Work array of length 6*NDATA-14.
       NL     - A counter.  (Input)
       IPL    - Work array of length 6*NDATA.
       NXOUT  - The number of elements in XOUT.  (Input)
       NYOUT  - The number of elements in YOUT.  (Input)
       XOUT   - Array of length NXOUT containing an increasing sequence
                of points.  (Input)
                These points will be the x-coordinates of a grid on
                which the interpolated surface is to be evaluated.
       YOUT   - Array of length NYOUT containing an increasing sequence
                of points.  (Input)
                These points will be the y-coordinates of a grid on
                which the interpolated surface is to be evaluated.
       NGP    - Work array of length 6*NDATA-14.
       IGP    - Work array of length NXOUT*NYOUT.

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_s7rf(Mfloat *xydata, Mint *nt, Mint ipt[], Mint *nl, Mint ipl[],
                   Mint *nxout, Mint *nyout, Mfloat xout[], Mfloat yout[],
                   Mint ngp[], Mint igp[])
#else
static void l_s7rf(xydata, nt, ipt, nl, ipl, nxout, nyout, xout,
        	  yout, ngp, igp)
	Mfloat          *xydata;
	Mint            *nt, ipt[], *nl, ipl[], *nxout, *nyout;
	Mfloat           xout[], yout[];
	Mint             ngp[], igp[];
#endif
{
#define XYDATA(I_,J_)   (xydata+(I_)*(2)+(J_))
	Mint             il0, il0t3, ilp1, ilp1t3, insd, ip1, ip2, ip3,
	                it0, it0t3, ixi, iximn, iximx, iyi, izi, jigp0,
	                jigp1, jigp1i, jngp0, jngp1, l, ngp0, ngp1, nl0,
	                nt0, nxi0, nxinyi, nyi0;
	Mfloat           s12, s21, s32, v12, v23,
	                v31, x1, x2, x3, xii, ximn, ximx, xmn, xmx, y1,
	                y2, y3, yii, yimn, yimx, ymn, ymx;

	/* STATEMENT FUNCTIONS */
#define SPDT(u1,v1,u2,v2,u3,v3)	(Mfloat)(((u1) - (u2))*((u3) - \
	 (u2)) + ((v1) - (v2))*((v3) - (v2)))
#define VPDT(u1,v1,u2,v2,u3,v3)	(Mfloat)(((u1) - (u3))*((v2) - \
	 (v3)) - ((v1) - (v3))*((u2) - (u3)))

	nt0 = *nt;
	nl0 = *nl;
	nxi0 = *nxout;
	nyi0 = *nyout;
	nxinyi = nxi0 * nyi0;
	ximn = imsl_f_min(xout[0], xout[nxi0 - 1]);
	ximx = imsl_f_max(xout[0], xout[nxi0 - 1]);
	yimn = imsl_f_min(yout[0], yout[nyi0 - 1]);
	yimx = imsl_f_max(yout[0], yout[nyi0 - 1]);
	/*
	 * Determines grid points inside the data area.
	 */
	jngp0 = 0;
	jngp1 = 2 * (nt0 + 2 * nl0) + 1;
	jigp0 = 0;
	jigp1 = nxinyi + 1;
	for (it0 = 1; it0 <= nt0; it0++) {
		ngp0 = 0;
		ngp1 = 0;
		it0t3 = it0 * 3;
		ip1 = ipt[it0t3 - 3];
		ip2 = ipt[it0t3 - 2];
		ip3 = ipt[it0t3 - 1];
		x1 = *XYDATA(ip1 - 1, 0);
		y1 = *XYDATA(ip1 - 1, 1);
		x2 = *XYDATA(ip2 - 1, 0);
		y2 = *XYDATA(ip2 - 1, 1);
		x3 = *XYDATA(ip3 - 1, 0);
		y3 = *XYDATA(ip3 - 1, 1);
		xmn = imsl_f_vmin(3, x1, x2, x3);
		xmx = imsl_f_vmax(3, x1, x2, x3);
		ymn = imsl_f_vmin(3, y1, y2, y3);
		ymx = imsl_f_vmax(3, y1, y2, y3);
		insd = 0;
		for (ixi = 1; ixi <= nxi0; ixi++) {
			if (xout[ixi - 1] < xmn || xout[ixi - 1] > xmx) {
				if (insd != 0) {
					iximx = ixi - 1;
					goto L_20;
				}
			} else if (insd != 1) {
				insd = 1;
				iximn = ixi;
			}
		}
		if (insd == 0)
			goto L_60;
		iximx = nxi0;
L_20:
		for (iyi = 1; iyi <= nyi0; iyi++) {
			yii = yout[iyi - 1];
			if (yii >= ymn && yii <= ymx) {
				for (ixi = iximn; ixi <= iximx; ixi++) {
					xii = xout[ixi - 1];
					l = 0;
					v12 = VPDT(x1, y1, x2, y2, xii, yii);
					if (v12 >= F_ZERO) {
						if (v12 == F_ZERO)
							l = 1;
						v23 = VPDT(x2, y2, x3, y3, xii, yii);
						if (v23 >= F_ZERO) {
							if (v23 == F_ZERO)
								l = 1;
							v31 = VPDT(x3, y3, x1, y1, xii, yii);
							if (v31 >= F_ZERO) {
								if (v31 == F_ZERO)
									l = 1;
								izi = nxi0 * (iyi - 1) + ixi;
								if (l != 1) {
									ngp0 += 1;
									jigp0 += 1;
									igp[jigp0 - 1] = izi;
								} else {
									if (jigp1 <= nxinyi) {
										for (jigp1i = jigp1; jigp1i <= nxinyi; jigp1i++) {
											if (izi == igp[jigp1i - 1])
												goto L_40;
										}
									}
									ngp1 += 1;
									jigp1 -= 1;
									igp[jigp1 - 1] = izi;
								}
							}
						}
					}
			L_40:
					;
				}
			}
		}
L_60:
		jngp0 += 1;
		ngp[jngp0 - 1] = ngp0;
		jngp1 -= 1;
		ngp[jngp1 - 1] = ngp1;
	}
	/*
	 * Determines grid points outside the data area. - IN semi-infinite
	 * rectangular area.
	 */
	for (il0 = 1; il0 <= nl0; il0++) {
		ngp0 = 0;
		ngp1 = 0;
		il0t3 = il0 * 3;
		ip1 = ipl[il0t3 - 3];
		ip2 = ipl[il0t3 - 2];
		x1 = *XYDATA(ip1 - 1, 0);
		y1 = *XYDATA(ip1 - 1, 1);
		x2 = *XYDATA(ip2 - 1, 0);
		y2 = *XYDATA(ip2 - 1, 1);
		xmn = ximn;
		xmx = ximx;
		ymn = yimn;
		ymx = yimx;
		if (y2 >= y1)
			xmn = imsl_f_min(x1, x2);
		if (y2 <= y1)
			xmx = imsl_f_max(x1, x2);
		if (x2 <= x1)
			ymn = imsl_f_min(y1, y2);
		if (x2 >= x1)
			ymx = imsl_f_max(y1, y2);
		insd = 0;
		for (ixi = 1; ixi <= nxi0; ixi++) {
			if (xout[ixi - 1] < xmn || xout[ixi - 1] > xmx) {
				if (insd != 0) {
					iximx = ixi - 1;
					goto L_90;
				}
			} else if (insd != 1) {
				insd = 1;
				iximn = ixi;
			}
		}
		if (insd == 0)
			goto L_130;
		iximx = nxi0;
L_90:
		for (iyi = 1; iyi <= nyi0; iyi++) {
			yii = yout[iyi - 1];
			if (yii >= ymn && yii <= ymx) {
				for (ixi = iximn; ixi <= iximx; ixi++) {
					xii = xout[ixi - 1];
					l = 0;
					v12 = VPDT(x1, y1, x2, y2, xii, yii);
					if (v12 <= F_ZERO) {
						if (v12 == F_ZERO)
							l = 1;
						s21 = SPDT(x2, y2, x1, y1, xii, yii);
						if (s21 >= F_ZERO) {
							if (s21 == F_ZERO)
								l = 1;
							s12 = SPDT(x1, y1, x2, y2, xii, yii);
							if (s12 >= F_ZERO) {
								if (s12 == F_ZERO)
									l = 1;
								izi = nxi0 * (iyi - 1) + ixi;
								if (l != 1) {
									ngp0 += 1;
									jigp0 += 1;
									igp[jigp0 - 1] = izi;
								} else {
									if (jigp1 <= nxinyi) {
										for (jigp1i = jigp1; jigp1i <= nxinyi; jigp1i++) {
											if (izi == igp[jigp1i - 1])
												goto L_110;
										}
									}
									ngp1 += 1;
									jigp1 -= 1;
									igp[jigp1 - 1] = izi;
								}
							}
						}
					}
			L_110:
					;
				}
			}
		}
L_130:
		jngp0 += 1;
		ngp[jngp0 - 1] = ngp0;
		jngp1 -= 1;
		ngp[jngp1 - 1] = ngp1;
		/* - In semi-infinite triangular area. */
		ngp0 = 0;
		ngp1 = 0;
		ilp1 = mod(il0, nl0) + 1;
		ilp1t3 = ilp1 * 3;
		ip3 = ipl[ilp1t3 - 2];
		x3 = *XYDATA(ip3 - 1, 0);
		y3 = *XYDATA(ip3 - 1, 1);
		xmn = ximn;
		xmx = ximx;
		ymn = yimn;
		ymx = yimx;
		if (y3 >= y2 && y2 >= y1)
			xmn = x2;
		if (y3 <= y2 && y2 <= y1)
			xmx = x2;
		if (x3 <= x2 && x2 <= x1)
			ymn = y2;
		if (x3 >= x2 && x2 >= x1)
			ymx = y2;
		insd = 0;
		for (ixi = 1; ixi <= nxi0; ixi++) {
			if (xout[ixi - 1] < xmn || xout[ixi - 1] > xmx) {
				if (insd != 0) {
					iximx = ixi - 1;
					goto L_150;
				}
			} else if (insd != 1) {
				insd = 1;
				iximn = ixi;
			}
		}
		if (insd == 0)
			goto L_190;
		iximx = nxi0;
L_150:
		for (iyi = 1; iyi <= nyi0; iyi++) {
			yii = yout[iyi - 1];
			if (yii >= ymn && yii <= ymx) {
				for (ixi = iximn; ixi <= iximx; ixi++) {
					xii = xout[ixi - 1];
					l = 0;
					s12 = SPDT(x1, y1, x2, y2, xii, yii);
					if (s12 <= F_ZERO) {
						if (s12 == F_ZERO)
							l = 1;
						s32 = SPDT(x3, y3, x2, y2, xii, yii);
						if (s32 <= F_ZERO) {
							if (s32 == F_ZERO)
								l = 1;
							izi = nxi0 * (iyi - 1) + ixi;
							if (l != 1) {
								ngp0 += 1;
								jigp0 += 1;
								igp[jigp0 - 1] = izi;
							} else {
								if (jigp1 <= nxinyi) {
									for (jigp1i = jigp1; jigp1i <= nxinyi; jigp1i++) {
										if (izi == igp[jigp1i - 1])
											goto L_170;
									}
								}
								ngp1 += 1;
								jigp1 -= 1;
								igp[jigp1 - 1] = izi;
							}
						}
					}
			L_170:
					;
				}
			}
		}
L_190:
		jngp0 += 1;
		ngp[jngp0 - 1] = ngp0;
		jngp1 -= 1;
		ngp[jngp1 - 1] = ngp1;
	}
	return;
}				/* end of function */
#undef	VPDT
#undef	SPDT
/* Structured by FOR_STRUCT, v0.2, on 08/23/90 at 14:46:00
    Options SET: fmt=t s=n
  -----------------------------------------------------------------------
    IMSL Name:  S4RF/DS4RF (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    August 12, 1986

    Purpose:    Compute a smooth bivariate interpolant to scattered data
                which is locally a quintic polynomial in two variables.

    Usage:      CALL S4RF (NDATA, XYDATA, FDATA, PD)

    Arguments:
       NDATA  - Number of data points.  (Input)
                NDATA must be at least 3.
       XYDATA - A 2 by NDATA Array containing the coordinates of
                the interpolation points.  (Input)
                These points must be distinct.  The x-coordinate
                of the Ith data point is stored in XYDATA(1,I) and the
                y-coordinate of the Ith data point is stored in
                XYDATA(2,I).
       FDATA  - Array of length NDATA containing the interpolation
                values.  (Input)
       PD     - Work array of size 5 by NDATA.

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1986 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_s4rf(Mint *ndata, Mfloat *xydata, Mfloat fdata[], Mfloat pd[])
#else
static void l_s4rf(ndata, xydata, fdata, pd)
	Mint            *ndata;
	Mfloat           *xydata, fdata[], pd[];
#endif
{
/*      Originally the argument pd[] was pd[][5]. But l_s2rf is sending
        wk[] so pd[][5] became pd[] and an appropriate #define was
        added to do the address arithmetic.
*/
#define PD(I_,J_)      (pd+(I_)*5 + (J_))
/*      Originally there was an array dimensioned as c(12,5)
        in the fortran subroutine. When imsl_l2qrr is called this
        array is received as a pointer to a float which caused a
        problem using Gnu C. (note the imsl_malloc call added)
*/
        Mfloat          *c = NULL;
#define C(I_,J_)	(c+(I_)*(12)+(J_))
#define XYDATA(I_,J_)   (xydata+(I_)*(2)+(J_))
	Mint             i,ii,  iwk[5], j, k, kbasis, na, npts;
	Mfloat           aclose[12], b[12], bb[5], big, bigx, bigy, d, eps,
	                r, reps, res[12], rmd, temp, tol, wkar2[10], wkarea[5],
	                xdd, ydd;
	static Mint      nc = 5;
	static Mint      ic = 12;


        c = (Mfloat *)imsl_malloc(60 * sizeof(*c));
	big = F_ZERO;
	reps = imsl_amach(4);
	for (i = 2; i <= *ndata; i++) {
		for (j = 1; j <= (i - 1); j++) {
			big = imsl_f_max(big, imsl_fi_power(*XYDATA(i - 1, 0) - *XYDATA(j - 1, 0), 2) +
					 imsl_fi_power(*XYDATA(i - 1, 1) - *XYDATA(j - 1, 1), 2));
		}
	}
	eps = 100.0 * sqrt(big) * reps;
	for (i = 1; i <= *ndata; i++) {
		/*
		 * Find radius which includes 10 neighbors
		 */
		sset(ic, big, aclose, 1);
		for (j = 1; j <= *ndata; j++) {
			if (i != j) {
				xdd = *XYDATA(j - 1, 0) - *XYDATA(i - 1, 0);
				ydd = *XYDATA(j - 1, 1) - *XYDATA(i - 1, 1);
				d = imsl_fi_power(xdd, 2) + imsl_fi_power(ydd, 2);
				if (d < aclose[ic - 1]) {
					aclose[ic - 1] = d;
					for (k = 1; k <= ic; k++) {
						if (aclose[k - 1] > aclose[ic - 1]) {
							temp = aclose[k - 1];
							aclose[k - 1] = aclose[ic - 1];
							aclose[ic - 1] = temp;
						}
					}
				}
			}
		}
		r = sqrt(aclose[ic - 1]);
		npts = 0;
		bigx = F_ZERO;
		bigy = F_ZERO;
		for (j = 1; j <= *ndata; j++) {
			if (i != j) {
				xdd = *XYDATA(j - 1, 0) - *XYDATA(i - 1, 0);
				ydd = *XYDATA(j - 1, 1) - *XYDATA(i - 1, 1);
				d = imsl_fi_power(xdd, 2) + imsl_fi_power(ydd, 2);
				if (d <= aclose[ic - 1]) {
					d = sqrt(d);
					rmd = r / d;
					bigx = imsl_f_max(bigx, fabs(xdd) * rmd);
					bigy = imsl_f_max(bigy, fabs(ydd) * rmd);
					npts += 1;
					*C(0, npts - 1) = rmd * xdd;
					*C(1, npts - 1) = rmd * ydd;
					*C(2, npts - 1) = *C(0, npts - 1) * xdd / F_TWO;
					*C(3, npts - 1) = *C(0, npts - 1) * ydd;
					*C(4, npts - 1) = *C(1, npts - 1) * ydd / F_TWO;
					b[npts - 1] = (fdata[j - 1] - fdata[i - 1]) * rmd;
					if (npts >= ic)
						goto L_60;
				}
			}
		}
L_60:
		na = nc;
		if (npts < 5)
			na = 2;
		kbasis = 0;
		tol = eps / r;
		/*
		 * Scale matrix to include 1st and 2nd columns.
		 */
		if (bigx > F_ONE)
			sscal(npts, bigx, C(0, 0), 1);
		if (bigy > F_ONE)
			sscal(npts, bigy, C(1, 0), 1);
		/*
		 * Least squares fit to neighbors using quadratic polynomial
		 */
#if defined(COMPUTER_VAXG)
                /* This was inserted because the code was so flakey on vaxg. */
                ii = 0;
#endif
		iset(na, 0, iwk, 1);
		imsl_l2qrr(&npts, &na, c, &ic, b, &tol, bb, res, &kbasis, c, wkarea,
			   iwk, wkar2);
		/* Rescale answers */
		if (bigx > F_ONE)
			bb[0] *= bigx;
		if (bigy > F_ONE)
			bb[1] *= bigy;

		if (imsl_n1rty(1) != 0)
			na = 0;
		for (j = 1; j <= 5; j++) {
/*			pd[i - 1][j - 1] = 0.0;  */
                                *PD(i - 1, j - 1) = F_ZERO;
			if (j <= na){
/*				pd[i - 1][j - 1] = bb[j - 1]; */
                                *PD(i - 1, j - 1) = bb[j - 1];
				}
		}
	}
        if (c != NULL) imsl_free(c);
	return;
}				/* end of function */
#undef  C
#undef  PD
