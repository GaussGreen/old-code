#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
static Mf_ppoly   *lv_pp;
static VA_LIST_HACK  PROTO(l_cub_spline_interp_e_cnd,(Mint, Mfloat[], Mfloat[], va_list));
static void     PROTO(l_c2per,(Mint*, Mfloat[], Mfloat[],Mfloat[], Mfloat*, Mfloat*, Mint[]));
                       
#ifdef ANSI
Mf_ppoly  *imsl_f_cub_spline_interp_e_cnd(Mint ndata, Mfloat xdata[], 
                    Mfloat fdata[], ...)
#else
Mf_ppoly  *imsl_f_cub_spline_interp_e_cnd(ndata,xdata,fdata,va_alist)
   Mint          ndata;
   Mfloat        xdata[];
   Mfloat        fdata[];
   va_dcl
#endif
{
    va_list     argptr;
    VA_START(argptr, fdata);

    E1PSH("imsl_f_cub_spline_interp_e_cnd","imsl_d_cub_spline_interp_e_cnd");

    lv_pp = NULL;
    IMSL_CALL(l_cub_spline_interp_e_cnd(ndata,xdata,fdata,argptr));
    va_end(argptr);
    E1POP("imsl_f_cub_spline_interp_e_cnd","imsl_d_cub_spline_interp_e_cnd");
    return lv_pp;
}


#ifdef ANSI
static VA_LIST_HACK   l_cub_spline_interp_e_cnd(Mint ndata, Mfloat xdata[],
                        Mfloat fdata[], va_list argptr)
#else
static VA_LIST_HACK   l_cub_spline_interp_e_cnd(ndata,xdata,fdata,argptr)
   Mint          ndata;
   Mfloat        xdata[];
   Mfloat        fdata[];
   va_list       argptr;
#endif
{
   Mint          tmp;
   Mint          code;
   Mint          i;
   Mint          arg_number         = 3;
   Mint          four               = 4;
   Mint          periodic_set       = 0;         /* DEFAULT VALUE */
   Mint          ileft              = 0;         /* DEFAULT VALUE */
   Mint          iright             = 0;         /* DEFAULT VALUE */
   Mfloat        dleft              = F_ZERO;       /* DEFAULT VALUE */
   Mfloat        dright             = F_ZERO;       /* DEFAULT VALUE */
   Mint          domain_dim         = 1;
   Mint          target_dim         = 1;
   Mint          *orders            = NULL;
   Mint          *num_breakpoints   = NULL;
   Mfloat        *coef_work         = NULL;
   Mfloat        *rwork             = NULL;
   Mint          *iwork             = NULL;
   Mint          free_the_structure = 0;
    code = 1;
    while (code > 0) {
        code = va_arg(argptr, Mint);
        arg_number++;
        switch (code) {
            case IMSL_LEFT:
                ileft = va_arg(argptr, Mint);
                arg_number++;
                dleft = (Mfloat)va_arg(argptr, Mdouble);
                arg_number++;
                break;
            case IMSL_RIGHT:
                iright = va_arg(argptr, Mint);
                arg_number++;
                dright =(Mfloat) va_arg(argptr, Mdouble);
                arg_number++;
                break;
            case IMSL_LEFT_ADR:
                ileft = *va_arg(argptr, Mint *);
                arg_number++;
                dleft = *(va_arg(argptr, Mfloat *));
                arg_number++;
                break;
            case IMSL_RIGHT_ADR:
                iright = *va_arg(argptr, Mint *);
                arg_number++;
                dright = *( va_arg(argptr, Mfloat *));
                arg_number++;
                break;
            case IMSL_PERIODIC:
                periodic_set = 1;
                break;
            case 0:
                break;
            default:
                imsl_e1sti (1, code);
                imsl_e1sti (2, arg_number);
                imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
                goto RETURN;
        }
    } 
   if(periodic_set == 1){
       
   orders = &four;
   tmp = ndata;
   num_breakpoints = &tmp;
                                   /* CREATE THE STRUCTURE */         
   lv_pp =   imsl_f_ppoly_create(domain_dim,target_dim,orders, num_breakpoints,0);
   if (imsl_n1rty(1)==4){
            imsl_e1mes(0,0," ");
            imsl_e1stl(1, "ndata");
            imsl_e1sti(1, ndata);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
       goto RETURN;
   }
        if (ndata < 4) {
                free_the_structure = 1;
                imsl_e1sti(1, ndata);
/*                imsl_ermes(5, 1, "The number of data points must be 4 or more while NDATA = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_AT_LEAST_4_PTS);
                goto PERIODIC_FREE_SPACE; 
        }
                                   /* GET THE WORKSPACE    */
                              /*NEED TO USE WORKSPACE FOR THE COEFFICIENTS*/
                               /* SINCE THE SPACE NEEDED BY C2INT IS LARGER
                                  THAN THE SPACE PROVIDED IN THE STRUCTURE
                                  IMSL_PPOLY. THIS IS CAUSED BY C2INT USING
                                  THE LAST 4 MEMORY SPACES OF CSCOEF FOR 
                                  WORKSPACE  */
   coef_work = (Mfloat *)imsl_malloc(4*(ndata)*sizeof(*coef_work));
   rwork = (Mfloat *)imsl_malloc(6*(ndata)*sizeof(*rwork));
   iwork = (Mint *)imsl_malloc(ndata*sizeof(*iwork));
   if ((iwork == NULL) ||(rwork == NULL) ||(coef_work == NULL) ){
            free_the_structure = 1;
            imsl_e1stl(1, "ndata");
            imsl_e1sti(1, ndata);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
            goto PERIODIC_FREE_SPACE; }

                                   /* CALL THE SPLINE ROUTINE */ 
   l_c2per(&ndata,xdata,fdata,lv_pp->breakpoints[0],coef_work,
        rwork,iwork);
   if (imsl_n1rty(1) >3) {
       free_the_structure = 1;
       goto PERIODIC_FREE_SPACE; 
   }
                               /* COPY THE COEFFICIENTS INTO THE STRUCTURE */
      for (i=0;i<lv_pp->num_coef[0];i++)
       lv_pp->coef[0][i]=coef_work[i];
     goto  PERIODIC_FREE_SPACE; 
   }
   /* THIS ELSE CLAUSE IS FOR THE CASE PERIODIC IS NOT SPECIFIED */
   else{
   orders = &four;
   tmp = ndata;
   num_breakpoints = &tmp;
                                   /* CREATE THE STRUCTURE */         
   lv_pp =   imsl_f_ppoly_create(domain_dim,target_dim,orders, num_breakpoints,0);
   if (imsl_n1rty(1)==4){
            imsl_e1mes(0,0," ");
            imsl_e1stl(1, "ndata");
            imsl_e1sti(1, ndata);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
            goto RETURN;
   }
                                   /* GET THE WORKSPACE    */
                                   /*NEED TO USE WORKSPACE FOR THE COEFFICIENTS
                                     SINCE THE SPACE NEEDED BY C2INT IS LARGER
                                     THAN THE SPACE PROVIDED IN THE STRUCTURE
                                     IMSL_PPOLY. THIS IS CAUSED BY C2INT USING
                                     THE LAST 4 MEMORY SPACES OF CSCOEF FOR 
                                     WORKSPACE  */
        if (ndata < 2) {
                free_the_structure = 1;
                imsl_e1sti(1, ndata);
/*                imsl_ermes(5, 1, "The number of data points must be 2 or more while NDATA = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_AT_LEAST_2_PTS);
                goto PERIODIC_FREE_SPACE; 
        }
   coef_work = (Mfloat *)imsl_malloc(4*ndata*sizeof(*coef_work));
   iwork = (Mint *)imsl_malloc(ndata*sizeof(*iwork));
   if ((iwork == NULL) ||( coef_work == NULL)){
            free_the_structure = 1;
            imsl_e1stl(1, "ndata");
            imsl_e1sti(1, ndata);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
            goto END_CND_FREE_SPACE; }
                                   /* CALL THE SPLINE ROUTINE */ 
   imsl_c2dec(&ndata,xdata,fdata,&ileft,&dleft,&iright,&dright,lv_pp->breakpoints[0],
              coef_work,iwork);
   if (imsl_n1rty(1) >3) {
       free_the_structure = 1;
       goto END_CND_FREE_SPACE; 
   }
                                  /* COPY THE COEFFICIENTS INTO THE STRUCTURE */
      for (i=0;i<lv_pp->num_coef[0];i++)
       lv_pp->coef[0][i]=coef_work[i];
   goto END_CND_FREE_SPACE;
   }
                                    /*  FREE THE WORKSPACE USED */
PERIODIC_FREE_SPACE:
   if ((free_the_structure == 1)&&(lv_pp!= NULL)) {
	imsl_free(lv_pp);
	lv_pp = NULL;
    }
   if ( iwork != NULL) imsl_free(iwork);
   if ( rwork != NULL) imsl_free(rwork);
   if ( coef_work != NULL) imsl_free(coef_work);
   goto RETURN;
END_CND_FREE_SPACE:
   if ((free_the_structure == 1)&&(lv_pp != NULL)) {
	imsl_free(lv_pp);
	lv_pp = NULL;
    }
   if ( iwork != NULL)  imsl_free(iwork);
   if ( coef_work != NULL)  imsl_free(coef_work);
   goto RETURN;
RETURN:
    return argptr;
}




/*  -----------------------------------------------------------------------
    IMSL Name:  C2PER/DC2PER (Single/Double precision version)

    Computer:   FORC/SINGLE

    Revised:    October 1, 1984

    Purpose:    Compute the cubic spline interpolant with periodic
                boundary conditions.

    Usage:      CALL C2PER (NDATA, XDATA, FDATA, BREAK, CSCOEF,
                            WORK, IPVT)

    Arguments:
       NDATA  - Number of data points.  (Input)
                NDATA must be at least 4.
       XDATA  - Array of length NDATA containing the data point
                abscissas.  (Input)
                The data point abscissas must be distinct.
       FDATA  - Array of length NDATA containing the data point
                ordinates.  (Input)
       BREAK  - Array of length NDATA containing the breakpoints for the
                piecewise cubic representation.  (Output)
       CSCOEF - Matrix of size 4 by NDATA containing the local
                coefficients of the cubic pieces.  (Output)
       WORK   - Work array of length 6*NDATA.
       IPVT   - Work array of length NDATA.

    Remark:
       Informational error
       Type Code
         3   1  The data set is not periodic, i.e., the function values
                at the smallest and largest XDATA points are not equal.
                The value at the smallest XDATA point is used.

    GAMS:       E1a

    Chapter:    MATH/LIBRARY Interpolation and Approximation

    Copyright:  1985 by IMSL, Inc.  All Rights Reserved.

    Warranty:   IMSL warrants only that IMSL testing has been applied
                to this code.  No other warranty, expressed or implied,
                is applicable.

  -----------------------------------------------------------------------
 */
#ifdef ANSI
static void l_c2per(Mint *ndata, Mfloat xdata[], Mfloat fdata[], Mfloat break_[],
           Mfloat *cscoef, Mfloat *work,  Mint ipvt[])
#else
static void l_c2per(ndata, xdata, fdata, break_, cscoef, work,  ipvt)
	Mint            *ndata;
	Mfloat           xdata[], fdata[], break_[], *cscoef, *work;
	Mint             ipvt[];
#endif
{
#define CSCOEF(I_,J_)	(cscoef+((I_)*(4)+(J_)))
#define WORK(I_,J_)	(work+((I_)*(6)+(J_)))
	Mint             i, im1, imj, ip1, ip2, j, n;
	Mfloat           a2, a4, b2, b3, b4, c2, c3, c4, d1, d2, d3, d4,
	                det1, det2, det3, det4, h[7], hh[6], ht1, ht2,
	                ht3, ht4, r;


	imsl_e1psh("L_C2PER ");
	/* CHECK ARGUMENT NDATA */

	if (*ndata < 4) {
		imsl_e1sti(1, *ndata);

/*		imsl_ermes(5, 1, "The number of data points must be 4 or more while NDATA = %(i1) is given.");
*/
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_AT_LEAST_4_PTS);
	}
	/* CHECK FOR ERRORS */
	if (imsl_n1rty(0) != 0)
		goto L_9000;
	/* COPY AND SORT INPUT DATA */
	imsl_c1sor(*ndata, xdata, fdata, break_, cscoef, 4, ipvt);
	if (imsl_n1rty(1) != 0)
		goto L_9000;

	n = *ndata - 1;
	if (*CSCOEF(0, 0) != *CSCOEF(*ndata - 1, 0)) {
		imsl_e1str(1, break_[0]);
		imsl_e1str(2, break_[*ndata - 1]);
		imsl_e1str(3, *CSCOEF(0, 0));
		imsl_e1str(4, *CSCOEF(*ndata - 1, 0));

/*		imsl_ermes(3, 1, "The function values at smallest and largest XDATA points, %(r1) and %(r2), are not equal.  They are %(r3) and %(r4),  respectively.  %(r3) is used.");
*/
                imsl_ermes(IMSL_WARNING, IMSL_NOT_PERIODIC);
		*CSCOEF(*ndata - 1, 0) = *CSCOEF(0, 0);
	}
	/*
	 * CALCULATE COEFFICIENTS OF MATRIX A WHERE A(I,J)=JTH SPLINE BASIS
	 * FUNCTION EVALUATED AT BREAK(I)
	 */
	for (i = 1; i <= n; i++) {
		for (j = 1; j <= 4; j++) {
			imj = mod(n + i + j - 4, n) + 1;
			h[j - 1] = break_[imj] - break_[imj - 1];
		}
		ht1 = h[2] + h[3];
		ht2 = h[0] + h[1];
		ht3 = h[1] + h[2];
		det1 = (ht1 + h[2]) * h[1] * ht2 + (ht2 + h[1]) * h[2] * ht1;
		*WORK(i - 1, 0) = F_ZERO;
		a2 = imsl_fi_power(h[0], 2) * (ht3 + h[3]) * ht3 / ht2 / det1;
		*WORK(i - 1, 1) = a2;
		*WORK(i - 1, 2) = F_ONE;
		a4 = imsl_fi_power(h[3], 2) * ht3 * (h[0] + ht3) / ht1 / det1;
		*WORK(i - 1, 3) = a4;
		*WORK(i - 1, 4) = F_ZERO;
		*WORK(i - 1, 5) = *CSCOEF(i - 1, 0);
	}
	*WORK(0, 0) = *WORK(n - 1, 3);
	*WORK(0, 4) = *WORK(0, 1);
	*WORK(n - 2, 0) = *WORK(n - 1, 1);
	*WORK(n - 2, 4) = *WORK(n - 2, 3);
	*WORK(n - 1, 0) = F_ONE;
	/*
	 * SOLVE LINEAR SYSTEM FOR COEFFICIENTS OF SPLINE BASIS FUNCTIONS
	 */
	for (i = 1; i <= (n - 1); i++) {

		r = *WORK(i - 1, 4) / *WORK(i - 1, 2);
		*WORK(n - 1, 0) += -r ** WORK(i - 1, 0);
		*WORK(n - 1, 5) += -r ** WORK(i - 1, 5);
		if (i == n - 1)
			goto L_30;
		*WORK(i, 4) += -r ** WORK(i, 1);
		r = *WORK(i - 1, 3) / *WORK(i - 1, 2);
		*WORK(i, 2) += -r ** WORK(i, 1);
		*WORK(i, 0) += -r ** WORK(i - 1, 0);
		*WORK(i, 5) += -r ** WORK(i - 1, 5);
L_30:
		;
	}
	*WORK(n - 1, 5) /= *WORK(n - 1, 0);
	*WORK(n - 2, 5) = (*WORK(n - 2, 5) - *WORK(n - 2, 0) ** WORK(n - 1, 5)) /
		*WORK(n - 2, 2);
	for (i = n - 2; i >= 1; i--) {
		*WORK(i - 1, 5) = (*WORK(i - 1, 5) - *WORK(i, 1) ** WORK(i, 5) -
		       *WORK(i - 1, 0) ** WORK(n - 1, 5)) / *WORK(i - 1, 2);
	}
	/*
	 * IN EACH INTERVAL CALCULATE COEFFICIENTS CSCOEF(J,I),J=3,5 FROM
	 * BASIS FUNCTION COEFFICIENTS
	 */
	for (i = 1; i <= n; i++) {
		im1 = mod(n + i - 2, n) + 1;
		ip1 = mod(n + i, n) + 1;
		ip2 = mod(n + i + 1, n) + 1;
		for (j = 1; j <= 7; j++) {
			imj = mod(n + i + j - 5, n) + 1;
			h[j - 1] = break_[imj] - break_[imj - 1];
		}
		for (j = 1; j <= 6; j++) {
			hh[j - 1] = h[j - 1] + h[j];
		}
		ht1 = h[3] + hh[3];
		ht2 = h[2] + hh[1];
		ht3 = h[4] + hh[4];
		ht4 = h[3] + hh[2];
		det1 = (h[2] + hh[2]) * h[1] * hh[0] + (h[1] + hh[0]) * h[2] * hh[2];
		det2 = ht1 * h[2] * hh[1] + ht2 * h[3] * hh[3];
		det3 = ht3 * h[3] * hh[2] + ht4 * h[4] * hh[4];
		det4 = (h[5] + hh[5]) * h[4] * hh[3] + (h[4] + hh[3]) * h[5] * hh[5];
		d1 = hh[4] * (hh[4] + h[6]) / h[3] / hh[3] / det4;
		b2 = F_THREE * h[2] * hh[3] * (hh[3] + h[5]) / hh[2] / det3;
		c2 = b2 / h[2];
		d2 = -(ht3 * ht4 + F_THREE * hh[2] * h[3] + h[2] * h[2] + h[4] * hh[4]) /
			hh[2] / h[3] / det3;
		b3 = F_THREE * (h[3] * hh[3] - h[2] * hh[1]) / det2;
		c3 = -F_THREE * (ht1 + ht2) / det2;
		d3 = (ht1 * ht2 + F_THREE * h[3] * hh[3] + h[4] * h[4] + h[2] * hh[1]) / hh[3] /
			h[3] / det2;
		b4 = -F_THREE * h[3] * hh[1] * (hh[0] + h[2]) / hh[2] / det1;
		c4 = -b4 / h[3];
		d4 = -c4 / h[3] / F_THREE;
		*CSCOEF(i - 1, 1) = *WORK(im1 - 1, 5) * b4 + *WORK(i - 1, 5) * b3 +
			*WORK(ip1 - 1, 5) * b2;
		*CSCOEF(i - 1, 2) = F_TWO * (*WORK(im1 - 1, 5) * c4 + *WORK(i - 1, 5) *
					   c3 + *WORK(ip1 - 1, 5) * c2);
		*CSCOEF(i - 1, 3) = F_SIX * (*WORK(im1 - 1, 5) * d4 + *WORK(i - 1, 5) *
		      d3 + *WORK(ip1 - 1, 5) * d2 + *WORK(ip2 - 1, 5) * d1);
	}

L_9000:
	;
	imsl_e1pop("L_C2PER ");

	return;
}				/* end of function */
