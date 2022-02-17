#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef ANSI
Mfloat imsl_f_spline_2d_integral(Mfloat a, Mfloat b, Mfloat c, Mfloat d,Mf_spline *pp)
#else
Mfloat imsl_f_spline_2d_integral(a,b,c,d,pp)
   Mfloat        a;
   Mfloat        b;
   Mfloat        c;
   Mfloat        d;
   Mf_spline    *pp;
#endif
{
   Mfloat        value=F_ZERO;
   Mfloat        temp_a, temp_b;
   Mfloat        temp_c, temp_d;
   Mint          temp_int;
   Mfloat        *wk1 = NULL;

#ifdef DOUBLE
   imsl_e1psh("imsl_d_spline_2d_integral");
#else
   imsl_e1psh("imsl_f_spline_2d_integral");
#endif
   temp_a= (Mfloat)a;
   temp_b= (Mfloat)b;
   temp_c= (Mfloat)c;
   temp_d= (Mfloat)d; 
        /* Check KXORD */
        if (pp->order[0] < 1) {
                imsl_e1sti(1, pp->order[0]);
 
                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_ORDER_X);
        }
        /* Check KYORD */
        if (pp->order[1] < 1) {
                imsl_e1sti(1, pp->order[1]);
                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_ORDER_Y);
        }
        if (imsl_n1rty(0) != 0)
                goto FREE_SPACE;
        /* Check NYCOEF */
        if (pp->num_coef[1] < pp->order[1]) {
                imsl_e1sti(1, pp->num_coef[1]);
                imsl_e1sti(2, pp->order[1]);

                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_COEFF_Y);
        }
        if (imsl_n1rty(0) != 0)
                goto FREE_SPACE;

   temp_int = 4*(imsl_i_max(pp->order[0],pp->order[1])+1)+pp->num_coef[1];
   wk1  = (Mfloat *)imsl_malloc(temp_int*sizeof(*wk1));

   if (wk1 == NULL) {
            imsl_e1stl(1, "x_order");
            imsl_e1sti(1, pp->order[0]);
            imsl_e1stl(1, "y_order");
            imsl_e1sti(1, pp->order[1]);
            imsl_e1stl(1, "nycoef");
            imsl_e1sti(1, pp->num_coef[1]);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_3);
            goto FREE_SPACE; }


   value = imsl_b22ig(&temp_a,&temp_b,&temp_c,&temp_d,&(pp->order[0]),&(pp->order[1]),
                      pp->knots[0],pp->knots[1],&(pp->num_coef[0]),&(pp->num_coef[1]),
                      pp->coef[0],wk1);
   if (imsl_n1rty(1) > 3) value = imsl_amach(6);

FREE_SPACE:
   if ( wk1 != NULL)         imsl_free(wk1);
#ifdef DOUBLE
   imsl_e1pop("imsl_d_spline_2d_integral");
#else
   imsl_e1pop("imsl_f_spline_2d_integral");
#endif
   return (value);
}

