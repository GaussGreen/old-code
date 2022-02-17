#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef ANSI
Mfloat imsl_f_spline_integral(Mfloat a, Mfloat b, Mf_spline *pp)
#else
Mfloat imsl_f_spline_integral(a,b,pp)
   Mfloat        a;
   Mfloat        b;
   Mf_spline    *pp;
#endif
{
   Mfloat        value = F_ZERO;
   Mfloat        temp_a;
   Mfloat        temp_b;
   Mfloat        *wk1;
   Mfloat        *wk2;
   Mfloat        *wk3;
   Mfloat        *wk4;

#ifdef DOUBLE
   imsl_e1psh("imsl_d_spline_integral");
#else
   imsl_e1psh("imsl_f_spline_integral");
#endif
   temp_a= (Mfloat)a;
   temp_b= (Mfloat)b;
   wk1  = (Mfloat *)imsl_malloc(((pp->order[0])+1)*sizeof(*wk1));
   wk2  = (Mfloat *)imsl_malloc(((pp->order[0])+1)*sizeof(*wk2));
   wk3  = (Mfloat *)imsl_malloc(((pp->order[0])+1)*sizeof(*wk3));
   wk4  = (Mfloat *)imsl_malloc(((pp->order[0])+1)*sizeof(*wk4));

   if ((wk4 == NULL)||(wk3 == NULL)||(wk2 == NULL)||(wk1 == NULL)) {
            imsl_e1stl(1, "order");
            imsl_e1sti(1, pp->order[0]);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
            goto FREE_SPACE; }


   value = imsl_b2itg(&temp_a,&temp_b,pp->order,pp->knots[0],pp->num_coef,
                      pp->coef[0], wk1,wk2,wk3,wk4);
   if (imsl_n1rty(1) > 3) value = imsl_amach(6);
FREE_SPACE:
   if ( wk4 != NULL)         imsl_free(wk4);
   if ( wk3 != NULL)         imsl_free(wk3);
   if ( wk2 != NULL)         imsl_free(wk2);
   if ( wk1 != NULL)         imsl_free(wk1);
#ifdef DOUBLE
   imsl_e1pop("imsl_d_spline_integral");
#else
   imsl_e1pop("imsl_f_spline_integral");
#endif
   return (value);
}

