#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
Mfloat imsl_f_cub_spline_integral(Mfloat a, Mfloat b, Mf_ppoly *pp)
#else
Mfloat imsl_f_cub_spline_integral(a,b,pp)
   Mfloat        a;
   Mfloat        b;
   Mf_ppoly    *pp;
#endif
{
   Mfloat        value;
   Mfloat        temp_a, temp_b;
   Mint          nintv;


#ifdef DOUBLE
   imsl_e1psh("imsl_d_cub_spline_integral");
#else
   imsl_e1psh("imsl_f_cub_spline_integral");
#endif
   temp_a= (Mfloat)a;
   temp_b= (Mfloat)b;
   nintv = (pp->num_breakpoints[0]) - 1;
   value = imsl_csitg(&temp_a,&temp_b,&nintv,pp->breakpoints[0],pp->coef[0]);
   if (imsl_n1rty(1) > 3) value = imsl_amach(6);

#ifdef DOUBLE
   imsl_e1pop("imsl_d_cub_spline_integral");
#else
   imsl_e1pop("imsl_f_cub_spline_integral");
#endif
   return (value);
}

