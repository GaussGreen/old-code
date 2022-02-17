#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static VA_LIST_HACK PROTO(l_cub_spline_value,(Mfloat temp_x, Mf_ppoly *pp, va_list argptr));
static Mfloat  value;
#ifdef ANSI
Mfloat imsl_f_cub_spline_value(Mfloat temp_x, Mf_ppoly *pp,...)
#else
Mfloat imsl_f_cub_spline_value(temp_x, pp, va_alist)
    Mfloat    temp_x;
    Mf_ppoly  *pp;
    va_dcl
#endif
{
    Mfloat    x;
    va_list     argptr;
    value = imsl_amach(6);
    VA_START(argptr,pp);

#ifdef DOUBLE
    imsl_e1psh("imsl_d_cub_spline_value");
#else
    imsl_e1psh("imsl_f_cub_spline_value");
#endif
    x = (Mfloat)temp_x;
    IMSL_CALL(l_cub_spline_value(x,pp,argptr));
    va_end(argptr);
#ifdef DOUBLE
    imsl_e1pop("imsl_d_cub_spline_value");
#else
    imsl_e1pop("imsl_f_cub_spline_value");
#endif
    return (value);

}

#ifdef ANSI
static VA_LIST_HACK l_cub_spline_value(Mfloat temp_x, Mf_ppoly *pp, va_list argptr)
#else
static VA_LIST_HACK l_cub_spline_value(temp_x,pp,argptr)
    Mfloat       temp_x;
    Mf_ppoly     *pp;
    va_list      argptr;
#endif
{
    Mfloat          x;
    Mint            arg_number = 2;
    Mint            code;
    Mint            i;
    Mint            num_intervals;
    Mint            derivative=0;
    Mfloat          *value_ptr=NULL;

    Mint	    length;
    Mfloat	    *input_vector = NULL;
    Mfloat	    **output_vector = NULL;
    Mfloat	    *output_vector_user = NULL;
    Mint	    grid = 0;
    Mint	    grid_user = 0;
    Mint	    *iwk = NULL;
    Mfloat	    *work1 = NULL;
    Mfloat	    *work2 = NULL;

    x = (Mfloat)temp_x;
    code = 1;
    while (code > 0) {
        code = va_arg(argptr, Mint);
        arg_number++;
        switch (code) {
            case IMSL_DERIV:
                derivative = va_arg(argptr, Mint);
                arg_number++;
                break;
	    case IMSL_GRID:
		length = va_arg(argptr, Mint);
		input_vector = va_arg(argptr, Mfloat*);
		output_vector = va_arg(argptr, Mfloat**);
		grid = 1;
		arg_number += 3;
		break;
	    case IMSL_GRID_USER:
		length = va_arg(argptr, Mint);
		input_vector = va_arg(argptr, Mfloat *);
		output_vector_user = va_arg(argptr, Mfloat *);
		grid_user = 1;
		arg_number += 3;
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
    if (imsl_n1rty(0)) goto RETURN;

    if (!grid && !grid_user) {
       num_intervals = (pp->num_breakpoints[0]) - 1;
       value =imsl_csder(&derivative,&x,&num_intervals,pp->breakpoints[0],pp->coef[0]);
       if (imsl_n1rty(1) > 3) value = imsl_amach(6);
     }

    if (grid) {
	iwk = (Mint *) imsl_malloc (length * sizeof(*iwk));
	work1 = (Mfloat *) imsl_malloc (length * sizeof(*work1));
	work2 = (Mfloat *) imsl_malloc (length * sizeof(*work2));
	*output_vector = (Mfloat *) imsl_malloc (length * sizeof(Mfloat));
	if ( (*output_vector == NULL) || (iwk == NULL) || (work1 == NULL)
		|| (work2 == NULL) ) {
		imsl_e1stl(1, "length");
		imsl_e1sti(1, length);
		imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
		goto FREE_SPACE;
	}
	num_intervals = (pp->num_breakpoints[0]) - 1;
	imsl_c21gd (&derivative, &length, input_vector, &num_intervals,
		pp->breakpoints[0], (Mfloat**) pp->coef[0], *output_vector, iwk, 
		work1, work2);
	if (imsl_n1rty(1) > 3) {
	    for (i = 0; i<length; i++) *(*(output_vector+i)) = imsl_amach(6);
	}
    }

    if (grid_user) {
        iwk = (Mint *) imsl_malloc (length * sizeof(*iwk));
        work1 = (Mfloat *) imsl_malloc (length * sizeof(*work1));
        work2 = (Mfloat *) imsl_malloc (length * sizeof(*work2));
        if ( (iwk == NULL) || (work1 == NULL) || (work2 == NULL) ) {
                imsl_e1stl(1, "length");
                imsl_e1sti(1, length);
                imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
                goto FREE_SPACE;
        }        
        num_intervals = (pp->num_breakpoints[0]) - 1;
        imsl_c21gd (&derivative, &length, input_vector, &num_intervals,
                pp->breakpoints[0], (Mfloat**)pp->coef[0], output_vector_user, iwk,
                work1, work2);
        if (imsl_n1rty(1) > 3) {
            for (i = 0; i<length; i++) *(output_vector_user +i) = imsl_amach(6);
        }
    }

FREE_SPACE:
;

RETURN:
    return (argptr);
}    
