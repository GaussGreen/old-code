#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
static VA_LIST_HACK PROTO( l_spline_value,(Mfloat x, Mf_spline *pp, va_list argptr));
static Mfloat  lv_value;

#ifdef ANSI
Mfloat imsl_f_spline_value(Mfloat tempx, Mf_spline *pp,...)
#else
Mfloat imsl_f_spline_value(tempx,pp,va_alist)
    Mfloat    tempx;
    Mf_spline *pp;
    va_dcl
#endif
{
    Mfloat      x;
    va_list     argptr;
    VA_START(argptr,pp);

    E1PSH("imsl_f_spline_value","imsl_d_spline_value");

    lv_value = imsl_amach(6);
    x = (Mfloat)tempx;

    IMSL_CALL(l_spline_value(x,pp,argptr));
    va_end(argptr);

    E1POP("imsl_f_spline_value","imsl_d_spline_value");
    return (lv_value);

}

#ifdef ANSI
static VA_LIST_HACK l_spline_value(Mfloat tempx, Mf_spline *pp, va_list argptr)
#else
static VA_LIST_HACK l_spline_value(tempx,pp,argptr)
    Mfloat       tempx;
    Mf_spline   *pp;
    va_list      argptr;
#endif
{
    Mint            arg_number = 2;
    Mint            code;
    Mint            i;
    Mint            wants_pointer = 0;
    Mint            derivative=0;
    Mfloat          x;
    Mfloat          *value_ptr=NULL;
    Mfloat          *wk1 = NULL;
    Mfloat          *wk2 = NULL;
    Mfloat          *wk3 = NULL;
    Mfloat          *wk4 = NULL;
    Mfloat          *wk5 = NULL;
    Mint	    *iwk = NULL;

    Mint            length;
    Mfloat          *input_vector = NULL;
    Mfloat          **output_vector = NULL;
    Mfloat          *output_vector_user = NULL;
    Mint            grid = 0;
    Mint            grid_user = 0;

    x = (Mfloat)tempx;
    code = 1;
    while (code > 0) {
        code = va_arg(argptr, int);
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
    if (pp->order[0] < 1) {
        imsl_e1sti(1, pp->order[0]);

        imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_BAD_ORDER);
        goto RETURN;
    }
    if (!grid && !grid_user) {
	wk1  = (Mfloat *)imsl_malloc(pp->order[0]*sizeof(*wk1));
	wk2  = (Mfloat *)imsl_malloc(pp->order[0]*sizeof(*wk2));
	wk3  = (Mfloat *)imsl_malloc(pp->order[0]*sizeof(*wk3));

	if ((wk1 == NULL)||(wk2 == NULL)||(wk3 == NULL)) {
            imsl_e1stl(1, "order");
            imsl_e1sti(1, pp->order[0]);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
            goto FREE_SPACE; 
	}

        lv_value =imsl_b2der(&derivative,&x,pp->order,pp->knots[0],
		pp->num_coef,pp->coef[0],wk1,wk2,wk3);
        if (imsl_n1rty(1) > 3) lv_value = imsl_amach(6);

    }

    if (grid) {
	wk1 = (Mfloat *) imsl_malloc (pp->order[0]*(pp->num_coef[0]-pp->order[0]+1)
		*sizeof(*wk1));
	wk2 = (Mfloat *) imsl_malloc ((pp->num_coef[0]-pp->order[0]+2) 
		*sizeof(*wk2));
	wk3 = (Mfloat *) imsl_malloc (length * sizeof(*wk3));
	wk4 = (Mfloat *) imsl_malloc (length * sizeof(*wk4));
	wk5 = (Mfloat *) imsl_malloc ((pp->order[0]+3)*pp->order[0] 
		*sizeof(*wk5));
	iwk = (Mint *) imsl_malloc (length * sizeof(*iwk));
	*output_vector = (Mfloat *) imsl_malloc (length * sizeof(Mfloat));

	if ( (wk1 == NULL) || (wk2 == NULL) || (wk3 == NULL) || 
		(wk4 == NULL) || (wk5 == NULL) || (iwk == NULL) ||
		(*output_vector == NULL) ) {
            imsl_e1stl(1, "order");
            imsl_e1sti(1, pp->order[0]);
            imsl_e1stl(2, "length");
            imsl_e1sti(2, length);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
            goto FREE_SPACE;
        }

        imsl_b21gd(&derivative, &length, input_vector, 
		pp->order, pp->knots[0], pp->num_coef, pp->coef[0],
                *output_vector, wk1, wk2, iwk, wk3, wk4, wk5);
        if (imsl_n1rty(1) > 3) {
		for (i=0; i<length; i++)
 		    *(*(output_vector+i)) = imsl_amach(6);
	}
    } 

    if (grid_user) {
        wk1 = (Mfloat *) imsl_malloc (pp->order[0]*(pp->num_coef[0]-pp->order[0]+1)
                *sizeof(*wk1));
        wk2 = (Mfloat *) imsl_malloc ((pp->num_coef[0]-pp->order[0]+2)
                *sizeof(*wk2));
        wk3 = (Mfloat *) imsl_malloc (length * sizeof(*wk3));
        wk4 = (Mfloat *) imsl_malloc (length * sizeof(*wk4));
        wk5 = (Mfloat *) imsl_malloc ((pp->order[0]+3)*pp->order[0]
                *sizeof(*wk5));
        iwk = (Mint *) imsl_malloc (length * sizeof(*iwk));
 
        if ( (wk1 == NULL) || (wk2 == NULL) || (wk3 == NULL) ||
                (wk4 == NULL) || (wk5 == NULL) || (iwk == NULL) ) {
            imsl_e1stl(1, "order");
            imsl_e1sti(1, pp->order[0]);
            imsl_e1stl(2, "length");
            imsl_e1sti(2, length);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
            goto FREE_SPACE;
        }
 
        imsl_b21gd(&derivative, &length, input_vector,
                pp->order, pp->knots[0], pp->num_coef, pp->coef[0],
                output_vector_user, wk1, wk2, iwk, wk3, wk4, wk5);
        if (imsl_n1rty(1) > 3) {
                for (i=0; i<length; i++)
                    output_vector_user[i] = imsl_amach(6);
        }
    }

 
FREE_SPACE:
   if ( wk1 != NULL)         imsl_free(wk1);
   if ( wk2 != NULL)         imsl_free(wk2);
   if ( wk3 != NULL)         imsl_free(wk3);
   if ( wk4 != NULL)         imsl_free(wk4);
   if ( wk5 != NULL)         imsl_free(wk5);
   if ( iwk != NULL)         imsl_free(iwk);
RETURN:
    return (argptr);
}    
