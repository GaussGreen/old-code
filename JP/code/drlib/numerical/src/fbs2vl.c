#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
static VA_LIST_HACK PROTO(l_spline_2d_value,(Mfloat x, Mfloat y, Mf_spline *pp, va_list argptr));
static Mfloat  lv_value;

#ifdef ANSI
Mfloat imsl_f_spline_2d_value(Mfloat tempx, Mfloat tempy, Mf_spline *pp,...)
#else
Mfloat imsl_f_spline_2d_value(tempx, tempy,pp,va_alist)
    Mfloat    tempx;
    Mfloat    tempy;
    Mf_spline *pp;
    va_dcl
#endif
{
    Mfloat      x;
    Mfloat      y;
    va_list     argptr;

    VA_START(argptr,pp);

    E1PSH("imsl_f_spline_2d_value", "imsl_d_spline_2d_value");
    x = (Mfloat)tempx;
    y = (Mfloat)tempy;
    lv_value = imsl_amach(6);

    IMSL_CALL(l_spline_2d_value(x,y,pp,argptr));
    va_end(argptr);

    E1POP("imsl_f_spline_2d_value", "imsl_d_spline_2d_value");
    return (lv_value);

}

#ifdef ANSI
static VA_LIST_HACK l_spline_2d_value(Mfloat tempx, Mfloat tempy,Mf_spline *pp, va_list argptr)
#else
static VA_LIST_HACK l_spline_2d_value(tempx,tempy,pp,argptr)
    Mfloat       tempx;
    Mfloat       tempy;
    Mf_spline   *pp;
    va_list      argptr;
#endif
{
    Mint            arg_number = 3;
    Mint            code;
    Mint            i;
    Mint            temp_int;
    Mint            wants_pointer = 0;
    Mint            derivative_x  =0;
    Mint            derivative_y  =0;
    Mfloat          x;
    Mfloat          y;
    Mfloat          *value_ptr= NULL;
    Mfloat          *wk1      = NULL;

    Mint	nx;
    Mint	ny;
    Mfloat	*xvec = NULL;
    Mfloat	*yvec = NULL;
    Mint	grid = 0;
    Mint	grid_user = 0;
    Mfloat	**output_vector = NULL;
    Mfloat	*output_vector_user = NULL;

    Mint	*leftx = NULL;
    Mint	*lefty = NULL;
    Mfloat	*a = NULL;
    Mfloat	*b = NULL;
    Mfloat	*dbiatx = NULL;
    Mfloat	*dbiaty = NULL;
    Mfloat	*bx = NULL;
    Mfloat	*by = NULL;
    Mint	kxord;
    Mint	kyord;

    x = (Mfloat)tempx;
    y = (Mfloat)tempy;
    code = 1;
    while (code > 0) {
        code = va_arg(argptr, int);
        arg_number++;
        switch (code) {
            case IMSL_DERIV:
                derivative_x = va_arg(argptr, Mint);
                arg_number++;
                derivative_y = va_arg(argptr, Mint);
                arg_number++;
                break;
            case IMSL_GRID:
                nx = va_arg(argptr, Mint);
                xvec = va_arg(argptr, Mfloat*);
                ny = va_arg(argptr, Mint);
                yvec = va_arg(argptr, Mfloat*);
                output_vector = va_arg(argptr, Mfloat**);
                grid = 1;
                arg_number += 3;
                break;
            case IMSL_GRID_USER:
                nx = va_arg(argptr, Mint);
                xvec = va_arg(argptr, Mfloat*);
                ny = va_arg(argptr, Mint);
                yvec = va_arg(argptr, Mfloat*);
                output_vector_user = va_arg(argptr, Mfloat*);
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

        /* Check X order */
        if (pp->order[0] < 1) {
                imsl_e1sti(1, pp->order[0]);

                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_ORDER_X);
        }
        /* Check Y order */
        if (pp->order[1] < 1) {
                imsl_e1sti(1, pp->order[1]);

                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_ORDER_Y);
        }
    if (imsl_n1rty(0) != 0) goto FREE_SPACE;

    if (!grid && !grid_user) {
        temp_int = 3*(imsl_i_max(pp->order[0],pp->order[1])) + pp->order[1];
        wk1  = (Mfloat *)imsl_malloc(temp_int*sizeof(*wk1));
        if (wk1 == NULL){
            imsl_e1stl(1, "x_order");
            imsl_e1sti(1, pp->order[0]);
            imsl_e1stl(2, "y_order");
            imsl_e1sti(2, pp->order[1]);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
            goto FREE_SPACE; 
	}
    }

    if (!grid && !grid_user) {
        lv_value =imsl_b22dr(&derivative_x,&derivative_y,&x,&y,
		&(pp->order[0]),&(pp->order[1]),pp->knots[0],pp->knots[1],
                &(pp->num_coef[0]),&(pp->num_coef[1]),pp->coef[0],wk1);
        if (imsl_n1rty(1) > 3) lv_value = imsl_amach(6);

    }


    if (grid || grid_user) {
        kxord = pp->order[0];
        kyord = pp->order[1];
        
        leftx = (Mint *) imsl_malloc (nx * sizeof(*leftx));
        lefty = (Mint *) imsl_malloc (ny * sizeof(*lefty));
        a = (Mfloat *) imsl_malloc (kxord*kxord * sizeof(*a));
        b = (Mfloat *) imsl_malloc (kyord*kyord * sizeof(*b));
        dbiatx = (Mfloat *) imsl_malloc (kxord*(derivative_x+1) *
                sizeof(*dbiatx));
        dbiaty = (Mfloat *) imsl_malloc (kyord*(derivative_y+1) *
                sizeof(*dbiaty));
        bx = (Mfloat *) imsl_malloc (kxord*nx * sizeof(*bx));
        by = (Mfloat *) imsl_malloc (kyord*ny * sizeof(*by));

        if ( (leftx == NULL) || (lefty == NULL) ||
                (a == NULL) || (b == NULL) ||
                (dbiatx == NULL) || (dbiaty == NULL) ||
                (bx == NULL) || (by == NULL) ) {
            imsl_e1stl(1, "nx");
            imsl_e1sti(1, nx);
            imsl_e1stl(2, "ny");
            imsl_e1sti(2, ny);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
            goto FREE_SPACE;
        }
    }
    if (grid) {
	*output_vector = (Mfloat *) imsl_malloc (nx*ny*sizeof(Mfloat));

	imsl_b22gd (&derivative_x, &derivative_y, &nx, xvec, &ny, yvec,
		&kxord, &kyord, pp->knots[0], pp->knots[1], 
		&(pp->num_coef[0]), &(pp->num_coef[1]), pp->coef[0],
		*output_vector, &nx, leftx, lefty, a, b, dbiatx, dbiaty,
		bx, by);
        imsl_f_m1ran(ny, nx, *output_vector, *output_vector);
    }

    if (grid_user) {
        imsl_b22gd (&derivative_x, &derivative_y, &nx, xvec, &ny, yvec,
                &kxord, &kyord, pp->knots[0], pp->knots[1],
                &(pp->num_coef[0]), &(pp->num_coef[1]), pp->coef[0],
                output_vector_user, &nx, leftx, lefty, a, b, dbiatx, dbiaty,
                bx, by);
        imsl_f_m1ran(ny, nx, output_vector_user, output_vector_user);
    }


FREE_SPACE:
    if ( wk1 != NULL)           imsl_free(wk1);
    if (leftx != NULL) 		imsl_free(leftx);
    if (lefty != NULL) 		imsl_free(lefty);
    if (a != NULL)		imsl_free(a);
    if (b != NULL)		imsl_free(b);
    if (dbiatx != NULL)		imsl_free(dbiatx);
    if (dbiaty != NULL)		imsl_free(dbiaty);
    if (bx != NULL)		imsl_free(bx);
    if (by != NULL)		imsl_free(by);

RETURN:
    return (argptr);
}
