#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

#ifdef ANSI
static VA_LIST_HACK l_vector_norm (Mint n, Mfloat *vector, va_list argptr);
#else
static VA_LIST_HACK l_vector_norm ();
#endif

static Mfloat lv_norm;
#ifdef ANSI
Mfloat     imsl_f_vector_norm (Mint n, Mfloat *vector, ...)
#else
Mfloat     imsl_f_vector_norm (n, vector, va_alist)
    Mint        n;
    Mfloat     *vector;
va_dcl
#endif
{
    va_list     argptr;
    VA_START (argptr, vector);
    E1PSH ("imsl_f_vector_norm", "imsl_d_vector_norm");
    lv_norm = -1.0;
    IMSL_CALL (l_vector_norm (n, vector, argptr));
    va_end (argptr);
    E1POP ("imsl_f_vector_norm", "imsl_d_vector_norm");
    return lv_norm;
}


#ifdef ANSI
static VA_LIST_HACK l_vector_norm (Mint n, Mfloat *vector, va_list argptr)
#else
static VA_LIST_HACK l_vector_norm (n, vector, argptr)
    Mint        n;
    Mfloat     *vector;
    va_list     argptr;
#endif
{
    Mfloat 	*second_vector = NULL;
    Mfloat	*diff_vector = NULL;
    Mint        code;
    Mint        arg_number = 2;
    Mint	one_norm = 0;
    Mint	two_norm = 1;
    Mint	two_vectors = 0;
    Mint	*index;
    Mint	i;

    code = 1;
    while (code > 0) {
	code = va_arg (argptr, Mint);
	arg_number++;
	switch (code) {
	case IMSL_ONE_NORM:
	    one_norm = 1;
	    two_norm = 0;
	    break;
	case IMSL_INF_NORM:
	    two_norm = 0;
	    index = va_arg(argptr, Mint *);
	    break;
	case IMSL_SECOND_VECTOR:
	    two_vectors = 1;
	    second_vector = va_arg(argptr, Mfloat *);
	    arg_number++;
	case 0:
	    break;
	default:
	    /* Argument number %(I2) is an unknown */
	    /* optional argument %(I1). */
	    imsl_e1sti (1, code);
	    imsl_e1sti (2, arg_number);
	    imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
	    break;
	}
    }

    if (n < 1) {
	imsl_e1sti (1, n);
	imsl_ermes (IMSL_TERMINAL, IMSL_LARGER_N_VALUE_NEEDED);
    }

    if (imsl_n1rty (0))
	goto RETURN;

    if (!two_vectors) {
    	if (one_norm) 
            lv_norm = imsl_sasum(n, vector, 1);
    	else if (two_norm)
	    lv_norm = imsl_snrm2(n, vector, 1);
   	else {
	    *index = imsl_isamax(n, vector, 1);
	    (*index)--;
	    lv_norm = fabs(vector[*index]);
	}
    }
    else {
	diff_vector = (Mfloat *) imsl_malloc(n*sizeof(*diff_vector));
	if (diff_vector == NULL) {
	    imsl_e1stl (1, "n");
	    imsl_e1sti (1, n);
	    imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
	    goto FREE_SPACE;
	}
	for (i=0; i<n; i++)
	    *(diff_vector+i) = *(vector+i) - *(second_vector+i);

        if (one_norm)
            lv_norm = imsl_sasum(n, diff_vector, 1);
        else if (two_norm)
            lv_norm = imsl_snrm2(n, diff_vector, 1);
        else {
            *index = imsl_isamax(n, diff_vector, 1);
	    (*index)--;
            lv_norm = fabs(diff_vector[*index]);
        }
    }
	
FREE_SPACE:
    if (diff_vector != NULL) imsl_free(diff_vector);
	
RETURN:
    if (imsl_n1rty (0) > 3) {
	lv_norm = imsl_amach(6);
    }
    return (argptr);
}
