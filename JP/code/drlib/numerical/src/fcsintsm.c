#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
static Mf_ppoly   *pp = NULL;
 

static VA_LIST_HACK  PROTO(l_cub_spline_smooth,(Mint ndata,Mfloat xdata[],
                                           Mfloat fdata[],va_list argptr));

#ifdef ANSI
Mf_ppoly  *imsl_f_cub_spline_smooth(Mint ndata, Mfloat xdata[], 
                    Mfloat fdata[], ...)
#else
Mf_ppoly  *imsl_f_cub_spline_smooth(ndata,xdata,fdata,va_alist)
   Mint          ndata;
   Mfloat        xdata[];
   Mfloat        fdata[];
   va_dcl
#endif
{
    va_list     argptr;
    VA_START(argptr, fdata);
#ifdef DOUBLE
   imsl_e1psh("imsl_d_cub_spline_smooth");
#else
   imsl_e1psh("imsl_f_cub_spline_smooth");
#endif
    IMSL_CALL(l_cub_spline_smooth(ndata,xdata,fdata,argptr));
    va_end(argptr);
#ifdef DOUBLE
   imsl_e1pop("imsl_d_cub_spline_smooth");
#else
   imsl_e1pop("imsl_f_cub_spline_smooth");
#endif
    return pp;
}


#ifdef ANSI
static VA_LIST_HACK   l_cub_spline_smooth(Mint ndata, Mfloat xdata[],
                        Mfloat fdata[], va_list argptr)
#else
static VA_LIST_HACK   l_cub_spline_smooth(ndata,xdata,fdata,argptr)
   Mint          ndata;
   Mfloat        xdata[];
   Mfloat        fdata[];
   va_list       argptr;
#endif
{
   Mint          arg_number = 3;
   Mint          *orders;
   Mint          *num_breakpoints;
   Mint          four = 4;
   Mint          tmp;
   Mint          code;
   Mint          i;
   Mint          domain_dim;
   Mint          target_dim;
   Mfloat        eps;
   Mfloat        temp_float;
   Mfloat        small;
   Mfloat        big;
   Mint          num_weights_zero   = 0;
   Mfloat        *coef_work  = NULL;
   Mfloat        *weights    = NULL;
   Mfloat        *new_weights= NULL;
   Mint          free_the_structure = 0;
   Mfloat        users_sigma;
   Mfloat        *wk_cssmh    = NULL;
   Mint          *iwk_cssmh   = NULL;
   Mfloat        *wk_csscv    = NULL;
   Mfloat        *sdwk_csscv  = NULL;
   Mint          *ipvt_csscv  = NULL;
   Mint          use_cssmh   = 0;
   Mint          weights_given= 0;
   Mint          iequal      = 0;
    code = 1;
    while (code > 0) {
        code = va_arg(argptr, Mint);
        arg_number++;
        switch (code) {
            case IMSL_WEIGHTS:
                weights= va_arg(argptr, Mfloat*);
                arg_number++;
                weights_given = 1;
                break;
            case IMSL_SMOOTHING_PAR:
                users_sigma= (Mfloat) va_arg(argptr, Mdouble);
                arg_number++;
                use_cssmh = 1;
                break;
            case IMSL_SMOOTHING_PAR_ADR:
                users_sigma=  *(va_arg(argptr, Mfloat *));
                arg_number++;
                use_cssmh = 1;
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
       
   domain_dim = 1;
   target_dim = 1;
   orders = &four;
   tmp = ndata;
   num_breakpoints = &tmp;
                                   /* CREATE THE STRUCTURE */         
   pp =   imsl_f_ppoly_create(domain_dim,target_dim,orders, num_breakpoints,0);
   if (imsl_n1rty(1)==4){
            imsl_e1mes(0,0," ");
            imsl_e1stl(1, "ndata");
            imsl_e1sti(1, ndata);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
       goto RETURN;
	  }

              /* Check NDATA */
        if (use_cssmh ==1){
        if (ndata < 2) {
                imsl_e1sti(1, ndata);
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_AT_LEAST_2_PTS);
                }
        }
        else{
          if (ndata <= 2) {
                imsl_e1sti(1, ndata);
                imsl_ermes(IMSL_TERMINAL, IMSL_NEED_AT_LEAST_3_PTS);
                }
        }
   if (imsl_n1rty(0) != 0) goto RETURN;

   if (weights_given == 0) {
       weights = (Mfloat*)imsl_malloc((ndata)*sizeof(*weights));
       if (weights != NULL) sset(ndata,F_ONE,weights,1);
       if (weights == NULL) {
            free_the_structure = 1;
            imsl_e1stl(1, "ndata");
            imsl_e1sti(1, ndata);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
            goto WEIGHTS_FREE_SPACE; }
     }


   if (use_cssmh == 0){
                                   /*GET THE WORKSPACE    */
                                   /*NEED TO USE WORKSPACE FOR THE COEFFICIENTS
                                     SINCE THE SPACE NEEDED BY C2SCV IS LARGER
                                     THAN THE SPACE PROVIDED IN THE STRUCTURE
                                     IMSL_PPOLY. THIS IS CAUSED BY C2INT USING
                                     THE LAST 4 MEMORY SPACES OF CSCOEF FOR 
                                     WORKSPACE  */
   coef_work  = (Mfloat*)imsl_malloc(4*(ndata)*sizeof(*coef_work));
   wk_csscv   = (Mfloat*)imsl_malloc(7*(ndata+2)*sizeof(*wk_csscv));
   sdwk_csscv = (Mfloat*)imsl_malloc(2*ndata*sizeof(*sdwk_csscv));
   ipvt_csscv = (Mint *)imsl_malloc(ndata*sizeof(*ipvt_csscv));
   if ((coef_work == NULL)||(wk_csscv == NULL)||(sdwk_csscv == NULL)||(ipvt_csscv == NULL)) {
            free_the_structure = 1;
            imsl_e1stl(1, "ndata");
            imsl_e1sti(1, ndata);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
            goto C2SCV_FREE_SPACE; }

                                   /* CALL THE SPLINE ROUTINE */
   iequal = 0; /* THIS ARGUMENT IS NO LONGER USED IN IMSL_C2SCV, BUT IT IS STILL PASSED */
   imsl_c2scv(&ndata, xdata, fdata, &iequal, pp->breakpoints[0], coef_work,
              weights, wk_csscv, sdwk_csscv, ipvt_csscv);
   if (imsl_n1rty(1) > 3) {
        free_the_structure = 1;
        goto C2SCV_FREE_SPACE;
   }
                               /* COPY THE COEFFICIENTS INTO THE STRUCTURE */
      for (i=0;i<pp->num_coef[0];i++)
       pp->coef[0][i]=coef_work[i];
      
      goto  C2SCV_FREE_SPACE; 
   
   }
   else{
                                   /*GET THE WORKSPACE    */
                                   /*NEED TO USE WORKSPACE FOR THE COEFFICIENTS
                                     SINCE THE SPACE NEEDED BY C2SCV IS LARGER
                                     THAN THE SPACE PROVIDED IN THE STRUCTURE
                                     IMSL_PPOLY. THIS IS CAUSED BY C2INT USING
                                     THE LAST 4 MEMORY SPACES OF CSCOEF FOR 
                                     WORKSPACE  */
   if (weights_given == 1){
       new_weights = (Mfloat*)imsl_malloc((ndata)*sizeof(*weights));
   }
   coef_work  = (Mfloat*)imsl_malloc(4*(ndata)*sizeof(*coef_work));
   wk_cssmh   = (Mfloat*)imsl_malloc((8*ndata+5)*sizeof(*wk_cssmh));
   iwk_cssmh = (Mint *)imsl_malloc(ndata*sizeof(*iwk_cssmh));
   if ((coef_work == NULL)||(wk_cssmh == NULL)||(iwk_cssmh == NULL)|| ((weights_given==1)&&(new_weights == NULL)) ) {
            free_the_structure = 1;
            imsl_e1stl(1, "ndata");
            imsl_e1sti(1, ndata);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
            goto C2SMH_FREE_SPACE; }

                                   /* ADJUST THE WEIGHTS FOR USE IN CSSMH */
                                   /* WE MUST INVERT AND SQUARE THE INPUT 
                                      WEIGHTS, AND ALLOW FOR ZERO VALUED 
                                      WEIGHTS*/

   if (weights_given == 1){
    small = imsl_amach(1);
    big  = imsl_amach(2);
    if (big*small > F_ONE) {
	eps = 1000.0*small;
    }
    else{
	eps = 1000.0/big;
    }
    temp_float = F_ONE/sqrt(eps);
    for (i=0;i<ndata;i++){
        if ( (weights[i] == F_ZERO)){
             num_weights_zero++;
             new_weights[i] = temp_float;
	}
        else if (weights[i] > F_ZERO){
	     new_weights[i] = F_ONE/sqrt(weights[i]);
	}

	else {
                        imsl_e1sti(1, i);
                        imsl_e1str(1, weights[i]);
                        imsl_e1stl(1,"X");
                        imsl_ermes(IMSL_FATAL, IMSL_NEGATIVE_WEIGHTS);
                        free_the_structure = 1;
                        goto C2SMH_FREE_SPACE;
	}
    }

    if (num_weights_zero == ndata){

                        imsl_ermes(IMSL_FATAL, IMSL_SPLINE_NO_POS_ELMNT);
                        free_the_structure = 1;
                        goto C2SMH_FREE_SPACE;
    }
   }

    if (weights_given == 1){
    imsl_c2smh(&ndata, xdata, fdata, new_weights, &users_sigma, pp->breakpoints[0],
               coef_work, wk_cssmh, iwk_cssmh);
    }
    else{
    imsl_c2smh(&ndata, xdata, fdata, weights, &users_sigma, pp->breakpoints[0],
               coef_work, wk_cssmh, iwk_cssmh);
	
    }
   if (imsl_n1rty(1) > 3) {
        free_the_structure = 1;
        goto C2SMH_FREE_SPACE;
   }
                               /* COPY THE COEFFICIENTS INTO THE STRUCTURE */
      for (i=0;i<pp->num_coef[0];i++)
       pp->coef[0][i]=coef_work[i];
     goto  C2SMH_FREE_SPACE; 
       
   }
WEIGHTS_FREE_SPACE:  
   if ((free_the_structure == 1)&&(pp != NULL)){
              imsl_free(pp);
              pp = NULL;
              }
   if ((weights != NULL)&&(weights_given == 0))  imsl_free(weights);
   goto RETURN;
C2SCV_FREE_SPACE:
   if ((free_the_structure == 1)&&(pp != NULL)){
              imsl_free(pp);
              pp = NULL;
              }
   if (coef_work != NULL)     imsl_free(coef_work);
   if ( wk_csscv != NULL)     imsl_free(wk_csscv);
   if ( sdwk_csscv != NULL)   imsl_free(sdwk_csscv);
   if ( ipvt_csscv != NULL)   imsl_free(ipvt_csscv);
   if ((weights != NULL)&&(weights_given == 0))  imsl_free(weights);
   goto RETURN;
C2SMH_FREE_SPACE:
   if ((free_the_structure == 1)&&(pp != NULL)){
              imsl_free(pp);
              pp = NULL;
              }
   if (coef_work != NULL)     			     imsl_free(coef_work);
   if ( wk_cssmh != NULL)                            imsl_free(wk_cssmh);
   if ( iwk_cssmh != NULL)                           imsl_free(iwk_cssmh);
   if ((weights != NULL)&&(weights_given == 0))      imsl_free(weights);
   if ((new_weights != NULL)&&(weights_given == 1))  imsl_free(new_weights);
RETURN:
    return argptr;
     }
