#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
static Mf_spline   *pp = NULL;

static VA_LIST_HACK  PROTO(l_spline_2d_interp,(Mint num_xdata,Mfloat xdata[],
                                           Mint num_ydata,Mfloat ydata[],
                                           Mfloat fdata[],va_list argptr));
#ifdef ANSI
Mf_spline  *imsl_f_spline_2d_interp(Mint num_xdata,Mfloat xdata[], Mint num_ydata,Mfloat ydata[],Mfloat fdata[], ...)
#else
Mf_spline  *imsl_f_spline_2d_interp(num_xdata,xdata,num_ydata,ydata,fdata,va_alist)
   Mint          num_xdata;
   Mfloat        xdata[];
   Mint          num_ydata;
   Mfloat        ydata[];
   Mfloat        fdata[];
   va_dcl
#endif
{
    va_list     argptr;
    VA_START(argptr, fdata);
#ifdef DOUBLE
   imsl_e1psh("imsl_d_spline_2d_interp");
#else
   imsl_e1psh("imsl_f_spline_2d_interp");
#endif
    pp = NULL;
    IMSL_CALL(l_spline_2d_interp(num_xdata,xdata,num_ydata,ydata,fdata,argptr));
    va_end(argptr);

#ifdef DOUBLE
   imsl_e1pop("imsl_d_spline_2d_interp");
#else
   imsl_e1pop("imsl_f_spline_2d_interp");
#endif
    return pp;
}




#ifdef ANSI
static VA_LIST_HACK   l_spline_2d_interp(Mint num_xdata,Mfloat xdata[],Mint num_ydata,Mfloat ydata[],
                                    Mfloat fdata[],va_list argptr)
#else
static VA_LIST_HACK   l_spline_2d_interp(num_xdata,xdata,num_ydata,ydata,fdata,argptr)
   Mint          num_xdata;
   Mfloat        xdata[];
   Mint          num_ydata;
   Mfloat        ydata[];
   Mfloat        fdata[];
   va_list       argptr;
#endif
{
   Mint          arg_number        = 5;
   Mint          four              = 4;
   Mint          code;
   Mint          users_xorder;
   Mint          users_yorder;
   Mfloat        *users_knots[2];
   Mint          users_fdata_col_dim;
   Mint          col_dim_fdata;
   Mint          domain_dim;
   Mint          target_dim;
   Mint          num_coefs[2];
   Mint          order[2];
   Mint          order_given        = 0;
   Mint          knots_given        = 0;
   Mint          col_dim_given      = 0;
   Mint          free_the_structure = 0;
   Mfloat        *wk                = NULL;
   Mint          *iwk               = NULL;
   Mint          temp_int1;
   Mint          temp_int2;
   Mint          temp_int3;
   Mint          temp_int4;
   Mint          *iwk_b2nak         = NULL;
   Mfloat        *xsrt              = NULL;
    code = 1;
    while (code > 0) {
        code = va_arg(argptr, Mint);
        arg_number++;
        switch (code) {
            case IMSL_ORDER:
                order_given = 1;
                users_xorder = va_arg(argptr, Mint);
                arg_number++;
                users_yorder = va_arg(argptr, Mint);
                arg_number++;
                break;
            case IMSL_KNOTS:
                knots_given = 1;
                users_knots[0] = va_arg(argptr,Mfloat*);
                arg_number++;
                users_knots[1] = va_arg(argptr,Mfloat*);
                arg_number++;
                break;
            case IMSL_FDATA_COL_DIM:
                col_dim_given = 1;
                users_fdata_col_dim = va_arg(argptr, Mint);
                arg_number++;
            case 0:
                break;
            default:
                imsl_e1sti (1, code);
                imsl_e1sti (2, arg_number);
                imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
                goto RETURN;
        }
    } 
       
   domain_dim = 2;
   target_dim = 1;
   
   if (order_given == 0){
    order[0] = four;
    order[1] = four;
   }
   else{
    order[0] = users_xorder;
    order[1] = users_yorder;
   }

   num_coefs[0] = num_xdata;
   num_coefs[1] = num_ydata;

   if (col_dim_given == 0){
     col_dim_fdata = num_ydata;
   }
   else{
     col_dim_fdata = users_fdata_col_dim;
   }

                                   /* CREATE THE STRUCTURE */      
   if (knots_given == 0)   
      pp =   imsl_f_spline_create(domain_dim,target_dim,order, num_coefs,0);
   else
      pp =   imsl_f_spline_create(domain_dim,target_dim,order, num_coefs,IMSL_KNOTS,users_knots,0);
   if (imsl_n1rty(1)== 4){
            imsl_e1mes(0,0," ");
            imsl_e1stl(1, "num_xdata");
            imsl_e1sti(1, num_xdata);
            imsl_e1stl(2, "num_ydata");
            imsl_e1sti(2, num_ydata);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
       goto RETURN;
   }
                                   /* CHECK X arguments */
        if (pp->order[0] < 1) {
                imsl_e1sti(1, pp->order[0]);
 
                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_ORDER_X);
        }
                                   /* CHECK KYORD */
        if (pp->order[1] < 1) {
                imsl_e1sti(1, pp->order[1]);
 
                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_ORDER_Y);
        }
        if (imsl_n1rty(0) != 0){
	        free_the_structure = 1;
                goto B22IN_FREE_SPACE;
	}
                                   /* CHECK NXDATA */
        if (num_xdata < pp->order[0]) {
                imsl_e1sti(1, num_xdata);
                imsl_e1sti(2, pp->order[0]);
                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_DATA_X);
        }
                                   /* CHECK NYDATA */
        if (num_ydata <  pp->order[1]) {
                imsl_e1sti(1, num_ydata);
                imsl_e1sti(2, pp->order[1]);
                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_DATA_Y);
        }
        if (imsl_n1rty(0) != 0){
	        free_the_structure = 1;
                goto B22IN_FREE_SPACE;
	}

                                   /* COMPUTE THE KNOTS IF THE USER DID NOT
                                      SUPPLY THEM */
  if (knots_given == 0){
                                   /* GET THE WORKSPACE NEEDED IN B2NAK    */ 
   temp_int1 = imsl_i_max((num_xdata), (num_ydata));
   xsrt      = (Mfloat *)imsl_malloc(temp_int1*sizeof(*xsrt));
   iwk_b2nak = (Mint *)imsl_malloc(temp_int1*sizeof(*iwk_b2nak));

   if ((iwk_b2nak == NULL)||(xsrt == NULL)) {
            free_the_structure = 1;
            imsl_e1mes(0,0," ");
            imsl_e1stl(1, "num_xdata");
            imsl_e1sti(1,  num_xdata);
            imsl_e1stl(2, "num_ydata");
            imsl_e1sti(2,  num_ydata);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
            goto B2NAK_FREE_SPACE; }
   imsl_b2nak(&num_xdata,xdata,&(pp->order[0]),pp->knots[0],xsrt,iwk_b2nak);
   imsl_b2nak(&num_ydata,ydata,&(pp->order[1]),pp->knots[1],xsrt,iwk_b2nak);
   if (imsl_n1rty(1)) {
      free_the_structure = 1;
      goto B2NAK_FREE_SPACE;
  }
  }
                                   /* GET THE WORKSPACE NEEDED IN B22IN    */
   temp_int1 = imsl_i_max(((2*order[0]-1)*num_xdata),((2*order[1]-1)*num_ydata));
   temp_int2 = imsl_i_max(((3*order[0]-2)*num_xdata),((3*order[1]-2)*num_ydata));
   temp_int3 = imsl_i_max((num_xdata), (num_ydata));
   temp_int4 = (num_xdata*num_ydata) + temp_int1 + temp_int2 + 2*temp_int3;
   wk  = (Mfloat *)imsl_malloc(temp_int4*sizeof(*wk));
   iwk = (Mint *)imsl_malloc(temp_int3*sizeof(*iwk));


   if ((iwk == NULL)||(wk==NULL)) {
            free_the_structure = 1;
            imsl_e1stl(1, "num_xdata");
            imsl_e1sti(1, num_xdata);
            imsl_e1stl(2, "num_ydata");
            imsl_e1sti(2, num_ydata);
            imsl_e1stl(3, "x_order");
            imsl_e1sti(3, order[0]);
            imsl_e1stl(4, "y_order");
            imsl_e1sti(4, order[1]);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_4);
            goto B22IN_FREE_SPACE; }

  /* Check COL_DIM_FDATA */
        if (col_dim_fdata < num_ydata) {
                free_the_structure = 1;
                imsl_e1sti(1, col_dim_fdata);
                imsl_e1sti(2, num_ydata);

                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_CD_FDATA);
                goto B22IN_FREE_SPACE;
        }
                                   /* CALL THE SPLINE ROUTINE */ 
   
   imsl_f_m1ran(num_xdata, col_dim_fdata, fdata, fdata);
   imsl_b22in(&num_xdata,xdata,&num_ydata,ydata,fdata,&num_xdata,&(pp->order[0]),&(pp->order[1]),
           pp->knots[0],pp->knots[1],pp->coef[0],wk,iwk);
   imsl_f_m1ran(col_dim_fdata,num_xdata, fdata, fdata);
                                   /* IF THE SPLINE COULD NOT BE
                                      COMPUTED, THEN FREE THE STRUCTURE 
                                         AND RETURN NULL */ 
   if (imsl_n1rty(1) > 3) free_the_structure = 1;
                                    /*  FREE THE WORKSPACE USED */
B22IN_FREE_SPACE:
   if ( wk != NULL)          imsl_free(wk);
   if ( iwk != NULL)         imsl_free(iwk);
B2NAK_FREE_SPACE:
   if (free_the_structure == 1) {
       if (pp != NULL)       imsl_free(pp);
       pp = NULL;
   }
   if (knots_given == 0){
   if ( iwk_b2nak != NULL)   imsl_free(iwk_b2nak);
   if ( xsrt != NULL)        imsl_free(xsrt);
   }
RETURN:
    return argptr;
}







