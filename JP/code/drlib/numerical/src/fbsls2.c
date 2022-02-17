#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
static Mf_spline   *pp = NULL;
static VA_LIST_HACK  PROTO(l_spline_2d_least_squares,(Mint num_xdata,Mfloat xdata[],
                                           Mint num_ydata,Mfloat ydata[],
                                           Mfloat fdata[],Mint x_spline_space_dim,
                                           Mint y_spline_space_dim, va_list argptr));

#ifdef ANSI
Mf_spline  *imsl_f_spline_2d_least_squares(Mint num_xdata,Mfloat xdata[], Mint num_ydata,Mfloat ydata[],Mfloat fdata[],
                                               Mint x_spline_space_dim, Mint y_spline_space_dim, ...)
#else
Mf_spline  *imsl_f_spline_2d_least_squares(num_xdata,xdata,num_ydata,ydata,fdata,x_spline_space_dim,y_spline_space_dim,va_alist)
   Mint          num_xdata;
   Mfloat        xdata[];
   Mint          num_ydata;
   Mfloat        ydata[];
   Mfloat        fdata[];
   Mint          x_spline_space_dim;
   Mint          y_spline_space_dim;
   va_dcl
#endif
{
    va_list     argptr;
    VA_START(argptr,y_spline_space_dim );
#ifdef DOUBLE
   imsl_e1psh("imsl_d_spline_2d_least_squares");
#else
   imsl_e1psh("imsl_f_spline_2d_least_squares");
#endif
    pp = NULL;
    IMSL_CALL(l_spline_2d_least_squares(num_xdata,xdata,num_ydata,ydata,fdata,x_spline_space_dim,y_spline_space_dim,argptr));
    va_end(argptr);
#ifdef DOUBLE
   imsl_e1pop("imsl_d_spline_2d_least_squares");
#else
   imsl_e1pop("imsl_f_spline_2d_least_squares");
#endif
    return pp;
}
#ifdef ANSI
static VA_LIST_HACK   l_spline_2d_least_squares(Mint num_xdata,Mfloat xdata[],Mint num_ydata,Mfloat ydata[],
                                    Mfloat fdata[],Mint x_spline_space_dim,
                                           Mint y_spline_space_dim,va_list argptr)
#else
static VA_LIST_HACK   l_spline_2d_least_squares(num_xdata,xdata,num_ydata,ydata,fdata,x_spline_space_dim,y_spline_space_dim,argptr)
   Mint          num_xdata;
   Mfloat        xdata[];
   Mint          num_ydata;
   Mfloat        ydata[];
   Mfloat        fdata[];
   Mint          x_spline_space_dim;
   Mint          y_spline_space_dim;
   va_list       argptr;
#endif
{
   Mint          arg_number = 7;
   Mint          four = 4;
   Mint          code;
   Mint          users_xorder;
   Mint          users_yorder;
   Mfloat        *users_knots[2];
   Mfloat        *weights[2];
   Mint          users_fdata_col_dim;
   Mint          i;
   Mint          j;
   Mint          domain_dim;
   Mint          target_dim;
   Mint          num_coefs[2];
   Mint          order[2];
   Mint          order_given   = 0;
   Mint          knots_given   = 0;
   Mint          col_dim_given = 0;
   Mint          weights_given = 0;
   Mint          sse_wanted    = 0;
   Mint          derivative_x  = 0;
   Mint          derivative_y  = 0;
   Mint          free_the_structure = 0;
   Mfloat        *sse_return   = NULL;
   Mfloat        sse;
   Mint          col_dim_fdata;
   Mfloat        *wk           = NULL;
   Mint          temp_int1;
   Mint          temp_int2;
   Mfloat        temp_float;
   Mfloat        x_small;
   Mfloat        x_big;
   Mfloat        x_interval_size;
   Mfloat        y_small;
   Mfloat        y_big;
   Mfloat        y_interval_size;
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
            case IMSL_WEIGHTS:
                weights_given = 1;
                weights[0] = va_arg(argptr,Mfloat*);
                arg_number++;
                weights[1] = va_arg(argptr,Mfloat*);
                arg_number++;
                break;
            case IMSL_SSE:
                sse_wanted    = 1;
                sse_return    = va_arg(argptr,Mfloat*);
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
                                  /* SET PARAMETERS */
                                  /* DOMAIN & RANGE DIMENSIONS */
   domain_dim = 2;
   target_dim = 1;
                                  /* ORDERS */
   if (order_given == 0){
    order[0] = four;
    order[1] = four;
   }
   else{
    order[0] = users_xorder;
    order[1] = users_yorder;
   }
                                  /* NUM_COEFS */
   num_coefs[0] = x_spline_space_dim;
   num_coefs[1] = y_spline_space_dim; 
                                  /* FDATA COLUMN DIMENSION*/
   if (col_dim_given == 0){
     col_dim_fdata = num_ydata;
   }
   else{
     col_dim_fdata = users_fdata_col_dim;
   }

                                  /* TEST INPUT PARAMETERS */
                                  /* CHECK NXDATA */
        if (num_xdata < 1) {
                imsl_e1stl(1, "X");
                imsl_e1sti(1, num_xdata);
 
                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_POS_DATA_PTS);
        }
                                  /* CHECK KXORD */
        if (order[0] < 1) {
                imsl_e1stl(1, "X");
                imsl_e1sti(1, order[0]);

                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_ORDER_POSI);
        }
        if (imsl_n1rty(0) != 0)
                goto RETURN;
                                  /* CHECK NXCOEF */
        if (num_coefs[0] < order[0]) {
                imsl_e1stl(1, "X");
                imsl_e1sti(1, num_coefs[0]);
                imsl_e1sti(2, order[0]);

                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_COEFF_XY);
        }
                                  /* CHECK NXCOEF */
        if (num_coefs[0] > num_xdata) {
                imsl_e1stl(1, "X");
                imsl_e1sti(1, num_coefs[0]);
                imsl_e1sti(2, num_xdata);

                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_MORE_COEF_REQ);
        }
        if (imsl_n1rty(0) != 0)
                goto RETURN;
                                  /* CHECK NYDATA */
        if (num_ydata < 1) {
                imsl_e1stl(1, "Y");
                imsl_e1sti(1, num_ydata);
 
                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_POS_DATA_PTS);
        }
                                  /* CHECK KYORD */
        if (order[1] < 1) {
                imsl_e1stl(1, "Y");
                imsl_e1sti(1, order[1]);

                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_ORDER_POSI);
        }
        if (imsl_n1rty(0) != 0)
                goto RETURN;
                                  /* CHECK NYCOEF */
        if (num_coefs[1] < order[1]) {
                imsl_e1stl(1, "Y");
                imsl_e1sti(1, num_coefs[1]);
                imsl_e1sti(2, order[1]);

                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_COEFF_XY);
        }
                                  /* CHECK NYCOEF */
        if (num_coefs[1] > num_ydata) {
                imsl_e1stl(1, "Y");
                imsl_e1sti(1, num_coefs[1]);
                imsl_e1sti(2, num_ydata);

                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_MORE_COEF_REQ);
        }
        if (imsl_n1rty(0) != 0)
                goto RETURN;
                                  /* SET WEIGHTS IF NEEDED */
   if (weights_given == 0){
       weights[0] = (Mfloat *)imsl_malloc(num_xdata*sizeof(Mfloat));
       weights[1] = (Mfloat *)imsl_malloc(num_ydata*sizeof(Mfloat));
   if ((weights[0] == NULL)||(weights[1] == NULL)) {
            free_the_structure = 1;
            imsl_e1stl(1, "num_xdata");
            imsl_e1sti(1, num_xdata);
            imsl_e1stl(2, "num_ydata");
            imsl_e1sti(2, num_ydata);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_2);
            goto B2LS2_FREE_SPACE; }
       for (i=0;i<num_xdata;i++){
           *(weights[0]+i)= F_ONE;
       }
       for (i=0;i<num_ydata;i++){
           *(weights[1]+i)= F_ONE;
       }
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
            goto B2LS2_FREE_SPACE; 
   }
                                   /* COMPUTE THE KNOTS IF THE USER DID NOT
                                      SUPPLY THEM.  THESE KNOTS ARE EQUALLY
                                      SPACED INT HE INTERVAL AND STACKED
                                      APPRPRIATELY AT THE ENDPOINTS. */
  if (knots_given == 0){
        x_big = xdata[0];
        x_small = xdata[0];
        for ( i=1;i<num_xdata;i++){
            if (xdata[i] < x_small){
                 x_small = xdata[i];}
            else if (xdata[i] > x_big){
                 x_big   = xdata[i];}
        }
        x_interval_size=fabs(x_big-x_small);

        y_big = ydata[0];
        y_small = ydata[0];
        for ( i=1;i<num_ydata;i++){
            if (ydata[i] < y_small){
                 y_small = ydata[i];}
            else if (ydata[i] > y_big){
                 y_big   = ydata[i];}
        }
        y_interval_size=fabs(y_big-y_small);

        /* Set up X knot sequence. */
        for (i = 1; i <= (num_coefs[0] - order[0] + 2); i++) {
                pp->knots[0][i + order[0] - 2] = x_small + x_interval_size*((Mfloat) (i - 1) / (Mfloat) (num_coefs[0] -
                                                                order[0] + 1));
        }
        /* Stack knots. */
        for (i = 1; i <= (order[0] - 1); i++) {
                pp->knots[0][i - 1] = pp->knots[0][order[0] - 1];
                pp->knots[0][i + num_coefs[0]] =pp->knots[0][num_coefs[0]];
        }
        /* Set up Y knot sequence. */
        for (i = 1; i <= (num_coefs[1] - order[1] + 2); i++) {
                pp->knots[1][i + order[1] - 2] = y_small + y_interval_size * ((Mfloat) (i - 1) / (Mfloat) (num_coefs[1] -
                                                                order[1] + 1));
        }
        /* Stack knots. */
        for (i = 1; i <= (order[1] - 1); i++) {
                pp->knots[1][i - 1] = pp->knots[1][order[1] - 1];
                pp->knots[1][i + num_coefs[1]] = pp->knots[1][num_coefs[1]];
        }
  }
                                   /* GET THE WORKSPACE NEEDED IN B2LS2    */
   temp_int1 = imsl_i_max(order[0],order[1]);
   temp_int2 = 3*temp_int1 + (x_spline_space_dim +1)*num_ydata +order[0]*x_spline_space_dim + order[1]*y_spline_space_dim;
   wk  = (Mfloat *)imsl_malloc(temp_int2*sizeof(*wk));
   if (wk==NULL) {
            free_the_structure = 1;
            imsl_e1stl(1, "num_xdata");
            imsl_e1sti(1, num_xdata);
            imsl_e1stl(2, "num_ydata");
            imsl_e1sti(2, num_ydata);
            imsl_e1stl(3, "x_order");
            imsl_e1sti(3, order[0]);
            imsl_e1stl(4, "y_order");
            imsl_e1sti(4, order[1]);
            imsl_e1stl(5,"x_spline_space_dim");
            imsl_e1sti(5, x_spline_space_dim);
            imsl_e1stl(6,"y_spline_space_dim");
            imsl_e1sti(6, y_spline_space_dim);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_6);
            goto B2LS2_FREE_SPACE; }

  /* CHECK COL_DIM_FDATA */
        if (col_dim_fdata < num_ydata) {
                free_the_structure = 1;
                imsl_e1sti(1, col_dim_fdata);
                imsl_e1sti(2, num_ydata);
                imsl_ermes(IMSL_TERMINAL, IMSL_SPLINE_CD_FDATA);
                goto B2LS2_FREE_SPACE;
        }



                                   /* CALL THE SPLINE ROUTINE */ 
    imsl_f_m1ran(num_xdata, col_dim_fdata, fdata, fdata);
    imsl_b2ls2(&num_xdata,xdata,&num_ydata,ydata,fdata,&num_xdata,&(pp->order[0]),&(pp->order[1]),
           pp->knots[0],pp->knots[1],&(pp->num_coef[0]),&(pp->num_coef[1]),weights[0],weights[1],pp->coef[0],wk);
                                   /* IF THE SPLINE COULD NOT BE
                                      COMPUTED, THEN FREE THE STRUCTURE 
                                      AND RETURN NULL */   
   if (imsl_n1rty(1)>3) {
             free_the_structure = 1;
             imsl_f_m1ran(col_dim_fdata,num_xdata, fdata, fdata);
             goto B2LS2_FREE_SPACE;
   }
   else{
             imsl_f_m1ran(col_dim_fdata,num_xdata, fdata, fdata);
   }
                                   /* IF THE USER WANTS THE
                                      SUM OF THE SQUARES OF THE ERRORS AT THE 
                                      ORIGINAL DATA POINTS, COMPUTE IT*/
   if (sse_wanted == 1){
   sse = F_ZERO;
       for (i= 0;i<num_xdata;i++){
           for (j=0;j<num_ydata;j++){
               temp_float = imsl_b22dr(&derivative_x,&derivative_y,&(xdata[i]),&(ydata[j]),&(pp->order[0]),&(pp->order[1]),pp->knots[0],pp->knots[1],
                                 &(pp->num_coef[0]),&(pp->num_coef[1]),pp->coef[0],wk);
               sse += (temp_float-*(fdata+i*col_dim_fdata +j))*(temp_float-*(fdata+i*col_dim_fdata +j));

           }
       }
    *sse_return = (sse);
   }
                                    /*  FREE THE WORKSPACE USED */
B2LS2_FREE_SPACE:
   if ( wk != NULL)                                  imsl_free(wk);
   if (free_the_structure == 1) {
       if (pp != NULL)                               imsl_free(pp);
       pp = NULL;
   }
   if ((weights_given == 0)&&(weights[0] != NULL))   imsl_free(weights[0]);
   if ((weights_given == 0)&&(weights[1] != NULL))   imsl_free(weights[1]);
RETURN:
    return argptr;
}

