#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef ANSI
void imsl_f_spline_print(Mf_spline *spline) 
#else
void imsl_f_spline_print(spline)  
  Mf_spline *spline;
#endif
{
   Mint          i;
   Mint          j;
   Mint          tot_coefs = 1;
   FILE          *nout;

   imsl_umach(2,&nout);
   imsl_umach(-3, &nout);

/*PRINT OUT THE CONTENTS OF THE STRUCTURE*/
   fprintf(nout,"\n\n PRINTOUT OF CONTENTS OF THE SPLINE STRUCTURE\n");
   fprintf(nout,"============================================ \n");
   fprintf(nout,"\n\n THE DOMAIN DIMENSION IS  %d \n",(spline->domain_dim));
   fprintf(nout," THE TARGET DIMENSION IS  %d \n\n",(spline->target_dim));
   for (i=0;i<spline->domain_dim;i++){
    fprintf(nout,"\n DOMAIN #%d \n",(i+1));
    fprintf(nout,"                 ORDER                   %d \n",(spline->order[i]));
    fprintf(nout,"                 NUMBER OF COEFFICIENTS  %d \n",(spline->num_coef[i]));
    fprintf(nout,"                 NUMBER OF KNOTS         %d \n",(spline->num_knots[i]));
    fprintf(nout,"                 THE KNOTS ARE:\n");
     for (j=0;j<(spline->num_knots[i]);j++)
    fprintf(nout,"                                       %4d.   %g\n",j+1, (spline->knots[i][j]));  
   }
    for (j=0;j<(spline->domain_dim);j++) tot_coefs *= spline->num_coef[j];
    fprintf(nout,"\n\n\n                 THE SETS OF COEFFICIENTS:\n\n");
     for (i=0;i<spline->target_dim;i++){
       fprintf(nout," *** TARGET  #%d ***\n\n",(i+1));
       fprintf(nout," COEF #       COEF \n");
       fprintf(nout," -----------------------\n");
       for (j=0;j<(tot_coefs);j++){
        fprintf(nout,"  %6d  %g \n",j+1, spline->coef[i][j]);
     }
     }
    
   fprintf(nout,"\n\n\n SOME ADDRESS INFORMATION \n\n");
   fprintf(nout," ADDRESS OF *spline            %d \n", (spline));
   fprintf(nout," ADDRESS OF spline->domain_dim %d \n",&(spline->domain_dim));
   fprintf(nout," ADDRESS OF spline->target_dim %d \n",&(spline->target_dim));
   fprintf(nout," ADDRESS OF spline->order      %d \n",&(spline->order));
   fprintf(nout," ADDRESS OF spline->num_coef   %d \n",&(spline->num_coef));
   fprintf(nout," ADDRESS OF spline->num_knots  %d \n",&(spline->num_knots));
   fprintf(nout," ADDRESS OF spline->knots      %d \n",&(spline->knots));
   fprintf(nout," ADDRESS OF spline->coef       %d \n\n",&(spline->coef));

   fprintf(nout," ADDRESS IN spline->order      %d \n",(spline->order));
   fprintf(nout," ADDRESS IN spline->num_coef   %d \n",(spline->num_coef));
   fprintf(nout," ADDRESS IN spline->num_knots  %d \n",(spline->num_knots));
   for (i=0;i<(spline->domain_dim);i++){
   fprintf(nout," ADDRESS IN spline->knots[%d]   %d \n",i,(spline->knots[i]));
   }
   for (i=0;i<(spline->target_dim);i++){
   fprintf(nout," ADDRESS IN spline->coef[%d]    %d \n",i,(spline->coef[i]));
   }
   return;

  }
