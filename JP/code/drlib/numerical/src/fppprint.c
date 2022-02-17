#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
#ifdef ANSI
void imsl_f_ppoly_print( Mf_ppoly *poly)
#else
void imsl_f_ppoly_print(poly)
  Mf_ppoly *poly;
#endif
{
   Mint          i;
   Mint          j;
   Mint          k;
   FILE          *nout;

   imsl_umach(2,&nout);
   imsl_umach(-3, &nout);

/*PRINT OUT THE CONTENTS OF THE STRUCTURE*/
   fprintf(nout,"\n\n PRINTOUT OF CONTENTS OF THE PPOLY STRUCTURE\n");
   fprintf(nout,"============================================ \n");
   fprintf(nout,"\n\n THE DOMAIN DIMENSION IS  %d \n",(poly->domain_dim));
   fprintf(nout," THE TARGET DIMENSION IS  %d \n\n",(poly->target_dim));
   for (i=0;i<poly->domain_dim;i++){
    fprintf(nout,"\n DOMAIN #%d \n",(i+1));
    fprintf(nout,"                 ORDER                   %d \n",(poly->order[i]));
    fprintf(nout,"                 NUMBER OF COEFFICIENTS  %d \n",(poly->num_coef[i]));
    fprintf(nout,"                 NUMBER OF BREAKPOINTS   %d \n",(poly->num_breakpoints[i]));
    fprintf(nout,"                 THE BREAKPOINTS ARE \n");
     for (j=0;j<(poly->num_breakpoints[i]);j++)
    fprintf(nout,"                                       %4d.   %g\n",j+1, (poly->breakpoints[i][j]));  
   }
    fprintf(nout,"\n\n\n                 THE SETS OF COEFFICIENTS:\n");
     for (i=0;i<poly->target_dim;i++){
    fprintf(nout," TARGET  #%d \n",(i+1));
     for (j=0;j<((poly->num_breakpoints[i])-1);j++){
      fprintf(nout,"      PP PIECE %4d",j+1);
      for (k=0;k<(poly->order[i]);k++)
#ifdef DOUBLE
       fprintf(nout,"  %30.12f ", poly->coef[i][j*(poly->order[i])+k]);
#else
       fprintf(nout,"  %14.6f ", poly->coef[i][j*(poly->order[i])+k]);
#endif
      fprintf(nout,"\n");
     }
     }
 

   fprintf(nout,"\n\n\n SOME ADDRESS INFORMATION \n\n");
   fprintf(nout," ADDRESS OF *poly                  %d \n", (poly));
   fprintf(nout," ADDRESS OF poly->domain_dim       %d \n",&(poly->domain_dim));
   fprintf(nout," ADDRESS OF poly->target_dim       %d \n",&(poly->target_dim));
   fprintf(nout," ADDRESS OF poly->order            %d \n",&(poly->order));
   fprintf(nout," ADDRESS OF poly->num_coef         %d \n",&(poly->num_coef));
   fprintf(nout," ADDRESS OF poly->num_breakpoints  %d \n",&(poly->num_breakpoints));
   fprintf(nout," ADDRESS OF poly->breakpoints      %d \n",&(poly->breakpoints));
   fprintf(nout," ADDRESS OF poly->coef             %d \n\n",&(poly->coef));

   fprintf(nout," ADDRESS IN poly->order            %d \n",(poly->order));
   fprintf(nout," ADDRESS IN poly->num_coef         %d \n",(poly->num_coef));
   fprintf(nout," ADDRESS IN poly->num_breakpoints  %d \n",(poly->num_breakpoints));
   for (i=0;i<(poly->domain_dim);i++){
   fprintf(nout," ADDRESS IN poly->breakpoints[%d]   %d \n",i,(poly->breakpoints[i]));
   }
   for (i=0;i<(poly->target_dim);i++){
   fprintf(nout," ADDRESS IN poly->coef[%d]          %d \n",i,(poly->coef[i]));
   }
   return;

  }
