#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif

static      Mf_ppoly   *lv_pp;

#ifdef ANSI
Mf_ppoly *imsl_f_cub_spline_interp(Mint ndata, Mfloat xdata[],Mfloat fdata[])
#else
Mf_ppoly *imsl_f_cub_spline_interp(ndata,xdata,fdata)
   Mint          ndata;
   Mfloat        *xdata;
   Mfloat        *fdata;
#endif
{
   Mint          *orders;
   Mint          *num_breakpoints;
   Mint          four = 4;
   Mint          tmp;
   Mint          *iwork;
   Mfloat        *coef_work;
   Mint          i;
   Mint          domain_dim;
   Mint          target_dim;
   Mint          free_the_structure = 0;

   E1PSH("imsl_f_cub_spline_interp","imsl_d_cub_spline_interp");
  
   domain_dim = 1;
   target_dim = 1;
   orders = &four;
   tmp = ndata;
   num_breakpoints = &tmp;
                                   /* CREATE THE STRUCTURE */         
   lv_pp =   imsl_f_ppoly_create(domain_dim,target_dim,orders, num_breakpoints,0);
   if (imsl_n1rty(1)==4){
            imsl_e1mes(0,0," ");
            imsl_e1stl(1, "ndata");
            imsl_e1sti(1, ndata);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
       goto RETURN;
   }
                                   /* GET THE WORKSPACE    */
                                   /*NEED TO USE WORKSPACE FOR THE COEFFICIENTS
                                     SINCE THE SPACE NEEDED BY C2INT IS LARGER
                                     THAN THE SPACE PROVIDED IN THE STRUCTURE
                                     IMSL_PPOLY. THIS IS CAUSED BY C2INT USING
                                     THE LAST 4 MEMORY SPACES OF CSCOEF FOR 
                                     WORKSPACE  */
   coef_work = (Mfloat *)imsl_malloc(4*ndata*sizeof(*coef_work));
   iwork = (Mint *)imsl_malloc(ndata*sizeof(*iwork));
   if (iwork == NULL) {
            free_the_structure = 1;
            imsl_e1stl(1, "ndata");
            imsl_e1sti(1, ndata);
            imsl_ermes (IMSL_TERMINAL, IMSL_OUT_OF_MEMORY_1);
            goto FREE_SPACE; }

   if (imsl_n1rty(1)) goto RETURN;
                                   /* CALL THE SPLINE ROUTINE */ 
   imsl_c2int(&ndata,xdata,fdata,lv_pp->breakpoints[0],coef_work,iwork);
                                  /* COPY THE COEFFICIENTS INTO THE STRUCTURE */
      for (i=0;i<lv_pp->num_coef[0];i++)
       lv_pp->coef[0][i]=coef_work[i];
                                  /*  FREE THE WORKSPACE USED */
FREE_SPACE:
   if (free_the_structure == 1) imsl_free(lv_pp);
   if (iwork     != NULL) imsl_free(iwork);
   if (coef_work != NULL) imsl_free(coef_work);
RETURN:

   E1POP("imsl_f_cub_spline_interp","imsl_d_cub_spline_interp");
   return lv_pp;
}
