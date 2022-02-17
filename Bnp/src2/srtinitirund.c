/* ==========================================================================
       FILENAME: SrtInitIRUnd.C                                        
     
       PURPOSE: Initialise an IR underlying and stores it in the underlying list 
    ========================================================================== */
/*        ========================================================================== 
  
	VISION Information Consulting
  
	Y2K Compliance
  
	Programmer   : Neil Glover
  
	Date         : 04/11/1998
  
	Fixes To     : SrtInitIrUnd
  
	Total Fixes  : 1
  
	===========================================================================*/
     
#include "srt_h_all.h"
#include "SrtAccess.h"
#include "grf_h_all.h" 
#ifdef PVMPI
#include "parallel.h"
#endif
    
char *SrtInitIRUnd(
		   char     *undName, 
		   char     *ycName, 
		   char     *model, 
		   int       volCrvRows,
		   int       volCrvCols,
		   double  **volCrvVals, 
		   int       tauCrvRows,
		   int       tauCrvCols,
		   double  **tauCrvVals,
		   double    beta, 
		   double    alpha,
		   double    gamma,
		   double    rho,
		   double    vovol,
		   double    etaOrBeta2,
		   double	 vasicek_init_cond,
		   int		 vasicek_mean_rev_level_rows,
		   int		 vasicek_mean_rev_level_cols,
		   double	 **vasicek_mean_rev_level_vals)
{
  Err				err;
  SrtUndListPtr   und_list;
  SrtMdlType		modelType;
  SrtMdlDim		modelDim; 
  TermStruct     *ts;
  long			today;
  SrtCurvePtr     crv;
  char           *ccy;
    	
#ifdef PVMPI
  if (!PRL_CONTEXT.bLocalComp) SendInitIRUndData(			
						 undName, 
						 ycName, 
						 model, 
						 volCrvRows,
						 volCrvCols,
						 volCrvVals, 
						 tauCrvRows,
						 tauCrvCols,
						 tauCrvVals,
						 beta, 
						 alpha,
						 gamma,
						 rho,
						 vovol,
						 etaOrBeta2);
#endif /*PVMPI*/	
  /* Get the Yield Curve */
  if ((crv = lookup_curve(ycName)) == NULL)
    {
      return serror("Fatal: (init_ir_und) Cannot find yield curve");
    }
  if (!ISCURVETYPE(crv,YIELD_CURVE))
    {
      return serror("Fatal: (init_ir_und) Must define a yield curve object");
    }
    
  /* Extract today and currency from this curve */
  today = get_today_from_curve(crv);
  ccy = get_curve_ccy(crv);   
    
  /* Interpret the model */
  if (err = srt_f_interp_model(model,&modelType,&modelDim))
    {
      return serror(err);
    }
    		
  /* Initialises the TermStruct of Volatility  */
  err = srt_f_init_IRM_TermStruct(
				  today,
				  volCrvVals,
				  volCrvCols,
				  volCrvRows,
				  tauCrvVals,
				  tauCrvCols,
				  tauCrvRows,
				  modelType, 
				  modelDim,
				  beta,
				  alpha,
				  gamma,
				  rho,
				  vovol,
				  etaOrBeta2,
				  vasicek_init_cond,
				  vasicek_mean_rev_level_rows,
				  vasicek_mean_rev_level_cols,
				  vasicek_mean_rev_level_vals,
				  &ts);
    	
  if (err)
    return err;
    
  /* Get the underlying list and check if not empty */
  und_list = get_underlying_list();
  if (und_list == NULL)
    return serror("No Underlying list defined: call SrtInit before");
    	
  /* Makes the Underlying Name UpperCase */
  strupper(undName);
    	
  /* Puts the Underlying in the Market List (updating its ticker) */
  err	= srt_f_addundtolist(
			     und_list,
			     undName,
			     "IR_UND",
			     ccy,
			     model,
			     ycName,
			     NULL,
				 NULL,
			     ts,
			     0.0);
  if (err)
    return serror("Fatal: (init_ir_und) Failed to add IR underlying");
    	
  /* Return a success message */
  return NULL;
    
} /* END Err SrtInitIRUnd(...) */
    

