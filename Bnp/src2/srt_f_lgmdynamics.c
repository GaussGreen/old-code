/* =================================================================================

   FILENAME: srt_f_lgmdynamics.c

   PURPOSE:  Give all the Equations required to descretise the LGM model:
				- drift of r
				- expectation of r at t+1 knowing t
				- local volatility of r
				- variance of r at t+1 knowing t (and same for phi)
   ================================================================================= */
   
#include "srt_h_all.h"
#include "srt_h_lgmdynamics.h"

/* --------------------------------------------------------------------------------
   FUNCNAME        :srt_f_lgmdrfatsam 
  
   DESCRIPTION     :calcs drift in the LGM model at stp,cur_sam (stp
			should correspond to cur_sam)
  
   MODIFIES        :drift_sam
   -------------------------------------------------------------------------------- */

void  srt_f_lgmdrfatsam
(
  SrtStpPtr   stp, 
  SrtSample  *cur_sam,
  SrtSample  *drift_sam,
  int         index 
)
{
  double mu;
  SrtIRMTmInf *tminf;

  if(!stp->next)return;

  tminf = (SrtIRMTmInf *)stp->tminf[index];

/* cur_sam and drift_sam are all fwd_sam's when this function is called, therefore
   they all should have index 0 in sam_get and samptr_get macros 
   K L Chau 6/6/95 */

  mu   =  -tminf->ev.onef.lambda * samptr_get(cur_sam,0,STATEVAR) 
							* stp->delta_t;
  
  samptr_get(drift_sam,0,STATEVAR) = samptr_get(cur_sam,0,STATEVAR) + mu;
  return;
}
/* ------------------------------------------------------------------------------- */
