/* ================================================================================
   
   FILENAME:      srt_f_chephiappx.c
	 
   PURPOSE:       Test a few approximations on the Phi variable for Cheyette

   ================================================================================ */

/* --------------------------------------------------------------------------------
   FUNCNAME        :srt_f_chelinrphi
   AUTHOR          :E.Auld
   MODIFIES        :sam->phi
   DESCRIPTION     :calculates the phi there would be at the time in stp
if the short rate had followed a linear path from today,
and if vol and lambda were constant at their initial values.
in this case, it is possible to solve the differential equation:
 
dphi(t)/dt = (sig^2 * [r(0) + (r(T)-r(0)t/T]^2 - 2*lambda*phi(t))dt
phi(0)= 0
   
   CALL            :

<%%END
  -------------------------------------------------------------------------------- */

#include "srt_h_all.h"
#include "math.h"

/*<%%STADEC*/
void   srt_f_chelinrphi
  (
  SrtStpPtr top, 
  SrtStpPtr stp, 
  SrtSample *sam      /* answer returned here */
  ) 
/*<%%ENDDEC*/
{
  double a,b,c;
  double sig,lam,T,r0,rT,onebyT,onebylam,dr,sig2;
  SrtIRMTmInf *tminf,*toptminf;

  toptminf = top->tminf[0];
  tminf = stp->tminf[0];

  sig  = toptminf->ev.onef.sig;
  lam  = toptminf->ev.onef.lambda;
  T  = stp->time;
  if(T<1.0e-6){samptr_get(sam,0,PHI)=0.0;return;}
  r0  = sam_get(toptminf->fwd_sam,0,SHORT_RATE);
  rT  = samptr_get(sam,0,SHORT_RATE);
  onebyT  = 1.0/T;
  onebylam = 1.0/lam;
  dr  = rT-r0;
  sig2  = sig*sig;

  a  = onebyT*dr;
  a  *= .5*a*onebylam*sig2;
  b  = sig2*r0*dr*onebyT*onebylam - a*onebylam;
  c  = sig2*dr*0.5*onebyT*onebylam*onebylam;
  c  = sig2*r0*r0*0.5*onebylam + c*(dr*.5*onebyT*onebylam - r0);

  samptr_get(sam,0,PHI) = a*T*T + b*T + c*(1 - exp(-2.0*lam*T));
  return;
}
