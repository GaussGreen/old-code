/******************************************************************************\
*                Copyright (c) 1995 PARIBAS Capital Markets Group              *
********************************************************************************
*                                                                              
*   SYSTEM          :   GRF
*                                                                              
*   MODULE NAME     :   GRF_F_LANG_PVRNG                                        
*
*   PURPOSE         :
*
*   AUTHOR          :
*
*   DATE            :
*
*   VERSION         :
*
*   DESCRIPTION     :   
*
*   FUNCTIONS USED  :   
*
*   PARAMETERS      :
*
*   RETURNS         :
*
*   
********************************************************************************
*                           Amendment History                                  *
********************************************************************************
*
*   AMENDED BY      :
*
*   DATE            :
*
*   VERSION         :
*
*   REASON          :
*
*   REQUEST/BUG NO  :
*
*   DESCRIPTION
*
\******************************************************************************/
 


/******************************************************************************\
*                               Import Include Files                           *
\******************************************************************************/
 

#include "grf_h_all.h"
#include "math.h"


/******************************************************************************\
*                           Public Function Definitions                        *
\******************************************************************************/


/******************************************************************************\

  DESCRIPTION     : initialize top assuming that it is a call to 
 pvrng(auxrng1,auxrng2,startdate).  As this is a financial function, we will
end up with only one node. (The nodes after top will be freed.)
 This assumes that auxrng1 contains dates, auxrng2 contains (known) cash flows;
 only the cashflows in auxrng2 corresponding to dates in auxrng1
 that are greater than startdate are included; pvrng is supposed to calculate
 the pv of those cashflows; hence it is a sum of discount factors.  In
 grfn_eval_event() a node of type COMLL_PVRNG  will be evaluated just like a 
level payment; the only difference is that the coverages will be different.

\******************************************************************************/


Err grfn_pvrange_setup  (   String      name,  
                            GrfnDeal    *gd, 
                            COMLL_PTR   top,
                            int         index   )
{
  int	auxrng1,auxrng2;
  Date	startdate;   
  int	numdates,i,indexfirstdate;

/** all three arguments should be deterministic;
  we know from earlier error checking that there are three arguments **/

  auxrng1		= DTOL(top->next->next->dval);
  auxrng2		= DTOL(top->next->dval);
  startdate     = DTOL(top->dval);

/** check that startdate is not before today **/
  if(startdate  < gd->event_dates[gd->I])
    return serror("%s: %d",GRERR_BAD_START,startdate);

/** check that auxrng ids are ok **/
  CheckAuxRange2(auxrng1,auxrng2,gd);
    
/** check that numbers in auxrng1 are increasing **/
  for(i=0;i<gd->auxlen[auxrng1]-1;i++){
    if(DTOL(gd->aux[auxrng1][i]) >= DTOL(gd->aux[auxrng1][i+1]))
      return serror("%s:dates in aux[%d] must increase",name,auxrng1);
  }
/** find num dates to be valued **/
  for(i=0;i<gd->auxlen[auxrng1];i++){
    if(DTOL(gd->aux[auxrng1][i])>startdate)break;
  }
  if(i<gd->auxlen[auxrng1])numdates=gd->auxlen[auxrng1]-i;
  else return serror("%s:nothing left to value",name);
  indexfirstdate= gd->auxlen[auxrng1]-numdates;   
   
/** make sure dates are not completely insane **/
  if(DTOL(gd->aux[auxrng1][gd->auxlen[auxrng1]-1]) > 50000)
    return serror("%s:aux[%d][%d]:silly date %d",name,auxrng1,
      gd->auxlen[auxrng1]-1,gd->aux[auxrng1][gd->auxlen[auxrng1]-1]); 
  
/** allocate date and cvg arrays **/
  top->cvg  = srt_calloc(numdates,sizeof(double));
  top->dfind     = srt_calloc(numdates+1,sizeof(Date));
  top->dfindlen = numdates+1;

/** Fill underlying **/
  top->ivec[0] = index;

/** fill date and cvg arrays **/
  top->dfind[0] = startdate;
  for(i=0;i<numdates;i++)
  {
    top->cvg[i] = gd->aux[auxrng2][indexfirstdate+i];
  }
  for(i=0;i<numdates;i++)
  {
    top->dfind[i+1]    = DTOL(gd->aux[auxrng1][indexfirstdate+i]);
  }

  return NULL;
}



/******************************************************************************\

    FUNCTION        :   grfn_dirtyprice_setup 
  
    DESCRIPTION     :   This function has not yet been tested.  
                        For that reason, it is not yet included in the GRFN 
                        language. The reason one might want this is really to 
                        compute the CLEAN price of a bond; for which we need 
                        an accrued interest function in the GRFN function as 
                        well.

 Initializes top assuming that it is a call to 
 drty(auxrng1,auxrng2,startdate,enddate).  As this is a financial function, 
we will end up with only one node. (The nodes after top will be freed.)
 This assumes that auxrng1 contains dates, auxrng2 contains (known) cash flows;
 only the cashflows in auxrng2 corresponding to dates in auxrng1
 that are greater than startdate and <= enddate are included; drty is supposed
 to calculate the pv of those cashflows (like the dirty price of a bond); 
hence it is a sum of discount factors.  It would be calculated in
grfn_eval_event just like a call to pvrng(...) see above. 

\******************************************************************************/


Err grfn_dirtyprice_setup(  String      name, 
                            GrfnDeal    *gd, 
                            COMLL_PTR   top, 
                            int         index)
{
    int     auxrng1, auxrng2, numdates, i, j;
    Ddate   startdate,enddate;  


    /* All three arguments should be deterministic. We know from earlier 
       error checking that there are three arguments 
    */

  auxrng1   = DTOL(top->next->next->next->dval);
  auxrng2   = DTOL(top->next->next->dval);
  startdate     = (top->next->dval);
  enddate   = (top->dval);

/** check that startdate is not before today **/
  if(DTOL(floor(startdate))     < gd->event_dates[gd->I])
    return serror("%s: %d",GRERR_BAD_START,startdate);

/** check that auxrng ids are ok **/
  CheckAuxRange2(auxrng1,auxrng2,gd);
    
/** check that numbers in auxrng1 are increasing **/
  for(i=0;i<gd->auxlen[auxrng1]-1;i++){
    if(DTOL(gd->aux[auxrng1][i]) >= DTOL(gd->aux[auxrng1][i+1]))
      return serror("%s:dates in aux[%d] must increase",name,auxrng1);
  }
/** find index of first date to be valued **/
  for(i=0;i<gd->auxlen[auxrng1];i++){
    if(DTOL(gd->aux[auxrng1][i])>startdate)break;
  }
/** find index of first date not to be valued **/
  for(j=i;j<gd->auxlen[auxrng1];j++){
    if(DTOL(gd->aux[auxrng1][j])>enddate)break;
  }
  
  numdates=j-i;
  if(numdates == 0){
    top->type = COMLL_REAL;
    top->dval = 0.0;
    return NULL;
  }
   
/** make sure dates are not completely insane **/
  if(DTOL(gd->aux[auxrng1][gd->auxlen[auxrng1]-1]) > 50000)
    return serror("%s:aux[%d][%d]:silly date %d",name,auxrng1,
      gd->auxlen[auxrng1]-1,gd->aux[auxrng1][gd->auxlen[auxrng1]-1]); 
  
/** allocate date and cvg arrays **/
  top->cvg  = srt_calloc(numdates+1,sizeof(double));
  top->dfind     = srt_calloc(numdates+1,sizeof(Date));
  top->dfindlen = numdates+1;


/** fill date and cvg arrays **/
  top->cvg[0] = 0.0;
  top->dfind[0] = (long) startdate;
  for(j=1;j<=numdates;j++){
    top->cvg[j] = gd->aux[auxrng2][i+j-1];
    top->dfind[j]    = DTOL(gd->aux[auxrng1][i+j-1]);
  }

  return NULL;
}
