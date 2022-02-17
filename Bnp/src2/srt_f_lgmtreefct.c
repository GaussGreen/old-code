/******************************************************************************/
/*                                                                            */
/*      SYSTEM:         SRT     SORT, Fixed Income 2020 Addins                */
/*      SUB_SYSTEM:     SWT     Swap Tools                                    */
/*                                                                            */
/*      MODULE NAME:    SRT_F_LGMTREEFCT.C                                    */
/*                                                                            */
/*      PURPOSE:        Functions to compute LGM tree 	                      */
/*                                                                            */
/*      AUTHORS:        E. Auld, K.L. Chau, J. Malhi, O. Van Eyseren          */
/*                                                                            */
/******************************************************************************/
/*                      Amendment History                                     */
/******************************************************************************/
/*                                                                            */
/*      AMENDED BY : Olivier Van Eyseren				      */
/*                                                                            */
/*      DATE:       August 1995                  		              */
/*                                                                            */
/*      REASON:     Compatibility with Two Factor Structures   		      */
/*                                                                            */
/******************************************************************************/
#include "math.h"
#include "srt_h_all.h"
/* #include "srt_h_lgmtreefct.h" */


/* ------------------------------------------------------------------------
   FUNCNAME        :lgm_calc_delta_x
   AUTHOR          :E. Auld
   DESCRIPTION     :calculates the variable step in an LGM trinomial tree
   MODIFIES	  :
   CALL            :
   INPUT		  : und_index refers to the underlying number in order to
		    collect the proper tminf (usually, this number is 0 :
		    first and only underlying)
  ------------------------------------------------------------------------- */

void LGM_calc_delta_x( SrtStpPtr top, double *delta_x, int und_index)
{
	double tmp,max_vol_dt=0,min_vol_dt,cum_var=0; 
	SrtStpPtr cur;
        SrtIRMTmInf *tminf;

	cur = top;

        tminf = (SrtIRMTmInf *)(cur->tminf[und_index]);
	min_vol_dt = tminf->ev.onef.sig2 * cur->delta_t;
	while(cur->next != NULL) {
		tminf = (SrtIRMTmInf *)cur->tminf[und_index];
		tmp = tminf->ev.onef.sig2 * cur->delta_t;
		cum_var += tmp;
		max_vol_dt = DMAX(max_vol_dt,tmp);
		min_vol_dt = DMIN(min_vol_dt,tmp);
        	cur=cur->next;
	} 
	tmp =  sqrt(cum_var/(double)cur->index)*DELTAOPTRATIO;
/*** the order of the following two tests is vital.  It was previously in the
opposite order and we considered that to be a bug E.Auld 17 Sep94 ***/
	max_vol_dt = sqrt(max_vol_dt);
	min_vol_dt = sqrt(min_vol_dt);
	
	if(tmp>DELTAMAXBOUND*min_vol_dt)
		tmp=DELTAMAXBOUND*min_vol_dt;
	if(tmp<DELTAMINBOUND*max_vol_dt)
		tmp=DELTAMINBOUND*max_vol_dt;

	*delta_x = tmp;
	
}

