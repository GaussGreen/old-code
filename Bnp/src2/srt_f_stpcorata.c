/* -------------------------------------------------------------------------

 	AUTHOR		: O Van Eyseren
	DATE		: Feb 23 1995
	FILE NAME	: srt_f_stpcorata.c
	PURPOSE		: attaches the right correlation/coefficients matrix to
				  the stps used (coeff for MC, correl for the tree)
   ------------------------------------------------------------------------- */

#include "srt_h_all.h"
#include "srt_h_stpcorata.h"

Err 	srt_f_attach_correl_to_stp(SrtStpPtr stp, SrtCorrLst *cls)
{
	Err err = NULL;
	SrtCorrLstVal *corrval;
	SrtLst *ls;
 
/* Makes a test on cls */
	if (!cls)
		return serror ("Passed a NULL SrtCorrLst in attach_correl_to_stp"); 
/* Sets stp to point to the first step (date == today) */
	stp = gototop(stp);

/* Sets ls to point to the first element in the list just after first step */
	ls = cls->head;
	while ((ls != NULL) 
			&& ( ((SrtCorrLstVal *)ls->element->val.pval)->time < stp->time) )
		ls = ls->next;
	if (ls == NULL)
		return serror("No data in local correlation matrix after today = %d",
											stp->date);

/* Now corrval->time >= stp->time */
/* Have stp->coeff point to corrval->coeff: do not free SrtCorrLst yet !! */
	while (stp != NULL)
	{
		corrval = (SrtCorrLstVal *)ls->element->val.pval;
	/* If stp->time<=corrval->time  or if no more element in the 
        SrtCorrLst, points to the current corval*/
		if ((stp->time <= corrval->time) || (ls->next == NULL ))
 		{	
			stp->coeff	= corrval->coeff;		
			stp->correl	= corrval->correl;
			stp = stp->next;
		}
		else
	/* If stp->time>corrval->time and still elements points to the next corrval (if exists)*/
		{
			ls = ls->next;
		}	
	}

/* Return a success message: a NULL pointer */	
	return err;	                                                
}