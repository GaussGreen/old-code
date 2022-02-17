/**********************************************************************
 *      Name: SrtGrfnMain.c                                           * 
 *  Function: Entry point to GRFN with raw data                       *
 * Copyright: (C) Paribas Capital Markets Ltd.                        *
 *--------------------------------------------------------------------*
 *    Author: Finbarr O'Sullivan                                      *
 *      Date: 18/10/95                                                *
 *--------------------------------------------------------------------*
 *    Inputs: Raw data from anywhere (Excel or 2020)                  *
 *   Returns:                                                         *
 *   Globals: Expects mkt and request list structures to exist        *
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 * 18/10/95 FOS     Created for SORT5-GRFN3 port to NT                *
 **********************************************************************/
#include "math.h"
#include "srt_h_all.h"
#include "SrtAccess.h"
#include "BGMEval.h"
#include "srt_h_allFx3F.h"

#define MAX_STP 3000

Err srt_f_set_GrfnCell(int tabRows, int tabCols, char ***tabStrings,
                   int **tabMask, GrfnCell ***tableau);
 

char *SrtGrfnMain(char *domestic,
	              int numParams, char **paramStrings, char **valueStrings, 
	              int numeventdates, long *eventdates,
                  long tableauRows, long tableauCols, char ***tableauStrings, int **tableauMask,
				  long auxWidth, long *auxLen, double **aux,
                  double *price, double *stdev, 
				  double **grfn_cells, double **pay_report)
{
Err            err = NULL;
int            status = 0;
SrtUndPtr      und;
SrtGrfnParam   grfnparam;
SrtIOStruct   *iolist;
GrfnCell     **tableau;
int            nUnderlyings = 0;
long          *evdates;

/* Set and Overwrite defaults with user defined parameters */
	if (err = srt_f_set_GrfnParams(numParams,paramStrings,
				valueStrings,&grfnparam))
	{
		return err;
	}

/* Gets the domestic underlying	*/
 	und = lookup_und(domestic); 
	if (!und)
	{
		return serror("Could not find % s underlying in market list", domestic);
	}

/*	New models	*/

	if (((SrtIrDesc*) (und->spec_desc))->mdl_type == DETERMINISTIC)
	{
		err = FirstValueDealInDeterministicIrModel (	numeventdates,
														eventdates,
														tableauRows,
														tableauCols,
														tableauStrings,
														tableauMask,
														auxWidth,
														auxLen,
														aux,
														domestic,
														&grfnparam,
														price);
		*stdev = 0.0;

		return err;	
	}
/*BGM*/
	if (((SrtIrDesc*) (und->spec_desc))->mdl_type == BGM)
	{
		err = FirstValueDealInBGMModel	(	numeventdates,
											eventdates,
											tableauRows,
											tableauCols,
											tableauStrings,
											tableauMask,
											auxWidth,
											auxLen,
											aux,
											domestic,
											&grfnparam,
											price,
											stdev);
		return err;
	}

/*	End of new models	*/

/* Construct GrfnCell structure from input string matrix */
	if (err = srt_f_set_GrfnCell(tableauRows,tableauCols,tableauStrings,
                            tableauMask,&tableau))
	{
		grfn_free_GrfnCellmatrix(tableau, tableauRows,tableauCols);
		return err;
	}
	
/* Gets the Main IO list */
	iolist = get_request_list();

/* Copies the input eventdates into a local array */
	evdates = (long *) srt_calloc( numeventdates, sizeof(long));
	memcpy (evdates, eventdates, numeventdates * sizeof(long));
	
/* Call to the Main grfn pricing function */
	err = srt_f_grfn(
				und,
				&grfnparam,
				numeventdates,
				&evdates,
				&tableauRows,
				&tableauCols,
				&tableau,
				0,
				NULL,
				auxWidth,
				auxLen,
				aux,
				iolist,
				grfn_cells, 
				pay_report);

/* Extract price from the "iolist" structure  */
	if (!err)
	{
	   srt_f_IOstructgetpremiumval(*iolist,price);
	   srt_f_IOstructgetstdevval(*iolist,stdev);
	}

/* Free tableau */
	grfn_free_GrfnCellmatrix(tableau, tableauRows,tableauCols);

/* Free the local copy of the event dates */
	srt_free(evdates);

	return err;
}

/* ============================================================= */
/*SrtGrfnMainExFrontier is the function called when the exercise frontier is optimized*/
/* ============================================================= */


char *SrtGrfnMainExFrontier(char *domestic,
	              int numParams, char **paramStrings, char **valueStrings, 
	              int numeventdates, long *eventdates,
                  long tableauRows, long tableauCols, char ***tableauStrings, int **tableauMask,
				  long auxWidth, long *auxLen, double **aux,
                  double *price, double *stdev,double *exfrontier, 
				  double **grfn_cells, double **pay_report)
{
Err            err = NULL;
int            status = 0;
SrtUndPtr      und;
SrtGrfnParam   grfnparam;
SrtIOStruct   *iolist;
GrfnCell     **tableau;
int            nUnderlyings = 0;
long          *evdates;

/* Set and Overwrite defaults with user defined parameters */
	if (err = srt_f_set_GrfnParams(numParams,paramStrings,
				valueStrings,&grfnparam))
	{
		return err;
	}

/* Gets the domestic underlying	*/
 	und = lookup_und(domestic); 
	if (!und)
	{
		return serror("Could not find % s underlying in market list", domestic);
	}

/*	New models	*/

	if (((SrtIrDesc*) (und->spec_desc))->mdl_type == DETERMINISTIC)
	{
		err = FirstValueDealInDeterministicIrModel (	numeventdates,
														eventdates,
														tableauRows,
														tableauCols,
														tableauStrings,
														tableauMask,
														auxWidth,
														auxLen,
														aux,
														domestic,
														&grfnparam,
														price);
		*stdev = 0.0;

		return err;	
	}

	/*BGM*/
	if (((SrtIrDesc*) (und->spec_desc))->mdl_type == BGM)
	{
		err = FirstValueBermudaDealOptimizationInBGMModel(	numeventdates,
															eventdates,
															tableauRows,
															tableauCols,
															tableauStrings,
															tableauMask,
															auxWidth,
															auxLen,
															aux,
															domestic,
															&grfnparam,
															price,
															stdev,
															exfrontier);
		return err;
	}

/*	End of new models	*/

/* Construct GrfnCell structure from input string matrix */
	if (err = srt_f_set_GrfnCell(tableauRows,tableauCols,tableauStrings,
                            tableauMask,&tableau))
	{
		grfn_free_GrfnCellmatrix(tableau, tableauRows,tableauCols);
		return err;
	}
	
/* Gets the Main IO list */
	iolist = get_request_list();

/* Copies the input eventdates into a local array */
	evdates = (long *) srt_calloc( numeventdates, sizeof(long));
	memcpy (evdates, eventdates, numeventdates * sizeof(long));
	
/* Call to the Main grfn pricing function */
	err = srt_f_grfn_ex_frontier(
				und,
				&grfnparam,
				numeventdates,
				&evdates,
				&tableauRows,
				&tableauCols,
				&tableau,
				0,
				NULL,
				auxWidth,
				auxLen,
				aux,
				iolist,
				grfn_cells, 
				pay_report,
				exfrontier);

/* Extract price from the "iolist" structure  */
	if (!err)
	{
	   srt_f_IOstructgetpremiumval(*iolist,price);
	   srt_f_IOstructgetstdevval(*iolist,stdev);
		
	}

/* Free tableau */
	grfn_free_GrfnCellmatrix(tableau, tableauRows,tableauCols);

/* Free the local copy of the event dates */
	srt_free(evdates);

	return err;
}
/* ============================================================= */


/* -----------------------------------------------------------
   From an array of Strings, construct an array of GrfnCell,
   with the right type and memory allocation for copied
   strings
   ------------------------------------------------------------ */
Err srt_f_set_GrfnCell(int tabRows, int tabCols, char ***tabStrings,
                   int **tabMask, GrfnCell ***tableau)
{
Err err = NULL;
GrfnCell **gc = NULL;
char tempStr[GRFN_DEF_ARGBUFSZ+1];
int i,j;

   if (!(gc = GrfnCellmatrix(tabRows,tabCols,0)))
   {
      return serror("Failed to allocate GrfnCell matrix");
   }

   for (i=0;i<tabRows;i++)
   {
      for (j=0;j<tabCols;j++)
	  {
	     switch(tabMask[i][j])
		 {
		    case GRFNBCELL:
			   (gc[i][j]).type = GRFNBCELL;
			   (gc[i][j]).str_alloced = SRT_NO;
			   break;
			case GRFNDCELL:
			   (gc[i][j]).type = GRFNDCELL;
			   (gc[i][j]).dval = strtod(tabStrings[i][j],NULL);
			   (gc[i][j]).str_alloced = SRT_NO;
			   break;
			case GRFNSCELL:
			   (gc[i][j]).type = GRFNSCELL;
			   strcpy(tempStr,tabStrings[i][j]);
			   (gc[i][j]).sval = new_string(tempStr);
			   (gc[i][j]).str_alloced = SRT_YES;
			   break;
			default:
			   break;
		 }
	  }
   }

   *tableau = gc;

   return NULL;

} /* END Err srt_f_SetGrfnCell(...) */

/* ============================================================= */



