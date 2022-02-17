/* ==========================================================
   FILENAME:     srt_h_grfn.h

   PURPOSE:     the Main call to Grfn
                and a few top level functions
   ========================================================== */
#ifndef SRT_H_GRFN
#define SRT_H_GRFN

/*
#include "srt_h_all.h"
#include "grf_h_pubtypes.h"
*/

/* ---------------- A Few Types Useful for Communication with Grfn ---------------- */

typedef void	*SrtDealDesc;

typedef void   (*EvalEventDfsFct) (
				GrfnEvent  *,
				SrtSample  *,
				SrtUndInfo *) ;

typedef Err    (*EvalEventFct)(
				GrfnEvent *,
				SrtSample *,
                SrtDealDesc,
				double *,
				EvalEventDfsFct, 
				SrtUndInfo *, 
				double *);

/* ---------------------------------------------------------------------------------- */



/* ----------------------------------------------------------------------------------
                  THE TOP LEVEL GRFN FUNCTION
   ---------------------------------------------------------------------------------- */

Err srt_f_grfn    (   
		SrtUndPtr      und,      
		SrtGrfnParam  *grfnparam,                     
		long           numeventdates,             
		Date         **eventdates,
		long          *nrows,
		long          *ncols,
		GrfnCell    ***sprdsht,
		long           numgrng,
		GrfnRng       *grng,
		long           auxwidth,
		long          *auxlen,
		double       **aux,
		void         *answer,       /* This is actually a pointer to the IoList */
		double       **cellcontents,  /* The results of the last path when using MC */
		double       **knownpayments); /* A report of known CF [1] and pay dates [1] */

/* ------------------------------------------------------------------------------------ */
/* For initial tableau construction: from strings to  GrfnCells */
 Err SrtSetGrfnCell(int tabRows, int tabCols, char ***tabStrings,
                   int **tabMask, GrfnCell ***tableau);

 /* grf_f_lngsymtab.c */
String grfn_symbol_helpstr(int i);

long grfn_symbol_list_size();



/* ------------------------- For GRFN evaluation --------------------------- */
Err srt_f_grfn_val_GrfnDeal(
						GrfnDeal     *gd, 
						SrtStpPtr     sptr, 
						SrtUndPtr     und, 
						SrtGrfnParam *grfnparam,
						void         *answer,
						SrtUndInfo   *und_info);

Err srt_f_grfn_val_historical_events(
			SrtStpPtr     *step, 
			GrfnDeal      *gd,
			SrtUndInfo    *und_info);

/* --------------------------- For Launching Grfn from a flat file --------- */
/* srt_f_grfn_flatfile.c */
Err srt_f_grfn_from_flat_file(
				String 		    filename,
				SrtGrfnParam   *grfnparam,
				SrtUndPtr 	    und,
				double 		   *answer);

/*---------------------------- The srt_f_grfn function used when optimizing an 
								exercise frontier for Bermudans-------------- */


Err srt_f_grfn_ex_frontier    (   
		SrtUndPtr       und,      
		SrtGrfnParam   *grfnparam,                     
		long            numeventdates,             
		Date          **eventdates,
		long           *nrows,
		long           *ncols,
		GrfnCell     ***sprdsht,
		long            numgrng,
		GrfnRng        *grng,
		long            auxwidth,
		long           *auxlen,
		double        **aux,
		void           *answer,       /* This is actually a pointer to the IoList */
		double        **cellcontents,  /* The results of the last path when using MC */
		double        **knownpayments, /* A report of known CF [2], paid on pay dates [0], according to event dates [1] ([0..*nrows-1])*/
		double         *exfrontier);





#endif
