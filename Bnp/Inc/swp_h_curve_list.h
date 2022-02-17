/* =========================================================================
   
   FILENAME:    swp_h_curve_list.h 

   PURPOSE:     Work with Curves in lists

   ========================================================================= */

#ifndef SWP_H_CURVE_LIST_H
#define SWP_H_CURVE_LIST_H

/* -------------------------------------------------------------------------
   A curve list is a usual double linked list
   ------------------------------------------------------------------------- */

typedef SrtList SrtCurveList, *SrtCurveListPtr ;

                                                  
/* -------------------------------------------------------------------------
   Functions to operate on the underlying list
   ------------------------------------------------------------------------- */

/* Allocate memory for the list where all the curves will be stored */
Err create_curve_list(String curve_list_name);

/* Destroy the list where all the curves are stored, as well as the curves */
Err destroy_all_curves();

/* An easy way to access the static pointer to the curves list */
SrtCurveListPtr get_curve_list(void);


/* -------------------------------------------------------------------------- 
              Functions to operate with a curve in a list
   -------------------------------------------------------------------------- */

/* Check if a list has go a curve */
SRT_Boolean swp_f_iscurvein_list(SrtCurveListPtr curve_list, String curve_name);

/* Get the Curve (structure) corresponding to the name in a list */
Err swp_f_getcurveinlist(SrtCurveListPtr curve_list,String curve_name, SrtCurvePtr *curve );

/* The function to pass to the srt_f_lstins function to free the SrtCurveDesc     */
Err srt_f_curvevalfree(void *crvdesc);

/* The same function, but with a SrtCurveDesc as an input */
Err srt_f_curvedescfree(SrtCurveDesc *crv_desc);

/* Add a curve to a list (overwrite the previous curve if it exists ) */
Err swp_f_addcurvetolist (
			SrtCurveListPtr  curve_list,
			String		     curve_name,       /* == curve_name in SrtCurveDesc */
 			String		     curve_label,      /* SPREAD, FORWARD, VOL, YC */
			String		     type_label,       /* is this a SRT YC, a MAP YC, a WES ?*/
			String		     curve_ccy,
			long             clcn_date,        /* Today == Calculation date */
			long             spot_date,        /* Spot Date */
			SrtCcyParam     *ccy_param,        /* Can be used to pass ccy_prm or ccy_str*/
			void 		    *crv_object       /* Still used for Forward Curves for the moment */
			);
			
/* ---------------------------------------------------------------------------------
                      FUNCTION TO LOOKUP A CURVE OR DESTROY IT                
   --------------------------------------------------------------------------------- */

/* Access a curve through its name */
SrtCurvePtr   swp_f_lookup_curve(String curve_name);      

#define   lookup_curve   swp_f_lookup_curve

/* Destroys a curve associated to a name */
Err swp_f_destroy_curve (String curve_name);


#endif

