/* =============================================================================
   FILENAME:   swp_f_curve_list.cxx

   PURPOSE:    functions to store and handle curves (zc        , spreads , vol
   ,...)
   =============================================================================
 */

/* -------------------------------------------------------------------------------
                              Include Statements
   -------------------------------------------------------------------------------
 */

#include "srt_h_types.h"
#include "swp_h_all.h"
#include "swp_h_curve_list.h"
#include "swp_h_curve_struct.h"

/* -------------------------------------------------------------------------------
                            Static Declarations
   -------------------------------------------------------------------------------
 */

/* ------ The Static Pointer used to refer to the Curve List ---------------- */
static SrtCurveListPtr _srt_curves_list = NULL;

/* -------------------------------------------------------------------------------
                   Functions to operate on the curve list
   -------------------------------------------------------------------------------
 */

/* Allocate memory for the list where all the curves will be stored */
Err create_curve_list(String curve_list_name) {
  Err err;

  if ((err = srt_f_lstcreate(&_srt_curves_list, curve_list_name)) != NULL)
    return serror("Could not create curves list %s", curve_list_name);

  return NULL;
}

/* -------------------------------------------------------------------------- */

/* Destroy the list where all the curves are stored        , as well as the
 * curves */
Err destroy_all_curves() {
  Err err;

  err = srt_f_lstfree(_srt_curves_list, SRT_YES);
  if (err)
    return serror("%s in destroy_all_curves", err);

  srt_free(_srt_curves_list);

  return NULL;
}

/* -------------------------------------------------------------------------- */

/* An easy way to access the static pointer to the curves list */
SrtCurveListPtr get_curve_list(void) { return _srt_curves_list; }

/* -------------------------------------------------------------------------- */

/* --------------------------------------------------------------------------
                    FUNCTIONS TO OPERATE ON A SINGLE CURVE IN A LIST
   -------------------------------------------------------------------------- */

/* Check if a list has go a curve */
SRT_Boolean swp_f_iscurveinlist(SrtCurveListPtr curve_list, String curve_name) {
  SRT_Boolean isthere = SRT_NO;

  isthere = srt_f_lsthas(*curve_list, curve_name, 0);

  return (isthere);
} /* END SRT_Boolean srt_f_iscurveinlist(...) */

/* -------------------------------------------------------------------------- */

/* Get the Curve (structure) corresponding to the name in a list */
Err swp_f_getcurveinlist(SrtCurveListPtr curve_list, String curve_name,
                         SrtCurvePtr *curve) {
  SrtListAtom *tmp;
  int found = 0;

  /* Check the curve is in the list */
  if (swp_f_iscurveinlist(curve_list, curve_name) == SRT_NO)
    return serror("Curve %s not in list", curve_name);

  /* Go to the first atom in the linked list */
  tmp = curve_list->head;

  /* Loop throught the elements of the list until one with the same name is
   * found */
  while (tmp != NULL) {
    if (strcmp(tmp->element->name, curve_name) == 0) {
      found = 1;
      break;
    }
    tmp = tmp->next;
  }

  /* If the list does not contain an element with the same name return an error
   */
  if (found == 0)
    return serror("Curve %s not found", curve_name);

  /* The curve is the pval of the SrtListObject in the double linked list */
  (*curve) = (SrtCurveObj *)tmp->element->val.pval;

  /* Return a success message */
  return NULL;

} /* END Err swp_f_getcurveinlist(...) */

/* -------------------------------------------------------------------------- */

/* The function to pass to the srt_f_lstins function to free the SrtCurveDesc */
Err srt_f_curvevalfree(void *crvdesc) {
  SrtCurveDesc *crv_desc = (SrtCurveDesc *)(crvdesc);
  Err err = NULL;

  err = srt_f_curvedescfree(crvdesc);

  return err;
}
/* ------------------------------------------------------------------------------
 */
/* The same function        , but with a SrtCurveDesc as an input */
Err srt_f_curvedescfree(SrtCurveDesc *curve_desc) {
  Err err = NULL;

  if (!curve_desc)
    return NULL;

  /* Free the Currency parameters if a Yield Curve */
  switch (curve_desc->curve_type) {
  case YIELD_CURVE:

    if (((SrtYCDesc *)curve_desc->curve_desc)->ccy_param)
      srt_free(((SrtYCDesc *)curve_desc->curve_desc)->ccy_param);
    break;

  case DVD_CURVE:
    /* Delete a SrtFwdObject and everything it contains */
    if (((SrtFwdDesc *)curve_desc->curve_desc)->dvd_obj)
      err = srt_f_dvdobj_delete(
          &(((SrtFwdDesc *)curve_desc->curve_desc)->dvd_obj));
    if (err)
      return err;

    break;
  default:
    break;
  }

  /* Free the remainder of the curve_desc */
  srt_free(curve_desc->curve_desc);
  srt_free(curve_desc);

  return NULL;

} /* END Err srt_f_curvedescfree(...) */

/* ------------------------------------------------------------------------------
 */

/* ---------------------------------------------------------------------------------
                  THE MAIN FUNCTION USED TO ADD A CURVE  TO A LIST
   ---------------------------------------------------------------------------------
 */
/* Add a curve to a list (overwrite the previous curve if it exists ) */

Err swp_f_addcurvetolist(
    SrtCurveListPtr curve_list,
    String curve_name,  /* == curve_name in SrtCurveDesc */
    String curve_label, /* SPREAD        , FORWARD        , VOL        , YC */
    String type_label,  /* is this a SRT YC        , a MAP YC        , a WES ?*/
    String curve_ccy, long clcn_date, long spot_date,
    SrtCcyParam *ccy_param, /* Can be used to pass ccy_prm or ccy_str*/
    void *crv_object        /* Still used for Forward Curves for the moment */
) {
  SrtCurveDesc *curve;
  SrtYCDesc *yc;
  SrtFwdDesc *dvd;
  SrtFwdDesc *repo;
  SrtErr err;
  SrtYieldCurveType yc_type;

  /* Create space for SrtCurveDesc */
  curve = (SrtCurveDesc *)srt_calloc(1, sizeof(SrtCurveDesc));

  /* Checks the curve type */
  err = srt_f_interp_curve(curve_label, &(curve->curve_type));
  if (err)
    return err;
  strcpy(curve->curve_label, curve_label);

  /* Sets curve name in the SrtCurveDesc: uppercase it  */
  strcpy(curve->curve_name, curve_name);
  /* OVE Since Westminster is  case and space sensitive        , do not do this
          strupper(curve->curve_name);
          strip_white_space(curve->curve_name);
  */
  /* Sets curve currency in the SrtCurveDesc  */
  strcpy(curve->curve_ccy, curve_ccy);

  /* Sets the ticker to one at the moment */
  /* it will be changed when we insert the object in the list */
  curve->curve_ticker = 1;

  /* Set up the elements in the curve object val */
  switch (curve->curve_type) {
  case YIELD_CURVE:
    yc = (SrtYCDesc *)srt_calloc(1, sizeof(SrtYCDesc));
    if (type_label != 0) {
      if (err = srt_f_interp_yc(type_label, &yc_type))
        return err;
    }
    yc->yc_type = yc_type;
    yc->ccy_param = ccy_param;
    yc->clcn_date = clcn_date;
    yc->spot_date = spot_date;
    curve->curve_desc = yc;
    break;

  case DVD_CURVE:
    dvd = (SrtFwdDesc *)srt_calloc(1, sizeof(SrtFwdDesc));
    dvd->dvd_obj = (SrtDvdObj *)crv_object;
    dvd->clcn_date = clcn_date;
    curve->curve_desc = dvd;
    break;

  case REPO_CURVE:
    repo = (SrtFwdDesc *)srt_calloc(1, sizeof(SrtFwdDesc));
    repo->repo_obj = (SrtRepoObj *)crv_object;
    repo->clcn_date = clcn_date;
    curve->curve_desc = repo;
    break;

  case CMT_CURVE:
    break;

  default:
    break;
  }

  /* Insert (overwrite & update ticker) the object in the curve list (with the
   * same name) */
  err = srt_f_lstins(curve_list, curve_name, 0.0, OBJ_PTR_CURVE, curve,
                     &srt_f_curvevalfree, &(curve->curve_ticker));
  if (err)
    return (serror("Error in initialising %s object", curve->curve_label));

  /* Return a success message */
  return NULL;

} /* END Err swp_f_addcurvetolist(...) */

/* ---------------------------------------------------------------------------------
 */

/* ---------------------------------------------------------------------------------
                      THE MAIN FUNCTION USED TO ACCESS A CURVE
   ---------------------------------------------------------------------------------
 */

/* Get the Curve corresponding to the name in the curve list used */
SrtCurvePtr swp_f_lookup_curve(String curve_name) {
  Err err;
  String copy_name;
  int len;
  SrtCurveListPtr curve_list;
  SrtCurvePtr curveptr = NULL;

  /* Get the full curve list */
  curve_list = get_curve_list();

  /* Make a copy of the string not to modify it when UpperCasing it */
  len = strlen(curve_name);
  copy_name = (char *)srt_malloc(sizeof(char) * (len + 1));
  strncpy(copy_name, curve_name, len);
  copy_name[len] = '\0';

  /* Remove the ticker from the name (and leave it in its case and spacing) */
  /*	strupper(copy_name);
          strip_white_space(copy_name);
  */
  rem_tick_string(copy_name, copy_name);

  /* Get the curve object in the list ( if it is there) */
  err = swp_f_getcurveinlist(curve_list, copy_name, &curveptr);

  /* Free the copy name */
  srt_free(copy_name);

  /* return whatever has been found */
  if (err) {
    return NULL;
  } else {
    return curveptr;
  }

} /* END SrtCurvePtr swp_f_lookup_curve(...) */

/* ------------------------------------------------------------------------- */

/* Destroy the curve corresponding to the CurveName in the curve list */
Err swp_f_destroy_curve(String curve_name) {
  SrtCurvePtr crvptr;
  Err err = NULL;

  SrtCurveListPtr curve_list;

  /* Get THE curve list */
  curve_list = get_curve_list();

  /* Get the curve in the list */
  crvptr = swp_f_lookup_curve(curve_name);

  /* If the curve does not exist: nothing to do */
  if (!crvptr) {
    return NULL;
  }

  /* Removes and frees the object from the list */
  err = srt_f_lstremobj(curve_list, curve_name, 0.0);
  if (err)
    return err;

  /* Return a success message */
  return NULL;

} /* END Err swp_f_destroy_curve (...) */

/* ========================================================================== */
