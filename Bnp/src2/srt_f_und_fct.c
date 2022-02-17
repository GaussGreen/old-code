/* ===============================================================================

   FILENAME	:	srt_f_und_fct.c

   PURPOSE:     Functions to extract information from underlyings

   ===============================================================================
 */

#include "srt_h_all.h"
#include "srt_h_und_fct.h"

/* temporary definition */
#define SRT_OK OK

/* -------------------------------------------------------------------------------
            FUNCTIONS TO EXTRACT INFORMATION FROM THE SrtUndPtr == SrtUndDesc
   -------------------------------------------------------------------------------
 */

/* -------------------------------------------------------------------------------
 */
Err get_underlying_mdltype(SrtUndPtr und, SrtMdlType *mdl_type) {
  Err err = NULL;

  if (!und)
    return serror("Empty underlying in get_und_mdltype");

  switch (und->underl_type) {
  case INTEREST_RATE_UND:
    *mdl_type = get_mdltype_from_irund(und);
    break;
  case EQUITY_UND:
    *mdl_type = get_mdltype_from_eqund(und);
    break;
  case FOREX_UND:
    *mdl_type = get_mdltype_from_fxund(und);
    break;
  case BOND_UND:
  default:
    *mdl_type = NONE;
    break;
  }
  return NULL;

} /* END Err get_und_mdltype(...) */

/* -------------------------------------------------------------------------------
 */

Err get_underlying_mdldim(SrtUndPtr und, SrtMdlDim *mdl_dim) {
  Err err = NULL;

  if (!und)
    return serror("Empty underlying in get_und_mdldim");

  switch (und->underl_type) {
  case INTEREST_RATE_UND:
    *mdl_dim = get_mdldim_from_irund(und);
    break;
  case EQUITY_UND:
    *mdl_dim = get_mdldim_from_eqund(und);
    break;
  case FOREX_UND:
    *mdl_dim = get_mdldim_from_fxund(und);
    break;
  case BOND_UND:
  default:
    *mdl_dim = NO_FAC;
    break;
  }
  return NULL;

} /* END Err get_und_mdldim(..) */

/* -------------------------------------------------------------------------------
 */

Err get_underlying_ts(SrtUndPtr und, TermStruct **ts) {
  Err err = NULL;

  if (!und)
    return serror("Empty underlying in get_underlying_ts");

  switch (und->underl_type) {
  case INTEREST_RATE_UND:
    *ts = get_ts_from_irund(und);
    break;
  case EQUITY_UND:
    *ts = get_ts_from_eqund(und);
    break;
  case FOREX_UND:
    *ts = get_ts_from_fxund(und);
    break;
  case BOND_UND:
    *ts = get_ts_from_bndund(und);
    break;
  default:
    *ts = NULL;
    return serror("No termstructure in underlying");
    break;
  }

  if (*ts == NULL)
    return serror("No termstructure in underlying");

  return NULL;

} /* END Err get_underlying_ts(...) */

/* -------------------------------------------------------------------------------
 */
/* A function to free the TS attached to an underlying */
Err free_underlying_ts(SrtUndPtr und) {
  Err err = NULL;

  if (!und)
    return NULL;

  switch (und->underl_type) {
  case INTEREST_RATE_UND:
    err = srt_f_free_IRM_TermStruct(&(get_ts_from_irund(und)));
    break;
  case EQUITY_UND:
    err = srt_f_free_EQ_TermStruct(&(get_ts_from_eqund(und)));
    break;
  case FOREX_UND:
    err = srt_f_free_FX_TermStruct(&(get_ts_from_fxund(und)));
    break;
  case BOND_UND:
    break;
  default:
    break;
  }

  /* For safety reasons  , attaches a NULL Term Struct */
  set_irund_ts(und, NULL);

  return err;

} /* Err free_underlying_ts(...) */

/* ------------------------------------------------------------------------- */

/* Get the Discount Curve Name attached to the underlying */
Err get_underlying_discname(SrtUndPtr und, String *yc_name) {
  Err err = NULL;

  if (!und)
    return serror("Empty underlying in get_underlying_ycname");

  switch (und->underl_type) {
  case INTEREST_RATE_UND:
    (*yc_name) = get_ycname_from_irund(und);
    break;
  case EQUITY_UND:
    (*yc_name) = get_discname_from_eqund(und);
    break;
  case FOREX_UND:
    (*yc_name) = get_discname_from_fxund(und);
    break;
  case BOND_UND:
    (*yc_name) = get_discname_from_bndund(und);
    break;
  default:
    (*yc_name) = NULL;
    return serror("No Yield Curve attached to underlying");
    break;
  }

  if ((*yc_name) == NULL)
    return serror("No Yield Curve attached to underlying");

  return NULL;
} /* END Err get_underlying_discname(...) */

/* ---------------------------------------------------------------------------
 */

/* Get the Discount Curve Name attached to the underlying */
String get_discname_from_underlying(SrtUndPtr und) {
  char *yc_name = NULL;
  Err err = NULL;

  err = get_underlying_discname(und, &yc_name);
  if (err)
    return NULL;
  else
    return yc_name;
} /* END String get_discname_from_underlying(...) */

/* ------------------------------------------------------------------------- */

/* Get the Growth Curve Name attached to the underlying */
Err get_underlying_growname(SrtUndPtr und, String *curve_name) {
  Err err = NULL;

  if (!und)
    return serror("Empty underlying in get_underlying_ycname");

  switch (und->underl_type) {
  case INTEREST_RATE_UND:
    (*curve_name) = get_ycname_from_irund(und);
    break;
  case EQUITY_UND:
    (*curve_name) = get_dvdname_from_eqund(und);
    break;
  case FOREX_UND:
    (*curve_name) = get_forname_from_fxund(und);
    break;
  case BOND_UND:
    (*curve_name) = get_reponame_from_bndund(und);
    break;
  default:
    (*curve_name) = NULL;
    return serror("No Growth Curve attached to underlying");
    break;
  }

  if ((*curve_name) == NULL)
    return serror("No Growth Curve attached to underlying");

  return NULL;

} /* END Err get_underlying_growname(...) */

/* ---------------------------------------------------------------------------
 */

/* Get the Discount Curve Name attached to the underlying */
String get_dividend_name_from_underlying(SrtUndPtr und) {
  char *curve_name = NULL;
  Err err = NULL;

  err = get_underlying_growname(und, &curve_name);
  if (err)
    return NULL;
  else
    return curve_name;

} /* END String get_dividend_name_from_underlying(...) */

String get_repo_name_from_underlying(SrtUndPtr und) {
  char *curve_name = NULL;
  Err err = NULL;

  curve_name = get_reponame_from_eqund(und);
  if (err)
    return NULL;
  else
    return curve_name;

} /* END String get_dividend_name_from_underlying(...) */

/* ---------------------------------------------------------------------------
 */
SRT_Boolean is_model_Cheyette_type(SrtMdlType mdl_type) {
  if ((mdl_type == LGM) || (mdl_type == CHEY) || (mdl_type == CHEY_BETA) ||
      (mdl_type == MIXED_BETA) || (mdl_type == LGM_STOCH_VOL) ||
      (mdl_type == CHEY_STOCH_VOL) || (mdl_type == CHEY_BETA_STOCH_VOL))
    return SRT_TRUE;
  else
    return SRT_FALSE;

} /* END SRT_Boolean is_model_Cheyette_type(...) */

/* Checks whether an interest rate underlying fits into the Cheyette framework
 */
SRT_Boolean is_irund_Cheyette_type(SrtUndPtr und) {
  SrtMdlType mdl_type;

  if (!ISUNDTYPE(und, INTEREST_RATE_UND))
    return SRT_FALSE;

  mdl_type = get_mdltype_from_irund(und);
  return is_model_Cheyette_type(mdl_type);

} /* END SRT_Boolean is_irund_Cheyette_type(...) */

/* ---------------------------------------------------------------------------
 */
SRT_Boolean is_model_New_type(SrtMdlType mdl_type) {
  if ((mdl_type == NEWLGM) || (mdl_type == NEWCHEYBETA))
    return SRT_TRUE;
  else
    return SRT_FALSE;

} /* END SRT_Boolean is_model_New_type(...) */

/* Checks whether an interest rate underlying fits into the Cheyette framework
 */
SRT_Boolean is_irund_New_type(SrtUndPtr und) {
  SrtMdlType mdl_type;

  if (!ISUNDTYPE(und, INTEREST_RATE_UND))
    return SRT_FALSE;

  mdl_type = get_mdltype_from_irund(und);
  return is_model_New_type(mdl_type);

} /* END SRT_Boolean is_irund_New_type(...) */

/* ---------------------------------------------------------------------------
 */

/* -------------------------------------------------------------------------
   Functions to extract Yield Curve Information through the Underlying
   ------------------------------------------------------------------------- */

/* Get the calculation date == today from the underlying through the YC */
Date get_today_from_underlying(SrtUndPtr und) {
  Err err;
  String yc_name;
  String und_domrate_str;
  SrtMdlType mdl_type;
  SrtUndPtr und_dom;
  Date today;
  SrtCurvePtr yldcrv;

  err = get_underlying_mdltype(und, &mdl_type);
  if (err)
    return -1;

  if ((mdl_type == FX_STOCH_RATES) || (mdl_type == FX_LGMSV) ||
      (mdl_type == FX_BETADLM)) {
    /* Store the associated curve in gd.domestic_und */
    und_domrate_str = get_domname_from_fxund(und);
    und_dom = lookup_und(und_domrate_str);
    if (!ISUNDTYPE(und_dom, INTEREST_RATE_UND))
      return -1;
    today = get_today_from_underlying(und_dom);
  } else if ((mdl_type == EQ_STOCH_RATES) ||
             (mdl_type == EQ_STOCH_RATES_SRVGS)) {
    /* Get the domestic underlying name and pointer to */
    und_domrate_str = get_discname_from_underlying(und);
    und_dom = lookup_und(und_domrate_str);
    if (!ISUNDTYPE(und_dom, INTEREST_RATE_UND))
      return -1;
    today = get_today_from_underlying(und_dom);
  } else {
    yc_name = get_discname_from_underlying(und);
    yldcrv = lookup_curve(yc_name);
    today = (Date)get_clcndate_from_curve(yldcrv);
  }

  return today;

} /* END Date get_today_from_underlying(...) */

/* ------------------------------------------------------------------------- */

/* Get the spot date from the underlying through the YC */
Date get_spotdate_from_underlying(SrtUndPtr und) {
  Err err;
  String yc_name;
  String und_domrate_str;
  SrtMdlType mdl_type;
  SrtUndPtr und_dom;
  Date spotdate;
  SrtCurvePtr yldcrv;

  err = get_underlying_mdltype(und, &mdl_type);
  if (err)
    return -1;

  if (mdl_type == FX_STOCH_RATES) {
    /* Store the associated curve in gd.domestic_und */
    und_domrate_str = get_domname_from_fxund(und);
    und_dom = lookup_und(und_domrate_str);
    if (!ISUNDTYPE(und_dom, INTEREST_RATE_UND))
      return -1;
    spotdate = get_spotdate_from_underlying(und_dom);
  } else if (mdl_type == EQ_STOCH_RATES) {
    /* Get the domestic underlying name and pointer to */
    und_domrate_str = get_discname_from_underlying(und);
    und_dom = lookup_und(und_domrate_str);
    if (!ISUNDTYPE(und_dom, INTEREST_RATE_UND))
      return -1;
    spotdate = get_spotdate_from_underlying(und_dom);
  } else {
    yc_name = get_discname_from_underlying(und);
    yldcrv = lookup_curve(yc_name);

    spotdate = (Date)get_spotdate_from_yldcrv(yldcrv);
  }
  return spotdate;

} /* END Date get_spotdate_from_und(...) */

/* ------------------------------------------------------------------------- */

/* Get the spot lag from the underlying through the YC */
int get_spotlag_from_underlying(SrtUndPtr und) {
  String yc_name;
  int spot_lag;
  String ccy_str;
  SrtCcyParam *ccy_param;
  SrtCurvePtr yldcrv;
  Err err;

  yc_name = get_discname_from_underlying(und);
  yldcrv = lookup_curve(yc_name);

  ccy_param = get_ccyparam_from_yldcrv(yldcrv);
  if (!ccy_param) {
    ccy_str = get_underlying_ccy(und);
    err = swp_f_get_CcyParam_from_CcyStr(ccy_str, &ccy_param);
    if (err)
      return 2;
  }
  spot_lag = ccy_param->spot_lag;

  return spot_lag;

} /* END int get_spotlag_from_underlying(...) */

/* -----------------------------------------------------------------------------
 */

/* Just for FX underlyings  , gets the currencies: dom and for */

Err get_fx_underlying_currencies(SrtUndPtr fxund, String *dom_ccy,
                                 String *for_ccy) {
  String dom_name;
  String for_name;
  SrtCrvPtr crv;
  SrtUndPtr irund;
  SrtMdlType mdl_type;

  /* Gets the domestic and foreign curves names from the fx underlying */
  dom_name = get_domname_from_fxund(fxund);
  for_name = get_forname_from_fxund(fxund);

  /* According to the model type  , the names correspond to underlyings or
   * curves */
  mdl_type = get_mdltype_from_fxund(fxund);
  if (mdl_type == FX_STOCH_RATES) {
    irund = lookup_und(dom_name);
    if (!irund)
      return serror("Underlying %s not found for fx currencies", dom_name);
    *dom_ccy = get_underlying_ccy(irund);
    irund = lookup_und(for_name);
    if (!irund)
      return serror("Underlying %s not found for fx currencies", for_name);
    *for_ccy = get_underlying_ccy(irund);
  } else if (mdl_type == FX_LGMSV) {
    irund = lookup_und(dom_name);
    if (!irund)
      return serror("Underlying %s not found for fx currencies", dom_name);
    *dom_ccy = get_underlying_ccy(irund);
    irund = lookup_und(for_name);
    if (!irund)
      return serror("Underlying %s not found for fx currencies", for_name);
    *for_ccy = get_underlying_ccy(irund);
  } else if (mdl_type == FX_BETADLM) {
    irund = lookup_und(dom_name);
    if (!irund)
      return serror("Underlying %s not found for fx currencies", dom_name);
    *dom_ccy = get_underlying_ccy(irund);
    irund = lookup_und(for_name);
    if (!irund)
      return serror("Underlying %s not found for fx currencies", for_name);
    *for_ccy = get_underlying_ccy(irund);
  } else {
    crv = lookup_curve(dom_name);
    if (!crv)
      return serror("Curve %s not found for fx currencies", dom_name);
    *dom_ccy = get_curve_ccy(crv);
    crv = lookup_curve(for_name);
    if (!crv)
      return serror("Curve %s not found for fx currencies", for_name);
    *for_ccy = get_curve_ccy(crv);
  }

  /* Return a success message */
  return NULL;
}

/* ======================================================================== */
/*  -----------------------------------------------------------------
  FUNCNAME        :srt_f_undaddbndund
  AUTHOR          :O. Van Eyseren
  DESCRIPTION     :adds a bond to the SrtundList
  MODIFIES	  :SrtundPtr
  CALL            :

 ---------------------------------------------------------------------*/

/*
SrtErr srt_f_undaddbndund(	SrtundPtr	und  ,
                                String		bond_name  ,
                                String		und_lbl  ,
                                String 		mdl_lbl  ,
                                SwapDP		sdp  ,
                                double		coupon  ,
                                double		clean  ,
                                double		yield  ,
                                int		settle_days  ,
                                String		swap_yc_name  ,
                                String		repo_yc_name  ,
                                TermStruct	*ts  ,
                                SrtBndAuxParam  aux)
{
  SrtUndDesc	*und ;
  SrtBndDesc	*undbnd ;
  SrtErr	err  ;


/* create space for SrtUndDesc

  und = (SrtUndDesc *) srt_calloc(1  ,sizeof(SrtUndDesc)) ;

/* set up the elements in the underlying object val

  err = srt_f_interp_under(und_lbl  , &(und->underl_type));
  if (err)
                return err;
  if (und->underl_type != BOND_UND)
        return(serror("Trying to use undaddBNDund for a %s underl"  ,und_lbl));
  strcpy(und->underl_lbl  ,und_lbl) ;



/* creation of the bond object

  if ( und->underl_type == BOND_UND )
  {
        undbnd = (SrtBndDesc *) srt_calloc(1  ,sizeof(SrtBndDesc)) ;
        if (mdl_lbl!=0)
        {
                err =srt_f_interp_model(mdl_lbl  , &(undbnd->mdl_type)
,&(undbnd->mdl_dim) ); if (err) return err; strcpy(undbnd->mdl_lbl  ,mdl_lbl) ;
        }
        else
        {
                undbnd->mdl_type = NONE ;
                strcpy(undbnd->mdl_lbl  ,"NONE") ;
        }

        undbnd->p		= sdp		;
        undbnd->coupon		= coupon	;
        undbnd->settle_days	= settle_days	;
        undbnd->ts		= ts		;
        undbnd->yield		= yield		;
        undbnd->clean		= clean		;
        strcpy(undbnd->swap_name  ,swap_yc_name) ;
        strcpy(undbnd->repo_name  ,repo_yc_name) ;
        undbnd->aux_par		= aux		;
        und->und_info		= undbnd ;
  }
  if ((err=srt_f_lstins(und  ,bond_name  ,0.0  ,OBJ_PTR_UND  ,und))
         !=NULL)
    return (serror("Error in initialising %s object"  ,und->underl_lbl)) ;

  return NULL ;
}

*/

/* ========================================================================= */

char *get_discname_from_fxund(SrtUndPtr undptr) {
  SrtUndPtr dom_ptr;
  char *dom_name = get_domname_from_fxund(undptr);

  if (get_mdltype_from_fxund(undptr) == FX_STOCH_RATES) {
    dom_ptr = lookup_und(dom_name);
    if (!dom_ptr) {
      return NULL;
    }

    return get_ycname_from_irund(dom_ptr);
  } else {
    return dom_name;
  }
}
