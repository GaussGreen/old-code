/* ============================================================
   FILENAME:      srt_f_grfn.c

   PURPOSE:       the main fucntion to call grfn : srt_f_grfn
                  and a couple of others (grfn_help...)
   ============================================================ */
#include "GRF_H_ALL.H>
#include "math.h"
#include "srt_h_all.h"
#include "srt_h_correlation_list.h"
#include "srt_h_grfn.h"
#include "srt_h_grfn_undinfo.h"
#include "srt_h_sample.h"

#define shift 0.01
/*#define GAMMACONST  0.02*/

/* ---------------------------------------------------------------------------
 */

/* Free all the memory used by Grfn : events  , dearls  , specific corr and und
 */
void free_the_world(SrtStpPtr sptr, GrfnDeal *gd, SrtUndInfo *und_info) {
  /* Frees grfn events*/
  if (sptr)
    srt_f_dettach_grfn_events_from_steps(sptr);

  /* Frees grfn deal*/
  if (gd)
    grfn_free_inGrfnDeal(gd);

  /* Frees the SrtStpPtr */
  if (sptr)
    free_list(sptr);

  /* Frees the SrtCorrLst (with the corr matrixes inside) in und_info */
  if (und_info->corr_ts)
    srt_f_corrlstdelete(&und_info->corr_ts);

  /* Frees the underlying list */
  grfn_free_und_name();

} /* END static void free_the_world(...) */

/*---------------------------------------------------------------------------*/

/* The following two functions try to recreate what used to be in
   grfn_init_event.  It recreates the loop structure in grfn_attach_event.
   See logic in grfn_attach_event in file grf_f_get_event.c */

static Err attach_ytat(double today, double amdt, GrfnEvent *ev, SrtUndPtr und,
                       int und_index) {
  Ddate *df_ddates;
  int i;
  SrtMdlType mdl_type;
  Err err;

  err = get_underlying_mdltype(und, &mdl_type);
  if (err)
    return err;

  if (ev->dflen[und_index] > 0) {
    ev->yp[und_index] = srt_calloc(ev->dflen[und_index], sizeof(YTatt_param));

    /* Copy the swap/fra/... dates in a format that is used by Y_T_at_t_param */
    df_ddates = srt_calloc(ev->dflen[und_index], sizeof(double));
    for (i = 0; i < ev->dflen[und_index]; i++)
      df_ddates[i] = (Ddate)ev->dfd[und_index][i];

    /* The parameters that we will need to rebuild a zero coupon in the future*/
    if (ev->d + amdt > today) {
      if (mdl_type == VASICEK) {
        Vasicek_Y_T_at_param((Ddate)ev->d + amdt, df_ddates,
                             ev->dflen[und_index], und, ev->yp[und_index]);

      } else {
        Y_T_at_t_param((Ddate)ev->d + amdt, df_ddates, ev->dflen[und_index],
                       und, ev->yp[und_index]);
      }
    } else /* The parameters that we will need to rebuild a zero coupon for a
              PAST payment in the future (with today) */
    {
      if (mdl_type == VASICEK) {
        Vasicek_Y_T_at_param((Ddate)today, df_ddates, ev->dflen[und_index], und,
                             ev->yp[und_index]);

      } else
        Y_T_at_t_param((Ddate)today, df_ddates, ev->dflen[und_index], und,
                       ev->yp[und_index]);
    }

    srt_free(df_ddates);
  }
  return NULL;
} /* static Err attach_ytat(...) */

/* ---------------------------------------------------------------------------------
 */

/* Attaches to GrfnDeal information necessary to rebuild zero-coupon (for ALL
 * UND)*/
Err attach_params_to_GrfnEvent(
    GrfnDeal *gd, SrtStpPtr sptr, SrtUndPtr und,
    /*			SRT_Boolean       end_of_day_flg  , */
    int und_index) {
  double amdt = 0.0;
  SrtStpPtr stp = NULL;
  int i;
  double x;
  Err err;

  /* Go to the first date */
  if (sptr)
    stp = gototop(sptr);

  i = 0;
  /* Move on until a simulation date (today) is found
          while (i < gd->first_unkn_index)
                  i++ ;
  */

  /* Correct the date with the end-of-day flag */
  /*
          if ((end_of_day_flg==SRT_YES) && (stp!=NULL))
                  stp = stp->next ;
  */

  /* Loop through all the time steps */
  while (stp) {
    /* Check there is an event attached to this date */
    if (stp->e != NULL) {
      /* Increment the event date to look at if step is an event date */
      if (stp->date == gd->event_dates[i]) {
        i++;
        amdt = 0.0;
      } else if (i > 0 && GrfnIsStatus(gd->rowstatus[i - 1], GRFNCSAMERICAN)) {
        if (stp->date <= 0) {
          x = stp->ddate - gd->event_dates[i - 1];
          amdt = DMAX(x - (double)DTOL(x), 0.0);
        }
      }

      /* Attaches at the time step all the info needed to compute zc (Y_T) */
      err = attach_ytat((double)gd->today, amdt, (GrfnEvent *)stp->e, und,
                        und_index);
      if (err)
        return err;
    }

    stp = stp->next;
  }

  return NULL;

} /* END static Err attach_params_to_GrfnEvent(...) */

/* ----------------------------------------------------------------------------------
 */

/* Parses the Grfn tableau  , and collect all the names and info of the
 * underlyings */
Err get_deal_und_info(GrfnDeal *gd, SrtUndInfo *und_info) {
  char **und_name;
  int i;
  Err err;

  /* Allocate memory for the array where all the und wil be stored */
  und_name = svector_size(0, MAXUNDERLYING - 1, SRTBUFSZ);
  if (und_name == NULL)
    return serror("Memory allocation error");

  /* There should at least be one underlying returned from the call  , the
   * domestic one */
  und_info->no_of_underlyings = 1;
  strcpy(und_info->und_data[0].und_name, (const char *)(gd->domestic_und));

  /* Parse all the cells row by row to create a list of underlyings  ,
          in the same order as stored in the static _grfn_und_list
     (grf_f_underlying) */
  err = grfn_list_deal_underlyings(und_name, gd);
  if (err) {
    free_svector_size(und_name, 0, MAXUNDERLYING - 1, SRTBUFSZ);
    return err;
  }

  /* Stores the names of all the remaining underlyings in the SrtUndInfo */
  for (i = 1; i < MAXUNDERLYING && strcmp(und_name[i], ""); i++) {
    strcpy(und_info->und_data[i].und_name, und_name[i]);
    und_info->no_of_underlyings++;
  }

  free_svector_size(und_name, 0, MAXUNDERLYING - 1, SRTBUFSZ);
  und_name = NULL;

  /* Set the Numeraire index by default to the P&L underlying (indexed 0) */
  und_info->numeraire_index = 0;

  /* Add the different underlying necessary for STOCH RATES FX QUANTOS if they
     exist  , and modify the Numeraire index if required */
  err = add_more_underlyings_to_und_info(und_info);
  if (err)
    return err;

  /* Sets the underlying info accordingly for Grfn discretisation */
  err = analyze_und_info(und_info);
  if (err)
    return err;

  return NULL;

} /* END static Err get_deal_und_info(...) */

/* ----------------------------------------------------------------------------------
 */

/* Find index of first event date in the tableau considered not to be history */

int first_unkn_date_index(GrfnDeal *gd, SRT_Boolean end_of_day_flag) {
  int i;
  Date d;

  d = (Date)((long)gd->today + (long)(end_of_day_flag ? 1 : 0));

  for (i = 0; i < gd->sslength; i++) {
    if (gd->event_dates[i] >= d)
      break;
  }
  return i;
}

/* ---------------------------------------------------------------------------
 */

/* --------------------------------------------------------------------------------

   FUNCNAME        :  srt_f_grfn

   DESCRIPTION     :  the main grfn function. pv of deal is placed in *answer  ,
                     value of cells in last path are placed in *cellcontents
                     if cellcontents is not equal to null.

  AMENDMENTS    :
  Reference     :
  Author        :   M Nadjm
  Date          :   30 Nov 1994
  Description   :   (1) Removed syntax checking of user-defined variables
                        (This was not implemented);

                    (2) Moved validation of input arguments (error checking)
                        to module grf_f_types.c;

                    (3) Removed sptr from the GRFN deal structure;

                    (4) Created SrtStpPtr sptr to replace the above;

                    (5) grf_copy_to_GrfnDeal returns Err;

                    (6) Introduced grf_dettach_GrfnDeal

  Author        :   K L Chau
  Date          :   19 Jan 1995
  Description   :   Removed call to Y_T_at_t_param from grfn_init_event  ,
                    and relocate it to the main loop (grfn()) here.

  Author        :   M Nadjm
  Date          :   27 Mar 1995
  Description   :   Added new function to deal with tableaux with
                    (a) Calclulation date after the last event date;
                    (b) All events occuring in the past!


  -------------------------------------------------------------------------------
*/

/* THE function to price using Grfn.... */
Err srt_f_grfn(
    SrtUndPtr und, SrtGrfnParam *grfnparam, long numeventdates,
    Date **eventdates, long *nrows, long *ncols, GrfnCell ***sprdsht,
    long numgrng, GrfnRng *grng, long auxwidth, long *auxlen, double **aux,
    void *answer,           /* This is actually a pointer to the IoList */
    double **cellcontents,  /* The results of the last path when using MC */
    double **knownpayments) /* A report of known CF [2]  , paid on pay dates [0]
                               , according to event dates [1] ([0..*nrows-1])*/
{
  GrfnDeal gd;
  SrtStpPtr sptr = NULL;
  Err err = NULL;
  int i, j;
  int first_unknown_index;
  String domestic_und;
  String disc_curve;
  SrtUndInfo und_info;
  Date today; /* Calculation date */
  long init_nrows;
  long init_ncols;

  /* -------------- STEP 1: CHECK ON INPUTS AND INITIALISATIONS --------- */

  /* If imm_exit is required  , return NULL */
  if (grfnparam->imm_exit == SRT_YES) {
    return NULL;
  }

  /* Check the size of dates vs number of rows in tableau */
  if (numeventdates != *nrows)
    return serror("# of Dates and Events do not match in Grfn: %d  , %d",
                  numeventdates, *nrows);

  /* Get the calculation day */
  today = get_today_from_underlying(und);

  /* Initialisation of und_info */
  memset(&und_info, 0, sizeof(SrtUndInfo));

  /* Removes all empty rows   , and checks the dates are increasing */
  init_nrows = *nrows;
  init_ncols = *ncols;
  err = grfn_check_tableau_and_dates(nrows, *ncols, (*sprdsht), (*eventdates));
  if (err)
    return err;

  /* If calculation date is > last event date  , then extend tableau: create a
   * phantom event having a zero cash-flow. */
  if ((today > (*eventdates)[*nrows - 1]) ||
      ((grfnparam->end_of_day_flg == SRT_YES) &&
       (today == (*eventdates)[*nrows - 1]))) {
    /* Extend the event dates to add date after today */
    (*eventdates) = (Date *)realloc((*eventdates), (*nrows + 1) * sizeof(Date));

    /* Add tomorrow as a date */
    (*eventdates)[*nrows] = today + 1;

    /* Extend the tableau into a larger one with a zero on the last row (nrows
     * will be incremented) */
    err = grfn_extend_tableau(sprdsht, nrows, *ncols);
    if (err) {
      free_the_world(sptr, &gd, &und_info);
      return err;
    }
  }

  /* Make copies of all the inputs and place them inside the GrfnDeal */
  err = grfn_copy_to_GrfnDeal(&gd, *nrows, *ncols, (*sprdsht), (*eventdates),
                              numgrng, grng, auxwidth, auxlen, aux);
  if (err) {
    free_the_world(sptr, &gd, &und_info);
    return err;
  }

  /* Set the Historical Fixing function in the GrfnDeal */
  swp_f_GetFixingFunc(&(gd.fixing_fct));

  /* Initialise to 0 the payments (and their dates) of the GrfnDeal (cf PAY
   * function) */
  gd.pay_dates = NULL;
  gd.event_pay_dates = NULL;
  gd.pay_amounts = NULL;
  gd.numpay = 0;

  /* Set today (extracted from the yield curve via the underlying) in the
   * GrfnDeal */
  gd.today = today;

  /* Store the domestic underlying name in gd.domestic_und  , with the
   * associated curve */
  domestic_und = get_underlying_name(und);
  strcpy((char *)gd.domestic_und, domestic_und);
  disc_curve = get_discname_from_underlying(und);
  strcpy((char *)gd.disc_curve, disc_curve);

  /* Find index of first event date considered not to be history */
  gd.first_unkn_index = first_unknown_index =
      first_unkn_date_index(&gd, grfnparam->end_of_day_flg);
  /* Set up the end-Of_day flags (for payments/exercise... and fixings) */
  gd.end_of_day_index = ((grfnparam->end_of_day_flg == SRT_YES) ? 1 : 0);
  gd.end_of_day_fixing = ((grfnparam->end_of_day_fixings == SRT_YES) ? 1 : 0);
  gd.end_of_day_payment = ((grfnparam->end_of_day_payment == SRT_YES) ? 1 : 0);

  /* Output files when debugging */
  if (grfnparam->prin_ok) {
    gd.outfileptr = fopen("grf_sam.dat", "w");
  } else {
    gd.outfileptr = NULL;
  }

  /* ---------- STEP 2: FIRST PARSING OF TABLEAU TO COLLECT UND INFO
   * ----------*/

  /* Get info concerning :number of underlyings  , which ones.. (parses the
   * deal) */
  err = get_deal_und_info(&gd, &und_info);
  if (err) {
    free_the_world(sptr, &gd, &und_info);
    return err;
  }

  /*  When two underlyings or more: use MONTE CARLO so far */
  if (und_info.no_of_underlyings >= 3) {
    grfnparam->force_mc = SRT_YES;
  }
  if ((und_info.jumping == SRT_YES) && (grfnparam->jumping == SRT_YES)) {
    und_info.jumping =
        SRT_YES; /* If the model and the user agree then go for it */
  } else {
    und_info.jumping = SRT_NO; /* If one doesn't agree don't do it */
  }

  /*
   When two underlyings  , if they are lgm/lgm or lgm/bs  , one wil soon be able
   to can use tree as well and then the test will be if (
   und_info.no_of_underlyings > 2 )
  */

  /* ------------------ STEP 3: TIME DISCRETISATION ----------------------- */

  /* Build all the time steps  , for all the dates  , past and future (according
   * to event dates  , grfn params...) */
  err = srt_f_make_time_step(gd.sslength /*-first_unknown_index*/,
                             &gd.event_dates[0] /*first_unknown_index]*/, und,
                             &sptr, grfnparam, &und_info);
  if (err) {
    free_the_world(sptr, &gd, &und_info);
    return err;
  }

  /* ---------- STEP 4: SECOND PARSING OF TABLEAU TO BUILD EVENTS ---------- */

  /* Transforms the tableau strings into a Grfn command representation (COMLL)
   */
  err = srt_f_attach_grfn_events_to_steps(&gd, sptr);
  if (err) {
    free_the_world(sptr, &gd, &und_info);
    return err;
  }

  /* Loop through underlyings to attach relevant interest rate information: df
   * , spreads  , ...*/
  for (i = 0; i < und_info.no_of_underlyings; i++) {
    und = lookup_und(und_info.und_data[i].und_name);
    if (err) {
      free_the_world(sptr, &gd, &und_info);
      return err;
    }

    /* Attach to GrfnDeal  information necessary to rebuild zero-coupon (for ALL
     * UNDS) */
    err = attach_params_to_GrfnEvent(
        &gd, sptr, und,
        /*				grfnparam->end_of_day_flg  , */
        i);
    if (err != NULL) {
      free_the_world(sptr, &gd, &und_info);
      return err;
    }
  }

  /* ----------- STEP 5: CHOICE OF THE SOLVING DIRECTION : FWD/BKWD --------- */

  /* new test for lsm */

  if (grfnparam->lsm == SRT_FORBACK)
    grfnparam->imp_type = IMPBACKWARDLATTICE;
  else {
    /* Set the Grfn direction depending on the results of the parsing on the
     * full tableau */
    if ((GrfnIsStatus(gd.sumstatus, GRFNCSFUTUREREF)) &&
        GrfnIsStatus(gd.sumstatus, GRFNCSPASTREF)) {
      err = serror(GRERR_METHOD_CONFLICT);
      free_the_world(sptr, &gd, &und_info);
      return err;
    }

    if (GrfnIsStatus(gd.sumstatus, GRFNCSPASTREF)) {
      gd.method = GRFNMONTECARLO;
    } else if (GrfnIsStatus(gd.sumstatus, GRFNCSFUTUREREF)) {
      gd.method = GRFNBACKWARDLATTICE;
    } else
      gd.method = GRFNNOMETHOD;

    /*  Choose srt implementation method (this should maybe be done elsewhere)
     */

    if ((gd.method == GRFNMONTECARLO) ||
        (((int)grfnparam->force_mc) && (gd.method != GRFNBACKWARDLATTICE))) {
      grfnparam->imp_type = IMPMONTECARLO;
    } else {
      grfnparam->imp_type = IMPBACKWARDLATTICE;
    }
  }

  /* ------------------ STEP 6: COMPUTATION OF THE PRICE  ----------------------
   */

  /* Treats all the historical events first ( up to today ) and deletes them */
  err = srt_f_grfn_val_historical_events(&sptr, &gd, &und_info);
  if (err) {
    free_the_world(sptr, &gd, &und_info);
    return err;
  }

  /* Calculate the price ( + other info ) of the GRFN deal (if not
   * skip_evaluation) */
  err = srt_f_grfn_val_GrfnDeal(&gd, sptr, und, grfnparam, answer, &und_info);

  /* Copy the results of the last MonteCarlo path (from gd.cells) */
  if ((!err) && (cellcontents)) {
    /* The tableau could have been extended or shrunk (but only for rows) */
    memset(&(cellcontents[0][0]), 0, init_nrows * init_ncols * sizeof(double));

    for (i = 0; i < ((init_nrows > *nrows) ? *nrows : init_nrows); i++) {
      for (j = 0; j < *ncols; j++)
        cellcontents[i][j] = gd.cells[i][j];
    }

    for (i = ((init_nrows > *nrows) ? *nrows : init_nrows); i < init_nrows;
         i++) {
      for (j = 0; j < *ncols; j++)
        cellcontents[i][j] = 0.0;
    }

  } /* END copy of last path values */

  /* Sets the KNOWN Payments report (pay dates  , event dates and cahsflow) as
   * stored in the GrfnDeal */
  if ((!err) && (knownpayments)) {
    /* The tableau could have been extended or shrunk (but only for rows) */
    memset(&(knownpayments[0][0]), 0, init_nrows * 3 * sizeof(double));

    for (i = 0; i < ((init_nrows > gd.numpay) ? gd.numpay : init_nrows); i++) {
      knownpayments[0][i] = gd.pay_dates[i];
      knownpayments[1][i] = gd.event_pay_dates[i];
      knownpayments[2][i] = gd.pay_amounts[i];
    }

  } /* END copy of known payment dates and cash flows */

  /* ------------------------- STEP 7: FREE EVERYTHING  ----------------------
   */

  /* Detach GRFN events and free GRFN deal */

  free_the_world(sptr, &gd, &und_info);

  /* Return the error if any  , or a success message */
  return err;

} /* END Err srt_f_grfn(...) */

/* ------------------------------------------------------------------------------------
 */

/* ------------------------------------------------------------------------------------
   FUNCNAME        :srt_f_grfn_val_GrfnDeal

   DESCRIPTION     :value a grfn deal initialized by grfn_init_GrfnDeal.
                    This function decides whether to call
                                MCCore or TreeCore depending on the value of
   grfnparam->imp_type  , which was set in the function grfn_init_GrfnDeal.  One
   thing it does which ought to be done by MCCore is remove nodes which have no
   events in the case of Monte Carlo and LGM (in this case we can use the
   jumping numeraire.) This function also sets the value of grfnparam->step_num
   and grfnparam->node_dim for the tree functions.  The first of these is
                                probably redundant because the tree function
   themselves will call the SrtStp function create_index(); the second of these
   is a legitimate piece of information that GRFN must somehow pass to the
                            valuation function.

                                        The *answer field is passed to the
   appropriate valuation function will set it.

   -------------------------------------------------------------------------*/

Err srt_f_grfn_val_GrfnDeal(
    GrfnDeal *gd, SrtStpPtr sptr, SrtUndPtr und, SrtGrfnParam *grfnparam,
    void *answer, /* This is actually a pointer to a IOStruct== IOList */
    SrtUndInfo *und_info) {
  Err err = NULL;
  long last_index;
  double premium;
  SrtMdlType mdl_type;
  SrtMdlDim mdl_dim;
  SRT_Boolean status = SRT_NO;

  /*  Set an initial value of 0.0 for the Premium in the Input/Output list */
  err = srt_f_IOstructsetpremium((SrtIOStruct *)answer, status, 0.0, "");
  if (err)
    return err;

  /*  Set an initial value of 0.0 for the StDev in the Input/Output list */
  err = srt_f_IOstructsetstdev((SrtIOStruct *)answer, status, 0.0, "");
  if (err)
    return err;

  /* If skip evaluation: do not evaluate the deal */
  if (grfnparam->skip_evaluation == SRT_NO) {
    /* There are time steps left  */
    if (sptr) {
      /* Creates indexes for the time steps (0  , 1  , 2  ,...) */
      last_index = gd->sptr_last_index = create_index(sptr);

      err = get_underlying_mdltype(und, &mdl_type);
      if (err)
        return err;

      err = get_underlying_mdldim(und, &mdl_dim);
      if (err)
        return err;

      /* Method is Forward: Monte-Carlo  **/
      if (grfnparam->imp_type == IMPMONTECARLO) {
        /* Initialise the global histogram structures */
        err = srt_f_ininewhist(grfnparam->num_MCarlo_paths);
        if (err)
          return (err);

        if (mdl_type == BDT)
          return (serror("Cannot use BDT type models in Monte Carlo yet"));

        /* Set the number of columns in the Grfn tableau for the Node dimension
         */
        grfnparam->node_dim = gd->sswidth;

        /* This will be executed for ( mdl_type ==LGM || mdl_type
         * ==EQD_BLACK_SCHOLES) */
        if (und_info->jumping) {
          /* With Jumping Numeraire  , no need to save steps with no evemts */
          sptr = rem_eventless_nodes(sptr, gd->today);
          gd->sptr_last_index = create_index(sptr);
        }

        /* Branching to the core of the Monte Carlo discretisation */
        /*In the Case of a callable product  , we use the optimization
          method which is called from MCCORE_ExFr*/

        if (grfnparam->exfrontier != 100) {
          err =
              MCCoreExFrontier(grfnparam, sptr, gd,
                               (EvalEventFct)grfn_eval_event, answer, und_info);
          if (err)
            return (err);
        }

        else {
          err = MCCore(grfnparam, sptr, gd, (EvalEventFct)grfn_eval_event,
                       answer, und_info);
          if (err)
            return (err);
        }

      } /* END if grfnparam->imp_type == IMPMONTECARLO */

      else
      /* Method is Backward: Finite Difference Method (Tree  , PDE  , McTree...)
       */
      {
        /* Set the number of time steps */
        grfnparam->step_num = last_index;

        /* Set the number of columns in the Grfn tableau for the Node dimension
         */
        grfnparam->node_dim = gd->sswidth;

        /* Branching to the core of the Tree/Backward discretisation */
        err = TreeCore(und, grfnparam, sptr, gd, (EvalEventFct)grfn_eval_event,
                       answer, und_info);
        if (err)
          return (err);
      } /* END branching to TreeCore */

      /* Get the PV from the IO list of the discretisation method (for futures
       * cash flows) */
      err = srt_f_IOstructgetpremiumval(*(SrtIOStruct *)answer, &premium);
      if (err)
        return (err);

      /* Add the PV of the Past Events to the current PV */
      premium += gd->pv_of_past;

      /* Set the Premium in the IO list */
      err =
          srt_f_IOstructsetpremium((SrtIOStruct *)answer, status, premium, "");
      if (err)
        return (err);

    } /* END if stpr */

    else

    {
      /* Since there are no time steps in the future: the PV is the PV of past
       */
      premium = gd->pv_of_past;

      err =
          srt_f_IOstructsetpremium((SrtIOStruct *)answer, status, premium, "");
      if (err)
        return (err);
    }

  } /* END if( grfnparam->skip_evaluation == SRT_NO ) */

  /* Return a success message */
  return err;

} /* END Err srt_f_grfn_val_GrfnDeal(...) */

/* ------------------------------------------------------------------------------------------
 */

/* This function immediately evaluates historical events and discrads them from
 * the steps */

Err srt_f_grfn_val_historical_events(SrtStpPtr *step, GrfnDeal *gd,
                                     SrtUndInfo *und_info) {
  Err err = NULL;
  SrtStpPtr stp;
  double cash_flow;

  /* Go to the first time step in the list */
  if (!step)
    return NULL;
  stp = gototop((*step));

  /* Sets the "is history" flag to true for evaluation */
  gd->is_history = SRT_YES;

  /* Loop on all the steps until the step corresponding to today (+ eoday)  ,
   * and deletes them one by one (but today)*/
  while (stp->ddate < gd->today + gd->end_of_day_index) {
    /* Evaluate the event attached to the current time step (we do not care
     * about cashflow: only PAY has an impact ) */
    err = grfn_eval_event((GrfnEvent *)(stp->e), NULL, /* State variables */
                          gd, NULL, /* Discounted Vector in tree          */
                          (EvalEventDfsFct)srt_f_calc_grfn_event_dfs,
                          und_info,    /* Underlying Information             */
                          &cash_flow); /* CashFlow calculated via the comll  */
    if (err)
      return err;

    /* Frees the event attached to the step  */
    grfn_free_GrfnEvent((GrfnEvent *)stp->e);
    stp->e = (SrtEvent)NULL;

    /* Frees the step (if before today) and moves on to the next one ( it should
     * exist: today has been created ) */
    if (stp->ddate < gd->today) {
      stp = free_node(stp);
      stp->prev = NULL;
    } else
    /* Just moves on to the next step */
    {
      stp = stp->next;
    }

  } /* END of loop on all time steps */

  /*	Return to today */
  stp = gototop(stp);
  /*	while ( ((long) (stp->ddate + 1.0e-08)) < gd->today) stp = stp->next; */
  if (((long)(stp->ddate + 1.0e-08)) != gd->today) {
    err = "No today in stp list";
    return err;
  }

  /* Makes sure all needed past variable information (va  ,vb...) is attached to
   * today's step */
  err = grfn_attach_past_info_to_future(gd, (GrfnEvent **)(&(stp->e)));
  if (err)
    return err;

  /* Get the input step to point to the first step on or after today's (+ end of
   * day) step */
  *step = gototop(stp);

  /* Sets the "is history" flag to false for future evaluation */
  gd->is_history = SRT_NO;

  return NULL;

} /* END Err srt_f_grfn_val_historical_events(...) */

/* =================================================================================
 */
/* =================================================================================
 */
/* =================================================================================
 */

Err srt_f_grfn_ex_frontier(
    SrtUndPtr und, SrtGrfnParam *grfnparam, long numeventdates,
    Date **eventdates, long *nrows, long *ncols, GrfnCell ***sprdsht,
    long numgrng, GrfnRng *grng, long auxwidth, long *auxlen, double **aux,
    SrtIOStruct *iolist,    /* This is actually a pointer to the IoList */
    double **cellcontents,  /* The results of the last path when using MC */
    double **knownpayments, /* A report of known CF [2]  , paid on pay dates [0]
                               , according to event dates [1] ([0..*nrows-1])*/
    double *exfrontier) {
  Err err = NULL;
  int numex; /*number of exercize dates  , equal to the length of the
                appropriate aux*/
  int exfr;  /*the index of the auxiliary containing the exercise frontier*/
  int n;
  int j;
  double *storeexfr;
  long num_paths;
  double price;
  double pricesh;
  double gradient, newgradient;
  long s;
  double gamma;
  double GAMMACONST;
  double realshift;
  /*initialization*/
  GAMMACONST = grfnparam->gammaexfrontier;
  exfr = grfnparam->exfrontier;
  numex = auxlen[exfr];
  num_paths = grfnparam->num_MCarlo_paths;
  grfnparam->num_MCarlo_paths = 5000;

  /* Memory allocation*/
  storeexfr = dvector(0, numex - 1);

  /* storage and no previous exercise condition*/
  if (grfnparam->recpay == SRT_RECEIVER) {
    for (j = 0; j < numex - 1; j++) {
      storeexfr[j] = aux[exfr][j];
      aux[exfr][j] = 0;
    }
  }

  else {
    for (j = 0; j < numex - 1; j++) {
      storeexfr[j] = aux[exfr][j];
      aux[exfr][j] = 999;
    }
  }

  for (n = numex - 2; n >= 0; n--) {

    gradient = 0.0;
    newgradient = 100;
    if (n == numex - 2)
      aux[exfr][n] = storeexfr[n];
    else
      aux[exfr][n] = aux[exfr][n + 1];

    s = 1;
    j = 0;

    while (j < 25) {

      if ((j == 0) && (n == numex - 2))
        grfnparam->exfrontierfirst_time = SRT_YES;

      grfnparam->exfrontiercounter = j;

      err = srt_f_grfn(und, grfnparam, numeventdates, eventdates, nrows, ncols,
                       sprdsht, 0, NULL, auxwidth, auxlen, aux, iolist,
                       cellcontents, knownpayments);

      if (!err) {
        srt_f_IOstructgetpremiumval(*iolist, &price);
      }

      grfnparam->exfrontierfirst_time = SRT_NO;

      realshift = shift * aux[exfr][n];
      aux[exfr][n] += shift;

      err = srt_f_grfn(und, grfnparam, numeventdates, eventdates, nrows, ncols,
                       sprdsht, 0, NULL, auxwidth, auxlen, aux, iolist,
                       cellcontents, knownpayments);

      if (!err) {
        srt_f_IOstructgetpremiumval(*iolist, &pricesh);
      }

      newgradient = (pricesh - price) / shift;

      aux[exfr][n] -= shift;

      if (gradient * newgradient <= 0)
        s++;
      gamma = (GAMMACONST - (GAMMACONST - 0.1 * GAMMACONST) / (n + 1)) / s;
      aux[exfr][n] += gamma * newgradient;
      if (grfnparam->recpay == SRT_RECEIVER) {
        if (aux[exfr][n] > aux[exfr][n + 1])
          aux[exfr][n] = aux[exfr][n + 1];
      } else {
        if (aux[exfr][n] < aux[exfr][n + 1])
          aux[exfr][n] = aux[exfr][n + 1];
      }
      gradient = newgradient;
      j++;
    }
  }

  grfnparam->exfrontierlast_time = SRT_YES;

  for (n = 0; n < numex; n++) {
    exfrontier[n] = aux[exfr][n];
  }

  grfnparam->num_MCarlo_paths = num_paths;

  err = srt_f_grfn(und, grfnparam, numeventdates, eventdates, nrows, ncols,
                   sprdsht, 0, NULL, auxwidth, auxlen, aux, iolist,
                   cellcontents, knownpayments);

  return err;
}
