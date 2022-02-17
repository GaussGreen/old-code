/* -------------------------------------------------------------------------

   FILE NAME	: srt_f_grfn_event.cxx

   PURPOSE		: Several useful functions to instantiate a call to
                  the Grfn langage (COMLL...)

   ------------------------------------------------------------------------- */
#include "grf_h_all.h"
#include "math.h"
#include "srt_h_all.h"

/* -------------------------------------------------------------------------

    FUNCTION        :   srt_f_attach_grfn_events_to_steps

    DESCRIPTION     :   Attach grfn events to steps in SrtStpLst inside gd.
                                  attach and delete Grfn events to a SrtStpPtr
                  (this used to be in grf_f_get_event.cxx)

                        This function parses events in the past.
                        If there are any events in the past        , it will
   augment the event at today with a COMLL command to restore the environment at
   the end history.

                        It then procedes by looping forward in time over the
                        nodes in step.

                        For each node of step:

                            o   Set values of all variables in the gd
                                (Grfn eval environment)

                            o   If a financial event is specified at that date ,
                                create event.


    IMPORTANT       :   This function assumes the following regarding step

                            o   First node corresponds to today
                            o   Node date field > 0 iff date is an event date
                            o   ddate is set for all dates
                            o   Nodes are in increasing order.

                            If step = NULL we assume all events are historical.


    Date            :   11 Jan 1995
    Description     :   Bug Fix for the following tableau:

                        Event Date  c[0]
                        18/01/95    1234

                        When end_of_day_flg is set to 0 (i.e. Not end of day)
                        correct value of 1234 is returned;
                        When end_of_day_flg is set to 1 (i.e. At end of day)
                        *** Crash *** Returned value should be ZERO.

                        Fix:    Check to see if step is NULL  (as returned by
                                srt_f_make_time_step(...) - if NULL we have
                                no event dates in the future!

   ------------------------------------------------------------------------- */

Err srt_f_attach_grfn_events_to_steps(GrfnDeal *gd, SrtStpPtr step) {
  /*..Local variables..*/
  SrtStpPtr top = NULL, bot = NULL, firstunknownstp = NULL, stp = NULL;
  Err error = NULL;

  /* Check that there are time steps: if not        , return */
  if (step != NULL) {
    top = gototop(step);
    bot = gotobot(step);
    stp = top;
  } else {
    return (NULL);
  }

  /* Sets today's and previous dates in the GrfnDeal */
  if (gd->first_unkn_index == 0) {
    /* The first unknown event date is in the future: set today */
    gd->nowdt = gd->today;
    gd->prvdt = gd->nowdt;
  } else {
    /* There were dates in the past */
    gd->nowdt = gd->event_dates[0];
    gd->prvdt = gd->nowdt;
  }

  /* Start at the top left hand corner of tableau (minus one row for potential
   * AM event )*/
  gd->I = -1;
  gd->J = 0;

  /* --------- PARSE HISTORICAL EVENTS FIRST --------- */

  gd->is_history = SRT_YES;

  /* While the next event date is before the first future (unknown) date */
  while (stp->ddate < gd->today + gd->end_of_day_index) {

    /* Creates the historical event        , taking fixings        , cvg into
     * account and evaluating what can be done */
    error = grfn_create_historical_event(
        gd, (GrfnEvent **)(&(stp->e)), stp->ddate,
        (stp->next == NULL) ? 0.0 : stp->next->ddate);
    if (error)
      return error;

    /* Move on to the next step ( where the next event will be attached ) */
    stp = stp->next;

    /* Checks the next step exists ; if not: return */
    if (!stp) {
      stp = top;
      return NULL;
    }
  }

  /* -------- CHECK THAT THE NEXT ROW INDEX IS THE FIRST UNKNOWN ONE
   * ----------- */

  if (gd->I + 1 != gd->first_unkn_index)
    return serror("Step dates mismatch in attach events to steps");
  else
    firstunknownstp = stp;

  /* --------- PARSE NON HISTORICAL EVENTS AFTER --------- */

  gd->is_history = SRT_NO;

  /* Loop through all the time steps        , one by one        , starting from
   * today  */
  while (stp) {
    /* Create and attach to the current time step the corresponding Grfn event
     */
    error =
        grfn_create_future_event(gd, (GrfnEvent **)(&(stp->e)), stp->ddate,
                                 (stp->next == NULL) ? 0.0 : stp->next->ddate);
    if (error)
      return error;

    /* Move on to the next discretisation date */
    stp = stp->next;

  } /* END while */

  /* Cut useless events off the end by finding last stp with real event */
  stp = bot;

  while (stp && stp->prev && stp->e == NULL) {
    stp = stp->prev;
  }

  if (stp && stp->next) {
    stp->next->prev = NULL;
    free_list(stp->next);
    stp->next = NULL;
  }

  return NULL;

} /* END Err srt_f_attach_grfn_events_to_steps(...)  */

/* ------------------------------------------------------------------------- */

/* The function that destroys all the events attched to the steps */
void srt_f_dettach_grfn_events_from_steps(SrtStpPtr sptr) {
  /* If we have a list(!) ... */

  if (sptr) {
    /* ... start at the top */
    sptr = gototop(sptr);

    /* While not at the end of the list ... */
    while (sptr) {
      /* ... if a financial event "e" is associated with the node ... */
      if (sptr->e) {
        /* ... free "e" */
        grfn_free_GrfnEvent((GrfnEvent *)sptr->e);
        sptr->e = (SrtEvent)NULL;
      } /* end if */

      /* Move to the next node */
      sptr = sptr->next;

    } /* end while */

  } /* end if */

} /* END void srt_f_dettach_grfn_events_from_steps (...) */

/* ------------------------------------------------------------------------- */

/* The function to pass to grfn_eval_event to calculate the DF at one simulation
 * point */
void srt_f_calc_grfn_event_dfs(GrfnEvent *event, SrtSample *sample,
                               SrtUndInfo *und_info) {
  int i, j;
  SrtUndPtr und;
  SrtMdlType mdl_type;
  SrtMdlDim mdl_dim;

  /* Loop through the underlyings in UndInfo one by one */
  for (i = 0; i < und_info->no_of_underlyings; i++) {
    /* Get the Underlying and its model */
    und = lookup_und(und_info->und_data[i].und_name);
    get_underlying_mdltype(und, &mdl_type);
    get_underlying_mdldim(und, &mdl_dim);

    if (mdl_type == VASICEK) {
      Vasicek_Y_T_at_t_compute(event->dflen[i], sample, event->yp[i],
                               event->df[i], i, mdl_dim, mdl_type);

    } else {
      /* Compute the log of ZeroCoupon prices using the Sample */
      Y_T_at_t_compute(event->dflen[i], sample, event->yp[i], event->df[i], i,
                       mdl_dim, mdl_type);
    }

    /* Now converts zero rates to discounting factors: */
    for (j = 0; j < event->dflen[i]; j++)
      event->df[i][j] = exp(-event->df[i][j]);
  }

} /* END void srt_f_calc_grfn_event_dfs(...) */

/* ------------------------------------------------------------------------- */
