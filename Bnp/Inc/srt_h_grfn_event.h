/* ==========================================================

   FILENAME:     srt_h_grfn_event.h

   PURPOSE:      functions to work with Grfn Events
                 once a model has been discretised

   ========================================================== */
#ifndef SRT_H_GRFN_EVENT_H
#define SRT_H_GRFN_EVENT_H

#include "grf_h_pubtypes.h"

/* Creates and attaches GrfnEvents to SrtStpPtr */
Err srt_f_attach_grfn_events_to_steps (
			GrfnDeal    *gd,
			SrtStpPtr   step   );

/* The function that destroys all the events attched to the steps */
void srt_f_dettach_grfn_events_from_steps (SrtStpPtr sptr);


/* The function to pass to grfn_eval_event to calculate the DF at one simulation point */
void srt_f_calc_grfn_event_dfs(
		GrfnEvent   *event,
		SrtSample   *sample,
		SrtUndInfo  *und_info);

/* This function computes the df using the domestic underlying for payments triggered from the past */
void srt_f_calc_grfn_hist_dom_dfs(
		GrfnEvent  *event,
		SrtSample  *sample,
		SrtUndInfo *hist_und_info);
#endif
