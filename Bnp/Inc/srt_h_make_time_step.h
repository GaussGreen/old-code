/* ========================================================================

  FILENAME:  srt_h_init_time_step.h

  PURPOSE:   When you want to make time steps for a discretisation using Grfn

  ========================================================================== */
#ifndef SRT_H_MAKE_TIME_STEP_H
#define SRT_H_MAKE_TIME_STEP_H

/*=========================================================*/
/*======discretize a model in time ========================*/
/*=========================================================*/

Err srt_f_make_time_step(long numeventdates, Date *eventdates, SrtUndPtr und,
                         SrtStpPtr *pointer_to_info, SrtGrfnParam *grfnparam,
                         SrtUndInfo *und_info);

#endif
