#ifndef SRT_H_TREBDT_H
#define SRT_H_TREBDT_H

Err srt_f_trebdt1f
	(
	SrtUndPtr und, /* underlying parameters */
	SrtGrfnParam   *grfnparam, /* model parameters */
	SrtStpPtr stp, /* discretization of deal in time, wi/ events attached*/
	GrfnDeal	   *gd,/* deal descriptor structure */
	EvalEventFct evalcf, /* cashflow evaluator */
	SrtIOStruct *iolist,
	SrtUndInfo *und_info
	) ;

Err srt_f_trebdt2f
	(
	SrtUndPtr und, /* underlying parameters */
	SrtGrfnParam   *grfnparam, /* model parameters */
	SrtStpPtr stp, /* discretization of deal in time, wi/ events attached*/
	GrfnDeal	   *gd,/* deal descriptor structure */
	EvalEventFct evalcf, /* cashflow evaluator */
	SrtIOStruct *iolist,
	SrtUndInfo *und_info
	) ;

#endif

