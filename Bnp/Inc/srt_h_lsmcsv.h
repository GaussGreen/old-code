/* ------------------------------------------------------------------------------
   AUTHOR: E. FOURNIE & al 

   DATE : Septembre 99

   FILENAME:    srt_h_lsmcsv.h

   PURPOSE:     

   MODIFICATION:

   ------------------------------------------------------------------------------ */

#ifndef SRT_H_LSMCSV_H
#define	SRT_H_LSMCSV_H

Err srt_f_csv_lsm(
	 SrtGrfnParam   *grfnparams,			/*grfn parameters*/
	 SrtStpPtr       stp,                   /* step pointer*/
	 GrfnDeal       *gd,                 /* GRFN deal description*/
	 EvalEventFct    evalcf,             /* cash-flow evaluation function*/
	 SrtIOStruct    *iolist,               /* list of requests*/ 
	 SrtUndInfo     *und_info           /* underlying info*/
	 );

#endif
/*--------------------------------- End of File -------------------------------------*/ 