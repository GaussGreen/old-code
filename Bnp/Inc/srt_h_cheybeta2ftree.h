/* ------------------------------------------------------------------------------
  AUTHOR: E. FOURNIE & al 

  DATE : MARS 98

  FILENAME: srt_h_cheybeta2dtree.h
 ----------------------------------------------------------------------------- */

#ifndef SRT_H_CHEYBETA2DTREE_H
#define	SRT_H_CHEYBETA2DTREE_H



Err srt_f_cheybeta2dtree (
	 SrtUndPtr       und,        /* underlying  pointer*/
	 SrtGrfnParam   *grfnparam,   /*grfn parameters*/
	 SrtStpPtr       stp,        /*step pointer*/
	 GrfnDeal	   *gd,         /*GRFN deal description*/
	 EvalEventFct    evalcf,     /*cash-flow evaluation function*/
	 SrtIOStruct    *iolist,    /*list of requests*/ 
	 SrtUndInfo     *und_info   /*underlying info*/
	 );


#endif
