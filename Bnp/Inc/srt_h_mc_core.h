/* ------------------------------------------------------------------------------ 
   FILENAME:      srt_h_mc_core.h

   PURPOSE:       THE Core function for a Monte Carlo discretisation 
   ------------------------------------------------------------------------------ */
#ifndef SRT_H_MC_CORE_H
#define SRT_H_MC_CORE_H

Err free_MCrandSrc();

Err MCCore	( 	
		SrtGrfnParam   *grfnparam, 
		SrtStpPtr       info, 
		GrfnDeal       *gd,
		EvalEventFct    evalcf, 
		void           *iolist,
		SrtUndInfo     *und_info
		) ;

Err MCOptimizeExFrontier(
						 SrtGrfnParam  *grfnparam, 
		SrtStpPtr      step, 
		GrfnDeal      *gd,
		EvalEventFct   evalcf, 
		void 		  *iolist,
		SrtUndInfo	  *und_info
		);

Err MCCoreExFrontier(
		SrtGrfnParam  *grfnparam, 
		SrtStpPtr      step, 
		GrfnDeal      *gd,
		EvalEventFct   evalcf, 
		void 		  *iolist,
		SrtUndInfo	  *und_info
		);


#endif