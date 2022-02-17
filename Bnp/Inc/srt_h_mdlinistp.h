/* ===========================================================================
   FILENAME:   srt_h_mdlinistp.h

   PURPOSE:    All the functions needed to attach at a time step all the
               informaton required for a discretisation that is not SIMULATION
			   dependent (sigma, taus, rho, ...) 
   
   FUNCTIONS:  - srt_f_irministp:      for any interest rate underlying
               - srt_f_lgministp:      specifically for LGM 
               - srt_f_betaetainistp:  for the EtaBeta model
               - srt_f_loginistp:      for BlackScholes type underlyings
               - srt_f_basicinistp:    for deterministic underlyings
		
   DESCRIPTION:  all the functions allocate space for a Srt...TmInf and
                 attach one (populated) to each time step
   =========================================================================== */

#ifndef SRT_H_MDLINISTP_H
#define SRT_H_MDLINISTP_H

/* ========================================================================

	IRM function: used for CHE and LGM 

   ======================================================================== */

Err srt_f_irministp(
	SrtStpPtr   stp,
	SrtUndPtr   und,
	int         und_index,
	SrtUndPtr   numeraire_und,
	SrtUndInfo   *und_info

)
;


/* ======================================================================

	LGM specific functions 

   ====================================================================== */


Err srt_f_lgministp(
	SrtStpPtr stp,
	TermStruct *ts,
	int index,
	SrtMdlDim mdl_dim,
	SrtUndPtr    numeraire_und,  /* THe Numeraire Used for Discounting (quanto's...) */
	SrtUndInfo   *und_info
		
)
;

/* ======================================================================

	ETABETA specific functions 

   ====================================================================== */


Err srt_f_betaetainistp(
	SrtStpPtr    stp,
	TermStruct   *ts,
	int          index,
	SrtMdlDim    mdl_dim);

/* ======================================================================

	Black Scholes specific functions 

   ====================================================================== */

Err srt_f_loginistp(
	SrtStpPtr   stp,
	SrtUndPtr   und,
	int         und_index,
	SrtUndPtr numeraire_und,
	SrtUndInfo   *und_info
)
;


/* ======================================================================

	Deterministic Model  specific functions 

   ====================================================================== */

Err srt_f_basicinistp(
	SrtStpPtr   stp,
	SrtUndPtr   und,
	int         und_index
)
;

/* -------------------------------------------------------------------------  
   When Dealing with Multi-Currencies, an extra initialisation is required:
   the Quanto adjustment in the drift
   ------------------------------------------------------------------------- */
  

Err srt_f_steps_set_quanto_adjustment(	
			SrtStpPtr     stp,
			SrtUndInfo   *und_info,
			int           und_index,
			SrtUndPtr     und,
			SrtUndPtr     numeraire_und);

/* Jumping Numeraire functions */

Err srt_f_fx_initstp(
					SrtStpPtr	stp,
					SrtUndInfo	  *und_info,
					SrtUndPtr	und,
					SrtUndPtr  dom_und,
					SrtUndPtr	for_und,					 
					int		und_index,
					SrtUndPtr numeraire_und);



#endif
