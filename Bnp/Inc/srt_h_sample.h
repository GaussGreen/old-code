/* =========================================================================
   
   FILENAME:    srt_h_sample.h 

   PURPOSE:     Functions to deal with several discretised unerlyings
                                 SrtSample

   ========================================================================= */

#ifndef SRT_H_SAMPLE_H
#define SRT_H_SAMPLE_H

/* -------- Maximum Number of Underlyings Authorised in a Grfn Tableau -------- */

#define   MAXUNDERLYING       100

/* ----------- Maximum Number of State Varaibles per Underlying --------------- */

#define   MAXSTATEVAR         10

/* --------------------------- State Variables Code ----------------------- */
#define   SHORT_RATE          0
#define   F_0_t		          SHORT_RATE         /* 0 */
#define   PHI                 1
#define   SPOT		          0
#define   STATEVAR	          MAXSTATEVAR - 1	 /* 9 */ 	
/* Do not try to define X: it will not compile */
#define   X1                  STATEVAR		     /* 9 */      
#define   X2                  MAXSTATEVAR - 2	 /* 8 */ 	
#define   PHI1		          PHI                /* 1 */ 
#define   PHI2		          PHI + 1            /* 2 */ 
#define   CROSSPHI	          PHI + 2            /* 3 */ 

#define   SIGMA               5     /* For stoch vol */
#define   BT                  6     /* money market account */
/* ------------------------- SrtSample ------------------------------------- */


/* ---------------- One Single Underlying: all the state variables --------- */
typedef struct {
  double		sv[MAXSTATEVAR]	;
} SrtUnderPt ;

/* ----------------- The SAMPLE: All the Underlyings ------------------------ */
typedef struct {
  
  SrtUnderPt     und[MAXUNDERLYING] ;

  long			 pathnum;

  double         numeraire;
  int			 numeraire_index;	/* discounting underlying index */
									/* it could be different from 0 */
									/* eg FX_STOCH_RATES */
} SrtSample ;


/* -------------------------------------------------------------------------- */

SrtErr srt_f_undptset(SrtUnderPt *undpt, int svarindex, double val) ;
SrtErr srt_f_undptget(SrtUnderPt *undpt, int svarindex, double *val);

SrtErr srt_f_getundptfromsam(	SrtSample *sam, int undindex, 
				SrtUnderPt *undpt ) ; 
SrtErr srt_f_samsetpath(SrtSample *sam, long pathnum)  ;
SrtErr srt_f_samgetpath(SrtSample *sam, long *pathnum) ;


#define sam_get(sam,index,field) (sam.und[index].sv[field])

#define samptr_get(samptr,index,field) ((samptr)->und[index].sv[field])


#endif