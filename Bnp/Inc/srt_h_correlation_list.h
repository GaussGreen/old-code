/* -------------------------------------------------------------------------

 	AUTHOR		: O Van Eyseren
	DATE		: Feb 23 1995
	FILE NAME	: srt_h_correlation_list.h
	PURPOSE		: use a SrtList to store a TermStructure of correlation matrixes
	              gives a range of functions that allow to use this list
				  properly
   ------------------------------------------------------------------------- */
#ifndef CORREL_LIST_H
#define CORREL_LIST_H

/* ----------------------------------------------------------------------- 
   The global SrtList used to store the TermStructure of correlation 
   ----------------------------------------------------------------------- */
typedef SrtListHdr SrtCorrLst, *SrtCorrLstPtr ;


/* ----------------------------------------------------------------------- 
   The pval of an element of the SrtCorrLst, used to store at each date
   the correlation matrix 
   ----------------------------------------------------------------------- */
typedef struct srtcorrlstval
{
	Ddate	date;
	double 	time;
	String	*und_names;   /* Sorted by alphabetical order */
	double 	**correl;     /* with und[1..n] vs und[1..n] */
	double 	**coeff;	  /* USed when correlating the underlyings*/
	int 	ncorrel;      /* = n_und*(n_und-1)/2: all but diagonal */
	int 	nund;
} SrtCorrLstVal;



/* ----------------------------------------------------------------------- 
   Returns the adress of the static used to store the correlation list.
   ----------------------------------------------------------------------- */
SrtCorrLstPtr srt_f_GetTheCorrelationList (void);

/* ----------------------------------------------------------------------- 
   Allocate space for the STATIC SrtCorrLstPtr corr_list, and gives it a name...
   ----------------------------------------------------------------------- */
Err create_correlation_list (String corr_list_name);

/* ----------------------------------------------------------------------- 
   Destroys the STATIC SrtCorrLstPtr corr_list
   ----------------------------------------------------------------------- */
Err destroy_correlation_list ();


/* ----------------------------------------------------------------------- */
/* Corresponding Functions related to the SrtCorrLst                       */


SrtErr srt_f_corrlstdelete(SrtCorrLstPtr *cls);


Err srt_f_corrlstcreate(SrtCorrLstPtr *cls, String list_name);


/* ----------------------------------------------------------------------- 
   Frees a SrtCorrLstVal (the pval of an element of the SrtCorrLst list
   ----------------------------------------------------------------------- */
Err srt_f_corrvalfree(void * corvalptr);



/* ----------------------------------------------------------------------- 
   Initialises (fills in) the STATIC SrtCorrLstPtr corr_list according to 
   data taken from spreadsheet...
   ----------------------------------------------------------------------- */
Err srt_f_init_Corr_TermStruct( 
					int       ndates, 
					int       ncorr, 
					double  **correl,
					double   *dates, 
					String  **und_names,
					SrtCorrLstPtr *the_corr_list);



/* ----------------------------------------------------------------------- 
   Initialises and fills in a local SrtCorrLst, extracting information
   from the STATIC SrtCorrLstPtr __corr_list according to 
   the underlyings used in the deal
   		und_names is a vector with the underlying names: [0..n_und]
   ----------------------------------------------------------------------- */
Err srt_f_make_deal_corrlist(
					String *und_names, 
					int n_und, 
					String list_name,
					SrtCorrLstPtr big_cls,
					SrtCorrLstPtr *deal_corr_list);


/* ----------------------------------------------------------------------- 
   Gets the correlation for two underlyings at a certain time (in yrs frm tdy), 
   according to the SrtCorrLstPtr cls that is passed as an argument. 
   There is no interpolation:
		the correlation is supposed to be constant between two dates
   ----------------------------------------------------------------------- */
Err srt_f_get_corr_from_CorrList(SrtCorrLstPtr cls, 	
	String	und_name1, String und_name2, double time, double *value);



#endif


