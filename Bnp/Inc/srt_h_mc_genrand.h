/* =====================================================================================
   FILENAME:   srt_m_mc_genrand.h
	
   PURPOSE:    random number generation and handling
   ===================================================================================== */

#ifndef SRT_H_MCGENRAND_H
#define SRT_H_MCGENRAND_H

/* The main structure to store random numbers */

typedef struct{
    double           ***r;        /* r[path][brownian][step]  */
    SrtMCSamType        sample_type;
    long                pl,       /* Path Lower index */
		                ph,       /* Path higher index */
						sl,       /* Step lower index */
						sh;       /* Step higher index */
	
	long                cur_path;   /* Index of the current path */
	
	long                lastindex;  /* Index of the last generated path */
	long 			    seed;    /* Initial seed for random number generation */
	
	long                nbr;       /* Number of Brownian motions used*/
    
	double            **correl;        /* Correlation matrix [0..nbr-1][0..nbr-1] */
    double            **coeff;        /* Linear coefficients derived from correl */
	SRT_Boolean             need_to_correl;	/* TO prevent the use of correl if not needed*/
    
}SrtMCRanSrc;

/* ----------------------------------------------------------------------- */



/* ----------------------------------------------------------------------- */

/*New function that does the same as the previous one but which does not take a SrtStpPtr as an input but only
the vector of times_at_steps  which is the only thing really needed*/

Err srt_f_New_gen_random(
				double     ***rand, 
				int           first_path, 
				int           last_path, 
				int           first_step, 
				int           last_step, 
				int           nbr, 
				SrtMCSamType  gen_method, 
				long         *seed, 
				double       *times_at_steps);

/* perform some sort of setup, allocation, whatever */
Err SrtMCRanSrc_start(
		SrtMCRanSrc  *mcrs, 
		long          pl, long ph,
		long          sl, long sh,
		long          nbr, 
		SrtMCSamType  gen_method, 
		long          seed,
		SrtStpPtr     top);

/* return a matrix of random gaussian samples, the dimension
	and computation method of which was determined by a 
	previous call to MCRSrc_start(), this function do once what SrtMCRanSrc_next
    does for each path */

/*New Function which is the same as the previous one only taking times_at_steps instead as a SrtStpPtr
as an input*/

Err SrtNewMCRanSrc_start(
		SrtMCRanSrc *mcrs, 
		long         pl,
		long         ph,
		long         sl,
		long         sh,
		long         nbr, 
		SrtMCSamType gen_method, 
		long         seed,
		double		*times_at_steps);


Err SrtMCRanSrc_init(SrtMCRanSrc *mcrs);

/* return a matrix of random gaussian samples, the dimension
	and computation method of which was determined by a 
	previous call to MCRSrc_start() and store the result in mrcs->cur_path*/

Err SrtMCRanSrc_next(SrtMCRanSrc *mcrs);

/* Correlates the numbers depending on SrtStpPtr*/
Err SrtMCRanSrc_correl(
		SrtMCRanSrc *mcrs, 
	SrtStpPtr	stp);

/* Undo whatever MCRSrc_init() did. */
Err SrtMCRanSrc_end(SrtMCRanSrc *mcrs); 

#endif
