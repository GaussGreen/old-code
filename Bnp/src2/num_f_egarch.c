/********************************************************************/
/* FILENAME :      num_f_egarch.c
/*																	*/
/* The following functions compute coefficients for a GARCH			*/
/* algorithm from a bunch of input data								*/
/*																	*/
/********************************************************************/
 

#include "num_h_allhdr.h"
#include "num_h_egarch.h"
#include "math.h"

#define NULL_S1		1.0e-4
#define FTOL		1.0e-8
#define MAX_TEMP	1000
#define MIN_TEMP	1.0e-8
#define DEC_FAC		10.00
#define P_INC		0.0001
#define NOT_FIXED	0


/* global variables */
static GarchData *	garcharray = NULL;
static double *		global_s2_array = NULL;
static double *		global_zeta_array = NULL;
static long			global_count;
static double		global_s2_init;
static double		global_beta_0;
static double		global_beta_1;
static double		global_beta_2;
static long *		global_fixed;

/* Sort GarchArray */
static int compare_garchdata (const void * elem1, const void * elem2)
{
GarchData * element1 = (GarchData *) elem1;
GarchData * element2 = (GarchData *) elem2;

	if ((*element1).date > (*element2).date) return  1;
	if ((*element1).date < (*element2).date) return -1;

	return 0;
}


/* To compute LL */
static double LL_computation (double * x) 
{
int i;
double LL0, LL1, LL2, LL;

	/* Initialisation of variables */
	i = 1;
	if (global_fixed[1] == NOT_FIXED)
	{
		global_s2_array[1] = x[i];
		i++;
	}
	else
		global_s2_array[1] = global_s2_init;

	if (global_fixed[2] == NOT_FIXED)
	{
		global_beta_0 = x[i];
		i++;
	}

	if (global_fixed[3] == NOT_FIXED)
	{
		global_beta_0 = x[i];
		i++;
	}

	if (global_fixed[4] == NOT_FIXED)
	{
		global_beta_0 = x[i];
		i++;
	}
	
	/* Fill the different array */
	if (garcharray[1].s1 == 0) garcharray[1].s1 = NULL_S1;
	if (global_s2_array[1] == 0) global_s2_array[1] = NULL_S1;
	global_zeta_array[1] = log(garcharray[1].rate/garcharray[0].rate)
						   /garcharray[1].s1/global_s2_array[1];	
	for (i = 2 ; i <= global_count ; i++)
	{
		global_s2_array[i] = global_s2_array[i-1]*exp(global_beta_0+global_beta_1*global_zeta_array[i-1]+global_beta_2*(fabs(global_zeta_array[i-1]) - sqrt(2.0)/SRT_PI));
		if (garcharray[i].s1 == 0) garcharray[i].s1 = NULL_S1;
		if (global_s2_array[i] == 0) global_s2_array[i] = NULL_S1;
		global_zeta_array[i] = log(garcharray[i].rate/garcharray[i-1].rate)
								/garcharray[i].s1/global_s2_array[i];	
	}
		

	/* MAXIMISATION OF F = MINIMISATION OF -F */
	/* we want a maximisation and simulated annealing do a minimisation so we consider -LL */

	/* computation of LL0 */
	LL0 = 0.5 * global_count * log ( 2 * SRT_PI);
	
	/* computation of LL1 */
	LL1 = 0;
	for (i = 1 ; i <= global_count ; i++)
	{
		LL1 += 2 * (log(fabs(global_s2_array[i]) + log(fabs(garcharray[i].s1))));
	}
	LL1 *= 0.5; 

	/* computation of LL2 */
	LL2 = 0;
	for (i = 1 ; i <= global_count ; i++)
	{
		LL2 += (garcharray[i].s1 != NULL_S1 ? pow(global_zeta_array[i],2) : 0);
	}
	LL2 *= 0.5; 

	LL = LL0 + LL1 + LL2;

	return (LL);
}

Err srt_f_egarchmethod(Ddate	date_chosen,	/* date chosen to compute the GARCH */
					   Ddate *	dates,			/* input dates */
					   double *	rates,			/* input rates */
					   double *	s1,				/* input s1 */
					   long		num_data,		/* number of input data */
					   long *	fixed,			/* array describing if param are fixed or not */
					   double * s2_init,		/* output & input Initialisation of S2 */
					   double * beta_0,			/* output & input beta_0 in the Garch algorithm */
					   double * beta_1,			/* output & input beta_0 in the Garch algorithm */
					   double * beta_2,			/* output & input beta_0 in the Garch algorithm */
					   double * s2_last,		/* output last value of S2 */
					   double * zeta_last,		/* output last value of Zeta */
					   long		num_iter		/* number of iterations for the algorithm */
					  )
{
Err err = NULL;
int i, temp;
int num_params = 0; /* number of parameters not fixed s2_init, beta_0, beta_1, beta_2 */
double ** p = NULL;
double * y = NULL;


	/* Sort the input data */
	garcharray = (GarchData *) calloc(num_data, sizeof(GarchData));
	if (!garcharray)
		return serror("Memory Allocation Error");

	for (i = 0 ; i < num_data ; i++)
	{
		garcharray[i].date = dates[i];
		garcharray[i].rate = rates[i];
		garcharray[i].s1   = s1[i];
	}

	qsort((void *)garcharray, (size_t) num_data, sizeof(GarchData), compare_garchdata);

	/* Get the number of dates before the chosen one */
	for (i = 0 ; i < num_data ; i++)
	{
		if (garcharray[i].date > date_chosen)
		{	
			global_count = --i;
			break;
		}
	}

	if (i >= num_data)
	{
		srt_free(garcharray);
		return serror("Dates input are not sufficient");
	}

	/* Allocation of memory for different array and initialisation */
	global_s2_array		= (double *) calloc (global_count+1, sizeof(double));
	global_zeta_array	= (double *) calloc (global_count+1, sizeof(double));
	if ((!global_s2_array) || (!global_zeta_array))
	{
		srt_free(global_s2_array);
		srt_free(garcharray);
		return serror("Memory Allocation Error");
	}

	global_fixed = fixed;
	global_s2_init = *s2_init;
	global_beta_0 = *beta_0;
	global_beta_1 = *beta_1;
	global_beta_2 = *beta_2;

	/*********************************/
	/* Simulated annealing algorithm */
	/*********************************/
	/* get the number of non fixed params */
	for (i= 1; i <= 4 ; i++)
		if (fixed[i] == NOT_FIXED) num_params++;
	if (num_params == 0)
		return serror ("Need at least an optimisation");

	p = dmatrix(1,num_params+1,1,num_params);
	y = dvector(1,num_params+1);
	if ((!p) || (!y))
	{
		if (p) free_dmatrix(p,1,num_params+1,1,num_params);
		srt_free(global_zeta_array);
		srt_free(global_s2_array);
		srt_free(garcharray);
		return serror("Memory Allocation Error");
	}

	/* fill p and y */
	for (i = 1 ; i <= num_params+1 ; i++)
	{
		temp = 1;
		if (fixed[1] == NOT_FIXED) { p[i][temp] = *s2_init; temp++; }
		if (fixed[2] == NOT_FIXED) { p[i][temp] = *beta_0; temp++; }
		if (fixed[3] == NOT_FIXED) { p[i][temp] = *beta_1; temp++; }
		if (fixed[4] == NOT_FIXED) { p[i][temp] = *beta_2; }
	}

	for (i = 2 ; i <= num_params+1 ; i++)
	{
		p[i][i-1] += P_INC;
	}

	for (i = 1 ; i <= num_params+1 ; i++)
	{
		y[i] = LL_computation(p[i]);
	}
	
	/* Algorithme */
	err = simulated_annealing(
		p,        
		y,
		num_params,
		FTOL,
		num_iter,
		LL_computation,
		MAX_TEMP,
		MIN_TEMP,
		DEC_FAC);

	if (err)
	{
		if (y) free_dvector(y,1,num_params+1);
		if (p) free_dmatrix(p,1,num_params+1,1,num_params);
		srt_free(global_zeta_array);
		srt_free(global_s2_array);
		srt_free(garcharray);
		return (err);
	}

	temp = 1;
	if (fixed[1] == NOT_FIXED) { *s2_init = p[1][temp]; temp++; }
	if (fixed[2] == NOT_FIXED) { *beta_0 = p[1][temp]; temp++; }
	if (fixed[3] == NOT_FIXED) { *beta_1 = p[1][temp]; temp++; }
	if (fixed[4] == NOT_FIXED) { *beta_2 = p[1][temp]; temp++; }

	/* Fill the different array */
	LL_computation(p[1]);
	
	*s2_last = global_s2_array[global_count]; 
	*zeta_last = global_zeta_array[global_count]; 
	
	if (y) free_dvector(y,1,num_params+1);
	if (p) free_dmatrix(p,1,num_params+1,1,num_params);
	srt_free(global_zeta_array);
	srt_free(global_s2_array);
	srt_free(garcharray);

	return NULL;
}

#undef NULL_S1
#undef FTOL
#undef MAX_TEMP
#undef MIN_TEMP
#undef DEC_FAC
#undef P_INC