/*******************************************************************************
**                      Include Files                               
*******************************************************************************/
#include		<math.h"
#include		<SPFNCTNS.H>
#include        "utallhdr.h"

/******************************************************************************/

#define SWITCH(a,b) {double *temp; temp=a; a=b; b=temp;} 

/*******************************************************************************
*                                                                           
* FUNCTION     	: srt_f_splocspr()                                          
*                                                                           
* PURPOSE      	: Calculates the local spread       
*                                                                           
* DESCRIPTION  	:                            	    
*                                                                           
* PARAMETERS                                   			         
*		INPUT	: premium		- premium on European option         	    
*              	: fwd_price    	- forward price of underlying price 	    
*              	: vol           - annual volatility  	      	  	    
*              	: mat           - maturity, in years      	      	    
*              	: disc          - discount factor      	      		    
*              	: call_put      - type of option: 0 call, 1 put       	    
*                                                                           
* RETURNS      	: strike_guess  - guess at implied strike                    
*                                                                           
*******************************************************************************/
static double reg_call (double spot, double strike)
{
	return (spot > strike ? spot - strike: 0.0);
}

static double reg_put (double spot, double strike)
{
	return (spot < strike ? - spot + strike: 0.0);
}

static double dig_call (double spot, double strike)
{
	return (spot > strike ? 1.0: 0.0);
}

static double dig_put (double spot, double strike)
{
	return (spot < strike ? 1.0: 0.0);
}

/************ Principal Function ******************/

double srt_f_splocspr  (double	fwd_cms_quanto_spread,	/* forward */
						double	maturity,				/* maturity */	
						double	strike,					/* strike */
						double	disc,					/* discount factor */
						int		num_spots,				/* number of spots */
						double	*spots,					/* array of spots */
						int		num_dates,				/* number of 'maturities' */
						double	*dates,					/* array of 'maturities' */				   
						double	**vols,					/* array of vols */					   					   
						SrtCallPutType  call_put,			/* Call or Put */
						int		num_steps,				/* number of steps */		
						int		type,					/* regular = 0, digital = 1 */					   
						int		option,					/* euro = 0, amer = 1 */   
						SrtGreekType greek)				/* Greek */
{
int i, j;
long step, node;
double *node1, *node2, date, spot, local_vol, maxvol, max_stdev, maxvar;
double proba, dt, sdt, disc_dt, intval;
double answer;

double (*PAYOFF)(double, double);

	if ((call_put == SRT_CALL) && (type == 0)) /* regular */
		PAYOFF = reg_call;
	else if ((call_put == SRT_PUT) && (type == 0)) /* regular */
		PAYOFF = reg_put;
	else if ((call_put == SRT_CALL) && (type == 1)) /* digital */
		PAYOFF = dig_call;
	else if ((call_put == SRT_PUT) && (type == 1)) /* digital */
		PAYOFF = dig_put;

	if (!(node1 = dvector (-num_steps, num_steps)))
		return MEMORY_ERR;
	
	if (!(node2 = dvector (-num_steps, num_steps)))
	{
		free_dvector(node1, -num_steps, num_steps);
		return MEMORY_ERR;
	}

	maxvol = 0.0;
	for (i=0; i<num_spots; i++) 
	{
		for (j=0; j<num_dates; j++)
		{
			if (vols[i][j] > maxvol) maxvol = vols[i][j];
		}
	}

	dt = maturity / num_steps;
	sdt = sqrt (dt);
	max_stdev = maxvol * sdt;
	maxvar = maxvol * maxvol;
	disc_dt = pow (disc, 1.0 / num_steps);
	
	spot = fwd_cms_quanto_spread + num_steps * max_stdev;

	for (node = num_steps; node >= -num_steps; node--)
	{
		node1[node] = PAYOFF (spot, strike);
		spot -= max_stdev;
	}

	date = maturity;
	for (step = num_steps-1; step >= 0; step --)
	{
		date -= dt;

		spot = fwd_cms_quanto_spread + step * max_stdev;


		for (node = step ; node >= -step ; node--)
		{
			local_vol = lin_interp_2d (spot, date, spots, dates, vols, num_spots, num_dates);
			proba = 0.50 * local_vol * local_vol / maxvar;
			node2[node] = disc_dt * 
						  ( proba * (node1[node + 1] + node1[node - 1])
							+ (1.0 - 2 * proba) * node1[node] );

			if (option == 1) /* American */
			{
				intval = PAYOFF (spot, strike);
				if (node2[node] <= intval) 
					node2[node] = intval;
			}

			spot -= max_stdev;
		}

		SWITCH (node1, node2);
	}

	if (greek == PREMIUM)
	{
		answer = node1[0];
	}
	else if ((greek == DELTA) || (greek == DELTA_FWD))
	{
		answer = disc_dt * 0.50 * (node2[1] - node2[-1]) / max_stdev;
	}
	else if (greek == GAMMA)
	{
		answer = disc_dt * (node2[1] + node2[-1] - 2 * node2[0]) / (maxvar * dt);
	}
	else
		answer = UNKNOWN_GREEK;

	free_dvector (node1, -num_steps, num_steps);
	free_dvector (node2, -num_steps, num_steps);

	return answer;
}

#undef SWITCH
