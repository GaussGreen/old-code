/*******************************************************************************
**                      Include Files                               
*******************************************************************************/
#include        "utallhdr.h"
#include        <OPFNCTNS.H>
#include		<math.h"

/*******************************************************************************
*                                                                          
* FUNCTION     	: srt_f_nono_bar                              
*                                                                          
* PURPOSE      	: Calculates the premium and the greeks for a scud options
*		  with multiple barriers 
*                                                                          
* DESCRIPTION  	: based on the nono-tree srt_f_nono_tree developed in SORT equity  
*	  	  by R. Carson 15/05/95 - originally in the file srt_f_nono_bar.c
*									    
* PARAMETERS   
*	INPUT : spot1		- spot value of 1st underlying 
*			spot2		- spot value of 2nd underlying   
*			fwd1 		- forward value of 1st underlying               
*			fwd2  		- forward value of 2nd underlying              
*			strike		- strike of option        
*			n_days		- numcer of days in tree           
*			vol1		- volatility of 1st underlying              
*			vol2		- volatility of 2nd underlying             
*			rho		- correlation between two underlyings              
*			barrier[]	- array of 4 barriers
*	 		disc_fac	- discount factor          
*			n_steps	- number of steps in tree         
*			call_put	- call/put          
*			info_wanted 	- info wanted (premium, greeks) 
*                                                                           
* RETURNS      	: premium, delta, gamma                                         
*                                                                           
*******************************************************************************/


/*******************************************************************************
*****  STRUCTURES 
*******************************************************************************/

typedef struct
{                    
	double	E_1; 
	double	E_S; 
	double	E_SS; 
	double	E_R; 
	double	E_RR; 
	double	E_RS;           
} MOMENTS;             


static double q2(double alpha,double a,double b)
{
double result;

	result  = exp(a*alpha + 0.5*b*b*alpha*alpha);
	result *= norm((a + alpha*b*b) / b);

	return(result);
}                                                       



/************************************************************************************************/
/*																								*/
/* IDEA		: Calculates the 1st 2 absorbed moments of a brownian motion (Barriers)				*/
/*																								*/
/************************************************************************************************/
   
static void absorbed_moment(double	spot,
							double 	barrier,
							double 	dt,             
							double 	vol,                       
							double 	rate,            
							double 	absorb_moments[2])
 		       	
{                        
double	m;  
double	mu;                                                                     
double	m1;  
double	m2;  
double	ft; 
double	fe; 
double	absorb_prob;  
double	absorb_moment[2];
long	up_down;

	up_down = ((spot > barrier) ? -1 : 1);

	vol	= DMIN(vol,1);
	vol	= DMAX(vol,0.0001);
          
	m  = log (barrier / spot);

	if (fabs(m) < 1e-6)             /* Send all paths to barrier, no vol.*/
	{
		absorb_moments[0] = barrier / spot;	
		absorb_moments[1] = absorb_moments[0] * absorb_moments[0] ;

		return;
	}
										/*       / Su^i+1	*/
	mu = (rate - vol * vol / 2);		/* L____/_    		*/
	m1 = (m - mu * dt) * up_down;		/*     /      		*/
	m2 = - (m + mu * dt) * up_down;		/*    /___ Su^i		*/
	ft = vol * sqrt(dt);				/*    \    			*/
	fe = exp (2 * mu * m / vol / vol);	/*     \   			*/
	absorb_prob = norm (m1 / ft)	   	/*      \  			*/
				  - fe * norm ( m2 /ft);/*       \ su^i-1	*/

	absorb_moment[0]  = barrier 
					    * (q2 (- up_down, m1,ft) - fe * q2(-up_down, m2, ft));   
	absorb_moment[0] /= absorb_prob;

	absorb_moments[0] = absorb_moment[0] * absorb_prob + barrier * (1 - absorb_prob);  
	absorb_moments[0] /= spot;                                          
                                               

	absorb_moment[1] =	barrier * barrier * (q2 (-2*up_down, m1, ft) 
						-  fe * q2 (-2*up_down, m2, ft))
						/ absorb_prob;  

	absorb_moments[1] = absorb_moment[1] * absorb_prob
						+ barrier * barrier * (1 - absorb_prob); 
    absorb_moments[1] /= spot * spot;                  

	if (absorb_moments[0]*absorb_moments[0] > absorb_moments[1])	/* Send all paths to barrier, no vol.*/
	{                                            
		absorb_moments[0] = barrier / spot;  
		absorb_moments[1] = absorb_moments[0] * absorb_moments[0] ;
	}
	
	return;
}



/* + + + + + + + + + + + + + + + + + + + + + + + + + + + */
static void	quick_live_val	(double		S[],
							 double		R[],
							 MOMENTS	*moments,
							 double		*val_u,    
							 double 	*val_m,    
							 double 	*val_d,     
							 double 	*new_val,
							 double		rho)
{       
double f0, fS, fSS, fr, frr, fSr;

double	E_F;                 
double *vals[3];
                                        
double	E_S;
double	E_SS;
double	E_R;
double	E_RR;
double	E_RS;
              
double	su;
double	sd;
double	sud;
double  f_uu,f_ud,f_du,f_dd;
double	ru;
double	rd;
double	rud;

                               
	E_S	 = moments->E_S;
	E_SS = moments->E_SS;
	E_R	 = moments->E_R;
	E_RR = moments->E_RR;
	E_RS = moments->E_RS;

	vals[0]	= val_d;
	vals[1]	= val_m;
	vals[2]	= val_u;
                             
	ru	= 1.0 / (R[2] - R[1]);
	rd	= 1.0 / (R[0] - R[1]);
	rud	= 1.0 / (R[2] - R[0]);

	su	= 1.0 / (S[2] - S[1]);
	sd	= 1.0 / (S[0] - S[1]);
	sud	= 1.0 / (S[2] - S[0]);
	
	f0	= vals [1][1];				
	fS	= 0;				
	fr	= 0;				
	fSS	= 0;				
	fSr	= 0;				
	frr	= 0;				

	if (sud != 1.0)
	{
		fSS	-=  vals [1][0] * sd * sud;
		fSS	+=  vals [1][1] * su * sd;
		fSS	+=  vals [1][2] * sud * su;
	}

	if (R[2] != R[0])
	{
		frr	-=  vals [0][1] * rd * rud;
		frr	+=  vals [1][1] * ru * rd;
		frr	+=  vals [2][1] * rud * ru;
	}		

	fS	 = (vals[1][0] - f0 - fSS / 2.0 * (S[0] - S[1]) * (S[0] - S[1]) )
		   * sd;
	fr	 = (vals[0][1] - f0 - frr / 2.0 * (R[0] - R[1]) * (R[0] - R[1]) )
		   * rd;

	if ( (S[2] != S[0]) && (R[2] != R[0]) )
	{
		f_uu = (vals[2][2] - vals[2][1] - vals[1][2] + vals[1][1]) 
				* ru * su ;
        f_du = (vals[0][2] - vals[0][1] - vals[1][2] + vals[1][1]) 
				* rd * su ;
        f_dd = (vals[0][0] - vals[0][1] - vals[1][0] + vals[1][1]) 
				* rd * sd ;
        f_ud = (vals[2][0] - vals[2][1] - vals[1][0] + vals[1][1]) 
				* ru * sd ;
		fSr = 0.25 *( (1 - rho) * (f_du + f_ud)  
			  + (1 + rho) * (f_uu + f_dd) ) ;
    }

	E_F	= f0 + fS * E_S + fr * E_R + fSS * E_SS  + frr * E_RR  + fSr * E_RS;

	*new_val = E_F;

	return;
}        


/* + + + + + + + + + + + + + + + + + + + + + + + + + + + */
static void	moments_stock_stock	(	double	fwd_S,
									double	nearest_S,
									double	fwd_S2,
									double	nearest_S2,
									double	vol_s,
									double	vol_s2,
									double	rho,
									double	dt,
									MOMENTS	*moments)
{
	moments->E_1	= 1;

	moments->E_S	= fwd_S - nearest_S;
	
	moments->E_SS	= fwd_S * fwd_S * exp(vol_s * vol_s * dt) 
					  - 2 * fwd_S * nearest_S 
					  + nearest_S * nearest_S;
	
	moments->E_R	= fwd_S2 - nearest_S2;

	moments->E_RR	= fwd_S2 * fwd_S2 *  exp(vol_s2 * vol_s2 * dt) 
					  - 2 * fwd_S2 * nearest_S2 
					  + nearest_S2 * nearest_S2;

	moments->E_RS	= fwd_S * fwd_S2 * exp (vol_s * vol_s2 * rho * dt)
					  - fwd_S  * nearest_S2
					  - fwd_S2 * nearest_S
					  + nearest_S * nearest_S2;
}


/* + + + + + + + + + + + + + + + + + + + + + + + + + + + */

static Err	init_nono_tree (double	**R,
		       				double	**S,
		       				double	**new_R,
							double	**new_S,
							double	**fwd_S,
							double	**mid_R,
							double	**gamma,
							double	**gamma_prime,
							double	**df_cve,
							double	**stoch_rate,
							double	**fwd_rate,
							double	***val,
							double	***new_val,
							long	n_steps)
{
	*mid_R			= dvector (0, n_steps+1);
	*gamma			= dvector (0, n_steps+1);
	*gamma_prime	= dvector (0, n_steps+1);
	*fwd_S			= dvector (0, n_steps+1);
	*df_cve			= dvector (0, n_steps+1);
	*stoch_rate		= dvector (0, n_steps+1);
	*fwd_rate		= dvector (0, n_steps+1);
	*R				= dvector(-n_steps-1, n_steps+1);
	*S				= dvector(-n_steps-1, n_steps+1);
	*new_R			= dvector(-n_steps-1, n_steps+1);
	*new_S			= dvector(-n_steps-1, n_steps+1);                                                        
	*val			= dmatrix(-n_steps-1, n_steps+1, -n_steps-1, n_steps+1);
	*new_val		= dmatrix(-n_steps-1, n_steps+1, -n_steps-1, n_steps+1);

	if  (  (*mid_R == NULL)
		||  (*gamma == NULL)
		||  (*gamma_prime == NULL)
		||  (*fwd_S == NULL)
		||  (*df_cve == NULL)
		||  (*stoch_rate == NULL)
		||  (*fwd_rate == NULL)
		||  (*R == NULL)
		||  (*S == NULL)
		||  (*new_R == NULL)
		||  (*new_S == NULL)
		||  (*val == NULL)
		||  (*new_val == NULL))
	{
		return serror ("Memory allocation failure");
	}

	return NULL;
}

static void free_all (	double	**R,
		       			double	**S,
		       			double	**new_R,
						double	**new_S,
						double	**fwd_S,
						double	**mid_R,
						double	**gamma,
						double	**gamma_prime,
						double	**df_cve,
						double	**stoch_rate,
						double	**fwd_rate,
						double	***val,
						double	***new_val,
						long	n_steps)
{
	free_dvector (*mid_R, 0, n_steps+1);
	free_dvector (*gamma, 0, n_steps+1);
	free_dvector (*gamma_prime, 0, n_steps+1);
	free_dvector (*fwd_S, 0, n_steps+1);
	free_dvector (*df_cve, 0, n_steps+1);
	free_dvector (*stoch_rate, 0, n_steps+1);
	free_dvector (*fwd_rate, 0, n_steps+1);
	free_dvector(*R, -n_steps-1, n_steps+1);
	free_dvector(*S, -n_steps-1, n_steps+1);
	free_dvector(*new_R, -n_steps-1, n_steps+1);
	free_dvector(*new_S, -n_steps-1, n_steps+1);                                                        
	free_dmatrix(*val, -n_steps-1, n_steps+1, -n_steps-1, n_steps+1);
	free_dmatrix(*new_val, -n_steps-1, n_steps+1, -n_steps-1, n_steps+1);
}



double 	srt_f_nono_bar (double		spot1,
						double		spot2,
						double		fwd1,
						double		fwd2,
						double		strike,   
						double		mat ,
						double		vol1,
						double		vol2,
						double		rho,
						double		barrier[],
						double		disc_fac,
						long		n_steps,
						SrtCallPutType	call_put,     
						SrtGreekType	greek)    
                       
{     
/* Declarations */
                                          
long	idx_S; 
long	idx_R;                    	
long	time_step;
                         
long	smooth_R;
long	smooth_S;
                      
long	barr11;
long	barr12;
long	barr21;
long	barr22;

double	*mid_R;     
double	*gamma;
double	*gamma_prime;

double	*S;
double	*R;
double	*swapper;  
double	*new_S;
double	*new_R;
                          
double	*fwd_S;   
double	*df_cve;  
double	*stoch_rate;  
double	*fwd_rate;   
  
double	**val;
double	**new_val; 
double	**swap_ptr;

double	temp_R[3];
double	temp_S[3];
double	absorb_moments[2];
                       
double	dlt_lat_up; 
double	dlt_lat_dn;       
double	delta_up;
double	delta_dn;

double  dlt_lat_up1;
double	dlt_lat_dn1;   
double	delta_up1;
double	delta_dn1;

double	current_time;

double	dt;
double	rt_dt;
                           
double	u_s;       
double	u_r;

double	ff1;
double	ff2;
double	fwd_r;
double	fwd_s;
double	DF;
                      
double	vol_R;
double	vol_S;

double	rate1;
double	rate2;

double  answer ;

MOMENTS	moments;
long	err_sts	= 0 /* = M_SUCCESS */ ;
long	n_days;   


	n_days = (long) mat * 365 ; 


	/* Memory Allocation */
	if ( init_nono_tree (&R, &S, &new_R, &new_S, &fwd_S,
						&mid_R, &gamma, &gamma_prime, &df_cve,
						&stoch_rate, &fwd_rate, &val, &new_val, n_steps)!= NULL)
	{
		free_all (&R, &S, &new_R, &new_S, &fwd_S,
				 &mid_R, &gamma, &gamma_prime, &df_cve,
				 &stoch_rate, &fwd_rate, &val, &new_val, n_steps);

		return MEMORY_ERR;
	}
              
	barr11 = (barrier[0] != 0.0);
	barr12 = (barrier[1] != 0.0);
	barr21 = (barrier[2] != 0.0);
	barr22 = (barrier[3] != 0.0);

	current_time	= 0;
	dt				= n_days / (double) (n_steps * 365.0);
	rt_dt			= sqrt( dt);
	u_r				= exp( vol2  * rt_dt * sqrt(3));
	u_s				= exp( vol1  * rt_dt * sqrt(3));
	DF  			= exp (log (disc_fac)/ (double)  n_steps);
	rate1			= log (fwd1 / spot1) / ((double) n_days / 365.0); 
	rate2			= log (fwd2 / spot2) / ((double) n_days / 365.0); 
	ff1				= exp (rate1 * dt);
	ff2				= exp (rate2 * dt);

    current_time	= (double) n_days / 365.0;					/* Loop should end here anyway */

	S[-n_steps - 1]	= fwd1 * pow (u_s, - (n_steps + 1));
	R[-n_steps - 1]	= fwd2 * pow (u_r, - (n_steps + 1));

	for (idx_R = -n_steps; idx_R <= n_steps + 1; idx_R++)
	{
			S[idx_R]	= S[idx_R - 1] * u_s;
		new_S[idx_R]	= S[idx_R];
			R[idx_R]	= R[idx_R - 1] * u_r;
		new_R[idx_R]	= R[idx_R];
	}
  
	for (idx_R = - (n_steps + 1); idx_R <= (n_steps + 1); idx_R++)
	{
		for (idx_S = -(n_steps + 1); idx_S <= (n_steps + 1); idx_S++)
		{
			switch (call_put)
			{
				case SRT_CALL:		/* call */
					val[idx_R][idx_S] = DMAX ( (S[idx_S] - strike) , 0.0);       
				break;

				case SRT_PUT:		/* put */
					val[idx_R][idx_S] = DMAX ( ( strike - S[idx_S] ) , 0.0);       
				break;

				/* + + + + + + + + + + + + + + + + + + + + + +	*/
				/*	possible to have SPREAD here etc   	        */
				/*  + + + + + + + + + + + + + + + + + + + + + +	*/
			}
		}
	}

	for (time_step	= n_steps - 1; time_step >= 0; time_step --)
	{                                                                 
		current_time -= dt;

		for (idx_R = -(time_step + 1); idx_R <= (time_step + 1); idx_R++)
		{
			new_S[idx_R]	= S[idx_R] / ff1;
       		new_R[idx_R]	= R[idx_R] / ff2;				
        }

		for (idx_R = -(time_step + 1); idx_R <= (time_step + 1); idx_R++)
		{                                                                                  
			/* + + + + + + + + + + + + + + + + + + + + + + + + + +	*/
			/*	          Loop on Spot to follow		            */
			/*  + + + + + + + + + + + + + + + + + + + + + + + + + +	*/

			if ((barr21) && (new_R[idx_R] < barrier[2]))
			{
				new_val[idx_R][idx_S] = 0.0;
				smooth_R = 0;
			}   
			else
			if ((barr22) && (R[idx_R] > barrier[3]))
			{
				new_val[idx_R][idx_S] = 0.0;
				smooth_R = 0;
			}   
			else
			if ((barr21) && (R[idx_R - 1] < barrier[2]))
			{
				absorbed_moment(new_R[idx_R]	,
								barrier[2]	,
								dt		,             
								vol2,
			                    rate2		,            
								absorb_moments	);

				fwd_r	= new_R[idx_R] * absorb_moments[0];
				vol_R	= sqrt (log (absorb_moments[1] / absorb_moments[0] / absorb_moments[0]) / dt);

				temp_R [0] = barrier[2];
				temp_R [1] = R[idx_R];
				temp_R [2] = R[idx_R + 1];

				for (idx_S = - time_step-2; idx_S <= (time_step + 2); idx_S++)
				{
					val[idx_R - 1][idx_S]  = 0;
                }

				smooth_R = 1;
            }
			else
			if ((barr22) && (R[idx_R + 1] > barrier[3]))
			{
				absorbed_moment(new_R[idx_R]	,
								barrier[3]	,
								dt		,             
								vol2,
			                    rate2		,            
								absorb_moments	);

				fwd_r	= new_R[idx_R] * absorb_moments[0];
				vol_R	= sqrt (log (absorb_moments[1] / absorb_moments[0] / absorb_moments[0]) / dt);

				temp_R [2] = barrier[3];
				temp_R [1] = R[idx_R];
				temp_R [0] = R[idx_R - 1];

				for (idx_S = - time_step-2; idx_S <= (time_step + 2); idx_S++)
				{
					val[idx_R + 1][idx_S]  = 0;
                }

				smooth_R = 1;
            }
			else
			{                          
				temp_R [0] = R[idx_R - 1];
				temp_R [1] = R[idx_R];
				temp_R [2] = R[idx_R + 1];

				fwd_r	= R[idx_R];
				vol_R	= vol2;
				smooth_R = 1;
            }

			if (smooth_R)
			{
				for (idx_S = -(time_step + 1); idx_S <= (time_step + 1); idx_S++)
				{                                                  
					/* + + + + + + + + + + + + + + + + + + + + + + + + + +	*/
					/*	      Recombination Algorithm on Spot		*/
					/*  + + + + + + + + + + + + + + + + + + + + + + + + + +	*/

					if ((barr11) && (new_S[idx_S] < barrier[0]))
					{
						new_val[idx_R][idx_S] = 0.0;
						smooth_S = 0;
					}   
					else
					if ((barr12) && (S[idx_S] > barrier[1]))
					{
						new_val[idx_R][idx_S] = 0.0;
						smooth_S = 0;
					}   
					else
					if ((barr11) && (S[idx_S - 1] < barrier[0]))
					{
						absorbed_moment(new_S[idx_S]	,
										barrier[0]	,
										dt		,             
										vol1		,
					                    rate1		,            
										absorb_moments	);
		
						fwd_s	= new_S[idx_S] * absorb_moments[0];
						vol_S	= sqrt (log (absorb_moments[1] / absorb_moments[0] / absorb_moments[0]) / dt);
		
						temp_S [0] = barrier[0];
						temp_S [1] = S[idx_S];
						temp_S [2] = S[idx_S + 1];
		
						val[idx_R + 1][idx_S - 1]  = 0;
						val[idx_R - 0][idx_S - 1]  = 0;
						val[idx_R - 1][idx_S - 1]  = 0;
		
						smooth_S = 1;
		            }
					else
					if ((barr12) && (S[idx_S + 1] > barrier[1]))
					{
						absorbed_moment(new_S[idx_S]	,
										barrier[1]	,
										dt		,             
										vol1		,
					                    rate1		,            
										absorb_moments	);
		
						fwd_s	= new_S[idx_S] * absorb_moments[0];
						vol_S	= sqrt (log (absorb_moments[1] / absorb_moments[0] / absorb_moments[0]) / dt);
		
						temp_S [2] = barrier[1];
						temp_S [1] = S[idx_S];
						temp_S [0] = S[idx_S - 1];
		
						val[idx_R - 1][idx_S + 1]  = 0;
						val[idx_R + 0][idx_S + 1]  = 0;
						val[idx_R + 1][idx_S + 1]  = 0;
		
						smooth_S = 1;
		            }
					else
					{                          
						temp_S [0] = S[idx_S - 1];
						temp_S [1] = S[idx_S];
						temp_S [2] = S[idx_S + 1];
		
						fwd_s	= S[idx_S];
						vol_S	= vol1;
						smooth_S = 1;
		            }

					if (smooth_S)
					{
						moments_stock_stock (	fwd_s,
												S[idx_S],  
												fwd_r,
												R[idx_R],
												vol_S,
												vol_R,
												rho,
												dt,
												&moments);
                    
						quick_live_val (&temp_S[0],
										&temp_R[0],
										&moments,
										&val[idx_R + 1][idx_S - 1],    
										&val[idx_R][idx_S - 1],    
										&val[idx_R - 1][idx_S - 1],    
										&new_val[idx_R][idx_S],
										rho);
  
						new_val[idx_R][idx_S] *= DF;
					}
				}
			}
		}

		swapper	= S;
			S	= new_S;
        new_S	= swapper;

		swapper	= R;
			R	= new_R;
        new_R	= swapper;

		swap_ptr = new_val;
		new_val	 = val;
		val		 = swap_ptr;
	}

	dlt_lat_up	=  u_s - 1.0         ;
	dlt_lat_dn	=  1.0 - 1.0 / u_s   ;
	delta_up 	= ( val[0][1] - val[0][0] )  
				  / ( dlt_lat_up * spot1 ) ;

	delta_dn 	= ( val[0][0] - val[0][-1] )  
				  / ( dlt_lat_dn * spot1 ) ;

	dlt_lat_up1 =  u_r - 1.0         ;
	dlt_lat_dn1	= 1.0 - 1.0 / u_r   ;

	delta_up1	= ( val[1][0] - val[0][0] )  
				  / ( dlt_lat_up1 * spot2 ) ;
	delta_dn1	= ( val[0][0] - val[0][-1] )  
				  / ( dlt_lat_dn1 * spot1 ) ;

	switch( greek )
	{
		case PREMIUM :
			answer = val[0][0];
		break ;

		case DELTAX :
       		answer = (delta_up + delta_dn) / 2.0;
		break ;

		case DELTAY :
       		answer = (delta_up1 + delta_dn1) / 2.0;
		break ;

		case GAMMAX :
			answer = (delta_up - delta_dn) / (dlt_lat_up + dlt_lat_dn) * 0.02 ;
		break ;

		case GAMMAY :
			answer = (delta_up1-delta_dn1)/ (dlt_lat_up1+dlt_lat_dn1) * 0.02 ;
		break ;

		default:
			answer = UNKNOWN_GREEK;
	}


	/* Free Memory */
	free_all (&R, &S, &new_R, &new_S, &fwd_S,
			 &mid_R, &gamma, &gamma_prime, &df_cve,
			 &stoch_rate, &fwd_rate, &val, &new_val, n_steps);

	return answer ;
}


