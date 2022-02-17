#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "math.h"

static double Z_Func(double x, double t1, double t2)
{
	return ((exp(x * t2) - exp(x * t1))/ x);
}

static double Phi_Func(double x, double T, double s, double t)
{
double result;

	result = (exp (-x * (T - t)) - exp (-x * (T - s))) / x;

	return result;
}

static double Etha_Func(double x, double T, double s, double t)
{
double result;

	result = (t - s - Phi_Func (x, T, s, t)) / x;

	return result;
}

static double B_Func(double x, double T, double s, double t)
{
double result;

	result = - (t * t - s * s) / 2.0 
		+ 1.0 / x * ((t - 1.0 / x) * exp (-x * (T - t)) 
		- (s - 1.0 / x) * exp(-x * (T - s)));

	return result;
}

static double Psi_Func(double x, double y, double T, double s, double t)
{
double result;

	result = 1.0 / (x * y) 
		* (t - s 
		- Phi_Func (x, T, s, t) 
		- Phi_Func (y, T, s, t) 
		+ Phi_Func (x + y, T, s, t));

	return result;
}

static double Psi2_Func(double x, double y, double Tx, double Ty, double s, double t)
{
double result;

	result = 1.0 / (x * y) 
		* (t - s 
		- Phi_Func(x, Tx, s, t) 
		- Phi_Func(y, Ty, s, t) 
		+ exp (-x * (Tx - Ty)) * Phi_Func (x + y, Ty, s, t));

	return result;
}


Err fill_mc_multi_init_corr(
					long		pay_date,
					double		pay_time,
					long		*date,
					double		*time,
					long		nb_dates,
					LINK_UND	link					
					)
{
double	var, expect, phi;
double	*sig_curve_dom;
double	lda_dom, sig_dom, fx_vol;
char*	dom_yc;
int		index_dom, index_for, dom_index, index_fx_quanto;
double	sig_und, sig_und2, sig_und_dom, sig_und_for;
double	sig_l_dom, sig_l_for, sig_l;
double	lda_und;
double	x_und, y_und;
double	QTexpect;

double	adj_pay, adj_quanto;
double	mat_pay, zc_pay, zc_dom, zc_for, pay_mat, mat;
double  ***do_covar = NULL;
int		i, k, und, l;
long	StartIndex, EndIndex;
double	start_date, end_date, start_mat, end_mat, T1, T2;
double	*correl = NULL,
		*correl_1 = NULL,
		*correl_2 = NULL,
		*correl_3 = NULL,
		**correl_fxfx = NULL;

Err		err = NULL;

	do_covar = f3tensor(0, nb_dates-1, 0, link->num_und-1, 0, link->num_und-1);
	correl = dvector(0, link->num_und-1);
	correl_1 = dvector(0, link->num_und-1);
	correl_2 = dvector(0, link->num_und-1);
	correl_3 = dvector(0, link->num_und-1);
	correl_fxfx = dmatrix(0, link->num_und-1, 0, 8);

	if (!do_covar || !correl || !correl_1 || !correl_2 || !correl_3 || !correl_fxfx)
	{
		err = "Memory allocation failure in fill_mc_multi_init_corr";
		goto FREE_RETURN;
	}

	/* Domestic informations */
	dom_index = link -> dom_index;
	sig_curve_dom = link -> sig_curve[dom_index];
	lda_dom = link -> lambda[dom_index];
	dom_yc = link -> yc[dom_index];
	
	mat_pay = (pay_date - date[0]) / 365.0;
	
	/* First we do the LGM Underlyings */
	for (und=0; und<link -> num_und; und++)
	{
		if (link -> type[und] < 2)
		{			
			phi = 0;
			lda_und = link -> lambda[und];
			
			start_date = date[0];
			start_mat = time[0];
			StartIndex = Get_Index (start_mat, link -> sig_dates, link -> nb_sig_dates);

			for (k=0; k<nb_dates-1; k++)
			{
				end_date = date[k+1];					
				end_mat = time[k+1];
				EndIndex = Get_Index (end_mat, link -> sig_dates, link -> nb_sig_dates);
				mat = (end_mat - start_mat);

				/* Initialisation */
				var = 0;
				expect = 0;		
				adj_pay = 0;
				adj_quanto = 0;
						
				for (l=0; l<link -> num_und; l++)
				{
					correl_1[l] = correl_2[l] = correl_3[l] = correl[l] = 0.0;								
				}
						
				link -> beta[und][k] = -1.0 / lda_und * (1.0 - exp(-lda_und * mat));
				
				if (und == link->dom_index)
				{
					zc_pay = swp_f_zr (start_date, pay_date, dom_yc);
					pay_mat = (pay_date - start_date) / 365.0;
					link -> dom_bond_pay[k] = exp(-zc_pay * pay_mat);
					link -> dom_beta_pay[k] = -1.0 / lda_dom * (1 - exp(-lda_dom * pay_mat));
				}														
				for (i = StartIndex; i < EndIndex + 1; i++)
				{
					if (i > StartIndex)
					{
						T1 = link -> sig_dates[i-1];
					}
					else
					{
						/* First part */
						T1 = start_mat;
					}
					

					if (i == EndIndex || StartIndex == EndIndex)
					{
						/* Last part */
						T2 = end_mat;
					}
					else
					{
						T2 = link -> sig_dates[i];
					}

					sig_und = sig_und2 = link -> sig_curve[und][i];
					sig_und2 *= sig_und2;			
					
					x_und = Z_Func(lda_und, T1, T2);
					y_und = Z_Func(2 * lda_und, T1, T2);
					
					expect += sig_und2 * (Phi_Func(lda_und, end_mat, T1, T2) 
											 -Phi_Func(2 * lda_und, end_mat, T1, T2));
					
					if (link -> type[und] == 0)
					{
						/* Q-Tpay adjustment */

						adj_pay += sig_und * link -> sig_curve[dom_index][i]
									* (x_und - exp(-lda_dom * mat_pay) * Z_Func(lda_dom + lda_und, T1, T2));
					}
					else	/* link -> type[und] == 1 */
					{
						adj_pay += sig_und * link -> sig_curve[dom_index][i]
									* (x_und - exp(-lda_dom * mat_pay) * Z_Func(lda_dom + lda_und, T1, T2))
									* link -> correlations[und][dom_index][i];

						/* Quanto Adjustment in the foreign case */

						adj_quanto += sig_und * link -> sig_curve[link -> fx_index[und]][i] * x_und
							* link -> correlations[und][link -> fx_index[und]][i];
					}
					
					var += sig_und2 * y_und;				
											
					/* Correlations */
					for (l=0; l<link -> num_und; l++)
					{
						if (!do_covar[k][und][l] &&  und != l)
						{
							switch (link -> type[l])
							{
								case 0:
								case 1:
								
								/* LGM1 / LGM2 */
								correl[l] += sig_und * link -> sig_curve[l][i] 
											* Z_Func(lda_und + link -> lambda[l], T1, T2)
											* link -> correlations[und][l][i];
								break;

								case 2:
								/* LGM / fx */
								sig_l_dom = link -> sig_curve[link -> dom_forex[l]][i];
								sig_l_for = link -> sig_curve[link -> for_forex[l]][i];

								correl_1[l] += sig_und * sig_l_dom * (									
										x_und - exp(-link -> lambda[link -> dom_forex[l]] * end_mat) 
										* Z_Func(lda_und + link -> lambda[link -> dom_forex[l]], T1, T2))
										* link -> correlations[und][link->dom_forex[l]][i];
												
								correl_2[l] += sig_und * sig_l_for * (									
										x_und - exp(-link -> lambda[link -> for_forex[l]] * end_mat) 
										* Z_Func(lda_und + link -> lambda[link -> for_forex[l]], T1, T2))
										* link -> correlations[und][link->for_forex[l]][i];
								 
								correl_3[l] += sig_und * link -> sig_curve[l][i] * x_und
										* link -> correlations[und][l][i];
								break;
							}
						}
					}
								
				}		
				phi += var;		

				/* Forward and Standard deviation */
				if (link -> type[und] == 0)
				{
					link -> fwd[und][k+1] = 1.0 / lda_und * (expect - exp(-lda_und * end_mat) * adj_pay)
							   +  link -> phi[und][k]  * exp(-lda_und * mat) 
								   * Phi_Func(-lda_und, start_mat, start_mat, end_mat);
				}
				else
				{
					link -> fwd[und][k+1] = 1.0 / lda_und * expect
						+ exp(-lda_und * end_mat) *
						(- 1.0 / lda_dom * adj_pay
						 - 1.0 * adj_quanto)
							   +  link -> phi[und][k]  * exp(-lda_und * mat) 
								   * Phi_Func(-lda_und, start_mat, start_mat, end_mat);
				}


				link -> std[und][k+1] = exp(-lda_und * end_mat) * sqrt(var);
				
				/* Covariance*/
				for (l=0; l<link -> num_und; l++)
				{
					if (!do_covar[k][und][l] && und != l)
					{
						switch (link -> type[l])
						{
							case 0:
							case 1:
							/* LGM1 / LGM2 */			
							link -> covariance[k+1][und][l] = 
								1.0 * correl[l] * exp(-(lda_und + link -> lambda[l]) * end_mat);
											
							break;		
							
							case 2:
							/* LGM / FX */
							link -> covariance[k+1][und][l] = 
								exp(-lda_und * end_mat) * 
									(1.0 * correl_1[l] / link -> lambda[link -> dom_forex[l]]
									- 1.0  * correl_2[l] / link -> lambda[link -> for_forex[l]]
									+ 1.0 * correl_3[l]);

							break;							
						}		
						link -> covariance[k+1][l][und] = link -> covariance[k+1][und][l];
						do_covar[k][und][l] = do_covar[k][l][und] = 1.0;
					}
				}
				link -> covariance[k+1][und][und] = link -> std[und][k+1] * link -> std[und][k+1];

				link -> phi[und][k+1] = phi * exp(-2 * lda_und * end_mat);

				start_date = end_date;
				start_mat = end_mat;
				StartIndex = EndIndex;
			}
						
			if (und == link->dom_index)
			{
				zc_pay = swp_f_zr (start_date, pay_date, dom_yc);
				pay_mat = (pay_date - start_date) / 365.0;
				link -> dom_bond_pay[k] = exp(-zc_pay * pay_mat);
				link -> dom_beta_pay[k] = -1.0 / lda_dom * (1 - exp(-lda_dom * pay_mat));
			}
		}
	}

	/* now we have to do the Fx cases */
	for (und=0; und<link -> num_und; und++)
	{
		if (link -> type[und] == 2)
		{									
			start_date = date[0];
			start_mat = time[0];
			StartIndex = Get_Index (start_mat, link -> sig_dates, link -> nb_sig_dates);

			index_dom = link -> dom_forex[und];
			index_for = link -> for_forex[und];

			if (index_dom != dom_index)
			{
				/* Need to find the fx between dom_index and index_dom */
				index_fx_quanto = link->fx_index[index_dom];
			}
			else
			{
				index_fx_quanto = -1;
			}

			for (k=0; k<nb_dates-1; k++)
			{
				end_date = date[k+1];					
				end_mat = time[k+1];
				EndIndex = Get_Index (end_mat, link -> sig_dates, link -> nb_sig_dates);
				mat = (end_mat - start_mat);

				/* Initialisation */
				var = 0;
				expect = 0;		
				adj_pay = 0;
				adj_quanto = 0;

				for (l=0; l<link -> num_und; l++)
				{
					correl_fxfx[l][0] = correl_fxfx[l][1] = correl_fxfx[l][2] = 0.0;
					correl_fxfx[l][3] = correl_fxfx[l][4] = correl_fxfx[l][5] = 0.0;
					correl_fxfx[l][6] = correl_fxfx[l][7] = correl_fxfx[l][8] = 0.0;
				}
						
				zc_dom = swp_f_zr (start_date, end_date, link -> yc[index_dom]);
				zc_for = swp_f_zr (start_date, end_date, link -> yc[index_for]);
									
				/* QTexpect of the log Fx */
				QTexpect = -mat * (zc_for-zc_dom) 
		    		 - 0.5 * (link -> beta[index_for][k] * link -> beta[index_for][k]  * link -> phi[index_for][k]
	    			 - link -> beta[index_dom][k] * link -> beta[index_dom][k]  * link -> phi[index_dom][k]);	
			
				/*	Implied Fx Vol */
				err = Fx3DtsImpliedVol_corr (	end_mat, start_mat, end_mat,
										link -> sig_dates,  link -> nb_sig_dates,
										link -> sig_curve[index_dom], link -> lambda[index_dom],
										link -> sig_curve[index_for], link -> lambda[index_for],
										link -> sig_dates,  link -> sig_curve[und],  link -> nb_sig_dates, 
										link -> sig_dates,
										link -> correlations[index_dom][index_for],
										link -> correlations[index_dom][und],
										link -> correlations[und][index_for],
										link -> nb_sig_dates,	&fx_vol);
				if (err)
					goto FREE_RETURN;

				/*	Calculate expectation of the log Fx under Q-Tfix */
				expect = QTexpect - 0.5 * fx_vol * fx_vol * mat;
				
				for (i = StartIndex; i < EndIndex + 1; i++)
				{
					if (i > StartIndex)
					{
						T1 = link -> sig_dates[i-1];
					}
					else
					{
						/* First part */
						T1 = start_mat;
					}
					
					if (i == EndIndex || StartIndex == EndIndex)
					{
						/* Last part */
						T2 = end_mat;
					}
					else
					{
						T2 = link -> sig_dates[i];
					}

					sig_und = sig_und2 = link -> sig_curve[und][i];
					sig_und2 *= sig_und2;
					sig_dom = link -> sig_curve[dom_index][i];
								
					/* Q-Tpay adjustment */
					
					adj_pay += link -> correlations[dom_index][und][i]
										* sig_dom * sig_und * (Etha_Func (lda_dom, end_mat, T1, T2)
															  -Etha_Func (lda_dom, mat_pay, T1, T2))
					   - link -> correlations[dom_index][index_for][i]
					   * sig_dom * link -> sig_curve[index_for][i] *
					   (Psi_Func (link -> lambda[index_for], lda_dom, end_mat, T1, T2) 
						 -Psi2_Func (link -> lambda[index_for], lda_dom, end_mat, mat_pay, T1, T2))
					   + link -> correlations[dom_index][index_dom][i]
					   * sig_dom * link -> sig_curve[index_dom][i] *
					   (Psi_Func (link -> lambda[index_dom], lda_dom, end_mat, T1, T2) 
						 -Psi2_Func (link -> lambda[index_dom], lda_dom, end_mat, mat_pay, T1, T2));  

					/* Quanto Adjustment in the foreign case */
					if (dom_index != index_dom)
					{	
						/* NOT VALID */
						/*
						adj_quanto += sig_und * link -> sig_curve[link -> fx_index[index_dom]][i] * (T2 - T1);
						*/
					}
																								
					/* Correlations */
					for (l=0; l<link -> num_und; l++)
					{
						if (!do_covar[k][und][l] &&  und != l)
						{
							switch (link -> type[l])
							{							
								case 2:
								
								sig_l_dom = link -> sig_curve[link -> dom_forex[l]][i];								
								sig_l_for = link -> sig_curve[link -> for_forex[l]][i];
								sig_l = link -> sig_curve[l][i];

								sig_und_dom = link -> sig_curve[link -> dom_forex[und]][i];								
								sig_und_for = link -> sig_curve[link -> for_forex[und]][i];

								/* Z_und with Z_l */
								correl_fxfx[l][0] += sig_und * sig_l * (T2 - T1)
									* link -> correlations[und][l][i];

								/* Z_und with Dom_l */
								correl_fxfx[l][1] += sig_und * sig_l_dom * (Etha_Func(link -> lambda[link -> dom_forex[l]], end_mat, T1, T2))
									* link -> correlations[und][link -> dom_forex[l]][i];

								/* Z_und with For_l */
								correl_fxfx[l][2] += sig_und * sig_l_for * (Etha_Func(link -> lambda[link -> for_forex[l]], end_mat, T1, T2))
									* link -> correlations[und][link -> for_forex[l]][i];

								/* Dom_und with Z_l */
								correl_fxfx[l][3] += sig_und_dom * sig_l * (Etha_Func(link -> lambda[link -> dom_forex[und]], end_mat, T1, T2))
									* link -> correlations[link -> dom_forex[und]][l][i];

								/* Dom_und with Dom_l */
								correl_fxfx[l][4] += sig_und_dom * sig_l_dom * (Psi_Func(link -> lambda[link -> dom_forex[und]], link -> lambda[link -> dom_forex[l]], end_mat, T1, T2))
									* link -> correlations[link -> dom_forex[und]][link -> dom_forex[l]][i];

								/* Dom_und with For_l */
								correl_fxfx[l][5] += sig_und_dom * sig_l_for * (Psi_Func(link -> lambda[link -> dom_forex[und]], link -> lambda[link -> for_forex[l]], end_mat, T1, T2))
									* link -> correlations[link -> dom_forex[und]][link -> for_forex[l]][i];

								/* For_und with Z_l */
								correl_fxfx[l][6] += sig_und_for * sig_l * (Etha_Func(link -> lambda[link -> for_forex[und]], end_mat, T1, T2))
									* link -> correlations[link -> for_forex[und]][l][i];

								/* For_und with Dom_l */
								correl_fxfx[l][7] += sig_und_for * sig_l_dom * (Psi_Func(link -> lambda[link -> for_forex[und]], link -> lambda[link -> dom_forex[l]], end_mat, T1, T2))
									* link -> correlations[link -> for_forex[und]][link -> dom_forex[l]][i];

								/* For_und with For_l */
								correl_fxfx[l][8] += sig_und_for * sig_l_for * (Psi_Func(link -> lambda[link -> for_forex[und]], link -> lambda[link -> for_forex[l]], end_mat, T1, T2))
									* link -> correlations[link -> for_forex[und]][link -> for_forex[l]][i];
								
								break;
							}
						}
					}
				}												

				/* Standard deviation */
				link -> std[und][k+1] = fx_vol * sqrt(mat);
				
				/* Covariance*/
				for (l=0; l<link -> num_und; l++)
				{
					if (!do_covar[k][und][l] &&  und != l)
					{
						switch (link -> type[l])
						{														
							case 2:

							link -> covariance[k+1][und][l] = 
				  1.0 * correl_fxfx[l][0]
				+ 1.0 * correl_fxfx[l][1]
				- 1.0 * correl_fxfx[l][2]
				+ 1.0 * correl_fxfx[l][3]
				+ 1.0 * correl_fxfx[l][4]
				- 1.0 * correl_fxfx[l][5]
				- 1.0 * correl_fxfx[l][6]
				- 1.0 * correl_fxfx[l][7]
				+ 1.0 * correl_fxfx[l][8];

							break;							
						}		
						link -> covariance[k+1][l][und] = link -> covariance[k+1][und][l];
						do_covar[k][und][l] = do_covar[k][l][und] = 1.0;
					}
				}

				link -> covariance[k+1][und][und] = link -> std[und][k+1] * link -> std[und][k+1];

				/* Quanto adjustment */
				if (index_fx_quanto >= 0)
				{
					adj_quanto = -link -> covariance[k+1][und][index_fx_quanto];
				}

				/* Forward and Standard deviation */
				link -> fwd[und][k+1] = expect + adj_pay + adj_quanto;
			
				start_date = end_date;
				start_mat = end_mat;
				StartIndex = EndIndex;
			}						
		}
	}

FREE_RETURN:

	if (do_covar)
		free_f3tensor(do_covar, 0, nb_dates-1, 0, link->num_und-1, 0, link->num_und-1);
	if (correl)
		free_dvector(correl, 0, link->num_und-1);
	if (correl_1)
		free_dvector(correl_1, 0, link->num_und-1);
	if (correl_2)
		free_dvector(correl_2, 0, link->num_und-1);
	if (correl_3)
		free_dvector(correl_3, 0, link->num_und-1);
	if (correl_fxfx)
		free_dmatrix(correl_fxfx, 0, link->num_und-1, 0, 8);
	
	return err;
}
