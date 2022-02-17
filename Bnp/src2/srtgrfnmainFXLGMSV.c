#include "math.h"
#include "srt_h_all.h"
#include "Fx3FCalib.h"
#include "Fx3FUtils.h"
#include "SrtAccess.h"
#include "BGMEval.h"
#include "LGMSVUtil.h"
#include "FXLGMSVGrfn.h"
#include "FXLGMSVMC.h"
#include "FXLGMSVUnd.h"
#include "srtgrfnmainFXLGMSV.h"


char *SrtGrfnFXLGMSVMC(
					char		*underlying,
					int			numeventdates,
					long		*eventdates,
					long		tableauRows,
					long		*tableauCols,
					char		***tableauStrings,
					int			**tableauMask,
					long		auxWidth,
					long		*auxLen,
					double		**aux,

					// for Optimisation of exercise boundary 
					int			do_optimisation,
					int			*optimise,
					MCEBPARAMS	params,
					long		*resRows,

					int			nstept,
					int			numpaths,				  
					int			*nb_prod,
					double		***prod_val)
{
	int					free_str = 0;
	FIRSTAllMkts		xStr;
	SrtGrfnParam		defParm;
	GRFNPARMMCFXLGMSV	grfn_prm;
	int					forback;
	int					flag = 0;
	long				nstp;	
						
	double				*evt_tms		= NULL,						
						*time			= NULL,
						*date			= NULL,

						*dom_sigma			= NULL,
						*dom_alpha			= NULL,
						*dom_rho			= NULL,
						*dom_rho2			= NULL,						
						*dom_lameps			= NULL,
						*dom_lvleps			= NULL,
						*dom_newalphaLGM	= NULL,
						*dom_newrhoLGM		= NULL,
						
						*dom_zcvol1_star	= NULL,
						*dom_zcvol2_star	= NULL,

						*dom_dff_star		= NULL,
						*dom_gam1_star		= NULL,
						*dom_gam2_star		= NULL,
						*dom_gam1_2_star	= NULL,
						*dom_gam2_2_star	= NULL,
						*dom_gam12_star		= NULL,

						*for_sigma			= NULL,
						*for_alpha			= NULL,
						*for_rho			= NULL,
						*for_rho2			= NULL,
						*for_lameps			= NULL,
						*for_lvleps			= NULL,
						*for_newalphaLGM	= NULL,
						*for_newrhoLGM		= NULL,
						
						*for_zcvol1_star	= NULL,
						*for_zcvol2_star	= NULL,

						*for_dff_star		= NULL,
						*for_gam1_star		= NULL,
						*for_gam2_star		= NULL,
						*for_gam1_2_star	= NULL,
						*for_gam2_2_star	= NULL,
						*for_gam12_star		= NULL;

	long				*evt_dts	= NULL;
	
	int					*is_event	= NULL;
	void				**void_prm	= NULL;
		
	long				index;
	double				final_mat;

	long				today, spot_date;
						
	int					num_col		= 0,
						max_num_df	= 0,
						num_evt		= 0,
						num_und		= 0;
						
	char				*domname,
						*forname;

	char				*dom_yc, *for_yc;

	SrtUndPtr			fx_und, dom_und, for_und;
	SrtUndPtr			*und_ptr	= NULL,
						und			= NULL;

	LGMSV_model			dom_model_, *dom_model = &dom_model_;
	LGMSV_model			for_model_, *for_model = &for_model_;

	int					i, j, k, l;

	int					fx_idx, dom_idx, for_idx;

	double				dom_df;
	double				temp1, temp2;
	
	Err					err = NULL;

	double				fx_spot;

	int					num_fx_vol;
	double				*fx_time=NULL;
	double				*fx_vol=NULL;

	int					num_rho;
	double				*rho_time=NULL;
	double				***CorrMatrix=NULL;

	double				*fx_sigma=NULL;
	double				***correlation=NULL;

	//	Initialisation of the LGMSV model 
	init_NULL_LGMSV_model(dom_model);
	init_NULL_LGMSV_model(for_model);

	//	Initialise the GRFN tableau 

	//	First, initialise the param struct 	
	err = srt_f_set_default_GrfnParams (&defParm);	
	defParm.force_mc = 0;

	err = FIRSTInitMktStruct(
							numeventdates,
							eventdates,
							tableauRows,
							*tableauCols,
							tableauStrings,
							tableauMask,
							auxWidth,
							auxLen,
							aux,
							underlying,
							&defParm,
							&forback,
							&xStr);
	if (err)
	{
		goto FREE_RETURN;
	}

	free_str = 1;

	/*	Now, lookup underlyings involved and their term structures */

	fx_und = lookup_und (underlying);
	
	if (!fx_und)
	{
		err = serror ("Couldn't find underlying named %s", underlying);
		goto FREE_RETURN;
	}
	
	today = get_today_from_underlying (fx_und);

	if (get_underlying_type (fx_und) != FOREX_UND)
	{
		err = serror ("Underlying %s is not of type FX", underlying);
		goto FREE_RETURN;
	}

	if (get_mdltype_from_fxund (fx_und) != FX_LGMSV)
	{
		err = serror ("Underlying %s is not of type FX Stoch Rates", underlying);
		goto FREE_RETURN;
	}

	fxlgmsv_get_fx_and_correl_ts(underlying, &fx_spot, &num_fx_vol, &fx_time, &fx_vol, &num_rho, &rho_time, &CorrMatrix);

	domname = get_domname_from_fxund (fx_und);
	dom_und = lookup_und (domname);
	if (!dom_und)
	{
		err = serror ("Couldn't find underlying named %s", domname);
		goto FREE_RETURN;
	}
	// Get the model 
	err = Get_LGMSV_model(	domname,
							dom_model);

	forname = get_forname_from_fxund (fx_und);
	for_und = lookup_und (forname);
	if (!for_und)
	{
		err = serror ("Couldn't find underlying named %s", forname);
		goto FREE_RETURN;
	}
	// Get the model 
	err = Get_LGMSV_model(	forname,
							for_model);

	// look for the today date 
	today = get_today_from_underlying (dom_und);
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	dom_yc = (char *) get_ycname_from_irund (dom_und);
	for_yc = (char *) get_ycname_from_irund (for_und);
	
	// Get number of columns 
	err = FIRSTGetNumColFromDeal(
									&xStr,
									&num_col);	
	if (err)
	{
		goto FREE_RETURN;
	}

	//	Get the maximum number of dfs required	
	err = FIRSTGetMaxNumDfFromDeal(
									&xStr,
									&max_num_df);
	if (err)
	{
		goto FREE_RETURN;
	}
	
	//	Next, get the time steps 
	err = FIRSTGetEvtDatesFromDeal(
									&xStr,
 									&num_evt,
 									&evt_dts,
									&evt_tms);
	if (err)
	{
		goto FREE_RETURN;
	}
	

	if (err) goto FREE_RETURN;
																					
	// discretise in time			

	nstp = num_evt;

	time = (double*) calloc (nstp, sizeof (double));
	if (!time)
	{
		err = "Memory allocation error (1) in SrtGrfnFXLGMSVMC";
		goto FREE_RETURN;
	}
	memcpy (time, xStr.tms, nstp * sizeof (double));

	final_mat = xStr.tms[nstp-1];
	i = 0;
	while (i < dom_model->iNbPWTime && dom_model->dPWTime[i] < final_mat)
	{
		num_f_add_number(&nstp, &time, dom_model->dPWTime[i]);
		i++;
	}
	i = 0;
	while (i < for_model->iNbPWTime && for_model->dPWTime[i] < final_mat)
	{
		num_f_add_number(&nstp, &time, for_model->dPWTime[i]);
		i++;
	}
	i = 0;
	while (i < num_fx_vol && fx_time[i] < final_mat)
	{
		num_f_add_number(&nstp, &time, fx_time[i]);
		i++;
	}
	i = 0;
	while (i < num_rho && rho_time[i] < final_mat)
	{
		num_f_add_number(&nstp, &time, rho_time[i]);
		i++;
	}
	
	num_f_sort_vector (nstp, time);

	//	Fill the time vector 
	err = fill_time_vector (	&time, 
								&nstp,
								0,
								NULL,
								0,
								NULL, 
								nstept);				
	if (err)
	{
		goto FREE_RETURN;
	}

	/*	Fill product structure */

	strupper (underlying);
	strip_white_space (underlying);
	strupper (domname);
	strip_white_space (domname);
	strupper (forname);
	strip_white_space (forname);
	for (i=0; i<xStr.num_und; i++)
	{
		strupper (xStr.und_data[i].und_name);
		strip_white_space (xStr.und_data[i].und_name);
	}

	fx_idx = -1;
	for (i=0; i<xStr.num_und; i++)
	{
		if (!strcmp (xStr.und_data[i].und_name, underlying))
		{
			fx_idx = i;
		}
	}
	if (fx_idx == -1)
	{
		err = "The Fx underlying is not present in the mdlcomm structure";
		goto FREE_RETURN;
	}

	dom_idx = -1;
	for (i=0; i<xStr.num_und; i++)
	{
		if (!strcmp (xStr.und_data[i].und_name, domname))
		{
			dom_idx = i;
		}
	}
	if (dom_idx == -1)
	{
		err = "The domestic underlying is not present in the mdlcomm structure";
		goto FREE_RETURN;
	}

	for_idx = -1;
	for (i=0; i<xStr.num_und; i++)
	{
		if (!strcmp (xStr.und_data[i].und_name, forname))
		{
			for_idx = i;
		}
	}
	if (for_idx == -1)
	{
		err = "The foreign underlying is not present in the mdlcomm structure";
		goto FREE_RETURN;
	}

	void_prm = (void**) calloc (nstp, sizeof (void*));
	is_event = (int *) calloc (nstp, sizeof (int));
	date = (double *) calloc (nstp, sizeof (double));

	dom_sigma = (double *) calloc (nstp, sizeof (double));	
	dom_alpha = (double *) calloc (nstp, sizeof (double));
	dom_rho = (double *) calloc (nstp, sizeof (double));	
	dom_rho2 = (double *) calloc (nstp, sizeof (double));
	dom_lameps = (double *) calloc (nstp, sizeof (double));
	dom_lvleps = (double *) calloc (nstp, sizeof (double));
	dom_newalphaLGM = (double *) calloc (nstp, sizeof (double));
	dom_newrhoLGM = (double *) calloc (nstp, sizeof (double));

	dom_zcvol1_star = (double *) calloc(nstp, sizeof(double));
	dom_zcvol2_star = (double *) calloc(nstp, sizeof(double));
	
	dom_dff_star = (double *) calloc(num_evt, sizeof(double));
	dom_gam1_star = (double *) calloc(num_evt, sizeof(double));
	dom_gam2_star = (double *) calloc(num_evt, sizeof(double));
	dom_gam1_2_star = (double *) calloc(num_evt, sizeof(double));
	dom_gam2_2_star = (double *) calloc(num_evt, sizeof(double));
	dom_gam12_star = (double *) calloc(num_evt, sizeof(double));
	
	for_sigma = (double *) calloc (nstp, sizeof (double));	
	for_alpha = (double *) calloc (nstp, sizeof (double));
	for_rho = (double *) calloc (nstp, sizeof (double));	
	for_rho2 = (double *) calloc (nstp, sizeof (double));
	for_lameps = (double *) calloc (nstp, sizeof (double));
	for_lvleps = (double *) calloc (nstp, sizeof (double));
	for_newalphaLGM = (double *) calloc (nstp, sizeof (double));
	for_newrhoLGM = (double *) calloc (nstp, sizeof (double));
	
	for_zcvol1_star = (double *) calloc(nstp, sizeof(double));
	for_zcvol2_star = (double *) calloc(nstp, sizeof(double));

	for_dff_star = (double *) calloc(num_evt, sizeof(double));
	for_gam1_star = (double *) calloc(num_evt, sizeof(double));
	for_gam2_star = (double *) calloc(num_evt, sizeof(double));
	for_gam1_2_star = (double *) calloc(num_evt, sizeof(double));
	for_gam2_2_star = (double *) calloc(num_evt, sizeof(double));
	for_gam12_star = (double *) calloc(num_evt, sizeof(double));

	fx_sigma = (double *) calloc(nstp, sizeof(double));
	correlation = f3tensor(0, nstp-1, 0, 6, 0, 6);

	if (!void_prm || !is_event || !date ||
		!dom_sigma || !dom_alpha || !dom_rho || !dom_rho2 || !dom_lameps || !dom_lvleps || !dom_newalphaLGM || !dom_newrhoLGM ||
		!fx_sigma || !correlation ||
		!dom_zcvol1_star || !dom_zcvol2_star || !for_zcvol1_star || !for_zcvol2_star ||
		!dom_dff_star || !dom_gam1_star || !dom_gam2_star || !dom_gam1_2_star || !dom_gam2_2_star || !dom_gam12_star ||
		!for_sigma || !for_alpha || !for_rho || !for_rho2 || !for_lameps || !for_lvleps || !for_newalphaLGM || !for_newrhoLGM ||
		!for_dff_star || !for_gam1_star || !for_gam2_star || !for_gam1_2_star || !for_gam2_2_star || !for_gam12_star)
	{
		err = "Memory allocation failure in SrtGrfnFXLGMSVMC";
		goto FREE_RETURN;
	}

	j = xStr.num_evt - 1;	

	for (i=nstp-1; i>=0; i--)
	{
		date[i] = ((long) (today + time[i] * 365.0 + 1.0E-08)) * 1.0;	
		

		// Calculate Tstart stuff which are relevant at step i
		if (i < nstp-1)
		{
			index = Get_Index(time[i+1], dom_model->dPWTime, dom_model->iNbPWTime);

			dom_sigma[i+1] = dom_model->dSigma[index];
			dom_newalphaLGM[i+1] = dom_model->dLGMAlpha[index];
			dom_newrhoLGM[i+1] = dom_model->dLGMRho[index];
			
			dom_alpha[i+1] = dom_model->dAlpha[index];
			dom_rho[i+1] = dom_model->dRho[index];
			dom_rho2[i+1] = dom_model->dRho2[index];
			dom_lameps[i+1] = dom_model->dLambdaEps[index];
			dom_lvleps[i+1] = dom_model->dLvlEps[index];
			
			// J.M.L 10 Mar 2004
			// dom_zcvol1_star represents the term (1 - exp(lbda1 * ( Tstar - t )) / lbda1 which appears in
			// the volatility of the bond B(t,Tstar). We can either code it directly:
			
			// dom_zcvol1_star[i+1] = (1.0 - exp(dom_model->dLambdaX * (dom_model->dTStar - time[i]))) / dom_model->dLambdaX;
			
			// but a better plan is to integrate the square of this expression which has the dimension of a volatility
			// and find an equivalent volatility for the bond which is valid between ti and ti+1
			
			dom_zcvol1_star[i+1] = -sqrt( 
				( time[i+1]-time[i] 
				+ 2*(exp(dom_model->dLambdaX * (dom_model->dTStar - time[i+1])) 
						- exp(dom_model->dLambdaX * (dom_model->dTStar - time[i])))/dom_model->dLambdaX
				+ 0.5*(exp(2*dom_model->dLambdaX*(dom_model->dTStar-time[i])) 
						- exp(2*dom_model->dLambdaX*(dom_model->dTStar-time[i+1])))/dom_model->dLambdaX
				) / (time[i+1] - time[i]) ) / dom_model->dLambdaX;

			// dom_zcvol2_star[i+1] = (1.0 - exp(dom_model->dLambdaX2 * (dom_model->dTStar - time[i]))) / dom_model->dLambdaX2;
			
			dom_zcvol2_star[i+1] = -sqrt( 
				( time[i+1]-time[i] 
				+ 2*(exp(dom_model->dLambdaX2*(dom_model->dTStar - time[i+1])) 
						- exp(dom_model->dLambdaX2*(dom_model->dTStar-time[i])))/dom_model->dLambdaX2
				+ 0.5*(exp(2*dom_model->dLambdaX2*(dom_model->dTStar-time[i])) 
						- exp(2*dom_model->dLambdaX2*(dom_model->dTStar-time[i+1])))/dom_model->dLambdaX2
				) / (time[i+1] - time[i]) ) / dom_model->dLambdaX2;

			//	for_zcvol1_star[i+1] = (1.0 - exp(for_model->dLambdaX * (for_model->dTStar - time[i]))) / for_model->dLambdaX;
			//	for_zcvol2_star[i+1] = (1.0 - exp(for_model->dLambdaX2 * (for_model->dTStar - time[i]))) / for_model->dLambdaX2;


			for_zcvol1_star[i+1] = -sqrt( 
				( time[i+1]-time[i] 
				+ 2*(exp(for_model->dLambdaX * (for_model->dTStar - time[i+1])) 
						- exp(for_model->dLambdaX * (for_model->dTStar - time[i])))/for_model->dLambdaX
				+ 0.5*(exp(2*for_model->dLambdaX*(for_model->dTStar-time[i])) 
						- exp(2*for_model->dLambdaX*(for_model->dTStar-time[i+1])))/for_model->dLambdaX
				) / (time[i+1] - time[i]) ) / for_model->dLambdaX;

			for_zcvol2_star[i+1] = -sqrt( 
				( time[i+1]-time[i] 
				+ 2*(exp(for_model->dLambdaX2*(for_model->dTStar - time[i+1])) 
						- exp(for_model->dLambdaX2*(for_model->dTStar-time[i])))/for_model->dLambdaX2
				+ 0.5*(exp(2*for_model->dLambdaX2*(for_model->dTStar-time[i])) 
						- exp(2*for_model->dLambdaX2*(for_model->dTStar-time[i+1])))/for_model->dLambdaX2
				) / (time[i+1] - time[i]) ) / for_model->dLambdaX2;

			index = Get_Index(time[i+1], for_model->dPWTime, for_model->iNbPWTime);

			for_sigma[i+1] = for_model->dSigma[index];
			for_newalphaLGM[i+1] = for_model->dLGMAlpha[index];
			for_newrhoLGM[i+1] = for_model->dLGMRho[index];

			for_alpha[i+1] = for_model->dAlpha[index];
			for_rho[i+1] = for_model->dRho[index];
			for_rho2[i+1] = for_model->dRho2[index];
			for_lameps[i+1] = for_model->dLambdaEps[index];
			for_lvleps[i+1] = for_model->dLvlEps[index];

			index = Get_Index(time[i+1], fx_time, num_fx_vol);
			fx_sigma[i+1] = fx_vol[index];

			index = Get_Index(time[i+1], rho_time, num_rho);
			for(l=0;l<7;++l)
			{
				for(k=0;k<7;++k)
				{
					correlation[i+1][l][k] = CorrMatrix[index][l][k];
				}
			}
		}

		// In case step i corresponds to an event date
		if (j >=0 && fabs(time[i] - evt_tms[j]) < 1.0E-08)
		{
			grfn_prm = malloc(sizeof(grfn_parm_mc_fxlgmsv));

			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j;
			grfn_prm->fx_idx = fx_idx;
			grfn_prm->dom_idx = dom_idx;
			grfn_prm->for_idx = for_idx;

			// DOMESTIC RECONSTRUCTION PRE-CALCULATIONS
			grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
			grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
			grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];

			// DF(t, T*) reconstruction 
			dom_dff_star[j] = swp_f_df (xStr.dts[j], dom_model->lTStarDate, dom_yc);
			dom_dff_star[j] = log(dom_dff_star[j]);
			dom_gam1_star[j] = (1.0 - exp(dom_model->dLambdaX * (dom_model->dTStar - evt_tms[j]))) / dom_model->dLambdaX;
			dom_gam2_star[j] = (1.0 - exp(dom_model->dLambdaX2 * (dom_model->dTStar - evt_tms[j]))) / dom_model->dLambdaX2;
			dom_gam12_star[j] = dom_gam1_star[j] * dom_gam2_star[j];
			dom_gam1_2_star[j] = 0.5 * dom_gam1_star[j] * dom_gam1_star[j];
			dom_gam2_2_star[j] = 0.5 * dom_gam2_star[j] * dom_gam2_star[j];

			if (grfn_prm->num_dom_df > 0 && grfn_prm->dom_idx != -1)
			{
				grfn_prm->dom_dff = dvector (0, grfn_prm->num_dom_df-1);
				grfn_prm->dom_gam1 = dvector (0, grfn_prm->num_dom_df-1);
				grfn_prm->dom_gam1_2 = dvector (0, grfn_prm->num_dom_df-1);
				grfn_prm->dom_gam2 = dvector (0, grfn_prm->num_dom_df-1);
				grfn_prm->dom_gam2_2 = dvector (0, grfn_prm->num_dom_df-1);
				grfn_prm->dom_gam12 = dvector (0, grfn_prm->num_dom_df-1);				

				if (!grfn_prm->dom_dff || !grfn_prm->dom_gam1 || !grfn_prm->dom_gam1_2 || !grfn_prm->dom_gam2 || !grfn_prm->dom_gam2_2 || !grfn_prm->dom_gam12)
				{
					err = "Memory allocation error in SrtGrfnFXLGMSVMC";
					goto FREE_RETURN;
				}

				for (k=0; k<grfn_prm->num_dom_df; k++)
				{
					grfn_prm->dom_dff[k] = swp_f_df (xStr.dts[j], grfn_prm->dom_df_dts[k], dom_yc);
					grfn_prm->dom_dff[k] = log(grfn_prm->dom_dff[k]);
					temp1 = (1.0 - exp(dom_model->dLambdaX * (dom_model->dTStar - evt_tms[j] - grfn_prm->dom_df_tms[k]))) / dom_model->dLambdaX;
					temp2 = (1.0 - exp(dom_model->dLambdaX2 * (dom_model->dTStar - evt_tms[j] - grfn_prm->dom_df_tms[k]))) / dom_model->dLambdaX2;
					grfn_prm->dom_gam1[k] = temp1 - dom_gam1_star[j];
					grfn_prm->dom_gam2[k] = temp2 - dom_gam2_star[j];					
					grfn_prm->dom_gam1_2[k] = 0.5 * temp1 * temp1 - dom_gam1_2_star[j];
					grfn_prm->dom_gam2_2[k] = 0.5 * temp2 * temp2 - dom_gam2_2_star[j];
					grfn_prm->dom_gam12[k] = temp1 * temp2 - dom_gam12_star[j];
				}

				grfn_prm->do_dom = 1;
			}
			else
			{
				grfn_prm->do_dom = 0;
			}
						

			// FX = DOMESTIC RECONSTRUCTION PRE-CALCULATIONS
			grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
			grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
			grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];

			if (grfn_prm->num_fx_df > 0 && grfn_prm->fx_idx != -1)
			{
				grfn_prm->fx_dff = dvector (0, grfn_prm->num_fx_df-1);
				grfn_prm->fx_gam1 = dvector (0, grfn_prm->num_fx_df-1);
				grfn_prm->fx_gam1_2 = dvector (0, grfn_prm->num_fx_df-1);
				grfn_prm->fx_gam2 = dvector (0, grfn_prm->num_fx_df-1);
				grfn_prm->fx_gam2_2 = dvector (0, grfn_prm->num_fx_df-1);
				grfn_prm->fx_gam12 = dvector (0, grfn_prm->num_fx_df-1);				

				if (!grfn_prm->fx_dff || !grfn_prm->fx_gam1 || !grfn_prm->fx_gam1_2 || !grfn_prm->fx_gam2 ||
					!grfn_prm->fx_gam2_2 || !grfn_prm->fx_gam12)
				{
					err = "Memory allocation error in SrtGrfnMainFXLGMSV";
					goto FREE_RETURN;
				}

				for (k=0; k<grfn_prm->num_fx_df; k++)
				{
					grfn_prm->fx_dff[k] = swp_f_df (xStr.dts[j], grfn_prm->fx_df_dts[k], dom_yc);
					grfn_prm->fx_dff[k] = log(grfn_prm->fx_dff[k]);
					temp1 = (1.0 - exp(dom_model->dLambdaX * (dom_model->dTStar - evt_tms[j] - grfn_prm->fx_df_tms[k]))) / dom_model->dLambdaX;
					temp2 = (1.0 - exp(dom_model->dLambdaX2 * (dom_model->dTStar - evt_tms[j] - grfn_prm->fx_df_tms[k]))) / dom_model->dLambdaX2;
					grfn_prm->fx_gam1[k] = temp1 - dom_gam1_star[j];
					grfn_prm->fx_gam2[k] = temp2 - dom_gam2_star[j];					
					grfn_prm->fx_gam1_2[k] = 0.5 * temp1 * temp1 - dom_gam1_2_star[j];
					grfn_prm->fx_gam2_2[k] = 0.5 * temp2 * temp2 - dom_gam2_2_star[j];
					grfn_prm->fx_gam12[k] = temp1 * temp2 - dom_gam12_star[j];
				}

				grfn_prm->do_fx = 1;
			}
			else
			{
				grfn_prm->do_fx = 0;
			}

			// FOREIGN RECONSTRUCTION PRE-CALCULATIONS
			grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
			grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
			grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];

			// DF(t, T*) reconstruction 
			for_dff_star[j] = swp_f_df (xStr.dts[j], for_model->lTStarDate, for_yc);
			for_dff_star[j] = log(for_dff_star[j]);
			for_gam1_star[j] = (1.0 - exp(for_model->dLambdaX * (for_model->dTStar - evt_tms[j]))) / for_model->dLambdaX;
			for_gam2_star[j] = (1.0 - exp(for_model->dLambdaX2 * (for_model->dTStar - evt_tms[j]))) / for_model->dLambdaX2;
			for_gam12_star[j] = for_gam1_star[j] * for_gam2_star[j];
			for_gam1_2_star[j] = 0.5 * for_gam1_star[j] * for_gam1_star[j];
			for_gam2_2_star[j] = 0.5 * for_gam2_star[j] * for_gam2_star[j];

			if (grfn_prm->num_for_df > 0 && grfn_prm->for_idx != -1)
			{
				grfn_prm->for_dff = dvector (0, grfn_prm->num_for_df-1);
				grfn_prm->for_gam1 = dvector (0, grfn_prm->num_for_df-1);
				grfn_prm->for_gam1_2 = dvector (0, grfn_prm->num_for_df-1);
				grfn_prm->for_gam2 = dvector (0, grfn_prm->num_for_df-1);
				grfn_prm->for_gam2_2 = dvector (0, grfn_prm->num_for_df-1);
				grfn_prm->for_gam12 = dvector (0, grfn_prm->num_for_df-1);				

				if (!grfn_prm->for_dff || !grfn_prm->for_gam1 || !grfn_prm->for_gam1_2 || !grfn_prm->for_gam2 || !grfn_prm->for_gam2_2 || !grfn_prm->for_gam12)
				{
					err = "Memory allocation error in SrtGrfnMainFXLGMSV";
					goto FREE_RETURN;
				}

				for (k=0; k<grfn_prm->num_for_df; k++)
				{
					grfn_prm->for_dff[k] = swp_f_df (xStr.dts[j], grfn_prm->for_df_dts[k], for_yc);
					grfn_prm->for_dff[k] = log(grfn_prm->for_dff[k]);
					temp1 = (1.0 - exp(for_model->dLambdaX * (for_model->dTStar - evt_tms[j] - grfn_prm->for_df_tms[k]))) / for_model->dLambdaX;
					temp2 = (1.0 - exp(for_model->dLambdaX2 * (for_model->dTStar - evt_tms[j] - grfn_prm->for_df_tms[k]))) / for_model->dLambdaX2;
					grfn_prm->for_gam1[k] = temp1 - for_gam1_star[j];
					grfn_prm->for_gam2[k] = temp2 - for_gam2_star[j];					
					grfn_prm->for_gam1_2[k] = 0.5 * temp1 * temp1 - for_gam1_2_star[j];
					grfn_prm->for_gam2_2[k] = 0.5 * temp2 * temp2 - for_gam2_2_star[j];
					grfn_prm->for_gam12[k] = temp1 * temp2 - for_gam12_star[j];
				}

				grfn_prm->do_for = 1;
			}
			else
			{
				grfn_prm->do_for = 0;
			}

			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;

			j--;
			while (j >= 0 && xStr.evt[j].evt == NULL)
			{
				j--;
			}
		}
		else if (j>=0 && xStr.am[j])
		{
			grfn_prm = malloc(sizeof(grfn_parm_mc_fxlgmsv));
			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j;
			grfn_prm->fx_idx = fx_idx;
			grfn_prm->dom_idx = dom_idx;
			grfn_prm->for_idx = for_idx;

			grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
			grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
			grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];

			grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
			grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
			grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];

			grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
			grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
			grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];

			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;
		}
		else
		{
			is_event[i] = 0;
			void_prm[i] = NULL;
		}
	}

	//	Eventually! call to function 
	
	*nb_prod = num_col;	
	*tableauCols = num_col;
	*resRows = max(num_col + 1, num_evt);

	if (do_optimisation)
	{		
		if (params->iMultiIndex)
		{
			params->iNbIndex = params->iMultiIndex;
		}
		else
		{
			params->iNbIndex = 1;
		}

		err = mceb_allocate_params(	params,
									num_evt);

		if (err) goto FREE_RETURN;

		*prod_val = dmatrix (0, *resRows - 1, 0, 2 + params->iNbIndex);
	}
	else
	{
		*prod_val = dmatrix (0, *nb_prod - 1, 0, 2);
	}
			
	err = fxlgmsv_mc_balsam(nstp,
							num_evt,
							time,
							date,
							
							numpaths,
							
							dom_model->dLambdaX,
							dom_model->dLambdaX2,
							
							dom_sigma,
							dom_newalphaLGM,
							dom_newrhoLGM,
							dom_alpha,
							dom_lameps,
							dom_lvleps,
							dom_rho,
							dom_rho2,

							dom_zcvol1_star,
							dom_zcvol2_star,
							
							dom_dff_star,
							dom_gam1_star,
							dom_gam2_star,
							dom_gam1_2_star,
							dom_gam2_2_star,
							dom_gam12_star,

							for_model->dLambdaX,
							for_model->dLambdaX2,
							
							for_sigma,
							for_newalphaLGM,
							for_newrhoLGM,
							for_alpha,
							for_lameps,
							for_lvleps,
							for_rho,
							for_rho2,
							
							for_zcvol1_star,
							for_zcvol2_star,

							for_dff_star,
							for_gam1_star,
							for_gam2_star,
							for_gam1_2_star,
							for_gam2_2_star,
							for_gam12_star,
							
							fx_spot * swp_f_df (today, for_model->lTStarDate, for_yc) / swp_f_df (today, dom_model->lTStarDate, dom_yc),
							fx_sigma,
							
							correlation,
							
							void_prm, 
							is_event,

							do_optimisation,
							optimise,
							params,
							
							NULL,
							
							payoff_fxlgmsv_mc,
							
							num_col,
							*prod_val);

	if (err) goto FREE_RETURN;

	if (do_optimisation)
	{
		/* Recopy Barrier / CoefLin for the moment */
		for (i=0; i<num_evt; i++)
		{
			(*prod_val)[i][2] = params->dBarrier[i];

			for (j=0; j<params->iNbIndex; j++)
			{
				(*prod_val)[i][3+j] = params->dCoefLin[i][j+1];
			}
		}
	}

	dom_df = swp_f_df (today, dom_model->lTStarDate, dom_yc);

	for (i=0; i<num_col; i++)
	{
		(*prod_val)[i][0] *= dom_df;
		(*prod_val)[i][1] *= dom_df;
	}

	if (do_optimisation)
	{
		(*prod_val)[num_col][0] *= dom_df;
		(*prod_val)[num_col][1] *= dom_df;
	}

	if (flag)
	{		
		for (i=0; i<num_evt-1; i++)
		{			
			(*prod_val)[i][2] = (*prod_val)[i+1][2];			
			
			if (do_optimisation)
			{
				for (k=0; k<params->iNbIndex; k++)
				{
					(*prod_val)[i][3+k] = (*prod_val)[i+1][3+k];
				}
			}
		}
	}
	
	//	Add PV of Past 
	(*prod_val)[num_col-1][0] += xStr.gd->pv_of_past;

		
FREE_RETURN:
	
	if (free_str)
	{
		FIRSTFreeUndFromDeal(
									num_und,
									&und_ptr
								   );

		FIRSTFreeEvtDatesFromDeal(
										nstp,
										&evt_dts,
										&evt_tms
										);

		FIRSTFreeMktStruct(&xStr);
	}

	if (void_prm)
	{
		for (i=0; i<nstp; i++)
		{
			if (void_prm[i])
			{
				grfn_prm = (GRFNPARMMCFXLGMSV) void_prm[i];
				free (grfn_prm);
			}
		}

		free (void_prm);
	}	

	if (time) free (time);
	if (is_event) free (is_event);
	if (date) free (date);	

	if (fx_sigma) free (fx_sigma);
	if (correlation) free_f3tensor(correlation, 0, nstp-1, 0, 6, 0, 6);

	if (dom_zcvol1_star) free (dom_zcvol1_star);
	if (dom_zcvol2_star) free (dom_zcvol2_star);
	if (for_zcvol1_star) free (for_zcvol1_star);
	if (for_zcvol2_star) free (for_zcvol2_star);

	if (dom_sigma) free (dom_sigma);
	if (dom_alpha) free (dom_alpha);
	if (dom_rho) free (dom_rho);
	if (dom_rho2) free (dom_rho2);
	if (dom_lameps) free (dom_lameps);
	if (dom_lvleps) free (dom_lvleps);
	if (dom_newalphaLGM) free (dom_newalphaLGM);
	if (dom_newrhoLGM) free (dom_newrhoLGM);

	if (dom_dff_star) free(dom_dff_star);
	if (dom_gam1_star) free(dom_gam1_star);
	if (dom_gam2_star) free(dom_gam2_star);
	if (dom_gam1_2_star) free(dom_gam1_2_star);
	if (dom_gam2_2_star) free(dom_gam2_2_star);
	if (dom_gam12_star) free(dom_gam12_star);

	if (for_sigma) free (for_sigma);
	if (for_alpha) free (for_alpha);
	if (for_rho) free (for_rho);
	if (for_rho2) free (for_rho2);
	if (for_lameps) free (for_lameps);
	if (for_lvleps) free (for_lvleps);
	if (for_newalphaLGM) free (for_newalphaLGM);
	if (for_newrhoLGM) free (for_newrhoLGM);

	if (for_dff_star) free(for_dff_star);
	if (for_gam1_star) free(for_gam1_star);
	if (for_gam2_star) free(for_gam2_star);
	if (for_gam1_2_star) free(for_gam1_2_star);
	if (for_gam2_2_star) free(for_gam2_2_star);
	if (for_gam12_star) free(for_gam12_star);

	if(fx_time)
	{
		free(fx_time);
		fx_time = NULL;
	}
	
	if(fx_vol)
	{
		free(fx_vol);
		fx_vol = NULL;
	}

	if(rho_time)
	{
		free(rho_time);
		rho_time = NULL;
	}

	if(CorrMatrix)
	{
		free_f3tensor(CorrMatrix, 0, num_rho-1, 0, 6, 0, 6);
		CorrMatrix = NULL;
	}

	free_LGMSV_model(dom_model);
	free_LGMSV_model(for_model);
	
	return err;
}


char *SrtGrfnQTOLGMSVMC(
					char		*underlying,
					int			numeventdates, 
					long		*eventdates,
					long		tableauRows,
					long		*tableauCols,
					char		***tableauStrings,
					int			**tableauMask,
					long		auxWidth,
					long		*auxLen,
					double		**aux,

					/* for Optimisation of exercise boundary */
					int			do_optimisation,
					int			*optimise,
					MCEBPARAMS	params,
					long		*resRows,

					int			nstept,
					int			numpaths,				  
					int			*nb_prod,
					double		***prod_val)
{
	Err			err = NULL;
	SrtUndPtr	pUndPtr;
	int			for_one2F;

	// Get the underlying through its name and check it exists
	// Check on the underlying type 
	pUndPtr  = lookup_und(underlying);
	if (!pUndPtr ) 
	{
		return serror("Underlying %s is not defined", underlying);				
	}

	if (get_underlying_type (pUndPtr) != FOREX_UND)
	{
		return serror("Underlying %s is not of type FX", underlying);
	}

	if (get_mdltype_from_irund (pUndPtr) != FX_LGMSV)
	{
		return serror("Underlying %s is not of type FXLGMSV", underlying);
	}

	err = qtolgmsv_check_und(underlying, &for_one2F);
	if(err)
	{
		return err;
	}

	if (for_one2F == 1)
	{
		err = SrtGrfnQTOLGMSV1FMC(underlying,
								numeventdates, 
								eventdates,
								tableauRows,
								tableauCols,
								tableauStrings,
								tableauMask,
								auxWidth,
								auxLen,
								aux,
								do_optimisation,
								optimise,
								params,
								resRows,
								nstept,
								numpaths,		
								nb_prod,
								prod_val);
	}
	else
	{
		err = SrtGrfnQTOLGMSV2FMC(underlying,
								numeventdates, 
								eventdates,
								tableauRows,
								tableauCols,
								tableauStrings,
								tableauMask,
								auxWidth,
								auxLen,
								aux,
								do_optimisation,
								optimise,
								params,
								resRows,
								nstept,
								numpaths,		
								nb_prod,
								prod_val);
	}

	return err;
}


char *SrtGrfnQTOLGMSV2FMC(
					char		*underlying,
					int			numeventdates,
					long		*eventdates,
					long		tableauRows,
					long		*tableauCols,
					char		***tableauStrings,
					int			**tableauMask,
					long		auxWidth,
					long		*auxLen,
					double		**aux,

					// for Optimisation of exercise boundary 
					int			do_optimisation,
					int			*optimise,
					MCEBPARAMS	params,
					long		*resRows,

					int			nstept,
					int			numpaths,				  
					int			*nb_prod,
					double		***prod_val)
{
	int					free_str = 0;
	FIRSTAllMkts		xStr;
	SrtGrfnParam		defParm;
	GRFNPARMMCFXLGMSV	grfn_prm;
	int					forback;
	int					flag = 0;
	long				nstp;	
						
	double				*evt_tms		= NULL,						
						*time			= NULL,
						*date			= NULL,

						*dom_sigma			= NULL,
						*dom_zcvol_star	= NULL,

						*dom_dff_star		= NULL,
						*dom_gam1_star		= NULL,
						*dom_gam1_2_star	= NULL,

						*for_sigma			= NULL,
						*for_alpha			= NULL,
						*for_rho			= NULL,
						*for_rho2			= NULL,
						*for_lameps			= NULL,
						*for_lvleps			= NULL,
						*for_newalphaLGM	= NULL,
						*for_newrhoLGM		= NULL,
						
						*for_zcvol1_star	= NULL,
						*for_zcvol2_star	= NULL,

						*for_dff_star		= NULL,
						*for_gam1_star		= NULL,
						*for_gam2_star		= NULL,
						*for_gam1_2_star	= NULL,
						*for_gam2_2_star	= NULL,
						*for_gam12_star		= NULL;

	long				*evt_dts	= NULL;
	
	int					*is_event	= NULL;
	void				**void_prm	= NULL;
		
	long				index;
	double				final_mat;

	long				today, spot_date;
						
	int					num_col		= 0,
						max_num_df	= 0,
						num_evt		= 0,
						num_und		= 0;
						
	char				*domname,
						*forname;

	char				*dom_yc, *for_yc;

	SrtUndPtr			fx_und, dom_und, for_und;
	SrtUndPtr			*und_ptr	= NULL,
						und			= NULL;

	LGMSV_model			dom_model_, *dom_model = &dom_model_;
	LGMSV_model			for_model_, *for_model = &for_model_;

	int					i, j, k, l;

	int					fx_idx, dom_idx, for_idx;

	double				dom_df;
	double				temp1, temp2;
	
	Err					err = NULL;

	double				fx_spot;

	int					num_fx_vol;
	double				*fx_time=NULL;
	double				*fx_vol=NULL;

	int					num_rho;
	double				*rho_time=NULL;
	double				***CorrMatrix=NULL;

	double				*fx_sigma=NULL;
	double				***correlation=NULL;

	//	Initialisation of the LGMSV model 
	init_NULL_LGMSV_model(dom_model);
	init_NULL_LGMSV_model(for_model);

	//	Initialise the GRFN tableau 

	//	First, initialise the param struct 	
	err = srt_f_set_default_GrfnParams (&defParm);	
	defParm.force_mc = 0;

	err = FIRSTInitMktStruct(
							numeventdates,
							eventdates,
							tableauRows,
							*tableauCols,
							tableauStrings,
							tableauMask,
							auxWidth,
							auxLen,
							aux,
							underlying,
							&defParm,
							&forback,
							&xStr);
	if (err)
	{
		goto FREE_RETURN;
	}

	free_str = 1;

	/*	Now, lookup underlyings involved and their term structures */

	fx_und = lookup_und (underlying);
	
	if (!fx_und)
	{
		err = serror ("Couldn't find underlying named %s", underlying);
		goto FREE_RETURN;
	}
	
	today = get_today_from_underlying (fx_und);

	if (get_underlying_type (fx_und) != FOREX_UND)
	{
		err = serror ("Underlying %s is not of type FX", underlying);
		goto FREE_RETURN;
	}

	if (get_mdltype_from_fxund (fx_und) != FX_LGMSV)
	{
		err = serror ("Underlying %s is not of type FX Stoch Rates", underlying);
		goto FREE_RETURN;
	}

	fxlgmsv_get_fx_and_correl_ts(underlying, &fx_spot, &num_fx_vol, &fx_time, &fx_vol, &num_rho, &rho_time, &CorrMatrix);

	domname = get_domname_from_fxund (fx_und);
	dom_und = lookup_und (domname);
	if (!dom_und)
	{
		err = serror ("Couldn't find underlying named %s", domname);
		goto FREE_RETURN;
	}
	// Get the model 
	err = Get_LGMSV_model(	domname,
							dom_model);

	forname = get_forname_from_fxund (fx_und);
	for_und = lookup_und (forname);
	if (!for_und)
	{
		err = serror ("Couldn't find underlying named %s", forname);
		goto FREE_RETURN;
	}
	// Get the model 
	err = Get_LGMSV_model(	forname,
							for_model);

	// look for the today date 
	today = get_today_from_underlying (dom_und);
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	dom_yc = (char *) get_ycname_from_irund (dom_und);
	for_yc = (char *) get_ycname_from_irund (for_und);
	
	// Get number of columns 
	err = FIRSTGetNumColFromDeal(
									&xStr,
									&num_col);	
	if (err)
	{
		goto FREE_RETURN;
	}

	//	Get the maximum number of dfs required	
	err = FIRSTGetMaxNumDfFromDeal(
									&xStr,
									&max_num_df);
	if (err)
	{
		goto FREE_RETURN;
	}
	
	//	Next, get the time steps 
	err = FIRSTGetEvtDatesFromDeal(
									&xStr,
 									&num_evt,
 									&evt_dts,
									&evt_tms);
	if (err)
	{
		goto FREE_RETURN;
	}
	

	if (err) goto FREE_RETURN;
																					
	// discretise in time			

	nstp = num_evt;

	time = (double*) calloc (nstp, sizeof (double));
	if (!time)
	{
		err = "Memory allocation error (1) in SrtGrfnQTOLGMSV2FMC";
		goto FREE_RETURN;
	}
	memcpy (time, xStr.tms, nstp * sizeof (double));

	final_mat = xStr.tms[nstp-1];
	i = 0;
	while (i < dom_model->iNbPWTime && dom_model->dPWTime[i] < final_mat)
	{
		num_f_add_number(&nstp, &time, dom_model->dPWTime[i]);
		i++;
	}
	i = 0;
	while (i < for_model->iNbPWTime && for_model->dPWTime[i] < final_mat)
	{
		num_f_add_number(&nstp, &time, for_model->dPWTime[i]);
		i++;
	}
	i = 0;
	while (i < num_fx_vol && fx_time[i] < final_mat)
	{
		num_f_add_number(&nstp, &time, fx_time[i]);
		i++;
	}
	i = 0;
	while (i < num_rho && rho_time[i] < final_mat)
	{
		num_f_add_number(&nstp, &time, rho_time[i]);
		i++;
	}
	
	num_f_sort_vector (nstp, time);

	//	Fill the time vector 
	err = fill_time_vector (	&time, 
								&nstp,
								0,
								NULL,
								0,
								NULL, 
								nstept);				
	if (err)
	{
		goto FREE_RETURN;
	}

	/*	Fill product structure */

	strupper (underlying);
	strip_white_space (underlying);
	strupper (domname);
	strip_white_space (domname);
	strupper (forname);
	strip_white_space (forname);
	for (i=0; i<xStr.num_und; i++)
	{
		strupper (xStr.und_data[i].und_name);
		strip_white_space (xStr.und_data[i].und_name);
	}

	fx_idx = -1;
	for (i=0; i<xStr.num_und; i++)
	{
		if (!strcmp (xStr.und_data[i].und_name, underlying))
		{
			fx_idx = i;
		}
	}
	if (fx_idx == -1)
	{
		err = "The Fx underlying is not present in the mdlcomm structure";
		goto FREE_RETURN;
	}

	dom_idx = -1;
	for (i=0; i<xStr.num_und; i++)
	{
		if (!strcmp (xStr.und_data[i].und_name, domname))
		{
			dom_idx = i;
		}
	}
	if (dom_idx == -1)
	{
		err = "The domestic underlying is not present in the mdlcomm structure";
		goto FREE_RETURN;
	}

	for_idx = -1;
	for (i=0; i<xStr.num_und; i++)
	{
		if (!strcmp (xStr.und_data[i].und_name, forname))
		{
			for_idx = i;
		}
	}
	if (for_idx == -1)
	{
		err = "The foreign underlying is not present in the mdlcomm structure";
		goto FREE_RETURN;
	}

	void_prm = (void**) calloc (nstp, sizeof (void*));
	is_event = (int *) calloc (nstp, sizeof (int));
	date = (double *) calloc (nstp, sizeof (double));

	dom_sigma = (double *) calloc (nstp, sizeof (double));	
	dom_zcvol_star = (double *) calloc(nstp, sizeof(double));
	
	dom_dff_star = (double *) calloc(num_evt, sizeof(double));
	dom_gam1_star = (double *) calloc(num_evt, sizeof(double));
	dom_gam1_2_star = (double *) calloc(num_evt, sizeof(double));
	
	for_sigma = (double *) calloc (nstp, sizeof (double));	
	for_alpha = (double *) calloc (nstp, sizeof (double));
	for_rho = (double *) calloc (nstp, sizeof (double));	
	for_rho2 = (double *) calloc (nstp, sizeof (double));
	for_lameps = (double *) calloc (nstp, sizeof (double));
	for_lvleps = (double *) calloc (nstp, sizeof (double));
	for_newalphaLGM = (double *) calloc (nstp, sizeof (double));
	for_newrhoLGM = (double *) calloc (nstp, sizeof (double));
	
	for_zcvol1_star = (double *) calloc(nstp, sizeof(double));
	for_zcvol2_star = (double *) calloc(nstp, sizeof(double));

	for_dff_star = (double *) calloc(num_evt, sizeof(double));
	for_gam1_star = (double *) calloc(num_evt, sizeof(double));
	for_gam2_star = (double *) calloc(num_evt, sizeof(double));
	for_gam1_2_star = (double *) calloc(num_evt, sizeof(double));
	for_gam2_2_star = (double *) calloc(num_evt, sizeof(double));
	for_gam12_star = (double *) calloc(num_evt, sizeof(double));

	fx_sigma = (double *) calloc(nstp, sizeof(double));
	
	// Original Call
	correlation = f3tensor(0, nstp-1, 0, 4, 0, 4);

	// TEST FXLGMSV
	// correlation = f3tensor(0, nstp-1, 0, 6, 0, 6);

	if (!void_prm || !is_event || !date ||
		!dom_sigma ||
		!fx_sigma || !correlation ||
		!dom_zcvol_star || !for_zcvol1_star || !for_zcvol2_star ||
		!dom_dff_star || !dom_gam1_star || !dom_gam1_2_star ||
		!for_sigma || !for_alpha || !for_rho || !for_rho2 || !for_lameps || !for_lvleps || !for_newalphaLGM || !for_newrhoLGM ||
		!for_dff_star || !for_gam1_star || !for_gam2_star || !for_gam1_2_star || !for_gam2_2_star || !for_gam12_star)
	{
		err = "Memory allocation failure in SrtGrfnQTOLGMSV2FMC";
		goto FREE_RETURN;
	}

	j = xStr.num_evt - 1;	

	for (i=nstp-1; i>=0; i--)
	{
		date[i] = ((long) (today + time[i] * 365.0 + 1.0E-08)) * 1.0;	
		
		if (i < nstp-1)
		{
			index = Get_Index(time[i+1], dom_model->dPWTime, dom_model->iNbPWTime);

			dom_sigma[i+1] = dom_model->dSigma[index];

			// J.M.L 10 Mar 2004
			// dom_zcvol_star represents the term (1 - exp(lbda1 * ( Tstar - t )) / lbda1 which appears in
			// the volatility of the bond B(t,Tstar). We can either code it directly:

			// Previous code:
			// dom_zcvol_star[i+1] = (1.0 - exp(dom_model->dLambdaX * (dom_model->dTStar - time[i]))) / dom_model->dLambdaX;
			
			// but a better plan is to integrate the square of this expression which has the dimension of a volatility
			// and find an equivalent volatility for the bond which is valid between ti and ti+1
			
			dom_zcvol_star[i+1] = -sqrt( 
				( time[i+1]-time[i] 
				+ 2*(exp(dom_model->dLambdaX * (dom_model->dTStar - time[i+1])) 
						- exp(dom_model->dLambdaX * (dom_model->dTStar - time[i])))/dom_model->dLambdaX
				+ 0.5*(exp(2*dom_model->dLambdaX*(dom_model->dTStar-time[i])) 
						- exp(2*dom_model->dLambdaX*(dom_model->dTStar-time[i+1])))/dom_model->dLambdaX
				) / (time[i+1] - time[i]) ) / dom_model->dLambdaX;

			// Previous code:
			// for_zcvol1_star[i+1] = (1.0 - exp(for_model->dLambdaX * (for_model->dTStar - time[i]))) / for_model->dLambdaX;
			// for_zcvol2_star[i+1] = (1.0 - exp(for_model->dLambdaX2 * (for_model->dTStar - time[i]))) / for_model->dLambdaX2;

			for_zcvol1_star[i+1] = -sqrt( 
				( time[i+1]-time[i] 
				+ 2*(exp(for_model->dLambdaX * (for_model->dTStar - time[i+1])) 
						- exp(for_model->dLambdaX * (for_model->dTStar - time[i])))/for_model->dLambdaX
				+ 0.5*(exp(2*for_model->dLambdaX*(for_model->dTStar-time[i])) 
						- exp(2*for_model->dLambdaX*(for_model->dTStar-time[i+1])))/for_model->dLambdaX
				) / (time[i+1] - time[i]) ) / for_model->dLambdaX;

			for_zcvol2_star[i+1] = -sqrt( 
				( time[i+1]-time[i] 
				+ 2*(exp(for_model->dLambdaX2*(for_model->dTStar - time[i+1])) 
						- exp(for_model->dLambdaX2*(for_model->dTStar-time[i])))/for_model->dLambdaX2
				+ 0.5*(exp(2*for_model->dLambdaX2*(for_model->dTStar-time[i])) 
						- exp(2*for_model->dLambdaX2*(for_model->dTStar-time[i+1])))/for_model->dLambdaX2
				) / (time[i+1] - time[i]) ) / for_model->dLambdaX2;

			index = Get_Index(time[i+1], for_model->dPWTime, for_model->iNbPWTime);

			for_sigma[i+1] = for_model->dSigma[index];
			for_newalphaLGM[i+1] = for_model->dLGMAlpha[index];
			for_newrhoLGM[i+1] = for_model->dLGMRho[index];
			
			for_alpha[i+1] = for_model->dAlpha[index];
			for_rho[i+1] = for_model->dRho[index];
			for_rho2[i+1] = for_model->dRho2[index];
			for_lameps[i+1] = for_model->dLambdaEps[index];
			for_lvleps[i+1] = for_model->dLvlEps[index];

			index = Get_Index(time[i+1], fx_time, num_fx_vol);
			fx_sigma[i+1] = fx_vol[index];

			index = Get_Index(time[i+1], rho_time, num_rho);

			
			correlation[i+1][0][0] = CorrMatrix[index][0][0];
			correlation[i+1][0][1] = CorrMatrix[index][0][3];
			correlation[i+1][0][2] = CorrMatrix[index][0][4];
			correlation[i+1][0][3] = CorrMatrix[index][0][5];
			correlation[i+1][0][4] = CorrMatrix[index][0][6];

			correlation[i+1][1][1] = CorrMatrix[index][3][3];
			correlation[i+1][1][2] = CorrMatrix[index][3][4];
			correlation[i+1][1][3] = CorrMatrix[index][3][5];
			correlation[i+1][1][4] = CorrMatrix[index][3][6];

			correlation[i+1][2][2] = CorrMatrix[index][4][4];
			correlation[i+1][2][3] = CorrMatrix[index][4][5];
			correlation[i+1][2][4] = CorrMatrix[index][4][6];

			correlation[i+1][3][3] = CorrMatrix[index][5][5];
			correlation[i+1][3][4] = CorrMatrix[index][5][6];

			correlation[i+1][4][4] = CorrMatrix[index][6][6];
			for(l=0;l<5;++l)
			{
				for(k=0;k<l;++k)
				{
					correlation[i+1][l][k] = correlation[i+1][k][l];
				}
			}

			// For testing purposes if all the correlation matrix is passed.
			/*
			for(l=0;l<7;++l)
			{
				for(k=0;k<7;++k)
				{
					correlation[i+1][l][k] = CorrMatrix[index][l][k];
				}
			}

			*/
		}

		if (j >=0 && fabs(time[i] - evt_tms[j]) < 1.0E-08)
		{
			grfn_prm = malloc(sizeof(grfn_parm_mc_fxlgmsv));

			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j;
			grfn_prm->fx_idx = fx_idx;
			grfn_prm->dom_idx = dom_idx;
			grfn_prm->for_idx = for_idx;

			// DOMESTIC RECONSTRUCTION PRE-CALCULATIONS
			grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
			grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
			grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];

			// DF(t, T*) reconstruction 
			dom_dff_star[j] = swp_f_df (xStr.dts[j], dom_model->lTStarDate, dom_yc);
			dom_dff_star[j] = log(dom_dff_star[j]);
			dom_gam1_star[j] = (1.0 - exp(dom_model->dLambdaX * (dom_model->dTStar - evt_tms[j]))) / dom_model->dLambdaX;
			dom_gam1_2_star[j] = 0.5 * dom_gam1_star[j] * dom_gam1_star[j];

			if (grfn_prm->num_dom_df > 0 && grfn_prm->dom_idx != -1)
			{
				grfn_prm->dom_dff = dvector (0, grfn_prm->num_dom_df-1);
				grfn_prm->dom_gam1 = dvector (0, grfn_prm->num_dom_df-1);
				grfn_prm->dom_gam1_2 = dvector (0, grfn_prm->num_dom_df-1);
				grfn_prm->dom_gam2 = dvector (0, grfn_prm->num_dom_df-1);
				grfn_prm->dom_gam2_2 = dvector (0, grfn_prm->num_dom_df-1);
				grfn_prm->dom_gam12 = dvector (0, grfn_prm->num_dom_df-1);				

				if (!grfn_prm->dom_dff || !grfn_prm->dom_gam1 || !grfn_prm->dom_gam1_2 || !grfn_prm->dom_gam2 || !grfn_prm->dom_gam2_2 || !grfn_prm->dom_gam12)
				{
					err = "Memory allocation error in SrtGrfnQTOLGMSV2FMC";
					goto FREE_RETURN;
				}

				for (k=0; k<grfn_prm->num_dom_df; k++)
				{
					grfn_prm->dom_dff[k] = swp_f_df (xStr.dts[j], grfn_prm->dom_df_dts[k], dom_yc);
					grfn_prm->dom_dff[k] = log(grfn_prm->dom_dff[k]);
					temp1 = (1.0 - exp(dom_model->dLambdaX * (dom_model->dTStar - evt_tms[j] - grfn_prm->dom_df_tms[k]))) / dom_model->dLambdaX;
					temp2 = (1.0 - exp(dom_model->dLambdaX2 * (dom_model->dTStar - evt_tms[j] - grfn_prm->dom_df_tms[k]))) / dom_model->dLambdaX2;
					grfn_prm->dom_gam1[k] = temp1 - dom_gam1_star[j];
					grfn_prm->dom_gam1_2[k] = 0.5 * temp1 * temp1 - dom_gam1_2_star[j];
				}

				grfn_prm->do_dom = 1;
			}
			else
			{
				grfn_prm->do_dom = 0;
			}

			// FX = DOMESTIC RECONSTRUCTION PRE-CALCULATIONS
			grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
			grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
			grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];

			if (grfn_prm->num_fx_df > 0 && grfn_prm->fx_idx != -1)
			{
				grfn_prm->fx_dff = dvector (0, grfn_prm->num_fx_df-1);
				grfn_prm->fx_gam1 = dvector (0, grfn_prm->num_fx_df-1);
				grfn_prm->fx_gam1_2 = dvector (0, grfn_prm->num_fx_df-1);
				grfn_prm->fx_gam2 = dvector (0, grfn_prm->num_fx_df-1);
				grfn_prm->fx_gam2_2 = dvector (0, grfn_prm->num_fx_df-1);
				grfn_prm->fx_gam12 = dvector (0, grfn_prm->num_fx_df-1);				

				if (!grfn_prm->fx_dff || !grfn_prm->fx_gam1 || !grfn_prm->fx_gam1_2 || !grfn_prm->fx_gam2 ||
					!grfn_prm->fx_gam2_2 || !grfn_prm->fx_gam12)
				{
					err = "Memory allocation error in SrtGrfnQTOLGMSV2FMC";
					goto FREE_RETURN;
				}

				for (k=0; k<grfn_prm->num_fx_df; k++)
				{
					grfn_prm->fx_dff[k] = swp_f_df (xStr.dts[j], grfn_prm->fx_df_dts[k], dom_yc);
					grfn_prm->fx_dff[k] = log(grfn_prm->fx_dff[k]);
					temp1 = (1.0 - exp(dom_model->dLambdaX * (dom_model->dTStar - evt_tms[j] - grfn_prm->fx_df_tms[k]))) / dom_model->dLambdaX;
					temp2 = (1.0 - exp(dom_model->dLambdaX2 * (dom_model->dTStar - evt_tms[j] - grfn_prm->fx_df_tms[k]))) / dom_model->dLambdaX2;
					grfn_prm->fx_gam1[k] = temp1 - dom_gam1_star[j];
					grfn_prm->fx_gam1_2[k] = 0.5 * temp1 * temp1 - dom_gam1_2_star[j];
				}

				grfn_prm->do_fx = 1;
			}
			else
			{
				grfn_prm->do_fx = 0;
			}

			// FOREIGN RECONSTRUCTION PRE-CALCULATIONS
			grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
			grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
			grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];

			// DF(t, T*) reconstruction 
			for_dff_star[j] = swp_f_df (xStr.dts[j], for_model->lTStarDate, for_yc);
			for_dff_star[j] = log(for_dff_star[j]);
			for_gam1_star[j] = (1.0 - exp(for_model->dLambdaX * (for_model->dTStar - evt_tms[j]))) / for_model->dLambdaX;
			for_gam2_star[j] = (1.0 - exp(for_model->dLambdaX2 * (for_model->dTStar - evt_tms[j]))) / for_model->dLambdaX2;
			for_gam12_star[j] = for_gam1_star[j] * for_gam2_star[j];
			for_gam1_2_star[j] = 0.5 * for_gam1_star[j] * for_gam1_star[j];
			for_gam2_2_star[j] = 0.5 * for_gam2_star[j] * for_gam2_star[j];

			if (grfn_prm->num_for_df > 0 && grfn_prm->for_idx != -1)
			{
				grfn_prm->for_dff = dvector (0, grfn_prm->num_for_df-1);
				grfn_prm->for_gam1 = dvector (0, grfn_prm->num_for_df-1);
				grfn_prm->for_gam1_2 = dvector (0, grfn_prm->num_for_df-1);
				grfn_prm->for_gam2 = dvector (0, grfn_prm->num_for_df-1);
				grfn_prm->for_gam2_2 = dvector (0, grfn_prm->num_for_df-1);
				grfn_prm->for_gam12 = dvector (0, grfn_prm->num_for_df-1);				

				if (!grfn_prm->for_dff || !grfn_prm->for_gam1 || !grfn_prm->for_gam1_2 || !grfn_prm->for_gam2 || !grfn_prm->for_gam2_2 || !grfn_prm->for_gam12)
				{
					err = "Memory allocation error in SrtGrfnQTOLGMSV2FMC";
					goto FREE_RETURN;
				}

				for (k=0; k<grfn_prm->num_for_df; k++)
				{
					grfn_prm->for_dff[k] = swp_f_df (xStr.dts[j], grfn_prm->for_df_dts[k], for_yc);
					grfn_prm->for_dff[k] = log(grfn_prm->for_dff[k]);
					temp1 = (1.0 - exp(for_model->dLambdaX * (for_model->dTStar - evt_tms[j] - grfn_prm->for_df_tms[k]))) / for_model->dLambdaX;
					temp2 = (1.0 - exp(for_model->dLambdaX2 * (for_model->dTStar - evt_tms[j] - grfn_prm->for_df_tms[k]))) / for_model->dLambdaX2;
					grfn_prm->for_gam1[k] = temp1 - for_gam1_star[j];
					grfn_prm->for_gam2[k] = temp2 - for_gam2_star[j];					
					grfn_prm->for_gam1_2[k] = 0.5 * temp1 * temp1 - for_gam1_2_star[j];
					grfn_prm->for_gam2_2[k] = 0.5 * temp2 * temp2 - for_gam2_2_star[j];
					grfn_prm->for_gam12[k] = temp1 * temp2 - for_gam12_star[j];
				}

				grfn_prm->do_for = 1;
			}
			else
			{
				grfn_prm->do_for = 0;
			}

			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;

			j--;
			while (j >= 0 && xStr.evt[j].evt == NULL)
			{
				j--;
			}
		}
		else if (j>=0 && xStr.am[j])
		{
			grfn_prm = malloc(sizeof(grfn_parm_mc_fxlgmsv));
			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j;
			grfn_prm->fx_idx = fx_idx;
			grfn_prm->dom_idx = dom_idx;
			grfn_prm->for_idx = for_idx;

			grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
			grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
			grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];

			grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
			grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
			grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];

			grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
			grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
			grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];

			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;
		}
		else
		{
			is_event[i] = 0;
			void_prm[i] = NULL;
		}
	}

	//	Eventually! call to function 
	
	*nb_prod = num_col;	
	*tableauCols = num_col;
	*resRows = max(num_col + 1, num_evt);

	if (do_optimisation)
	{		
		if (params->iMultiIndex)
		{
			params->iNbIndex = params->iMultiIndex;
		}
		else
		{
			params->iNbIndex = 1;
		}

		err = mceb_allocate_params(	params,
									num_evt);

		if (err) goto FREE_RETURN;

		*prod_val = dmatrix (0, *resRows - 1, 0, 2 + params->iNbIndex);
	}
	else
	{
		*prod_val = dmatrix (0, *nb_prod - 1, 0, 2);
	}

	err = qtolgmsv2f_mc_balsam(nstp,
							num_evt,
							time,
							date,
							
							numpaths,
							
							dom_model->dLambdaX,
							
							dom_sigma,
							dom_zcvol_star,
							
							dom_dff_star,
							dom_gam1_star,
							dom_gam1_2_star,

							for_model->dLambdaX,
							for_model->dLambdaX2,
							
							for_sigma,
							for_newalphaLGM,
							for_newrhoLGM,
							for_alpha,
							for_lameps,
							for_lvleps,
							for_rho,
							for_rho2,
							
							for_zcvol1_star,
							for_zcvol2_star,

							for_dff_star,
							for_gam1_star,
							for_gam2_star,
							for_gam1_2_star,
							for_gam2_2_star,
							for_gam12_star,
							
							fx_sigma,
							
							correlation,
							
							void_prm, 
							is_event,

							do_optimisation,
							optimise,
							params,
							
							NULL,
							
							payoff_qtolgmsv2f_mc,
							
							num_col,
							*prod_val);

	if (err) goto FREE_RETURN;

	if (do_optimisation)
	{
		/* Recopy Barrier / CoefLin for the moment */
		for (i=0; i<num_evt; i++)
		{
			(*prod_val)[i][2] = params->dBarrier[i];

			for (j=0; j<params->iNbIndex; j++)
			{
				(*prod_val)[i][3+j] = params->dCoefLin[i][j+1];
			}
		}
	}

	dom_df = swp_f_df (today, dom_model->lTStarDate, dom_yc);

	for (i=0; i<num_col; i++)
	{
		(*prod_val)[i][0] *= dom_df;
		(*prod_val)[i][1] *= dom_df;
	}

	if (do_optimisation)
	{
		(*prod_val)[num_col][0] *= dom_df;
		(*prod_val)[num_col][1] *= dom_df;
	}

	if (flag)
	{		
		for (i=0; i<num_evt-1; i++)
		{			
			(*prod_val)[i][2] = (*prod_val)[i+1][2];			
			
			if (do_optimisation)
			{
				for (k=0; k<params->iNbIndex; k++)
				{
					(*prod_val)[i][3+k] = (*prod_val)[i+1][3+k];
				}
			}
		}
	}
	
	//	Add PV of Past 
	(*prod_val)[num_col-1][0] += xStr.gd->pv_of_past;

		
FREE_RETURN:
	
	if (free_str)
	{
		FIRSTFreeUndFromDeal(
									num_und,
									&und_ptr
								   );

		FIRSTFreeEvtDatesFromDeal(
										nstp,
										&evt_dts,
										&evt_tms
										);

		FIRSTFreeMktStruct(&xStr);
	}

	if (void_prm)
	{
		for (i=0; i<nstp; i++)
		{
			if (void_prm[i])
			{
				grfn_prm = (GRFNPARMMCFXLGMSV) void_prm[i];
				free (grfn_prm);
			}
		}

		free (void_prm);
	}	

	if (time) free (time);
	if (is_event) free (is_event);
	if (date) free (date);	

	if (fx_sigma) free (fx_sigma);

	// TEST FXLGMSV
	// if (correlation) free_f3tensor(correlation, 0, nstp-1, 0, 6, 0, 6);

	// Original Call
	if (correlation) free_f3tensor(correlation, 0, nstp-1, 0, 4, 0, 4);

	if (dom_zcvol_star) free (dom_zcvol_star);
	if (for_zcvol1_star) free (for_zcvol1_star);
	if (for_zcvol2_star) free (for_zcvol2_star);

	if (dom_sigma) free (dom_sigma);

	if (dom_dff_star) free(dom_dff_star);
	if (dom_gam1_star) free(dom_gam1_star);
	if (dom_gam1_2_star) free(dom_gam1_2_star);

	if (for_sigma) free (for_sigma);
	if (for_alpha) free (for_alpha);
	if (for_rho) free (for_rho);
	if (for_rho2) free (for_rho2);
	if (for_lameps) free (for_lameps);
	if (for_lvleps) free (for_lvleps);
	if (for_newalphaLGM) free (for_newalphaLGM);
	if (for_newrhoLGM) free (for_newrhoLGM);

	if (for_dff_star) free(for_dff_star);
	if (for_gam1_star) free(for_gam1_star);
	if (for_gam2_star) free(for_gam2_star);
	if (for_gam1_2_star) free(for_gam1_2_star);
	if (for_gam2_2_star) free(for_gam2_2_star);
	if (for_gam12_star) free(for_gam12_star);

	if(fx_time)
	{
		free(fx_time);
		fx_time = NULL;
	}
	
	if(fx_vol)
	{
		free(fx_vol);
		fx_vol = NULL;
	}

	if(rho_time)
	{
		free(rho_time);
		rho_time = NULL;
	}

	if(CorrMatrix)
	{
		free_f3tensor(CorrMatrix, 0, num_rho-1, 0, 6, 0, 6);
		CorrMatrix = NULL;
	}

	free_LGMSV_model(dom_model);
	free_LGMSV_model(for_model);

	mceb_free_params(params);
	
	return err;
}



char *SrtGrfnQTOLGMSV1FMC(
					char		*underlying,
					int			numeventdates,
					long		*eventdates,
					long		tableauRows,
					long		*tableauCols,
					char		***tableauStrings,
					int			**tableauMask,
					long		auxWidth,
					long		*auxLen,
					double		**aux,

					// for Optimisation of exercise boundary 
					int			do_optimisation,
					int			*optimise,
					MCEBPARAMS	params,
					long		*resRows,

					int			nstept,
					int			numpaths,				  
					int			*nb_prod,
					double		***prod_val)
{
	int					free_str = 0;
	FIRSTAllMkts		xStr;
	SrtGrfnParam		defParm;
	GRFNPARMMCFXLGMSV	grfn_prm;
	int					forback;
	int					flag = 0;
	long				nstp;	
						
	double				*evt_tms		= NULL,						
						*time			= NULL,
						*date			= NULL,

						*dom_sigma			= NULL,
						*dom_zcvol_star	= NULL,

						*dom_dff_star		= NULL,
						*dom_gam1_star		= NULL,
						*dom_gam1_2_star	= NULL,

						*for_sigma			= NULL,
						*for_alpha			= NULL,
						*for_rho			= NULL,
						*for_lameps			= NULL,
						*for_lvleps			= NULL,
						*for_newalphaLGM	= NULL,
						*for_newrhoLGM		= NULL,
						
						*for_zcvol_star	= NULL,

						*for_dff_star		= NULL,
						*for_gam1_star		= NULL,
						*for_gam1_2_star	= NULL;

	long				*evt_dts	= NULL;
	
	int					*is_event	= NULL;
	void				**void_prm	= NULL;
		
	long				index;
	double				final_mat;

	long				today, spot_date;
						
	int					num_col		= 0,
						max_num_df	= 0,
						num_evt		= 0,
						num_und		= 0;
						
	char				*domname,
						*forname;

	char				*dom_yc, *for_yc;

	SrtUndPtr			fx_und, dom_und, for_und;
	SrtUndPtr			*und_ptr	= NULL,
						und			= NULL;

	LGMSV_model			dom_model_, *dom_model = &dom_model_;
	LGMSV_model			for_model_, *for_model = &for_model_;

	int					i, j, k, l;

	int					fx_idx, dom_idx, for_idx;

	double				dom_df;
	double				temp1, temp2;
	
	Err					err = NULL;

	double				fx_spot;

	int					num_fx_vol;
	double				*fx_time=NULL;
	double				*fx_vol=NULL;

	int					num_rho;
	double				*rho_time=NULL;
	double				***CorrMatrix=NULL;

	double				*fx_sigma=NULL;
	double				***correlation=NULL;

	//	Initialisation of the LGMSV model 
	init_NULL_LGMSV_model(dom_model);
	init_NULL_LGMSV_model(for_model);

	//	Initialise the GRFN tableau 

	//	First, initialise the param struct 	
	err = srt_f_set_default_GrfnParams (&defParm);	
	defParm.force_mc = 0;

	err = FIRSTInitMktStruct(
							numeventdates,
							eventdates,
							tableauRows,
							*tableauCols,
							tableauStrings,
							tableauMask,
							auxWidth,
							auxLen,
							aux,
							underlying,
							&defParm,
							&forback,
							&xStr);
	if (err)
	{
		goto FREE_RETURN;
	}

	free_str = 1;

	/*	Now, lookup underlyings involved and their term structures */

	fx_und = lookup_und (underlying);
	
	if (!fx_und)
	{
		err = serror ("Couldn't find underlying named %s", underlying);
		goto FREE_RETURN;
	}
	
	today = get_today_from_underlying (fx_und);

	if (get_underlying_type (fx_und) != FOREX_UND)
	{
		err = serror ("Underlying %s is not of type FX", underlying);
		goto FREE_RETURN;
	}

	if (get_mdltype_from_fxund (fx_und) != FX_LGMSV)
	{
		err = serror ("Underlying %s is not of type FX Stoch Rates", underlying);
		goto FREE_RETURN;
	}

	fxlgmsv_get_fx_and_correl_ts(underlying, &fx_spot, &num_fx_vol, &fx_time, &fx_vol, &num_rho, &rho_time, &CorrMatrix);

	domname = get_domname_from_fxund (fx_und);
	dom_und = lookup_und (domname);
	if (!dom_und)
	{
		err = serror ("Couldn't find underlying named %s", domname);
		goto FREE_RETURN;
	}
	// Get the model 
	err = Get_LGMSV_model(	domname,
							dom_model);

	forname = get_forname_from_fxund (fx_und);
	for_und = lookup_und (forname);
	if (!for_und)
	{
		err = serror ("Couldn't find underlying named %s", forname);
		goto FREE_RETURN;
	}
	// Get the model 
	err = Get_LGMSV_model(	forname,
							for_model);

	// look for the today date 
	today = get_today_from_underlying (dom_und);
	spot_date = add_unit (today, 2, SRT_BDAY, MODIFIED_SUCCEEDING);

	dom_yc = (char *) get_ycname_from_irund (dom_und);
	for_yc = (char *) get_ycname_from_irund (for_und);
	
	// Get number of columns 
	err = FIRSTGetNumColFromDeal(
									&xStr,
									&num_col);	
	if (err)
	{
		goto FREE_RETURN;
	}

	//	Get the maximum number of dfs required	
	err = FIRSTGetMaxNumDfFromDeal(
									&xStr,
									&max_num_df);
	if (err)
	{
		goto FREE_RETURN;
	}
	
	//	Next, get the time steps 
	err = FIRSTGetEvtDatesFromDeal(
									&xStr,
 									&num_evt,
 									&evt_dts,
									&evt_tms);
	if (err)
	{
		goto FREE_RETURN;
	}
	

	if (err) goto FREE_RETURN;
																					
	// discretise in time			

	nstp = num_evt;

	time = (double*) calloc (nstp, sizeof (double));
	if (!time)
	{
		err = "Memory allocation error (1) SrtGrfnQTOLGMSV1FMC";
		goto FREE_RETURN;
	}
	memcpy (time, xStr.tms, nstp * sizeof (double));

	final_mat = xStr.tms[nstp-1];
	i = 0;
	while (i < dom_model->iNbPWTime && dom_model->dPWTime[i] < final_mat)
	{
		num_f_add_number(&nstp, &time, dom_model->dPWTime[i]);
		i++;
	}
	i = 0;
	while (i < for_model->iNbPWTime && for_model->dPWTime[i] < final_mat)
	{
		num_f_add_number(&nstp, &time, for_model->dPWTime[i]);
		i++;
	}
	i = 0;
	while (i < num_fx_vol && fx_time[i] < final_mat)
	{
		num_f_add_number(&nstp, &time, fx_time[i]);
		i++;
	}
	i = 0;
	while (i < num_rho && rho_time[i] < final_mat)
	{
		num_f_add_number(&nstp, &time, rho_time[i]);
		i++;
	}
	
	num_f_sort_vector (nstp, time);

	//	Fill the time vector 
	err = fill_time_vector (	&time, 
								&nstp,
								0,
								NULL,
								0,
								NULL, 
								nstept);				
	if (err)
	{
		goto FREE_RETURN;
	}

	/*	Fill product structure */

	strupper (underlying);
	strip_white_space (underlying);
	strupper (domname);
	strip_white_space (domname);
	strupper (forname);
	strip_white_space (forname);
	for (i=0; i<xStr.num_und; i++)
	{
		strupper (xStr.und_data[i].und_name);
		strip_white_space (xStr.und_data[i].und_name);
	}

	fx_idx = -1;
	for (i=0; i<xStr.num_und; i++)
	{
		if (!strcmp (xStr.und_data[i].und_name, underlying))
		{
			fx_idx = i;
		}
	}
	if (fx_idx == -1)
	{
		err = "The Fx underlying is not present in the mdlcomm structure";
		goto FREE_RETURN;
	}

	dom_idx = -1;
	for (i=0; i<xStr.num_und; i++)
	{
		if (!strcmp (xStr.und_data[i].und_name, domname))
		{
			dom_idx = i;
		}
	}
	if (dom_idx == -1)
	{
		err = "The domestic underlying is not present in the mdlcomm structure";
		goto FREE_RETURN;
	}

	for_idx = -1;
	for (i=0; i<xStr.num_und; i++)
	{
		if (!strcmp (xStr.und_data[i].und_name, forname))
		{
			for_idx = i;
		}
	}
	if (for_idx == -1)
	{
		err = "The foreign underlying is not present in the mdlcomm structure";
		goto FREE_RETURN;
	}

	void_prm = (void**) calloc (nstp, sizeof (void*));
	is_event = (int *) calloc (nstp, sizeof (int));
	date = (double *) calloc (nstp, sizeof (double));

	dom_sigma = (double *) calloc (nstp, sizeof (double));	
	dom_zcvol_star = (double *) calloc(nstp, sizeof(double));
	
	dom_dff_star = (double *) calloc(num_evt, sizeof(double));
	dom_gam1_star = (double *) calloc(num_evt, sizeof(double));
	dom_gam1_2_star = (double *) calloc(num_evt, sizeof(double));
	
	for_sigma = (double *) calloc (nstp, sizeof (double));	
	for_alpha = (double *) calloc (nstp, sizeof (double));
	for_rho = (double *) calloc (nstp, sizeof (double));	
	for_lameps = (double *) calloc (nstp, sizeof (double));
	for_lvleps = (double *) calloc (nstp, sizeof (double));
	for_newalphaLGM = (double *) calloc (nstp, sizeof (double));
	for_newrhoLGM = (double *) calloc (nstp, sizeof (double));
	
	for_zcvol_star = (double *) calloc(nstp, sizeof(double));

	for_dff_star = (double *) calloc(num_evt, sizeof(double));
	for_gam1_star = (double *) calloc(num_evt, sizeof(double));
	for_gam1_2_star = (double *) calloc(num_evt, sizeof(double));

	fx_sigma = (double *) calloc(nstp, sizeof(double));
	correlation = f3tensor(0, nstp-1, 0, 3, 0, 3);

	if (!void_prm || !is_event || !date ||
		!dom_sigma ||
		!fx_sigma || !correlation ||
		!dom_zcvol_star || !for_zcvol_star ||
		!dom_dff_star || !dom_gam1_star || !dom_gam1_2_star ||
		!for_sigma || !for_alpha || !for_rho || !for_lameps || !for_lvleps || !for_newalphaLGM || !for_newrhoLGM ||
		!for_dff_star || !for_gam1_star || !for_gam1_2_star )
	{
		err = "Memory allocation failure in SrtGrfnQTOLGMSV1FMC";
		goto FREE_RETURN;
	}

	j = xStr.num_evt - 1;	

	for (i=nstp-1; i>=0; i--)
	{
		date[i] = ((long) (today + time[i] * 365.0 + 1.0E-08)) * 1.0;	
		
		if (i < nstp-1)
		{
			index = Get_Index(time[i+1], dom_model->dPWTime, dom_model->iNbPWTime);

			dom_sigma[i+1] = dom_model->dSigma[index];

			// J.M.L 10 Mar 2004
			// dom_zcvol_star represents the term (1 - exp(lbda1 * ( Tstar - t )) / lbda1 which appears in
			// the volatility of the bond B(t,Tstar). We can either code it directly:

			// Previous code:
			// dom_zcvol_star[i+1] = (1.0 - exp(dom_model->dLambdaX * (dom_model->dTStar - time[i]))) / dom_model->dLambdaX;
			
			// but a better plan is to integrate the square of this expression which has the dimension of a volatility
			// and find an equivalent volatility for the bond which is valid between ti and ti+1
			
			dom_zcvol_star[i+1] = -sqrt( 
				( time[i+1]-time[i] 
				+ 2*(exp(dom_model->dLambdaX * (dom_model->dTStar - time[i+1])) 
						- exp(dom_model->dLambdaX * (dom_model->dTStar - time[i])))/dom_model->dLambdaX
				+ 0.5*(exp(2*dom_model->dLambdaX*(dom_model->dTStar-time[i])) 
						- exp(2*dom_model->dLambdaX*(dom_model->dTStar-time[i+1])))/dom_model->dLambdaX
				) / (time[i+1] - time[i]) ) / dom_model->dLambdaX;

			// Previous code:
			// for_zcvol_star[i+1] = (1.0 - exp(for_model->dLambdaX * (for_model->dTStar - time[i]))) / for_model->dLambdaX;

			for_zcvol_star[i+1] = -sqrt( 
				( time[i+1]-time[i] 
				+ 2*(exp(for_model->dLambdaX * (for_model->dTStar - time[i+1])) 
						- exp(for_model->dLambdaX * (for_model->dTStar - time[i])))/for_model->dLambdaX
				+ 0.5*(exp(2*for_model->dLambdaX*(for_model->dTStar-time[i])) 
						- exp(2*for_model->dLambdaX*(for_model->dTStar-time[i+1])))/for_model->dLambdaX
				) / (time[i+1] - time[i]) ) / for_model->dLambdaX;

			index = Get_Index(time[i+1], for_model->dPWTime, for_model->iNbPWTime);

			for_sigma[i+1] = for_model->dSigma[index];
			for_newalphaLGM[i+1] = for_model->dLGMAlpha[index];
			for_newrhoLGM[i+1] = for_model->dLGMRho[index];
			
			for_alpha[i+1] = for_model->dAlpha[index];
			for_rho[i+1] = for_model->dRho[index];
			for_lameps[i+1] = for_model->dLambdaEps[index];
			for_lvleps[i+1] = for_model->dLvlEps[index];

			index = Get_Index(time[i+1], fx_time, num_fx_vol);
			fx_sigma[i+1] = fx_vol[index];

			index = Get_Index(time[i+1], rho_time, num_rho);

			correlation[i+1][0][0] = CorrMatrix[index][0][0];
			correlation[i+1][0][1] = CorrMatrix[index][0][3];
			correlation[i+1][0][2] = CorrMatrix[index][0][5];
			correlation[i+1][0][3] = CorrMatrix[index][0][6];

			correlation[i+1][1][1] = CorrMatrix[index][3][3];
			correlation[i+1][1][2] = CorrMatrix[index][3][5];
			correlation[i+1][1][3] = CorrMatrix[index][3][6];

			correlation[i+1][2][2] = CorrMatrix[index][5][5];
			correlation[i+1][2][3] = CorrMatrix[index][5][6];

			correlation[i+1][3][3] = CorrMatrix[index][6][6];
			for(l=0;l<4;++l)
			{
				for(k=0;k<l;++k)
				{
					correlation[i+1][l][k] = correlation[i+1][k][l];
				}
			}
		}

		if (j >=0 && fabs(time[i] - evt_tms[j]) < 1.0E-08)
		{
			grfn_prm = malloc(sizeof(grfn_parm_mc_fxlgmsv));

			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j;
			grfn_prm->fx_idx = fx_idx;
			grfn_prm->dom_idx = dom_idx;
			grfn_prm->for_idx = for_idx;

			grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
			grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
			grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];

			// DF(t, T*) reconstruction 
			dom_dff_star[j] = swp_f_df (xStr.dts[j], dom_model->lTStarDate, dom_yc);
			dom_dff_star[j] = log(dom_dff_star[j]);
			dom_gam1_star[j] = (1.0 - exp(dom_model->dLambdaX * (dom_model->dTStar - evt_tms[j]))) / dom_model->dLambdaX;
			dom_gam1_2_star[j] = 0.5 * dom_gam1_star[j] * dom_gam1_star[j];

			if (grfn_prm->num_dom_df > 0 && grfn_prm->dom_idx != -1)
			{
				grfn_prm->dom_dff = dvector (0, grfn_prm->num_dom_df-1);
				grfn_prm->dom_gam1 = dvector (0, grfn_prm->num_dom_df-1);
				grfn_prm->dom_gam1_2 = dvector (0, grfn_prm->num_dom_df-1);
				grfn_prm->dom_gam2 = dvector (0, grfn_prm->num_dom_df-1);
				grfn_prm->dom_gam2_2 = dvector (0, grfn_prm->num_dom_df-1);
				grfn_prm->dom_gam12 = dvector (0, grfn_prm->num_dom_df-1);				

				if (!grfn_prm->dom_dff || !grfn_prm->dom_gam1 || !grfn_prm->dom_gam1_2 || !grfn_prm->dom_gam2 || !grfn_prm->dom_gam2_2 || !grfn_prm->dom_gam12)
				{
					err = "Memory allocation error in SrtGrfnQTOLGMSV1FMC";
					goto FREE_RETURN;
				}

				for (k=0; k<grfn_prm->num_dom_df; k++)
				{
					grfn_prm->dom_dff[k] = swp_f_df (xStr.dts[j], grfn_prm->dom_df_dts[k], dom_yc);
					grfn_prm->dom_dff[k] = log(grfn_prm->dom_dff[k]);
					temp1 = (1.0 - exp(dom_model->dLambdaX * (dom_model->dTStar - evt_tms[j] - grfn_prm->dom_df_tms[k]))) / dom_model->dLambdaX;
					temp2 = (1.0 - exp(dom_model->dLambdaX2 * (dom_model->dTStar - evt_tms[j] - grfn_prm->dom_df_tms[k]))) / dom_model->dLambdaX2;
					grfn_prm->dom_gam1[k] = temp1 - dom_gam1_star[j];
					grfn_prm->dom_gam1_2[k] = 0.5 * temp1 * temp1 - dom_gam1_2_star[j];
				}

				grfn_prm->do_dom = 1;
			}
			else
			{
				grfn_prm->do_dom = 0;
			}

		
			grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
			grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
			grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];

			if (grfn_prm->num_fx_df > 0 && grfn_prm->fx_idx != -1)
			{
				grfn_prm->fx_dff = dvector (0, grfn_prm->num_fx_df-1);
				grfn_prm->fx_gam1 = dvector (0, grfn_prm->num_fx_df-1);
				grfn_prm->fx_gam1_2 = dvector (0, grfn_prm->num_fx_df-1);
				grfn_prm->fx_gam2 = dvector (0, grfn_prm->num_fx_df-1);
				grfn_prm->fx_gam2_2 = dvector (0, grfn_prm->num_fx_df-1);
				grfn_prm->fx_gam12 = dvector (0, grfn_prm->num_fx_df-1);				

				if (!grfn_prm->fx_dff || !grfn_prm->fx_gam1 || !grfn_prm->fx_gam1_2 || !grfn_prm->fx_gam2 ||
					!grfn_prm->fx_gam2_2 || !grfn_prm->fx_gam12)
				{
					err = "Memory allocation error in SrtGrfnQTOLGMSV1FMC";
					goto FREE_RETURN;
				}

				for (k=0; k<grfn_prm->num_fx_df; k++)
				{
					grfn_prm->fx_dff[k] = swp_f_df (xStr.dts[j], grfn_prm->fx_df_dts[k], dom_yc);
					grfn_prm->fx_dff[k] = log(grfn_prm->fx_dff[k]);
					temp1 = (1.0 - exp(dom_model->dLambdaX * (dom_model->dTStar - evt_tms[j] - grfn_prm->fx_df_tms[k]))) / dom_model->dLambdaX;
					temp2 = (1.0 - exp(dom_model->dLambdaX2 * (dom_model->dTStar - evt_tms[j] - grfn_prm->fx_df_tms[k]))) / dom_model->dLambdaX2;
					grfn_prm->fx_gam1[k] = temp1 - dom_gam1_star[j];
					grfn_prm->fx_gam1_2[k] = 0.5 * temp1 * temp1 - dom_gam1_2_star[j];
				}
				grfn_prm->do_fx = 1;
			}
			else
			{
				grfn_prm->do_fx = 0;
			}

			grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
			grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
			grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];

			// DF(t, T*) reconstruction 
			for_dff_star[j] = swp_f_df (xStr.dts[j], for_model->lTStarDate, for_yc);
			for_dff_star[j] = log(for_dff_star[j]);
			for_gam1_star[j] = (1.0 - exp(for_model->dLambdaX * (for_model->dTStar - evt_tms[j]))) / for_model->dLambdaX;
			for_gam1_2_star[j] = 0.5 * for_gam1_star[j] * for_gam1_star[j];

			if (grfn_prm->num_for_df > 0 && grfn_prm->for_idx != -1)
			{
				grfn_prm->for_dff = dvector (0, grfn_prm->num_for_df-1);
				grfn_prm->for_gam1 = dvector (0, grfn_prm->num_for_df-1);
				grfn_prm->for_gam1_2 = dvector (0, grfn_prm->num_for_df-1);
				grfn_prm->for_gam2 = dvector (0, grfn_prm->num_for_df-1);
				grfn_prm->for_gam2_2 = dvector (0, grfn_prm->num_for_df-1);
				grfn_prm->for_gam12 = dvector (0, grfn_prm->num_for_df-1);				

				if (!grfn_prm->for_dff || !grfn_prm->for_gam1 || !grfn_prm->for_gam1_2 || !grfn_prm->for_gam2 || !grfn_prm->for_gam2_2 || !grfn_prm->for_gam12)
				{
					err = "Memory allocation error in SrtGrfnQTOLGMSV1FMC";
					goto FREE_RETURN;
				}

				for (k=0; k<grfn_prm->num_for_df; k++)
				{
					grfn_prm->for_dff[k] = swp_f_df (xStr.dts[j], grfn_prm->for_df_dts[k], for_yc);
					grfn_prm->for_dff[k] = log(grfn_prm->for_dff[k]);
					temp1 = (1.0 - exp(for_model->dLambdaX * (for_model->dTStar - evt_tms[j] - grfn_prm->for_df_tms[k]))) / for_model->dLambdaX;
					temp2 = (1.0 - exp(for_model->dLambdaX2 * (for_model->dTStar - evt_tms[j] - grfn_prm->for_df_tms[k]))) / for_model->dLambdaX2;
					grfn_prm->for_gam1[k] = temp1 - for_gam1_star[j];
					grfn_prm->for_gam1_2[k] = 0.5 * temp1 * temp1 - for_gam1_2_star[j];
				}
				grfn_prm->do_for = 1;
			}
			else
			{
				grfn_prm->do_for = 0;
			}

			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;

			j--;
			while (j >= 0 && xStr.evt[j].evt == NULL)
			{
				j--;
			}
		}
		else if (j>=0 && xStr.am[j])
		{
			grfn_prm = malloc(sizeof(grfn_parm_mc_fxlgmsv));
			grfn_prm->global = &xStr;
			grfn_prm->local = xStr.evt + j;
			grfn_prm->fx_idx = fx_idx;
			grfn_prm->dom_idx = dom_idx;
			grfn_prm->for_idx = for_idx;

			grfn_prm->num_dom_df = xStr.evt[j].evt->dflen[dom_idx];
			grfn_prm->dom_df_tms = xStr.evt[j].evt->dft[dom_idx];
			grfn_prm->dom_df_dts = xStr.evt[j].evt->dfd[dom_idx];

			grfn_prm->num_fx_df = xStr.evt[j].evt->dflen[fx_idx];
			grfn_prm->fx_df_tms = xStr.evt[j].evt->dft[fx_idx];
			grfn_prm->fx_df_dts = xStr.evt[j].evt->dfd[fx_idx];

			grfn_prm->num_for_df = xStr.evt[j].evt->dflen[for_idx];
			grfn_prm->for_df_tms = xStr.evt[j].evt->dft[for_idx];
			grfn_prm->for_df_dts = xStr.evt[j].evt->dfd[for_idx];

			is_event[i] = 1;
			void_prm[i] = (void*) grfn_prm;
		}
		else
		{
			is_event[i] = 0;
			void_prm[i] = NULL;
		}
	}

	//	Eventually! call to function 
	
	*nb_prod = num_col;	
	*tableauCols = num_col;
	*resRows = max(num_col + 1, num_evt);

	if (do_optimisation)
	{		
		if (params->iMultiIndex)
		{
			params->iNbIndex = params->iMultiIndex;
		}
		else
		{
			params->iNbIndex = 1;
		}

		err = mceb_allocate_params(	params,
									num_evt);

		if (err) goto FREE_RETURN;

		*prod_val = dmatrix (0, *resRows - 1, 0, 2 + params->iNbIndex);
	}
	else
	{
		*prod_val = dmatrix (0, *nb_prod - 1, 0, 2);
	}
			
	err = qtolgmsv1f_mc_balsam(nstp,
							num_evt,
							time,
							date,
							
							numpaths,
							
							dom_model->dLambdaX,
							
							dom_sigma,
							dom_zcvol_star,
							
							dom_dff_star,
							dom_gam1_star,
							dom_gam1_2_star,

							for_model->dLambdaX,
							
							for_sigma,
							for_alpha,
							for_lameps,
							for_lvleps,
							for_rho,
							
							for_zcvol_star,

							for_dff_star,
							for_gam1_star,
							for_gam1_2_star,
							
							fx_sigma,
							
							correlation,
							
							void_prm, 
							is_event,

							do_optimisation,
							optimise,
							params,
							
							NULL,
							
							payoff_qtolgmsv1f_mc,
							
							num_col,
							*prod_val);
	if (err) goto FREE_RETURN;

	if (do_optimisation)
	{
		/* Recopy Barrier / CoefLin for the moment */
		for (i=0; i<num_evt; i++)
		{
			(*prod_val)[i][2] = params->dBarrier[i];

			for (j=0; j<params->iNbIndex; j++)
			{
				(*prod_val)[i][3+j] = params->dCoefLin[i][j+1];
			}
		}
	}

	dom_df = swp_f_df (today, dom_model->lTStarDate, dom_yc);

	for (i=0; i<num_col; i++)
	{
		(*prod_val)[i][0] *= dom_df;
		(*prod_val)[i][1] *= dom_df;
	}

	if (do_optimisation)
	{
		(*prod_val)[num_col][0] *= dom_df;
		(*prod_val)[num_col][1] *= dom_df;
	}

	if (flag)
	{		
		for (i=0; i<num_evt-1; i++)
		{			
			(*prod_val)[i][2] = (*prod_val)[i+1][2];			
			
			if (do_optimisation)
			{
				for (k=0; k<params->iNbIndex; k++)
				{
					(*prod_val)[i][3+k] = (*prod_val)[i+1][3+k];
				}
			}
		}
	}
	
	//	Add PV of Past 
	(*prod_val)[num_col-1][0] += xStr.gd->pv_of_past;

		
FREE_RETURN:
	
	if (free_str)
	{
		FIRSTFreeUndFromDeal(
									num_und,
									&und_ptr
								   );

		FIRSTFreeEvtDatesFromDeal(
										nstp,
										&evt_dts,
										&evt_tms
										);

		FIRSTFreeMktStruct(&xStr);
	}

	if (void_prm)
	{
		for (i=0; i<nstp; i++)
		{
			if (void_prm[i])
			{
				grfn_prm = (GRFNPARMMCFXLGMSV) void_prm[i];
				free (grfn_prm);
			}
		}

		free (void_prm);
	}	

	if (time) free (time);
	if (is_event) free (is_event);
	if (date) free (date);	

	if (fx_sigma) free (fx_sigma);
	if (correlation) free_f3tensor(correlation, 0, nstp-1, 0, 3, 0, 3);

	if (dom_zcvol_star) free (dom_zcvol_star);
	if (for_zcvol_star) free (for_zcvol_star);

	if (dom_sigma) free (dom_sigma);

	if (dom_dff_star) free(dom_dff_star);
	if (dom_gam1_star) free(dom_gam1_star);
	if (dom_gam1_2_star) free(dom_gam1_2_star);

	if (for_sigma) free (for_sigma);
	if (for_alpha) free (for_alpha);
	if (for_rho) free (for_rho);
	if (for_lameps) free (for_lameps);
	if (for_lvleps) free (for_lvleps);

	if (for_dff_star) free(for_dff_star);
	if (for_gam1_star) free(for_gam1_star);
	if (for_gam1_2_star) free(for_gam1_2_star);

	if(fx_time)
	{
		free(fx_time);
		fx_time = NULL;
	}
	
	if(fx_vol)
	{
		free(fx_vol);
		fx_vol = NULL;
	}

	if(rho_time)
	{
		free(rho_time);
		rho_time = NULL;
	}

	if(CorrMatrix)
	{
		free_f3tensor(CorrMatrix, 0, num_rho-1, 0, 6, 0, 6);
		CorrMatrix = NULL;
	}

	free_LGMSV_model(dom_model);
	free_LGMSV_model(for_model);
	
	return err;
}

