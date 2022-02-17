#include "swp_h_all.h"
#include "srt_h_all.h"
#include "opfnctns.h"
#include "srtaccess.h"
#include "math.h"
#include "BGMUtils.h"
#include "swp_h_vol.h"

//******************************************************************************
/* Closed from pricing of a future in LGM using a calibrated LGM               */
//******************************************************************************

Err srt_f_price_future_lgm (
			long			  today,
			long              *fut_dates, 
			long              num_fut_dates,
			char              *yc_name,
			char              *ref_rate_code,
			char              *calibrated_underlying,
			int				  result_type,
			double            *fut_price)
{
Err			       err = NULL;
int				   i;
Date *			   end_futdates = NULL;
double			   df_start, df_end;
double             spread;
SrtCompounding     float_compounding;
SrtBasisCode       float_basis;
char			   und_name[256];
double*			   price = NULL;
SrtUndPtr		   und=NULL;
double fwd;

// Check if the underlying is empty
	und = lookup_und(calibrated_underlying);
	if(!und)
	{
		err = "srt_f_price_future_lgm cannot find the underlying";
		goto e1;
	}

/* Get the reference rate details (compounding and basis) */
	err = swp_f_get_ref_rate_details(ref_rate_code, &float_basis, &float_compounding);
	if (err)
		goto e1;

// Memory allocation (useful for calibration
	end_futdates   = (Date *) calloc(num_fut_dates, sizeof(Date));

//Compute the end dates of the futures
	for (i = 0 ; i < num_fut_dates ; i++)
	{
	// End dates of each future
		end_futdates[i] = add_unit(fut_dates[i], (int)(12/float_compounding), 
			SRT_MONTH, MODIFIED_SUCCEEDING);
	}

// ------------------------------------------------------------------------------------
//The Pricing of futures according to the closed form for LGM
//------------------------------------------------------------------------------------- 

// Underlying name
	strcpy(und_name, calibrated_underlying);


		for (i = 0 ; i < num_fut_dates ; i++)
		{	
			err = srt_f_lgmcashfuture(	fut_dates[i], end_futdates[i],
										float_basis, und_name, 
										&fut_price[i]);
			if(err) goto e1;

			spread = swp_f_spread( fut_dates[i], end_futdates[i], ref_rate_code );
			fut_price[i] -= spread;

			//Result Type is conv adj
			if(result_type !=0)
			{
				
				//cash LIBOR
				df_start = swp_f_df(today, fut_dates[i], yc_name);
				df_end = swp_f_df(today, end_futdates[i], yc_name);
				fwd = (df_start / df_end - 1.0) / coverage(fut_dates[i], end_futdates[i], float_basis);
				fwd += spread;

				fut_price[i] = 1.0 - fut_price[i] - fwd;
			}
		}

		//Free memory
e1:	
	if (end_futdates) free(end_futdates);
	return err;

}// END Err srt_f_price_future_lgm(...)






//******************************************************************************
/* Function to determine the price of a future Given Calibration Instruments */
//******************************************************************************
Err srt_f_autocal_future_input_instr (
			long              *fut_dates, 
			long              num_fut_dates,
			long			  ninstruments,
			char			  **type_str,
			long			  *start_dates,
			long			  *end_dates,
			char			  **compd_str,
			char 			  **basis_str,
			char			  **recpay_str,
			char			  **refratecode,
			char              *mdl_name, 
			char              *yc_name,
			char			  *volcurvename,
			char              *ref_rate_code,
			double            tau, 
			double            alpha, 
			double            gamma, 
			double            rho,
			double            *fut_price,
			double			  *instr_price,
			int				  bumpinst)
{
Err			       err = NULL;
int				   i;
Date *			   end_futdates = NULL;
SrtBasisCode *	   basis = NULL;
SrtCompounding *   compounding = NULL;
double             spread;
double			   lognorm = 1.0;
SrtCompounding     float_compounding;
SrtBasisCode       float_basis;

char			   und_name[256];

SrtMdlType		   mdl_typ;
SrtMdlDim		   mdl_dim;

long	            volCrvRows ;
long                volCrvCols ;
double            **volCrvVals;

long	            tauCrvRows ;
long                tauCrvCols ;
double            **tauCrvVals;
long			   spotlag,today;
long				fixingdate;
double				fwd,vol,mat;
char **			   grfnparams_str = NULL;
char **			   grfnvalue_str = NULL;
double*			   price = NULL;
long			   numgrfnparams;
long			  numpaths = 100;
double			*strike = NULL;
double			*bondstrike = NULL;
SrtCrvPtr		  yldcrv;

yldcrv = lookup_curve(yc_name);
spotlag = get_spotlag_from_curve(yldcrv);
today = get_today_from_curve(yldcrv);


	/* Memory allocation (useful for calibration */
	end_futdates   = (Date *) calloc(num_fut_dates, sizeof(Date));
	basis       = (SrtBasisCode *) calloc (num_fut_dates, sizeof(SrtBasisCode));
	compounding = (SrtCompounding *) calloc (num_fut_dates, sizeof(SrtCompounding));
	price       = dvector(0, ninstruments-1);
	strike		= dvector(0,ninstruments-1);
	bondstrike	= dvector(0,ninstruments-1);

	if ((!end_dates) )
	{
		err = "Memory Allocation Failure";
		goto e1;
	}

/* Get the reference rate details (copmpounding and basis) */
	err = swp_f_get_ref_rate_details(ref_rate_code, &float_basis, &float_compounding);
	if (err)
		goto e1;

/* Fill arrays needed for the calibration (each instrument is a single caplet) */
	for (i = 0 ; i < num_fut_dates ; i++)
	{
	/* End dates of each caplet */
		end_futdates[i] = add_unit(fut_dates[i], (int)(12/float_compounding), 
			SRT_MONTH, MODIFIED_SUCCEEDING);
		
	/* Compounding of each caplet */
		compounding[i] = float_compounding;
		

	/* Basis of each caplet */
		basis[i] = float_basis;
		



	}/* END loop on instruments needed for calibration */


/* -------------------------------------------------------
	Get The relevant information for the Calibration
   ------------------------------------------------------- */

	for (i=0;i<ninstruments;i++)
	{
		if (strcmp(type_str[i],"SWAPTION")==0)
		{
			err = swp_f_ForwardRate(start_dates[i],
									end_dates[i],
									compd_str[i],
									basis_str[i],
									yc_name,
									refratecode[i],
									&fwd);
			
			err = swp_f_vol(volcurvename,start_dates[i],end_dates[i],fwd,&vol,&lognorm); 

			if (i+1 == bumpinst) vol+=0.01; 

			err = swp_f_Swaption(start_dates[i], end_dates[i],compd_str[i], basis_str[i], vol, fwd, recpay_str[i], refratecode[i],yc_name, "PREMIUM", "LOGNORMAL",&(price[i]));

		}
		else
		{
			err = srt_BGM_FwdRate(start_dates[i], end_dates[i], refratecode[i], yc_name, 1, &fwd);
			
			err = swp_f_vol(volcurvename,start_dates[i],end_dates[i],fwd,&vol,&lognorm); 

			if (i+1 == bumpinst) vol+=0.01; 

			fixingdate = add_unit(start_dates[i],-spotlag,SRT_BDAY, MODIFIED_SUCCEEDING);
			
			mat = coverage(today, fixingdate, BASIS_ACT_365);

			if(strcmp(recpay_str[i],"PAY") == 0)
			{
				price[i] = coverage(start_dates[i],end_dates[i],basis[i])*swp_f_df(today, end_dates[i], yc_name)*srt_f_optblksch(fwd, fwd, vol, mat, 1, SRT_CALL, SRT_PREMIUM);
			}
			else
			{
				price[i] = coverage(start_dates[i],end_dates[i],basis[i])*swp_f_df(today, end_dates[i], yc_name)*srt_f_optblksch(fwd, fwd, vol, mat, 1, SRT_PUT, SRT_PREMIUM);
			}

			instr_price[i]=price[i];

		}

		strike[i] = fwd;
		bondstrike[i] = 1.0;
	}




/* -------------------------------------------------------
   Create the associated interest rate underlying
   ------------------------------------------------------- */

/* Get the model type and model dimension */
	if (err = srt_f_interp_model(mdl_name, &mdl_typ, &mdl_dim))
	{
		goto e1;
	}

/* Underlying name */
	strcpy(und_name, "Future_ir_und");

/* Fill volatility */
	volCrvRows = 1;
	volCrvCols = 1;
	volCrvVals = dmatrix(0, volCrvCols-1, 0, volCrvRows-1);
	if (!volCrvVals)
	{
		err = "Memory Allocation Failure";
		goto e1;
	}

	if (mdl_typ == LGM)
	{
		volCrvVals[0][0] =  0.01;
	}	
	else if (mdl_typ == CHEY)
	{
		volCrvVals[0][0] =  0.20;
	}

/* Fill Tau */
	tauCrvRows = 1;
	tauCrvCols = 1;
	tauCrvVals = dmatrix(0, tauCrvCols-1, 0, tauCrvRows-1);
	if (!tauCrvVals)
	{
		err = "Memory Allocation Failure";
		goto e1;
	}
	tauCrvVals[0][0] =  8.0;
	
/* Initialise underlying */
	err = SrtInitIRUnd(	und_name, yc_name, mdl_name,
						volCrvRows, volCrvCols, volCrvVals,
						tauCrvRows, tauCrvCols, tauCrvVals,
						0.0, alpha, gamma, rho, 0.0, 0.0,
						0.0,0,0,NULL); 
	
	
	if (err)
	{
		goto e1;
	}


/* -------------------------------------
   The calibration itself
   ------------------------------------- */
/* Model parameters */
	numgrfnparams = 5;
	grfnparams_str = (char **) svector_size(0, numgrfnparams-1, 32);
	grfnvalue_str = (char **) svector_size(0, numgrfnparams-1, 32);
	if ((!grfnparams_str) || (!grfnvalue_str))
	{
		err = "Memory Allocation Failure";
		goto e1;
	}

	strcpy(grfnparams_str[0], "FORCEMC");
	if /*(*/(mdl_typ == CHEY)/* && (mdl_dim == TWO_FAC))*/
		strcpy(grfnvalue_str[0], "YES");
	else
		strcpy(grfnvalue_str[0], "NO");

	strcpy(grfnparams_str[1], "SAMPLETYPE");
	strcpy(grfnvalue_str[1], "BALANTISAM");

	strcpy(grfnparams_str[2], "MINNODE");
	strcpy(grfnvalue_str[2], "1");

	strcpy(grfnparams_str[3], "MAXTIME");
	strcpy(grfnvalue_str[3], "0.0384615");

	strcpy(grfnparams_str[4], "NUMPATH");
	sprintf(grfnvalue_str[4], "%ld", numpaths);

/* Call to the calibration routine: bootstrap */
	err = SrtBootstrap(numgrfnparams, grfnparams_str, grfnvalue_str, (int*) &ninstruments,
                   start_dates, end_dates, compd_str, basis_str,
		           strike, bondstrike, type_str, recpay_str,
			       refratecode,
					price,  price, tau, 
			       alpha, gamma, rho,
				   und_name);
	if (err)
	{
		goto e1;
	}

/* -------------------------------------------------------------------------------------
   The Pricing of futures according to the model (closed form for LGM, tableau for Chey)
   ------------------------------------------------------------------------------------- */
	if (mdl_typ == LGM)
	{
		for (i = 0 ; i < num_fut_dates ; i++)
		{	
			if (err = srt_f_lgmcashfuture(	fut_dates[i], end_futdates[i],
										basis[i], und_name, 
										&fut_price[i]))
			{
				goto e1;
			}
			spread = swp_f_spread( fut_dates[i], end_futdates[i], ref_rate_code );
			fut_price[i] -= spread;
		}
	}
	else 
	if (mdl_typ == CHEY)
	{
		for (i = 0 ; i < num_fut_dates ; i++)
		{	
			if (err = srt_f_cheycashfuture(	fut_dates[i], end_futdates[i],
										"MM", und_name,
										numpaths,
										&fut_price[i]))
			{
				goto e1;
			}
			spread = swp_f_spread( fut_dates[i], end_futdates[i], ref_rate_code );
			fut_price[i] -= spread;
		}
	}

/* Destroy the underlying */
/*	err = srt_f_destroy_und(und_name);*/

/* Free memory */	

e1:	if (grfnparams_str) free_svector_size(grfnparams_str, 0, numgrfnparams-1, 32);
	if (grfnvalue_str) free_svector_size(grfnvalue_str, 0, numgrfnparams-1, 32);
	if (end_futdates) free(end_futdates);
	if (basis) free(basis);
	if (price) free_dvector(price, 0, ninstruments-1);
	if (strike) free_dvector(strike,0,ninstruments-1);
	if (bondstrike)	free_dvector(bondstrike,0,ninstruments-1);
	
	return err;

}/* END Err srt_f_autocal_future(...) */ 

/* --------------------------------------------------------------------------- */
/* Function to do an autocalibration to determine the price of a future */
Err srt_f_autocal_future (
			long              *fut_dates, 
			long              num_fut_dates,
			Err			(*GetVol)(double dStart, double dEnd, double dStrike, 
						double dForward, double dSpread, double *pdBsVol),
							  /* function that returns the volatility of the reference swptns */
		    String            logNormStr,
			char              *mdl_name, 
			char              *yc_name,
			char              *ref_rate_code,
			double            tau, 
			double            alpha, 
			double            gamma, 
			double            rho,
			long              numpaths,
			double            *fut_price)
{
Err			       err = NULL;
int				   i;
Date *			   end_dates = NULL;
char **			   compd_str = NULL;
char *			   basis_temp = NULL;
char *			   compounding_temp = NULL;
char **			   basis_str = NULL;
SrtBasisCode *	   basis = NULL;
SrtCompounding *   compounding = NULL;
double *		   strike = NULL;
double *		   bondstrike = NULL;
char **			   type_str = NULL;
char **			   recpay_str = NULL;
double *		   price = NULL;
char **            refratecode = NULL;
double             spread;

SrtCurvePtr		   yldcrv = NULL;
double			   today;
char 			   ccy_str[4]="";
SrtCompounding     float_compounding;
SrtBasisCode       float_basis;

char			   und_name[256];

SrtMdlType		   mdl_typ;
SrtMdlDim		   mdl_dim;

long	            volCrvRows ;
long                volCrvCols ;
double            **volCrvVals;

long	            tauCrvRows ;
long                tauCrvCols ;
double            **tauCrvVals;

char **			   grfnparams_str = NULL;
char **			   grfnvalue_str = NULL;
long			   numgrfnparams;

/* Get the curve */
	yldcrv = lookup_curve(yc_name);
	if (!yldcrv)
	{
		return serror("Can't find market");
	}
	
/* Get today */
	today	= get_clcndate_from_curve(yldcrv);
	

/* -------------------------------------------------------------
   Create the caplets associated with each future date 
   ------------------------------------------------------------- */
	/* Memory allocation (useful for calibration */
	end_dates   = (Date *) calloc(num_fut_dates, sizeof(Date));
	compd_str   = (char **) svector_size(0, num_fut_dates-1, 3);
 	basis_str   = (char **) svector_size(0, num_fut_dates-1, 8);
	basis       = (SrtBasisCode *) calloc (num_fut_dates, sizeof(SrtBasisCode));
	compounding = (SrtCompounding *) calloc (num_fut_dates, sizeof(SrtCompounding));
	strike      = dvector(0, num_fut_dates-1);
	bondstrike  = dvector(0, num_fut_dates-1);
	type_str    = (char **) svector_size(0, num_fut_dates-1, 9);
	recpay_str  = (char **) svector_size(0, num_fut_dates-1, 4);
	price       = dvector(0, num_fut_dates-1);
	refratecode = (char **) svector_size(0, num_fut_dates-1, 16);
	
	if ((!end_dates) || (!compd_str) || (!basis_str)
		 || (!basis) || (!compounding) || (!strike) || (!bondstrike)
		|| (!type_str) || (!recpay_str) || (!price) ||  (!refratecode) )
	{
		err = "Memory Allocation Failure";
		goto e1;
	}

/* Get the reference rate details (compounding and basis) */
	err = swp_f_get_ref_rate_details(ref_rate_code, &float_basis, &float_compounding);
	if (err)
		goto e1;

/* Fill arrays needed for the calibration (each instrument is a single caplet) */
	for (i = 0 ; i < num_fut_dates ; i++)
	{
	/* End dates of each caplet */
		end_dates[i] = add_unit(fut_dates[i], (int)(12/float_compounding), 
			SRT_MONTH, MODIFIED_SUCCEEDING);
		
	/* Compounding of each caplet */
		strcpy(compd_str[i], "Q");
		err = translate_compounding(&compounding_temp, compounding[i]);
		if (err)
			goto e1;

	/* Basis of each caplet */
		basis[i] = float_basis;
		err = translate_basis(&basis_temp, basis[i]);
		if (err)
			goto e1;
		strncpy(basis_str[i], basis_temp, 7);
		basis_str[i][7] = '\0';

	/* The strike of each caplet is set at the money for the reference rate */
		strike[i] = swp_f_fra((Ddate)fut_dates[i], (Ddate)end_dates[i], basis[i], 
						yc_name, ref_rate_code); 

	/* The bond strike is set to one (caplets) */
		bondstrike[i] = 1.0;

	/* Set the option type: CAPFLOOR */
		strcpy(type_str[i], "CAPFLOOR");

	/* Set the receiver or payer type: cap */
		strcpy(recpay_str[i], "CAP");

	/* Set the reference rate code for each caplet */
		strcpy(refratecode[i], ref_rate_code);

	/* Get the PRICE of each caplet */
		err = swp_f_CapFloor( fut_dates[i], end_dates[i], strike[i], GetVol, 
				recpay_str[i], refratecode[i], yc_name, "PREMIUM", logNormStr,
				&price[i]);
		if (err)
			goto e1;		
	

	}/* END loop on instruments needed for calibration */


/* -------------------------------------------------------
   Create the associated interest rate underlying
   ------------------------------------------------------- */

/* Get the model type and model dimension */
	if (err = srt_f_interp_model(mdl_name, &mdl_typ, &mdl_dim))
	{
		goto e1;
	}

/* Underlying name */
	strcpy(und_name, "Future_ir_und");

/* Fill volatility */
	volCrvRows = 1;
	volCrvCols = 1;
	volCrvVals = dmatrix(0, volCrvCols-1, 0, volCrvRows-1);
	if (!volCrvVals)
	{
		err = "Memory Allocation Failure";
		goto e1;
	}

	if (mdl_typ == LGM)
	{
		volCrvVals[0][0] =  0.01;
	}	
	else if (mdl_typ == CHEY)
	{
		volCrvVals[0][0] =  0.20;
	}

/* Fill Tau */
	tauCrvRows = 1;
	tauCrvCols = 1;
	tauCrvVals = dmatrix(0, tauCrvCols-1, 0, tauCrvRows-1);
	if (!tauCrvVals)
	{
		err = "Memory Allocation Failure";
		goto e1;
	}
	tauCrvVals[0][0] =  8.0;
	
/* Initialise underlying */
	err = SrtInitIRUnd(	und_name, yc_name, mdl_name,
						volCrvRows, volCrvCols, volCrvVals,
						tauCrvRows, tauCrvCols, tauCrvVals,
						0.0, alpha, gamma, rho, 0.0, 0.0,
						0.0,0,0,NULL); /* vasicek params */
	
	
	if (err)
	{
		goto e1;
	}


/* -------------------------------------
   The calibration itself
   ------------------------------------- */
/* Model parameters */
	numgrfnparams = 5;
	grfnparams_str = (char **) svector_size(0, numgrfnparams-1, 32);
	grfnvalue_str = (char **) svector_size(0, numgrfnparams-1, 32);
	if ((!grfnparams_str) || (!grfnvalue_str))
	{
		err = "Memory Allocation Failure";
		goto e1;
	}

	strcpy(grfnparams_str[0], "FORCEMC");
	if /*(*/(mdl_typ == CHEY)/* && (mdl_dim == TWO_FAC))*/
		strcpy(grfnvalue_str[0], "YES");
	else
		strcpy(grfnvalue_str[0], "NO");

	strcpy(grfnparams_str[1], "SAMPLETYPE");
	strcpy(grfnvalue_str[1], "BALANTISAM");

	strcpy(grfnparams_str[2], "MINNODE");
	strcpy(grfnvalue_str[2], "1");

	strcpy(grfnparams_str[3], "MAXTIME");
	strcpy(grfnvalue_str[3], "0.0384615");

	strcpy(grfnparams_str[4], "NUMPATH");
	sprintf(grfnvalue_str[4], "%ld", numpaths);

/* Call to the calibration routine: bootstrap */
	err = SrtBootstrap(numgrfnparams, grfnparams_str, grfnvalue_str, (int *)&num_fut_dates,
                   fut_dates, end_dates, compd_str, basis_str,
		           strike, bondstrike, type_str, recpay_str,
			       refratecode,
					price,  price, tau, 
			       alpha, gamma, rho,
				   und_name);
	if (err)
	{
		goto e1;
	}

/* -------------------------------------------------------------------------------------
   The Pricing of futures according to the model (closed form for LGM, tableau for Chey)
   ------------------------------------------------------------------------------------- */
	if (mdl_typ == LGM)
	{
		for (i = 0 ; i < num_fut_dates ; i++)
		{	
			if (err = srt_f_lgmcashfuture(	fut_dates[i], end_dates[i],
										basis[i], und_name, 
										&fut_price[i]))
			{
				goto e1;
			}
			spread = swp_f_spread( fut_dates[i], end_dates[i], ref_rate_code );
			fut_price[i] -= spread;
		}
	}
	else 
	if (mdl_typ == CHEY)
	{
		for (i = 0 ; i < num_fut_dates ; i++)
		{	
			if (err = srt_f_cheycashfuture(	fut_dates[i], end_dates[i],
										basis_str[i], und_name,
										numpaths,
										&fut_price[i]))
			{
				goto e1;
			}
			spread = swp_f_spread( fut_dates[i], end_dates[i], ref_rate_code );
			fut_price[i] -= spread;
		}
	}

/* Destroy the underlying */
	err = srt_f_destroy_und(und_name);

/* Free memory */	

e1:	if (grfnparams_str) free_svector_size(grfnparams_str, 0, numgrfnparams-1, 32);
	if (grfnvalue_str) free_svector_size(grfnvalue_str, 0, numgrfnparams-1, 32);
	if (end_dates) free(end_dates);
	if (compd_str) free_svector_size(compd_str, 0, num_fut_dates-1, 3);
	if (basis_str) free_svector_size(basis_str, 0, num_fut_dates-1, 8);
	if (basis) free(basis);
	if (strike) free_dvector(strike, 0, num_fut_dates-1);
	if (bondstrike) free_dvector(bondstrike, 0, num_fut_dates-1);
	if (type_str) free_svector_size(type_str, 0, num_fut_dates-1, 9);
	if (recpay_str) free_svector_size(recpay_str, 0, num_fut_dates-1, 4);
	if (price) free_dvector(price, 0, num_fut_dates-1);
	
	return err;

}/* END Err srt_f_autocal_future(...) */ 


/* --------------------------------------------------------------------------- */

/* Function to calculate the convexity biais between future and forward */
/* in a LGM model */
Err srt_f_lgmcashfuture (
			Ddate          start_date, 
			Ddate          end_date,
			SrtBasisCode   basis, 
			char          *und_name,
			double        *future_price)
{
SrtErr        err;
SrtUndPtr     und;
TermStruct	  *ts;
SrtMdlType    mdltype;
SrtMdlDim     mdldim;
String        yc_name;
double        G, H;
SrtTFTSMat    Two_fac_H;
double        cvg, disc_fac;
double        fixing_time;
Ddate         fixing_date, today; 
Ddate         start_end[2];
double        phi[2][2];
YTatt_param   YT[2];
double        zc_exp, zc_var, expectation;
int           spot_lag;

	/* Check if it is LGM and get the dimension */
	und = lookup_und(und_name);
	if (!ISUNDTYPE(und,INTEREST_RATE_UND))
		return serror ("%s should be an interest rate underlying",und_name); 

	err = get_underlying_mdltype(und, &mdltype);
	if (err)
		return (err);

	if (mdltype != LGM)
		return serror("Can compute future only for LGM");

	err = get_underlying_mdldim(und, &mdldim);
	if (err)
		return (err);

	today = get_today_from_underlying(und);
	spot_lag = get_spotlag_from_underlying(und);
	
	/* Compute the future only if it is available */
	if (start_date < today)
	{
		*future_price = 0.0;
		return NULL;
	}

	/* Get the expectation of 1/B(T1, T2) under measure beta */
	fixing_date		= (Ddate) add_unit((Date)start_date, -spot_lag, SRT_BDAY, SUCCEEDING);
	fixing_time		= (double) ( fixing_date - today) * YEARS_IN_DAY;

	/* get reconstruction coefficients */
	start_end[0] = start_date;
	start_end[1] = end_date;

	err = Y_T_at_t_param (fixing_date, start_end, 2, und, YT);
	if (err)
	{
		return serror (err);
	}

	err = get_underlying_ts (und, &ts);
	if (err)
		return (err);

	err = srt_f_lgm_phifunc	(und, (Ddate) fixing_date, phi);
	if (err)
	{
		return serror (err);
	}

	if (mdldim == ONE_FAC)
	{
		G_H_func (fixing_time, ts, &G, &H);
		
		zc_exp =	(YT[1].phi_coeff[0][0] - YT[0].phi_coeff[0][0])
					* phi[0][0]
					+ (YT[1].x_coeff[0] - YT[0].x_coeff[0])
					* H;

		zc_var =	phi[0][0]
					* (YT[1].x_coeff[0] - YT[0].x_coeff[0])		  		
					* (YT[1].x_coeff[0] - YT[0].x_coeff[0]); 
	}
	else if (mdldim == TWO_FAC)
	{
		err = get_2f_H_funcs (fixing_time, ts, &Two_fac_H);

		zc_exp =   (YT[1].phi_coeff[0][0] - YT[0].phi_coeff[0][0])
				   * phi[0][0]
				   + (YT[1].phi_coeff[1][1] - YT[0].phi_coeff[1][1])
				   * phi[1][1]
				   + (YT[1].phi_coeff[0][1] - YT[0].phi_coeff[0][1])
				   * phi[0][1]
				   + (YT[1].phi_coeff[1][0] - YT[0].phi_coeff[1][0])
				   * phi[1][0]
				   + (YT[1].x_coeff[0] - YT[0].x_coeff[0])
				   * (Two_fac_H[0][0] + Two_fac_H[0][1])
				   + (YT[1].x_coeff[1] - YT[0].x_coeff[1])
				   * (Two_fac_H[1][1] + Two_fac_H[1][0]);
		
		zc_var =   phi[0][0] 
				   * (YT[1].x_coeff[0] - YT[0].x_coeff[0])		  
				   * (YT[1].x_coeff[0] - YT[0].x_coeff[0])
				   + phi[1][1] 
				   * (YT[1].x_coeff[1] - YT[0].x_coeff[1])		  
				   * (YT[1].x_coeff[1] - YT[0].x_coeff[1])
				   + phi[0][1] 
				   * (YT[1].x_coeff[0] - YT[0].x_coeff[0])		  
				   * (YT[1].x_coeff[1] - YT[0].x_coeff[1])
				   + phi[1][0] 
				   * (YT[1].x_coeff[1] - YT[0].x_coeff[1])		  
				   * (YT[1].x_coeff[0] - YT[0].x_coeff[0]);
	}

	expectation = exp(zc_exp + zc_var / 2);
	
	/* coverage */
	cvg = coverage((Date)start_date, (Date)end_date, basis);

	/* Get the yield curve associated with the underlying */
	yc_name = get_ycname_from_irund(und);

	/* discount factor */
	disc_fac = swp_f_df(start_date, end_date, yc_name); 

	*future_price = 1.0 - (expectation / disc_fac - 1.0) / cvg;

	return err;
}

/* ----------------------------------------------------------------------------- */

/* Function to calculate the convexity biais between future and forward */
/* in a CHEYETTE model */
/* To do it we build a GRFN tableau */
SrtErr srt_f_cheycashfuture (Ddate start_date, Ddate end_date,
						char * basis, char * und_name,
						long numpaths,
						double * future_price)
{
int			i;
SrtErr		err = NULL;
SrtUndPtr   und;
SrtMdlType	mdltype;
SrtMdlDim	mdldim;

Ddate		fixing_date;
int			numeventdates;
long *		eventdates = NULL;

int			tableaurows;
int			tableaucols;
char ***	tableaustrings = NULL;
int **		tableaumask = NULL;

char **		grfnparams_str = NULL;
char **		grfnvalue_str = NULL;
long		numgrfnparams;

const char * histo_name = "HISTO";
double *	histo_x = NULL;
double *	histo_y = NULL;
long		x_size, y_size;
double **	histo_values = NULL;

int         spot_lag;
double		price, stdev;	

	/* Check if it is CHEY and get the dimension */
	und = lookup_und(und_name);
	if (!ISUNDTYPE(und,INTEREST_RATE_UND))
		return serror ("%s should be an interest rate underlying",und_name); 

	err = get_underlying_mdltype(und, &mdltype);
	if (err)
		return (err);

	if (mdltype != CHEY)
		return serror("Can compute future only for CHEY");

	err = get_underlying_mdldim(und, &mdldim);
	if (err)
		return (err);
	spot_lag = get_spotlag_from_underlying(und);

	/**************************************************************/
	/* To price a future in a cheyette model we do a GRFN tableau */
	/**************************************************************/
	/* model parameters */
	numgrfnparams = 5;
	grfnparams_str = (char **) svector_size(0, numgrfnparams-1, 32);
	grfnvalue_str = (char **) svector_size(0, numgrfnparams-1, 32);
	if ((!grfnparams_str) || (!grfnvalue_str))
	{
		err = "Memory Allocation Failure";
		goto e1;
	}

	strcpy(grfnparams_str[0], "FORCEMC");
	strcpy(grfnvalue_str[0], "YES");

	strcpy(grfnparams_str[1], "SAMPLETYPE");
	strcpy(grfnvalue_str[1], "BALANTISAM");

	strcpy(grfnparams_str[2], "MINNODE");
	strcpy(grfnvalue_str[2], "1");

	strcpy(grfnparams_str[3], "MAXTIME");
	strcpy(grfnvalue_str[3], "0.0384615");

	strcpy(grfnparams_str[4], "NUMPATH");
	sprintf(grfnvalue_str[4], "%ld", numpaths);

	/* Dates */
	numeventdates = 1;
	eventdates = lngvector(0,numeventdates-1);
	if (!eventdates)
	{
		err = "Memory Allocation Failure";
		goto e2;
	}

	fixing_date		= (Ddate) add_unit((Date)start_date, -spot_lag, SRT_BDAY, SUCCEEDING);
	eventdates[0]	= (Date) fixing_date;

	/* Tableau */
	tableaurows = 1;
	tableaucols = 1;
	tableaustrings = smatrix_size(0,0,0,0,GRFN_DEF_ARGBUFSZ);
	tableaumask = imatrix(0,0,0,0);
	if ((!tableaustrings) || (!tableaumask))
	{
		err = "Memory Allocation Failure";
		goto e3;
	}

	sprintf(tableaustrings[0][0], "FRA(%.0f,%.0f,\"%s\")|hist(\"%s\")",
									start_date, end_date, basis, histo_name); 
	tableaumask[0][0] = GRFNSCELL;

	/* Call to GRFN */
	err = SrtGrfnMain(	und_name,
						numgrfnparams, grfnparams_str, grfnvalue_str,
						numeventdates, eventdates,
						tableaurows, tableaucols, tableaustrings, tableaumask,
						0, NULL, NULL, 
						&price, &stdev, NULL, NULL);

	if (err)
    {
		goto e3;
	}


/* Get the result from the histogram */
	err = srt_f_makehistogram(	(char *)histo_name,0,"",0,
								&histo_x, &x_size,&histo_y, &y_size, &histo_values);
	if (err)
		goto e4;

/* Do the mean from the histogram */
	*future_price = 0;
	for (i = 0 ;  i < x_size ; i++)
	{
		*future_price += histo_values[i][0];
	}

	*future_price /= x_size;
	*future_price = 1.0 - (*future_price);

e4:	if (histo_x) free_dvector(histo_x, 0, x_size - 1);
	if (histo_y) free_dvector(histo_y, 0, y_size - 1);
	if (histo_values) free_dmatrix (histo_values, 0, x_size - 1, 0, y_size - 1);


e3:	if (tableaustrings) free_smatrix_size(tableaustrings,0,0,0,0,GRFN_DEF_ARGBUFSZ);
	if (tableaumask) free_imatrix(tableaumask,0,0,0,0);

e2:	if (eventdates) free_lngvector(eventdates,0,numeventdates-1);

e1:	if (grfnparams_str) free_svector_size(grfnparams_str, 0, numgrfnparams-1, 32);
	if (grfnvalue_str) free_svector_size(grfnvalue_str,0, numgrfnparams-1, 32);

	return err;
}