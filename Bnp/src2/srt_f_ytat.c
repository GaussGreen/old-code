/* ====================================================================================
   FILENAME:  srt_f_ytat.c

   PURPOSE:   The parameters needed at a date t to compute a DF maturity T
              in any model (interest rate OR NOT)
   ==================================================================================== */

#include "math.h"
#include "srt_h_all.h"
#include "srt_h_ytat.h"
#include "srt_h_mdl_types.h"
#include "srt_h_powermodel.h"



/*----------------------------------------------------------------------------*/

/* ----------------------------------------------------------------------------

  FUNCNAME        : Y_T_at_t_param
  
DESCRIPTION     :   Y_T_at_t_param computes 3 components in eqn. A.18, 
   Cheyette(1991) which enter Y_T_at_T_compute which computes the actual 
   zero-coupon yield of maturity T_date at time t.

   NOTE: Since the actual data structure for "sample" will be a 2-dim array of
         structures, containing r, phi, Z, 
         one needs to pass one row of random numbers corresponding to one 
         sample path of the simulation to this routine, which will
         o/p one sample path (i.e. one vector) of the relevant state variables.
   
   variables: t=date acting as origin for YTM Y_T_at_t with timetomat T ;
              fwd_rate_t:= fwd rate for time t, today;
              Tmat:= pointer to array of relevant maturity dates
                     for which yields-to-maturity will be computed.
              dim_Tdate:= dimension of above array.  

            Major changes for Two Factor Model
		    note that now this coeff are computed to return 
		    (T-t) * Y(t,T) = - ln B(t,T) and not Y(t,T) anymore
---------------------------------------------------------------------------*/

Err Vasicek_Y_T_at_param(Ddate		fix_date, 
						Ddate		*pay_dates,
						long		dim_pay_date,
						SrtUndPtr	und,
						YTatt_param *YT_sam)
{
	Err			err = NULL;
	TermStruct  *ts;
	long		 i;
	double		 F_t,G_t,H_fix_t,Psi_t,S_t,O_t,
				 vasicek_mean_sr_t,vasicek_mean_int_sr_t;
	double		 vasicek_mean_int_sr_pay_T,G_pay_T,H_pay_T;
	double       fix_t,pay_T;
	String       yc_name;
	SrtCurvePtr  yldcrv;
	Date         today;

	/* get the yc name */
	yc_name  = get_discname_from_underlying(und) ;
	yldcrv = lookup_curve(yc_name) ;
	if (!yldcrv)
		return serror("Yield Curve %s is not defined",yc_name);

	/* extract today from the curve */
	today = get_clcndate_from_curve(yldcrv);
	fix_t = (double)(fix_date - today) * YEARS_IN_DAY;

	/* get the underlying term struct */
	err = get_underlying_ts(und, &ts) ;
	if (err) return err;

	/* compute the variable that are only t_date dependent */
	F_t = F_func(fix_t, ts);
	G_H_func(fix_t,ts,&G_t,&H_fix_t);
	Psi_t = Psi_func(fix_t, ts);
	S_t = S_func(fix_t,ts);
	O_t = O_func(fix_t,ts);

	err = srt_f_get_vasicek_mean_sr(fix_t,
									ts,
									&vasicek_mean_sr_t);
	if(err) return err;

	err = srt_f_get_vasicek_mean_int_sr(fix_t,
										ts,
										&vasicek_mean_int_sr_t);
	if(err) return err;

	/* Compute the part of the coefficient that is T dependent */
	for(i=0; i< dim_pay_date; i++)
	{
		if(pay_dates[i] > fix_date)
		{ 
			pay_T = (double)(pay_dates[i] - today) * YEARS_IN_DAY;
			G_H_func(pay_T,ts,&G_pay_T,&H_pay_T);

			YT_sam[i].vasicek_sr_coeff = ( Psi_func(pay_T,ts) - Psi_t ) / F_t;
			
			err = srt_f_get_vasicek_mean_int_sr(pay_T,
												ts,
												&vasicek_mean_int_sr_pay_T);
			if(err) return err;

			YT_sam[i].vasicek_mean_coeff = (vasicek_mean_int_sr_pay_T-vasicek_mean_int_sr_t);

			YT_sam[i].vasicek_mean_coeff -= (Psi_func(pay_T,ts) - Psi_t)*vasicek_mean_sr_t;

			YT_sam[i].vasicek_var_coeff = Psi_func(pay_T,ts)*Psi_func(pay_T,ts)*(G_pay_T-G_t);
			YT_sam[i].vasicek_var_coeff -= 2*Psi_func(pay_T,ts)*(S_func(pay_T,ts)-S_t);
			YT_sam[i].vasicek_var_coeff += (G_pay_T*Psi_func(pay_T,ts)*Psi_func(pay_T,ts)-2*O_func(pay_T,ts));
			YT_sam[i].vasicek_var_coeff -= (G_t*Psi_t*Psi_t-2*O_t);

		}
	}

	return err;

}

Err Y_T_at_t_param( 
		Ddate         t_date, /* Simulation date */
		Ddate        *T_date, /* Discount Factor expiration dates (array) */ 
		long          dim_Tdate, 
		SrtUndPtr     und,
		YTatt_param  *YT_sam )
{   
long          i;
double        a0,a1;
double        F_t,Psi_t;
SrtTFTSVec    twofac_F_t, twofac_Psi_t, twofac_Psi_T;
Date          today;
double        t_in_y;
double        T_in_y;
double        df;
TermStruct    *ts                 = NULL ;

String        yc_name ;
SrtMdlDim     mdl_dim;
SrtMdlType    mdl_type;
Err           err                 = NULL;
SrtCurvePtr   yldcrv;

double		  **M;
double		  s,theta;
double		  beta,eta;
double        lambda_t,lambda_T;
double        zeta_t,zeta_T;
double        A_t, A_T;
double        H_t, H_T;

	
/* Get the underlying model and type */
	err = get_underlying_mdldim(und, &mdl_dim);
	if (err)
		return err;
	err = get_underlying_mdltype(und, &mdl_type);
	if (err)
		return err;
	
/* Get the name of the discount curve associated to the underlying */
	yc_name  = get_discname_from_underlying(und) ;
	yldcrv = lookup_curve(yc_name) ;
	if (!yldcrv)
		return serror("Yield Curve %s is not defined",yc_name);

/* Extract today from the curve */
	today = get_clcndate_from_curve(yldcrv);
	t_in_y = (double)(t_date - today) * YEARS_IN_DAY;

/* Compute the forward zero-rates from the yield curve for all the zero-coupons */
	for(i=0; i< dim_Tdate; i++)
	{
	/* If the zero-coupn date is STRICTLY above today */
		if(T_date[i] > t_date)
		{ 
			T_in_y = (double)(T_date[i] - today) * YEARS_IN_DAY;

		/* Compute the log of forward zero-coupon: (T-t) * Y(0,t,T) */	
			df = swp_f_df(t_date, T_date[i], yc_name);
			if (df == SRT_DF_ERROR)
				return serror("Could not compute df(%f,%f) for %s", t_date, T_date[i], yc_name);  
			
			YT_sam[i].fwd_zero_rate = -log(df);
		}
	
		else
	
	/* The zero-coupon rate is set to 0.0 */
		YT_sam[i].fwd_zero_rate = 0.0;
	
	} /* END for(i=0; i< dim_Tdate; i++) */
	
/* If we have an interest rate model: more has to be done */
	if (ISUNDTYPE(und,INTEREST_RATE_UND))
	{
	/* Get the underlying Term Structure of volatility */
		err = get_underlying_ts(und, &ts) ;
		if (err)
			return err;
		
	/* Computes and stores what is only t (simulation date) dependent */
		if ( is_irund_Cheyette_type(und) )
		{
			switch (mdl_dim)
			{
				case ONE_FAC:
					F_t = F_func(t_in_y, ts);
					Psi_t = Psi_func(t_in_y, ts);
				break;

				case TWO_FAC:
					err = get_2f_F_funcs(t_in_y, ts, &twofac_F_t);	
					err = get_2f_Psi_funcs(t_in_y, ts, &twofac_Psi_t);	
					if (err)
						return err;
				break;
			}
		}

		else
	/* For the ETABETA model */
		if (mdl_type == ETABETA)
		{
			M		 = find_M_eta_beta(t_in_y,ts);	
			beta     = find_beta(t_in_y,ts);  
			eta      = find_eta(t_in_y,ts);  
		/* Lambda is the factor function (multiplied by X in the numeraire) */
			lambda_t = Psi_func(t_in_y,ts);
		/* Zeta is the integral of the square volatility */
			zeta_t   = Zeta_func(t_in_y,ts);
			
			/* A(t) = M(0,0,t) */
			s        = beta*beta*zeta_t;
			theta    = lambda_t/beta;
			A_t      = M_eta_beta_func(s,theta,eta,M); 
		
		}
		else
	/* NEW MODEL */
		if (is_model_New_type(mdl_type))
		{
			H_t = find_struct_interp( t_in_y, H, ts);
		}

	/* Compute the part of the coefficient that is T dependent */
		for(i=0; i< dim_Tdate; i++)
		{
			if(T_date[i] > t_date)
			{ 
				T_in_y = (double)(T_date[i] - today) * YEARS_IN_DAY;
				
			/* Compute the coefficients in front of the state variables  */	
				if ( is_irund_Cheyette_type(und) )
				{
					switch (mdl_dim)
					{
						case ONE_FAC:
							a0 = ( Psi_func(T_in_y,ts) - Psi_t ) / F_t;
							
							YT_sam[i].x_coeff[0]       = a0; 		 
							YT_sam[i].phi_coeff[0][0]  =  0.5 * a0 *a0 ; 
						break;
					
						case TWO_FAC:
							err = get_2f_Psi_funcs(T_in_y, ts, &twofac_Psi_T);	
							a0 = (twofac_Psi_T[0] - twofac_Psi_t[0]) / twofac_F_t[0] ;
							a1 = (twofac_Psi_T[1] - twofac_Psi_t[1]) / twofac_F_t[1] ;
							
							YT_sam[i].x_coeff[0]       = a0;	               
							YT_sam[i].phi_coeff[0][0]  = 0.5 * a0 * a0;  
							YT_sam[i].phi_coeff[0][1]  = 0.5 * a0 * a1;   
							YT_sam[i].x_coeff[1]       = a1;	
							YT_sam[i].phi_coeff[1][1]  = 0.5 * a1 * a1; 
							YT_sam[i].phi_coeff[1][0]  = 0.5 * a1 * a0;  
						break;
					}/* END switch statement */	
				}
				else
				
				if ( mdl_type == ETABETA )
				{
					lambda_T = Psi_func(T_in_y,ts);
			/*		G_H_func(T_in_y,ts,&lambda_T,&H_temp); */  
			/* looks like Lambda_func and G_H are same */
					zeta_T   = Zeta_func(T_in_y,ts);
				
					/* A(T) = M(0,0,T) */
					s = beta*beta*zeta_T;
					theta = lambda_T/beta;
					A_T = M_eta_beta_func(s,theta,eta,M);

					YT_sam[i].beta				  = beta;
					YT_sam[i].eta				  = eta;
					YT_sam[i].lambda_t			  = lambda_t;
					YT_sam[i].lambda_T			  = lambda_T;
					YT_sam[i].zeta_t			  = zeta_t;
					YT_sam[i].zeta_T			  = zeta_T;
					YT_sam[i].A_t				  = A_t;
					YT_sam[i].A_T                 = A_T;
					YT_sam[i].M		              = M;
				}
				else
				if (is_model_New_type(mdl_type))
				{
					H_T = find_struct_interp( T_in_y, H, ts);
					YT_sam[i].x_coeff[0] = - H_T + H_t;
					YT_sam[i].phi_coeff[0][0] = 0.5 * (H_T*H_T - H_t*H_t);
				}

				if (err)
					break;
			
			} /* END if T_date[i] > t */
			else
			{
				YT_sam[i].x_coeff[0] = YT_sam[i].x_coeff[1] = 0.0;
				YT_sam[i].phi_coeff[0][0] = YT_sam[i].phi_coeff[0][1] =0.0;
				YT_sam[i].phi_coeff[1][0] = YT_sam[i].phi_coeff[1][1] =0.0;

			/* Eta Beta */
				YT_sam[i].A_t= 0.0;
				YT_sam[i].A_T= 0.0;
				YT_sam[i].M= NULL;
				YT_sam[i].zeta_t= 0.0;
				YT_sam[i].zeta_T= 0.0;
				YT_sam[i].lambda_t= 0.0;
				YT_sam[i].lambda_T= 0.0;

			}
		
		}/* end for i loop on zero-coupon dates */
	
	} /* END if (ISUNTYPE(und,INTEREST_RATE) ) */	
	
/* Return a success message */
	return err;       

} /* END function Y_T_at_t_param */

Err Vasicek_Y_T_at_t_compute(long           dim_pay_date, 
							SrtSample		*sam, 
							YTatt_param		*YT_sam, 
							double			*YTt_yield, 
							int				index ,
							SrtMdlDim		mdl_dim,
							SrtMdlType		mdl_type )
{
	long i;

	for(i = 0; i < dim_pay_date; ++i)
	{
		YTt_yield[i] = YT_sam[i].vasicek_sr_coeff*samptr_get(sam,index,SHORT_RATE);
		YTt_yield[i] += YT_sam[i].vasicek_mean_coeff;
		YTt_yield[i] += 0.5*YT_sam[i].vasicek_var_coeff;
	}

	return NULL;
	
}

/*--------------------------------------------------------------------------

  FUNCNAME        :"Y_T_at_t_compute" 

  DESCRIPTION     :
   "Y_T_at_t_compute" calculates Y_T_on_t,i.e. the relevant zero-coupon 
   yield-to-maturities as defined by eqn. A.18, Cheyette(1991).

   NOTE: Since the actual data structure for "sample" will be a 2-dim array of
         structures, containing r, phi, Z, one needs to pass one row 
         of random numbers corresponding to one 
         sample path of the simulation to this routine, which will
         o/p one sample path (i.e. one vector) of the relevant state variables.
   
   variables: t=date acting as origin for YTM Y_T_at_t with timetomat T ;
              fwd_rate_t:= fwd rate for time t, today;
              T_date:= pointer to array of relevant maturity dates
                      for which yields-to-maturity will be computed.
              dim_Tdate:= dimension of T_date vector.  

  Description     : Computes (T-t)*Y(t,T) instead of previous Y(t,T)
		    major changes for two factor model...

   CAUTION: here, STATEVAR (resp X1,X2) should correspond to:
		 r -f(0,t) for the one factor model
		 ri  for the two factor model with r = f(o,t) + r1 + r2
-----------------------------------------------------------------------------*/

Err Y_T_at_t_compute( 
		long           dim_Tdate, 
		SrtSample     *sam, 
		YTatt_param   *YT_sam, 
		double        *YTt_yield, 
		int            index ,
		SrtMdlDim      mdl_dim,
		SrtMdlType     mdl_type )
{ 
long     i;
double	 M_stheta;
double   s,theta;		
Err      err = NULL;
                      
	for(i = 0; i < dim_Tdate; ++i)
	{
	/* For any underlying, there is at least the forward zero-coupon price */
		YTt_yield[i] =    YT_sam[i].fwd_zero_rate;
		
	/* If no sample has been passed, no simulation : continue */
		if (sam == NULL)
			continue;

	/* For Cheyette type underlying, there is more to do  */
		if (is_model_Cheyette_type(mdl_type) )
		{
			YTt_yield[i] += YT_sam[i].x_coeff[0] * samptr_get(sam,index,STATEVAR);
			YTt_yield[i] += YT_sam[i].phi_coeff[0][0]*samptr_get(sam,index,PHI) ;
			
			if (mdl_dim == TWO_FAC)
			{
				YTt_yield[i] += YT_sam[i].x_coeff[1] * samptr_get(sam,index,X2)  
					+ YT_sam[i].phi_coeff[1][1]*samptr_get(sam,index,PHI2)
				+ 2*YT_sam[i].phi_coeff[0][1]*samptr_get(sam,index,CROSSPHI) ;
			}
		}
		else
	/* For EtaBeta type underlying, there is more to do */
		if (mdl_type == ETABETA)
		{
		/* Compute s,theta for M(s,theta) reference */
			s = YT_sam[i].beta*YT_sam[i].beta*(YT_sam[i].zeta_T-YT_sam[i].zeta_t)/
				pow(1.+YT_sam[i].beta*samptr_get(sam,index,STATEVAR),2.*(1.-YT_sam[i].eta));
		
			if (YT_sam[i].beta == 0.)
				theta = 0.;
			else
				theta = YT_sam[i].lambda_T*(1.+YT_sam[i].beta*samptr_get(sam,index,STATEVAR))/YT_sam[i].beta;

			M_stheta = M_eta_beta_func(s,theta,YT_sam[i].eta,YT_sam[i].M);   
		

		/* Compute YTt_yield[i] from M_tT,A_t,A_T,lambda_t */
			YTt_yield[i] += (YT_sam[i].lambda_T - YT_sam[i].lambda_t)*samptr_get(sam,index,STATEVAR); 
			YTt_yield[i] += YT_sam[i].A_T - YT_sam[i].A_t;  /* A(T) - A(t) */
			YTt_yield[i] -= M_stheta;  

		}
		else
	/* NEW TYPE LGM MODEL more simple .... */
		if (is_model_New_type(mdl_type))
		{
			YTt_yield[i] += YT_sam[i].x_coeff[0] * samptr_get(sam,index,STATEVAR);
			YTt_yield[i] += YT_sam[i].phi_coeff[0][0] * samptr_get(sam,index,PHI);
		}
	
	} /* end for i loop on all zero-coupons */
	
/* Return a success message */
	return NULL;

} /* END function Y_T_at_t_compute */ 

/* -------------------------------------------------------------------------------- */


/* --------------------------------------------------------------------------------
  SOME UTILITIES FOR DEALING WITH TERM STRUCTURE OF VOL WHEN
  USING SrtRtFnc STRUCTURES
  E.Auld Jan 96
  -------------------------------------------------------------------------------- */
/*
 * set the Y_T_at_t_params needed to price a rate described
 * by a SrtRtFnc structure.  Will only allocate space if has not already
 * been allocated.
 */

SrtErr srt_f_rtinityp(SrtRtFnc rt, SrtUndPtr und, Ddate obsdate)
{
  YTatt_param *yp;
  long dim_Tdate;
  Ddate *T_date;
  SrtErr err;

/*
 * get dates
 */
  T_date = (Ddate*) srt_f_rtgetdt(rt);
  dim_Tdate = (long)srt_f_rtgetlen(rt);
  if(DTOL(obsdate) > DTOL(T_date[0]))
  {
    return serror("rate observed after it is set");
  }

/* 
 * get YTatt_params; they might not exist yet
 */
  yp = (YTatt_param *) srt_f_rtgetyp(rt);
  if(yp == NULL)
  {
    yp = (YTatt_param *)srt_calloc(dim_Tdate,sizeof(YTatt_param));
    srt_f_rtsetyp(rt,(void*)yp);
  }

/*
 * initialize them
 */

  err = Y_T_at_t_param(obsdate, T_date, dim_Tdate, und, yp);
  if(err) return err;
  return NULL;
}


/*
 * eval the Y_T_at_t_params contained in a SrtRtFnc structure;
 * and set the dfs.  If df field is blank, will allocate it.
 */

SrtErr srt_f_rtevalyp(
		SrtRtFnc      rt, 
		SrtSample     *sam, 
		int           index,  
		SrtMdlDim     mdl_dim, 
		SrtMdlType    mdl_type)
{
  int i;
  double *df;
  long dim_Tdate;
  YTatt_param *yp;
  SrtErr err = (SrtErr)NULL;

  dim_Tdate = (long)srt_f_rtgetlen(rt);
  yp = (YTatt_param *) srt_f_rtgetyp(rt);
  if(yp == NULL)
  {
    return serror("yield param not set!");
  }
  df = (double*)srt_f_rtgetdf(rt);
  if(df == NULL)
  {
    df = (double*)srt_calloc(dim_Tdate,sizeof(double));
    srt_f_rtsetdf(rt,(void*)df);
  }
  err = Y_T_at_t_compute(dim_Tdate, sam, yp, df, index , mdl_dim, mdl_type);
  for(i=0;i<dim_Tdate;i++)
  {
    df[i] = exp(-df[i]);
  }
  return err;
}
