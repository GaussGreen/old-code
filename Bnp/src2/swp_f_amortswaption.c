/* ===================================================================================
   FILENAME:      swp_f_amortswaption.c
   
   PURPOSE:       Computes price of a european option on an amortized swap
				in a co-initial Swap market model context
   =================================================================================== */

#pragma warning(disable:4786)	// Disable long name warnings

#include "math.h"
#include "num_h_allhdr.h"
#include "swp_h_all.h"
#include "swp_h_swaption.h"
#include "swp_h_swap_pricing.h"
#include "swp_h_swap_simple.h"
#include "swp_h_fwd_swaption.h"
#include "swp_h_amortswaption.h"
#include "AmortSwapUtils.h"
#include "AmortMidatCalib.h"
#include "swp_h_vol.h"
#include "opHeston.h"
#include "opfnctns.h"
#include "crtdbg.h"

Err ShiftedLogCopulaSimul(double **corr_mtx, long NPaths, long nFactorDim, int *pCutDim,double **res)
{
	int i, j,k;
	
	double *pdGauss = _alloca(sizeof(double)*nFactorDim);
	double **UnifSample=dmatrix(0,NPaths,0,nFactorDim-1);
	double **sqrt_corr_mtx=dmatrix(0,nFactorDim-1,0,nFactorDim-1);
	double **_corr_mtx = dmatrix(0,nFactorDim-1,0,nFactorDim-1);
	if( !(UnifSample&&sqrt_corr_mtx&&_corr_mtx))
		return "Memory allocation failed in AmortizedSwaptionShiftedLog";

	for(j=0;j<nFactorDim;++j)
		for(k=0;k<nFactorDim;++k)
			_corr_mtx[j][k]=corr_mtx[pCutDim[j]][pCutDim[k]];
	
	PositiveMatrix(_corr_mtx, nFactorDim);

	//////////  Step1: performs Chol decomposition of the correlation matrix corr_mtx
	CholDec	(&CholAlg, _corr_mtx, sqrt_corr_mtx, nFactorDim);

	//////////  Step2: Generates the matrix of r.v. 
	GetSobolMatrix(NPaths,nFactorDim,UnifSample);
		
	//////////  Step3: Generates the matrix of r.v. compatible with the Copula and the marginals
	for (i=0;i<NPaths;i++)
	{	
		/// fill in independent Gaussians
		const double *pdUnifSample = UnifSample[i];
		for (j=0; j<nFactorDim; j++)
		{
			//GaussSample[i][j][0]=inv_cumnorm_fast(UnifSample[i][j]);
			pdGauss[j]=inv_cumnorm_fast(pdUnifSample[j]);
		}

		// Correlate independent Caussians
		for (j=0; j<nFactorDim; j++)
		{
			//res[i][j]=0.0;
			res[i][j]=0.;
			for (k=0; k<nFactorDim;k++)
			{
				//res[i][j]+=sqrt_corr_mtx[j][k]*GaussSample[i][k][0];
				res[i][j] += sqrt_corr_mtx[j][k]*pdGauss[k];
			}
		}
	}

	free_dmatrix(UnifSample,0,NPaths-1,0,nFactorDim-1);
	free_dmatrix(sqrt_corr_mtx,0,nFactorDim-1,0,nFactorDim-1);
	free_dmatrix(_corr_mtx,0,nFactorDim-1,0,nFactorDim-1);
	return 0;
}

void ShiftedLogCopulaSimul_NotOpt(double **corr_mtx, long NPaths, long Dim, 
						   double *Fwd, double *Shifts, 
						   double *Vols, double **res)
{
	long i,idum=-8935807;
	int j,k;
	double	***GaussSample=NULL,**UnifSample=NULL,
			**sqrt_corr_mtx=NULL;
	double sign;
				
		//////////  Memory allocation
		GaussSample=dcube(0,NPaths,0, Dim-1,0,0);
		sqrt_corr_mtx=dmatrix(0,Dim-1,0,Dim-1);

		PositiveMatrix(corr_mtx, Dim);

		//////////  Step1: performs Chol decomposition of the correlation matrix corr_mtx
		CholDec	(&CholAlg,
				 corr_mtx,       
				 sqrt_corr_mtx,
				 Dim);

		//////////  Step2: Generates the matrix of r.v. 
		UnifSample=dmatrix(0,NPaths,0,Dim-1);
		GetSobolMatrix(NPaths,Dim,UnifSample);
			
		//////////  Step3: Generates the matrix of r.v. compatible with the Copula and the marginals
		for (i=0;i<NPaths;i++)
		{	
			for (j=0; j<Dim; j++)
			{
				GaussSample[i][j][0]=inv_cumnorm_fast(UnifSample[i][j]);
			}
			
			for (j=0; j<Dim; j++)
			{
				res[i][j]=0.0;
				for (k=0; k<Dim ;k++)
				{
					res[i][j]+=sqrt_corr_mtx[j][k]*GaussSample[i][k][0];
				}

				sign = 1;
				if(Vols[j] < 0)
				{
					sign = -1;
				}
//				res[i][j] = sign * ( (sign * Fwd[j] + Shifts[j]) * exp(-0.5 * Vols[j] * Vols[j] + Vols[j] * res[i][j]) - Shifts[j] );
				res[i][j] = ( (Fwd[j] + Shifts[j]) * exp(-0.5 * Vols[j] * Vols[j] + Vols[j] * res[i][j]) ) - Shifts[j];
			}
		}

	free_dmatrix(UnifSample,0,NPaths,0,Dim-1);
	free_dcube(GaussSample,0,NPaths,0,Dim-1,0,0);
	free_dmatrix(sqrt_corr_mtx,0,Dim-1,0,Dim-1);
}


//-------------------------------------------------------------------------------------
//---Function to compute the payoff of an amortized swap ------------------------------
//-------------------------------------------------------------------------------------

Err	GetAmortizedSwaptionPayoff_ForMAD(double  **dCopulaCube,
							 SrtCompounding srtFreq,
							 SrtBasisCode srtBasis,
			                 int     iNCoInitSwaps,
							 double  *dLvls,
							 double  *dStrikes,
							 double  *dNotionals,
							 int     iMC,
							 double  *dPayoff,
							 double	*controlvar,
							 int flag)
{
	Err err = NULL;
	int i;

	*dPayoff = 0;
	for(i=0;i<iNCoInitSwaps;++i)
	{
		*dPayoff += dNotionals[i] * dLvls[i] *(dCopulaCube[iMC][i]-dStrikes[i]) + controlvar[i];
	}

	return err;
}

//-------------------------------------------------------------------------------------
//---End Of Function to compute the payoff of an amortized swap ----------------------
//-------------------------------------------------------------------------------------



//--------------------------------------------------------------------------------------
//------------------Utilities Functions for Shifted Log Model---------------------------
//--------------------------------------------------------------------------------------

Err  GetShiftedLogDensity_Lvl_Measure_ForMAD(
						long today,
						char	*cYCname,
						char	*cVCname,
						char	*cRefRname,
						SrtCompounding srtFreq,
						SrtBasisCode srtBasis,
						long start,
						long end,
						double dStrike,
						int	   iCalibshift, //0 => Calibrate the Shift parameter, 
											//1 => Use the input Shift
						double *dShift,
						double *dVol
						)
{
	double *dVec = NULL, *dX=NULL,*dY=NULL, dThreshold=1.e-3;
	Err err = NULL;
	double epsilon = 1.0e-4;
	double dFwdSwap, dFwdCash, dSpread, dLvl, ex_time;
//	double sigma, alpha, beta, rho;
	char *cBasis, *cFreq;
	double power, atmvol;
	double shift, shifted_vol, equi_beta;
	double k1, k2, vol1, vol2, price1, price2;
	
	err = translate_basis(&cBasis, srtBasis);
	err = translate_compounding(&cFreq, srtFreq);

	err = swp_f_ForwardRate(start, end, cFreq, cBasis, cYCname, cRefRname, &dFwdSwap);
	dFwdCash = swp_f_swapcashrate( start, end, srtBasis, srtFreq, cYCname, cRefRname);
	err = swp_f_LevelPayment(start, end, cFreq, cBasis, cYCname, cRefRname, &dLvl);
	ex_time = (add_unit (start, - 2, SRT_BDAY, MODIFIED_SUCCEEDING) - today) * YEARS_IN_DAY;

	dSpread = dFwdSwap - dFwdCash;

/*	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&sigma, &power, SABR_BETAVOL );
	if (err) return err;
	if(sigma==0)
	{
		err = "Need a SABR Market";
		return err;
	}
	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&alpha, &power, SABR_ALPHA );
	if (err) return err;
	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&beta, &power, SABR_BETA );
	if (err) return err;
	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&rho, &power, SABR_RHO );
	if (err) return err;

	err = srt_f_optsarbvol(
				dFwdSwap,
				dFwdSwap, 
				ex_time, 
				sigma,  
				alpha,
				beta,
				rho,
				SRT_BETAVOL,
				SRT_LOGNORMAL,
				&atmvol);
*/

	err = swp_f_vol(cVCname, start, end, dFwdSwap, &atmvol, &power);
	if (err) return err;

	if(power<1)
	{
		atmvol = atmvol / dFwdSwap;
	}

	k1 = dFwdSwap;
	k2 = DMIN(dFwdSwap * (1 + 4*sqrt(ex_time)*atmvol), DMAX(dFwdSwap * exp(-4*sqrt(ex_time)*atmvol) , dStrike) );
	if( fabs(k2 - k1) < 0.001 )
	{
		k1 = k2 + 0.001;
	}
	
	err = swp_f_vol(cVCname, start, end, k1, &vol1, &power);
	err = swp_f_vol(cVCname, start, end, k2, &vol2, &power);
	if(err) return err;

	if(power<1)
	{
		price1 = srt_f_optblknrm(dFwdSwap, k1, vol1, ex_time, 1, SRT_CALL, SRT_PREMIUM);
		price2 = srt_f_optblknrm(dFwdSwap, k2, vol2, ex_time, 1, SRT_CALL, SRT_PREMIUM);

		err = srt_f_optimpvol(price1, dFwdSwap, k1, ex_time, 1, SRT_CALL, SRT_LOGNORMAL, &vol1);
		err = srt_f_optimpvol(price2, dFwdSwap, k2, ex_time, 1, SRT_CALL, SRT_LOGNORMAL, &vol2);
		if(err) return err;
	}

/*	err = srt_f_optsarbvol(
				dFwdSwap,
				k1, 
				ex_time, 
				sigma,  
				alpha,
				beta,
				rho,
				SRT_BETAVOL,
				SRT_LOGNORMAL,
				&vol1);
	
	err = srt_f_optsarbvol(
				dFwdSwap,
				k2, 
				ex_time, 
				sigma,  
				alpha,
				beta,
				rho,
				SRT_BETAVOL,
				SRT_LOGNORMAL,
				&vol2);
*/

	if(iCalibshift==0)
	{
		err = find_shifted_log_params(
					ex_time,
					dFwdSwap,
					k1,
					vol1,
					k2,
					vol2,
					&shift,
					&shifted_vol,
					&equi_beta);
		*dShift = shift;
	}
	else
	{
		price2 = srt_f_optblksch(dFwdSwap,
									k2,
									vol2,
									ex_time,
									1.0,
									SRT_CALL, 
									PREMIUM);

		err = srt_f_optimpvol(price2,
							*dShift+dFwdSwap, 
							*dShift+k2,
							ex_time, 
							1.0, 
							SRT_CALL,
							SRT_LOGNORMAL,
							&shifted_vol);
	}

	*dVol = shifted_vol * sqrt(ex_time);

	return err;
}


Err       Fwd_Swaption_Get_Cum_ShiftedLogDensity_ForMAD(Date lToday,
								 char *cYCname,
								 char *cVCname,
						         char *cRefRname,
								 SrtCompounding srtFreq,
								 SrtBasisCode srtBasis,

								 int iNDates,
								 Date *lStartDates,
								 Date *lEndDates,

								 double *dStrikes,

								 int	iCalibShift,
								 double Shift,
								 double *dShifts,
								 double *dVols)
{
	Err err = NULL;
	int i;
	double dMat,dCov;
	char *cBasis, *cFreq;
	double dDf;
	double shift;
	double vol;

	dMat=(lStartDates[0]-lToday)/365.; 
	
//////////////////////////////////////////
/// translate SORT types into strings
//////////////////////////////////////////

	err = translate_basis(&cBasis, srtBasis);
	err = translate_compounding(&cFreq, srtFreq);

	if (cFreq == "ANNUAL") {
		dCov=1.0;
	} else if (cFreq == "SEMIANNUAL") {
		dCov=0.5;
	} else if (cFreq == "QUARTERLY") {
		dCov=0.25;
	} else if (cFreq == "MONTHLY") {
		dCov=1./12.;
	}

//////////////////////////////////////////////////////////
////////////////// generates the cumulative densities
//////////////////////////////////////////////////////////	

	dDf=swp_f_df(lToday, lStartDates[0], cYCname);
	for (i=0;i<iNDates;i++)
	{
		shift = Shift;
		err = GetShiftedLogDensity_Lvl_Measure_ForMAD(lToday,
										cYCname,cVCname,cRefRname,
										srtFreq,srtBasis,
										lStartDates[0],
										lEndDates[i],
										dStrikes[i],
										iCalibShift,
										&shift,
										&vol);
		dShifts[i] = shift;
		dVols[i] = vol;
	}

	return err;

}	

//--------------------------------------------------------------------------------------
//----------------End Of Utilities Functions for Shifted Log Model----------------------
//--------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------
//---------------------Shifted Log Model : Function for MAD-----------------------------
//--------------------------------------------------------------------------------------
Err		AmortizedSwaptionShiftedLogForMAD(
						 char *cYCname,
						 char *cVCname,
						 char *cRefRname,

						 SrtCompounding srtFreq,
						 SrtBasisCode srtBasis,

						 double exer_fee,
						 
						 long	*lFixStartDates,
						 long	*lFixEndDates,
						 long	lNFixDates,
						 double *dFixNotionals,
						 double *dFixRates,

						 long	*lFloatStartDates,
						 long	*lFloatEndDates,
						 long	lNFloatDates,
						 double *dFloatNotionals,
						 double *dFloatMargins,
						 double *dFloatSpreads,

						 SrtCallPutType srtCallPut,

						 int xlNumSim,

						 int iCalibShift,
						 double Shift,
						 int UseVol,

						 double **Correl,
						 double *dPrice
			            )
{
	///// Objects and variables declaration
	Err err = NULL;
	SrtCurvePtr	pCurve = lookup_curve (cYCname);
	const Date lToday = (Date) get_today_from_curve (pCurve);
	const Date lSpotLag = (Date) get_spotlag_from_curve(pCurve);	
	double dIV, mean;

#ifdef _DEBUG
	double debug = 0.;
#endif //_DEBUG
	
	const int flag = (srtCallPut == SRT_CALL)?1:-1;
	int i, j;
  
	// fixed leg 
	double *dCoInitNotionals = _alloca(sizeof(double)*lNFixDates);
	double *dLvls = _alloca(sizeof(double)*lNFixDates);
	double *dFwdSwapRates = _alloca(sizeof(double)*lNFixDates);
	double *dFixCoverages =_alloca(sizeof(double)*lNFixDates);
	double *dShifts = _alloca(sizeof(double)*lNFixDates);
	double *dVols = _alloca(sizeof(double)*lNFixDates);
	double *dStrikes = _alloca(sizeof(double)*lNFixDates);

	// funding leg
	double *dFloatCoverages=_alloca(sizeof(double)*lNFloatDates);
	
	SrtCompounding srtFloatFreq;
	SrtBasisCode srtFloatBasis;
	char *cBasis, *cFreq;

	// determine the # of factors, i.e Brownian motions
	const int nForwardDim = lNFixDates; // # of forwards
	const int nFactorMax = 3; // max # of factors
	const int nFactorDim = nForwardDim>nFactorMax?nFactorMax:nForwardDim; // # of factors
	const double dIntDim = ((double)(nForwardDim-1)/(nFactorDim));
	int *pCutDim = _alloca(sizeof(int)*nFactorDim);
	
	// for simulation 
	double **dCopulaCube= 0;
	double *dPath = 0;

	/// Each bucket X(T)=X(0) +A*N^2 + B*N + C
	double *pdA = memset(_alloca(sizeof(double)*nFactorDim),0,sizeof(double)*nFactorDim);
	double *pdB = memset(_alloca(sizeof(double)*nFactorDim),0,sizeof(double)*nFactorDim);
	double *pdC = memset(_alloca(sizeof(double)*nFactorDim),0,sizeof(double)*nFactorDim);

	if(err = translate_basis(&cBasis, srtBasis)) return err;
	if(err = err = translate_compounding(&cFreq, srtFreq)) return err;
	if(err = err = swp_f_get_ref_rate_details(cRefRname, &srtFloatBasis, &srtFloatFreq)) return err;
	
	for(i=0;i<lNFixDates;++i)
		dFixCoverages[i] = coverage(lFixStartDates[i], lFixEndDates[i], srtBasis);
	
	for(i=0;i<lNFloatDates;++i)
		dFloatCoverages[i] = coverage(lFloatStartDates[i], lFloatEndDates[i], srtFloatBasis);
	
	/// produce the replicating portfolio
	err = ConvertAmortSwapWithMarginsInCoInitSwapPortfolioForMAD2(
								 cYCname,
								 cVCname,
								 cRefRname,

								 srtCallPut,
								 exer_fee,

								 srtFreq,
								 srtBasis,
								 
								 lNFixDates,
								 lFixStartDates,
								 lFixEndDates,
								 dFixCoverages,
								 dFixNotionals,
								 dFixRates,

								 lNFloatDates,
								 lFloatStartDates,
								 lFloatEndDates,
								 dFloatCoverages,
								 dFloatNotionals,
								 dFloatSpreads,
								 dFloatMargins,

								 dStrikes,
								 dCoInitNotionals,

								 dLvls,
								 dFwdSwapRates,
								 
								 UseVol);
	if (err)
	{
		return err;
	}

	/// compute instrinsic value
	dIV = 0.;
	for(i=0;i<nForwardDim;++i)
		dIV +=dCoInitNotionals[i] * dLvls[i]*(dFwdSwapRates[i]-dStrikes[i]);

	if( add_unit (lToday, 2, SRT_BDAY, MODIFIED_SUCCEEDING) >= lFixStartDates[0])
	{
		*dPrice = DMAX(0.0, flag * dIV);
		return 0;
	}


//---------------------------------------------------------------------------------
//-----------------Retrieves the densities ----------------------------------------
//---------------------------------------------------------------------------------
	err = Fwd_Swaption_Get_Cum_ShiftedLogDensity_ForMAD(
								 lToday,
								 cYCname,
								 cVCname,
								 cRefRname,
								 srtFreq,
								 srtBasis,

								 lNFixDates,
                                 lFixStartDates,
                                 lFixEndDates,

								 dStrikes,

								 iCalibShift,
								 Shift,
								 dShifts,
								 dVols);	
	if (err)
	{
		return err;
	}

	// compute the bucketing of swap rates
	for(j=0;j<nFactorDim;++j)
		pCutDim[j]=(int)(dIntDim*(j+1));
	_ASSERTE(nFactorDim>0);
	_ASSERTE(pCutDim[nFactorDim-1]==nForwardDim-1);

	// simulate normal deviates
	dCopulaCube= dmatrix(0,xlNumSim-1,0,nFactorDim-1);
	dPath = dvector(0,xlNumSim-1);
	if(!dCopulaCube||!dPath) 
		return "Memory allocation failed in AmortizedSwaptionShiftedLog";

	err = ShiftedLogCopulaSimul(Correl, xlNumSim, nFactorDim, pCutDim, dCopulaCube);
	if(err)
	{
		free_dvector(dPath,0,xlNumSim-1);
		free_dmatrix(dCopulaCube,0,xlNumSim-1,0,nFactorDim-1);
		return err;
	}
	
	*dPrice = 0.;

#ifdef _DEBUG
	debug = 0.;
#endif //_DEBUG
	
	// recompuate constants for each bucket
	for(i=0,j=0;j<nFactorDim;++j)
	{
		for(;i<=pCutDim[j];++i)	
		{
			const double dHalfVolSqrd = 0.5*dVols[i]*dVols[i];
			const double dWLFS = dCoInitNotionals[i] * dLvls[i]*(dFwdSwapRates[i]+dShifts[i]);
			const double dWLFSE = dWLFS*exp(-dHalfVolSqrd);
			
			pdA[j] += dWLFSE*dHalfVolSqrd;
			pdB[j] += dWLFSE*dVols[i];
			pdC[j] += dWLFSE-dWLFS;
		}
	}

	mean = 0.;
	for (i=0; i<xlNumSim; i++)
	{
		const double *pdRand=dCopulaCube[i];
		dPath[i]=0.;

		for(j=0;j<nFactorDim;++j)
		{
			dPath[i] += pdA[j]*pdRand[j]*pdRand[j]+pdB[j]*pdRand[j]+pdC[j];
		}
		mean +=dPath[i];
	}
	
	mean /= xlNumSim;
	mean -=dIV;// take the differerence 
	
	//-------MonteCarlo----------------------------------------------------------------
	for (i=0; i<xlNumSim; i++)
	{
		dPath[i]-=mean;

#ifdef _DEBUG			
		debug += dPath[i];
#endif //_DEBUG
		
		*dPrice += DMAX(0,flag * dPath[i]) ;
	}

#ifdef _DEBUG
	debug/=(xlNumSim);
	_ASSERTE(fabs(debug-dIV)<1.e-14);
#endif //_DEBUG

	*dPrice /= xlNumSim;

	free_dvector(dPath,0,xlNumSim-1);
	free_dmatrix(dCopulaCube,0,xlNumSim-1,0,nFactorDim-1);

	return err;
}
//--------------------------------------------------------------------------------------
//----------------Shifted Log Model : End Of Function For MAD---------------------------
//--------------------------------------------------------------------------------------
// previous version of the european amortizing simulator
Err		AmortizedSwaptionShiftedLogForMAD_NotOpt(
						 char *cYCname,
						 char *cVCname,
						 char *cRefRname,

						 SrtCompounding srtFreq,
						 SrtBasisCode srtBasis,

						 double exer_fee,
						 
						 long	*lFixStartDates,
						 long	*lFixEndDates,
						 long	lNFixDates,
						 double *dFixNotionals,
						 double *dFixRates,

						 long	*lFloatStartDates,
						 long	*lFloatEndDates,
						 long	lNFloatDates,
						 double *dFloatNotionals,
						 double *dFloatMargins,
						 double *dFloatSpreads,

						 SrtCallPutType srtCallPut,

						 int xlNumSim,

						 int iCalibShift,
						 double Shift,
						 int UseVol,

						 double **Correl,
						 double *dPrice
			            )
{
	
	///// Objects and variables declaration
	Err err = NULL;

	SrtCurvePtr	pCurve = lookup_curve (cYCname);

	Date lToday,lSpotLag;

	double **dCopulaCube=NULL,**dCopulaCubeShiftedLog=NULL,*dMean_v=NULL,*dAvg=NULL;

	double **dPPDx=NULL,**dPPDy=NULL,**dCPDy=NULL;
	double **dPPDShiftedLogx=NULL;

	int i, j;
	
	int flag;

	int nMonth;

	double dPayoff;

	Date *lCoInitSwapsEndDates=NULL;
	double *dCoInitStrikes=NULL;
	double *dCoInitNotionals=NULL;

	SrtCompounding srtFloatFreq;
	SrtBasisCode srtFloatBasis;

	double *dShifts=NULL;
	double *dVols=NULL;

	double *dStrikes=NULL;

	char *cBasis, *cFreq;

	double *dFixCoverages=NULL;
	double *dFloatCoverages=NULL;

//---------------------Variable for variate control----------------------
	double *mean=NULL;
	double *dLvls=NULL;
	double *dFwdSwapRates=NULL;
	double dCov, dFwdSwap, dFwdLvl, iv;

	long	*lFixStartDatesOut=NULL;
	long	*lFixEndDatesOut=NULL;
	long	*lFloatStartDatesOut=NULL;
	long	*lFloatEndDatesOut=NULL;

	lFixStartDatesOut = lngvector(0,lNFixDates-1);
	lFixEndDatesOut = lngvector(0,lNFixDates-1);

	lFloatStartDatesOut = lngvector(0,lNFloatDates-1);
	lFloatEndDatesOut = lngvector(0,lNFloatDates-1);

/*	err = amortMidat_modify_dates(
					lNFixDates,
					lFixStartDates,
					lFixEndDates,
					
					lNFloatDates,
					lFloatStartDates,
					lFloatEndDates,
					
					lFixStartDatesOut,
					lFixEndDatesOut,
					
					lFloatStartDatesOut,
					lFloatEndDatesOut);
*/
	err = translate_basis(&cBasis, srtBasis);
	err = translate_compounding(&cFreq, srtFreq);

	err = swp_f_get_ref_rate_details(cRefRname, &srtFloatBasis, &srtFloatFreq);
	/////////////////////////////////////////////////////////////////////////
	///// defines the pointer to the curve and retrieves today and spot dates
	/////////////////////////////////////////////////////////////////////////

	pCurve = lookup_curve (cYCname);
	lToday = (Date) get_today_from_curve (pCurve);
	lSpotLag = (Date) get_spotlag_from_curve(pCurve);	

	dFixCoverages = dvector(0,lNFixDates-1);
	for(i=0;i<lNFixDates;++i)
	{
//		dFixCoverages[i] = coverage(lFixStartDatesOut[i], lFixEndDatesOut[i], srtBasis);
		dFixCoverages[i] = coverage(lFixStartDates[i], lFixEndDates[i], srtBasis);
	}
	dFloatCoverages = dvector(0,lNFloatDates-1);
	for(i=0;i<lNFloatDates;++i)
	{
//		dFloatCoverages[i] = coverage(lFloatStartDatesOut[i], lFloatEndDatesOut[i], srtFloatBasis);
		dFloatCoverages[i] = coverage(lFloatStartDates[i], lFloatEndDates[i], srtFloatBasis);
	}


//-----------------------------------------------------------------------
//------------------------Calibration Strikes----------------------------
//-----------------------------------------------------------------------
	dStrikes = dvector(0,lNFixDates-1);
	if (!dStrikes)
	{
		err = "Memory allocation failed in AmortizedSwaptionShiftedLog";
		goto FREE_RETURN;
	}

//---------------------------------------------------------------------------------
//--------Convert AMort Swap in Coinit Swap Portfolio------------------------------
//---------------------------------------------------------------------------------
	dCoInitNotionals = dvector(0,lNFixDates-1);
	dLvls = dvector(0,lNFixDates-1);
	dFwdSwapRates = dvector(0,lNFixDates-1);
	err = ConvertAmortSwapWithMarginsInCoInitSwapPortfolioForMAD2(
								 cYCname,
								 cVCname,
								 cRefRname,

								 srtCallPut,
								 exer_fee,

								 srtFreq,
								 srtBasis,
								 
								 lNFixDates,
//								 lFixStartDatesOut,
//								 lFixEndDatesOut,
								 lFixStartDates,
								 lFixEndDates,
								 dFixCoverages,
								 dFixNotionals,
								 dFixRates,

								 lNFloatDates,
//								 lFloatStartDatesOut,
//								 lFloatEndDatesOut,
								 lFloatStartDates,
								 lFloatEndDates,
								 dFloatCoverages,
								 dFloatNotionals,
								 dFloatSpreads,
								 dFloatMargins,

								 dStrikes,
								 dCoInitNotionals,

								 dLvls,
								 dFwdSwapRates,
								 
								 UseVol);
	if (err)
	{
		goto FREE_RETURN;
	}

	flag = -1;
    if (srtCallPut == SRT_CALL) flag=1;

	if (cFreq == "ANNUAL") {
		dCov=1.0;
		nMonth = 12;

	} else if (cFreq == "SEMIANNUAL") {
		dCov=0.5;
		nMonth = 6;
	} else if (cFreq == "QUARTERLY") {
		dCov=0.25;
		nMonth = 3;
	} else if (cFreq == "MONTHLY") {
		dCov=1./12.;
		nMonth = 1;
	}

	mean = dvector(0,lNFixDates-1);
	iv = 0.0;
	for(j=0;j<lNFixDates;++j)
	{
		mean[j] = 0.0;
		iv += dCoInitNotionals[j] * dLvls[j] * ( dFwdSwapRates[j] - dStrikes[j]);
	}

//	if( add_unit (lToday, 2, SRT_BDAY, MODIFIED_SUCCEEDING) >= lFixStartDatesOut[0])
	if( add_unit (lToday, 2, SRT_BDAY, MODIFIED_SUCCEEDING) >= lFixStartDates[0])
	{
		*dPrice = DMAX(0.0, flag * iv);
		goto FREE_RETURN;
	}


//---------------------------------------------------------------------------------
//-----------------Retrieves the densities ----------------------------------------
//---------------------------------------------------------------------------------
	dShifts = dvector(0,lNFixDates-1);
	dVols = dvector(0,lNFixDates-1);
	err = Fwd_Swaption_Get_Cum_ShiftedLogDensity_ForMAD(
								 lToday,
								 cYCname,
								 cVCname,
								 cRefRname,
								 srtFreq,
								 srtBasis,

								 lNFixDates,
//                               lFixStartDatesOut,
//                               lFixEndDatesOut,
                                 lFixStartDates,
                                 lFixEndDates,

								 dStrikes,

								 iCalibShift,
								 Shift,
								 dShifts,
								 dVols);	
	if (err)
	{
		goto FREE_RETURN;
	}

	
//---------------------------------------------------------------------------------
//--------Copula Retrieves the Matrix of Joint distributions-----------------------
//---------------------------------------------------------------------------------
	dMean_v=dvector(0,lNFixDates-1);
	dAvg=dvector(0,lNFixDates-1);
	for (i=0; i<lNFixDates; i++) {
		dAvg[i]=0.0;
		dMean_v[i]=0.0;  
	}

	dCopulaCube = dmatrix(0,xlNumSim,0,lNFixDates-1);

	ShiftedLogCopulaSimul_NotOpt(Correl, xlNumSim, lNFixDates, dFwdSwapRates, dShifts, dVols, dCopulaCube);

	for (i=3; i<xlNumSim; i++)
	{
		for(j=0;j<lNFixDates;++j)
		{
			mean[j] += dCoInitNotionals[j] * dLvls[j] *(dCopulaCube[i][j]-dStrikes[j]);
		}
	}

	for(j=0;j<lNFixDates;++j)
	{
		dFwdSwap = dFwdSwapRates[j];
		dFwdLvl = dLvls[j];
		mean[j] = dCoInitNotionals[j] * dFwdLvl * ( dFwdSwap - dStrikes[j]) - mean[j] / (xlNumSim-2);
	}


	*dPrice=0.0;

	//-------MonteCarlo----------------------------------------------------------------
	for (i=3; i<xlNumSim; i++)
	{
		err = GetAmortizedSwaptionPayoff_ForMAD(dCopulaCube,
											srtFreq,
											srtBasis,
											lNFixDates,
											dLvls,
											dStrikes,
											dCoInitNotionals,
											i,
											&dPayoff,
											mean,
											flag);

		*dPrice += (DMAX(0,flag * dPayoff) - *dPrice)/(i-2);
	}


	//-----------------Frees the memory and return-------------------------------------

FREE_RETURN:

	if(dFixCoverages) free_dvector(dFixCoverages, 0,lNFixDates-1);
	if(dFloatCoverages) free_dvector(dFloatCoverages, 0,lNFloatDates-1);

	if(dStrikes) free_dvector(dStrikes, 0, lNFixDates-1);

	if(dShifts) free_dvector(dShifts, 0, lNFixDates-1);
	if(dVols) free_dvector(dVols, 0, lNFixDates-1);

	if(dMean_v) free_dvector(dMean_v,0,lNFixDates-1);
	if(dAvg) free_dvector(dAvg,0,lNFixDates-1);
	if(dCopulaCube) free_dmatrix(dCopulaCube,0,xlNumSim,0,lNFixDates-1);

	if(dCoInitStrikes) free(dCoInitStrikes);
	if(dCoInitNotionals) free_dvector(dCoInitNotionals, 0, lNFixDates-1);

	if(mean) free_dvector(mean, 0,lNFixDates-1);
	if(dLvls) free_dvector(dLvls, 0,lNFixDates-1);
	if(dFwdSwapRates) free_dvector(dFwdSwapRates, 0,lNFixDates-1);

//	if(lFixStartDatesOut) free_lngvector(lFixStartDatesOut, 0,lNFixDates);
//	if(lFixEndDatesOut) free_lngvector(lFixEndDatesOut, 0,lNFixDates);

//	if(lFloatStartDatesOut) free_lngvector(lFloatStartDatesOut, 0,lNFloatDates);
//	if(lFloatEndDatesOut) free_lngvector(lFloatEndDatesOut, 0,lNFloatDates);


	return err;
}


//--------------------------------------------------------------------------------------
//---------------------Shifted Log Model : Main Function--------------------------------
//--------------------------------------------------------------------------------------
Err		AmortizedSwaptionShiftedLog(
						 char *cYCname,
						 char *cVCname,
						 char *cRefRname,

						 long xlStartDate,
						 long xlEndDate,

						 SrtCompounding srtFreq,
						 SrtBasisCode srtBasis,

						 double exer_fee,
						 
						 long	lNFixNot,
						 double *dFixNotionals,
						 double *dFixRates,

						 long	lNFloatNot,
						 double *dFloatNotionals,
						 double *dMargins,

						 SrtCallPutType srtCallPut,

						 int xlStudDegree,
						 int xlNumSim,
						 double xlUpBound,
						 int xlnClasses,

						 int UseVol,

						 double **Correl,
						 double *dPrice
			            )
{
	
	///// Objects and variables declaration
	Err err = NULL;

	SrtCurvePtr	pCurve = lookup_curve (cYCname);

	SwapDP Swap;

	Date lToday,lSpotLag;

	Date *lPayDates=NULL,*lStartDates=NULL,*lEndDates=NULL;
	double *dCoverages=NULL;
	int iNDates, iNPayDates;

	double **dCopulaCube=NULL,**dCopulaCubeShiftedLog=NULL,*dMean_v=NULL,*dAvg=NULL;

	double **dPPDx=NULL,**dPPDy=NULL,**dCPDy=NULL;
	double **dPPDShiftedLogx=NULL;

	int i, j;
	
	int nClass=1<<xlnClasses;
	
	int flag;

	int nMonth;

	double dPayoff;
	
	Date *lCpnDates=NULL;
	double *dCpns=NULL;

	Date *lCoInitSwapsEndDates=NULL;
	double *dCoInitStrikes=NULL;
	double *dCoInitNotionals=NULL;

	double *dSpreads=NULL;

	double *dShifts=NULL;

	double *dVols=NULL;

	double *dStrikes=NULL;

	char *cBasis, *cFreq;

	int shortstub;

	//---Variable for variate control--------------------------------------------------
	double *mean=NULL;
	double *dLvls=NULL;
	double *dFwdSwapRate=NULL;
	double dCov, dFwdSwap, dFwdLvl, iv;

	err = translate_basis(&cBasis, srtBasis);
	err = translate_compounding(&cFreq, srtFreq);
	/////////////////////////////////////////////////////////////////////////
	///// defines the pointer to the curve and retrieves today and spot dates
	/////////////////////////////////////////////////////////////////////////

	pCurve = lookup_curve (cYCname);
	///// if (!pCurve) throw XSort("Fatal: (FwdSwaption) Yield Curve Object not found");
	lToday = (Date) get_today_from_curve (pCurve);
	lSpotLag = (Date) get_spotlag_from_curve(pCurve);	
	
	///////////////////////////////////////////////////////////////
	///// defines the relevant swap objects
	///////////////////////////////////////////////////////////////
    err = swp_f_setSwapDP(xlStartDate, xlEndDate, srtFreq, srtBasis, &Swap);
	if (err)
	{
		goto FREE_RETURN;
	}
	
	err = swp_f_make_FixedLegDatesAndCoverages(&Swap,lToday,&lPayDates,&iNPayDates,
								  			   &lStartDates,&lEndDates,
								  			   &dCoverages,&iNDates);
	if (err)
	{
		goto FREE_RETURN;
	}

//-----------------------------------------------------------------------
//------------------------Compute Calibration Strikes--------------------
//-----------------------------------------------------------------------
	dStrikes = dvector(0,iNDates-1);
	if (!dStrikes)
	{
		err = "Memory allocation failed in AmortizedSwaptionShiftedLog";
		goto FREE_RETURN;
	}
//---------------------------------------------------------------------------------
//---------------------End Of Compute Calibration Strikes--------------------------
//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
//--------Convert AMort Swap in Coinit Swap Portfolio------------------------------
//---------------------------------------------------------------------------------
	dCoInitNotionals = dvector(0,iNDates-1);
	dLvls = dvector(0,lNFixNot-1);
	dFwdSwapRate = dvector(0,lNFixNot-1);
	err = ConvertAmortSwapWithMarginsInCoInitSwapPortfolio2(
								 cYCname,
								 cVCname,
								 cRefRname,
								 srtCallPut,
								 xlStartDate,
								 xlEndDate,
								 srtFreq,
								 srtBasis,
								 exer_fee,
								 lNFixNot,
								 dFixNotionals,
								 dFixRates,
								 lNFloatNot,
								 dFloatNotionals,
								 dMargins,
								 dStrikes,
								 dCoInitNotionals,
								 dLvls,
								 dFwdSwapRate,
								 UseVol
								 );

	if (err)
	{
		goto FREE_RETURN;
	}

	flag = -1;
    if (srtCallPut == SRT_CALL) flag=1;

	if (cFreq == "ANNUAL") {
		dCov=1.0;
		nMonth = 12;

	} else if (cFreq == "SEMIANNUAL") {
		dCov=0.5;
		nMonth = 6;
	} else if (cFreq == "QUARTERLY") {
		dCov=0.25;
		nMonth = 3;
	} else if (cFreq == "MONTHLY") {
		dCov=1./12.;
		nMonth = 1;
	}

	shortstub = 0;
	if( (iNDates>1) 
			&&
	( (int)(12*(lEndDates[0]-lStartDates[0])/365.0 + 0.5)!=(int)(12*(lEndDates[1]-lStartDates[1])/365.0+0.5) ) )
	{
		shortstub = 1;
	}

	mean = dvector(0,lNFixNot-1);
	iv = 0.0;
	for(j=0;j<lNFixNot;++j)
	{
		mean[j] = 0.0;
		iv += dCoInitNotionals[j] * dLvls[j] * ( dFwdSwapRate[j] - dStrikes[j]);
	}

	if( add_unit (lToday, 2, SRT_BDAY, MODIFIED_SUCCEEDING) >= xlStartDate)
	{
		*dPrice = DMAX(0.0, flag * iv);
		goto FREE_RETURN;
	}


//---------------------------------------------------------------------------------
//-----------------Retrieves the densities ----------------------------------------
//---------------------------------------------------------------------------------
	dShifts = dvector(0,iNDates-1);
	dVols = dvector(0,iNDates-1);
	err = Fwd_Swaption_Get_Cum_ShiftedLogDensity_ForMAD(
								 lToday,
								 cYCname,
								 cVCname,
								 cRefRname,
								 srtFreq,
								 srtBasis,

								 iNPayDates-1,
                                 lStartDates,
                                 lEndDates,

								 dStrikes,

								 0,
								 0,
								 dShifts,
								 dVols);	
	if (err)
	{
		goto FREE_RETURN;
	}

//---------------------------------------------------------------------------------
//--------Copula Retrieves the Matrix of Joint distributions-----------------------
//---------------------------------------------------------------------------------
	dMean_v=dvector(0,iNDates-1);
	dAvg=dvector(0,iNDates-1);
	for (i=0; i<iNDates; i++) {
		dAvg[i]=0.0;
		dMean_v[i]=0.0;  
	}

	dCopulaCube = dmatrix(0,xlNumSim,0,iNDates-1);

	ShiftedLogCopulaSimul_NotOpt(Correl, xlNumSim, iNDates, dFwdSwapRate, dShifts, dVols, dCopulaCube);

	for (i=3; i<xlNumSim; i++)
	{
		for(j=0;j<iNDates;++j)
		{
			mean[j] += dCoInitNotionals[j] * dLvls[j] *(dCopulaCube[i][j]-dStrikes[j]);
		}
	}

	for(j=0;j<iNDates;++j)
	{
		dFwdSwap = dFwdSwapRate[j];
		dFwdLvl = dLvls[j];
		mean[j] = dCoInitNotionals[j] * dFwdLvl * ( dFwdSwap - dStrikes[j]) - mean[j] / (xlNumSim-2);
	}


	*dPrice=0.0;

	//-------MonteCarlo----------------------------------------------------------------
	for (i=3; i<xlNumSim; i++)
	{
		err = GetAmortizedSwaptionPayoff_ForMAD(dCopulaCube,
											srtFreq,
											srtBasis,
											lNFixNot,
											dLvls,
											dStrikes,
											dCoInitNotionals,
											i,
											&dPayoff,
											mean,
											flag);

		*dPrice += (DMAX(0,flag * dPayoff) - *dPrice)/(i-2);
	}


	//-----------------Frees the memory and return-------------------------------------

FREE_RETURN:
	if(dStrikes) free_dvector(dStrikes, 0, iNDates-1);

	if(dShifts) free_dvector(dShifts, 0, iNDates-1);
	if(dVols) free_dvector(dVols, 0, iNDates-1);

	if(dMean_v) free_dvector(dMean_v,0,iNDates-1);
	if(dAvg) free_dvector(dAvg,0,iNDates-1);
	if(dCopulaCube) free_dmatrix(dCopulaCube,0,xlNumSim,0,iNDates-1);

	if(dCoInitStrikes) free(dCoInitStrikes);
	if(dCoInitNotionals) free_dvector(dCoInitNotionals, 0, iNDates-1);

	if(mean) free_dvector(mean, 0,iNDates-1);
	if(dLvls) free_dvector(dLvls, 0,iNDates-1);
	if(dFwdSwapRate) free_dvector(dFwdSwapRate, 0,iNDates-1);

	return err;
}
//--------------------------------------------------------------------------------------
//----------------Shifted Log Model : End Of Main Function------------------------------
//--------------------------------------------------------------------------------------

Err		PriceAmortizedSwaptionShiftedLog(
						 char *cYCname,
						 char *cVCname,
						 char *cRefRname,

						 long xlStartDate,
						 long xlEndDate,

						 SrtCompounding srtFreq,
						 SrtBasisCode srtBasis,
						 
						 long	lNFixNot,
						 double *dFixNotionals,
						 long	lNFloatNot,
						 double *dFloatNotionals,
						 double *dFixRates,

						 SrtCallPutType srtCallPut,

						 int xlStudDegree,
						 int xlNumSim,
						 double xlUpBound,
						 int xlnClasses,
						 double **Correl,
						 double *dPrice,
						 int iPayoffMethod,
						 int iDistribMethod
			            )
{
	Err err=NULL;

	err = "Function PriceAmortizedSwaptionShiftedLog no more supported, please use AmortizedSwaptionShiftedLogForMAD";

	return err;
}

Err		PriceAmortizedSwaptionSABR(
						 char *cYCname,
						 char *cVCname,
						 char *cRefRname,

						 long xlStartDate,
						 long xlEndDate,

						 SrtCompounding srtFreq,
						 SrtBasisCode srtBasis,
						 
						 long	lNFixNot,
						 double *dFixNotionals,
						 long	lNFloatNot,
						 double *dFloatNotionals,
						 double *dFixRates,

						 SrtCallPutType srtCallPut,

						 int xlStudDegree,
						 int xlNumSim,
						 double xlUpBound,
						 int xlnClasses,
						 double **Correl,
						 double *dPrice,
						 int iPayoffMethod,
						 int iDistribMethod
			            )
{
	Err err=NULL;

	err = "Function PriceAmortizedSwaptionSABR no more supported, please use AmortizedSwaptionShiftedLogForMAD";

	return err;
}


//--------------------------------------------------------------------------------------
//----------------I don't like to delete code... ---------------------------------------
//----------------Below, among others pricing function with Heston model ---------------
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

/*
void CopulaSimul(
						double **corr_mtx,   
						const long p,
						const int d, 
						double **xa,
						double **ya,
						const long n_pts,
						double **res
					   )
{
	// NB d is the total numbers of dimensions: d=dim space

	long i,idum=-8935807;
	int j,k;
	double	***GaussSample=NULL,**UnifSample=NULL,
			**sqrt_corr_mtx=NULL,
			IntermRes,q;
				
		//////////  Memory allocation
		GaussSample=dcube(0,p,0,d-1,0,0);
		sqrt_corr_mtx=dmatrix(0,d-1,0,d-1);

		PositiveMatrix(corr_mtx, d);

		//////////  Step1: performs Chol decomposition of the correlation matrix corr_mtx
		CholDec	(&CholAlg,
				 corr_mtx,       
				 sqrt_corr_mtx,
				 d);

		//////////  Step2: Generates the matrix of r.v. 
		UnifSample=dmatrix(0,p,0,d-1);
		GetSobolMatrix(p,d,UnifSample);
			
		//////////  Step3: Generates the matrix of r.v. compatible with the Copula and the marginals
		for (i=0;i<p;i++)
		{	
			for (j=0; j<d; j++)
			{
				GaussSample[i][j][0]=inv_cumnorm_fast(UnifSample[i][j]);
			}
			
			for (j=0; j<d; j++)
			{
				res[i][j]=0.0;
				for (k=0; k<d ;k++)
				{
					res[i][j]+=sqrt_corr_mtx[j][k]*GaussSample[i][k][0];
				}
				res[i][j] = norm_accurate(res[i][j]);
			}

			for (j=0; j<d; j++)
			{
				IntermRes=res[i][j];
				res[i][j]=interp_columns(ya, 
										xa, 
										n_pts, 
										IntermRes, 
										0, 
										&q,
										j);
			}

		}
	free_dmatrix(UnifSample,0,p,0,d-1);
				
	free_dcube(GaussSample,0,p,0,d-1,0,0);
	free_dmatrix(sqrt_corr_mtx,0,d-1,0,d-1);

}
*/

/*
double amortSwap_PV(long lToday, char *cYCname, int iNCpns, long *lCpnDates, double *dCpns)
{
	int i;
	double PV;

	PV = 0;
	for(i=0;i<iNCpns;++i)
	{
		PV = PV + dCpns[i] * swp_f_df(lToday, lCpnDates[i], cYCname);
	}

	return PV;
}

*/

/*
//----------------------------------------------------------------
//---5) Coinit Swaps, level cash & no stripping ------------------
//----------------------------------------------------------------
Err		GetAmortizedSwaptionPayoff5(double  **dCopulaCube,
							 SrtCompounding srtFreq,
							 SrtBasisCode srtBasis,
							 double	 *dCoverages,
			                 int     iNCoInitSwaps,
							 double  *dStrikes,
							 double  *dNotionals,
							 int     iMC,
							 double  *dPayoff,
							 int flag)
{
	Err err = NULL;
	double dLevelCash;
	int i;
	double dCov;
	char *cBasis, *cFreq;

	err = translate_basis(&cBasis, srtBasis);
	err = translate_compounding(&cFreq, srtFreq);

	if (cFreq == "ANNUAL") {
		dCov=1.0;
	} else if (cFreq == "SEMIANNUAL") {
		dCov=0.5;
	} else if (cFreq == "QUARTERLY") {
		dCov=0.25;
	} else if (cFreq == "MONTHLY") {
		dCov=1./12.;
	}


//	*dPayoff = 0;
//	dLevelCash = 0;
//	df = 1.0;
//	for(i=0;i<iNCoInitSwaps;++i)
//	{
//		df = df / (1 + dCoverages[i] * dCopulaCube[iMC][i]);
//		dLevelCash = dLevelCash + dCoverages[i] * df;
//		*dPayoff += dNotionals[i] * dLevelCash *(dCopulaCube[iMC][i]-dStrikes[i]);
//	}


	*dPayoff = 0;
	for(i=0;i<iNCoInitSwaps;++i)
	{
		dLevelCash =(1.-1./pow(1+dCopulaCube[iMC][i]*dCov,i+1))/dCopulaCube[iMC][i];
		*dPayoff += dNotionals[i] * dLevelCash *(dCopulaCube[iMC][i]-dStrikes[i]);
	}

	return err;
}


Err		GetAmortizedSwaptionPayoff5_c(double  **dCopulaCube,
							 SrtCompounding srtFreq,
							 SrtBasisCode srtBasis,
							 double	 *dCoverages,
			                 int     iNCoInitSwaps,
							 double  *dLvls,
							 double  *dStrikes,
							 double  *dNotionals,
							 int     iMC,
							 double  *dPayoff,
							 double	*controlvar,
							 int flag)
{
	Err err = NULL;
	double dLevelCash;
	int i;
	double dCov;
	char *cBasis, *cFreq;

	err = translate_basis(&cBasis, srtBasis);
	err = translate_compounding(&cFreq, srtFreq);

	if (cFreq == "ANNUAL") {
		dCov=1.0;
	} else if (cFreq == "SEMIANNUAL") {
		dCov=0.5;
	} else if (cFreq == "QUARTERLY") {
		dCov=0.25;
	} else if (cFreq == "MONTHLY") {
		dCov=1./12.;
	}

	*dPayoff = 0;
	for(i=0;i<iNCoInitSwaps;++i)
	{
		dLevelCash =(1.-1./pow(1+dCopulaCube[iMC][i]*dCov,i+1))/dCopulaCube[iMC][i];
		dLevelCash =dLvls[i];
		*dPayoff += dNotionals[i] * dLevelCash *(dCopulaCube[iMC][i]-dStrikes[i]) + controlvar[i];
	}

	return err;
}
*/

/*

Err  GetShiftedLogDensity_Lvl_Measure6(
						long today,
						int		nClass,
						char	*cYCname,
						char	*cVCname,
						char	*cRefRname,
						SrtCompounding srtFreq,
						SrtBasisCode srtBasis,
						long start,
						long end,
						long nTen,
						double dCov,
						double dStrike,
						double dDf,
						double UpBound,
						double **dPPDx,
						double **dPPDy,
						double **dCPDy,
						double *dShift
						)
{
	double Df_x, LvlC_x, LvlCFwd;
	double *dVec = NULL, *dX=NULL,*dY=NULL, dThreshold=1.e-3;
	int i;
	double deriv1, deriv2;
	Err err = NULL;
	double epsilon = 1.0e-4;
	double temp, temp_plus, temp_moins;
	double dK;
	double dFwdSwap, dFwdCash, dSpread, dLvl, ex_time;
	double sigma, alpha, beta, rho;
	char *cBasis, *cFreq;
	double power, atmvol;
	double shift, shifted_vol, equi_beta;
	double k1, k2, vol1, vol2;
	SrtCallPutType callput;
	double sgn;
	
	err = translate_basis(&cBasis, srtBasis);
	err = translate_compounding(&cFreq, srtFreq);

	err = swp_f_ForwardRate(start, end, cFreq, cBasis, cYCname, cRefRname, &dFwdSwap);
	dFwdCash = swp_f_swapcashrate( start, end, srtBasis, srtFreq, cYCname, cRefRname);
	err = swp_f_LevelPayment(start, end, cFreq, cBasis, cYCname, cRefRname, &dLvl);
	ex_time = (add_unit (start, - 2, SRT_BDAY, MODIFIED_SUCCEEDING) - today) * YEARS_IN_DAY;

	dSpread = dFwdSwap - dFwdCash;

	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&sigma, &power, SABR_BETAVOL );
	if (err) return err;
	if(sigma==0)
	{
		err = "Need a SABR Market";
		return err;
	}
	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&alpha, &power, SABR_ALPHA );
	if (err) return err;
	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&beta, &power, SABR_BETA );
	if (err) return err;
	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&rho, &power, SABR_RHO );
	if (err) return err;

	err = srt_f_optsarbvol(
				dFwdSwap,
				dFwdSwap, 
				ex_time, 
				sigma,  
				alpha,
				beta,
				rho,
				SRT_BETAVOL,
				SRT_LOGNORMAL,
				&atmvol);

	k1 = dFwdSwap;
	k2 = DMAX(dFwdSwap*exp(-3*sqrt(ex_time)*atmvol), dStrike);
	if( fabs(k2 - k1) < 0.001 )
	{
		k1 = k2 + 0.001;
	}
	err = srt_f_optsarbvol(
				dFwdSwap,
				k1, 
				ex_time, 
				sigma,  
				alpha,
				beta,
				rho,
				SRT_BETAVOL,
				SRT_LOGNORMAL,
				&vol1);

	err = srt_f_optsarbvol(
				dFwdSwap,
				k2, 
				ex_time, 
				sigma,  
				alpha,
				beta,
				rho,
				SRT_BETAVOL,
				SRT_LOGNORMAL,
				&vol2);

	err = find_shifted_log_params(
					ex_time,
					dFwdSwap,
					k1,
					vol1,
					k2,
					vol2,
					&shift,
					&shifted_vol,
					&equi_beta);

//	if(fabs(shift) > 10000)
//	{
//		shift = 10000;
//		price2 = srt_f_optblksch(dFwdSwap, k2, vol2, ex_time, 1, SRT_CALL, SRT_PREMIUM);
//		err = srt_f_optimpvol(price2, shift+dFwdSwap, shift+k2, ex_time, 1, SRT_CALL, SRT_LOGNORMAL, &shifted_vol);
//	}

	callput = SRT_CALL;
	sgn = 1.0;
	if(shifted_vol < 0)
	{
		sgn = -1.0;
		callput = SRT_PUT;
		shifted_vol = - shifted_vol;
		shift = - shift;
	}

	*dShift = shift;

	///// gets the two boundaries
	if(sgn > 0)
	{
		dPPDx[nTen][0] = DMAX(-UpBound, -shift + 0.0001);
		dPPDx[nTen][2*nClass] = UpBound;
	}
	else
	{
		dPPDx[nTen][0] = -UpBound;
		dPPDx[nTen][2*nClass] = DMIN(shift - 0.0001, UpBound);
	}

	dK = (dPPDx[nTen][2*nClass]-dPPDx[nTen][0])/(2*nClass);

	for (i=1;i<2*nClass;i++)
	{
		dPPDx[nTen][i]=dPPDx[nTen][0] + i * dK;
	}

	///// runs through the points 
	for (i=0;i<=2*nClass;i++)
	{
		dCPDy[nTen][i]=0.0;
		dPPDy[nTen][i]=0.0;
	}

	////i=0
	temp = srt_f_optblksch(shift + sgn * dFwdSwap,
								shift + sgn * dPPDx[nTen][0],
								shifted_vol,
								ex_time,
								1.0,
								callput, 
								PREMIUM);

	temp_plus = srt_f_optblksch(shift + sgn * dFwdSwap,
								shift + sgn * dPPDx[nTen][1],
								shifted_vol,
								ex_time,
								1.0,
								callput, 
								PREMIUM);

	LvlC_x = (1.-1./pow(1+dPPDx[nTen][0]*dCov,nTen+1))/dPPDx[nTen][0];
	LvlCFwd = (1.-1./pow(1+dFwdSwap*dCov,nTen+1))/dFwdSwap;
	Df_x = 1;//./(1+dPPDx[nTen][0]*dCov);

	deriv1 = 1.0;
	deriv2 = DMIN(deriv1,(temp - temp_plus)/dK);

	dPPDy[nTen][0]=DMAX(0,(deriv1 - deriv2)/dK)/Df_x;//*LvlCFwd/LvlC_x;
	dCPDy[nTen][0]=0.0;
	
	for (i=1;i<2*nClass;i++)
	{
		temp_moins = temp;
		temp = temp_plus;

		temp_plus = srt_f_optblksch(shift + sgn * dFwdSwap,
								shift + sgn * dPPDx[nTen][i+1],
								shifted_vol,
								ex_time,
								1.0,
								callput, 
								PREMIUM);

		LvlC_x = (1.-1./pow(1+dPPDx[nTen][i]*dCov,nTen+1))/dPPDx[nTen][i];
		Df_x = 1;//./(1+dPPDx[nTen][i]*dCov);

		deriv1 = DMIN(1.0,(temp_moins - temp)/dK);
		deriv2 = DMIN(deriv1,(temp - temp_plus)/dK);

		dPPDy[nTen][i]=DMAX(0,(deriv1 - deriv2)/dK)/Df_x;//*LvlCFwd/LvlC_x;

		dCPDy[nTen][i]=dCPDy[nTen][i-1]+(dPPDy[nTen][i]+dPPDy[nTen][i-1])/2.*dK;
	}

	temp_moins = temp;
	temp = temp_plus;

	LvlC_x = (1.-1./pow(1+dPPDx[nTen][2*nClass]*dCov,nTen+1))/dPPDx[nTen][2*nClass];
	Df_x = 1;//./(1+dPPDx[nTen][2*nClass]*dCov);

	deriv1 = (temp_moins - temp)/(dPPDx[nTen][2*nClass]-dPPDx[nTen][2*nClass-1]);
	deriv2 = 0.0;

	dPPDy[nTen][2*nClass]=DMAX(0,(deriv1 - deriv2)/dK)/Df_x;//*LvlCFwd/LvlC_x;

	dCPDy[nTen][2*nClass]=dCPDy[nTen][2*nClass-1]
		+(dPPDy[nTen][2*nClass]+dPPDy[nTen][2*nClass-1])/2.*dK;

	for (i=0;i<=2*nClass;i++) 
	{
		dCPDy[nTen][i]/=dCPDy[nTen][2*nClass];
	}

	return err;

}
*/

/*
Err       Fwd_Swaption_Get_Cum_ShiftedLogDensity(Date lToday,
								 char *cYCname,
								 char *cVCname,
						         char *cRefRname,
								 SrtCompounding srtFreq,
								 SrtBasisCode srtBasis,

								 int iNPayDatesLong,
								 Date *lPayDatesLong,
								 double *dCoveragesLong,

								 double *dStrikes,

								 double *dShifts,

								 double dUpBound,
								 int nClass,

								 double **dPPDx,
                                 double **dPPDy,
							     double **dCPDy
								 )

{
	Err err = NULL;
	int i;
	double dMat,dCov;
	char *cBasis, *cFreq;
//	double shiftHeston;
	double dDf;
	double shift;

	dMat=(lPayDatesLong[0]-lToday)/365.; 
//		dDisc=swp_f_df(lToday, lPayDatesLong[0], cYCname);
	
//////////////////////////////////////////
/// translate SORT types into strings
//////////////////////////////////////////

	err = translate_basis(&cBasis, srtBasis);
	err = translate_compounding(&cFreq, srtFreq);

	if (cFreq == "ANNUAL") {
		dCov=1.0;
	} else if (cFreq == "SEMIANNUAL") {
		dCov=0.5;
	} else if (cFreq == "QUARTERLY") {
		dCov=0.25;
	} else if (cFreq == "MONTHLY") {
		dCov=1./12.;
	}

//////////////////////////////////////////////////////////
////////////////// generates the cumulative densities
//////////////////////////////////////////////////////////	

	dDf=swp_f_df(lToday, lPayDatesLong[0], cYCname);
	for (i=1;i<iNPayDatesLong;i++)
	{
		err = GetShiftedLogDensity_Lvl_Measure6(lToday, nClass,
										cYCname,cVCname,cRefRname,
										srtFreq,srtBasis,
										lPayDatesLong[0],
										lPayDatesLong[i],
										i-1,
										dCov,
										dStrikes[i-1],
										dDf,
										dUpBound,
										dPPDx,
										dPPDy,
										dCPDy,
										&shift);
		dShifts[i-1] = shift;
	}

	return err;

}	
*/

/*
//--------------------------------------------------------------------------------------
//----------------Shifted Log Model : Old Function--------------------------------------
//--------------------------------------------------------------------------------------
Err		AmortizedSwaptionShiftedLogOld(
						 char *cYCname,
						 char *cVCname,
						 char *cRefRname,

						 long xlStartDate,
						 long xlEndDate,

						 SrtCompounding srtFreq,
						 SrtBasisCode srtBasis,

						 double exer_fee,
						 
						 long	lNFixNot,
						 double *dFixNotionals,
						 double *dFixRates,

						 long	lNFloatNot,
						 double *dFloatNotionals,
						 double *dMargins,

						 SrtCallPutType srtCallPut,

						 int xlStudDegree,
						 int xlNumSim,
						 double xlUpBound,
						 int xlnClasses,

						 int UseVol,

						 double **Correl,
						 double *dPrice
			            )
{
	
	///// Objects and variables declaration
	Err err = NULL;

	SrtCurvePtr	pCurve = lookup_curve (cYCname);

	SwapDP Swap;

	Date lToday,lSpotLag;

	Date *lPayDates=NULL,*lStartDates=NULL,*lEndDates=NULL;
	double *dCoverages=NULL;
	int iNDates, iNPayDates;

	double **dCopulaCube=NULL,**dCopulaCubeShiftedLog=NULL,*dMean_v=NULL,*dAvg=NULL;

	double **dPPDx=NULL,**dPPDy=NULL,**dCPDy=NULL;
	double **dPPDShiftedLogx=NULL;

	int i, j;
	
	int nClass=1<<xlnClasses;
	
	int flag;

	int nMonth;

	double dPayoff;
	
	Date *lCpnDates=NULL;
	double *dCpns=NULL;

	Date *lCoInitSwapsEndDates=NULL;
	double *dCoInitStrikes=NULL;
	double *dCoInitNotionals=NULL;

	double *dSpreads=NULL;

	double *dShifts=NULL;

	double *dStrikes=NULL;

	char *cBasis, *cFreq;

	int shortstub;

	//---Variable for variate control--------------------------------------------------
	double *mean=NULL;
	double *dLvls=NULL;
	double *dFwdSwapRate=NULL;
	double dCov, dFwdSwap, dFwdLvl, iv;

	err = translate_basis(&cBasis, srtBasis);
	err = translate_compounding(&cFreq, srtFreq);
	/////////////////////////////////////////////////////////////////////////
	///// defines the pointer to the curve and retrieves today and spot dates
	/////////////////////////////////////////////////////////////////////////

	pCurve = lookup_curve (cYCname);
	///// if (!pCurve) throw XSort("Fatal: (FwdSwaption) Yield Curve Object not found");
	lToday = (Date) get_today_from_curve (pCurve);
	lSpotLag = (Date) get_spotlag_from_curve(pCurve);	
	
	///////////////////////////////////////////////////////////////
	///// defines the relevant swap objects
	///////////////////////////////////////////////////////////////
    err = swp_f_setSwapDP(xlStartDate, xlEndDate, srtFreq, srtBasis, &Swap);
	if (err)
	{
		goto FREE_RETURN;
	}
	
	err = swp_f_make_FixedLegDatesAndCoverages(&Swap,lToday,&lPayDates,&iNPayDates,
								  			   &lStartDates,&lEndDates,
								  			   &dCoverages,&iNDates);
	if (err)
	{
		goto FREE_RETURN;
	}

//-----------------------------------------------------------------------
//------------------------Compute Calibration Strikes--------------------
//-----------------------------------------------------------------------
	dStrikes = dvector(0,iNDates-1);
	if (!dStrikes)
	{
		err = "Memory allocation failed in AmortizedSwaptionShiftedLog";
		goto FREE_RETURN;
	}
//---------------------------------------------------------------------------------
//---------------------End Of Compute Calibration Strikes--------------------------
//---------------------------------------------------------------------------------

//---------------------------------------------------------------------------------
//--------Convert AMort Swap in Coinit Swap Portfolio------------------------------
//---------------------------------------------------------------------------------
	dCoInitNotionals = dvector(0,iNDates-1);
	dLvls = dvector(0,lNFixNot-1);
	dFwdSwapRate = dvector(0,lNFixNot-1);
	err = ConvertAmortSwapWithMarginsInCoInitSwapPortfolio2(
								 cYCname,
								 cVCname,
								 cRefRname,
								 srtCallPut,
								 xlStartDate,
								 xlEndDate,
								 srtFreq,
								 srtBasis,
								 exer_fee,
								 lNFixNot,
								 dFixNotionals,
								 dFixRates,
								 lNFloatNot,
								 dFloatNotionals,
								 dMargins,
								 dStrikes,
								 dCoInitNotionals,
								 dLvls,
								 dFwdSwapRate,
								 UseVol
								 );

	if (err)
	{
		goto FREE_RETURN;
	}

	flag = -1;
    if (srtCallPut == SRT_CALL) flag=1;

	if (cFreq == "ANNUAL") {
		dCov=1.0;
		nMonth = 12;

	} else if (cFreq == "SEMIANNUAL") {
		dCov=0.5;
		nMonth = 6;
	} else if (cFreq == "QUARTERLY") {
		dCov=0.25;
		nMonth = 3;
	} else if (cFreq == "MONTHLY") {
		dCov=1./12.;
		nMonth = 1;
	}

	shortstub = 0;
	if( (iNDates>1) 
			&&
	( (int)(12*(lEndDates[0]-lStartDates[0])/365.0 + 0.5)!=(int)(12*(lEndDates[1]-lStartDates[1])/365.0+0.5) ) )
	{
		shortstub = 1;
	}

	mean = dvector(0,lNFixNot-1);
	iv = 0.0;
	for(j=0;j<lNFixNot;++j)
	{
		mean[j] = 0.0;
		iv += dCoInitNotionals[j] * dLvls[j] * ( dFwdSwapRate[j] - dStrikes[j]);
	}

	if( add_unit (lToday, 2, SRT_BDAY, MODIFIED_SUCCEEDING) >= xlStartDate)
	{
		*dPrice = DMAX(0.0, flag * iv);
		goto FREE_RETURN;
	}


//---------------------------------------------------------------------------------
//-----------------Retrieves the densities ----------------------------------------
//---------------------------------------------------------------------------------
	dShifts = dvector(0,iNDates-1);
	dPPDx=dmatrix(0,iNDates-1,0,2*nClass);
	dPPDy=dmatrix(0,iNDates-1,0,2*nClass);
	dCPDy=dmatrix(0,iNDates-1,0,2*nClass);
	err = Fwd_Swaption_Get_Cum_ShiftedLogDensity(
								 lToday,
								 cYCname,
								 cVCname,
								 cRefRname,
								 srtFreq,
								 srtBasis,

								 iNPayDates,
                                 lPayDates,
								 dCoverages,

								 dStrikes,

								 dShifts,

								 xlUpBound,
								 nClass,

								 dPPDx,
								 dPPDy,
								 dCPDy
								);	
	if (err)
	{
		goto FREE_RETURN;
	}

//---------------------------------------------------------------------------------
//--------Conversion into shifted Log before Copula--------------------------------
//---------------------------------------------------------------------------------
	dPPDShiftedLogx=dmatrix(0,iNDates-1,0,2*nClass);
	for(i=0;i<iNDates;++i)
	{
		for(j=0;j<2*nClass+1;++j)
		{
			dPPDShiftedLogx[i][j] = log(dShifts[i] + dPPDx[i][j]);
		}
	}

	
//---------------------------------------------------------------------------------
//--------Copula Retrieves the Matrix of Joint distributions-----------------------
//---------------------------------------------------------------------------------
	dMean_v=dvector(0,iNDates-1);
	dAvg=dvector(0,iNDates-1);
	for (i=0; i<iNDates; i++) {
		dAvg[i]=0.0;
		dMean_v[i]=0.0;  
	}

	dCopulaCube = dmatrix(0,xlNumSim,0,iNDates-1);
	dCopulaCubeShiftedLog = dmatrix(0,xlNumSim,0,iNDates-1);

//	GetStudentCplDev (xlStudDegree, dMean_v, Correl, xlNumSim, xlStudDegree+iNDates,dPPDShiftedLogx,
//						  dCPDy,2*nClass+1,0,RANDOM_GAUSS,dCopulaCubeShiftedLog);

	CopulaSimul(Correl,xlNumSim,iNDates,dPPDShiftedLogx,dCPDy,2*nClass+1,dCopulaCubeShiftedLog);

	//--------Get the swap rates back from their shifted log---------------------------
	for(i=0;i<iNDates;++i)
	{
		for(j=0;j<xlNumSim+1;++j)
		{
			dCopulaCube[j][i] = exp(dCopulaCubeShiftedLog[j][i]) - dShifts[i];
		}
	}

	for (i=3; i<xlNumSim; i++)
	{
		for(j=0;j<lNFixNot;++j)
		{
			mean[j] += dCoInitNotionals[j] * dLvls[j] *(dCopulaCube[i][j]-dStrikes[j]);
		}
	}

	for(j=0;j<lNFixNot;++j)
	{
//		err = swp_f_ForwardRate(lStartDates[0], lEndDates[j], cFreq, cBasis, cYCname, cRefRname, &dFwdSwap);
		dFwdSwap = dFwdSwapRate[j];
		dFwdLvl = dLvls[j];
		mean[j] = dCoInitNotionals[j] * dFwdLvl * ( dFwdSwap - dStrikes[j]) - mean[j] / (xlNumSim-2);
	}


	*dPrice=0.0;

	//-------MonteCarlo----------------------------------------------------------------
	for (i=3; i<xlNumSim; i++)
	{
		err = GetAmortizedSwaptionPayoff5_c(dCopulaCube,
											srtFreq,
											srtBasis,
											dCoverages,
											lNFixNot,
											dLvls,
											dStrikes,
											dCoInitNotionals,
											i,
											&dPayoff,
											mean,
											flag);

		*dPrice += (DMAX(0,flag * dPayoff) - *dPrice)/(i-2);
	}


	//-----------------Frees the memory and return-------------------------------------

FREE_RETURN:
	if(dStrikes) free_dvector(dStrikes, 0, iNDates-1);

	if(dShifts) free_dvector(dShifts, 0, iNDates-1);

	if(dPPDx) free_dmatrix(dPPDx,0,iNDates-1,0,2*nClass);
	if(dPPDy) free_dmatrix(dPPDy,0,iNDates-1,0,2*nClass);
	if(dCPDy) free_dmatrix(dCPDy,0,iNDates-1,0,2*nClass);

	if(dPPDShiftedLogx) free_dmatrix(dPPDShiftedLogx,0,iNDates-1,0,2*nClass);

	if(dMean_v) free_dvector(dMean_v,0,iNDates-1);
	if(dAvg) free_dvector(dAvg,0,iNDates-1);
	if(dCopulaCube) free_dmatrix(dCopulaCube,0,xlNumSim,0,iNDates-1);
	if(dCopulaCubeShiftedLog) free_dmatrix(dCopulaCubeShiftedLog,0,xlNumSim,0,iNDates-1);

	if(dCoverages) free(dCoverages);
	if(lStartDates) free(lStartDates);
	if(lEndDates) free(lEndDates);
	if(lPayDates) free(lPayDates);
	
	if(lCoInitSwapsEndDates) free(lCoInitSwapsEndDates);
	if(dCoInitStrikes) free(dCoInitStrikes);
	
	if(dCoInitNotionals) free_dvector(dCoInitNotionals, 0, iNDates-1);

	if(mean) free_dvector(mean, 0,lNFixNot-1);
	if(dLvls) free_dvector(dLvls, 0,lNFixNot-1);
	if(dFwdSwapRate) free_dvector(dFwdSwapRate, 0,lNFixNot-1);

	if(lCpnDates) free(lCpnDates);
	if(dCpns) free(dCpns);

	if(dSpreads) free(dSpreads);

	return err;
}
//--------------------------------------------------------------------------------------
//----------------Shifted Log Model : End Of Old Function------------------------------
//--------------------------------------------------------------------------------------
*/

/*
//--------------------------------------------------------------------------------------
//----------- Convert Amortized Swap in CoInitial Swaps Portfolio (2) ------------------
//--------------------------------------------------------------------------------------
Err		ConvertAmortSwapInCoInitalSwapsPortfolio2(
								 char *cYCname,
								 char *cRefRname,
								 long StartDate,
								 long EndDate,
								 SrtCompounding srtFreq,
								 SrtBasisCode srtBasis,
								 long	lNFixNot,
								 double *dFixNotionals,
								 double *dFixRates,
								 double  *dStrikes,
								 long	*lNCoInitSwaps,
								 double  **dAmounts)
{
	Err err = NULL;
	int i;
	SrtCurvePtr	pCurve = lookup_curve (cYCname);
	SwapDP Swap;
	Date lToday;

	double sum;

	long iNFixPayDates, iNFixDates;
	long *lFixPayDates=NULL, *lFixStartDates=NULL, *lFixEndDates=NULL;
	double *dFixCoverages=NULL;

	pCurve = lookup_curve (cYCname);
	lToday = (Date) get_today_from_curve (pCurve);

    err = swp_f_setSwapDP(StartDate, EndDate, srtFreq, srtBasis, &Swap);
	if(err) goto FREE_RETURN;
														
	err = swp_f_make_FixedLegDatesAndCoverages(&Swap,lToday,&lFixPayDates,&iNFixPayDates,
								  			   &lFixStartDates,&lFixEndDates,
								  			   &dFixCoverages,&iNFixDates);
	if(err) goto FREE_RETURN;

	if(iNFixDates!=lNFixNot)
	{
		err = "Amortized Notionals : Wrong Dimensions";
		goto FREE_RETURN;
	}

	*dAmounts = (double*) calloc (iNFixDates, sizeof(double));

	(*dAmounts)[iNFixDates-1] = dFixNotionals[iNFixDates-1] 
								* (1+dFixCoverages[iNFixDates-1]*dFixRates[iNFixDates-1])
								/ (1+dFixCoverages[iNFixDates-1]*dStrikes[iNFixDates-1]);

	sum = dStrikes[iNFixDates-1] * (*dAmounts)[iNFixDates-1];
	for(i=1;i<iNFixDates;++i)
	{
		(*dAmounts)[iNFixDates-i-1] = ( dFixNotionals[iNFixDates-i-1] - dFixNotionals[iNFixDates-i] 
										+ dFixNotionals[iNFixDates-i-1] 
											* dFixCoverages[iNFixDates-i-1] 
											* dFixRates[iNFixDates-i-1] 
										- dFixCoverages[iNFixDates-i-1] * sum )
										/ (1 + dFixCoverages[iNFixDates-i-1] * dStrikes[iNFixDates-i-1] );
		sum += dStrikes[iNFixDates-i-1] * (*dAmounts)[iNFixDates-i-1];
	}

	*lNCoInitSwaps = iNFixDates;

FREE_RETURN :

	if(lFixPayDates) free(lFixPayDates);
	if(lFixStartDates) free(lFixStartDates);
	if(lFixEndDates) free(lFixEndDates);
	if(dFixCoverages) free(dFixCoverages);

	return err;
}

//--------------------------------------------------------------------------------------
//------------End Of Convert Amortized Swap in CoInitial Swaps Portfolio (2)------------
//--------------------------------------------------------------------------------------
*/

/*
double PV_IRR(char *cYCname, long lToday, int iIndexStart, int iNPayDates, long *lPayDates, double *dCpns, double R)
{
	int i;
	double PV;

	PV=0.0;
	for(i=iIndexStart;i<iNPayDates;++i)
	{
		PV+= dCpns[i]
			* (swp_f_df(lToday, lPayDates[i], cYCname) / swp_f_df(lToday, lPayDates[0], cYCname) )
			* exp(-R*(lPayDates[i]-lPayDates[0])/365.0);
	}
	return PV;
}


Err		ComputeIRRStrikes(char *cYCname,
						  long lToday,
						 char *cRefRname,
						 long StartDate,
						 long EndDate,
						 SrtCompounding srtFreq,
						 SrtBasisCode srtBasis,
						 long	lNFixNot,
						 double *dFixNotionals,
						 long	lNFloatNot,
						 double *dFloatNotionals,
						 double *dFixRates,
						 double *dStrikes)
{
	SwapDP Swap;
	long iNFixPayDates, iNFixDates;
	long *lFixPayDates=NULL, *lFixStartDates=NULL, *lFixEndDates=NULL;
	double *dFixCoverages=NULL;

	long iNFloatPayDates, iNFloatDates;
	long *lFloatFixingDates=NULL, *lFloatPayDates=NULL, *lFloatStartDates=NULL, *lFloatEndDates=NULL;
	double *dFloatCoverages=NULL;
	double *dFloatSpreads=NULL;
	int i, FixFloatMult;
	double *dCpns=NULL;
	double *dCpns_Fix=NULL;
	double *dCpns_Float=NULL;
	long *lCpnDates=NULL;
	double R, R2;
	double PV, PV2;
	double PV_Float, PV_Fix;

	int Compt;

	Err err=NULL;

    err = swp_f_setSwapDP(StartDate, EndDate, srtFreq, srtBasis, &Swap);
	if(err) goto FREE_RETURN;
														
	err = swp_f_make_FixedLegDatesAndCoverages(&Swap,lToday,&lFixPayDates,&iNFixPayDates,
								  			   &lFixStartDates,&lFixEndDates,
								  			   &dFixCoverages,&iNFixDates);
	if(err) goto FREE_RETURN;


	err = swp_f_make_FloatLegDatesCoveragesAndSpreads(&Swap, lToday, cRefRname, 
												&lFloatPayDates,
												&iNFloatPayDates, 
												&lFloatFixingDates, 
												&lFloatStartDates,
												&lFloatEndDates,
												&dFloatCoverages, 
												&dFloatSpreads, 
												&iNFloatDates);
	if(err) goto FREE_RETURN;

	if((iNFloatDates!=lNFloatNot)||(iNFixDates!=lNFixNot))
	{
		err = "Amortized Notionals : Wrong Dimensions";
		goto FREE_RETURN;
	}

	dCpns = dvector(0,iNFloatPayDates-1);
	lCpnDates = lvector(0,iNFloatPayDates-1);
	
	dCpns_Float = dvector(0,iNFloatPayDates-1);
	dCpns_Fix = dvector(0,iNFixPayDates-1);

	FixFloatMult = (int)((double)(iNFloatDates)/iNFixDates+0.5);

	lCpnDates[0] = lFloatPayDates[0];
	dCpns[0] = -dFloatNotionals[0];		
	dCpns_Float[0] = -1.0;
	dCpns_Fix[0] = 0.;
	for(i=1;i<iNFloatPayDates-1;++i)
	{
		lCpnDates[i] = lFloatPayDates[i];

		dCpns[i] = - (dFloatNotionals[i] - dFloatNotionals[i-1]
						+ dFloatNotionals[i-1] 
							* dFloatCoverages[i-1]
							* dFloatSpreads[i-1]);

		dCpns_Float[i] = - dFloatCoverages[i-1] * dFloatSpreads[i-1];

		if(lFixPayDates[(int)(i/FixFloatMult)]==lFloatPayDates[i])
		{
			dCpns[i] = dCpns[i] + dFixCoverages[(int)(i/FixFloatMult)-1]*dFixRates[(int)(i/FixFloatMult)-1] * dFixNotionals[(int)(i/FixFloatMult)-1];

			dCpns_Fix[(int)(i/FixFloatMult)] = dFixCoverages[(int)(i/FixFloatMult)-1];
		}
	}

	lCpnDates[iNFloatPayDates-1] = lFloatPayDates[iNFloatPayDates-1];

	dCpns[iNFloatPayDates-1] = dFloatNotionals[iNFloatPayDates-2]
				- dFloatNotionals[iNFloatPayDates-2]
					* dFloatCoverages[iNFloatPayDates-2]
					* dFloatSpreads[iNFloatPayDates-2];

	dCpns[iNFloatPayDates-1] = dCpns[iNFloatPayDates-1] 
							+ dFixRates[iNFixPayDates-2] 
								* dFixCoverages[iNFixPayDates-2]
								* dFixNotionals[iNFixPayDates-2];

	dCpns_Float[iNFloatPayDates-1] = - dFloatCoverages[iNFloatPayDates-2]
										* dFloatSpreads[iNFloatPayDates-2];

	dCpns_Fix[iNFixPayDates-1] = dFixCoverages[iNFixPayDates-2];

	R = 0.05;
	PV = PV_IRR(cYCname, lToday, 0, iNFloatPayDates, lFloatPayDates, dCpns, R);

	Compt = 0;
	while((fabs(PV)>1e-12)&&(Compt<15))
	{
		R2 = R + 0.0001;
		PV2 = PV_IRR(cYCname, lToday, 0, iNFloatPayDates, lFloatPayDates, dCpns, R2);

		R = R - PV * (R2 - R) / (PV2 - PV);
		PV = PV_IRR(cYCname, lToday, 0, iNFloatPayDates, lFloatPayDates, dCpns, R);

		Compt = Compt + 1;
	}

	for(i=1;i<iNFixPayDates;++i)
	{
		dCpns_Float[i*FixFloatMult] = dCpns_Float[i*FixFloatMult] + 1;

		PV_Float = PV_IRR(cYCname, lToday, 0, i*FixFloatMult+1, lFloatPayDates, dCpns_Float, R);
		PV_Fix = PV_IRR(cYCname, lToday, 0, i+1, lFixPayDates, dCpns_Fix, R);

		dStrikes[i-1] = -PV_IRR(cYCname, lToday, 0, i*FixFloatMult+1, lFloatPayDates, dCpns_Float, R);
		dStrikes[i-1] = dStrikes[i-1] / PV_IRR(cYCname, lToday, 0, i+1, lFixPayDates, dCpns_Fix, R);

		dCpns_Float[i*FixFloatMult] = dCpns_Float[i*FixFloatMult] - 1;
	}


FREE_RETURN :

	if(dCpns) free_dvector(dCpns, 0,iNFloatPayDates-1);
	if(lCpnDates) free_lvector(lCpnDates, 0,iNFloatPayDates-1);

	if(dCpns_Float) free_dvector(dCpns_Float, 0,iNFloatPayDates-1);
	if(dCpns_Fix) free_dvector(dCpns_Fix, 0,iNFixPayDates-1);

	if(lFixPayDates) free(lFixPayDates);
	if(lFixStartDates) free(lFixStartDates);
	if(lFixEndDates) free(lFixEndDates);
	if(dFixCoverages) free(dFixCoverages);

	if(lFloatPayDates) free(lFloatPayDates);
	if(lFloatFixingDates) free(lFloatFixingDates);
	if(lFloatStartDates) free(lFloatStartDates);
	if(lFloatEndDates) free(lFloatEndDates);
	if(dFloatCoverages) free(dFloatCoverages);
	if(dFloatSpreads) free(dFloatSpreads);
	
	return err;

}
*/

/*
//--------------------------------------------------------------------------------------
//---------------------Shifted Log Model : Main Function--------------------------------
//--------------------------------------------------------------------------------------
Err		PriceAmortizedSwaptionShiftedLog(
						 char *cYCname,
						 char *cVCname,
						 char *cRefRname,

						 long xlStartDate,
						 long xlEndDate,

						 SrtCompounding srtFreq,
						 SrtBasisCode srtBasis,
						 
						 long	lNFixNot,
						 double *dFixNotionals,
						 long	lNFloatNot,
						 double *dFloatNotionals,
						 double *dFixRates,

						 SrtCallPutType srtCallPut,

						 int xlStudDegree,
						 int xlNumSim,
						 double xlUpBound,
						 int xlnClasses,
						 double **Correl,
						 double *dPrice,
						 int iPayoffMethod,
						 int iDistribMethod
			            )
{
	
	///// Objects and variables declaration
	Err err = NULL;

	SrtCurvePtr	pCurve = lookup_curve (cYCname);

	SwapDP Swap;

	Date lToday,lSpotLag;

	Date *lPayDates=NULL,*lStartDates=NULL,*lEndDates=NULL;
	double *dCoverages=NULL;
	int iNDates, iNPayDates;

	double **dCopulaCube=NULL,**dCopulaCubeShiftedLog=NULL,*dMean_v=NULL,*dAvg=NULL;

	double **dPPDx=NULL,**dPPDy=NULL,**dCPDy=NULL;
	double **dPPDShiftedLogx=NULL;

	int i, j;
	
	int nClass=1<<xlnClasses;
	
	int flag;

	double dPayoff;
	
//	int iNCpns;
	Date *lCpnDates=NULL;
	double *dCpns=NULL;

	int iNCoInitSwaps;
	Date *lCoInitSwapsEndDates=NULL;
	double *dCoInitStrikes=NULL;
	double *dCoInitNotionals=NULL;

//	int NSpreads;
	double *dSpreads=NULL;

	double *dShifts=NULL;

	double *dStrikes=NULL;

//	double dSwapPV;
	double df_start;

	char *cBasis, *cFreq;

	//---Variable for variate control--------------------------------------------------
	double *mean=NULL;
	double *dLvls=NULL;
	double dCov, dFwdSwap, dFwdLvl;


	err = translate_basis(&cBasis, srtBasis);
	err = translate_compounding(&cFreq, srtFreq);
	/////////////////////////////////////////////////////////////////////////
	///// defines the pointer to the curve and retrieves today and spot dates
	/////////////////////////////////////////////////////////////////////////

	pCurve = lookup_curve (cYCname);
	///// if (!pCurve) throw XSort("Fatal: (FwdSwaption) Yield Curve Object not found");
	lToday = (Date) get_today_from_curve (pCurve);
	lSpotLag = (Date) get_spotlag_from_curve(pCurve);	
	
	///////////////////////////////////////////////////////////////
	///// defines the relevant swap objects
	///////////////////////////////////////////////////////////////
    err = swp_f_setSwapDP(xlStartDate, xlEndDate, srtFreq, srtBasis, &Swap);
	if (err)
	{
		goto FREE_RETURN;
	}
	
	err = swp_f_make_FixedLegDatesAndCoverages(&Swap,lToday,&lPayDates,&iNPayDates,
								  			   &lStartDates,&lEndDates,
								  			   &dCoverages,&iNDates);
	if (err)
	{
		goto FREE_RETURN;
	}

//-----------------------------------------------------------------------
//------------------------Compute Calibration Strikes--------------------
//-----------------------------------------------------------------------
	dStrikes = dvector(0,iNDates-1);
	err = ComputeIRRStrikes(cYCname, lToday, cRefRname,
						 xlStartDate, xlEndDate,
						 srtFreq,srtBasis,
						 lNFixNot,dFixNotionals,
						 lNFloatNot,dFloatNotionals,dFixRates,
						 dStrikes);
	if (err)
	{
		goto FREE_RETURN;
	}
//-----------------------------------------------------------------------
//---------------------End Of Compute Calibration Strikes----------------
//-----------------------------------------------------------------------


	//-----------------Retrieves the densities ----------------------------------------
	dShifts = dvector(0,iNDates-1);
	dPPDx=dmatrix(0,iNDates-1,0,2*nClass);
	dPPDy=dmatrix(0,iNDates-1,0,2*nClass);
	dCPDy=dmatrix(0,iNDates-1,0,2*nClass);
	err = Fwd_Swaption_Get_Cum_ShiftedLogDensity(
								 lToday,
								 cYCname,
								 cVCname,
								 cRefRname,
								 srtFreq,
								 srtBasis,

								 iNPayDates,
                                 lPayDates,
								 dCoverages,

								 dStrikes,

								 dShifts,

								 xlUpBound,
								 nClass,

								 dPPDx,
								 dPPDy,
								 dCPDy
								);	
	if (err)
	{
		goto FREE_RETURN;
	}

	//--------Conversion into shifted Log befor Copula---------------------------------
	dPPDShiftedLogx=dmatrix(0,iNDates-1,0,2*nClass);
	for(i=0;i<iNDates;++i)
	{
		for(j=0;j<2*nClass+1;++j)
		{
			dPPDShiftedLogx[i][j] = log(dShifts[i] + dPPDx[i][j]);
		}
	}

	
	//--------Copula Retrieves the Matrix of Joint distributions-----------------------
	dMean_v=dvector(0,iNDates-1);
	dAvg=dvector(0,iNDates-1);
	for (i=0; i<iNDates; i++) {
		dAvg[i]=0.0;
		dMean_v[i]=0.0;  
	}

	dCopulaCube = dmatrix(0,xlNumSim,0,iNDates-1);
	dCopulaCubeShiftedLog = dmatrix(0,xlNumSim,0,iNDates-1);

//	GetStudentCplDev (xlStudDegree,dMean_v,Correl,xlNumSim,xlStudDegree+iNDates,dPPDShiftedLogx,
//						  dCPDy,2*nClass+1,0,RANDOM_GAUSS,dCopulaCubeShiftedLog);

	CopulaSimul(Correl,xlNumSim,iNDates,dPPDShiftedLogx,dCPDy,2*nClass+1,dCopulaCubeShiftedLog);

	//--------Get the swap rates back from their shifted log---------------------------
	for(i=0;i<iNDates;++i)
	{
		for(j=0;j<xlNumSim+1;++j)
		{
			dCopulaCube[j][i] = exp(dCopulaCubeShiftedLog[j][i]) - dShifts[i];
		}
	}


	//--------Convert AMort Swap in Coinit Swap Portfolio------------------------------
	err=ConvertAmortSwapInCoInitalSwapsPortfolio2(
							cYCname,
							cRefRname,
							xlStartDate,
							xlEndDate,
							srtFreq,
							srtBasis,
							lNFixNot,
							dFixNotionals,
							dFixRates,
							dStrikes,
							&iNCoInitSwaps,
//							&lCoInitSwapsEndDates,
//							&dCoInitStrikes,
							&dCoInitNotionals);
	if (err)
	{
		goto FREE_RETURN;
	}

	*dPrice=0.0;

	//-------MonteCarlo----------------------------------------------------------------
	if (cFreq == "ANNUAL") {
		dCov=1.0;
	} else if (cFreq == "SEMIANNUAL") {
		dCov=0.5;
	} else if (cFreq == "QUARTERLY") {
		dCov=0.25;
	} else if (cFreq == "MONTHLY") {
		dCov=1./12.;
	}

	mean = dvector(0,iNCoInitSwaps-1);
	dLvls = dvector(0,iNCoInitSwaps-1);
	for(j=0;j<iNCoInitSwaps;++j)
	{
		mean[j] = 0.0;
		err = swp_f_LevelPayment(lStartDates[0], lEndDates[j], cFreq, cBasis, cYCname, cRefRname, &(dLvls[j]));
	}

	for (i=3; i<xlNumSim; i++)
	{
		for(j=0;j<iNCoInitSwaps;++j)
		{
//			err = swp_f_LevelPayment(lStartDates[0], lEndDates[j], cFreq, cBasis, cYCname, cRefRname, &dFwdLvl);
//			dLevelCash =(1.-1./pow(1+dCopulaCube[i][j]*dCov,j+1))/dCopulaCube[i][j];
			mean[j] += dCoInitNotionals[j] * dLvls[j] *(dCopulaCube[i][j]-dStrikes[j]);
		}
	}

	for(j=0;j<iNCoInitSwaps;++j)
	{
		err = swp_f_ForwardRate(lStartDates[0], lEndDates[j], cFreq, cBasis, cYCname, cRefRname, &dFwdSwap);
//		dFwdLvl = (1.-1./pow(1+dFwdSwap*dCov,j+1))/dFwdSwap;
		dFwdLvl = dLvls[j];
//		err = swp_f_LevelPayment(lStartDates[0], lEndDates[j], cFreq, cBasis, cYCname, cRefRname, &dFwdLvl);
		mean[j] = dCoInitNotionals[j] * dFwdLvl * ( dFwdSwap - dStrikes[j]) - mean[j] / (xlNumSim-2);
	}

	
//	dSwapPV = amortSwap_PV(lToday, cYCname, iNCpns, lCpnDates, dCpns);
	df_start = swp_f_df(lToday, xlStartDate, cYCname);
	flag = -1;
    if (srtCallPut == SRT_CALL) flag=1;
	for (i=3; i<xlNumSim; i++)
	{
//		if(iPayoffMethod==5)
//		{

			err = GetAmortizedSwaptionPayoff5_c(dCopulaCube,
											srtFreq,
											srtBasis,
											dCoverages,
											iNCoInitSwaps,
											dLvls,
											dStrikes,
											dCoInitNotionals,
											i,
											&dPayoff,
											mean,
											flag);
		*dPrice += (DMAX(0,flag * dPayoff) - *dPrice)/(i-2);
//		*dPrice += (flag * dPayoff - *dPrice)/(i-2);
	}

//	*dPrice = df_start* (*dPrice);

	//-----------------Frees the memory and return-------------------------------------

FREE_RETURN:
	if(dStrikes) free_dvector(dStrikes, 0, iNDates-1);

	if(dShifts) free_dvector(dShifts, 0, iNDates-1);

	if(dPPDx) free_dmatrix(dPPDx,0,iNDates-1,0,2*nClass);
	if(dPPDy) free_dmatrix(dPPDy,0,iNDates-1,0,2*nClass);
	if(dCPDy) free_dmatrix(dCPDy,0,iNDates-1,0,2*nClass);

	if(dPPDShiftedLogx) free_dmatrix(dPPDShiftedLogx,0,iNDates-1,0,2*nClass);

	if(dMean_v) free_dvector(dMean_v,0,iNDates-1);
	if(dAvg) free_dvector(dAvg,0,iNDates-1);
	if(dCopulaCube) free_dmatrix(dCopulaCube,0,xlNumSim,0,iNDates-1);
	if(dCopulaCubeShiftedLog) free_dmatrix(dCopulaCubeShiftedLog,0,xlNumSim,0,iNDates-1);

	if(dCoverages) free(dCoverages);
	if(lStartDates) free(lStartDates);
	if(lEndDates) free(lEndDates);
	if(lPayDates) free(lPayDates);
	
	if(lCoInitSwapsEndDates) free(lCoInitSwapsEndDates);
	if(dCoInitStrikes) free(dCoInitStrikes);
	if(dCoInitNotionals) free(dCoInitNotionals);

	if(lCpnDates) free(lCpnDates);
	if(dCpns) free(dCpns);

	if(dSpreads) free(dSpreads);

	return err;
}
//--------------------------------------------------------------------------------------
//----------------Shifted Log Model : End Of Main Function------------------------------
//--------------------------------------------------------------------------------------
*/



/*
//-----------------------------------------------------------------------------------
//        Function that transforms smile from the swap measure into Fwd measure   
//-----------------------------------------------------------------------------------

/*

Err  GetDensity_Lvl_Measure9_b(
							int    nClass,
							double dMat,
							double dFwd,
							double dSigmaH,
							double dAlphaH,
							double dSigmaIftyH,
							double dMeanRevH,
							double dBetaH,
							double dRhoH,
							double dDisc,
							double dUpBound,
							int    nSteps,
							int    nTen,
							double dCov,
							double **dPPDx,
							double **dPPDy,
							double **dCPDy
							)
{

	double dEqBSVol, dStdDev, LvlC_x, LvlCFwd;
	double *dVec = NULL, *dX=NULL,*dY=NULL, dThreshold=1.e-3;
	int i;
	double dShiftH;
	double deriv1, deriv2;
	Err err = NULL;
	double epsilon = 1.0e-4;
	double *strike=NULL, *strike_plus=NULL, *strike_moins=NULL;
	double *temp=NULL, *temp_plus=NULL, *temp_moins=NULL;
	double dK;
	double price_moins;
	double price;
	double price_plus;
	double *strikeVect=NULL;
	double *priceVect=NULL;

//	dShiftH = HESTON_CONST*(1.-dBetaH)/dBetaH;
	dShiftH = dFwd*(1.-dBetaH)/dBetaH;

	dEqBSVol=(dFwd+dShiftH)*dSigmaH;
	dStdDev=10.;

	strikeVect=dvector(1,2*nClass+1);
	priceVect=dvector(1,2*nClass+1);

	strike=dvector(1,1);
	strike_plus=dvector(1,1);
	strike_moins=dvector(1,1);

	temp=dvector(1,1);
	temp_plus=dvector(1,1);
	temp_moins=dvector(1,1);

	///// gets the two boundaries

	dPPDx[nTen][0] = -dShiftH + 0.0001;

	dPPDx[nTen][2*nClass] = 1.0;//sqrt(dMat);//dFwd + dStdDev*dEqBSVol*sqrt(dMat);

	dK = (dPPDx[nTen][2*nClass]-dPPDx[nTen][0])/(2*nClass);

	for (i=1;i<2*nClass;i++)
	{
		dPPDx[nTen][i]=dPPDx[nTen][0] + i * dK;
	}

	///// runs through the points 
	for (i=0;i<=2*nClass;i++)
	{
		strikeVect[i+1] = dPPDx[nTen][i];
		dCPDy[nTen][i]=0.0;
		dPPDy[nTen][i]=0.0;
	}

	HestonPrice(dFwd,strikeVect,2*nClass+1,dMat,
		dSigmaH,dAlphaH,dSigmaIftyH,dMeanRevH,dBetaH,dRhoH,
		1.0,dUpBound,SRT_CALL,PREMIUM,SRT_TRUE,2,nSteps,priceVect);

	////i=0
	price = priceVect[1];
	price_plus = priceVect[2];

	LvlC_x = (1.-1./pow(1+dPPDx[nTen][0]*dCov,nTen+1))/dPPDx[nTen][0];
	LvlCFwd = (1.-1./pow(1+dFwd*dCov,nTen+1))/dFwd;

	deriv1 = 1.0;
	deriv2 = DMIN(1.0,(price-price_plus)/dK);

	dPPDy[nTen][0]=DMAX(0,(deriv1 - deriv2))/dK;///(dK*LvlC_x));

	dCPDy[nTen][0]=0.0;
	
	for (i=1;i<2*nClass;i++)
	{
		price_moins = price;
		price = price_plus;
		price_plus = priceVect[i+2];

		LvlC_x = (1.-1./pow(1+dPPDx[nTen][i]*dCov,nTen+1))/dPPDx[nTen][i];

		deriv1 = DMIN(1.0,(price_moins - price)/dK);
		deriv2 = DMIN(1.0,(price - price_plus)/dK);

		dPPDy[nTen][i]=LvlCFwd * DMAX(0,(deriv1 - deriv2))/dK;//(dK*LvlC_x));

		dCPDy[nTen][i]=dCPDy[nTen][i-1]
			+(dPPDy[nTen][i]+dPPDy[nTen][i-1])/2.*dK;
	}

	price_moins = price;
	price = price_plus;

	LvlC_x = (1.-1./pow(1+dPPDx[nTen][2*nClass]*dCov,nTen+1))/dPPDx[nTen][2*nClass];

	deriv1 = 0.0;
	deriv2 = (price - price_moins)/dK;

	dPPDy[nTen][2*nClass]=DMAX(0,(deriv1 - deriv2))/dK;///(dK*LvlC_x));

	dCPDy[nTen][2*nClass]=dCPDy[nTen][2*nClass-1]
		+(dPPDy[nTen][2*nClass]+dPPDy[nTen][2*nClass-1])/2.*dK;

	for (i=0;i<=2*nClass;i++) 
	{
		dCPDy[nTen][i]/=dCPDy[nTen][2*nClass];
	}

	return err;

}
*/

/*
Err       Fwd_Swaption_Get_Cum_HestonDensity(Date lToday,
								 char *cYCname,
						         char *cRefRname,
								 SrtCompounding srtFreq,
								 SrtBasisCode srtBasis,

								 int iNPayDatesLong,
								 Date *lPayDatesLong,
								 double *dCoveragesLong,

								 double *SigmaHeston,
								 double *AlphaHeston,
								 double *LambdaHeston,
								 double *BetaHeston,
								 double *RhoHeston,

								 int iNPts,
								 double dUpBound,
								 int nClass,

								 double **dPPDx,
                                 double **dPPDy,
							     double **dCPDy,

								 int iDistribMethod
								 )

{
	Err err = NULL;
	int i;
	double dMat,dFwd,dLevel,dCov;
	char *cBasis, *cFreq;

	dMat=(lPayDatesLong[0]-lToday-2)/365.; 
	
//////////////////////////////////////////
/// translate SORT types into strings
//////////////////////////////////////////

	err = translate_basis(&cBasis, srtBasis);
	err = translate_compounding(&cFreq, srtFreq);

	if (cFreq == "ANNUAL") {
		dCov=1.0;
	} else if (cFreq == "SEMIANNUAL") {
		dCov=0.5;
	} else if (cFreq == "QUARTERLY") {
		dCov=0.25;
	} else if (cFreq == "MONTHLY") {
		dCov=1./12.;
	}

//////////////////////////////////////////////////////////
////////////////// generates the cumulative densities
//////////////////////////////////////////////////////////	

	for (i=1;i<iNPayDatesLong;i++) {

//      computes the Fwd swap rate

		err = swp_f_ForwardRate(lPayDatesLong[0],lPayDatesLong[i],cFreq,cBasis,cYCname,cRefRname,&dFwd);
		err = swp_f_LevelPayment(lPayDatesLong[0],lPayDatesLong[i],cFreq,cBasis,cYCname,cRefRname,&dLevel);

		err = GetDensity_Lvl_Measure9_b(nClass,dMat,dFwd,SigmaHeston[i-1],AlphaHeston[i-1],1.0,LambdaHeston[i-1],
			                 BetaHeston[i-1],RhoHeston[i-1],dLevel,dUpBound,iNPts,i-1,dCov,dPPDx,dPPDy,dCPDy);
	}

	return err;

}	
*/

/*

//--------------------------------------------------------------------------------------
//-------------------Utilities Functions for Heston Model-------------------------------
//--------------------------------------------------------------------------------------
Err  GetHestonParameters(long lToday,
						 char	*cYCname,
						 char	*cVCname,
						 long lStartDate, long lEndDate, 
						 SrtCompounding srtFreq,
						 SrtBasisCode srtBasis,
						 char	*cRefRname,
						 double UpperBound,
						 int nSteps,
						 int iIntegerType,
						 double nStdDev,
						 double **HestonSigma,
						 double **HestonAlpha,
						 double **HestonBeta,
						 double **HestonRho,
						 double **HestonLambda)
{
	//---Swap Details------------------------------------------------------------------
	SwapDP Swap;
	char *cBasis, *cFreq;
	Date *lPayDates=NULL,*lStartDates=NULL,*lEndDates=NULL;
	double *dCoverages=NULL;
	int iNDates, iNPayDates;

	int i;
	double SABRsigma, SABRalpha, SABRbeta, SABRrho;
	double HESTONsigma, HESTONalpha, HESTONbeta, HESTONrho;
	double HESTONsigmainfinity, HESTONlambda;

	long start, end;

	double dFwdSwap, dFwdCash, dSpread, dLvl, ex_time;
	double power;
//	double shift;

	Err err=NULL;

	//----------Convert Basis and Freq in SRTcode--------------------------------------
	err = translate_basis(&cBasis, srtBasis);
	err = translate_compounding(&cFreq, srtFreq);


	//-----Init and Fill SWAPDP--------------------------------------------------------
    err = swp_f_setSwapDP(lStartDate, lEndDate, srtFreq, srtBasis, &Swap);
	if(err)
	{
		goto FREE_RETURN;
	}
	err = swp_f_make_FixedLegDatesAndCoverages(&Swap,lToday,&lPayDates,&iNPayDates,
								  			   &lStartDates,&lEndDates,
								  			   &dCoverages,&iNDates);	
	if(err)
	{
		goto FREE_RETURN;
	}

	*HestonSigma = (double*) calloc (iNDates, sizeof (double));
	*HestonAlpha = (double*) calloc (iNDates, sizeof (double));
	*HestonBeta = (double*) calloc (iNDates, sizeof (double));
	*HestonRho = (double*) calloc (iNDates, sizeof (double));
	*HestonLambda = (double*) calloc (iNDates, sizeof (double));

	for(i=0;i<iNDates;++i)
	{
		start = lStartDate;
		end = lEndDates[i];

		err = swp_f_ForwardRate(start, end, cFreq, cBasis, cYCname, cRefRname, &dFwdSwap);
		dFwdCash = swp_f_swapcashrate( start, end, srtBasis, srtFreq, cYCname, cRefRname);
		err = swp_f_LevelPayment(start, end, cFreq, cBasis, cYCname, cRefRname, &dLvl);
		ex_time = (add_unit (start, - 2, SRT_BDAY, MODIFIED_SUCCEEDING) - lToday) * YEARS_IN_DAY;

		dSpread = dFwdSwap - dFwdCash;

		err = swp_f_SABRvol( cVCname, start, end, 0.05,
				&SABRsigma, &power, SABR_BETAVOL );
		if (err) return err;
		err = swp_f_SABRvol( cVCname, start, end, 0.05,
				&SABRalpha, &power, SABR_ALPHA );
		if (err) return err;
		err = swp_f_SABRvol( cVCname, start, end, 0.05,
				&SABRbeta, &power, SABR_BETA );
		if (err) return err;
		err = swp_f_SABRvol( cVCname, start, end, 0.05,
				&SABRrho, &power, SABR_RHO );
		if (err) return err;

		HESTONsigma = SABRsigma;
		HESTONsigmainfinity = SABRsigma;
		HESTONalpha = SABRalpha;
		HESTONbeta = SABRbeta;
		HESTONrho = SABRrho;
		HESTONlambda = 0.05;

		err = CalibrateHestonToSABR(dFwdSwap,
									ex_time,
									SABRsigma,
									SABRalpha,
									SABRbeta,
									SABRrho,
									&HESTONsigma,
									&HESTONalpha,
									&HESTONsigmainfinity,
									&HESTONlambda,
									&HESTONbeta,
									&HESTONrho,
									UpperBound,
									nSteps,
									iIntegerType,
									nStdDev,
									SRT_TRUE,
									CHI2_MIN);

		(*HestonSigma)[i] = HESTONsigma;	
		(*HestonAlpha)[i] = HESTONalpha;	
		(*HestonBeta)[i] = HESTONbeta;	
		(*HestonRho)[i] = HESTONrho;	
		(*HestonLambda)[i] = HESTONlambda;	

	}

FREE_RETURN:

	if(dCoverages) free(dCoverages);
	if(lStartDates) free(lStartDates);
	if(lEndDates) free(lEndDates);
	if(lPayDates) free(lPayDates);


	return err;
	
}

*/

//--------------------------------------------------------------------------------------
//------------End Of Utilities Functions for Heston Model-------------------------------
//--------------------------------------------------------------------------------------


/*
//--------------------------------------------------------------------------------------
//--------------------------Heston Model : Main Function--------------------------------
//--------------------------------------------------------------------------------------
Err		PriceAmortizedSwaptionHeston(
						 char *cYCname,
						 char *cVCname,
						 char *cRefRname,
						 long xlStartDate,
						 long xlEndDate,
						 SrtCompounding srtFreq,
						 SrtBasisCode srtBasis,
						 long	lNFixNot,
						 double *dFixNotionals,
						 long	lNFloatNot,
						 double *dFloatNotionals,
						 double *dFixRates,
						 SrtCallPutType srtCallPut,
						 double xlnStdDev,
						 int xlStudDegree,
						 int xlNumSim,
						 char *xlTypeVol,
						 int xlNPts,
						 double xlUpBound,
						 int xlnClasses,
						 double **Correl,
						 double *dPrice,
						 int iPayoffMethod,
						 int iDistribMethod
			            )
{
	//----------Declarations-----------------------------------------------------------
	//---Err containing result---------------------------------------------------------
	Err err = NULL;
	
	//---Curve-------------------------------------------------------------------------
	SrtCurvePtr	pCurve = lookup_curve (cYCname);
	
	//---SwapDP for schedule generation------------------------------------------------
	SwapDP Swap;

	//---For Dates---------------------------------------------------------------------
	Date lToday,lSpotLag;

	//---Swap Details------------------------------------------------------------------
	char *cBasis, *cFreq;
	Date *lPayDates=NULL,*lStartDates=NULL,*lEndDates=NULL;
	double *dCoverages=NULL;
	int iNDates, iNPayDates;

	//---Monte Carlo Paths-------------------------------------------------------------
	double **dCopulaCube=NULL,**dCopulaCubeShiftedLog=NULL,*dMean_v=NULL,*dAvg=NULL;

	//---Distribution Discretization---------------------------------------------------
	double **dPPDx=NULL,**dPPDy=NULL,**dCPDy=NULL;
	double **dPPDShiftedLogx=NULL;
	
	//---Miscelaneous Integers---------------------------------------------------------
	int i, j;
	int nClass=1<<xlnClasses;
	int flag;

	//---Simulated Payoff--------------------------------------------------------------
	double dPayoff;

	//---Amortized Swap as a Bond------------------------------------------------------
//	int iNCpns;
	Date *lCpnDates=NULL;
	double *dCpns=NULL;

	//---Amortized Swap as a Coinit Swaps Portfolio------------------------------------
	int iNCoInitSwaps;
	Date *lCoInitSwapsEndDates=NULL;
	double *dCoInitStrikes=NULL;
	double *dCoInitNotionals=NULL;

	//---Fwd Swaps Spreads-------------------------------------------------------------
//	int NSpreads;
	double *dSpreads=NULL;

	//---Decomposition Strikes---------------------------------------------------------
	double *dStrikes=NULL;

	//---temporary variable for shifts-------------------------------------------------
	double shift;

	//---Heston Parameters-------------------------------------------------------------
	double *SigmaHeston=NULL;
	double *AlphaHeston=NULL;
	double *LambdaHeston=NULL;
	double *BetaHeston=NULL;
	double *RhoHeston=NULL;

	//---Variable for variate control--------------------------------------------------
	double *mean=NULL;
	double *dLvls=NULL;
	double dFwdSwap;
	double dCov;//, dFwdSwap, dFwdLvl;

	//--------End Of Declarations------------------------------------------------------


	//----------Convert Basis and Freq in SRTcode--------------------------------------
	err = translate_basis(&cBasis, srtBasis);
	err = translate_compounding(&cFreq, srtFreq);

	//-------defines the pointer to the curve and retrieves today and spot dates-------
	pCurve = lookup_curve (cYCname);
	lToday = (Date) get_today_from_curve (pCurve);
	lSpotLag = (Date) get_spotlag_from_curve(pCurve);
	
	//-----Init and Fill SWAPDP--------------------------------------------------------
    err = swp_f_setSwapDP(xlStartDate, xlEndDate, srtFreq, srtBasis, &Swap);
	if(err)
	{
		goto FREE_RETURN;
	}
	err = swp_f_make_FixedLegDatesAndCoverages(&Swap,lToday,&lPayDates,&iNPayDates,
								  			   &lStartDates,&lEndDates,
								  			   &dCoverages,&iNDates);	
	if(err)
	{
		goto FREE_RETURN;
	}

	err = GetHestonParameters(lToday,
								cYCname,
								cVCname,
								xlStartDate, 
								xlEndDate, 
								srtFreq,
								srtBasis,
								cRefRname,
								xlUpBound,
								xlNPts,
								2,
								xlnStdDev,
								&SigmaHeston,
								&AlphaHeston,
								&BetaHeston,
								&RhoHeston,
								&LambdaHeston);


	//--------Retrieves the densities--------------------------------------------------
	dPPDx=dmatrix(0,iNDates-1,0,2*nClass);
	dPPDy=dmatrix(0,iNDates-1,0,2*nClass);
	dCPDy=dmatrix(0,iNDates-1,0,2*nClass);
	err = Fwd_Swaption_Get_Cum_HestonDensity(
								 lToday,
								 cYCname,
								 cRefRname,
								 srtFreq,
								 srtBasis,

								 iNPayDates,
                                 lPayDates,
								 dCoverages,

								 SigmaHeston,
								 AlphaHeston,
								 LambdaHeston,
								 BetaHeston,
								 RhoHeston,

								 xlNPts,
								 xlUpBound,
								 nClass,

								 dPPDx,
								 dPPDy,
								 dCPDy,

								 iDistribMethod
								);	

	//--------Conversion into shifted Log befor Copula---------------------------------
	dPPDShiftedLogx=dmatrix(0,iNDates-1,0,2*nClass);
	for(i=0;i<iNDates;++i)
	{
		err = swp_f_ForwardRate(lStartDates[0], lEndDates[i], cFreq, cBasis, cYCname, cRefRname, &dFwdSwap);
		for(j=0;j<2*nClass+1;++j)
		{
//			shift = HESTON_CONST*(1.-BetaHeston[i])/BetaHeston[i];
			shift = dFwdSwap*(1.-BetaHeston[i])/BetaHeston[i];
			dPPDShiftedLogx[i][j] = log(shift + dPPDx[i][j]);
		}
	}


	//--------Copula Retrieves the Matrix of Joint distributions-----------------------
	dMean_v=dvector(0,iNDates-1);
	dAvg=dvector(0,iNDates-1);
	for (i=0; i<iNDates; i++) {
		dAvg[i]=0.0;
		dMean_v[i]=0.0;  
	}

	dCopulaCube = dmatrix(0,xlNumSim,0,iNDates-1);
	dCopulaCubeShiftedLog = dmatrix(0,xlNumSim,0,iNDates-1);

//	GetStudentCplDev (xlStudDegree,dMean_v,Correl,xlNumSim,xlStudDegree+iNDates,dPPDShiftedLogx,
//						  dCPDy,2*nClass+1,0,RANDOM_GAUSS,dCopulaCubeShiftedLog);

	CopulaSimul(Correl,xlNumSim,iNDates,dPPDShiftedLogx,dCPDy,2*nClass+1,dCopulaCubeShiftedLog);

	//--------Get the swap rates back from their shifted log---------------------------
	for(i=0;i<iNDates;++i)
	{
		err = swp_f_ForwardRate(lStartDates[0], lEndDates[i], cFreq, cBasis, cYCname, cRefRname, &dFwdSwap);
		for(j=0;j<xlNumSim+1;++j)
		{
//			shift = HESTON_CONST*(1.-BetaHeston[i])/BetaHeston[i];
			shift = dFwdSwap*(1.-BetaHeston[i])/BetaHeston[i];
			dCopulaCube[j][i] = exp(dCopulaCubeShiftedLog[j][i]) - shift;
		}
	}

	//--------Convert AMort Swap in Coinit Swap Portfolio------------------------------
	dStrikes = dvector(0,iNDates-1);
	err = ComputeIRRStrikes(cYCname, lToday, cRefRname,
						 xlStartDate, xlEndDate,
						 srtFreq,srtBasis,
						 lNFixNot,dFixNotionals,
						 lNFloatNot,dFloatNotionals,dFixRates,
						 dStrikes);
	if (err)
	{
		goto FREE_RETURN;
	}
	err=ConvertAmortSwapInCoInitalSwapsPortfolio2(
							cYCname,
							cRefRname,
							xlStartDate,
							xlEndDate,
							srtFreq,
							srtBasis,
							lNFixNot,
							dFixNotionals,
							dFixRates,
							dStrikes,
							&iNCoInitSwaps,
//							&lCoInitSwapsEndDates,
//							&dCoInitStrikes,
							&dCoInitNotionals);
	if (err)
	{
		goto FREE_RETURN;
	}


	//-------MonteCarlo----------------------------------------------------------------
	if (cFreq == "ANNUAL") {
		dCov=1.0;
	} else if (cFreq == "SEMIANNUAL") {
		dCov=0.5;
	} else if (cFreq == "QUARTERLY") {
		dCov=0.25;
	} else if (cFreq == "MONTHLY") {
		dCov=1./12.;
	}

	mean = dvector(0,iNCoInitSwaps-1);
	dLvls = dvector(0,iNCoInitSwaps-1);
	for(j=0;j<iNCoInitSwaps;++j)
	{
		mean[j] = 0.0;
		err = swp_f_LevelPayment(lStartDates[0], lEndDates[j], cFreq, cBasis, cYCname, cRefRname, &(dLvls[j]));
	}

	*dPrice=0.0;
	flag = -1;
    if (srtCallPut == SRT_CALL) flag=1;
	for (i=3; i<xlNumSim; i++)
	{
			err = GetAmortizedSwaptionPayoff5_c(dCopulaCube,
											srtFreq,
											srtBasis,
											dCoverages,
											iNCoInitSwaps,
											dLvls,
											dStrikes,
											dCoInitNotionals,
											i,
											&dPayoff,
											mean,
											flag);

		*dPrice += (DMAX(0,flag * dPayoff) - *dPrice)/(i-2);
//		*dPrice += (flag * dPayoff - *dPrice)/(i-2);
	}

	//-----------------Frees the memory and return-----------------------------
FREE_RETURN:

	if(dPPDx) free_dmatrix(dPPDx,0,iNDates-1,0,2*nClass);
	if(dPPDy) free_dmatrix(dPPDy,0,iNDates-1,0,2*nClass);
	if(dCPDy) free_dmatrix(dCPDy,0,iNDates-1,0,2*nClass);

	if(dPPDShiftedLogx) free_dmatrix(dPPDShiftedLogx,0,iNDates-1,0,2*nClass);

	if(dMean_v) free_dvector(dMean_v,0,iNDates-1);
	if(dAvg) free_dvector(dAvg,0,iNDates-1);
	if(dCopulaCube) free_dmatrix(dCopulaCube,0,xlNumSim,0,iNDates-1);
	if(dCopulaCubeShiftedLog) free_dmatrix(dCopulaCubeShiftedLog,0,xlNumSim,0,iNDates-1);

	if(dCoverages) free(dCoverages);
	if(lStartDates) free(lStartDates);
	if(lEndDates) free(lEndDates);
	if(lPayDates) free(lPayDates);

	if(lCoInitSwapsEndDates) free(lCoInitSwapsEndDates);
	if(dCoInitStrikes) free(dCoInitStrikes);
	if(dCoInitNotionals) free(dCoInitNotionals);

	if(lCpnDates) free(lCpnDates);
	if(lCpnDates) free(dCpns);

	if(dSpreads) free(dSpreads);

	if(mean) free_dvector(mean, 0,iNCoInitSwaps-1);
	if(dLvls) free_dvector(dLvls, 0,iNCoInitSwaps-1);
	if(dStrikes) free_dvector(dStrikes, 0,iNDates-1);

	return err;
}
//--------------------------------------------------------------------------------------
//----------------------Heston Model : End Of Main Function-----------------------------
//--------------------------------------------------------------------------------------
*/

/*
//--------------------------------------------------------------------------------------
//--------------------------Heston Model : Main Function--------------------------------
//--------------------------------------------------------------------------------------
Err		AmortizedSwaptionHeston(
						 char *cYCname,
						 char *cVCname,
						 char *cRefRname,
						 
						 long xlStartDate,
						 long xlEndDate,
						 
						 SrtCompounding srtFreq,
						 SrtBasisCode srtBasis,

						 double exer_fee,
						 
						 long	lNFixNot,
						 double *dFixNotionals,
						 double *dFixRates,

						 long	lNFloatNot,
						 double *dFloatNotionals,
						 double *dMargins,
						 
						 SrtCallPutType srtCallPut,
						 
						 double xlnStdDev,
						 int xlStudDegree,
						 int xlNumSim,
						 int xlNPts,
						 double xlUpBound,
						 int xlnClasses,

						 int UseVol,
						 
						 double **Correl,
						 double *dPrice
			            )
{
	//----------Declarations-----------------------------------------------------------
	//---Err containing result---------------------------------------------------------
	Err err = NULL;
	
	//---Curve-------------------------------------------------------------------------
	SrtCurvePtr	pCurve = lookup_curve (cYCname);
	
	//---SwapDP for schedule generation------------------------------------------------
	SwapDP Swap;

	//---For Dates---------------------------------------------------------------------
	Date lToday,lSpotLag;

	//---Swap Details------------------------------------------------------------------
	char *cBasis, *cFreq;
	Date *lPayDates=NULL,*lStartDates=NULL,*lEndDates=NULL;
	double *dCoverages=NULL;
	int iNDates, iNPayDates;

	//---Monte Carlo Paths-------------------------------------------------------------
	double **dCopulaCube=NULL,**dCopulaCubeShiftedLog=NULL,*dMean_v=NULL,*dAvg=NULL;

	//---Distribution Discretization---------------------------------------------------
	double **dPPDx=NULL,**dPPDy=NULL,**dCPDy=NULL;
	double **dPPDShiftedLogx=NULL;
	
	//---Miscelaneous Integers---------------------------------------------------------
	int i, j;
	int nClass=1<<xlnClasses;
	int flag;

	//---Simulated Payoff--------------------------------------------------------------
	double dPayoff;

	//---Amortized Swap as a Bond------------------------------------------------------
	Date *lCpnDates=NULL;
	double *dCpns=NULL;

	//---Amortized Swap as a Coinit Swaps Portfolio------------------------------------
	Date *lCoInitSwapsEndDates=NULL;
	double *dCoInitStrikes=NULL;
	double *dCoInitNotionals=NULL;

	//---Fwd Swaps Spreads-------------------------------------------------------------
	double *dSpreads=NULL;

	//---Decomposition Strikes---------------------------------------------------------
	double *dStrikes=NULL;

	//---temporary variable for shifts-------------------------------------------------
	double shift;

	//---Heston Parameters-------------------------------------------------------------
	double *SigmaHeston=NULL;
	double *AlphaHeston=NULL;
	double *LambdaHeston=NULL;
	double *BetaHeston=NULL;
	double *RhoHeston=NULL;

	//---Variable for variate control--------------------------------------------------
	double *mean=NULL;
	double *dLvls=NULL;
	double dCov, dFwdSwap;//, dFwdLvl;

	//--------End Of Declarations------------------------------------------------------


	//----------Convert Basis and Freq in SRTcode--------------------------------------
	err = translate_basis(&cBasis, srtBasis);
	err = translate_compounding(&cFreq, srtFreq);

	//-------defines the pointer to the curve and retrieves today and spot dates-------
	pCurve = lookup_curve (cYCname);
	lToday = (Date) get_today_from_curve (pCurve);
	lSpotLag = (Date) get_spotlag_from_curve(pCurve);
	
	//-----Init and Fill SWAPDP--------------------------------------------------------
    err = swp_f_setSwapDP(xlStartDate, xlEndDate, srtFreq, srtBasis, &Swap);
	if(err)
	{
		goto FREE_RETURN;
	}
	err = swp_f_make_FixedLegDatesAndCoverages(&Swap,lToday,&lPayDates,&iNPayDates,
								  			   &lStartDates,&lEndDates,
								  			   &dCoverages,&iNDates);	
	if(err)
	{
		goto FREE_RETURN;
	}

	if(iNDates!=lNFixNot)
	{
		err = "Fix Notionals : Wrong Dimension";
		goto FREE_RETURN;
	}

	err = GetHestonParameters(lToday,
								cYCname,
								cVCname,
								xlStartDate, 
								xlEndDate, 
								srtFreq,
								srtBasis,
								cRefRname,
								xlUpBound,
								xlNPts,
								2,
								xlnStdDev,
								&SigmaHeston,
								&AlphaHeston,
								&BetaHeston,
								&RhoHeston,
								&LambdaHeston);


	//--------Retrieves the densities--------------------------------------------------
	dPPDx=dmatrix(0,iNDates-1,0,2*nClass);
	dPPDy=dmatrix(0,iNDates-1,0,2*nClass);
	dCPDy=dmatrix(0,iNDates-1,0,2*nClass);
	err = Fwd_Swaption_Get_Cum_HestonDensity(
								 lToday,
								 cYCname,
								 cRefRname,
								 srtFreq,
								 srtBasis,

								 iNPayDates,
                                 lPayDates,
								 dCoverages,

								 SigmaHeston,
								 AlphaHeston,
								 LambdaHeston,
								 BetaHeston,
								 RhoHeston,

								 xlNPts,
								 xlUpBound,
								 nClass,

								 dPPDx,
								 dPPDy,
								 dCPDy,

								 9
								);	

	//--------Conversion into shifted Log befor Copula---------------------------------
	dPPDShiftedLogx=dmatrix(0,iNDates-1,0,2*nClass);
	for(i=0;i<iNDates;++i)
	{
		err = swp_f_ForwardRate(lStartDates[0], lEndDates[i], cFreq, cBasis, cYCname, cRefRname, &dFwdSwap);
		for(j=0;j<2*nClass+1;++j)
		{
//			shift = HESTON_CONST*(1.-BetaHeston[i])/BetaHeston[i];
			shift = dFwdSwap*(1.-BetaHeston[i])/BetaHeston[i];
			dPPDShiftedLogx[i][j] = log(shift + dPPDx[i][j]);
		}
	}


	//--------Copula Retrieves the Matrix of Joint distributions-----------------------
	dMean_v=dvector(0,iNDates-1);
	dAvg=dvector(0,iNDates-1);
	for (i=0; i<iNDates; i++) {
		dAvg[i]=0.0;
		dMean_v[i]=0.0;  
	}

	dCopulaCube = dmatrix(0,xlNumSim,0,iNDates-1);
	dCopulaCubeShiftedLog = dmatrix(0,xlNumSim,0,iNDates-1);

//	GetStudentCplDev (xlStudDegree,dMean_v,Correl,xlNumSim,xlStudDegree+iNDates,dPPDShiftedLogx,
//						  dCPDy,2*nClass+1,0,RANDOM_GAUSS,dCopulaCubeShiftedLog);

	CopulaSimul(Correl,xlNumSim,iNDates,dPPDShiftedLogx,dCPDy,2*nClass+1,dCopulaCubeShiftedLog);

	//--------Get the swap rates back from their shifted log---------------------------
	for(i=0;i<iNDates;++i)
	{
		err = swp_f_ForwardRate(lStartDates[0], lEndDates[i], cFreq, cBasis, cYCname, cRefRname, &dFwdSwap);
		for(j=0;j<xlNumSim+1;++j)
		{
//			shift = HESTON_CONST*(1.-BetaHeston[i])/BetaHeston[i];
			shift = dFwdSwap*(1.-BetaHeston[i])/BetaHeston[i];
			dCopulaCube[j][i] = exp(dCopulaCubeShiftedLog[j][i]) - shift;
		}
	}

	//--------Convert AMort Swap in Coinit Swap Portfolio------------------------------
	dStrikes = dvector(0,iNDates-1);
	dCoInitNotionals = dvector(0,iNDates-1);
	err = ConvertAmortSwapWithMarginsInCoInitSwapPortfolio(
								 cYCname,
								 cVCname,
								 cRefRname,
								 srtCallPut,
								 xlStartDate,
								 xlEndDate,
								 srtFreq,
								 srtBasis,
								 exer_fee,
								 lNFixNot,
								 dFixNotionals,
								 dFixRates,
								 lNFloatNot,
								 dFloatNotionals,
								 dMargins,
								 dStrikes,
								 dCoInitNotionals,
								 UseVol);
	if (err)
	{
		goto FREE_RETURN;
	}


	//-------MonteCarlo----------------------------------------------------------------
	if (cFreq == "ANNUAL") {
		dCov=1.0;
	} else if (cFreq == "SEMIANNUAL") {
		dCov=0.5;
	} else if (cFreq == "QUARTERLY") {
		dCov=0.25;
	} else if (cFreq == "MONTHLY") {
		dCov=1./12.;
	}

	mean = dvector(0,lNFixNot-1);
	dLvls = dvector(0,lNFixNot-1);
	for(j=0;j<lNFixNot;++j)
	{
		mean[j] = 0.0;
		err = swp_f_LevelPayment(lStartDates[0], lEndDates[j], cFreq, cBasis, cYCname, cRefRname, &(dLvls[j]));
	}

	*dPrice=0.0;
	flag = -1;
    if (srtCallPut == SRT_CALL) flag=1;
	for (i=3; i<xlNumSim; i++)
	{
		err = GetAmortizedSwaptionPayoff5_c(dCopulaCube,
											srtFreq,
											srtBasis,
											dCoverages,
											lNFixNot,
											dLvls,
											dStrikes,
											dCoInitNotionals,
											i,
											&dPayoff,
											mean,
											flag);

		*dPrice += (DMAX(0,flag * dPayoff) - *dPrice)/(i-2);
	}

	//-----------------Frees the memory and return-----------------------------
FREE_RETURN:

	if(dPPDx) free_dmatrix(dPPDx,0,iNDates-1,0,2*nClass);
	if(dPPDy) free_dmatrix(dPPDy,0,iNDates-1,0,2*nClass);
	if(dCPDy) free_dmatrix(dCPDy,0,iNDates-1,0,2*nClass);

	if(dPPDShiftedLogx) free_dmatrix(dPPDShiftedLogx,0,iNDates-1,0,2*nClass);

	if(dMean_v) free_dvector(dMean_v,0,iNDates-1);
	if(dAvg) free_dvector(dAvg,0,iNDates-1);
	if(dCopulaCube) free_dmatrix(dCopulaCube,0,xlNumSim,0,iNDates-1);
	if(dCopulaCubeShiftedLog) free_dmatrix(dCopulaCubeShiftedLog,0,xlNumSim,0,iNDates-1);

	if(dCoverages) free(dCoverages);
	if(lStartDates) free(lStartDates);
	if(lEndDates) free(lEndDates);
	if(lPayDates) free(lPayDates);

	if(lCoInitSwapsEndDates) free(lCoInitSwapsEndDates);
	if(dCoInitStrikes) free(dCoInitStrikes);

	if(dCoInitNotionals) free_dvector(dCoInitNotionals,0,iNDates-1);

	if(lCpnDates) free(lCpnDates);
	if(lCpnDates) free(dCpns);

	if(mean) free_dvector(mean, 0,lNFixNot-1);
	if(dLvls) free_dvector(dLvls, 0,lNFixNot-1);

	if(dSpreads) free(dSpreads);

	if(dStrikes) free_dvector(dStrikes, 0,iNDates-1);

	if(SigmaHeston) free(SigmaHeston);
	if(AlphaHeston) free(AlphaHeston);
	if(BetaHeston) free(BetaHeston);
	if(RhoHeston) free(RhoHeston);
	if(LambdaHeston) free(LambdaHeston);
								
	return err;
}
//--------------------------------------------------------------------------------------
//----------------------Heston Model : End Of Main Function-----------------------------
//--------------------------------------------------------------------------------------
*/

//--------------------------------------------------------------------------------------
//-----------------------SABR Model : Utilities Functions-------------------------------
//--------------------------------------------------------------------------------------

/*
//------------------------------------------------------------
//----------regular mesh and finite difference within mesh ---
//------------------------------------------------------------
Err  GetSABRDensity_Lvl_Measure(
						long today,
						int		nClass,
						char	*cYCname,
						char	*cVCname,
						char	*cRefRname,
						SrtCompounding srtFreq,
						SrtBasisCode srtBasis,
						long start,
						long end,
						long nTen,
						double dCov,
						double dDf,
						double UpBound,
						double **dPPDx,
						double **dPPDy,
						double **dCPDy
						)
{
	double LvlC_x;
	double *dVec = NULL, *dX=NULL,*dY=NULL, dThreshold=1.e-3;
	int i;
	double deriv1, deriv2;
	Err err = NULL;
	double epsilon = 1.0e-4;
	double temp, temp_plus, temp_moins;
	double dK;
	double dFwdSwap, dFwdCash, dSpread, dLvl, ex_time;
	double sigma, alpha, beta, rho;
	char *cBasis, *cFreq;
	double power;
	
	err = translate_basis(&cBasis, srtBasis);
	err = translate_compounding(&cFreq, srtFreq);

	err = swp_f_ForwardRate(start, end, cFreq, cBasis, cYCname, cRefRname, &dFwdSwap);
	dFwdCash = swp_f_swapcashrate( start, end, srtBasis, srtFreq, cYCname, cRefRname);
	err = swp_f_LevelPayment(start, end, cFreq, cBasis, cYCname, cRefRname, &dLvl);
	ex_time = (add_unit (start, - 2, SRT_BDAY, MODIFIED_SUCCEEDING) - today) * YEARS_IN_DAY;

	dSpread = dFwdSwap - dFwdCash;
	
	///// gets the two boundaries
	dPPDx[nTen][0]=0.0001;
	dPPDx[nTen][2*nClass]=UpBound;

	dK = (dPPDx[nTen][2*nClass]-dPPDx[nTen][0])/(2*nClass);

	for (i=1;i<2*nClass;i++)
	{
		dPPDx[nTen][i]=dPPDx[nTen][0] + i * dK;
	}

	///// runs through the points 
	for (i=0;i<=2*nClass;i++)
	{
		dCPDy[nTen][i]=0.0;
		dPPDy[nTen][i]=0.0;
	}

	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&sigma, &power, SABR_BETAVOL );
	if (err) return err;
	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&alpha, &power, SABR_ALPHA );
	if (err) return err;
	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&beta, &power, SABR_BETA );
	if (err) return err;
	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&rho, &power, SABR_RHO );
	if (err) return err;


	////i=0

	temp = srt_f_optblkschbetastochquick(
								dFwdSwap,
								dPPDx[nTen][0],
								ex_time,
								sigma,
								alpha,
								beta,
								rho,
								1.0,
								SRT_BETAVOL, 
								SRT_LOGNORMAL, 
								SRT_CALL, 
								PREMIUM);

	temp_plus = srt_f_optblkschbetastochquick(
								dFwdSwap,
								dPPDx[nTen][1],
								ex_time,
								sigma,
								alpha,
								beta,
								rho,
								1.0,
								SRT_BETAVOL, 
								SRT_LOGNORMAL, 
								SRT_CALL, 
								PREMIUM);

	LvlC_x = dDf * (1.-1./pow(1+dPPDx[nTen][0]*dCov,nTen+1))/dPPDx[nTen][0];

	deriv1 = 1.0;
	deriv2 = DMIN(deriv1,(temp - temp_plus)/dK);

	dPPDy[nTen][0]=DMAX(0,(deriv1 - deriv2)/dK)*dLvl/LvlC_x;
	dCPDy[nTen][0]=0.0;
	
	for (i=1;i<2*nClass;i++)
	{
		temp_moins = temp;
		temp = temp_plus;

		temp_plus = srt_f_optblkschbetastochquick(
								dFwdSwap,
								dPPDx[nTen][i+1],
								ex_time,
								sigma,
								alpha,
								beta,
								rho,
								1.0,
								SRT_BETAVOL, 
								SRT_LOGNORMAL, 
								SRT_CALL, 
								PREMIUM);

		LvlC_x = dDf * (1.-1./pow(1+dPPDx[nTen][i]*dCov,nTen+1))/dPPDx[nTen][i];

		deriv1 = DMIN(1.0,(temp_moins - temp)/dK);
		deriv2 = DMIN(deriv1,(temp - temp_plus)/dK);

		dPPDy[nTen][i]=DMAX(0,(deriv1 - deriv2)/dK)*dLvl/LvlC_x;

		dCPDy[nTen][i]=dCPDy[nTen][i-1]
			+(dPPDy[nTen][i]+dPPDy[nTen][i-1])/2.*dK;
	}

	temp_moins = temp;
	temp = temp_plus;

	LvlC_x = dDf * (1.-1./pow(1+dPPDx[nTen][2*nClass]*dCov,nTen+1))/dPPDx[nTen][2*nClass];

	deriv1 = (temp_moins - temp)/dK;
	deriv2 = 0.0;

	dPPDy[nTen][2*nClass]=DMAX(0,(deriv1 - deriv2)/dK)*dLvl/LvlC_x;

	dCPDy[nTen][2*nClass]=dCPDy[nTen][2*nClass-1]
		+(dPPDy[nTen][2*nClass]+dPPDy[nTen][2*nClass-1])/2.*dK;

	for (i=0;i<=2*nClass;i++) 
	{
		dCPDy[nTen][i]/=dCPDy[nTen][2*nClass];
	}

	return err;

}
*/

/*

Err  GetSABRDensity_Lvl_Measure2(
						long today,
						int		nClass,
						char	*cYCname,
						char	*cVCname,
						char	*cRefRname,
						SrtCompounding srtFreq,
						SrtBasisCode srtBasis,
						long start,
						long end,
						long nTen,
						double dCov,
						double dDf,
						double UpBound,
						double **dPPDx,
						double **dPPDy,
						double **dCPDy
						)
{
	double LvlC_x;
	double *dVec = NULL, *dX=NULL,*dY=NULL, dThreshold=1.e-3;
	int i;
	double deriv1, deriv2;
	Err err = NULL;
	double epsilon = 1.0e-4;
	double temp, temp_plus, temp_moins;
	double dK;
	double dFwdSwap, dFwdCash, dSpread, dLvl, ex_time;
	double sigma, alpha, beta, rho;
	char *cBasis, *cFreq;
	double power;
	long kfail;
	double PrMin, PrAdj;
	
	err = translate_basis(&cBasis, srtBasis);
	err = translate_compounding(&cFreq, srtFreq);

	err = swp_f_ForwardRate(start, end, cFreq, cBasis, cYCname, cRefRname, &dFwdSwap);
	dFwdCash = swp_f_swapcashrate( start, end, srtBasis, srtFreq, cYCname, cRefRname);
	err = swp_f_LevelPayment(start, end, cFreq, cBasis, cYCname, cRefRname, &dLvl);
	ex_time = (add_unit (start, - 2, SRT_BDAY, MODIFIED_SUCCEEDING) - today) * YEARS_IN_DAY;

	dSpread = dFwdSwap - dFwdCash;
	
	///// gets the two boundaries
	dPPDx[nTen][0]=0.0001;
	dPPDx[nTen][2*nClass]=UpBound;

	dK = (dPPDx[nTen][2*nClass]-dPPDx[nTen][0])/(2*nClass);

	for (i=1;i<2*nClass;i++)
	{
		dPPDx[nTen][i]=dPPDx[nTen][0] + i * dK;
	}

	///// runs through the points 
	for (i=0;i<=2*nClass;i++)
	{
		dCPDy[nTen][i]=0.0;
		dPPDy[nTen][i]=0.0;
	}

	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&sigma, &power, SABR_BETAVOL );
	if (err) return err;
	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&alpha, &power, SABR_ALPHA );
	if (err) return err;
	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&beta, &power, SABR_BETA );
	if (err) return err;
	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&rho, &power, SABR_RHO );
	if (err) return err;


	////i=2N

	temp_moins = srt_f_optblkschbetastochquick(
								dFwdSwap,
								dPPDx[nTen][2*nClass],
								ex_time,
								sigma,
								alpha,
								beta,
								rho,
								1.0,
								SRT_BETAVOL, 
								SRT_LOGNORMAL, 
								SRT_CALL, 
								PREMIUM);

	temp = temp_moins;
	temp_plus = temp_moins;

	LvlC_x = dDf * (1.-1./pow(1+dPPDx[nTen][0]*dCov,nTen+1))/dPPDx[nTen][0];

	deriv1 = (temp_moins - temp)/dK;
	deriv2 = (temp - temp_plus)/dK;

	dPPDy[nTen][2*nClass]=(deriv1 - deriv2)/dK*dLvl/LvlC_x;
	dCPDy[nTen][2*nClass]=1.0;
	
	kfail = 0;
	for (i=2*nClass;i>0;--i)
	{
		temp_plus = temp;
		temp = temp_moins;

		temp_moins = srt_f_optblkschbetastochquick(
								dFwdSwap,
								dPPDx[nTen][i-1],
								ex_time,
								sigma,
								alpha,
								beta,
								rho,
								1.0,
								SRT_BETAVOL, 
								SRT_LOGNORMAL, 
								SRT_CALL, 
								PREMIUM);

		LvlC_x = dDf * (1.-1./pow(1+dPPDx[nTen][i]*dCov,nTen+1))/dPPDx[nTen][i];

		deriv1 = (temp_moins - temp)/dK;
		deriv2 = (temp - temp_plus)/dK;

		dPPDy[nTen][i-1]=(deriv1 - deriv2)/dK*dLvl/LvlC_x;

		dCPDy[nTen][i-1]=DMAX(0,dCPDy[nTen][i]-(dPPDy[nTen][i]+dPPDy[nTen][i-1])/2.*dK);

		if(dCPDy[nTen][i-1]>dCPDy[nTen][i])
		{
			kfail = i;
			PrMin = dCPDy[nTen][i];
			i=0;
		}
	}


//-----------------------------------------------------------------------------
//---------------------------SABR Failed ?-------------------------------------
//-----------------------------------------------------------------------------
	while(kfail>0)
	{
		alpha = alpha * 0.9;
		temp_moins = srt_f_optblkschbetastochquick(
								dFwdSwap,
								dPPDx[nTen][kfail],
								ex_time,
								sigma,
								alpha,
								beta,
								rho,
								1.0,
								SRT_BETAVOL, 
								SRT_LOGNORMAL, 
								SRT_CALL, 
								PREMIUM);

		temp = srt_f_optblkschbetastochquick(
								dFwdSwap,
								dPPDx[nTen][kfail+1],
								ex_time,
								sigma,
								alpha,
								beta,
								rho,
								1.0,
								SRT_BETAVOL, 
								SRT_LOGNORMAL, 
								SRT_CALL, 
								PREMIUM);

		temp_plus = srt_f_optblkschbetastochquick(
								dFwdSwap,
								dPPDx[nTen][kfail+2],
								ex_time,
								sigma,
								alpha,
								beta,
								rho,
								1.0,
								SRT_BETAVOL, 
								SRT_LOGNORMAL, 
								SRT_CALL, 
								PREMIUM);

		LvlC_x = dDf * (1.-1./pow(1+dPPDx[nTen][kfail]*dCov,nTen+1))/dPPDx[nTen][kfail];

		deriv1 = (temp_moins - temp)/dK;
		deriv2 = (temp - temp_plus)/dK;

		PrAdj = PrMin/deriv2;
	
		for (i=kfail;i>0;--i)
		{
			kfail = 0;
			temp_plus = temp;
			temp = temp_moins;

			temp_moins = srt_f_optblkschbetastochquick(
								dFwdSwap,
								dPPDx[nTen][i-1],
								ex_time,
								sigma,
								alpha,
								beta,
								rho,
								1.0,
								SRT_BETAVOL, 
								SRT_LOGNORMAL, 
								SRT_CALL, 
								PREMIUM);

			LvlC_x = dDf * (1.-1./pow(1+dPPDx[nTen][i]*dCov,nTen+1))/dPPDx[nTen][i];

			deriv1 = (temp_moins - temp)/dK;
			deriv2 = (temp - temp_plus)/dK;

			dPPDy[nTen][i-1]=(deriv1 - deriv2)/dK*dLvl/LvlC_x;

			dCPDy[nTen][i-1]=DMAX(0,dCPDy[nTen][i]-(dPPDy[nTen][i]+dPPDy[nTen][i-1])/2.*dK/PrAdj);
			if(dCPDy[nTen][i-1]>dCPDy[nTen][i])
			{
				kfail = i;
				PrMin = dCPDy[nTen][i];
				i=0;
			}

		}

	}

	return err;

}

*/

/*

Err  GetSABRDensity_Lvl_Measure3(
						long today,
						int		nClass,
						char	*cYCname,
						char	*cVCname,
						char	*cRefRname,
						SrtCompounding srtFreq,
						SrtBasisCode srtBasis,
						long start,
						long end,
						long nTen,
						double dCov,
						double dDf,
						double UpBound,
						double **dPPDx,
						double **dPPDy,
						double **dCPDy
						)
{
	double LvlC_x;
	double *dVec = NULL, *dX=NULL,*dY=NULL, dThreshold=1.e-3;
	int i;
	double deriv1, deriv2;
	Err err = NULL;
	double epsilon = 1.0e-4;
	double temp, temp_plus, temp_moins;
	double dK;
	double dFwdSwap, dFwdCash, dSpread, dLvl, ex_time;
	double sigma, alpha, beta, rho;
	char *cBasis, *cFreq;
	double power;
	long kfail;
	double PrMin;
	double x_PrMin, discriminent, shifted_vol, shift;
	
	err = translate_basis(&cBasis, srtBasis);
	err = translate_compounding(&cFreq, srtFreq);

	err = swp_f_ForwardRate(start, end, cFreq, cBasis, cYCname, cRefRname, &dFwdSwap);
	dFwdCash = swp_f_swapcashrate( start, end, srtBasis, srtFreq, cYCname, cRefRname);
	err = swp_f_LevelPayment(start, end, cFreq, cBasis, cYCname, cRefRname, &dLvl);
	ex_time = (add_unit (start, - 2, SRT_BDAY, MODIFIED_SUCCEEDING) - today) * YEARS_IN_DAY;

	dSpread = dFwdSwap - dFwdCash;
	
	///// gets the two boundaries
	dPPDx[nTen][0]=0.0001;
	dPPDx[nTen][2*nClass]=UpBound;

	dK = (dPPDx[nTen][2*nClass]-dPPDx[nTen][0])/(2*nClass);

	for (i=1;i<2*nClass;i++)
	{
		dPPDx[nTen][i]=dPPDx[nTen][0] + i * dK;
	}

	///// runs through the points 
	for (i=0;i<=2*nClass;i++)
	{
		dCPDy[nTen][i]=0.0;
		dPPDy[nTen][i]=0.0;
	}

	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&sigma, &power, SABR_BETAVOL );
	if (err) return err;
	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&alpha, &power, SABR_ALPHA );
	if (err) return err;
	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&beta, &power, SABR_BETA );
	if (err) return err;
	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&rho, &power, SABR_RHO );
	if (err) return err;


	////i=2N

	temp_moins = srt_f_optblkschbetastochquick(
								dFwdSwap,
								dPPDx[nTen][2*nClass],
								ex_time,
								sigma,
								alpha,
								beta,
								rho,
								1.0,
								SRT_BETAVOL, 
								SRT_LOGNORMAL, 
								SRT_CALL, 
								PREMIUM);

	temp = temp_moins;
	temp_plus = temp_moins;

	LvlC_x = dDf * (1.-1./pow(1+dPPDx[nTen][0]*dCov,nTen+1))/dPPDx[nTen][0];

	deriv1 = (temp_moins - temp)/dK;
	deriv2 = (temp - temp_plus)/dK;

	dPPDy[nTen][2*nClass]=(deriv1 - deriv2)/dK*dLvl/LvlC_x;
	dCPDy[nTen][2*nClass]=1.0;
	
	kfail = 0;
	for (i=2*nClass;i>0;--i)
	{
		temp_plus = temp;
		temp = temp_moins;

		temp_moins = srt_f_optblkschbetastochquick(
								dFwdSwap,
								dPPDx[nTen][i-1],
								ex_time,
								sigma,
								alpha,
								beta,
								rho,
								1.0,
								SRT_BETAVOL, 
								SRT_LOGNORMAL, 
								SRT_CALL, 
								PREMIUM);

		LvlC_x = dDf * (1.-1./pow(1+dPPDx[nTen][i]*dCov,nTen+1))/dPPDx[nTen][i];

		deriv1 = (temp_moins - temp)/dK;
		deriv2 = (temp - temp_plus)/dK;

		dPPDy[nTen][i-1]=(deriv1 - deriv2)/dK*dLvl/LvlC_x;

		dCPDy[nTen][i-1]=DMAX(0,dCPDy[nTen][i]-(dPPDy[nTen][i]+dPPDy[nTen][i-1])/2.*dK);

		if(dCPDy[nTen][i-1]>dCPDy[nTen][i])
		{
			kfail = i;
			PrMin = dCPDy[nTen][i];
			x_PrMin = inv_cumnorm_fast(PrMin);
//			shift = HESTON_CONST*(1.-beta)/beta;
			shift = 0.0;
			discriminent = x_PrMin * x_PrMin * ex_time - 2 * ex_time * log((shift+dPPDx[nTen][i])/(shift+dFwdSwap));
			shifted_vol = (x_PrMin*sqrt(ex_time) + sqrt(discriminent))/ex_time;
			i=0;
		}
	}


//-----------------------------------------------------------------------------
//---------------------------SABR Failed ?-------------------------------------
//-----------------------------------------------------------------------------
	if(kfail>0)
	{
		temp_moins = srt_f_optblksch(
								shift+dFwdSwap,
								shift+dPPDx[nTen][kfail],
								shifted_vol,
								ex_time,
								1.0,
								SRT_CALL, 
								PREMIUM);

		temp = srt_f_optblksch(
								shift+dFwdSwap,
								shift+dPPDx[nTen][kfail+1],
								shifted_vol,
								ex_time,
								1.0,
								SRT_CALL, 
								PREMIUM);

		for (i=kfail;i>0;--i)
		{
			kfail = 0;
			temp_plus = temp;
			temp = temp_moins;

			temp_moins = srt_f_optblksch(
								shift+dFwdSwap,
								shift+dPPDx[nTen][i-1],
								shifted_vol,
								ex_time,
								1.0,
								SRT_CALL, 
								PREMIUM);

			LvlC_x = dDf * (1.-1./pow(1+dPPDx[nTen][i]*dCov,nTen+1))/dPPDx[nTen][i];

			deriv1 = (temp_moins - temp)/dK;
			deriv2 = (temp - temp_plus)/dK;

			dPPDy[nTen][i-1]=(deriv1 - deriv2)/dK*dLvl/LvlC_x;

			dCPDy[nTen][i-1]=DMAX(0,dCPDy[nTen][i]-(dPPDy[nTen][i]+dPPDy[nTen][i-1])/2.*dK);
		}
	}

	return err;

}
*/

/*
Err  impliedATMSigmaBeta(double atm_premium,
					  double Fwd,
					  double ex_time,
					  double firstguess_betavol,
					  double alpha,
					  double beta,
					  double rho,
					  double disc,
					  double *sigmabeta)
{
	double betavol, betavol2;
	double atm_model, atm_model2;
	Err err=NULL;

	betavol = firstguess_betavol;
	atm_model = srt_f_optblkschbetastochquick(
								Fwd,
								Fwd,
								ex_time,
								betavol,
								alpha,
								beta,
								rho,
								disc,
								SRT_BETAVOL, 
								SRT_LOGNORMAL, 
								SRT_CALL, 
								PREMIUM);

	while(fabs(atm_model - atm_premium)>1e-7)
	{
		betavol2 = betavol + 1e-5;
		atm_model2 = srt_f_optblkschbetastochquick(
								Fwd,
								Fwd,
								ex_time,
								betavol2,
								alpha,
								beta,
								rho,
								disc,
								SRT_BETAVOL, 
								SRT_LOGNORMAL, 
								SRT_CALL, 
								PREMIUM);

		betavol=betavol+(atm_premium-atm_model)*(betavol2-betavol)/(atm_model2-atm_model);
		atm_model = srt_f_optblkschbetastochquick(
								Fwd,
								Fwd,
								ex_time,
								betavol,
								alpha,
								beta,
								rho,
								disc,
								SRT_BETAVOL, 
								SRT_LOGNORMAL, 
								SRT_CALL, 
								PREMIUM);

	}

	*sigmabeta = betavol;

	return err;

}
*/

/*

Err  GetSABRDensity_Lvl_Measure4(
						long today,
						int		nClass,
						char	*cYCname,
						char	*cVCname,
						char	*cRefRname,
						SrtCompounding srtFreq,
						SrtBasisCode srtBasis,
						long start,
						long end,
						long nTen,
						double dCov,
						double dDf,
						double UpBound,
						double **dPPDx,
						double **dPPDy,
						double **dCPDy,
						double *dShifts,
						double *dMeans
						)
{
	double LvlC_x;
	double *dVec = NULL, *dX=NULL,*dY=NULL, dThreshold=1.e-3;
	int i;
	double deriv1, deriv2;
	Err err = NULL;
	double epsilon = 1.0e-4;
	double temp, temp_plus, temp_moins;
	double dK;
	double dFwdSwap, dFwdCash, dSpread, dLvl, ex_time;
	double sigma, alpha, beta, rho;
	char *cBasis, *cFreq;
	double power;
	double shift, shifted_vol, atm_premium;
	
	err = translate_basis(&cBasis, srtBasis);
	err = translate_compounding(&cFreq, srtFreq);

	err = swp_f_ForwardRate(start, end, cFreq, cBasis, cYCname, cRefRname, &dFwdSwap);
	dFwdCash = swp_f_swapcashrate( start, end, srtBasis, srtFreq, cYCname, cRefRname);
	err = swp_f_LevelPayment(start, end, cFreq, cBasis, cYCname, cRefRname, &dLvl);
	ex_time = (add_unit (start, - 2, SRT_BDAY, MODIFIED_SUCCEEDING) - today) * YEARS_IN_DAY;

	dSpread = dFwdSwap - dFwdCash;


	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&sigma, &power, SABR_BETAVOL );
	if (err) return err;
	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&alpha, &power, SABR_ALPHA );
	if (err) return err;
	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&beta, &power, SABR_BETA );
	if (err) return err;
	err = swp_f_SABRvol( cVCname, start, end, 0.05,
			&rho, &power, SABR_RHO );
	if (err) return err;

//	shift = HESTON_CONST*(1.-beta)/beta;
	shift = dFwdSwap*(1.-beta)/beta;

	*dShifts = shift;

	atm_premium = srt_f_optblkschbetastochquick(
								dFwdSwap,
								dFwdSwap,
								ex_time,
								sigma,
								alpha,
								beta,
								rho,
								1.0,
								SRT_BETAVOL, 
								SRT_LOGNORMAL, 
								SRT_CALL, 
								PREMIUM);

	err = srt_f_optimpvol(atm_premium,
							shift+dFwdSwap, 
							shift+dFwdSwap,
							ex_time, 
							1.0, 
							SRT_CALL,
							SRT_LOGNORMAL,
							&shifted_vol);

	err = impliedATMSigmaBeta(atm_premium,
								shift+dFwdSwap, 
								ex_time,
								shifted_vol,
								alpha,
								1.0,
								rho,
								1.0,
								&shifted_vol);
	
	///// gets the two boundaries
	dPPDx[nTen][0]=-shift + 0.01;
	dPPDx[nTen][2*nClass]=UpBound;

	dK = (dPPDx[nTen][2*nClass]-dPPDx[nTen][0])/(2*nClass);

	for (i=1;i<2*nClass;i++)
	{
		dPPDx[nTen][i]=dPPDx[nTen][0] + i * dK;
	}

	///// runs through the points 
	for (i=0;i<=2*nClass;i++)
	{
		dCPDy[nTen][i]=0.0;
		dPPDy[nTen][i]=0.0;
	}


	////i=0
	temp_moins = srt_f_optblkschbetastochquick(
								shift+dFwdSwap,
								DMAX(0.0001, shift+dPPDx[nTen][0] - dK),
								ex_time,
								shifted_vol,
								alpha,
								1.0,
								rho,
								1.0,
								SRT_BETAVOL, 
								SRT_LOGNORMAL, 
								SRT_CALL, 
								PREMIUM);

	temp = srt_f_optblkschbetastochquick(
								shift+dFwdSwap,
								shift+dPPDx[nTen][0],
								ex_time,
								shifted_vol,
								alpha,
								1.0,
								rho,
								1.0,
								SRT_BETAVOL, 
								SRT_LOGNORMAL, 
								SRT_CALL, 
								PREMIUM);

	temp_plus = srt_f_optblkschbetastochquick(
								shift+dFwdSwap,
								shift+dPPDx[nTen][1],
								ex_time,
								shifted_vol,
								alpha,
								1.0,
								rho,
								1.0,
								SRT_BETAVOL, 
								SRT_LOGNORMAL, 
								SRT_CALL, 
								PREMIUM);

	LvlC_x = dDf * (1.-1./pow(1+dPPDx[nTen][0]*dCov,nTen+1))/dPPDx[nTen][0];

//	deriv1 = 1.0;
	deriv1 = DMIN(1.0,(temp_moins - temp)/dK);
	deriv2 = DMIN(deriv1,(temp - temp_plus)/dK);

	dPPDy[nTen][0]=DMAX(0,(deriv1 - deriv2)/dK)*dLvl/LvlC_x;
	dCPDy[nTen][0]=0.0;
	
	for (i=1;i<2*nClass;i++)
	{
		temp_moins = temp;
		temp = temp_plus;

		temp_plus = srt_f_optblkschbetastochquick(
								shift+dFwdSwap,
								shift+dPPDx[nTen][i+1],
								ex_time,
								shifted_vol,
								alpha,
								1.0,
								rho,
								1.0,
								SRT_BETAVOL, 
								SRT_LOGNORMAL, 
								SRT_CALL, 
								PREMIUM);

		LvlC_x = dDf * (1.-1./pow(1+dPPDx[nTen][i]*dCov,nTen+1))/dPPDx[nTen][i];

		deriv1 = DMIN(1.0,(temp_moins - temp)/dK);
		deriv2 = DMIN(deriv1,(temp - temp_plus)/dK);

		dPPDy[nTen][i]=DMAX(0,(deriv1 - deriv2)/dK)*dLvl/LvlC_x;

		dCPDy[nTen][i]=dCPDy[nTen][i-1]
			+(dPPDy[nTen][i]+dPPDy[nTen][i-1])/2.*dK;
	}

	temp_moins = temp;
	temp = temp_plus;

	LvlC_x = dDf * (1.-1./pow(1+dPPDx[nTen][2*nClass]*dCov,nTen+1))/dPPDx[nTen][2*nClass];

	deriv1 = (temp_moins - temp)/dK;
//	deriv2 = 0.0;
	deriv2 = (temp - temp_plus)/dK;

	dPPDy[nTen][2*nClass]=DMAX(0,(deriv1 - deriv2)/dK)*dLvl/LvlC_x;

	dCPDy[nTen][2*nClass]=dCPDy[nTen][2*nClass-1]
		+(dPPDy[nTen][2*nClass]+dPPDy[nTen][2*nClass-1])/2.*dK;

	*dMeans = 0;
	for (i=0;i<=2*nClass;i++) 
	{
		dCPDy[nTen][i]/=dCPDy[nTen][2*nClass];
		*dMeans += (1-dCPDy[nTen][i])*dK;
	}
	*dMeans = *dMeans - shift;

	for (i=0;i<=2*nClass;i++) 
	{
		dCPDy[nTen][i]/=dCPDy[nTen][2*nClass];
		*dMeans += (1-dCPDy[nTen][i])*dK;
	}

	return err;

}
*/

/*
//--------------------------------------------------------------------------------------
//-----------Function that transforms smile from the swap measure into Fwd measure------
//--------------------------------------------------------------------------------------
Err       Fwd_Swaption_Get_Cum_SABRDensity(Date lToday,
								 char *cYCname,
								 char *cVCname,
						         char *cRefRname,
								 SrtCompounding srtFreq,
								 SrtBasisCode srtBasis,

								 int iNPayDatesLong,
								 Date *lPayDatesLong,
								 double *dCoveragesLong,

								 double dUpBound,
								 int nClass,

								 double **dPPDx,
                                 double **dPPDy,
							     double **dCPDy,

								 int iDistribMethod,

								 double *dStrikes,

								 double *dShifts,

								 double *dMeans
								 )

{
	Err err = NULL;
	int i;
	double dMat,dCov;
	char *cBasis, *cFreq;
	double shift, mean;
	double dDf;

	dMat=(lPayDatesLong[0]-lToday)/365.; 
//		dDisc=swp_f_df(lToday, lPayDatesLong[0], cYCname);
	
//////////////////////////////////////////
/// translate SORT types into strings
//////////////////////////////////////////

	err = translate_basis(&cBasis, srtBasis);
	err = translate_compounding(&cFreq, srtFreq);

	if (cFreq == "ANNUAL") {
		dCov=1.0;
	} else if (cFreq == "SEMIANNUAL") {
		dCov=0.5;
	} else if (cFreq == "QUARTERLY") {
		dCov=0.25;
	} else if (cFreq == "MONTHLY") {
		dCov=1./12.;
	}

//////////////////////////////////////////////////////////
////////////////// generates the cumulative densities
//////////////////////////////////////////////////////////	

	dDf=swp_f_df(lToday, lPayDatesLong[0], cYCname);
	for (i=1;i<iNPayDatesLong;i++)
	{
		err = GetSABRDensity_Lvl_Measure4(lToday, nClass,
										cYCname,cVCname,cRefRname,
										srtFreq,srtBasis,
										lPayDatesLong[0],
										lPayDatesLong[i],
										i-1,
										dCov,
										dDf,
										dUpBound,
										dPPDx,
										dPPDy,
										dCPDy,
										&shift,
										&mean);

		dShifts[i-1] = shift;
		dMeans[i-1] = mean;
	}

	return err;

}	
//--------------------------------------------------------------------------------------
//-------End Of Function that transforms smile from the swap measure into Fwd measure---
//--------------------------------------------------------------------------------------
*/

/*
//--------------------------------------------------------------------------------------
//--------------------------SABR Model : Main Function----------------------------------
//--------------------------------------------------------------------------------------
Err		PriceAmortizedSwaptionSABR(
						 char *cYCname,
						 char *cVCname,
						 char *cRefRname,

						 long xlStartDate,
						 long xlEndDate,

						 SrtCompounding srtFreq,
						 SrtBasisCode srtBasis,
						 
						 long	lNFixNot,
						 double *dFixNotionals,
						 long	lNFloatNot,
						 double *dFloatNotionals,
						 double *dFixRates,

						 SrtCallPutType srtCallPut,

						 int xlStudDegree,
						 int xlNumSim,
						 double xlUpBound,
						 int xlnClasses,
						 double **Correl,
						 double *dPrice,
						 int iPayoffMethod,
						 int iDistribMethod
			            )
{
	//----------Declarations-----------------------------------------------------------
	//---Err containing result---------------------------------------------------------
	Err err = NULL;

	SrtCurvePtr	pCurve = lookup_curve (cYCname);

	Date lToday,lSpotLag;

	SwapDP Swap;
	Date *lPayDates=NULL,*lStartDates=NULL,*lEndDates=NULL;
	double *dCoverages=NULL;
	int iNDates, iNPayDates;

	double **dPPDx=NULL,**dPPDy=NULL,**dCPDy=NULL;
	double **dPPDShiftedLogx=NULL;
	double **dCopulaCube=NULL,**dCopulaCubeShiftedLog=NULL,*dMean_v=NULL,*dAvg=NULL;

	int i, j, nClass=1<<xlnClasses, flag;

	double dPayoff;

//	int iNCpns;
	Date *lCpnDates=NULL;
	double *dCpns=NULL;

	int iNCoInitSwaps;
	Date *lCoInitSwapsEndDates=NULL;
	double *dCoInitStrikes=NULL;
	double *dCoInitNotionals=NULL;

//	int NSpreads;
	double *dSpreads=NULL;
	
	double *dShifts=NULL;
	double *dStrikes=NULL;
	double *dMeans=NULL;
	
	double df_start;//, dSwapPV;

	char *cBasis, *cFreq;

	err = translate_basis(&cBasis, srtBasis);
	if(err)
	{
		goto FREE_RETURN;
	}
	err = translate_compounding(&cFreq, srtFreq);
	if(err)
	{
		goto FREE_RETURN;
	}

	//------defines the pointer to the curve and retrieves today and spot dates------------
	pCurve = lookup_curve (cYCname);
	///// if (!pCurve) throw XSort("Fatal: (FwdSwaption) Yield Curve Object not found");
	//------Get Today and SpotLag----------------------------------------------------------
	lToday = (Date) get_today_from_curve (pCurve);
	lSpotLag = (Date) get_spotlag_from_curve(pCurve);	
	
	
	//------Inits and fills swapDP---------------------------------------------------------
    err = swp_f_setSwapDP(xlStartDate, xlEndDate, srtFreq, srtBasis, &Swap);
	if(err)
	{
		goto FREE_RETURN;
	}
														
	err = swp_f_make_FixedLegDatesAndCoverages(&Swap,lToday,&lPayDates,&iNPayDates,
								  			   &lStartDates,&lEndDates,
								  			   &dCoverages,&iNDates);
	if(err)
	{
		goto FREE_RETURN;
	}
	
//-----------------------------------------------------------------------
//------------------------Compute Calibration Strikes--------------------
//-----------------------------------------------------------------------
	dStrikes = dvector(0,iNDates-1);
	err = ComputeIRRStrikes(cYCname, lToday, cRefRname,
						 xlStartDate, xlEndDate,
						 srtFreq,srtBasis,
						 lNFixNot,dFixNotionals,
						 lNFloatNot,dFloatNotionals,dFixRates,
						 dStrikes);
	if (err)
	{
		goto FREE_RETURN;
	}

	//-----------Retrieves the densities---------------------------------------------------
	dPPDx=dmatrix(0,iNDates-1,0,2*nClass);
	dPPDy=dmatrix(0,iNDates-1,0,2*nClass);
	dCPDy=dmatrix(0,iNDates-1,0,2*nClass);

	dShifts = dvector(0,iNDates-1);
	dMeans = dvector(0,iNDates-1);
	err = Fwd_Swaption_Get_Cum_SABRDensity(
								 lToday,
								 cYCname,
								 cVCname,
								 cRefRname,
								 srtFreq,
								 srtBasis,

								 iNPayDates,
                                 lPayDates,
								 dCoverages,

								 xlUpBound,
								 nClass,

								 dPPDx,
								 dPPDy,
								 dCPDy,

								 iDistribMethod,

								 dStrikes,

								 dShifts,

								 dMeans
								);
	if(err)
	{
		goto FREE_RETURN;
	}





	//--------Conversion into shifted Log befor Copula---------------------------------
	dPPDShiftedLogx=dmatrix(0,iNDates-1,0,2*nClass);
	for(i=0;i<iNDates;++i)
	{
		for(j=0;j<2*nClass+1;++j)
		{
			dPPDShiftedLogx[i][j] = log(dShifts[i] + dPPDx[i][j]);
		}
	}

	
	//--------Copula Retrieves the Matrix of Joint distributions-----------------------
	dMean_v=dvector(0,iNDates-1);
	dAvg=dvector(0,iNDates-1);
	for (i=0; i<iNDates; i++) {
		dAvg[i]=0.0;
		dMean_v[i]=0.0;  
	}

	dCopulaCube = dmatrix(0,xlNumSim,0,iNDates-1);
	dCopulaCubeShiftedLog = dmatrix(0,xlNumSim,0,iNDates-1);

//	GetStudentCplDev (xlStudDegree,dMean_v,Correl,xlNumSim,xlStudDegree+iNDates,dPPDShiftedLogx,
//						  dCPDy,2*nClass+1,0,RANDOM_GAUSS,dCopulaCubeShiftedLog);

	CopulaSimul(Correl,xlNumSim,iNDates,dPPDShiftedLogx,dCPDy,2*nClass+1,dCopulaCubeShiftedLog);

	//--------Get the swap rates back from their shifted log---------------------------
	for(i=0;i<iNDates;++i)
	{
		for(j=0;j<xlNumSim+1;++j)
		{
			dCopulaCube[j][i] = exp(dCopulaCubeShiftedLog[j][i]) - dShifts[i];
		}
	}


	//--------Convert AMort Swap in Coinit Swap Portfolio------------------------------
	err=ConvertAmortSwapInCoInitalSwapsPortfolio2(
							cYCname,
							cRefRname,
							xlStartDate,
							xlEndDate,
							srtFreq,
							srtBasis,
							lNFixNot,
							dFixNotionals,
							dFixRates,
							dStrikes,
							&iNCoInitSwaps,
//							&lCoInitSwapsEndDates,
//							&dCoInitStrikes,
							&dCoInitNotionals);
	if (err)
	{
		goto FREE_RETURN;
	}


	*dPrice=0.0;

	//-------MonteCarlo----------------------------------------------------------------
//	dSwapPV = amortSwap_PV(lToday, cYCname, iNCpns, lCpnDates, dCpns);
	df_start = swp_f_df(lToday, xlStartDate, cYCname);
	flag = -1;
    if (srtCallPut == SRT_CALL) flag=1;
	for (i=3; i<xlNumSim; i++)
	{
//		if(iPayoffMethod==5)
//		{
			err = GetAmortizedSwaptionPayoff5(dCopulaCube,
											srtFreq,
											srtBasis,
											dCoverages,
											iNCoInitSwaps,
//											dCoInitStrikes,
											dStrikes,
											dCoInitNotionals,
											i,
											&dPayoff,
											flag);
			if(err)
			{
				goto FREE_RETURN;
			}

//		*dPrice += (dSwapPV + df_start * dPayoff);

//		*dPrice += df_start * DMAX(0, flag * dPayoff);
//		*dPrice += df_start * DMAX(0, flag * dPayoff) - flag * (dSwapPV + df_start * dPayoff);

		*dPrice += (DMAX(0,flag * dPayoff) - *dPrice)/(i-2);
//		*dPrice += df_start * (flag * dPayoff - *dPrice)/(i-2);
//		*dPrice += df_start * flag * dPayoff + (df_start * dSimuSwapPV - dSwapPV);
	}

	*dPrice = df_start * (*dPrice);


//----------------Frees the memory and return-----------------------
FREE_RETURN:

	if(dPPDx) free_dmatrix(dPPDx,0,iNDates-1,0,2*nClass);
	if(dPPDy) free_dmatrix(dPPDy,0,iNDates-1,0,2*nClass);
	if(dCPDy) free_dmatrix(dCPDy,0,iNDates-1,0,2*nClass);

	if(dPPDShiftedLogx) free_dmatrix(dPPDShiftedLogx,0,iNDates-1,0,2*nClass);

	if(dMean_v) free_dvector(dMean_v,0,iNDates-1);
	if(dAvg) free_dvector(dAvg,0,iNDates-1);
	if(dCopulaCube) free_dmatrix(dCopulaCube,0,xlNumSim,0,iNDates-1);
	if(dCopulaCubeShiftedLog) free_dmatrix(dCopulaCubeShiftedLog,0,xlNumSim,0,iNDates-1);

	if(dCoverages) free(dCoverages);
	if(lStartDates) free(lStartDates);
	if(lEndDates) free(lEndDates);
	if(lPayDates) free(lPayDates);

	if(lCoInitSwapsEndDates) free(lCoInitSwapsEndDates);
	if(dCoInitStrikes) free(dCoInitStrikes);
	if(dCoInitNotionals) free(dCoInitNotionals);

	if(dStrikes) free_dvector(dStrikes, 0, iNDates-1);
	if(dShifts) free_dvector(dShifts, 0, iNDates-1);
	if(dMeans) free_dvector(dMeans, 0, iNDates-1);
	if(dSpreads) free(dSpreads);

	if(lCpnDates) free(lCpnDates);
	if(dCpns) free(dCpns);

	return err;
}
//--------------------------------------------------------------------------------------
//---------------------SABR Model : End Of Main Function--------------------------------
//--------------------------------------------------------------------------------------
*/