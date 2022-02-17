
#include "bsconvadjustrep.h"
#include "bsmodel.h"
#include "correlmanager.h"
#include "irindex.h"
#include "fromto.h"
#include "gaussian.h"
#include "gpbase/datestrip.h"

using ARM::ARM_DateStrip;


ARM_CMSVolInterSpline::ARM_CMSVolInterSpline(const ARM_Vector& x, const ARM_Vector& y)
{
	if(x.size() < 3)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_OPERATION, "ARM_CMSVolInterSpline : need at least 3 points to build spline");
	if(x.size() != y.size())
		throw Exception(__LINE__, __FILE__, ERR_INVALID_OPERATION, "ARM_CMSVolInterSpline : abscisses and ordinates should have same size");

	itsX = x;
	itsY = y;

	buildSpline();

	itsIsInit = true;
}

ARM_CMSVolInterSpline::ARM_CMSVolInterSpline(const ARM_CMSVolInterSpline& rhs) : itsX(rhs.itsX),
	itsY(rhs.itsY), itsMK(rhs.itsMK), itsAK(rhs.itsAK), itsBK(rhs.itsBK),
	itsIsInit(rhs.itsIsInit)
{
}

void ARM_CMSVolInterSpline::buildSpline()
{
	int k, size = itsX.size();

	double * fk = new double [size];
	double * hk = new double [size];

	for(k = 0; k < size - 1; k++) hk[k] = itsX[k+1] - itsX[k];
	for(k = 1; k < size - 1; k++) fk[k] = 6. * ((itsY[k+1] - itsY[k]) / hk[k] - (itsY[k] - itsY[k-1]) / hk[k-1]);

	findMK(fk, hk);

	itsAK.Resize(size);
	itsBK.Resize(size);

	for(k = 0; k < size-1; k++) itsAK[k] = (itsY[k+1] - itsY[k]) / hk[k] - hk[k] * (itsMK[k+1] - itsMK[k]) / 6.;

	for(k = 0; k < size; k++) itsBK[k] = itsY[k] - itsMK[k] * hk[k] * hk[k] / 6.;

	delete [] hk;
	delete [] fk;
}

void ARM_CMSVolInterSpline::findMK(double * fk, double * hk)
{
	int k, n = itsX.size()-1;

	itsMK.Resize(n+1);

	double * beta = new double [n];
	double * gprime = new double [n];

	beta[n-1]	= 2. * (hk[n-1] + hk[n-2]);
	gprime[n-1] = fk[n-1];

	for(k = n-2; k > 0; k--)
	{
		beta[k]		= 2. * (hk[k] + hk[k-1]) - hk[k] * hk[k-1] / beta[k+1];
		gprime[k]	= fk[k] - hk[k] * gprime[k+1] / beta[k+1];		
	}

	itsMK[1] = gprime[1] / beta[1];
	for(k = 2; k < n; k++) itsMK[k] = (gprime[k] - hk[k-1] * itsMK[k-1]) / beta[k];

	itsMK[0] = 0.;
	itsMK[n] = 0.;

	delete [] beta;
	delete [] gprime;
}

ARM_BSConvAdjustRep::ARM_BSConvAdjustRep(ARM_Model * UsedModel, const ARM_VolLInterpol& swoptVolCurve, 
										 const ARM_Vector& stddev, const ARM_VolLInterpol& MR, int NbPointsForRepliq,
										 bool fullRepliq, double upperProba, double lowerProba) : ARM_BSConvAdjust()
{
	SetName(ARM_BSCONVADJUSTREP);

	if(UsedModel == NULL)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_OPERATION,
                             " ARM_BSConvAdjustRep : need a model for replication");
	}

	itsUsedModel = UsedModel;

	itsCvxCMSVol = swoptVolCurve;

	int i, k, nbLines = swoptVolCurve.GetVolatilities()->GetNumLines(), nbCols = swoptVolCurve.GetVolatilities()->GetNumCols();

	for(i = 0; i < nbLines; i++)
		for(k = 0; k < nbCols; k++)
			itsCvxCMSVol.GetVolatilities()->Elt(i,k) = -1.;

	if(stddev.size() == 0)
	{
		itsStdDev.Resize(7);
		
		itsStdDev[0] = -5;
		itsStdDev[1] = -2;
		itsStdDev[2] = -1;
		itsStdDev[3] = 0;
		itsStdDev[4] = 1;
		itsStdDev[5] = 2;
		itsStdDev[6] = 3;

		itsIdxATM = 3;
	}
	else
	{
		if(stddev.size() < 3)
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_OPERATION,
                            " ARM_BSConvAdjustRep : std dev for interpolation should have at least 3 points");
		}

		itsStdDev = stddev;
		if(itsStdDev.find(0) == -1) itsStdDev.push_back(0.);
		itsStdDev.Sort();
		itsIdxATM = itsStdDev.find(0.);
	}

	itsCMSVol.resize(itsStdDev.size());

	for(k = 0; k < itsStdDev.size(); k++)
	{
		itsCMSVol[k] = itsCvxCMSVol;
	}

	itsNbPointsForRepliq = NbPointsForRepliq;
	itsMR = MR;
	itsFullRepliq = fullRepliq;
	itsUpperProba = upperProba;
	itsLowerProba = lowerProba;
}

ARM_BSConvAdjustRep::ARM_BSConvAdjustRep(const ARM_BSConvAdjustRep& rhs) : ARM_BSConvAdjust(rhs),
	itsNbPointsForRepliq(rhs.itsNbPointsForRepliq), itsIdxATM(rhs.itsIdxATM), 
	itsCvxCMSVol(rhs.itsCvxCMSVol), itsStdDev(rhs.itsStdDev), 
	itsCMSVol(rhs.itsCMSVol), itsMR(rhs.itsMR), itsFullRepliq(rhs.itsFullRepliq),
	itsUpperProba(rhs.itsUpperProba), itsLowerProba(rhs.itsLowerProba)
{
	SetName(ARM_BSCONVADJUSTREP);

	itsUsedModel = rhs.itsUsedModel;
}

ARM_BSConvAdjustRep& ARM_BSConvAdjustRep::operator = (const ARM_BSConvAdjustRep& rhs)
{
	if(&rhs != this)
	{
		this->~ARM_BSConvAdjustRep();

		new (this) ARM_BSConvAdjustRep (rhs);
	}

	return *this;
}

ARM_BSConvAdjustRep::~ARM_BSConvAdjustRep()
{
}

double ARM_BSConvAdjustRep::NaturalAdjstCMS(ARM_Model * Model, const NaturalAdjustData& Input, 
											StoreFwdRateInfo *  StoreInfo)
{
	double Tenor = ROUND((Input.End-Input.Start)/K_YEAR_LEN);
	int tInf = itsCvxCMSVol.GetStrikes()->LookupPrevIndex(Tenor);
	int tSup = fabs((*(itsCvxCMSVol.GetStrikes()))[tInf] - Tenor) < K_DOUBLE_TOL ? tInf : tInf == itsCvxCMSVol.GetStrikes()->size() - 1 ? tInf : tInf+1;

	if((*itsCvxCMSVol.GetVolatilities()).Elt(0,tInf) < 0.)
	{
		doColReplication(tInf);
	}

	if((*itsCvxCMSVol.GetVolatilities()).Elt(0,tSup) < 0.)
	{
		doColReplication(tSup);
	}

	double asof= Model->GetStartDate().GetJulian();
	double resetLag = (Input.Reset-asof)/K_YEAR_LEN;
    double swapVol	= itsCvxCMSVol.ComputeVol(resetLag, Tenor);

	int dayCount= Input.Daycount;

	double Term_se  = Tenor;
 	double theta	= 1.0/Input.Frequency;
	double N		= ROUND((Input.End-Input.Start)/K_YEAR_LEN)*Input.Frequency;
	double rawSwap	= Input.RawSwap/100.0;
   	int YieldDecompFreq = Input.YieldDecompFreq;
	double Margin		= Input.Margin/100.0;
	double payLagAdj	= Input.PayLagAdj;


    double lambda=1.0-(theta*N*rawSwap)/((1.0+theta*rawSwap)*(pow((1.0+theta*rawSwap),N)-1.0));
    double exp1=exp(swapVol*swapVol*resetLag)-1.0;
    double ki = lambda*exp1;

	if (YieldDecompFreq > 1 || Margin != 0.0)
	{
		double convlag = (1+ki);
        
		if(GetSUMMITFormulaUsed()==K_CONV_ADJ_EXP)// 0
        {
            convlag*=payLagAdj;
        }

		double cmsadust = rawSwap * convlag;
		double swapRateDcpSpr = FromRateToRate((rawSwap+Margin)*100.0, 1.0,
                                      K_COMP_ANNUAL, YieldDecompFreq)/100.0;

		double swapRateDcpSprAdjust = FromRateToRate((rawSwap*convlag+Margin)*100.0, 1.0,
                                      K_COMP_ANNUAL, YieldDecompFreq)/100.0;

		double xi = 0.0;
		if (YieldDecompFreq > 1)
		{
			xi = 0.5*pow(swapVol*(cmsadust+Margin),2.0);
			xi *=(1.0/YieldDecompFreq-1)*pow(1.0+cmsadust+Margin,1.0/YieldDecompFreq-2.0)*resetLag;
			xi /= (YieldDecompFreq*(pow(1.0+cmsadust+Margin,1.0/YieldDecompFreq)-1.0));
		}

		ki = swapRateDcpSprAdjust/swapRateDcpSpr*exp(xi)/payLagAdj-1.0;
	}

	return ki+1.0;	
}

double ARM_BSConvAdjustRep::ComputeCMSVolatility(ARM_Model * Model, const NaturalAdjustData& Input, 
												 StoreFwdRateInfo * info)
{
	double Tenor = ROUND(Input.End-Input.Start)/K_YEAR_LEN;

	if(itsFullRepliq == true)
	{
		double asof = itsUsedModel->GetStartDate().GetJulian();
		ARM_Currency * ccy = itsUsedModel->GetZeroCurve()->GetCurrencyUnit();

		ARM_Date StartDate(Input.Reset);
		ARM_Date EndDate(StartDate);

		Tenor < 1. ? EndDate.AddMonths(Tenor*12) : EndDate.AddYears((int)Tenor);

		double zcPay = itsUsedModel->GetZeroCurve()->DiscountPrice(StartDate);

		// Réplication du forward
		CMSReplicator func(itsUsedModel, StartDate, EndDate, ccy->GetFixedDayCount(), ccy->GetFixedPayFreq(), 1, itsMR.ComputeVol((Input.Reset - asof)/K_YEAR_LEN,Tenor), itsUpperProba, itsLowerProba);
		double strike = (fabs(Input.Strike) < K_DOUBLE_TOL || fabs(Input.Fwd) < K_DOUBLE_TOL) ? func.GetFwdRate() : Input.Strike / 100.;
		double fwd = (fabs(Input.Strike) < K_DOUBLE_TOL || fabs(Input.Fwd) < K_DOUBLE_TOL) ? func.GetFwdRate() : Input.Fwd / 100.;
		double xlim = func.limitFromStrike(strike);

		if(strike > fwd)
		{
			func.setCoeffs(1., -strike);
			func.setType(1);
			double opt = LegendreInt(func, xlim, func.GetUBound(), itsNbPointsForRepliq) + func.extraTerms(xlim, strike);
			
			return ImpCMSVol(fwd, (Input.Reset - asof) / K_YEAR_LEN, strike, 1, opt/zcPay);
		}
		else
		{
			func.setCoeffs(-1., strike);
			func.setType(-1);
			double opt = LegendreInt(func, func.GetLBound(), xlim, itsNbPointsForRepliq) + func.extraTerms(xlim, strike);
			
			return ImpCMSVol(fwd, (Input.Reset - asof) / K_YEAR_LEN, strike, -1, opt/zcPay);
		}
	}

	int tInf = itsCvxCMSVol.GetStrikes()->LookupPrevIndex(Tenor);
	int tSup = tInf == itsCvxCMSVol.GetStrikes()->size() - 1 ? tInf : tInf+1;

	if((*itsCvxCMSVol.GetVolatilities()).Elt(0,tInf) < 0.)
	{
		doColReplication(tInf);
	}

	if((*itsCvxCMSVol.GetVolatilities()).Elt(0,tSup) < 0.)
	{
		doColReplication(tSup);
	}

	double asof= Model->GetStartDate().GetJulian();
	double resetLag = (Input.Reset-asof)/K_YEAR_LEN;

	if(fabs(Input.Strike) < K_DOUBLE_TOL || fabs(Input.Fwd) < K_DOUBLE_TOL)
		return itsCMSVol[itsIdxATM].ComputeVol(resetLag, Tenor);

	ARM_Vector vars(itsStdDev.size());

	for(int k = 0; k < itsStdDev.size(); k++)
	{
		vars[k] = SQR(itsCMSVol[k].ComputeVol(resetLag, Tenor)) * resetLag;
	}

	double stddev = log(Input.Strike/Input.Fwd) / sqrt(vars[itsIdxATM]);

	// interpolation sur les variances
	ARM_CMSVolInterSpline spline(itsStdDev, vars);

	double var = spline(stddev);

	return sqrt(var / resetLag);
}

void ARM_BSConvAdjustRep::doColReplication(int col)
{
	double asof = itsUsedModel->GetStartDate().GetJulian();
	double Tenor = (*itsCvxCMSVol.GetStrikes())[col];
	int size = itsCvxCMSVol.GetExpiryTerms()->size();

	ARM_Currency * ccy = itsUsedModel->GetZeroCurve()->GetCurrencyUnit();
	double theta = 1. / (double)ccy->GetFixedPayFreq();
	double N = ROUND(Tenor)*ccy->GetFixedPayFreq();

	for(int k = 0; k < size; k++)
	{
		
		double xlim, call_part, put_part, opt, strike, volATM, mat;
		mat = (*itsCvxCMSVol.GetExpiryTerms())[k];

		ARM_Date StartDate(asof + mat*K_YEAR_LEN);
		ARM_Date EndDate(StartDate);

		Tenor < 1. ? EndDate.AddMonths(Tenor*12) : EndDate.AddYears((int)Tenor);

		double zcPay = itsUsedModel->GetZeroCurve()->DiscountPrice(StartDate);

		// Réplication du forward
		CMSReplicator func(itsUsedModel, StartDate, EndDate, ccy->GetFixedDayCount(), ccy->GetFixedPayFreq(), 1, itsMR.ComputeVol(mat, Tenor), itsUpperProba, itsLowerProba);
		xlim = func.limitFromStrike(func.GetFwdRate());
		func.setCoeffs(1.,0.);
		func.setType(1);
		call_part = LegendreInt(func, xlim, func.GetUBound(), itsNbPointsForRepliq) + func.extraTerms(xlim, func.GetFwdRate());
		func.setType(-1);
		put_part = LegendreInt(func, func.GetLBound(), xlim, itsNbPointsForRepliq) + func.extraTerms(xlim, func.GetFwdRate());

		double cmsByReplic = (call_part + put_part) / zcPay;

		// la vol équivalente pour l'ajustement de convexité
		EquivConvAdjVol func2(func.GetFwdRate(), theta, N, mat);
		double best;
		double equivol = brentSolve(func2, cmsByReplic, 1e-6,1.,1e-8, 100, 0, &best);

		itsCvxCMSVol.GetVolatilities()->Elt(k, col) = equivol;

		if(itsFullRepliq == true) continue;

		// Réplication de la vol monnaie
		func.setCoeffs(-1., cmsByReplic);
		xlim = func.limitFromStrike(cmsByReplic);
		func.setType(-1);
		opt = LegendreInt(func, func.GetLBound(), xlim, itsNbPointsForRepliq) + func.extraTerms(xlim, cmsByReplic);

		volATM = ImpCMSVol(cmsByReplic, mat, cmsByReplic, -1, opt/zcPay);

		itsCMSVol[itsIdxATM].GetVolatilities()->Elt(k, col) = volATM;

		// Réplication des options strikées
		for(int i = 0; i < (int)itsStdDev.size(); i++)
		{
			if(i == itsIdxATM) continue;

			strike = cmsByReplic * exp(itsStdDev[i] * volATM * sqrt(mat));
			xlim = func.limitFromStrike(strike);

			if(i < itsIdxATM)
			{
				func.setCoeffs(-1., strike);
				func.setType(-1);
				opt = LegendreInt(func, func.GetLBound(), xlim, itsNbPointsForRepliq) + func.extraTerms(xlim, strike);

				itsCMSVol[i].GetVolatilities()->Elt(k, col) = ImpCMSVol(cmsByReplic, mat, strike, -1, opt/zcPay);
			}
			else if(i > itsIdxATM)
			{
				func.setCoeffs(1., -strike);
				func.setType(1);
				opt = LegendreInt(func, xlim, func.GetUBound(),itsNbPointsForRepliq) + func.extraTerms(xlim, strike);

				itsCMSVol[i].GetVolatilities()->Elt(k, col) = ImpCMSVol(cmsByReplic, mat, strike, 1, opt/zcPay);
			}
		}
	}
}

ARM_BSConvAdjustRep::CMSReplicator::CMSReplicator(ARM_Model * model, ARM_Date& StartDate, ARM_Date& EndDate, 
								   int basis, int freq, int CallPut, double mr,
								   double upperProba, double lowerProba, double strike)
{
	
	ARM_DateStrip datestrip(StartDate, EndDate, freq, basis, 
		model->GetDomesticCcy(), K_MOD_FOLLOWING, K_UNADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, freq,
		0, model->GetDomesticCcy(), K_ADVANCE, K_ARREARS, false);

	int k, size = datestrip.size();

	double zcPay = model->GetZeroCurve()->DiscountPrice(StartDate);
	double asof = model->GetStartDate().GetJulian();

	double dt = (StartDate.GetJulian() - asof)/K_YEAR_LEN;

	itsDayCount.Resize(size);
	itsZCFwd.Resize(size);
	itsvarZC.Resize(size);
	itsracVarZC.Resize(size);
	itsTheta.Resize(size);

	itsAnnuity = 0.;
	for(k = 0; k < size; k++)
	{
		itsDayCount[k] = (*(datestrip.GetInterestTerms()))[k];

		itsZCFwd[k] = model->GetZeroCurve()->DiscountPrice(((*(datestrip.GetFlowEndDates()))[k]-asof)/K_YEAR_LEN) / zcPay;
		itsAnnuity += itsDayCount[k] * itsZCFwd[k] * zcPay;
		itsTheta[k] = ((*(datestrip.GetFlowEndDates()))[k] - StartDate.GetJulian()) / K_YEAR_LEN;
	}

	itsFwdRate = (1. - itsZCFwd[size-1]) / (itsAnnuity / zcPay);
	itsStartDate = (StartDate.GetJulian() - asof)/K_YEAR_LEN;
	itsEndDate = (EndDate.GetJulian() - asof)/K_YEAR_LEN;
	itsResetDate = itsStartDate;
	itsUsedModel = model;
	itsMR = mr;

	double K = fabs(strike) < K_DOUBLE_TOL ? itsFwdRate : strike;
	double mktPrice = itsUsedModel->EuroSwaption(itsResetDate, itsStartDate, itsEndDate, -1, itsFwdRate, K, itsAnnuity); 

	double vol = GetHWVol(mktPrice / zcPay, K);

	for(k = 0; k < size; k++)
	{
		itsvarZC[k] = SQR(vol * itsTheta[k] * KSI(mr * itsTheta[k])) * itsResetDate * KSI(mr*itsResetDate);
		itsracVarZC[k] = sqrt(itsvarZC[k]);
	}

	double uBound = BoundFromProba(upperProba);
	double lBound = BoundFromProba(lowerProba);

	itsUBound = uBound;
	itsLBound = lBound;
}

ARM_BSConvAdjustRep::CMSReplicator::ReplicationCoeffs::ReplicationCoeffs(double x, const ARM_Vector& zcFwd, 
																		 const ARM_Vector& varZC, const ARM_Vector& racVarZC, 
																		 const ARM_Vector& dcf, double a, double b)
{
	int size = varZC.size();

	double zc, dZc, d2Zc;
	double annuity = 0., dAnnuity = 0., d2Annuity = 0.;

	int k;
	for(k = 0 ; k < size ; k++)
	{
		zc = dcf[k] * zcFwd[k] * exp( - 0.5 * varZC[k] - racVarZC[k] * x);
		annuity += zc;
		dAnnuity -= racVarZC[k]*zc;
		d2Annuity += varZC[k]*zc;
	}
	k--;
	zc = zcFwd[k] * exp( - 0.5 * varZC[k] - racVarZC[k] * x);
	dZc = - racVarZC[k] * zc;
	d2Zc = varZC[k] * zc;

	double swapRate, dSwapRate, d2SwapRate;
	double beta, dBeta, d2Beta;

	swapRate = (1. - zc) / annuity;
	dSwapRate = swapRate * ( - dZc / (1. - zc) - dAnnuity / annuity );
	d2SwapRate = dSwapRate * ( - dZc / (1. - zc) - dAnnuity / annuity ) 
		+ swapRate * (	- d2Zc / (1. - zc) - dZc * dZc / (1. - zc) / (1. - zc) 
						- d2Annuity / annuity + dAnnuity * dAnnuity / annuity / annuity ) ;

	beta = (a * swapRate + b) / annuity;
	dBeta = a * dSwapRate / annuity - (a * swapRate + b) * dAnnuity / annuity / annuity ;
	d2Beta	= a * d2SwapRate / annuity - a * dSwapRate * dAnnuity / annuity / annuity 
			- (a * dSwapRate * dAnnuity + (a * swapRate + b) * d2Annuity) / annuity / annuity 
				+ 2.* (a * swapRate + b) * dAnnuity * dAnnuity / annuity / annuity / annuity;

	itsSwapRate			= swapRate;
	itsDerivSwapRate	= dSwapRate;
	itsBeta				= beta;
	itsRatio			= dBeta / dSwapRate;
	itsDerivRatio		= d2Beta / dSwapRate - dBeta * d2SwapRate / dSwapRate / dSwapRate; 
}

double ARM_BSConvAdjustRep::CMSReplicator::GetHWVol(double target, double strike)
{
	double hwvol = 0.;

	ARM_BSConvAdjustRep::CMSReplicator::HWSwopt func(this, itsMR, itsResetDate, strike, -1);
	double best;

	hwvol = brentSolve(func, target, 1e-8, 0.1,1e-8,100,0,&best);

	return hwvol;
}

double ARM_BSConvAdjustRep::CMSReplicator::operator () (double x) const
{
	CMSReplicator::ReplicationCoeffs coeffs(x, itsZCFwd, itsvarZC, itsracVarZC, itsDayCount, itsAlpha, itsBeta);

	return coeffs.itsDerivRatio * itsUsedModel->EuroSwaption(itsResetDate, itsStartDate, itsEndDate, itsCallPut, itsFwdRate, coeffs.itsSwapRate, itsAnnuity); 
}

double ARM_BSConvAdjustRep::CMSReplicator::extraTerms(double xlim, double strike) const
{
	CMSReplicator::ReplicationCoeffs coeffs(xlim, itsZCFwd, itsvarZC, itsracVarZC, itsDayCount, itsAlpha, itsBeta);
	
	double opt = itsUsedModel->EuroSwaption(itsResetDate, itsStartDate, itsEndDate, itsCallPut, itsFwdRate, strike, itsAnnuity); 

	double bin = itsUsedModel->EuroSwaption(itsResetDate, itsStartDate, itsEndDate, itsCallPut, itsFwdRate, strike + 0.0001, itsAnnuity)
			   - itsUsedModel->EuroSwaption(itsResetDate, itsStartDate, itsEndDate, itsCallPut, itsFwdRate, strike - 0.0001, itsAnnuity); 

	bin = fabs(bin) / 2. / 0.0001;

	return coeffs.itsBeta * bin + itsCallPut * coeffs.itsRatio * opt;
}

double ARM_BSConvAdjustRep::CMSReplicator::limitFromStrike(double strike)
{
	CMSReplicator::SwapFromFactor func(this, &itsvarZC, &itsracVarZC);

	double best;
	return brentSolve(func, strike, itsLBound, itsUBound, 1e-8, 100, 0, &best); 
}

double ARM_BSConvAdjustRep::CMSReplicator::SwapFromFactor::operator () (double x) const
{
	double annuity = 0., zc;
	int k, size = _this->itsZCFwd.size();

	for(k = 0 ; k < size ; k++)
	{
		zc = _this->itsZCFwd[k] * exp( - 0.5 * (*itsVars)[k] - (*itsracVars)[k] * x);
		annuity += zc * _this->itsDayCount[k];
	}

	return (1. - zc) / annuity;
}

ARM_BSConvAdjustRep::CMSReplicator::HWSwopt::HWSwopt(CMSReplicator * replicator, double MR, double Mat, double K, int CallPut)
{
	_this = replicator;
	itsMR = MR;
	itsMat = Mat;
	itsStrike = K;
	itsCallPut = CallPut;

	int k, size = _this->itsZCFwd.size();
	itsFlows.Resize(size);
	itsVars.Resize(size);
	itsracVars.Resize(size);

	for(k = 0; k < size; k++)
		itsFlows[k] = - itsStrike * _this->itsDayCount[k];
	itsFlows[size-1] -= 1.;
}

double ARM_BSConvAdjustRep::CMSReplicator::HWSwopt::operator () (double x) const
{
	double result = 0.;
	int k, size = _this->itsZCFwd.size();

	for(k = 0; k < size; k++)
	{
		itsVars[k] = SQR(x * _this->itsTheta[k] * KSI(itsMR * _this->itsTheta[k])) * itsMat * KSI(itsMR * itsMat);
		itsracVars[k] = sqrt(itsVars[k]);
	}

	ARM_BSConvAdjustRep::CMSReplicator::SwapFromFactor func(_this, &itsVars, &itsracVars);
	double best;
	double xlim = brentSolve(func, itsStrike, _this->itsLBound, _this->itsUBound, 1e-8, 100, 0, &best);

	for(k = 0; k < size; k++)
		result += itsCallPut * _this->itsZCFwd[k] * itsFlows[k] * cdfNormal(- itsCallPut * (xlim + itsracVars[k]));

	result += itsCallPut * cdfNormal(- itsCallPut * xlim);

	return result;
}

double ARM_BSConvAdjustRep::CMSReplicator::BoundFromProba(double proba)
{
	CMSReplicator::EquivStrikeForProba func(itsResetDate, itsStartDate, itsEndDate, itsFwdRate, itsAnnuity, itsUsedModel);

	double best;
	double strike = brentSolve(func, proba, 1e-6, 1., 1e-6, 100, 0, &best);

	itsUBound = 20.;
	itsLBound = -20.;

	return limitFromStrike(strike);
}

double ARM_BSConvAdjustRep::ImpCMSVol(double fwd, double mat, double strike, int CallPut, double price)
{
	EquivCMSVol func(fwd, mat, strike, CallPut);
	double best;
	return brentSolve(func, price, 1e-6,2.,1e-8,100,0,&best);
}