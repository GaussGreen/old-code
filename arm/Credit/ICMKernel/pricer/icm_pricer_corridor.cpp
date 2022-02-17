#include "firsttoinc.h"

#include "ICMKernel/pricer/icm_pricer_corridor.h"

#include "ICMKernel/glob/icm_maths.h"
#include "ICMKernel/util/icm_utils.h"
#include "ICMKernel/glob/icm_enums.h"
#include "ICMKernel/crv/icm_fixing_curve.h"
#include "ICMKernel/util/icm_integrator.h"
#include "ICMKernel/inst/icm_corridorleg.h"
#include "ICMKernel/crv/icm_defaultcurve.h"
#include "ARMKernel/crv/volflat.h"
#include "ICMKernel/glob/icm_correlation.h"
#include "ICMKernel/inst/icm_credit_index.h"


double ICM_Pricer_Corridor::ComputePrice(qCMPMETH measure)
{
	// bidouille temporaire 
	if ( measure == qCMPPRICE) measure = qCMPFEELEGPV;
	if (getFlg(measure)) return getValue(measure);
	ICM_CorridorLeg*	CorridorLeg=NULL;
	CorridorLeg	=	((ICM_CorridorLeg*) GetSecurity());
	ARM_Vector* resetsched = CorridorLeg->GenResetSched();
	double result = CptPrice_Corridor(measure);
	return (result);
}
void ICM_Pricer_Corridor::Set(ARM_Security *sec, ARM_Object *mod, const ICM_Parameters& parameters , const ARM_Date& valdate )
{
	ICM_Pricer::Set(sec,mod,parameters,&valdate);
	ICM_CorridorLeg* leg = dynamic_cast<ICM_CorridorLeg*>(sec);

	if (leg) {
		string name;
		string issuer = leg->GetCreditIndex()->GetLabels()[0];
		string iIRIndex = leg->GetIRIndex()->GetIndexName();

		name=GetZeroCurveName(leg->GetCurrencyUnit()->GetCcyName(),GetAsOfDate());
		itsZeroCurve = (ARM_ZeroCurve*)((ICM_MktDataMng*)mod)->find(name);
		if (!itsZeroCurve) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Corridor:No zero curve "<<name); 

		name = GetDefaultCurveName(issuer,GetAsOfDate());
		itsDefCurveIndex = (ICM_DefaultCurve*)((ICM_MktDataMng*)mod)->find(name);
		if (!itsDefCurveIndex) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Corridor:No def curve "<<name); 
		
		// fixing curve IR
		name = GetFixingCurveName(iIRIndex, leg->GetIRIndex()->GetCurrencyUnit()->GetCcyName());
		itsFixingCurveIR = (ICM_Fixing_Curve*) ((ICM_MktDataMng*)mod)->find(name);
		if (!itsFixingCurveIR) 
		//	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Corridor : No Fixing curve "<<name); 
		itsFixingCurveIR = new 	ICM_Fixing_Curve;	

		// fixing curve Credit 
		name = GetFixingCurveName(issuer, leg->GetCreditIndex()->GetCurrencyUnit()->GetCcyName());
		itsFixingCurveCredit = (ICM_Fixing_Curve*) ((ICM_MktDataMng*)mod)->find(name);
		if (!itsFixingCurveCredit) 
			itsFixingCurveCredit = new ICM_Fixing_Curve;
		//	ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Corridor : No Fixing curve "<<name); 

		// Vol taux.
		name = GetVolCurveName(iIRIndex,leg->GetIRIndex()->GetCurrencyUnit()->GetCcyName(),GetAsOfDate());
		itsIRVolatility = (ARM_VolCurve*)((ICM_MktDataMng*)mod)->find(name);
		if (!itsIRVolatility) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Corridor:No IR vol "<<name); 

		// index Vol (Cube)
		name = issuer;
		name = GetVolCurveName(name,leg->GetCreditIndex()->GetCurrencyUnit()->GetCcyName(),GetAsOfDate());
		itsIdxVolatility = (ARM_VolCurve*)((ICM_MktDataMng*)mod)->find(name);
		if (!itsIdxVolatility) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Corridor:No Credit  vol "<<name); 
		
		// Correlation
		name = leg->GetSingleName();
		name = GetCorrelationName(name, iIRIndex, issuer, GetAsOfDate());
		itsCorrelation = (ICM_Correlation*)((ICM_MktDataMng*)mod)->find(name);
		if (!itsCorrelation) {
			string name2 = GetCorrelationName(leg->GetSingleName(),"", "", GetAsOfDate());
			itsCorrelation = (ICM_Correlation*)((ICM_MktDataMng*)mod)->find(name2);
			if (!itsCorrelation) {
				ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Corridor:No Correl "<<name<<" nor "<< name2); 
			}
		}

	}
	else {
		ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Pricer_Corridor::Set : security is not a corridorLeg") ; 
	}
}
double ICM_Pricer_Corridor::Accrued(void) 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
}
double ICM_Pricer_Corridor::FeeLegPV(void) 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
}
double ICM_Pricer_Corridor::DefLegPV(void) 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
}
double ICM_Pricer_Corridor::Compute_Fwd_Spread(const class ARM_Date &,const class ARM_Date &, double& dur)
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
}
double ICM_Pricer_Corridor::ComputeDuration(void) 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
}
double ICM_Pricer_Corridor::ComputeSpread(const double &) 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
}
double ICM_Pricer_Corridor::ComputeImpliedVol(const double &) 
{
	ICMTHROW(ERR_INVALID_ARGUMENT,"Not Implemented") ; 
}

// -------------------------------------------------------------------
// Floating Leg
// -------------------------------------------------------------------

double	ICM_Pricer_Corridor::CorridorSubPeriodPVFloat(double ResetDate,
												double BarrierDown,
												double BarrierUp,
												double VolIrCurve,
												double FwdStartDate,
												double FwdEndDate,
												double MaturityIdxYF,
												double FixingCredit,
												double& voldown,
												double& volup,
												double& FwdSpreadcvx)
{
	double	AsOfDate_Julian	=	GetAsOfDate().GetJulian();

	double	FwdSpread1 = 0.;
	double	YFResetDate	=	(ResetDate-GetAsOfDate().GetJulian())/365.;
	double	ConvexifiedFwdSpread1 = 0.;
	double	ConvexitySpread1= 0.;

	double VolStrikeDown =0.;
	double VolStrikeUp =0.;
	double VolStrikeDown_Epsilon =0.;
	double VolStrikeUp_Epsilon =0.;
	double VolATM =0.;


	ARM_Date	Date_Exer(ResetDate);
	ARM_Date	MatuIndex;

	double	T_FwdStart_yf=0.,T_MatuIndex_yf=0.,Epsilon = 1.e-2;

	// ------------------------------------------------
	// STEP 1: GET INDEX FORWARD
	// ------------------------------------------------

	//Récupération Vol Tx + Correl pour ajustement Float
	//Correl
	double correl = itsCorrelation->GetCorrelation("NONE","NONE");

	if (YFResetDate<=0.)
		{YFResetDate=0.001;}

	T_FwdStart_yf		=	MAX((FwdStartDate - AsOfDate_Julian)/365.,0.);
	T_MatuIndex_yf		=	MAX((FwdEndDate - AsOfDate_Julian)/365.,0.);

	if (FixingCredit == CREDIT_DEFAULT_VALUE)
	{
	FwdSpread1			=	itsDefCurveIndex->FwdSpread((ARM_Date)FwdStartDate, (ARM_Date)FwdEndDate);  //itsDefCurveIndex->FwdSpread_AsIndex(Date_Exer, MatuIndex, Flat_RPV01_1, RPV01_1);

	// Recupération Vol ATM (ajustement convexite)
	if (YFResetDate>=0.){
		VolATM = itsIdxVolatility->ComputeVolatility(YFResetDate,0.,MaturityIdxYF)/100.;}

	//Vol
	ARM_VolFlat v1((ARM_Date)AsOfDate_Julian,VolATM*100.);
	if (T_MatuIndex_yf>=0.){
		ConvexitySpread1		=	itsDefCurveIndex->AjustConvexity(T_FwdStart_yf, T_MatuIndex_yf, FwdSpread1, &v1);}
		ConvexifiedFwdSpread1	=	FwdSpread1 + ConvexitySpread1;
		ConvexifiedFwdSpread1 *=exp(correl*VolATM*VolIrCurve*YFResetDate);
		// Recupération Vol Smile (pricing digitales)
		if (YFResetDate>=0.){
			voldown = itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierDown-ConvexifiedFwdSpread1,MaturityIdxYF)/100.;
			volup = itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierUp-ConvexifiedFwdSpread1,MaturityIdxYF)/100.;
			VolStrikeDown	= itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierDown-ConvexifiedFwdSpread1-Epsilon/2.,MaturityIdxYF)/100.;
			VolStrikeUp	=	itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierUp-ConvexifiedFwdSpread1-Epsilon/2.,MaturityIdxYF)/100.;
			VolStrikeDown_Epsilon	=	itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierDown-ConvexifiedFwdSpread1+Epsilon/2.,MaturityIdxYF)/100.;
			VolStrikeUp_Epsilon	=	itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierUp-ConvexifiedFwdSpread1+Epsilon/2.,MaturityIdxYF)/100.;
		}
	}
	else
	{
		ConvexifiedFwdSpread1=FixingCredit;
		// Recupération Vol Smile (pricing digitales)
		if (YFResetDate>=0.){
			voldown = itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierDown-ConvexifiedFwdSpread1,MaturityIdxYF)/100.;
			volup = itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierUp-ConvexifiedFwdSpread1,MaturityIdxYF)/100.;
			VolStrikeDown	= itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierDown-ConvexifiedFwdSpread1-Epsilon/2.,MaturityIdxYF)/100.;
			VolStrikeUp	=	itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierUp-ConvexifiedFwdSpread1-Epsilon/2.,MaturityIdxYF)/100.;
			VolStrikeDown_Epsilon	=	itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierDown-ConvexifiedFwdSpread1+Epsilon/2.,MaturityIdxYF)/100.;
			VolStrikeUp_Epsilon	=	itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierUp-ConvexifiedFwdSpread1+Epsilon/2.,MaturityIdxYF)/100.;
		}
	}

	FwdSpreadcvx = ConvexifiedFwdSpread1;

	// ------------------------------------------------
	// STEP 5: PRICE
	// ------------------------------------------------

	double D_KA = 0.,D_KB = 0.;
	double C_K = 0.;
	double C_K_Epsilon = 0.;
	double output=0.,NPV=0.;
	
	double spot = ConvexifiedFwdSpread1;

	NPV = CorridorPriceBS(YFResetDate, /*double& yf*/
					   FwdSpreadcvx, /*double& spot,*/
					   BarrierDown, /*double& K1,*/
					   BarrierUp, /*double& K2,*/
					   VolStrikeDown, /*double& volK1,*/
					   VolStrikeDown_Epsilon,/*double& vol_epsilonk1,*/
					   VolStrikeUp, /*double& volK2,*/
					   VolStrikeUp_Epsilon, /*double& vol_epsilonk2,*/
					   Epsilon/*double& epsilon*/);

	return (NPV);
}
 
// -------------------------------------------------------------------
// Fix Leg
// -------------------------------------------------------------------

double	ICM_Pricer_Corridor::CorridorSubPeriodPV(double ResetDate,
												  double BarrierDown,
												  double BarrierUp,
												  double FwdStartDate,
												  double FwdEndDate,
												  double MaturityIdxYF,
												  double FixingCredit,
												  double& voldown,
												  double& volup,
												  double& FwdSpreadcvx)
{
	double	AsOfDate_Julian	=	GetAsOfDate().GetJulian();
	double	FwdSpread1=0.;
	double	YFResetDate	=	(ResetDate-GetAsOfDate().GetJulian())/365.;
	ARM_Date	Date_Exer(ResetDate);
	ARM_Date	MatuIndex;
	double	T_FwdStart_yf=0.,T_MatuIndex_yf=0.,Epsilon = 1.e-2;
	double	ConvexifiedFwdSpread1;
	double	ConvexitySpread1;
	double VolStrikeDown =0.;
	double VolStrikeUp =0.;
	double VolStrikeDown_Epsilon =0.;
	double VolStrikeUp_Epsilon =0.;
	double VolATM =0.;
	
	if (YFResetDate<=0.)
		{YFResetDate=0.001;}

	T_FwdStart_yf		=	MAX((FwdStartDate - AsOfDate_Julian)/365.,0.);
	T_MatuIndex_yf		=	MAX((FwdEndDate - AsOfDate_Julian)/365.,0.);

	if (FixingCredit == CREDIT_DEFAULT_VALUE)
	{
		FwdSpread1			=	itsDefCurveIndex->FwdSpread((ARM_Date)FwdStartDate, (ARM_Date)FwdEndDate);
		// Recupération Vol ATM (ajustement convexite)
		if (YFResetDate>=0.){
		VolATM = itsIdxVolatility->ComputeVolatility(YFResetDate,0.,MaturityIdxYF)/100.;}
		ARM_VolFlat v1((ARM_Date)AsOfDate_Julian,VolATM*100.);
		if (T_MatuIndex_yf>=0.){
			ConvexitySpread1		=	itsDefCurveIndex->AjustConvexity(T_FwdStart_yf, T_MatuIndex_yf, FwdSpread1, &v1);}
			ConvexifiedFwdSpread1	=	FwdSpread1 + ConvexitySpread1;
			if (YFResetDate>=0.)
			{
				voldown = itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierDown-ConvexifiedFwdSpread1,MaturityIdxYF)/100.;
				volup = itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierUp-ConvexifiedFwdSpread1,MaturityIdxYF)/100.;
				VolStrikeDown	= itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierDown-ConvexifiedFwdSpread1-Epsilon/2.,MaturityIdxYF)/100.;
				VolStrikeUp	=	itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierUp-ConvexifiedFwdSpread1-Epsilon/2.,MaturityIdxYF)/100.;
				VolStrikeDown_Epsilon	=	itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierDown-ConvexifiedFwdSpread1+Epsilon/2.,MaturityIdxYF)/100.;
				VolStrikeUp_Epsilon	=	itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierUp-ConvexifiedFwdSpread1+Epsilon/2.,MaturityIdxYF)/100.;
			}
	}
	else 
	{
		ConvexifiedFwdSpread1 = FixingCredit;
		if (YFResetDate>=0.)
		{
			voldown = itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierDown-ConvexifiedFwdSpread1,MaturityIdxYF)/100.;
			volup = itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierUp-ConvexifiedFwdSpread1,MaturityIdxYF)/100.;
			VolStrikeDown	= itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierDown-ConvexifiedFwdSpread1-Epsilon/2.,MaturityIdxYF)/100.;
			VolStrikeUp	=	itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierUp-ConvexifiedFwdSpread1-Epsilon/2.,MaturityIdxYF)/100.;
			VolStrikeDown_Epsilon	=	itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierDown-ConvexifiedFwdSpread1+Epsilon/2.,MaturityIdxYF)/100.;
			VolStrikeUp_Epsilon	=	itsIdxVolatility->ComputeVolatility(YFResetDate, BarrierUp-ConvexifiedFwdSpread1+Epsilon/2.,MaturityIdxYF)/100.;
		}

	}

	FwdSpreadcvx = ConvexifiedFwdSpread1;

	// ------------------------------------------------
	// STEP 5: PRICE
	// ------------------------------------------------

	double D_KA = 0.,D_KB = 0.;
	double C_K = 0.;
	double C_K_Epsilon = 0.;
	double output=0.,NPV=0.;
	
	double spot = ConvexifiedFwdSpread1;

	NPV = CorridorPriceBS(YFResetDate, /*double& yf*/
					   FwdSpreadcvx,
					   BarrierDown, /*double& K1,*/
					   BarrierUp, /*double& K2,*/
					   VolStrikeDown, /*double& volK1,*/
					   VolStrikeDown_Epsilon,/*double& vol_epsilonk1,*/
					   VolStrikeUp, /*double& volK2,*/
					   VolStrikeUp_Epsilon, /*double& vol_epsilonk2,*/
					   Epsilon/*double& epsilon*/);

	return (NPV);
}

// ---------------------------------------------------------------------
// corridor pricing
// ---------------------------------------------------------------------

double	ICM_Pricer_Corridor::CptPrice_Corridor(qCMPMETH measure)
{
	ICM_CorridorLeg* CorridorLeg = (ICM_CorridorLeg*) GetSecurity();
	ICM_Security* security = CorridorLeg->GetCreditInfos();
	security->ComputeYF(GetAsOfDate());
	ARM_Vector* ResetSched= CorridorLeg->GenResetSched();
	vector<double> SubResetSched;
	int	i =0,j=0,NbFlows=0;
	NbFlows = security->GetAccStartDates().GetSize();	
	double	CurrentFlow=0.,NPV	=	0.0;
	double	Duration = 0.0;
	double PV01 = 0.,ProbaRangeFix=0.,SumProbaRangeFix=0.,SumProbaRangeFloat=0.;
	double ProbaRangeFloat=0.;
	double startdate =0.,enddate=0.,paydate=0.,ir_fwdstart=0.,ir_fwdend=0.,ir_reset=0.;
	double rate=0.;	
	ARM_Date Start_;
	ARM_Date Matu_;
	ARM_Date FwdStartDate;
	ARM_Date FwdEndDate;
	double voldown=0.;
	double volup=0.;
	double FwdSpreadcvx=0.;
	double VolIrCurve=0.;
	CorridorLeg->GetResetProbas().clear();
	CorridorLeg->GetResetProbas().resize(NbFlows );
	CorridorLeg->GetVolK1().clear();
	CorridorLeg->GetVolK1().resize(NbFlows);
	CorridorLeg->GetVolK2().clear();
	CorridorLeg->GetVolK2().resize(NbFlows);
	CorridorLeg->GetFwdSpreadAdj().clear();
	CorridorLeg->GetFwdSpreadAdj().resize(NbFlows);
	CorridorLeg->GetFwdSpreadAdj_Float().clear();
	CorridorLeg->GetFwdSpreadAdj_Float().resize(NbFlows);
	CorridorLeg->GetResetProbas_Float().clear();
	CorridorLeg->GetResetProbas_Float().resize(NbFlows);
	double accrued = 0.;
	int subsize = 0;

	ARM_IRIndex* ir_index = CorridorLeg->GetIRIndex();
	ICM_Credit_Index* idx_index = CorridorLeg->GetCreditIndex(); 


	double LevFloat = CorridorLeg->GetLeverageFloatIdx();
	double deltaTPeriod=0.,yt_pay=0.,DP=0.,fwd=0.,deltaT=0.,notional=0.;
	int iPayment0 = security->PaymentIndex(GetAsOfDate()) ;	// index of first flow to be paid (-1 if no flows) 
	int iAccrued = security->PeriodIndex(GetAsOfDate()) ;	// index of the period including AsOf (-1 if no period) 
	double KUp = CorridorLeg->GetRefValue_kup()->CptReferenceValue(startdate);
	double KDown = CorridorLeg->GetRefValue_kdw()->CptReferenceValue(startdate);
	double fixingCredit = CREDIT_DEFAULT_VALUE;

#ifdef _DEBUG
		FILE* pFile = NULL;
		pFile = fopen("c:\\temp\\PricingCORRILeg.txt", "w");
#endif
	//case before asOfDate just for View
	for (i=0; i<iPayment0; i++){
		subsize = CorridorLeg->GetSubResetSchedule()[i].size();
		ARM_Vector vGetResetProbas(subsize);
		ARM_Vector vFwdSpreadAdj(subsize);
		for (j=0; j<subsize; j++)
		{
			double reset = CorridorLeg->GetSubResetSchedule()[i][j];
			vFwdSpreadAdj[j] = 0;
			if (itsFixingCurveCredit->hasValue(ARM_Date(reset))) 
					fixingCredit = itsFixingCurveCredit->getValue(ARM_Date(reset));
			else 
				fixingCredit = 0.;
			vFwdSpreadAdj[j] = fixingCredit;
			 
			
			if (KDown < fixingCredit && fixingCredit < KUp )
				vGetResetProbas[j] = 1;
			else
				vGetResetProbas[j] = 0;
		}
		CorridorLeg->GetResetProbas()[i] = vGetResetProbas;
		CorridorLeg->GetResetProbas_Float()[i] = vGetResetProbas;
		CorridorLeg->GetFwdSpreadAdj()[i] = vFwdSpreadAdj;
		CorridorLeg->GetFwdSpreadAdj_Float()[i] = vFwdSpreadAdj;

	}

	fixingCredit = CREDIT_DEFAULT_VALUE;
	for (i=iPayment0; i<NbFlows; i++)
	{
		fwd =0.; // reinitialisation
		startdate = security->GetAccStartDates().Elt(i);
		enddate = security->GetAccEndDates().Elt(i);
		paydate = security->GetPayDates().Elt(i);
		ir_fwdstart = security->GetFwdStartDates().Elt(i);
		ir_fwdend = security->GetFwdEndDates().Elt(i);
		ir_reset = security->GetResetDates().Elt(i);
		double dAsOf = GetAsOfDate().GetJulian();

		if (ir_reset< dAsOf){
			// getting the ir fixing
			if (itsFixingCurveIR->hasValue(ARM_Date(ir_reset)))
				fwd = itsFixingCurveIR->getValue(ARM_Date(ir_reset));
		}
		double yfir_reset = security->GetYFResetDates().Elt(i);
		VolIrCurve=0.;

		if (yfir_reset>=0.){
			VolIrCurve = itsIRVolatility->ComputeVolatility(yfir_reset,0.,ir_index->GetYearTerm());
			VolIrCurve/=100.;
		}
		deltaTPeriod = CountYears(CorridorLeg->GetDayCount(),startdate,enddate);
		rate = CorridorLeg->GetSpreads()->CptReferenceValue(startdate);
		yt_pay= (paydate-GetAsOfDate().GetJulian())/365.;
		if (yt_pay<0)	{ yt_pay = 0.;}
		DP = itsZeroCurve->DiscountPrice(yt_pay);
		if (ir_fwdend>dAsOf && fwd ==0)
		{
				double date = MAX(GetAsOfDate().GetJulian(),ir_fwdstart);
				fwd = itsZeroCurve->ForwardYield((ARM_Date)date,(ARM_Date)ir_fwdend,-1,CorridorLeg->GetDayCount());
		}
		notional = security->GetNotionals().Elt(i);	
		//Reinit Proba periode
		SumProbaRangeFix=0.,SumProbaRangeFloat=0.;
		subsize = CorridorLeg->GetSubResetSchedule()[i].size();

		ARM_Vector vGetResetProbas(subsize);
		ARM_Vector vVolK1(subsize);
		ARM_Vector vVolK2(subsize);
		ARM_Vector vFwdSpreadAdj(subsize);
		ARM_Vector vFwdSpreadAdj_Float(subsize);
		ARM_Vector vResetProbas_Float(subsize);

		for (j=0; j<subsize; j++)
		{
			voldown=volup=FwdSpreadcvx=0.;
			//double substart = CorridorLeg->GetSubStart()[i][j];
			//double subend = CorridorLeg->GetSubEnd()[i][j];
			double subfwdstart = CorridorLeg->GetSubFwdStart()[i][j];
			double subfwdend = CorridorLeg->GetSubFwdEnd()[i][j];
			double reset = CorridorLeg->GetSubResetSchedule()[i][j];
			//double subpaydate = CorridorLeg->GetSubResetSchedule()[i][j];
			double fixingCredit = CREDIT_DEFAULT_VALUE;

			if (reset< dAsOf){
				// getting the Credit fixing
				if (itsFixingCurveCredit->hasValue(ARM_Date(reset)))
					fixingCredit = itsFixingCurveCredit->getValue(ARM_Date(reset)); // fixingCredit must be in %
			}

#ifdef _DEBUG
			fprintf(pFile, "(i,j) = (%d,%d) \t fixing credit %f  \n", i, j, fixingCredit);
#endif
			if (rate)
			{
		  		ProbaRangeFix =	CorridorSubPeriodPV(reset,
								CorridorLeg->GetRefValue_kdw()->CptReferenceValue(startdate), 
								CorridorLeg->GetRefValue_kup()->CptReferenceValue(startdate),
								subfwdstart,subfwdend,idx_index->GetYearTerm(),fixingCredit,
								voldown,volup,FwdSpreadcvx);

			vGetResetProbas[j] = (ProbaRangeFix);
			vVolK1[j] = (voldown);
			vVolK2[j] = (volup);
			vFwdSpreadAdj[j] = (FwdSpreadcvx);
			}

			if (LevFloat)
			{
				ProbaRangeFloat =	CorridorSubPeriodPVFloat(reset,
								CorridorLeg->GetRefValue_kdw()->CptReferenceValue(startdate), 
								CorridorLeg->GetRefValue_kup()->CptReferenceValue(startdate),
								VolIrCurve,subfwdstart,subfwdend,idx_index->GetYearTerm(),fixingCredit,
								voldown,volup,FwdSpreadcvx);

				if (rate==0.) vVolK1[j] = (voldown);
				if (rate==0.) vVolK2[j] = (volup);
				vFwdSpreadAdj_Float[j] = (FwdSpreadcvx);
				vResetProbas_Float [j] = (ProbaRangeFloat);
			}

			/*
			deltaT = CountYears(CorridorLeg->GetDayCount(),substart,subend);
			SumProbaRangeFix += deltaT*ProbaRangeFix;
			SumProbaRangeFloat += deltaT*ProbaRangeFloat;
			*/

			SumProbaRangeFix += ProbaRangeFix;
			SumProbaRangeFloat += ProbaRangeFloat;

		}	
			
			CorridorLeg->GetResetProbas()[i] = vGetResetProbas;
			CorridorLeg->GetVolK1()[i] = vVolK1;
			CorridorLeg->GetVolK2()[i] = vVolK2;
			CorridorLeg->GetFwdSpreadAdj()[i] = vFwdSpreadAdj;
			CorridorLeg->GetFwdSpreadAdj_Float()[i] = vFwdSpreadAdj_Float;
			CorridorLeg->GetResetProbas_Float()[i] = vResetProbas_Float;

		//gestion de l'accrued

		Duration += deltaTPeriod*DP*SumProbaRangeFix/subsize;

		CurrentFlow	= deltaTPeriod*DP*notional*(rate*SumProbaRangeFix+LevFloat*fwd*SumProbaRangeFloat)/(j*100);

#ifdef _DEBUG

		fprintf(pFile, "i = %d \t  \n", i); 
		fprintf(pFile, "\t rate= %lf \t  SumProbaRangeFix = %lf \n",rate , SumProbaRangeFix); 
		fprintf(pFile, "\t LevFloat= %lf \t  fwd = %lf \t Duration = %lf \n\n",LevFloat , fwd, Duration); 
#endif
		security->GetCouponRates()[i] = (rate*SumProbaRangeFix+fwd*SumProbaRangeFloat)/j; // *100 because /100 dans ICM_Security::FullCouponRate.
		security->GetFwdSpreads()[i] = fwd;

		if (i==iAccrued)
		{
			accrued = security->AccruedCoupon(GetAsOfDate());
			
		}
		
		NPV	+=	CurrentFlow;
	}

#ifdef _DEBUG
	if (pFile) fclose(pFile);
#endif
	// --------------------------
	// Implied Spread Calculation
	// --------------------------
	// to discrimine whether it is Risky or not

	double	Spread = 0;
	double FloatLeg = 0;
	if (Duration)
		Spread	=	10000.0 * FloatLeg	/	Duration;		// in bps.
	else
		Spread	=	0.0;

	if (ResetSched) delete ResetSched;
	//this->GetRcvOrPay()
	NPV*= ((double)CorridorLeg->GetRcvOrPay());
	setValue(qCMPFEELEGPV,NPV);
	setValue(qCMPSPREAD,Spread);
	setValue(qCMPDEFLEGPV,FloatLeg);
	setValue(qCMPDURATION,Duration);
	SetAccrued(accrued);
	return getValue(measure);
}


void ICM_Pricer_Corridor::View(char* id, FILE* ficOut) {

	FILE* fOut;
	char  fOutName[200];

	if ( ficOut == NULL )
	{
	ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
	
	   (void) unlink(fOutName);

       fOut = fopen(fOutName, "w"); 
    }
	else
	{
	fOut = ficOut;
	} 
	fprintf(fOut, "\n ZC Curve \n");
	itsZeroCurve->View(id, fOut);
	fprintf(fOut, "\n\n Def Index Curve \n");
	itsDefCurveIndex->View(id, fOut);
	fprintf(fOut, "\n\n IR Vol \n");
	itsIRVolatility->View(id, fOut);
	fprintf(fOut, "\n\n Index Vol \n");
	itsIdxVolatility->View(id, fOut);
	fprintf(fOut, "\n\n Correlation Curve \n");
	itsCorrelation->View(id, fOut);
	fprintf(fOut, "\n\n Credit FixingCurve Curve \n");
	itsFixingCurveCredit->View(id, fOut);
	fprintf(fOut, "\n\n IR FixingCurve Curve \n");
	itsFixingCurveIR->View(id, fOut);

	ICM_Pricer_MktMng::View(id, fOut);

	if ( ficOut == NULL )
	{
		fclose(fOut);
	}

}
