/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file cppi_options.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
/*
#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include <cmath>
#include "gpclosedforms/cppi_options.h"
#include "gpclosedforms/vanilla_bs.h"
#include "zerocurv.h"
#include "gpbase/datemanip.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelparamtype.h"
#include "gpinfra/surfacemodelparam.h"
#include "gpinfra/gensecurity.h"
#include "gpmodels/sabr_eq.h"
#include "gpcalib/vanillamepi.h"

*/


#include <glob/firsttoinc.h>
/*
#include "ARM_local_gp_calculators.h"
#include "ARM_local_glob.h"
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_local_wrapper.h"
#include "ARM_local_class.h"
*/

/// gpbase

#include <gpbase/datestripcombiner.h>
#include <gpbase/curve.h>
#include <gpbase/curvetypedef.h>
#include <gpbase/interpolator.h>
#include <gpbase/stringmanip.h>
#include <gpbase/datemanip.h>

#include <gpbase/typedef.h>

#include "gpbase/env.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/autocleaner.h"
#include "gpbase/singleton.h"
#include "gpbase/datestrip.h"
#include "gpbase/stringconvert.h"
#include "gpbase/surface.h"

/// gpinfra
#include <gpinfra/mktdatamanagerrep.h>
#include <gpinfra/gensecurity.h>
#include <gpinfra/pricingmodel.h>
#include <gpinfra/dealdescription.h>
#include <gpinfra/pricingadviser.h>
#include <gpinfra/modelnrefcall.h>
#include <gpinfra/gramnode.h>
#include <gpinfra/modelparam.h>
#include <gpinfra/modelparams.h>
#include "gpinfra/modelnamemap.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/curvemodelparam.h"
#include <gpinfra/surfacemodelparam.h>

#include "gpinfra/genpricer.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/numerairefactory.h"
#include "gpinfra/pricerinfo.h"
#include "gpinfra/gramfunctorargdict.h"
#include "gpinfra/cstmanager.h"



/// gpcalib
#include <gpcalib/calibmethod.h>
#include <gpcalib/vanillamepi.h>

/// gpcalculators
#include <gpcalculators/argconvdefault.h>
#include <gpcalculators/gencalculator.h>
#include <gpcalculators/crfcalculator.h>
#include <gpcalculators/tarncalculator.h>
#include <gpcalculators/tarncalculatorsnowball.h>
#include <gpcalculators/maturitycapcalculator.h>
#include <gpcalculators/captioncalculator.h>
#include <gpcalculators/prdccalculator.h>
#include <gpcalculators/callablesnowballcalculator.h>
#include <gpcalculators/csocalculator.h>
#include <gpcalculators/bermudaswaptioncalculator.h>


/// gpmodels
#include <gpmodels/HW1F.h>
#include <gpmodels/HW2F.h>
#include <gpmodels/SFRM.h>
#include <gpmodels/HybridBasisFwdIR.h>
#include "gpmodels/MultiAssets.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/SABR_Eq.h"

/// gpnumlib
#include "gpnumlib/randomgenfactory.h"
#include "gpnumlib/compositegen.h"
#include "gpnumlib/antitheticgen.h"

/// gpnummethods
#include "gpnummethods/mcmethod.h"
#include "gpnummethods/finummethod.h"

/// gphelp
#include <gphelp/crmcookies.h>

/// ARM Kernel
#include "portfolio.h"
#include "powrev.h"
#include "model.h"
#include <inst/forex.h>
#include <inst/swaption.h>


#include "gpmodels/sabr_eq.h"
#include "gpcalib/vanillamepi.h"
#include "gpclosedforms/cppi_options.h"
#include "gpclosedforms/vanilla_bs.h"

//#include <ARM\libarm_frometk\arm_local_parsexml_util.h>

//// using the namespace directive to access ARM object!
using ARM::ARM_DateStrip;
using ARM::ARM_DateStripPtr;
using ARM::ARM_DateStripVector;
using ARM::ARM_DateStripCombiner;
using ARM::ARM_GenCalculator;
using ARM::ARM_CRFCalculator;
using ARM::ARM_TARNCalculator;
using ARM::ARM_TARNCalculatorSnowBall;
using ARM::ARM_MaturityCapCalculator;
using ARM::ARM_CaptionCalculator;
using ARM::ARM_PRDCCalculator;
using ARM::ARM_CallableSnowBallCalculator;
using ARM::ARM_CSOCalculator;
using ARM::ARM_BermudaSwaptionCalculator;
using ARM::ARM_MarketData_ManagerRep;
using ARM::ARM_GenSecurity;
using ARM::ARM_PricingModel;
using ARM::ARM_CalibMethod;
using ARM::ARM_GenSecurityPtr;
using ARM::ARM_PricingModelPtr;
using ARM::ARM_CalibMethodPtr;
using ARM::ARM_MarketData_ManagerRepPtr;
using ARM::ARM_CRMCookies;
using ARM::ARM_Curve;
using ARM::ARM_HullWhite1F;
using ARM::ARM_HullWhite2F;
using ARM::ARM_SFRM;
using ARM::ARM_HybridBasisFwdIR;
using ARM::ARM_GP_Vector;
using ARM::ARM_GP_Matrix;
using ARM::stringGetUpper;
using ARM::ARM_ObjectVector;
using ARM::ARM_BoolVector;
using ARM::ARM_ModelParam;
using ARM::ARM_VanillaMepi;
using ARM::ConvertXLDateToJulian;






CC_BEGIN_NAMESPACE(ARM)

// MC default pricing values
const double MC_NB_INTER_STEPS		= 1;


////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///												BinomialTree
///
///
///////////////////////////////////////////////////////////////////////////////////////////////////

BinomialTree::BinomialTree(int n, int size,double T0)
	:itsSize(n),path_state_size(size), itsUnderlyingContainer(NULL), 
		itsPathStateContainer(NULL), itsPathOptionContainer(NULL),T(T0)
{
	itsUnderlyingContainer=	new vector<double>(itsSize*itsSize);
	itsPathStateContainer=	new vector<double>(path_state_size*itsSize*itsSize);
	itsPathOptionContainer=	new vector<double>(path_state_size*itsSize*itsSize);
}

BinomialTree::~BinomialTree()
	{
		delete itsUnderlyingContainer;
		delete itsPathStateContainer;
		delete itsPathOptionContainer;
		itsUnderlyingContainer = NULL;
		itsPathStateContainer = NULL;
		itsPathOptionContainer = NULL;
	}


double BinomialTree::Price()
	{
		/// Backward pass
		int i,j,k;
		double Si,Siplus1_Up,Siplus1_Down,Di,Diplus1,PathVariable0,PathVariableUp,PathVariableDown;
		double DownProbability=1.-UpProbability;
		Diplus1=1.0;
		Di=1.0/CapitalizationFactor;
		double current_time;
		for(i=itsSize-2;i>=0;i--)
		{
			current_time=(i+1)*T/(itsSize-1.0);
			Diplus1=Di;
			Di=Di/CapitalizationFactor;
			for(j=0;j<=i;j++)
			{
				Si=get_underlyingstate(i,j);
				Siplus1_Up=get_underlyingstate(i+1,j+1);
				Siplus1_Down=get_underlyingstate(i+1,j);
				for(k=0;k<path_state_size;k++)
				{
					PathVariable0=		get_pathstate(i,j,k);
					PathVariableUp=		PortfolioValueNextStep(current_time,PathVariable0, Si, Siplus1_Up,Di,Diplus1);
					PathVariableDown=	PortfolioValueNextStep(current_time,PathVariable0, Si, Siplus1_Down,Di,Diplus1);
					set_pathoption(i,j,k,(UpProbability*interpole_pathoption(i+1,j+1,PathVariableUp)+DownProbability*interpole_pathoption(i+1,j,PathVariableDown))/CapitalizationFactor);
				}
			}
		}
		return get_pathoption(0,0,0);
	}


double BinomialTree::interpole_pathoption(int i,int j,double v)
	 {
		int k=0;
		if(v<=get_pathstate(i,j,0)) return get_pathoption(i,j,0);
		while((k<path_state_size)&&(v>=get_pathstate(i,j,k))) k++;
		if(k>=path_state_size) return get_pathoption(i,j,path_state_size-1);
		double dp=get_pathstate(i,j,k)-get_pathstate(i,j,k-1);
		if(dp!=0.)
		{
			return get_pathoption(i,j,k-1)+(get_pathoption(i,j,k)-get_pathoption(i,j,k-1))*(v-get_pathstate(i,j,k-1))/dp;
		} 
		return	(get_pathoption(i,j,k)+get_pathoption(i,j,k-1))/2.;
	 }

////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///												BinomialTreeMepiVanillaOption
///
///
///////////////////////////////////////////////////////////////////////////////////////////////////

BinomialTreeMepiVanillaOption::BinomialTreeMepiVanillaOption(int n0, int nz0,double p0,double K0,double Vol0,double r0,double BorrowingSpread0,double YearlyFees0,double T0,
			double minExp,double MaximumExposure0,double riskFactor,double g0,double g,double minport,double maxport,int AveragingPeriodNb0,int CallOrPut0):
		BinomialTree(n0,nz0,T0),K(K0),r(r0),BorrowingSpread(BorrowingSpread0),MinimumExposure(minExp),MinimumGaranty(g),YearlyFees(YearlyFees0),
			PATHSTATE_LIMITDOWN(minport),PATHSTATE_LIMITUP(maxport),CallOrPut(CallOrPut0),RiskFactor(riskFactor),AveragingPeriodNb(AveragingPeriodNb0),InitialMinimumGaranty(g0),
		MaximumExposure(MaximumExposure0)
		{
			S=100.;
			FeesPerPeriod=YearlyFees*T/(itsSize-1.0);
			VolRatio=sqrt(1.0-2.0/3.0*(AveragingPeriodNb-1.0)/(itsSize-1.0));
			VanillaOptionInit(Vol0*VolRatio,p0);
		}
		
		
		
		void BinomialTreeMepiVanillaOption::VanillaOptionInit(double Vol0, double initialPortfolio)
		{
			InitialPortfolio=initialPortfolio;
			Vol=Vol0;
			int i,j,k;
			u=exp((r-0.5*Vol*Vol*VolRatio*VolRatio)*T/(itsSize-1)+Vol*VolRatio*sqrt(T/(itsSize-1)));
			d=exp((r-0.5*Vol*Vol*VolRatio*VolRatio)*T/(itsSize-1)-Vol*VolRatio*sqrt(T/(itsSize-1)));
			CapitalizationFactor=exp(r*T/(itsSize-1));
			BorrowingSpreadingFactor=exp(BorrowingSpread*T/(itsSize-1));
			UpProbability=0.5;
			double St=S*pow(d,itsSize-1);
			double uOverd=u/d;
			double Invd=1/d;
			
			/// Filling of underlying prices
			for(j=0;j<itsSize;j++)
			{
				set_underlyingstate(itsSize-1,j,St);
				St*=uOverd;
			}
			for(i=itsSize-2;i>=0;i--)
			{
				for(j=0;j<=i;j++)
				{
					set_underlyingstate(i,j,Invd*get_underlyingstate(i+1,j));
				}
			}
			
			/// Filling of the portfolio prices and the nb of portfolio prices : Forward pass
			/// the purpose of the forward pass is to compute the limits of the path state for every i,j 
			/// then to fill the intermediate value in an arithmetic way (not geometric to handle the case where we have negative values for the portfolio)
			/// 
			double Pup,Pdown,x0,dx,Di,Diplus1;
			Di=pow(CapitalizationFactor,-(itsSize-1));
			Diplus1=Di*CapitalizationFactor;
			double current_time;
			for(k=path_state_size-1;k>=0;k--)
			{
				unsafe_set_pathstate(0,0,k,InitialPortfolio);
			}	
			for(i=0;i<itsSize-1;i++)
			{
				current_time=(i+1)*T/(itsSize-1.0);
				Di=Diplus1;
				Diplus1=Di*CapitalizationFactor;
				for(j=0;j<=i;j++)
				{
					x0=get_pathstate(i,j,0);
					if(path_state_size>1) dx=(get_pathstate(i,j,path_state_size-1)-x0)/(path_state_size-1.0);
					/*				for(k=1;k<path_state_size-1;k++)
					{
					set_pathstate(i,j,k,x0);
					}
					*/
					for(k=path_state_size-1;k>=1;k--)
					{
						set_pathstate(i,j,k,x0+k*dx);
					}
					for(k=0;k<path_state_size;k++)
					{
						Pup=min(PATHSTATE_LIMITUP*InitialPortfolio,PortfolioValueNextStep(current_time,get_pathstate(i,j,k),get_underlyingstate(i,j),get_underlyingstate(i+1,j+1),Di,Diplus1));
						Pdown=max(PATHSTATE_LIMITDOWN*InitialPortfolio,PortfolioValueNextStep(current_time,get_pathstate(i,j,k),get_underlyingstate(i,j),get_underlyingstate(i+1,j),Di,Diplus1));
						add_pathstate(i+1,j+1,Pup);
						add_pathstate(i+1,j,Pdown);
					}
				}
			}
			for(j=0;j<=itsSize-1;j++)
			{
				x0=get_pathstate(itsSize-1,j,0);
				if(path_state_size>1) dx=(get_pathstate(itsSize-1,j,path_state_size-1)-x0)/(path_state_size-1.0);
				/*			for(k=1;k<path_state_size-1;k++)
				{
				set_pathstate(itsSize-1,j,k,x0);
				}
				*/			for(k=path_state_size-1;k>=1;k--)
				{
					set_pathstate(itsSize-1,j,k,x0+k*dx);
				}
			}
			/// Filling of the last date column of pathoption prices : 
			if(CallOrPut==K_CALL) 
			{  
				for(j=0;j<itsSize;j++) 
				{
					for(k=0;k<path_state_size;k++)
					{
						set_pathoption(itsSize-1,j,k,max(0.0,get_pathstate(itsSize-1,j,k)-K));
					}
				}
			}
			else
			{
				for(j=1;j<itsSize;j++) 
				{
					for(k=0;k<path_state_size;k++)
					{
						set_pathoption(itsSize-1,j,k, max(0.0,K-get_pathstate(itsSize-1,j,k)));
					}
				}
			}
		}
		

	
	void BinomialTreeMepiVanillaOption::add_pathstate(int i,int j,double v){
		if(v<get_pathstate(i,j,0))					set_pathstate(i,j,0,v);
		if(v>get_pathstate(i,j,path_state_size-1))	set_pathstate(i,j,path_state_size-1,v);
	}
	
	double BinomialTreeMepiVanillaOption::PortfolioValueNextStep(double currenttime,double PortfolioPrecedingValue, double S0, double S1,double D0,double D1)
	{
		if(InitialMinimumGaranty==0)
		{
			double exposed_value=max(MinimumExposure*PortfolioPrecedingValue,RiskFactor*(PortfolioPrecedingValue-MinimumGaranty*InitialPortfolio*D0));
			if (exposed_value<PortfolioPrecedingValue)
			{
				return exposed_value*S1/S0+D1/D0*(PortfolioPrecedingValue-exposed_value)-FeesPerPeriod;
			}
			else
			{
				return exposed_value*S1/S0+BorrowingSpreadingFactor*D1/D0*(PortfolioPrecedingValue-exposed_value)-FeesPerPeriod;
			}
		}
		else
		{
			double exposed_value=min(InitialPortfolio*MaximumExposure,
				max(MinimumExposure*PortfolioPrecedingValue,RiskFactor*(PortfolioPrecedingValue-
				(InitialMinimumGaranty+(MinimumGaranty-InitialMinimumGaranty)*currenttime/T)*InitialPortfolio)));
			if (exposed_value<PortfolioPrecedingValue)
			{
				return exposed_value*S1/S0+D1/D0*(PortfolioPrecedingValue-exposed_value)-FeesPerPeriod;
			}
			else
			{
				return exposed_value*S1/S0+BorrowingSpreadingFactor*D1/D0*(PortfolioPrecedingValue-exposed_value)-FeesPerPeriod;
			}
		}
	}



double BinomialTreeMepiVanillaOption::InitialExposition()
{
	if(InitialMinimumGaranty==0)
	{
		return max(MinimumExposure,RiskFactor*(1.0-MinimumGaranty*exp(-r*T)));
		
	}
	else
	{
		return min(MaximumExposure,max(MinimumExposure,RiskFactor*(1.0-InitialMinimumGaranty)));
		
	}
	
}


////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///												BinomialTreeMepiVanillaOptionSV
///
///
///////////////////////////////////////////////////////////////////////////////////////////////////

double BinomialTreeMepiVanillaOptionSV::PriceSV()
	{
		if(LegendreNb>1)
		{
			double mub=log(Vol*Vol*VolRatio*VolRatio*T)+(VolDrift-VolOfVol*VolOfVol/2.0)*T;
			double sigmab=ARM_NumericConstants::ARM_TWO_DIVIDED_BY_SQRT_3*VolOfVol*VolRatio*sqrt(T);
			ReducedGaussHermite_Coefficients c(LegendreNb);
			double Sum=0;
			double x;
			for(int i=0;i<LegendreNb;i++){
				x=exp(sigmab*c.get_point(i)+mub);
				VanillaOptionInit(sqrt(x/T),InitialPortfolio);
				/// change the volatility according to the integrated volatility
				Sum+= Price()*c.get_weight(i);
			}
			return Sum;
		}
		else
		{	
			return Price();
		}
	}


double BinomialTreeMepiVanillaOptionSV::PriceSV_delta()
	{
		double price,shifted_price,x;
			double InitialPortfolio_Save=InitialPortfolio;
			double InitialPortfolio_Shifted=InitialPortfolio*(1.0+0.01);
		if(LegendreNb>1)
		{
			double mub=log(Vol*Vol*VolRatio*VolRatio*T)+(VolDrift-VolOfVol*VolOfVol/2.0)*T;
			double sigmab=ARM_NumericConstants::ARM_TWO_DIVIDED_BY_SQRT_3*Vol*VolRatio*sqrt(T);
			GaussLegendre_Coefficients c(LegendreNb,1e-10,exp(mub+14.0*sigmab*sigmab));
			double Sum=0;
		
			for(int i=0;i<LegendreNb;i++){
				x=c.get_point(i);
				VanillaOptionInit(sqrt(x/T),InitialPortfolio_Save);
				price=Price();
				VanillaOptionInit(sqrt(x/T),InitialPortfolio_Shifted);
				shifted_price=Price();
				/// change the volatility according to the integrated volatility
				Sum+= (shifted_price-price)*IntegratedVolSquarredDensity(x,mub,sigmab)*c.get_weight(i);
			}
			InitialPortfolio=InitialPortfolio_Save;
			return Sum*100.;
		}
		else
		{	
			x=Vol*Vol*VolRatio*VolRatio*T;
			VanillaOptionInit(sqrt(x/T),InitialPortfolio_Save);
			price=Price();
			VanillaOptionInit(sqrt(x/T),InitialPortfolio_Shifted);
			shifted_price=Price();
			InitialPortfolio=InitialPortfolio_Save;
			return (shifted_price-price)*100.;
		}
	}


////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///												Global Functions 
///
///
///////////////////////////////////////////////////////////////////////////////////////////////////


double mepi_VanillaOption_STOBS(double date,double p0,double K, double T, ARM_ZeroCurve* ZcCurve,string CurveName,
						   double br,double YearlyFees, double cashSpread,ARM_SABR_Eq* SABR_Model,string EquityName,
						  double minExp,double maxExp,double riskFac,double g0,double g,
						  double minport,double maxport,int AveragingPeriodNb,double alreadyAsianed,
						  int CallOrPut,int n,int nz,int nbLegendre)
{	
	double r;
	
	if(T==0) 
	{
		r=0;
	}
	else
	{
		r= -log(ZcCurve->DiscountPrice(T))/T;
	}
		
	double vol			= SABR_Model->GetModelParams()->GetModelParam( ARM_ModelParamType::Alpha ).GetValue(0);
	double VolDrift		=0;
	double VolOFVol		= SABR_Model->GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).GetValue(0);
		
	BinomialTreeMepiVanillaOptionSV b(n,nz,nbLegendre,p0, K,vol,VolDrift,VolOFVol,r,br,YearlyFees,T,minExp,maxExp,riskFac,g0,g,minport,maxport,AveragingPeriodNb,CallOrPut) ;
	return b.PriceSV();	
}

double mepi_VanillaOption_STOBS_delta(double date,double p0,double K, double T, ARM_ZeroCurve* ZcCurve,string CurveName,
						   double br,double YearlyFees, double cashSpread,ARM_SABR_Eq* SABR_Model,string EquityName,
						  double minExp,double maxExp,double riskFac,double g0,double g,
						  double minport,double maxport,int AveragingPeriodNb,double alreadyAsianed,
						  int CallOrPut,int n,int nz,int nbLegendre)
{	
	double r;
	
	if(T==0) 
	{
		r=0;
	}
	else
	{
		r= -log(ZcCurve->DiscountPrice(T))/T;
	}
		
	double vol			= SABR_Model->GetModelParams()->GetModelParam( ARM_ModelParamType::Alpha ).GetValue(0);
	double VolDrift		=0;
	double VolOFVol		= SABR_Model->GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).GetValue(0);
		
	BinomialTreeMepiVanillaOptionSV b(n,nz,nbLegendre,p0, K,vol,VolDrift,VolOFVol,r,br,YearlyFees,T,minExp,maxExp,riskFac,g0,g,minport,maxport,AveragingPeriodNb,CallOrPut) ;
	return b.PriceSV_delta();	
}

double mepi_VanillaOption_SABR(double date,double p0,double K, double Maturity_date, ARM_ZeroCurve* ZcCurve,string CurveName,
						   double br,double YearlyFees, double cashSpread,ARM_PricingModel* Model,string EquityName,
						  double minExp,double maxExp,double riskFac,double g0,double g,
						  double minport,double maxport,int AveragingPeriodNb,double alreadyAsianed,int resetfreq)
{			
	ARM_GenSecurity * gensec = NULL;
	double startDate=date;
	double endDate=Maturity_date;
	double resetFreq=resetfreq;
	double maxBorrow=maxExp-1.;
	double protectionCurveStart=g0;
	double protectionCurveEnd=g;
	double startingPortfolio=p0;
	double startingCash=0.0;
	double minInvested=minExp;
	double leverageCost=br;
	double asianDatesNb=AveragingPeriodNb;

	/// fix fix  :only pour verifier les parametres
	/*
		ARM_CountedPtr<ARM_ZeroCurve> zerocurvePtr;
		ARM_PricingModelPtr EquityModelPtr;
		ARM_PricingModelPtr  DiscountingModelPtr;
	ARM_MultiAssetsModel* multiassetPtr=dynamic_cast<ARM_MultiAssetsModel*>(Model);
	if (multiassetPtr)
	{
		ARM_ModelNameMap* modelMapPtr=multiassetPtr->GetModelMap();
		
		ARM_ModelNameMap::iterator iterEquityModel =(*modelMapPtr)[EquityName];
		EquityModelPtr=(*iterEquityModel).Model();
		
		ARM_ModelNameMap::iterator iterDiscountingModel =(*modelMapPtr)[CurveName];
		DiscountingModelPtr=(*iterDiscountingModel).Model();
		zerocurvePtr=DiscountingModelPtr->GetZeroCurve();
		
		
	}
	else
	{
		zerocurvePtr=Model->GetZeroCurve();
		ARM_PricingModelPtr EquityModelPtr1(Model);
		EquityModelPtr=EquityModelPtr1;
	}
	ARM_SABR_Eq* SABR_Model= dynamic_cast<ARM_SABR_Eq*>(&(*EquityModelPtr));

	size_t modelNb		= SABR_Model->GetModelNb();
	double rho			= SABR_Model->GetModelParams()->GetModelParam( ARM_ModelParamType::Correlation ).GetValue(0);
	double beta			= SABR_Model->GetModelParams()->GetModelParam( ARM_ModelParamType::Beta ).GetValue(0);
	double volOfVol		= SABR_Model->GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).GetValue(0);
	double alpha		= SABR_Model->GetModelParams()->GetModelParam( ARM_ModelParamType::Alpha).GetValue(0);
	double f= (dynamic_cast<ARM_ModelParamsSABR_EqFx*>(SABR_Model->GetModelParams()))->GetSpot();
	ARM_MCMethod* nmc=dynamic_cast<ARM_MCMethod*>(&(*(SABR_Model->GetNumMethod())));
	ARM_VectorPtr nbMCPtr=nmc->GetBuckets();
	*/

	/// end of fix fix

	
	ARM_VanillaMepi vm(CurveName, EquityName,
		startDate, endDate,
		resetFreq, riskFac, K, maxBorrow, protectionCurveStart, protectionCurveEnd, startingPortfolio,
		startingCash, minInvested, leverageCost, cashSpread, YearlyFees, alreadyAsianed, asianDatesNb );
	
	vm.CreateAndSetGenSec(CurveName, startDate);
	gensec = static_cast<ARM_GenSecurity *> ( vm.GetGenSecurity()->Clone() );

	ARM_GenPricer* genPricer = new ARM_GenPricer( gensec,Model);	
	double price = genPricer->Price();
	delete gensec;
	delete genPricer;
	return price;
	
}




double mepi_VanillaOption_SABR_delta(double date,double p0,double K, double TDate, ARM_ZeroCurve* ZcCurve,string CurveName,
						   double br,double YearlyFees, double cashSpread,ARM_PricingModel* SABR_Model,string EquityName,
						  double minExp,double maxExp,double riskFac,double g0,double g,
						  double minport,double maxport,int AveragingPeriodNb,double alreadyAsianed,double resetfreq,double shift)
{		
	ARM_GenSecurity * gensec = NULL;
	ARM_GenSecurity * gensec1 = NULL;
	double startDate=date;
	double endDate=TDate;
	double resetFreq=resetfreq;
	double maxBorrow=maxExp;
	double protectionCurveStart=g0;
	double protectionCurveEnd=g;
	double startingCash=0.0;
	double minInvested=minExp;
	double leverageCost=br;
	double asianDatesNb=AveragingPeriodNb;

/*	double portfolioshift=p0*shift;	
	startingPortfolio=p0-portfolioshift/2.0;
	ARM_VanillaMepi vm(CurveName, EquityName,
		startDate, endDate,
		resetFreq, riskFac, K, maxBorrow, protectionCurveStart, protectionCurveEnd, startingPortfolio,
		startingCash, minInvested, leverageCost, cashSpread, YearlyFees, alreadyAsianed, asianDatesNb );
	
	vm.CreateAndSetGenSec(CurveName, startDate);
	gensec = static_cast<ARM_GenSecurity *> ( vm.GetGenSecurity()->Clone() );
	ARM_GenPricer* genPricer = new ARM_GenPricer( gensec,SABR_Model);	
	double price0 = genPricer->Price();


	startingPortfolio=p0+portfolioshift/2.0;
	ARM_VanillaMepi vm1(CurveName, EquityName,
		startDate, endDate,
		resetFreq, riskFac, K, maxBorrow, protectionCurveStart, protectionCurveEnd, startingPortfolio,
		startingCash, minInvested, leverageCost, cashSpread, YearlyFees, alreadyAsianed, asianDatesNb );

	vm1.CreateAndSetGenSec(CurveName, startDate);
	gensec1 = static_cast<ARM_GenSecurity *> ( vm1.GetGenSecurity()->Clone() );
	ARM_GenPricer* genPricer1 = new ARM_GenPricer( gensec1,SABR_Model);	
	double price1 = genPricer1->Price();
	delete gensec;
	delete gensec1;
	delete genPricer;
	delete genPricer1;
	return(price1-price0)/portfolioshift;
	*/

	ARM_VanillaMepiDelta vm(CurveName, EquityName,
		startDate, endDate,
		resetFreq, riskFac, K, maxBorrow, protectionCurveStart, protectionCurveEnd, p0,
		startingCash, minInvested, leverageCost, cashSpread, YearlyFees, alreadyAsianed, asianDatesNb, shift);
	
	vm.CreateAndSetGenSec(CurveName, startDate);
	gensec = static_cast<ARM_GenSecurity *> ( vm.GetGenSecurity()->Clone() );

	ARM_GenPricer* genPricer = new ARM_GenPricer( gensec,SABR_Model);	
	double price = genPricer->Price();
	delete gensec;
	delete genPricer;
	return price;
	
}
double mepi_VanillaOption_delta(double p0,double K, double T, double r, double br,double YearlyFees,double vol,double VolDrift,double VolOFVol,
								double minExp,double maxExp,double riskFac,double g0,double g,double minport,double maxport,int AveragingPeriodNb,int CallOrPut,int n,int nz,int nbLegendre)
{
	BinomialTreeMepiVanillaOptionSV b(n,nz,nbLegendre,p0, K,vol,VolDrift,VolOFVol,r,br,YearlyFees,T,minExp,maxExp,riskFac,g0,g,minport,maxport,AveragingPeriodNb,CallOrPut) ;
	return b.PriceSV_delta();

}

double IntegratedVolSquarredDensity(double x,double mub,double sigmab)
	{
		/// the process is dsigma/sigma=VolDrift*dt + VolOfVol*dw2(t)
		/// so we do the normal approximation to compute the variance, then we use the normal / lognormal equivalence to have a true density 
			
			double arg=(mub-log(x))/sigmab;	
			return exp(-arg*arg/2.0)/(ARM_NumericConstants::ARM_SQRT_2_PI*sigmab*x);
	}



double Export_Fund_VanillaOption(double f,double K,double T,double r,double sig,double VolDrift,double VolOfVol,double callput,int LegendreNb)
{
	if(LegendreNb>1)
		{
			
			double mub=log(sig*sig*T)+(VolDrift-VolOfVol*VolOfVol/2.0)*T;
			double sigmab=ARM_NumericConstants::ARM_TWO_DIVIDED_BY_SQRT_3*VolOfVol*sqrt(T);
			ReducedGaussHermite_Coefficients c(LegendreNb);
			double Sum=0;
			double x,Vol;
			int i;
			for(i=0;i<LegendreNb;i++){
				x=exp(sigmab*c.get_point(i)+mub);
				Vol=sqrt(x/T);
				/// change the volatility according to the integrated volatility
				Sum+= BlackSholes_Formula(f,Vol*sqrt(T),exp(-r*T),K,callput)*c.get_weight(i);
			}
			return Sum;
		}
		else
		{	
			return  BlackSholes_Formula(f,sig*sqrt(T),exp(-r*T),K,callput);
		}
}



CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
