/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file MarketHybridModel.cpp
 *
 *  \brief 
 *
 *	\author  J-M Prié
 *	\version 1.0
 *	\date March 2006
 */


/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/MarketHybridModel.h"

/// gpmodel
#include "gpmodels/MarketIRModel.h"
#include "gpmodels/BS_ModelParams.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/ModelParamsQ1f.h"
#include "gpmodels/Q1f.h"

/// gpbase
#include "gpbase/gpmatrix.h"

/// gpinfra
#include "gpinfra/mktdatamanagerrep.h"
#include "gpinfra/pricingstates.h"

/// gpclosedform
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/spreadoption_lognormal_formula.h"
#include "gpclosedforms/spreadoption_lognormal_interface.h"
#include "gpclosedforms/vanilla_normal.h"
#include "gpclosedforms/spreadoption_normal_formula.h"
#include "gpclosedforms/spreadoption_normal_interface.h"

/// kernel
#include "mod/bsmodel.h"
//#include "crv/volcurv.h"

CC_BEGIN_NAMESPACE( ARM )

/// Default MDM key names
const string MKT_IR_KEY_NAME		= "MKT_IR_MODEL";
const string MKT_EQFX_KEY_NAME		= "MKT_EQFX_MODEL";
const string MKT_HYB_DATAS_KEY_NAME	= "MKT_HYB_DATAS";

////////////////////////////////////////////////////
///	Class  : ARM_MarketHybridModel
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_MarketHybridModel::ARM_MarketHybridModel(const ARM_Date& asOfDate,const ARM_ObjectVector& mktDatas,const ARM_StringVector& mktKeys)
:	ARM_AnalyticIRModel(), itsIsMarketVNS(true), itsIsLogNorRates(true), itsIrStrikeStatus(0),
	itsEqFxRate(0.0),itsEqFxStrike(0.0),itsIrRate(0.0),itsIrStrike(0.0),
	itsRefEqFxOptionIdx(-1)
{
	/// Analyse input markets & models to keep predefined order and
	/// fill market data manager
	ARM_StringVector mdmKeys(0);
	size_t i,nbDatas=mktDatas.size();
	ARM_IntVector idx(nbDatas);
	for(i=0;i<nbDatas;++i)
		idx[i]=i;

	ARM_MarketIRModel* irModel;

	/// Create an empty MdM in the model
	ARM_MarketData_ManagerRep mdm(asOfDate);
	SetMktDataManager(mdm,mdmKeys);

	/// Find Market IR Model. If not found VNS evaluation will rely on
	/// precomputed hybrid datas
	for(i=0;i<idx.size();++i)
	{
		if((irModel=dynamic_cast<ARM_MarketIRModel*>(mktDatas[idx[i]])))
		{

			///... then fill it
			mdmKeys.push_back(idx[i]<mktKeys.size() ? mktKeys[idx[i]] : MKT_IR_KEY_NAME);
			GetMktDataManager()->RegisterData(mdmKeys[mdmKeys.size()-1],irModel);
			idx.erase(idx.begin()+i);
			break;
		}
	}
	itsIsMarketVNS = (i<idx.size());
	if(!itsIsMarketVNS)
		mdmKeys.push_back(MKT_IR_KEY_NAME); /// keep keys order


	/// Find Market EQ/FX Model
	for(i=0;i<idx.size();++i)
	{
		if(dynamic_cast<ARM_BSModel*>(mktDatas[idx[i]]))
		{
			mdmKeys.push_back(idx[i]<mktKeys.size() ? mktKeys[idx[i]] : MKT_EQFX_KEY_NAME);
			GetMktDataManager()->RegisterData(mdmKeys[mdmKeys.size()-1],mktDatas[idx[i]]);
			idx.erase(idx.begin()+i);
			break;
		}
	}
	if(i>=idx.size())
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : market EQ/FX model is missing !");


	/// Find hybrid datas
	for(i=0;i<idx.size();++i)
	{
		if(dynamic_cast< ARM_GP_T_Vector< ARM_VectorPtr > *>(mktDatas[idx[i]]))
		{
			mdmKeys.push_back(idx[i]<mktKeys.size() ? mktKeys[idx[i]] : MKT_HYB_DATAS_KEY_NAME);
			GetMktDataManager()->RegisterData(mdmKeys[mdmKeys.size()-1],mktDatas[idx[i]]);

			/// Update MdM keys
			SetKeys(mdmKeys);

			return;
		}
	}

	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : market hybrid datas are missing !");
}


////////////////////////////////////////////////////
///	Class   : ARM_MarketHybridModel
///	Routines: Copy Constructor
///	Returns :
///	Action  :
////////////////////////////////////////////////////
ARM_MarketHybridModel::ARM_MarketHybridModel(const ARM_MarketHybridModel& rhs )
: ARM_AnalyticIRModel(), itsIsMarketVNS(rhs.itsIsMarketVNS), itsIsLogNorRates(rhs.itsIsLogNorRates),
  itsIrStrikeStatus(rhs.itsIrStrikeStatus),itsRefEqFxOptionIdx(rhs.itsRefEqFxOptionIdx),
  itsEqFxRate(rhs.itsEqFxRate),itsEqFxStrike(rhs.itsEqFxStrike),
  itsIrRate(rhs.itsIrRate),itsIrStrike(rhs.itsIrStrike),
  itsDomZcModel( CreateClonedPtr(&*(rhs.itsDomZcModel)) ),
  itsForZcModel( CreateClonedPtr(&*(rhs.itsForZcModel)) )
{}


////////////////////////////////////////////////////
///	Class   : ARM_MarketHybridModel
///	Routines: SetHybridDatas
///	Returns : nothing
///	Action  : Register hybrid datas (correlations & vols) in
///			  its market data manager
////////////////////////////////////////////////////
void ARM_MarketHybridModel::SetHybridDatas(ARM_GP_T_Vector< ARM_VectorPtr >& hybridDatas)
{
	GetMktDataManager()->RegisterData(GetKeys()[HybridDatasKey],&hybridDatas);
}


////////////////////////////////////////////////////
///	Class   : ARM_MarketHybridModel
///	Routines: GetSettlementCalendar
///	Returns : string
///	Action  : default implementation
////////////////////////////////////////////////////
string ARM_MarketHybridModel::GetSettlementCalendar(const string& modelName) const
{
    ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[MarketEqFxModelKey]) );
	return ARM_PricingFunctionEquity::GetSettlementCalendar(fxBSModel->GetYCModel()->GetZeroCurve(),fxBSModel->GetDividend());
}


////////////////////////////////////////////////////
///	Class   : ARM_MarketHybridModel
///	Routines: GetSettlementGap
///	Returns : double
///	Action  : default implementation
////////////////////////////////////////////////////
double ARM_MarketHybridModel::GetSettlementGap(const string& modelName) const
{
    ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[MarketEqFxModelKey]) );
	return ARM_PricingFunctionEquity::GetSettlementGap(fxBSModel->GetYCModel()->GetZeroCurve(),fxBSModel->GetDividend());
}

////////////////////////////////////////////////////
///	Class  : ARM_MarketHybridModel
///	Routine: HybridCallVectorial
///	Returns: price a hybrid call that is a spread
///			 between a forward equity or FX strip
///			 and an IR variable notional swap (VNS) :
///
///			 Call = Et[max(EqFxStrip(T) - (FloatLeg(T)-Strike*FixO1dom(T)) , 0)] 
///				  = FixO1dom(T) * Et[max(VNXRate(T) - (VNSRate(T) - Strike) , 0)] 
///				  = FixO1dom(T) * Et[max(VNXRate(T) - VNSRate(T) - (-Strike) , 0)] 
///				  where VNXRate(T) = EqFxStrip(T)/FixO1dom(T) & VNSRate(T) = FloatLeg(T)/FixO1dom(T)
///
///			 Put  = Et[max((FloatLeg(T)-Strike*FixO1dom(T) - EqFxStrip(T)) , 0)]
///				  = FixO1dom(T) * Et[max((-Strike) -(VNXRate(T) - VNSRate(T)), 0)]
///
///			 If Market VNS is activated, VNS rate volatility is a market data and
///			 correlation with FX rate is read from hybrid datas structure
///
///			 If Market VNS flag is not activated, VNS rate volatility & correlations
///			 are read from hybrid datas structure
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MarketHybridModel::HybridCallVectorial(
	const string& modelName,
    double evalTime,
	double expiryTime,
	int callPut,
	const std::vector<double>& strikesPerState,

	/// Strip of forwards FX (or equity)
	const std::vector<double>& fxExpiryTimes,
	const std::vector<double>& fxSettlementTimes,
	const std::vector<double>& fxPayTimes,
	const std::vector<double>& fxNotionals,

	/// IR Swap
	double swapResetTime,
	const std::vector<double>& fixNotionals,
	const std::vector<double>& floatNotionals,
	double floatStartTime,
	double floatEndTime,
	const std::vector<double>& floatResetTimes,
	const std::vector<double>& floatStartTimes,
	const std::vector<double>& floatEndTimes,
	const std::vector<double>& floatIntTerms,
	const std::vector<double>& fixPayTimes,
	const std::vector<double>& fixPayPeriods,

	const ARM_PricingStatesPtr& states,
    ARM_PricingContext* context) const
{

	/// Only spot evaluation allowed
	if(evalTime>0.0 || states->size()!=1 || strikesPerState.size()!=1)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : hybrid call evaluation only available at spot date");

	size_t i,nbFwds=fxExpiryTimes.size();
	if(fxSettlementTimes.size()!=nbFwds || fxPayTimes.size()!=nbFwds || fxNotionals.size()!=nbFwds )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : inconsistency in fx strip schedule");

	if(expiryTime>swapResetTime)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : call expiry must be prior to swap reset date");

	for(i=0;i<nbFwds;++i)
	{
		if(expiryTime>fxExpiryTimes[i])
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : call expiry must be prior to any FX reset dates");
	}

	/// Restore hybrid datas at right expiry date
    ARM_GP_T_Vector<ARM_VectorPtr> * hybridDatasStruct = dynamic_cast< ARM_GP_T_Vector<ARM_VectorPtr> * >( GetMktDataManager()->GetData(GetKeys()[HybridDatasKey]) );
	ARM_VectorPtr hybridDatas;
	if(!hybridDatasStruct)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : hybrid datas structure is missing");

	size_t nbHybridDatas=hybridDatasStruct->size();
	for(i=0;i<nbHybridDatas;++i)
	{
		hybridDatas = (*hybridDatasStruct)[i];
		if(hybridDatas->size() < Time+1)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : time is missing in the hybrid datas structure");
		if((*hybridDatas)[Time] <= expiryTime+K_NEW_DOUBLE_TOL)
			break;
	}
	if(hybridDatas->size() < NbDatas)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : hybrid datas structure doesn't contain all required datas");

	double domForCor,domFxCor,forFxCor,irFxSpreadCor,domO1Vol,forO1Vol,irVol,fxVol;
	domForCor		=(*hybridDatas)[DomForCor];		// integrated correlation between fwd dom fix O1 (=FixO1dom(t)/Bdom(t,T)) & fwd for O1 (=vol of O1for(t)/Bfor(t,T))
	domFxCor		=(*hybridDatas)[DomFxCor];		// integrated correlation between 1st fwd FX (=X(t,T)) & fwd dom fix O1 
	forFxCor		=(*hybridDatas)[ForFxCor];		// integrated correlation between 1st fwd FX & fwd for O1 
	irFxSpreadCor	=(*hybridDatas)[IrFxSpreadCor];	// integrated correlation between VNX rate (=EqFxStrip(t)/FixO1dom(t)) & VNS rate (=FloatLeg(t)/FixO1dom(t))

	domO1Vol		=(*hybridDatas)[DomO1LogNorVol];// lognormal fwd dom fix O1 volatility
	forO1Vol		=(*hybridDatas)[ForO1LogNorVol];// lognormal fwd for O1 volatility
	irVol			=(*hybridDatas)[IrNorVol];		// normal VNS rate volatility


	/// Only FX hybrid allowed at the moment
    ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[MarketEqFxModelKey]) );
	ARM_ZeroCurve* domZcCurve = fxBSModel->GetYCModel()->GetZeroCurve();
	ARM_ZeroCurve* forZcCurve = fxBSModel->GetDividend();
	if( !strcmp(domZcCurve->GetCurrencyUnit()->GetCcyName(), forZcCurve->GetCurrencyUnit()->GetCcyName()) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : hybrid call evaluation only available for FX/IR hybrids");

	double strike = strikesPerState[0];

	/// Compute FX strip forward value
	double fwdFx,fxStripPrice=0.0;
	double zcPay,zcDomSet,zcForSet,fxSpot = fxBSModel->GetSpot();
	double flow,forNumeraire=0.0;
	for(i=0;i<nbFwds;++i)
	{
		fwdFx = ARM_ModelParams_Fx::Forward(domZcCurve,forZcCurve,fxSpot,fxSettlementTimes[i]);
		zcPay = domZcCurve->DiscountPrice(fxPayTimes[i]/K_YEAR_LEN);
		flow  = zcPay*fxNotionals[i];
		fxStripPrice += fwdFx*flow;

		zcForSet = forZcCurve->DiscountPrice(fxSettlementTimes[i]/K_YEAR_LEN);
		if(fxSettlementTimes[i] != fxPayTimes[i])
		{
			zcDomSet	= domZcCurve->DiscountPrice(fxSettlementTimes[i]/K_YEAR_LEN);
			flow		= zcForSet/zcDomSet*flow;
		}
		else
			flow		= zcForSet*fxNotionals[i];

		forNumeraire	+= flow;
	}

	/// Compute IR swap floating leg forward value
	size_t nbFloatFlows=floatNotionals.size();
	double floatLegPrice=0.0;
	double zcEnd,zcStart = domZcCurve->DiscountPrice(floatStartTimes[0]/K_YEAR_LEN);
	for(i=0;i<nbFloatFlows;++i)
	{
		zcEnd			= domZcCurve->DiscountPrice(floatEndTimes[i]/K_YEAR_LEN);
		floatLegPrice	+= (zcStart-zcEnd)*floatNotionals[i];
		zcStart			= zcEnd;
	}

	/// Compute reference numeraire (= domestic fixed leg O1)
	size_t nbFixFlows=fixNotionals.size();
	double domFixNumeraire=0.0;
	for(i=0;i<nbFixFlows;++i)
	{
		zcPay			= domZcCurve->DiscountPrice(fixPayTimes[i]/K_YEAR_LEN);
		domFixNumeraire	+= zcPay*fixPayPeriods[i]*fixNotionals[i];
	}

	/// Compute Zcs at settlemenet date of the 1st option (previously is was
	/// the settle date of the hybrid option expiry)
	double expiryYf=expiryTime/K_YEAR_LEN;
	double refFxOptionSetTime;

	if(itsRefEqFxOptionIdx<0)
	{
		/// Reference Fx option expiries at the hybrid option expiry
		ARM_Date expirySetDate( ARM_ModelParams_Fx::ComputeSettlementDate(domZcCurve,forZcCurve,domZcCurve->GetAsOfDate()+expiryTime) );
		refFxOptionSetTime = expirySetDate.GetJulian()-domZcCurve->GetAsOfDate().GetJulian();
	}
	else
		/// Reference Fx option is defined by the index itsRefEqFxOptionIdx
		refFxOptionSetTime = fxSettlementTimes[itsRefEqFxOptionIdx < nbFwds ? itsRefEqFxOptionIdx : nbFwds-1];

	zcDomSet = domZcCurve->DiscountPrice(refFxOptionSetTime/K_YEAR_LEN);
	zcForSet = forZcCurve->DiscountPrice(refFxOptionSetTime/K_YEAR_LEN);

	double fxSwapRate	= fxStripPrice/domFixNumeraire;
	double irSwapRate	= floatLegPrice/domFixNumeraire;

	double irSwapStrike = fxSwapRate + strike;
	if(itsIsLogNorRates && irSwapStrike < 0.0005) // 5bp
	{
		/// Save strike limitation status
		if(irSwapStrike < -0.0025) // -25bp
			++itsIrStrikeStatus;

		/// Limit strike to allow a lognormal implied vol
		irSwapStrike = 0.0005;
	}

	int swapCallPut = K_PUT;	// Implicit VNS receiver swaption

	if(itsIsMarketVNS)
	{
		/// Review the normal IR volatility assuming a market VNS (=> option on IR swap basket)

		/// Price variable notional swaption to get
		/// its forward value, normal volatility and numeraire
		ARM_MarketIRModel* mktIRModel = static_cast< ARM_MarketIRModel* >( GetMktDataManager()->GetData(GetKeys()[MarketIrModelKey]) );
		ARM_GP_Matrix strikes(1,nbFixFlows,irSwapStrike);
		bool isConstantNotional=false;
		bool isConstantSpread=true;
		bool isConstantStrike=true;
		ARM_VectorPtr vnsPrice = mktIRModel->VanillaSwaption(modelName,evalTime,
			expiryTime,fixNotionals,floatNotionals,
			floatStartTime,floatEndTime,floatResetTimes,floatStartTimes,floatEndTimes,floatIntTerms,
			fixPayTimes,fixPayPeriods,strikes,swapCallPut,
			states,
			isConstantNotional,isConstantSpread,isConstantStrike);
		irVol = mktIRModel->GetVnsBasketVol();
	}
	if(itsIsLogNorRates)
	{
		/// Compute a normal price then convert to its implied lognormal vol
		double norPrice = VanillaOption_N(irSwapRate,irVol,irSwapStrike,expiryYf,swapCallPut);
		double initNorVol = irVol/irSwapRate;
		irVol = VanillaImpliedVol_BS (irSwapRate,irSwapStrike,expiryYf,norPrice,swapCallPut,&initNorVol);
	}


	/// Compute lognormal volatility of reference FX at its equivalent strike
	/// from 0 to hybrid option expiry assuming a constant local volatility
	/// If reference FX is the option that expiries at hybrid option expiry (i.e. itsRefEqFxOptionIdx=-1)
	/// there is no assumption because this FX volatility is market data
	ARM_FXVolCurve* fxBSVol = dynamic_cast<ARM_FXVolCurve*>(fxBSModel->GetVolatility());
	if(!fxBSVol)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : hybrid call evaluation doesn't support the input FX Vol");
	double fxStrike = (floatLegPrice-strike*domFixNumeraire)*zcForSet/(forNumeraire*zcDomSet);
	fwdFx = ARM_ModelParams_Fx::Forward(domZcCurve,forZcCurve,fxSpot,refFxOptionSetTime);
	double fxLogNorVol = 0.01*fxBSVol->computeVol(fxStrike/fwdFx,expiryYf);


	/// Compute lognormal variable notional fx swaption volatility
	/// VNXRate(t) = X(t,T) * (O1for(t)/Bfor(t,T)) / (FixO1dom(t)/Bdom(t,T))
	/// + lognormal assumption for X(t,T), O1for(t)/Bfor(t,T) and FixO1dom(t)/Bdom(t,T)
	/// Conversion to equivalent normal vol is possible
	double fxSwapLogNorVol = domO1Vol*(domO1Vol - 2*forO1Vol*domForCor) + forO1Vol*forO1Vol
							 + fxLogNorVol*( fxLogNorVol + 2*(forO1Vol*forFxCor - domO1Vol*domFxCor) );
	if(fxSwapLogNorVol<K_NEW_DOUBLE_TOL)
		fxSwapLogNorVol = 0.0;
	else
		fxSwapLogNorVol = sqrt(fxSwapLogNorVol);

	fxVol = fxSwapLogNorVol;
	if(!itsIsLogNorRates)
	{
		/// Compute a lognormal price then convert to its normal vol
		double fxSwapStrike	= irSwapRate - strike;
		double logNorPrice = BS(fxSwapRate,fxSwapStrike,expiryYf,fxVol,K_CALL);
		double initLogNorVol = fxSwapRate*fxVol;
		fxVol = VanillaImpliedVol_N(fxSwapRate,logNorPrice,fxSwapStrike,expiryYf,K_CALL,&initLogNorVol);
	}


	/// Compute IR/FX spread option : the closed form formula is
	/// coded for Call=S1-S2-K and  Put=K-(S1-S2)
	double spreadStrike = - strike;
	double optionType = ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION;
	int n = 64; // Nb pts in Gauss-Legendre
	double hybridPrice;
	if(itsIsLogNorRates)
		hybridPrice = domFixNumeraire * Export_LogNormal_SpreadOption(fxSwapRate,irSwapRate,
			fxVol,irVol,irFxSpreadCor,spreadStrike,expiryYf,callPut,optionType,n);
	else
		hybridPrice = domFixNumeraire * Export_Normal_SpreadOption(fxSwapRate,irSwapRate,
			fxVol,irVol,irFxSpreadCor,spreadStrike,expiryYf,callPut,optionType);

	/// Save computed rates & strikes
	itsEqFxRate		= fxSwapRate;
	itsEqFxStrike	= fxStrike;
	itsIrRate		= irSwapRate;
	itsIrStrike		= irSwapStrike;

	return ARM_VectorPtr(new std::vector<double>(1,hybridPrice));
}

////////////////////////////////////////////////////
///	Class  : ARM_MarketHybridModel
///	Routine: ImpliedVol
///	Returns: int
///	Action : tells the type of the model
////////////////////////////////////////////////////
double ARM_MarketHybridModel::ImpliedVol(const ARM_VanillaArg& arg) const
{
    return DefaultImpliedVol(arg);
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

