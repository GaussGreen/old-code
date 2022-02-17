/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 */

/// this header comes first as it includes some preprocessor constants!

/// gpcalib
#include "gpcalib/basket.h"
#include "gpcalib/basketsensi.h"
#include "gpcalib/stripper.h"
#include "gpcalib/vanillapricer.h"
#include "gpcalib/vanillaspreadoption.h"
#include "gpcalib/vanillaswaption.h"
#include "gpcalib/kerneltogp.h"
#include "gpcalib/argconvdefault.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/datestrip.h"
#include "gpbase/gpmatrixlinalg.h"


/// gpmodels
#include "gpmodels/marketirmodel.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingstates.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_normal.h"

// kernel
//#include <inst/spreadoption.h>
//#include <inst/swap.h>
//#include <inst/swaption.h>
#include <util/refvalue.h>
//#include <inst/fixleg.h>


CC_BEGIN_NAMESPACE( ARM )

const double DEFAULT_PRICE			= 1.0e+100;
const double DEFAULT_PRECISION		= 1.0;
const double DEFAULT_WEIGHT			= 1.0;

const double SENSI_CMS				= 0.000001;

const double NON_CALL_FEE			= 1.0e15;



////////////////////////////////////////////////////
///	Class  : ARM_BasketCalib
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_BasketCalib::ARM_BasketCalib(ARM_DateStrip* ds,const ARM_ReferenceValue& notionalProfile,const ARM_ReferenceValue& feesProfile,double side,BasketType bt, BasketStrike strike)
:	itsExerDateStrip(static_cast<ARM_DateStrip*>(ds->Clone())),
	itsNotional(notionalProfile),
	itsCallFees(feesProfile),
	itsSide(side),
	itsBasketType(bt),
	itsBasketStrike(strike),
	itsPortfolio(ARM_StdPortfolioPtr(NULL)),
	itsBasketCoefs(0)
{

}


////////////////////////////////////////////////////
///	Class  : ARM_BasketCalib
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_BasketCalib::ARM_BasketCalib(const ARM_BasketCalib& rhs)
:	itsExerDateStrip(static_cast<ARM_DateStrip*>(rhs.itsExerDateStrip->Clone())),
	itsNotional(rhs.itsNotional),
	itsCallFees(rhs.itsCallFees),
	itsSide(rhs.itsSide),
	itsBasketType(rhs.itsBasketType),
	itsBasketStrike(rhs.itsBasketStrike),
	itsPortfolio(NULL),//(rhs.itsPortfolio!=ARM_StdPortfolioPtr(NULL) ?  (ARM_StdPortfolio*)rhs.itsPortfolio->Clone() : NULL)),
	itsBasketCoefs(rhs.itsBasketCoefs)
{
}	



////////////////////////////////////////////////////
///	Class  : ARM_BasketCalib
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_BasketCalib::~ARM_BasketCalib()
{
	for (int i=0;i<itsStructLeg.size();i++)
	{
		delete itsStructLeg[i];
	}
	if (itsExerDateStrip)
		delete itsExerDateStrip;
}

////////////////////////////////////////////////////
///	Class  : ARM_BasketCalib
///	Routine: Price
///	Returns: void
///	Action : computes price with mkmo
////////////////////////////////////////////////////
void ARM_BasketCalib::Price(ARM_PricingModel* pMkmo)
{
	//double asOfDate = 0;// pMkmo->GetZeroCurve()->GetAsOfDateJul();
	//ARM_MarketIRModel* pCastModel = dynamic_cast<ARM_MarketIRModel*> (pMkmo);

	//for(int i=0;i<itsPortfolio->GetSize();++i)
	//{
	//	ARM_VanillaSwaptionArg* swaptionGP = static_cast<ARM_VanillaSwaptionArg*>(ARM_ConverterFromKernel::ConvertSecuritytoArgObject(itsPortfolio->GetAsset(i), asOfDate));
	//	double price = swaptionGP->Price(pMkmo);
	//	itsPortfolio->SetPrice (price, i);
	//	std::vector<double> coefs = pCastModel->GetBasketCoefs();
	//	if (i==0)
	//		itsBasketCoefs.resize(coefs.size(),itsPortfolio->GetSize());
	//	for (int j=0;j<coefs.size();j++)
	//		itsBasketCoefs(j,i)=coefs(j);
	//}
}

////////////////////////////////////////////////////
///	Class  : ARM_BasketCalib
///	Routine: Compute
///	Returns: void
///	Action : computes sensis of securities and equivalent VNS
////////////////////////////////////////////////////
void ARM_BasketCalib::Compute(vector<ARM_Security*> securities, vector<ARM_Model*> models, const std::vector<double>& weights)
{
	/*double mktprice;
	size_t size = securities.size();

	std::vector<double>* resetDates = itsExerDateStrip->GetResetDates();
	std::vector<double>* startDates = itsExerDateStrip->GetFlowStartDates();
	
	if (size>0)
	{
		itsStructLeg.resize(size);
		
		for(size_t i=0;i<size;i++)
		{
			ARM_BSModel* BSModel = dynamic_cast< ARM_BSModel* >(models[i]);
			if (BSModel)
			{
				ARM_SpreadOption* spreadOption = dynamic_cast<ARM_SpreadOption*>(securities[i]);
				if (spreadOption)
				{
					spreadOption->SetModel(BSModel);
					mktprice	= spreadOption->ComputePrice();
					if (spreadOption->IsCorridorSpread())
						itsStructLeg[i] = new ARM_BasketSensiCSO(models[i]->GetZeroCurve(),spreadOption,weights[i],mktprice);
					else if (spreadOption->IsSpreadOption())
						itsStructLeg[i] = new ARM_BasketSensiSO(models[i]->GetZeroCurve(),spreadOption,weights[i],mktprice);

				}
				else
				{
					ARM_SwapLeg* leg = dynamic_cast<ARM_SwapLeg*>(securities[i]);
					if (leg)
					{
						leg->SetModel(BSModel);
						mktprice	= leg->ComputePrice();
						if (leg->IsFloatLeg())
							itsStructLeg[i] = new ARM_BasketSensiFunding(models[i]->GetZeroCurve(),leg,weights[i],mktprice);
						else if (leg->IsFixedLeg())
							itsStructLeg[i] = new ARM_BasketSensiFixLeg(models[i]->GetZeroCurve(),leg,weights[i],mktprice);
						else
							ARM_THROW( ERR_INVALID_ARGUMENT, "Basket : swapleg neither float nor fixed" );

					}
					else
						ARM_THROW( ERR_INVALID_ARGUMENT, "Basket : only soleg, floatleg and fixleg supported" );
				}
			}
			else
				ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BasketCalib::Compute : check model type" );
		}
		
		itsEndDate = itsStructLeg[0]->GetEndDate();
		itsRealEndDate = itsStructLeg[0]->GetRealEndDate();
		for(i=0;i<size;i++)
		{
			if (itsStructLeg[i]->GetEndDate()>itsEndDate)
				itsEndDate = itsStructLeg[i]->GetEndDate();
			if (itsStructLeg[i]->GetRealEndDate()>itsRealEndDate)
				itsRealEndDate = itsStructLeg[i]->GetRealEndDate();

			itsStructLeg[i]->FillCpnIndex(*resetDates,*startDates);
		}

		CreateAndSetPortfolio(models[0]->GetZeroCurve());
	}*/
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BasketCalib
///	Routine: CreateSwaptionPortfolio
///	Returns: a portfolio
///	Action : create the list of diagonal swaptions
/////////////////////////////////////////////////////////////////
//void ARM_BasketCalib::CreateAndSetPortfolio(ARM_ZeroCurve* curve)
//{
//	//double asOfDate	  = curve->GetAsOfDateJul();
//
//	//std::vector<double>& resetDates = itsExerDateStrip->GetResetDates();
//	//std::vector<double>& startDates = itsExerDateStrip->GetFlowStartDates();
//
//	//
//	//int i, nbFlows;
//
//	//list< ARM_Security* > swaptionList;
//
//	//nbFlows = startDates->size();
//
//	//// End Date of all the calibrated swaptions
//	//ARM_Date endDate = GetEndDate();
//
//	//double fees;
//	//
//	//for (i = 0; i < nbFlows; ++i)
//	//{
//	//	ARM_Date startDate ((*startDates)[i]);
//
//	//	fees = itsCallFees.Interpolate((*resetDates)[i]);
//
//	//	if (((*resetDates)[i] > asOfDate) && (fees < NON_CALL_FEE))
//	//	{
//	//		ARM_Swaption* swaption = NULL;
//
//	//		swaption = CreateVarNotionalSwaptionAtExer(curve,i,fees);
//	//			
//	//		if (swaption)
//	//			swaptionList.push_back(swaption);
//	//	}
//	//}
//
//	//itsPortfolio = ARM_StdPortfolioPtr(new ARM_StdPortfolio(swaptionList));
//	//for(i=0;i<itsPortfolio->size();++i)
//	//{ 
//	//	itsPortfolio->SetPrecision(DEFAULT_PRECISION,i);
// //       itsPortfolio->SetWeight(DEFAULT_WEIGHT,i);
//	//	itsPortfolio->SetPrice(DEFAULT_PRICE,i);
//	//}
//}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BasketCalib
///	Routine: CreateVarNotionalSwaptionAtExer
///	Returns: void
///	Action : create the variable notional swaption
///			 
/////////////////////////////////////////////////////////////////
ARM_Swaption* ARM_BasketCalib::CreateVarNotionalSwaptionAtExer (ARM_ZeroCurve* curve,int exerIdx, double fee)
{
	//int i,j,k;
	//
	//double asOfDate	  = curve->GetAsOfDateJul();
	//ARM_Currency* ccy = curve->GetCurrencyUnit();

	//std::vector<double>& exerDates = GetExerDateStrip()->GetResetDates();	
	//double ExerciseDate = (*exerDates)[exerIdx];	

	///// get defaults swaption params for currency
	//int spotDays			 = ccy->GetSpotDays(); 
	//ARM_INDEX_TYPE indexType = ccy->GetVanillaIndexType();
	//char* resetCalendar		 = ccy->GetResetCalName(indexType);
	//char* payCalendar		 = ccy->GetPayCalName(indexType);
	//int stdFixFreq			 = ccy->GetFixedPayFreq();
	//int stdFixDayCount		 = ccy->GetFixedDayCount();
	//int resetGap			 = - spotDays;
	//int fwdRule				 = K_MOD_FOLLOWING;	// for forward dates
 //   int intRule				 = K_ADJUSTED;			// for interest dates
 //   int stubRule			 = K_SHORTSTART;
	//int resetTiming			 = K_ADVANCE;
	//int payTiming			 = K_ARREARS;
	//	
	//// compute start date
	//ARM_Date startDate(ExerciseDate);
	//startDate.GapBusinessDay(spotDays,resetCalendar);
	///// find (unadjusted) end date for variable notional swaption
	//ARM_Date realEndDate = GetRealEndDate();
	//ARM_Date endDate(startDate), prevEndDate;
	//while (endDate<realEndDate)
	//{
	//	prevEndDate = endDate;
	//	endDate.AddMonths( 12/stdFixFreq ) ;
	//}
	///*if (fabs(prevEndDate.GetJulian()-realEndDate.GetJulian()) < fabs(endDate.GetJulian()-realEndDate.GetJulian()) )
	//	endDate = prevEndDate;*/


	///// create date strip
	//ARM_DateStrip dateStrip  (	startDate,endDate,stdFixFreq,stdFixDayCount,
	//							resetCalendar,
	//							fwdRule,intRule,stubRule,
	//							resetGap,stdFixFreq, GETDEFAULTVALUE, payCalendar, resetTiming, payTiming);	
	//std::vector<double>& resetDates	 = dateStrip.GetResetDates();
	//std::vector<double>& startDates	 = dateStrip.GetFlowStartDates();
	//std::vector<double>& endDates		 = dateStrip.GetFlowEndDates();
	//std::vector<double>& payDates		 = dateStrip.GetPaymentDates();
	//std::vector<double>& interestTerms = dateStrip.GetInterestTerms();
	//int sizeSched = startDates->size();
	//double startDf = curve->DiscountPrice(((*startDates)[0]-asOfDate)/K_YEAR_LEN);
	//
	///// instanciate stripper and compute start df + fwd swap rates to feed stripper
	//ARM_STRIPPER stripper(asOfDate, ARM_DateStripPtr(new ARM_DateStrip(dateStrip)));
	//stripper.setStartDf(startDf);

	//std::vector<double> forwards(sizeSched);
	//double level=0.,df;

	//for (j=0; j<sizeSched; j++)
	//{
	//	df = curve->DiscountPrice(((*endDates)[j]-asOfDate)/K_YEAR_LEN);
	//	level += df * (*interestTerms)[j];
	//	forwards[j] = (startDf - df) / level;
	//	stripper.setSwapRate (j, forwards[j]);
	//}
	//stripper.strip();

	///// create notional
	//ARM_Vector fixAbs(sizeSched,1.);
	//ARM_Vector fixNotios(sizeSched,1.);
	//double fixNotional;
	//bool AllZero = FillNotionals(fixAbs,fixNotios,*endDates,*payDates,fixNotional);
	//if (AllZero) 
	//	return NULL;

	///// initial values for dfs and swaprates
	//for (k = 0; k < itsStructLeg.size() ; k++)
	//	for (i = itsStructLeg[k]->GetCpnIndex(exerIdx); i<itsStructLeg[k]->GetCpnSize(); i++)
	//		itsStructLeg[k]->Init(i,stripper);

	///// compute underlyingPv
	//double underlyingPv (0.0);
	//double feeRompu(0.0);

	//for (k = 0; k < itsStructLeg.size() ; k++)
	//{
	//	for (i = itsStructLeg[k]->GetCpnIndex(exerIdx); i<itsStructLeg[k]->GetCpnSize(); i++)
	//	{
	//		underlyingPv	+= itsStructLeg[k]->Notional(i) * itsStructLeg[k]->GetCpn(i) * itsStructLeg[k]->GetDF(i) * itsStructLeg[k]->GetW();
	//	}
	//	feeRompu		+= itsStructLeg[k]->GetFee(exerIdx) * itsStructLeg[k]->GetW();
	//}

	//underlyingPv -= ( feeRompu + fee ) * curve->DiscountPrice((ExerciseDate-asOfDate)/K_YEAR_LEN);
 //
	//double dfFee0 =  stripper.df(ExerciseDate-asOfDate);
	//double Numeraire0=0., Numeraire;
	//for (j=0; j<sizeSched; j++)
	//	Numeraire0 += (*interestTerms)[j] * (fixNotios[j]/fixNotional) * stripper.df((*endDates)[j]-asOfDate) ;

	//std::vector<double> NumerizedLiborFlow0 (sizeSched, 0.0);
	//for (j=0; j<sizeSched; j++)
	//	NumerizedLiborFlow0[j] =  (stripper.df((*startDates)[j]-asOfDate) -  stripper.df((*endDates)[j]-asOfDate)) / Numeraire0;
	//
	///// compute deal sensitivities w.r.t fwd swap rates
	//std::vector<double> sensitivities (sizeSched, 0.0);
	//std::vector<double> invNumeraireSensis (sizeSched, 0.0);
	//ARM_GP_Matrix toInvert(sizeSched, sizeSched);
	//
	//double weight = (itsBasketType==SIMPLE?0.:1.);
	//	
	///// main loop on bumped forwards
	//for (j=0; j<sizeSched; j++)
	//{
	//	/// shift fwd swap rate #j
	//	stripper.setSwapRate(j, forwards[j] + SENSI_CMS);
	//	stripper.strip();

	//	/// (1) for structured leg
	//	for (k = 0; k < itsStructLeg.size() ; k++)
	//	{
	//		for (i = itsStructLeg[k]->GetCpnIndex(exerIdx); i<itsStructLeg[k]->GetCpnSize(); i++)
	//		{
	//			df = stripper.df(itsStructLeg[k]->GetPay(i)-asOfDate);
	//			sensitivities[j] += itsStructLeg[k]->Notional(i) *itsStructLeg[k]->GetCpn(i) * (df - itsStructLeg[k]->GetDF(i)) / SENSI_CMS * itsStructLeg[k]->GetW();
	//			sensitivities[j] += itsStructLeg[k]->Notional(i) *itsStructLeg[k]->GetDCpn(i,stripper) * itsStructLeg[k]->GetDF(i) * weight / SENSI_CMS * itsStructLeg[k]->GetW();
	//		}
	//	}

	//	/// take fee into account
	//	df = stripper.df(ExerciseDate-asOfDate);
	//	sensitivities[j] -= ( feeRompu + fee ) * (df - dfFee0);
	//

	//	/// (2) for 1/Numeraire
	//	Numeraire = 0.0;
	//	for (size_t k=0; k<sizeSched; k++)
	//	{	
	//		df = stripper.df((*endDates)[k]-asOfDate);
	//		Numeraire += (*interestTerms)[k] * (fixNotios[k]/fixNotional) * df;
	//	}

	//	for (k=0; k<sizeSched; k++)
	//	{
	//		toInvert(j,k)  =  (stripper.df((*startDates)[k]-asOfDate) -  stripper.df((*endDates)[k]-asOfDate)) / Numeraire;
	//		toInvert(j,k) -= NumerizedLiborFlow0[k];
	//		toInvert(j,k) /= SENSI_CMS;
	//		toInvert(j,k) *= itsSide;
	//	}
	//	
	//	invNumeraireSensis[j] = (1./Numeraire - 1./Numeraire0) / SENSI_CMS;

	//	stripper.setSwapRate(j, forwards[j]);
	//}

	///// Reset stripper
	//stripper.strip();

	/////
	///// build matrix to invert
	///// for the moment its contains the sensis of the Underlying w.r.t swap rates
	///// now, we want it to contain the sensi of Underlying/Numeraire w.r.t. swap rates
	//for (j=0; j<sizeSched; j++)
	//{
	//	sensitivities[j] *= 1./Numeraire0;
	//	sensitivities[j] += underlyingPv * invNumeraireSensis[j];
	//}
	//
	///// invert matrix
	//LinSolve(&toInvert, &sensitivities);
	//ARM_Vector floatNotios (sizeSched,1.);
	//for (j=0; j<sizeSched; j++)
	//	floatNotios[j] = sensitivities[j];
	//
	//	
	///// compute strike
	//double strike=0;
	//double varFloatLeg = 0.0;

	//for (j=0; j<sizeSched; j++)
	//{			
	//	varFloatLeg += floatNotios[j] * (stripper.df((*startDates)[j]-asOfDate) - stripper.df((*endDates)[j]-asOfDate));
	//}

	//if (itsBasketStrike == ATM)
	//	underlyingPv = 0.0;

	//strike = (itsSide * varFloatLeg - underlyingPv) / (Numeraire0 * fixNotional) ;


	/////  vector need to be cloned (deleted in ARM_ReferenceValue destructor)
	//ARM_ReferenceValue floatNotionalRefVal ((ARM_Vector*)fixAbs.Clone(), (ARM_Vector*)floatNotios.Clone());
	//floatNotionalRefVal.SetCalcMethod(K_STEPUP_RIGHT);

	//	///  vectors need to be cloned (deleted in ARM_ReferenceValue destructor)
	//ARM_ReferenceValue fixNotionalRefVal((ARM_Vector*)fixAbs.Clone(), (ARM_Vector*)fixNotios.Clone());
	//fixNotionalRefVal.SetCalcMethod(K_STEPUP_RIGHT);
	//	
	//// index type: same frequency as fixed leg
	//ARM_IRIndex irIndex (indexType, stdFixFreq, stdFixFreq, ccy);
	//irIndex.SetTerm(stdFixFreq);
	//irIndex.SetYearTerm(1.0/stdFixFreq);

	//ARM_SwapLeg armFloatLeg( startDate, 
	//						 endDate, 
	//						 &irIndex, 
	//						 itsSide==1.?K_PAY:K_RCV,
	//						 0.0, 
	//						 K_SHORTSTART, 
	//						 K_COMP_PROP,
	//						 ccy,
	//						 ccy->GetLiborIndexDayCount());

	///// floatNotionalRefVal is cloned in SetAmount
	//armFloatLeg.SetAmount(&floatNotionalRefVal);


	//// create fixed leg with good strike
	//ARM_FixLeg armFixLeg ( startDate,
	//					   endDate, 
	//					   strike * 100.0,
	//					   itsSide==1.?K_RCV:K_PAY,
	//					   stdFixFreq,
	//					   stdFixDayCount,
	//					   K_COMP_PROP,
	//					   K_ARREARS, 
	//					   K_ADJUSTED,
	//					   K_SHORTSTART,
	//					   ccy);
	//
	///// fixNotionalRefVal is cloned in SetAmount
	//armFixLeg.SetAmount(&fixNotionalRefVal);
	//
	//
	///// create swap
	//ARM_Swap swap(&armFixLeg, &armFloatLeg);
	//ARM_Date expiryDate (ExerciseDate);
	//		
	/// create swaptionj
	ARM_Swaption* swaption = NULL;//new ARM_Swaption((&swap), K_RCV,K_EUROPEAN, strike * 100.0, expiryDate);

	return swaption;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BasketCalib
///	Routine: FillNotionals
///	Returns: bool
///	Action : computes notionals from inputed reference value, with sanity checks
/////////////////////////////////////////////////////////////////
bool ARM_BasketCalib::FillNotionals(ARM_Vector& fixAbs,ARM_Vector& fixNotios,const std::vector<double>& endDates,const std::vector<double>& payDates,double& fixNotional)
{
	int j;
	int sizeSched=fixAbs.size();
/// check that notionals are all > 0 (as snv requires so)
	ARM_Vector* notio = itsNotional.GetDiscreteValues();
	for (j=0; j<notio->size(); j++)
	{
		if (notio->Elt(j)<0)
			ARM_THROW( ERR_INVALID_ARGUMENT, "CreateVarNotionalSwaptionAtExer : negative notional" );
	}


	// compute fixIndex so that we are as close as possible from product end date
	size_t fixIndex (0);

	if (endDates[sizeSched-1]<GetEndDate().GetJulian())
		ARM_THROW( ERR_INVALID_ARGUMENT, "CreateVarNotionalSwaptionAtExer : unresolved date problem" );

	while (endDates[fixIndex]<GetEndDate().GetJulian())
		fixIndex ++;


	if (	(fixIndex>0)
		&&  ( endDates[fixIndex] - GetEndDate().GetJulian() >= GetEndDate().GetJulian()  - endDates[fixIndex-1] )  )
		fixIndex = fixIndex - 1;
	
	for (j=0; j<sizeSched; j++)
	{
		fixAbs[j] = endDates[j];
		
		if (j>fixIndex)
			fixNotios[j] = 0.0;
		else
		{	/// rather approximative ...
			fixNotios[j] = itsNotional.Interpolate(payDates[j]); 
			
			/// take 0 if <0 (with lin interp, result can be sometimes <0 for small notio values)
			if (fixNotios[j]<0)
				fixNotios[j] = 0.0;
		}
	}

		
	/// if fixNotios are all 0, exclude this swaption of calib set
	bool AllZero = true;
	for (j=0; j<sizeSched; j++)
	{
		if (fixNotios[j])
		{
			AllZero = false;
			break;
		}
	}

	/// Notional used for normalization
	for (j=0; j<sizeSched; j++)
	{
		if (fixNotios[j])
		{
			fixNotional = fixNotios[j];
			break;
		}
	}

	return AllZero;
}

////////////////////////////////////////////////////
///	Class  : ARM_BasketCalib
///	Routine: Clone method
///	Returns: a new copy of this
///	Action : Copy this
////////////////////////////////////////////////////
ARM_Object* ARM_BasketCalib::Clone() const
{
	return new ARM_BasketCalib( *this );
}

////////////////////////////////////////////////////
///	Class  : ARM_BasketCalib
///	Routine: View
///	Returns: void
///	Action : standard view
////////////////////////////////////////////////////
void ARM_BasketCalib::View(char* id, FILE* ficOut) const
{
    FILE* fOut;
    char fOutName[200];
	
    if ( ficOut == NULL )
    {
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w");
    }
    else
		fOut = ficOut;   

    fprintf(fOut, "\n INFOS FOR BASKET CALIB\n" );
	fprintf(fOut, "%s", toString( " | \t" ).c_str() ); 
    fprintf(fOut, " | \n" );
    fprintf(fOut, "\n Basket weights computed by Market Model \n" );
	fprintf(fOut, " | \n" );
    fprintf(fOut, "%s", itsBasketCoefs.toString( " | \t" ).c_str() ); 
    fprintf(fOut, " | \n" );

    
    /*if(itsPortfolio!=ARM_StdPortfolioPtr(NULL))
		itsPortfolio->View(id,fOut);*/

    fprintf(fOut, " | \n" );

    fprintf(fOut, " ======> END OF INFOS BASKET CALIB <================== \n\n" );
    
    if ( ficOut == NULL )
       fclose(fOut);
}


////////////////////////////////////////////////////
///	Class   : ARM_BasketCalib
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_BasketCalib::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Basket Calibration\n";
    os << indent << "Method : \t" << ARM_ArgConvReverse_BasketCalibrationType.GetString(itsBasketType)<< " \n";
	os << indent << "Strike : \t" << ARM_ArgConvReverse_BasketCalibrationStrike.GetString(itsBasketStrike)<< " \n";
	
	return os.str();	
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

