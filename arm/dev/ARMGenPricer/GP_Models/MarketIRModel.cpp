/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file MarketIRModel.cpp
 *
 *  \brief class for analytic pricing using ARM_PricingModel's market data manager
 *	\author  A. Chaix
 *	\version 1.0
 *	\date July 2005
 */

/// remove identified warning
#include "gpbase/removeidentifiedwarning.h"

#include "gpmodels/marketirmodel.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/mktdatamanagerrep.h"
#include "gpinfra/nullmodelparams.h"

/// gpbase
#include "gpbase/surface.h"
#include "gpbase/datestrip.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/vanilla_normal.h"

/// gpmodels

/// kernel
#include <inst/swaption.h>
#include <crv/volint.h>
#include <crv/volcube.h>




CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_MarketIRModel
///	Routine: Constructor
///	Returns: 
///	Action : builds the ARM_MarketIRModel
////////////////////////////////////////////////////

ARM_MarketIRModel::ARM_MarketIRModel(const ARM_MarketData_ManagerRep& mktDataManager,
									 const ARM_StringVector& mdmKeys, 
									 VnsPricingMethod method,
									  double moyenessLevel)
:	
	ARM_AnalyticIRModel( ARM_ZeroCurvePtr(new ARM_ZeroCurve())/*fake zero curve*/,	ARM_NullModelParams	()/*fake model params*/ ),
		itsVnsPricingMethod (method),
		itsMoyenessLevel(moyenessLevel)
{
	itsVnsBasketCorrels.reserve(30,30);
	SetMktDataManager(mktDataManager, mdmKeys);

	// default choice for mkt data
	//  o take as zero curve the first zero curve appearing in mkt data manager
	//  o taka as BS model the first BS modelo appearing in mkt data manager
	// use SetZeroCurveKey and SetBsModelKey methods for modification
	size_t i;
	ARM_ZeroCurve* zcCurve;
	for (i=0; i<GetKeys().size(); i++)
	{
		if ( zcCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[i])) )
		{
			itsZeroCurveKey = GetKeys()[i];
			break;
		}
		if (i == GetKeys().size() - 1)
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MarketIRModel: market data manager is required to contain at least one zero curve" );
			
	}
	for (i=0; i<GetKeys().size(); i++)
	{
		if ( dynamic_cast<ARM_BSModel*>(GetMktDataManager()->GetData(GetKeys()[i])) )
		{
			itsBsModelKey = GetKeys()[i];
			break;
		}
		if (i == GetKeys().size() - 1)
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MarketIRModel: market data manager is required to contain at least one BS model" );
			
	}

	// don't forget to set ZeroCurve Properly
	SetZeroCurve( ARM_ZeroCurvePtr((ARM_ZeroCurve*)zcCurve->Clone()) );
	SetModelName(zcCurve->GetCurrencyUnit()->GetCcyName());
}


////////////////////////////////////////////////////
///	Class   : ARM_MarketIRModel
///	Routines: Copy constructor
///	Returns :
///	Action  : 
////////////////////////////////////////////////////
ARM_MarketIRModel::ARM_MarketIRModel(const ARM_MarketIRModel& rhs)
:	ARM_AnalyticIRModel( rhs ),
	itsZeroCurveKey (rhs.itsZeroCurveKey),
	itsBsModelKey	(rhs.itsBsModelKey),
	itsVnsPricingMethod (rhs.itsVnsPricingMethod),
	itsMoyenessLevel(rhs.itsMoyenessLevel),
	itsVnsBasketCorrels(rhs.itsVnsBasketCorrels)
{}

////////////////////////////////////////////////////
///	Class  : ARM_MarketIRModel
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_MarketIRModel::~ARM_MarketIRModel()
{}


////////////////////////////////////////////////////
///	Class  : ARM_MarketIRModel
///	Routine: VanillaCaplet
///	Returns: ARM_VectorPtr
///	Action : computes the price of a caplet
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MarketIRModel::VanillaCaplet(
	const string& curveName, 
	double evalTime,
	double payTime, 
	double period,
    double payNotional,
	double fwdResetTime, 
	double fwdStartTime,
    double fwdEndTime,
	double fwdPeriod,
	const ARM_GP_Vector& strikesPerState,
    int capFloor,
	const ARM_PricingStatesPtr& states) const
{
	ARM_THROW(ERR_INVALID_ARGUMENT, "ARM_MarketIRModel::VanillaCaplet : not implemented" );
	return ARM_VectorPtr(NULL);
}	

////////////////////////////////////////////////////
///	Class   : none / local function
///	Routines: ComputeDiscountFactorsFromSwapRates
///	Returns : void
///	Action  : Stripping of fwd swap rates into DF
///			  Used for variable notional swaption
///			  price
////////////////////////////////////////////////////
void ComputeDiscountFactorsFromSwapRates (	const ARM_GP_Vector& fixPayPeriods,
											const ARM_GP_Vector& swapRates,
											double startDf, 
											/// result
											ARM_GP_Vector& dfs)
											   
{	size_t size = swapRates.size();
	if (size != dfs.size()) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MarketIRModel / Variable Notional Swaption:  bad size" );

	double Sratio, df;

	for (size_t i(0); i<size; i++)
	{
		if (i==0)
			dfs[i] = startDf / ( 1.0 + fixPayPeriods[i] * swapRates[i] ) ;
		else
		{
			Sratio = swapRates[i] / swapRates[i-1];
			df  = Sratio * dfs[i-1] + (1.0-Sratio) * startDf;
			df /= 1.0 + fixPayPeriods[i] * swapRates[i];
			dfs[i] = df;
		}
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_AnalyticIRModel
///	Routines: VanillaSwaption
///	Returns : ARM_VectorPtr
///	Action  : Vanilla Swaption : standard + variable notional
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MarketIRModel::VanillaSwaption(
		const string& curveName,
		double evalTime,
		double swapResetTime,
		const ARM_GP_Vector& fixNotional,
		const ARM_GP_Vector& floatNotional,
		double floatStartTime,
		double floatEndTime,
		const ARM_GP_Vector& floatResetTimes,
		const ARM_GP_Vector& floatStartTimes,
		const ARM_GP_Vector& floatEndTimes,
		const ARM_GP_Vector& floatIntTerms,	
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
		const ARM_PricingStatesPtr& states,
		bool isConstantNotional,
		bool isConstantSpread,
		bool isConstantStrike) const
{
	ARM_ZeroCurve* curve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(itsZeroCurveKey));
	ARM_BSModel* bsModel = dynamic_cast<ARM_BSModel*>  (GetMktDataManager()->GetData(itsBsModelKey));

	if (!curve || !bsModel)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MarketIRModel: problem with market data manager" );
	
	ARM_VolCurve* correlCube  = bsModel->GetCorrelationMatrix();

	if (!correlCube)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MarketIRModel: invalid correl cube" );

		
	/// some validations...
	if (!isConstantSpread)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MarketIRModel: variable spread not supported)" );

	if (!isConstantStrike)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MarketIRModel: variable strike not supported)" );

	if( !GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
        /// We need to compute the floating leg by forward method and no more by double notional
        /// but we have not at the moment all floating leg datas => throw an error
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "Swaption pricing not implemented for differents discount & fixing curves" );

	if (fixPayTimes.size() != fixPayPeriods.size() || fixPayTimes.size() != fixNotional.size())
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MarketIRModel / Variable Notional Swaption: invalid sizes." );

	if (states->size() != 1)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MarketIRModel / Variable Notional Swaption:  multi state not supported" );

	if (strikesPerState.rows() != 1)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MarketIRModel / Variable Notional Swaption:  strike per state not supported" );

	if (strikesPerState.cols() != fixPayTimes.size()) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MarketIRModel / Variable Notional Swaption:  strikesPerState: bad size" );
	
	if ( fabs(floatEndTime - fixPayTimes[fixPayTimes.size()-1])>0.001)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MarketIRModel / Variable Notional Swaption: float and fixed end dates are not matching" );

	if (floatNotional.size() != fixNotional.size())
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MarketIRModel / Variable Notional Swaption: float an fix legs are required to be of same frequency." );


	size_t size, i, j;
	size = fixPayTimes.size();
	
	double NPV, Level, Result, Df, DfStart,
		   BasketVol, BasketForward, Correl, 
		   Tmp, nStdev, Strike_i, Expiry, SqrtTe, 
		   Shift, DfCoefStart, X, X0, Numeraire, Numeraire0;

	ARM_GP_Vector DfCoefs (size+1), Dfs(size);
	itsVnsStrikes.resize(size);
	itsVnsForwards.resize(size);
	itsVnsAtmNormVols.resize(size);
	itsVnsNormVols.resize(size);
	itsVnsBasketCoefs.resize(size);
	itsVnsStartDates.resize(size);
	itsVnsEndDates.resize(size);
	
	/// for toString display
	itsVnsFloatNotionals = floatNotional;
	itsVnsStrike		 = strikesPerState(0,0);

	for (i = 0; i<size; i++)
		itsVnsEndDates[i] = GetDateFromTime(fixPayTimes[i]);
	
	itsVnsStartDates[0] = GetDateFromTime(floatStartTime);
	for (i = 1; i<size; i++)
		itsVnsStartDates[i] = GetDateFromTime(fixPayTimes[i-1]);;
	
		
	DfCoefStart = floatNotional[0];

	for (i = 0; i<size-1; i++)
		DfCoefs[i] = floatNotional[i+1] - floatNotional[i];
	DfCoefs[size - 1] = - floatNotional[size-1];

				
	// Compute Level of swap with variable notional
	Level	 = 0.0;
	DfStart = curve->DiscountPrice(floatStartTime/K_YEAR_LEN);
		
	for (i=0; i<size; i++)
	{
		Df = curve->DiscountPrice(fixPayTimes[i]/K_YEAR_LEN);;
		Level       += Df * fixPayPeriods[i];
		itsVnsForwards[i] = (DfStart - Df) / Level;
	}
	
	
	// compute X0
	X0 = 0.0;
	ComputeDiscountFactorsFromSwapRates(fixPayPeriods, itsVnsForwards, DfStart, Dfs);	

	X0 += DfCoefStart * DfStart;
	for (j=0; j<size; j++)
		X0 += DfCoefs[j] * Dfs[j];
	
	// Modif : Var notio swaption = float leg notios & fixed leg notios
	Numeraire0 = 0.0;
	
	for (j=0; j<size; j++)
	{
		if (fixNotional[j]<0.0)
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MarketIRModel / Variable Notional Swaption: Fixed leg is required to have only positive notionals");
			
		Numeraire0 += fixPayPeriods[j] * fixNotional[j] * Dfs[j];
	}

	if (Numeraire0 < 1e-12)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MarketIRModel / Variable Notional Swaption: At least one fixed leg notionals should be > 0");

	itsVnsNumeraire = Numeraire0;
		
	X0 /= Numeraire0;

	Shift = 0.00001; // 0.1 bp
	
	for (i=0; i<size; i++)
	{
		// shift forward swap rate
		itsVnsForwards[i] += Shift;
		
		ComputeDiscountFactorsFromSwapRates(fixPayPeriods, itsVnsForwards, DfStart, Dfs);
			
		// compute X
		X  = 0.0;
		X += DfCoefStart * DfStart;
		for (j=0; j<size; j++)
			X += DfCoefs[j] * Dfs[j];

		Numeraire = 0.0;
		
		for (j=0; j<size; j++)
			Numeraire += fixPayPeriods[j] * fixNotional[j] * Dfs[j];
			
		X /= Numeraire;

		// compute VolCoeff
		itsVnsBasketCoefs[i] =  (X - X0) / Shift;
		
		// RAZ
		itsVnsForwards[i] -= Shift;
	}
	
	
	// fwd of variable-notional swap rate
	BasketForward = X0;
	
	/// for display in toString method
	itsVnsForward = X0;

		
	// Time to expiry
	Expiry = swapResetTime / K_YEAR_LEN;
	SqrtTe = sqrt(Expiry);
		
	for (i=0; i<size; i++)
	{
		if (fabs(itsVnsBasketCoefs[i])>1e-15)
		{
			// ATM -> on fait fwd * bsvol à la place de normal vol pour gagner du temps (evite la conversion vol lognorm -> vol norm)
			double tenor  = (fixPayTimes[i] - floatStartTime) / K_YEAR_LEN;
			double vol;

			// compute lognormal atm vol
			vol = bsModel->ComputeVol(Expiry, tenor, 100. * itsVnsForwards[i], 100. * itsVnsForwards[i], K_LIBOR);
			vol *= 0.01;

			itsVnsAtmNormVols[i] = itsVnsForwards[i] * vol ;
		}
		else
			itsVnsAtmNormVols[i] = 0.0;
	}
		
	// compute variance of basket
	BasketVol = 0.0;

	if (itsVnsPricingMethod == MONEYNESS) 
	{
		for (i=0; i<size; i++)
		{
			if (fabs(itsVnsBasketCoefs[i])>1e-15)
			{
				// diag terms
				Tmp = itsVnsBasketCoefs[i] * itsVnsAtmNormVols[i];
				BasketVol += Tmp * Tmp;
				double tenor_i = (fixPayTimes[i] - floatStartTime) / K_YEAR_LEN;

				// non diag terms
				for (j=0; j<i; j++)
				{
					if (fabs(itsVnsBasketCoefs[j])>1e-15)
					{					
						double tenor_j = (fixPayTimes[j] - floatStartTime) / K_YEAR_LEN;
						Correl  = correlCube->ComputeCorrelByExpiry(Expiry, tenor_i, tenor_j);
						Correl *= 0.01;
						if (Correl > 1)			Correl =  1;
						else if (Correl < -1)	Correl = -1;
						
							
						BasketVol += 2.0 * Correl 
										 * itsVnsBasketCoefs[i] * itsVnsAtmNormVols[i]
										 * itsVnsBasketCoefs[j] * itsVnsAtmNormVols[j];
					}
				}
			}
		}

		BasketVol = sqrt (BasketVol);
		if (BasketVol<1e-15) 
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MarketIRModel / Variable Notional Swaption: unresolved problem");
	}

	// compute strikes & build itsVnsNormVols
	double Strike = strikesPerState(0,0);
						
	for (i=0; i<size; i++)
	{
		if (fabs(itsVnsBasketCoefs[i])>1e-15)
		{							
			if (itsVnsPricingMethod == MONEYNESS)
			{
				nStdev = (Strike - BasketForward) /(BasketVol * SqrtTe);
				if( nStdev > 6.0) nStdev = 6.0;
				if( nStdev <- 6.0) nStdev = -6.0;

				if (itsVnsBasketCoefs[i]<0)
					nStdev = -nStdev; 
									
				// compute strike
				Strike_i = itsVnsForwards[i] + itsMoyenessLevel*nStdev * itsVnsAtmNormVols[i] * SqrtTe;
				if (Strike_i<0.0010) Strike_i = 0.0010;
			}
			else if (itsVnsPricingMethod == ATM)
			{
				Strike_i = itsVnsForwards[i];
			}
			else
				ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MarketIRModel / Variable Notional Swaption: invalid VNS Pricing Method");


			itsVnsStrikes[i] = Strike_i;

			double tenor  = (fixPayTimes[i] - floatStartTime) / K_YEAR_LEN;
			double vol;

			// compute lognormal vol @ strike
			
			vol = bsModel->ComputeVol(Expiry, tenor, 100. * itsVnsForwards[i], 100. * Strike_i, K_LIBOR);
			vol *= 0.01;

			// convert into normal vol
			int callPut = (itsVnsForwards[i]>Strike_i) ? -1 : 1 ; // consider OTM option
			double bs = BS(itsVnsForwards[i], Strike_i, Expiry, vol, callPut);
			itsVnsNormVols[i] = VanillaImpliedVol_N(itsVnsForwards[i], bs, Strike_i, Expiry, callPut,&itsVnsAtmNormVols[i]);
	 
		}
		else 
		{
			itsVnsNormVols[i] = 0.0;
		}
	}


	// re-compute basket vol with computed itsVnsNormVols
	BasketVol = 0.0;
	itsVnsBasketCorrels.resize(size,size);
	for (i=0; i<size; i++)
	{
		if (fabs(itsVnsBasketCoefs[i])>1e-15)
		{
			// diag terms
			Tmp = itsVnsBasketCoefs[i] * itsVnsNormVols[i];
			BasketVol += Tmp * Tmp;

			double tenor_i = (fixPayTimes[i] - floatStartTime) / K_YEAR_LEN;

			itsVnsBasketCorrels(i,i)=1;

			// non diag terms
			for (j=0; j<i; j++)
			{
				if (fabs(itsVnsBasketCoefs[j])>1e-15)
				{

					double tenor_j = (fixPayTimes[j] - floatStartTime) / K_YEAR_LEN;
					Correl  = correlCube->ComputeCorrelByExpiry(Expiry, tenor_i, tenor_j);
					Correl *= 0.01;
					if (Correl > 1)			Correl =  1;
					else if (Correl < -1)	Correl = -1;
										
					BasketVol += 2.0 * Correl 
									 * itsVnsBasketCoefs[i] * itsVnsNormVols[i]
									 * itsVnsBasketCoefs[j] * itsVnsNormVols[j];

					/// for display
					itsVnsBasketCorrels(i,j)=Correl;
					itsVnsBasketCorrels(j,i)=Correl;
				}
			}
		}
	}
	BasketVol = sqrt (BasketVol);
	
	if (BasketVol<1e-15) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MarketIRModel / Variable Notional Swaption: unresolved problem");

	
	/// for display
	itsVnsBasketVol = BasketVol;

	// on appelle la formule normale
	Result = VanillaOption_N(BasketForward, BasketVol, Strike, Expiry, callPut);
		
	NPV = Numeraire0 * Result ;
	
	return ARM_VectorPtr(new ARM_GP_Vector(1,NPV));
}


    
////////////////////////////////////////////////////
///	Class   : ARM_MarketIRModel
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_MarketIRModel::Clone() const
{
	return new ARM_MarketIRModel(*this);
}

////////////////////////////////////////////////////
///	Class   : ARM_MarketIRModel
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_MarketIRModel::toString(const string& indent,const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "IR Market Model \n";
    os << indent << "---------------------\n";
    os << ARM_PricingModel::toString(indent);



	if (itsVnsStartDates.size())
	{
		os << "\n\n";

		os << "======================================================\n";
		os << "Details on last priced Variable Notional Swaption\n";
		os << "======================================================\n";

		ARM_Date vnsStartDate = itsVnsStartDates[0];
		ARM_Date vnsEndDate   = itsVnsEndDates[itsVnsEndDates.size()-1];


		os << "Start Date    : " << vnsStartDate.toString() << "\n";
		os << "End Date      : " << vnsEndDate.toString()   << "\n";
		os << "Strike (%)    : " << CC_NS(std,setw)(4)   << CC_NS(std,setprecision)(4) << itsVnsStrike    * 100.0 <<  "\n";
		os << "Forward (%)   : " << CC_NS(std,setw)(4)   << CC_NS(std,setprecision)(4) << itsVnsForward   * 100.0 <<  "\n";
		os << "BasketVol (%) : " << CC_NS(std,setw)(4)   << CC_NS(std,setprecision)(4) << itsVnsBasketVol * 100.0 <<  "\n";

		os << "\n\nStart Date\tEnd Date\tFwd\tStrike\tCoef\tAtmVol\tVolUsed\tCorrels\n";
		os << "---------------------------------------------------------------------------------------------\n";

		for (size_t i = 0; i < itsVnsEndDates.size(); ++i )
		{
			
			os  <<  itsVnsStartDates[0].toString()						<< "\t" 
				<<  itsVnsEndDates[i].toString()						<< "\t" 
				<<  CC_NS(std,setw)(4)   << CC_NS(std,setprecision)(4) 
				<<  100.0 * itsVnsForwards[i]							<< "\t" 
				<<  100.0 * itsVnsStrikes[i]							<< "\t" 
				<<  CC_NS(std,setw)(3)   << CC_NS(std,setprecision)(3) 
				<<  0.001 * floor(1000.0*itsVnsBasketCoefs[i])			<< "\t" 
				<<  100.0 * itsVnsAtmNormVols[i]						<< "\t"
				<<  100.0 * itsVnsNormVols[i]							<< "\t";
			for (size_t j = 0; j < itsVnsEndDates.size(); ++j )
			{
				os	<<	100 * itsVnsBasketCorrels(i,j) << "\t";
			}
			os << "\n";

		}

		os << "\n\nFloat Leg Notionals \n";
		os <<    "------------------------\n";
		os << itsVnsFloatNotionals.toString()<< "\n";

	}
	
    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_MarketIRModel
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_MarketIRModel::ValidateModelParams(const ARM_ModelParams& params) const
{	
	return true;
}



CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

