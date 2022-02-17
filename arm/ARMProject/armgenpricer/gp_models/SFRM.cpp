/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file SFRM.cpp
 *  \brief SFRM model implements the SFRM model
 *	the analytic part and the specialised part on numerical
 *	method
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 */


/// this header comes first as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

/// gpmodels
#include "gpmodels/SFRM.h"
#include "gpmodels/ModelParamsSFRM.h"
#include "gpmodels/VanillaSwaptionArgSFRM.h"
#include "gpmodels/VanillaCapArgSFRM.h"
#include "gpmodels/VanillaDigitalArgSFRM.h"
#include "gpmodels/VanillaCorridorLegArgSFRM.h"
#include "gpmodels/SFRMPricingStates.h"
#include "gpmodels/VanillaSpreadOptionletArgSFRM.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/autocleaner.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/gpmatrixtriangular.h"
#include "gpbase/comparisonfunctor.h"
#include "gpbase/datestrip.h"
#include "gpbase/stringmanip.h"
#include "gpbase/cloneutilityfunc.h"
#include "gpbase/utilityport.h"
#include "gpbase/curve.h"
#include "gpbase/eventviewerfwd.h"
#include "gpbase/gpmatrixlinalg.h"
#include "gpbase/interpolatorvector.h"
#include "gpbase/globalconstant.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/timeinfo.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/discretisationscheme.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/modelparam.h"
#include "gpmodels/HW2F.h"

#include "gpinfra/gensecurity.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/curvemodelparam.h"

/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/modelfitter.h"
#include "gpcalib/kerneltogp.h"
#include "gpcalib/vanillaarg.h"
#include "gpcalib/vanillaspreadoption.h"

/// gpclosedforms
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/normal.h"

/// gpnummethods
#include "gpnummethods/treemethod.h"

/// gpnumlib
#include "gpnumlib/normalinvcum.h"
#include "gpnumlib/gaussiananalytics.h"

/// kernel
#include <inst/irindex.h>
#include <inst/portfolio.h>
#include <inst/swaption.h>
#include <inst/corridorleg.h>
#include <inst/armdigital.h>
#include <inst/spreadoption.h>
#include <inst/sumopt.h>
/// kernel
//#include <util/irindex.h>

/// STL
#include <vector>
CC_USING_NS(std,vector)
#include <memory>
CC_USING_NS(std,auto_ptr)
#include <iomanip>
#include <functional>
#include <algorithm>
#include <new>
#include <new.h>



CC_BEGIN_NAMESPACE( ARM )

const bool ARM_SFRM::ARM_SFRM_STD_FWD = false;
const bool ARM_SFRM::ARM_SFRM_AUX_FWD = true;

const int INTERMEDIATESTEPS = 2; // 2 steps per year 

const int ProbaChangePolyDegre = 6;
const double FRMVOL_LASTLAG_THRESHOLD = 10;

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_SFRM::ARM_SFRM( const ARM_ZeroCurvePtr& zc, const ARM_ModelParamsSFRM& params, ARM_Portfolio* ConvPort, bool diffusionDrift )
:	ARM_PricingModelIR(zc, &params ), 
	itsFwdResetTimes(NULL),
	itsFwdStartTimes(NULL), 
	itsFwdEndTimes(NULL),
	itsFwdValues(NULL),
	itsFwdIntTerms(NULL),
	itsFwdShifts(NULL),
	itsFwdStatus(0),
	itsIsFixStartDate(false),
	itsIsFixEndDate(false),
	itsFactors(NULL),
    itsCurrentArg(NULL),
	itsFwdsCorrelVarCoeffs(NULL),
	itsFwdsVarCoeffs(NULL),
	itsFwdDrifts(NULL),
#ifdef __GP_SFRM_STOCH_TERMS
	itsUpdateVolSwapvolFRA(false),
	itsStochasticTerms(NULL)
#else
	itsUpdateVolSwapvolFRA(false)
#endif
{
	/// validation to see that the model params are compatible!
	ValidateModelParams(params);
	size_t factorsNb= params.FactorCount();
	itsFactors		= ARM_VectorPtr( new std::vector<double>(factorsNb) );
	ConvertToShiftorBetaParam(*ConvPort);
    CC_ARM_SETNAME(ARM_SFRMMODEL);
	ResetDFMap();

// FIXMEFRED: mig.vc8 (25/05/2007 15:52:36): pointer on function
	if (diffusionDrift)
		itsTreeStatesToModelStatesFunc = &ARM::ARM_SFRM::TreeStatesToModelStatesNonParamDrift;
	else
		itsTreeStatesToModelStatesFunc = &ARM::ARM_SFRM::TreeStatesToModelStatesDriftMarkovian;
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_SFRM::ARM_SFRM(const ARM_SFRM& rhs)
:	ARM_PricingModelIR(rhs), 	
	itsFwdResetTimes(rhs.itsFwdResetTimes),
	itsFwdStartTimes(rhs.itsFwdStartTimes),
	itsFwdEndTimes(rhs.itsFwdEndTimes),
	itsFwdValues(rhs.itsFwdValues),
	itsFwdIntTerms(rhs.itsFwdIntTerms),
	itsFwdShifts(rhs.itsFwdShifts),
	itsFwdStatus(rhs.itsFwdStatus),
	itsIsFixStartDate(rhs.itsIsFixStartDate),
	itsFixStartDate(rhs.itsFixStartDate),
	itsIsFixEndDate(rhs.itsIsFixEndDate),
	itsFixEndDate(rhs.itsFixEndDate),
	itsFactors(rhs.itsFactors),
	itsFwdsVarCoeffs(rhs.itsFwdsVarCoeffs),
	itsFwdDrifts(rhs.itsFwdDrifts),
	itsDFMap(rhs.itsDFMap),
#ifdef __GP_SFRM_STOCH_TERMS
	itsUpdateVolSwapvolFRA(rhs.itsUpdateVolSwapvolFRA),
	itsStochasticTerms(NULL),
#else	
	itsUpdateVolSwapvolFRA(rhs.itsUpdateVolSwapvolFRA),
#endif
	itsTreeStatesToModelStatesFunc(rhs.itsTreeStatesToModelStatesFunc)
{
	DuplicateCloneablePointorVectorInPlace<ARM_GP_TriangularMatrix>(rhs.itsFwdsCorrelVarCoeffs, itsFwdsCorrelVarCoeffs);
    itsCurrentArg = rhs.itsCurrentArg ? (ARM_VanillaArg*)rhs.itsCurrentArg->Clone() :  NULL;
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_SFRM::~ARM_SFRM()
{
    delete itsCurrentArg;
    itsCurrentArg = NULL;
#ifdef __GP_SFRM_STOCH_TERMS
	delete itsStochasticTerms;
	itsStochasticTerms = NULL;
#endif

	DeletePointorVector<ARM_GP_TriangularMatrix>( itsFwdsCorrelVarCoeffs );
	
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_SFRM& ARM_SFRM::operator=(const ARM_SFRM& rhs)
{
	if(this != &rhs)
	{
		ARM_PricingModelIR::operator=(rhs);
		itsDFMap				= rhs.itsDFMap;
		itsFwdResetTimes		= rhs.itsFwdResetTimes;
		itsFwdStartTimes		= rhs.itsFwdStartTimes;
		itsFwdEndTimes			= rhs.itsFwdEndTimes;
		itsFwdValues			= rhs.itsFwdValues;
		itsFwdIntTerms			= rhs.itsFwdIntTerms;
		itsFwdShifts			= rhs.itsFwdShifts;
		itsFwdDrifts			= rhs.itsFwdDrifts;
		itsFwdStatus			= rhs.itsFwdStatus;
		itsIsFixStartDate		= rhs.itsIsFixStartDate;
		itsFixStartDate			= rhs.itsFixStartDate;
		itsIsFixEndDate			= rhs.itsIsFixEndDate;
		itsFixEndDate			= rhs.itsFixEndDate;
		itsFactors				= rhs.itsFactors;
		itsFwdsVarCoeffs		= rhs.itsFwdsVarCoeffs;
		itsUpdateVolSwapvolFRA	= rhs.itsUpdateVolSwapvolFRA;
	    itsCurrentArg			= rhs.itsCurrentArg ? (ARM_VanillaArg*)rhs.itsCurrentArg->Clone() :  NULL;

		DeletePointorVector<ARM_GP_TriangularMatrix>( itsFwdsCorrelVarCoeffs );
		DuplicateCloneablePointorVectorInPlace<ARM_GP_TriangularMatrix>(rhs.itsFwdsCorrelVarCoeffs, itsFwdsCorrelVarCoeffs);

#ifdef __GP_SFRM_STOCH_TERMS
		delete itsStochasticTerms;
		itsStochasticTerms		= NULL;
#endif
		itsTreeStatesToModelStatesFunc=rhs.itsTreeStatesToModelStatesFunc;
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: ResetDFMap
///	Returns: 
///	Action : Reset the DF Map
////////////////////////////////////////////////////
void ARM_SFRM::ResetDFMap () const
{
	/// Erase saved functions if any
    itsDFMap = map<double, ARM_VectorPtr>();
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: SetZeroCurve
///	Returns: void
///	Action : To set ZeroCurve and postprocess model
////////////////////////////////////////////////////
void ARM_SFRM::SetZeroCurve( const ARM_ZeroCurvePtr& zc)
{
	if( zc != GetZeroCurve())
	{
		((ARM_ModelParamsSFRM*) GetModelParams())->PostProcessing();
		ResetDFMap();
	}
	ARM_PricingModel::SetZeroCurve(zc);
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: toString
///	Returns: this
///	Action : view of the object
////////////////////////////////////////////////////

string ARM_SFRM::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "=======> SFRM Model " << FactorCount() << "F <======\n";
    os << indent << "------------------------------\n\n";
    os << indent << ARM_PricingModel::toString(indent);

	if (itsFwdResetTimes != ARM_VectorPtr(NULL)
		&& itsFwdStartTimes != ARM_VectorPtr(NULL)
		&& itsFwdEndTimes != ARM_VectorPtr(NULL)
		&& itsFwdIntTerms != ARM_VectorPtr(NULL))
	{
		os << "\n";
		os << indent << " =======> SFRM Scheduler <====== \n";
		os << indent << "Reset Dates\tStartDates\tEndDates\t";
		os << setw(11) << "Fwd\t";
		os << setw(11) << "Term\t";
		os << "Status\n";

		int nbFwds = itsFwdResetTimes->size();

		double asOfDate	= GetAsOfDate().GetJulian();

		os << setiosflags(ios::left);

		for (int i = 0; i < nbFwds; ++i)
		{
			ARM_Date resetDate(asOfDate+(*itsFwdResetTimes)[i]);
			ARM_Date startDate(asOfDate+(*itsFwdStartTimes)[i]);
			ARM_Date endDate(asOfDate+(*itsFwdEndTimes)[i]);

			os	<< indent << resetDate.toString()	<< "\t"
				<< startDate.toString() << "\t"
				<< endDate.toString() << "\t";

			double fwdValue = 0;

			if (itsFwdValues != ARM_VectorPtr(NULL))
				fwdValue = (*itsFwdValues)[i];

			os	<< setfill('0') <<  setw(11) << setprecision(8) << fwdValue << "\t";
			os	<< setfill('0') << setw(11) << setprecision(8) << (*itsFwdIntTerms)[i] << "\t";

			if (itsFwdStatus[i] == ARM_SFRM::ARM_SFRM_STD_FWD)
				os << "STD" << endl;
			else
				os << "AUX" << endl;
		}
	}

	return os.str();	
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: ValidateModelParams
///	Returns: true/false
///	Action : Check the consistency of the model
///          parameters
////////////////////////////////////////////////////
bool ARM_SFRM::ValidateModelParams(const ARM_ModelParams& params) const
{
	const ARM_ModelParamsSFRM* modelParamsSFRM = dynamic_cast<const ARM_ModelParamsSFRM*>(&params);
	if( !modelParamsSFRM )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_ModelParamsSFRM" );

	/// test that the irIndex has the same currency as the curve
	string curveCcy = GetZeroCurve()->GetCurrencyUnit()->GetCcyName();
	string indexCcy = ((ARM_ModelParamsSFRM&) params).GetIRIndex()->GetCurrencyUnit()->GetCcyName();

	if( !StrIEqual(curveCcy,indexCcy) )
	{
		CC_Ostringstream os;
		os << "ir index ccy " << indexCcy << " while curve ccy " << curveCcy;
		throw Exception(__LINE__, __FILE__,ERR_INVALID_ARGUMENT,os.str());
	}

	return true;
}


////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: Init
///	Returns: 
///	Action : Default initialisation of the model and the
///          associated numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_SFRM::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
	size_t factorsNb=FactorCount();
	//itsFactors		= ARM_VectorPtr( new std::vector<double>(factorsNb) );

    int nbEvents=timeInfos.size();
    bool isSpotUse = nbEvents == 0 || (nbEvents==1 && timeInfos[0]->GetEventTime() <= K_NEW_DOUBLE_TOL);
	ResetDFMap();
	
    if(!isSpotUse)
    {
        // Numerical method and numeraire are needed to go on
        ARM_NumMethodPtr numMethod=GetNumMethod();
		if( numMethod == ARM_NumMethodPtr(NULL))
		{
			CC_Ostringstream os;
			os << ARM_USERNAME << ": numerical method not set in the SFRM mode!";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}

        ARM_NumerairePtr numeraire=GetNumeraire();
        if( numeraire == ARM_NumerairePtr(NULL) )
		{
			CC_Ostringstream os;
			os << ARM_USERNAME << ": numeraire not set in the SFRM mode!";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
            
		}

        if(	   numeraire->GetType() != ARM_Numeraire::TerminalZc
			&& numeraire->GetType() != ARM_Numeraire::TerminalEventZc
			&& numeraire->GetType() != ARM_Numeraire::RollingEvent
			&& numeraire->GetType() != ARM_Numeraire::RollingPayment)
        {
			CC_Ostringstream os;
			os  << ARM_USERNAME << ": only TerminalZc, TerminalEventZc, RollingEvent, RollingPayment numeraire"
				<< " supported by SFRM mode at the moment!";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
				
        }
		
        /// creates the model schedule (smart pointor for exception safety!)
		ARM_DiscretisationScheme* discretisationScheme;
				
		if(numMethod->GetPricingDirection() ==  ARM_NumMethod::GP_FWDLOOKING || numMethod->GetPricingDirection() == ARM_NumMethod::GP_FWDBCKWDLOOKING)
        {
			discretisationScheme = new ARM_EventAndModelTime();
		}
		else if(numMethod->GetPricingDirection() ==  ARM_NumMethod::GP_BCKWDLOOKING)
        {
			discretisationScheme = new ARM_EventTime();	
			CC_NS(std,auto_ptr)<std::vector<double>> pModelTimes( (*this).ComputeModelTimes( timeInfos ) );

        }
		CC_NS(std,auto_ptr)<std::vector<double>> ptimeSteps( discretisationScheme->ModelTimesFromTimesInfo(timeInfos, *this) );

		delete discretisationScheme;

        //// Initialise the numeraire
		numeraire->Init(*(GetDiscountFunctor()),payModelName,timeInfos, ComputeNumeraireTimes( timeInfos ) );

		/// Set the basic schedule in the numerical method and...
		numMethod->SetTimeSteps(*ptimeSteps);

		// The first induct time is the first reset date of the generic security
		double firstInductTime = timeInfos[0]->GetEventTime();

		/// ...initialise it
		return numMethod->Init(*this,firstInductTime);
    }
    else
    {
        // Compute a single model states set to (0.0,...,0.0)
        int nbDir = FactorCount();
        ARM_PricingStatesPtr initStates(new ARM_PricingStates(1,nbDir,0));
		initStates->EnqueuePricingStatesContext( ARM_PricingStatesContextPtr( new ARM_SFRMPricingStatesContext () ) );
        for(int i=0;i<nbDir;++i)
            initStates->SetModelState(0,i,0.0);
		
        return initStates;
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: ComputeNumeraireTimes
///	Returns: 
///	Action : compute the corresonding numeraire times
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SFRM::ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const
{
	size_t nbEvents = timeInfos.size();

	switch( GetNumeraire()->GetType() )
	{
	case ARM_Numeraire::TerminalZc:
		{
			double lastPayTime=0;
			/// Set the very last payment
			for(int i=0;i<nbEvents;i++)
			{
				const std::vector<double>& lastPayTimes = timeInfos[i]->GetPayTimes();
				if(lastPayTime<lastPayTimes[lastPayTimes.size()-1])
				lastPayTime	= lastPayTimes[lastPayTimes.size()-1];
			}
			size_t payIndex = CC_NS(std,lower_bound)(itsFwdEndTimes->begin(), itsFwdEndTimes->end(), lastPayTime-K_NEW_DOUBLE_TOL) - itsFwdEndTimes->begin();

			// test if the forward is auxiliary
			while( payIndex < itsFwdResetTimes->size() && itsFwdStatus[payIndex] == ARM_SFRM::ARM_SFRM_AUX_FWD )
				++payIndex;

			/// test if the fwd is dead?
			while( payIndex < itsFwdResetTimes->size() && (*itsFwdResetTimes)[payIndex] < timeInfos[nbEvents-1]->GetEventTime() )
				++payIndex;
			
			/// if it is a non std fwd, goes to the nest std fwd... we are sure the last fwd is std!
			while( payIndex < itsFwdStatus.size() && itsFwdStatus[payIndex] == ARM_SFRM::ARM_SFRM_AUX_FWD )
				++payIndex;

			if( payIndex >= itsFwdEndTimes->size() )
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					ARM_USERNAME + ": pb with model times! there is no std fwd at the end!" );
			return static_cast<ARM_VectorPtr>(new std::vector<double>(1,(*itsFwdEndTimes)[payIndex] ));
		}
		break;
	case ARM_Numeraire::TerminalEventZc:
		{
			/// Set the very last payment
			double lastEventTime = timeInfos[nbEvents-1]->GetEventTime();
			size_t resetIndex = CC_NS(std,lower_bound)(itsFwdResetTimes->begin(), itsFwdResetTimes->end(), lastEventTime ) - itsFwdResetTimes->begin();

			// test if the forward is auxiliary
			while( resetIndex < itsFwdResetTimes->size() && itsFwdStatus[resetIndex] == ARM_SFRM::ARM_SFRM_AUX_FWD )
				++resetIndex;
			
			if( resetIndex >= itsFwdEndTimes->size() )
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					ARM_USERNAME + ": pb with model times! last reset before last event Time" );
			return static_cast<ARM_VectorPtr>(new std::vector<double>(1,(*itsFwdEndTimes)[resetIndex] ));
		}
		break;
	case ARM_Numeraire::RollingPayment:
		{
			std::vector<double>& paymentDates = new std::vector<double>;
			for (size_t i = 0; i < nbEvents; ++i)
			{
				const std::vector<double>& payTimes = timeInfos[i]->GetPayTimes();
				double lastPayTime	= CC_Max<double>(payTimes[payTimes.size()-1],(*itsFwdEndTimes)[i]);
				size_t payIndex = CC_NS(std,lower_bound)(itsFwdEndTimes->begin(), itsFwdEndTimes->end(), lastPayTime ) - itsFwdEndTimes->begin();

				while( payIndex < itsFwdResetTimes->size() && itsFwdStatus[payIndex] == ARM_SFRM::ARM_SFRM_AUX_FWD )
					++payIndex;

				paymentDates->push_back((*itsFwdEndTimes)[payIndex]);
			}
			return static_cast<ARM_VectorPtr>(paymentDates);
		}
		break;
	case ARM_Numeraire::RollingEvent:
		{
			std::vector<double>& paymentDates = new std::vector<double>;
			for (size_t i = 0; i < nbEvents; ++i)
			{
				double lastEventTime = timeInfos[i]->GetEventTime();
				size_t resetIndex = CC_NS(std,lower_bound)(itsFwdResetTimes->begin(), itsFwdResetTimes->end(), lastEventTime ) - itsFwdResetTimes->begin();

				while ((resetIndex <= itsFwdStatus.size()) && (itsFwdStatus[resetIndex] != ARM_SFRM_STD_FWD))
					++resetIndex;

				paymentDates->push_back((*itsFwdEndTimes)[resetIndex]);
			}
			return static_cast<ARM_VectorPtr>(paymentDates);
		}
		break;
	default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			ARM_USERNAME + ": unknown numeraire type, allowed TerminalZc, TerminalEventZc" );
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: Libor
///	Returns : a vector of libor values
///	Action  : Libor computation
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SFRM::Libor( 
		const string& curveName, 
		double evalTime,
		double fwdStartTime,
		double fwdEndTime,
		double period,
		double fwdResetTime,
		double payTime,
		const ARM_PricingStatesPtr& states) const
{
    ARM_VectorPtr ZcStart= GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdStartTime,states);
    ARM_VectorPtr ZcEnd  = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTime,states);

	/// handle libor with no convexity correction
	if (	(	evalTime <= K_NEW_DOUBLE_TOL 
	 	||	states   == ARM_PricingStatesPtr(NULL) ) 
		&&  (-ARM_GlobalConstant::ARM_SEVENDAYS_LAG <= fwdEndTime - payTime && 
		fwdEndTime - payTime <= ARM_GlobalConstant::ARM_SEVENDAYS_LAG) )
    {
	    ARM_VectorPtr Zct = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTime,states);
		size_t payoffSize = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
		double libor = ((*ZcStart)[0]/(*ZcEnd)[0]-1.0)/period;
		return ARM_VectorPtr( new std::vector<double>(payoffSize,libor) );
    }

    int i,j,nbStates=states->size();
    ARM_VectorPtr values( new std::vector<double>(nbStates) );
    if(-ARM_GlobalConstant::ARM_SEVENDAYS_LAG <= fwdEndTime - payTime &&
		fwdEndTime - payTime <= ARM_GlobalConstant::ARM_SEVENDAYS_LAG)
    {
        /// No convexity (exact if pay = end !)
        for(i=0;i<nbStates;++i)
            (*values)[i]=((*ZcStart)[i]/(*ZcEnd)[i]-1.0)/period;
    }
    else
    {
		double var  = ((ARM_ModelParamsSFRM*) GetModelParams() )->IntegratedVariance(evalTime,fwdResetTime,fwdResetTime);
		double M    = ( (ARM_ModelParamsSFRM*) GetModelParams() )->ShiftValue(fwdResetTime);
        double F,FM;
        if(payTime <= fwdStartTime + ARM_GlobalConstant::ARM_SEVENDAYS_LAG)
        {
            /// Use the in-arrears correction (exact if pay = start !)
            double payConvex = exp(var) - 1.0;
            for(i=0;i<nbStates;++i)
            {
                F = ((*ZcStart)[i]/(*ZcEnd)[i]-1.0)/period;
                FM = F + M;
                (*values)[i] = F + period*FM*FM*payConvex/(1+period*F);
            }
        }
        else
        {
            /// Case start < pay <= end or end < pay
            /// Use the payment lag correction with the additive approximation
            /// It converges to the exact in-arrears correction if pay=start
            /// and obviously to no conversion if pay=end

            ARM_VectorPtr ZcPay = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,payTime,states);

            bool isPayBeforeEnd = (payTime < fwdEndTime);

            // Compute the ccovariance between L(t,Start,End) & L(t,Start,Pay)
		    ARM_VectorPtr coVarSP(new std::vector<double>(nbStates));

			ARM_IRIndex* pIndex	= ((ARM_ModelParamsSFRM*)GetModelParams())->GetIRIndex();
			double indexPeriod = 1.0/pIndex->GetResetFrequency();

			int nbFwds = (payTime-fwdStartTime)/K_YEAR_LEN/indexPeriod;

			double startFwd, endFwd, fwdPeriod, shift, cov, fwdVal;

			ARM_VectorPtr ZcFwdStart;
			ARM_VectorPtr ZcFwdEnd;

			for (i = 0; i <= nbFwds; ++i)
			{
				startFwd = fwdStartTime+i*indexPeriod*K_YEAR_LEN;
				if (i < nbFwds)
					endFwd = startFwd+indexPeriod*K_YEAR_LEN;
				else
					endFwd = payTime;

				if (fabs(endFwd-startFwd) > K_DOUBLE_TOL)
				{
					shift = ( (ARM_ModelParamsSFRM*) GetModelParams() )->ShiftValue(startFwd);

					cov = ((ARM_ModelParamsSFRM*)GetModelParams())->IntegratedCovariance(evalTime,fwdStartTime,fwdStartTime,startFwd);;

					ZcFwdStart=GetDiscountFunctor()->DiscountFactor(curveName,evalTime,startFwd,states);
					ZcFwdEnd=GetDiscountFunctor()->DiscountFactor(curveName,evalTime,endFwd,states);

					fwdPeriod = (endFwd-startFwd)/360;

					for (j = 0; j < nbStates; ++j)
					{
						fwdVal= ((*ZcFwdStart)[j]/(*ZcFwdEnd)[j]-1.0)/fwdPeriod;
						(*coVarSP)[j] += fwdPeriod*(fwdVal+shift)/(1+fwdPeriod*fwdVal)*cov;
					}
				}
			}



            double MSP = M; // assuming L(t,Start,Pay) same shift as L(t,Start,End)
            double deltaSP = period*(payTime - fwdStartTime)/(fwdEndTime-fwdStartTime);

		    double MEP  = M; // assuming L(t,Pay,End) or L(t,End,Pay) same shift as L(t,Start,End)
            double deltaEP,varEP=0.0;
            if(isPayBeforeEnd)
                deltaEP = period*(fwdEndTime-payTime)/(fwdEndTime-fwdStartTime);
            else
            {
				ARM_IRIndex* pIndex	= ((ARM_ModelParamsSFRM*)GetModelParams())->GetIRIndex();
                deltaEP = period*(payTime-fwdEndTime)/(fwdEndTime-fwdStartTime);
                varEP = ((ARM_ModelParamsSFRM*) GetModelParams() )->IntegratedCovariance(evalTime,fwdResetTime,fwdEndTime,payTime,payTime);
            }

            double FSP,FMSP; // L(t,Start,Pay) & shifted libor 
            double FEP,FMEP; // L(t,Pay,End) or L(t,End,Pay) & shifted libor
            double coVar,volCoefEP,adjustEP=0.0;
            for(i=0;i<nbStates;++i)
            {
                F = ((*ZcStart)[i]/(*ZcEnd)[i]-1.0)/period;
                FM = F + M;

                FSP = ((*ZcStart)[i]/(*ZcPay)[i]-1.0)/deltaSP;
                FMSP = FSP + MSP;

                if(isPayBeforeEnd)
                {
                    FEP = ((*ZcPay)[i]/(*ZcEnd)[i]-1.0) / deltaEP;
                    volCoefEP = deltaEP*(FEP+MEP)/(1 + deltaEP*FEP);
                }
                else
                {
                    FEP = ((*ZcEnd)[i]/(*ZcPay)[i]-1.0) / deltaEP;
                    FMEP = FEP + MEP;
                    volCoefEP = - deltaEP * FMEP / (1 + deltaEP*FEP);
                    adjustEP = - volCoefEP * FMEP * (exp(varEP)-1.0);
                    FEP += adjustEP;
                    adjustEP *= F/FMEP;
                }
                FMEP = FEP + MEP;

                coVar = ( period*FM/(1+period*F)*var - (*coVarSP)[i] ) / volCoefEP;

                (*values)[i] = F + volCoefEP * ( FM * (exp(coVar) - 1.0) + adjustEP );
            }
        }
    }

    return values;
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: VanillaCaplet
///	Returns: a vector of Caplet(t,L(R,S),K,S-E)
///	Action : Closed form formula for standard
///          caplet/floorlet (i.e. on libor resetting
///          in advance, paying in arrears)
///				after testing, we found out that the 
///				pricing using ARM_Caplet was inefficient!
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SFRM::VanillaCaplet(
		const string& curveName, 
		double evalTime,
		double payTime, 
		double period,
        double payNotional,
		double fwdResetTime, 
		double fwdStartTime,
        double fwdEndTime,
		double fwdPeriod,
        const std::vector<double>& strikesPerState,
        int capFloor,
		const ARM_PricingStatesPtr& states) const
{
	/// handle the case of dummy states!
	if( states == ARM_PricingStatesPtr(NULL) )
// FIXMEFRED: mig.vc8 (25/05/2007 15:54:21):cast
		return static_cast<ARM_VectorPtr>(new std::vector<double>(1,0.0));
		
	ARM_VectorPtr zcPay,libor;
	size_t i,nbStates;

	if( NULL == itsCurrentArg )
	{
		/// Does computation in place!
		zcPay	= GetDiscountFunctor()->DiscountFactor(curveName,evalTime,payTime,states);
        if(-2 <= fwdEndTime - payTime && fwdEndTime - payTime <= 2)
        {
            /// No convexity : libor = (DFStart/DFEnd - 1.0)/ period
		    libor	= GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdStartTime,states);
		    ARM_VectorPtr zcFwdEnd;

		    if( fwdEndTime == payTime && GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
			    zcFwdEnd = zcPay;
            else
			    zcFwdEnd = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTime,states);

		    for(i=0;i<zcPay->size(); ++i)
			    (*libor)[i] = ((*libor)[i]/(*zcFwdEnd)[i]-1.0)/fwdPeriod;
        }
        else
            /// Compute an adjusted forward rate
            libor = Libor(curveName,evalTime,fwdStartTime,fwdEndTime,fwdPeriod,fwdResetTime,payTime,states);
	}
	else
	{
		size_t modelNb = GetModelNb();
		size_t index = states->GetModelState(0,modelNb);
		libor = ((ARM_VanillaCapArgSFRM&) *itsCurrentArg).itsLibors[index];
		zcPay = ((ARM_VanillaCapArgSFRM&) *itsCurrentArg).itsZCPays[index];
	}
	
	nbStates		= zcPay->size();

#if defined( __GP_STRICT_VALIDATION )
	if( strikesPerState.size() != nbStates )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " strike vector size != states size" );
#endif


	ARM_VectorPtr caplet(new std::vector<double>(nbStates));

	/// do we just need to compute the intrinsic value!
	if( fwdResetTime <= K_NEW_DOUBLE_TOL || evalTime >= fwdResetTime-K_NEW_DOUBLE_TOL )
	{
		for(i=0;i<nbStates;i++)
			(*caplet)[i]=payNotional*period*CC_Max<double>(capFloor*(*libor)[i]-capFloor*strikesPerState[i],0) *(*zcPay)[i];
	}
	else
	{
		const double shift	= ( (ARM_ModelParamsSFRM*) GetModelParams() )->ShiftValue(fwdResetTime);
		const double stdDev	= ((ARM_ModelParamsSFRM*) GetModelParams() )->IntegratedVol(evalTime,fwdResetTime,fwdResetTime);

		for(i=0;i<nbStates;i++)
			(*caplet)[i]=payNotional*period*BlackSholes_Formula((*libor)[i]+shift,stdDev,(*zcPay)[i],strikesPerState[i]+shift,capFloor);
	}

	return caplet;
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: VanillaDigital
///	Returns: a vector of Digital(t,L(R,S),K,S-E)
///	Action : Closed form formula for standard
///          digital caplet/floorlet (i.e. on libor resetting
///          in advance, paying in arrears)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SFRM::VanillaDigital(
		const string& curveName, 
		double evalTime,
		double payTime, 
		double period,
        double payNotional,
		double fwdResetTime, 
		double fwdStartTime,
        double fwdEndTime,
        double fwdPeriod,
        const std::vector<double>& strikesPerState,
        int capFloor,
		const ARM_PricingStatesPtr& states) const
{
	/// handle the case of dummy states!
	if( states == ARM_PricingStatesPtr(NULL) )
// FIXMEFRED: mig.vc8 (30/05/2007 15:55:57):cast
		return static_cast<ARM_VectorPtr>(new std::vector<double>(1,0.0));

	/// test case of negative time!
	fwdResetTime = fwdResetTime < 0.0? 0.0 : fwdResetTime;
	ARM_VectorPtr zcPay,libor;
	double shift = ( (ARM_ModelParamsSFRM*) GetModelParams() )->ShiftValue(fwdResetTime);

	size_t i,nbStates;
	if( NULL == itsCurrentArg )
	{
		/// libor = (DFStart/DFEnd - 1.0)/ period
		/// does computation in place!
		libor	= GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdStartTime,states);
		zcPay	= GetDiscountFunctor()->DiscountFactor(curveName,evalTime,payTime,states);
		ARM_VectorPtr zcFwdEnd;

		if(fwdEndTime == payTime && GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
			zcFwdEnd = zcPay;
        else
			zcFwdEnd = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTime,states);

		for(i=0;i<zcPay->size(); ++i)
			(*libor)[i] = ((*libor)[i]/(*zcFwdEnd)[i]-1.0)/fwdPeriod;
	}
	else
	{
		size_t modelNb = GetModelNb();
		size_t index = states->GetModelState(0,modelNb);
		libor = ((ARM_VanillaDigitalArgSFRM&) *itsCurrentArg).itsLibors[index];
		zcPay = ((ARM_VanillaDigitalArgSFRM&) *itsCurrentArg).itsZCPays[index];
	}
	
	nbStates		= zcPay->size();
#if defined( __GP_STRICT_VALIDATION )
	if( strikesPerState.size() != nbStates )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " strike vector size != states size" );
#endif

	double stdDev	= ((ARM_ModelParamsSFRM*) GetModelParams() )->IntegratedVol(evalTime,fwdResetTime,fwdResetTime);
	ARM_VectorPtr digitalCaplet(new std::vector<double>(nbStates));
	double d2;
	for(i=0;i<nbStates;i++)
	{
		d2 = log(((*libor)[i]+shift)/(strikesPerState[i]+shift))/stdDev-.5*stdDev;
        (*digitalCaplet)[i]	= payNotional*period*(*zcPay)[i]*cdfNormal(capFloor*d2);
	}
	return digitalCaplet;
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: VanillaCorridorlet
///	Returns: a vector of Corridor let(t,L(R,S),K,S-E)
///	Action : Closed form formula for standard
///          Corridor caplet/floorlet (i.e. on libor resetting
///          in advance, paying in arrears)
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SFRM::VanillaCorridorlet(
		const   string& curveName, 
		double  evalTime,
        double  payTime,
        double  resetTime,
        double  startTime,
        double  endTime,
        int     fwdPaymentType, 
        double  fwdPaymentPeriod,
        const std::vector<double>& RefIdxResettimes,
        const std::vector<double>& RefIdxStarttimes,
        const std::vector<double>& RefIdxEndtimes,
        const std::vector<double>& RefFwdPeriods,
        const std::vector<double>& RefIndexWeight,
        double  couponMargin,
        const vector<const std::vector<double>*> downBarrierPerState,
        const vector<const std::vector<double>*> upBarrierPerState,
        double  payNotional,
        int     capFloor,
        const   ARM_PricingStatesPtr& states) const 
{
	/// handle the case of dummy states!
	if( states == ARM_PricingStatesPtr(NULL) )
// FIXMEFRED: mig.vc8 (30/05/2007 15:58:05):cast
		return static_cast<ARM_VectorPtr>(new std::vector<double>(1,0.0));
    size_t nbStates = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
	
	/// FIX FIX
	size_t modelNb = 0;
    size_t index    = states->GetModelState(0,modelNb);

     /// test case of negative time!
	resetTime = resetTime < 0.0? 0.0 : resetTime;

    ARM_VectorPtr corridor(new std::vector<double>(nbStates));
    std::vector<double> ExerciseProbability(nbStates,0.0); 
    ARM_VectorPtr libor,zcPay;

    if(itsCurrentArg )
    {
		libor = ((ARM_VanillaCorridorLegArgSFRM*)itsCurrentArg)->itsLibors[index];
		zcPay = ((ARM_VanillaCorridorLegArgSFRM*)itsCurrentArg)->itsZCPays[index];
    }
    else
	{
        zcPay	= GetFixingFunctor()->DiscountFactor(curveName,evalTime,payTime,states);
        ARM_VectorPtr zcEnd = ( (endTime == payTime) && ( GetFixingFunctor()->IsSameModel( *GetDiscountFunctor() ) ) )? zcPay :
                                                        GetFixingFunctor()->DiscountFactor(curveName,evalTime,endTime,states);    
        if(fwdPaymentType != IDXFIXED)
        {
            libor   = GetFixingFunctor()->DiscountFactor(curveName,evalTime,startTime,states);	
            for(int i=0 ; i < nbStates; ++i)
                (*libor)[i] = ((*libor)[i]/(*zcEnd)[i]-1.0)/fwdPaymentPeriod + couponMargin;
        }
        else
        {
            libor = ARM_VectorPtr( new std::vector<double>(nbStates, couponMargin ) );
        }
    }
    
    double shift	= ((ARM_ModelParamsSFRM*) GetModelParams())->ShiftValue(resetTime);

    double d2;
    size_t i,k;
    for( k = 0; k < RefIdxResettimes.size(); ++k)
    {
        double ref_resetTime = RefIdxResettimes[k];
        double ref_startTime = RefIdxStarttimes[k];
        double ref_endTime   = RefIdxEndtimes[k];

        const double stdDevFwd  = ((ARM_ModelParamsSFRM*) GetModelParams() )->IntegratedVol(evalTime,ref_resetTime,ref_resetTime);
        const double ref_shift	= ( (ARM_ModelParamsSFRM*) GetModelParams() )->ShiftValue(ref_resetTime);

        ARM_VectorPtr ref_Fwd;
        if(itsCurrentArg)
        {
            ref_Fwd = ARM_VectorPtr((std::vector<double>&)((((ARM_VanillaCorridorLegArgSFRM*)itsCurrentArg)->itsRef_Libors[index])[k])->Clone());
            if (fabs(ref_endTime - payTime) > 5.0)
            {
                const ARM_VectorPtr   ref_DfFwd = (((ARM_VanillaCorridorLegArgSFRM*)itsCurrentArg)->itsRef_DfFwds[index])[k];
                const ARM_VectorPtr   ref_zcEnd = (((ARM_VanillaCorridorLegArgSFRM*)itsCurrentArg)->itsRef_ZCEnds[index])[k];

				/// WARNING ASSUMPTION: correlation is set to 100% 
                for(i=0 ; i < nbStates; ++i)
                    (*ref_Fwd)[i] *= (1.0 + (*ref_DfFwd)[i]* stdDevFwd*(ref_shift+(*ref_zcEnd)[i]) * stdDevFwd);
            }
        }
        else
        {
            ref_Fwd           = DiscountFactor(curveName,evalTime,ref_startTime,states);
            const ARM_VectorPtr ref_zcEnd   = DiscountFactor(curveName,evalTime,ref_endTime,states);
            
            for(i=0 ; i < nbStates; ++i)
                (*ref_Fwd)[i] = ((*ref_Fwd)[i]/(*ref_zcEnd)[i]-1.0)/RefFwdPeriods[k];     
        
            /// adjustment between payment date for Foward Libor
            if (fabs(ref_endTime - payTime) > 5.0)
            {
				 //const double stdDevFwd_AtAsof  = ((ARM_ModelParamsSFRM*) GetModelParams() )->IntegratedVol(0.0,ref_resetTime,ref_resetTime);

                if (ref_endTime > payTime)
                {
                    double periodPayEnd = CountYears(KACTUAL_360, payTime, ref_endTime);               
                    for(i=0 ; i < nbStates; ++i)
                    {
                        double DfFwd     = (*ref_zcEnd)[i]/(*zcPay)[i];
                        double ref_zcPay = (1.0/DfFwd -1.0)/periodPayEnd;
                        (*ref_Fwd)[i] *= (1.0 + periodPayEnd * DfFwd* stdDevFwd*(ref_shift+ref_zcPay) * stdDevFwd); 
                    }
                }
                else
                {
                    double periodPayEnd = CountYears(KACTUAL_360, ref_endTime,payTime);
                    for(i=0 ; i < nbStates; ++i)
                    {
                        double DfFwd    = (*ref_zcEnd)[i]/(*zcPay)[i];
                        double ref_zcPay= (DfFwd -1.0)/periodPayEnd;
                        (*ref_Fwd)[i] *= (1.0 - periodPayEnd * DfFwd* stdDevFwd*(ref_shift+ref_zcPay) * stdDevFwd);
                    }
                }
            }
        }
          
        ARM_VectorPtr DownProbability(new std::vector<double>(nbStates,1.0));
        ARM_VectorPtr UpProbability(new std::vector<double>(nbStates,1.0));
        for(i=0 ; i < nbStates; ++i)
        {
            /// Warning: To be coherent with kernel code
            if((*(downBarrierPerState[i]))[k] > 0.0001)
            {
                d2 = log(((*ref_Fwd)[i]+ ref_shift)/( (*(downBarrierPerState[i]))[k] +ref_shift))/stdDevFwd - 0.5*stdDevFwd;
                (*DownProbability)[i] = cdfNormal(d2);
            }

			if(((*(upBarrierPerState[i]))[k] > 0.0001) && ((*(upBarrierPerState[i]))[k] >(*(downBarrierPerState[i]))[k]))
			{
				d2 = log(((*ref_Fwd)[i]+ ref_shift)/( (*(upBarrierPerState[i]))[k] +ref_shift))/stdDevFwd - 0.5*stdDevFwd;
				(*UpProbability)[i] = cdfNormal(d2);
			}

        }
        
        if(fwdPaymentType != IDXFIXED)
        {
            resetTime = resetTime > ref_resetTime ? ref_resetTime : resetTime;
            
			/// ATM convexity correction
            double VolPay   = ((ARM_ModelParamsSFRM*) GetModelParams() )->IntegratedVol(evalTime,MIN(ref_startTime,resetTime),resetTime);
            double VolStart = ((ARM_ModelParamsSFRM*) GetModelParams() )->IntegratedVol(evalTime,MIN(ref_startTime,resetTime),ref_startTime);
            
			/// WARNING ASSUMPTION: correlation is set to 100% 
            double proba=exp(VolPay*VolStart);
            
			/// in the case where the rate is not a CMS, we use a change of probability
            for(i=0 ; i < nbStates; ++i)
            {
                double ref_Fwdadjust = proba*((*ref_Fwd)[i] + ref_shift) - ref_shift;
                
                /// Warning: To be coherent with kernel code
                if((*(downBarrierPerState[i]))[k] > 0.0001)
                {
                    d2 = log((ref_Fwdadjust+ ref_shift)/( (*(downBarrierPerState[i]))[k] +ref_shift))/stdDevFwd - 0.5*stdDevFwd;
                    (*DownProbability)[i] = (cdfNormal(d2) *((*libor)[i] + shift - couponMargin) 
                                             + (*DownProbability)[i] * (couponMargin - shift))/(*libor)[i];
                }
				if(((*(upBarrierPerState[i]))[k] > 0.0001) && ((*(upBarrierPerState[i]))[k] >=(*(downBarrierPerState[i]))[k]))
				{

					d2 = log((ref_Fwdadjust+ ref_shift)/( (*(upBarrierPerState[i]))[k] +ref_shift))/stdDevFwd - 0.5*stdDevFwd;
					(*UpProbability)[i] = (cdfNormal(d2) *((*libor)[i] + shift - couponMargin)  
										   + (*UpProbability)[i] * (couponMargin - shift))/(*libor)[i];
				}
            }
        }
        
        for(i=0 ; i < nbStates; ++i)
            ExerciseProbability[i] += ((*DownProbability)[i] - (*UpProbability)[i])* RefIndexWeight[k];
    }

    for(i=0; i<nbStates; i++)
        (*corridor)[i]	= payNotional*capFloor*(*zcPay)[i]*(*libor)[i]*ExerciseProbability[i];

    return corridor;
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: VanillaCorridorlet Restrikeable
///	Returns: a vector of Corridor let(t,L(R,S),K,S-E)
///	Action : Closed form formula for standard
///          Corridor caplet/floorlet (i.e. on libor resetting
///          in advance, paying in arrears) + numerical integration
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SFRM::VanillaCorridorletRestrikeable(
		const   string& curveName, 
		double  evalTime,
        double  payTime,
        double  resetTime,
        double  startTime,
        double  endTime,
        int     fwdPaymentType, 
        double  fwdPaymentPeriod,
        const std::vector<double>& RefIdxResettimes,
        const std::vector<double>& RefIdxStarttimes,
        const std::vector<double>& RefIdxEndtimes,
        const std::vector<double>& RefFwdPeriods,
        const std::vector<double>& RefIndexWeight,
        double  couponMargin,
        const vector<const std::vector<double>*> downBarrierPerState,
        const vector<const std::vector<double>*> upBarrierPerState,
        double  payNotional,
        int     capFloor,
        const   ARM_PricingStatesPtr& states) const 
{
	/// handle the case of dummy states!
	if( states == ARM_PricingStatesPtr(NULL) )
// FIXMEFRED: mig.vc8 (30/05/2007 15:58:14):cast
	return static_cast<ARM_VectorPtr>(new std::vector<double>(1,0.0));
    size_t nbStates = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
	
	/// FIX FIX
	size_t modelNb = 0;
    size_t index    = states->GetModelState(0,modelNb);

     /// test case of negative time!
	resetTime = resetTime < 0.0? 0.0 : resetTime;

	/////// Start Amine
	//// Law of F at Reset Date Under QTpay
	double m = 1;
	double var = 0.04;
	double stdDev = 0.2;

	//// Numerical integration
	int nbQuadPoints = 100;
	// Numerical integration. We consider the inyervall [m-5stdDev,m+5stdDev]
	GaussLegendre_Coefficients Quadrature(nbQuadPoints,m - 5 * stdDev,m + 5 * stdDev);
	
	double t,w,value_i;
    double result=0.0;
    for(int i=0;i<nbQuadPoints;++i)
    {
        t   = Quadrature.get_point(i);
        w   = Quadrature.get_weight(i);
		value_i = 1.0;
        result += w * value_i;
    }

	ARM_VectorPtr corridor(new std::vector<double>(nbStates,result));
    return corridor;
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: VanillaSwaption
///	Returns: a vector of Swaption(t,S(R,Ti),K)
///	Action : 1 factor like closed form formula for standard
///          swaption (i.e. on standard swap with
///          a "double notional" evaluation of its
///          floating leg)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SFRM::VanillaSwaption(
	const string& curveName,
	double evalTime,
	double swapResetTime,
	const std::vector<double>& fixNotional,
	const std::vector<double>& floatNotional,
	double floatStartTime,
	double floatEndTime,
	const std::vector<double>& floatResetTimes,
	const std::vector<double>& floatStartTimes,
	const std::vector<double>& floatEndTimes,
	const std::vector<double>& floatIntTerms,
	const std::vector<double>& fixPayTimes,
	const std::vector<double>& fixPayPeriods,
    const ARM_GP_Matrix& strikesPerState,
    int callPut,
	const ARM_PricingStatesPtr& states,
	bool isConstantNotional,
	bool isConstantSpread,
	bool isConstantStrike) const
{
	double swapNotional = fixNotional[0];

	/// handle the case of dummy states!
	if( states == ARM_PricingStatesPtr(NULL) )
// FIXMEFRED: mig.vc8 (30/05/2007 15:58:27):cast
		return static_cast<ARM_VectorPtr>(new std::vector<double>(1,0.0));

	if (!(isConstantSpread&&isConstantStrike))
				ARM_THROW( ERR_INVALID_ARGUMENT, "The Model can not price a swaption with variable notional, Spread or Strike!" );

    if( !GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
    {
        /// We need to compute the floating leg by forward method and no more by double notional
        /// but we have not at the moment all floating leg datas => throw an error
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "Swaption pricing not implemented for different discount & fixing curves" );
    }
    if(itsUpdateVolSwapvolFRA)
    {
        /// computes the averageShift
	    ((ARM_VanillaSwaptionArgSFRM*)itsCurrentArg)->SetAverageShift( 
			((ARM_ModelParamsSFRM*) GetModelParams())->AverageShift(*((ARM_VanillaSwaptionArgSFRM*)itsCurrentArg)->GetFloatResetTimes()) );
        /// Computes  Mu
        ((ARM_VanillaSwaptionArgSFRM*)itsCurrentArg)->SetMu( 
			((ARM_ModelParamsSFRM*) GetModelParams())->ComputeVolSwapvolFRA(*(ARM_VanillaSwaptionArgSFRM*)itsCurrentArg,*this,isConstantNotional) );
    }

	ARM_VanillaSwaptionArgSFRM* arg;
	// This flag check if arg is a local variable
	bool isArgLocal = false;
	if( NULL == itsCurrentArg )
	{
		arg = GetVanillaSwaptionArg( curveName, evalTime, swapResetTime,swapNotional,
			fixNotional, floatNotional, floatStartTime, floatEndTime,floatResetTimes,floatStartTimes,floatEndTimes,floatIntTerms,
			fixPayTimes,fixPayPeriods, strikesPerState, callPut, states,isConstantNotional );
		isArgLocal = true;
	}
	else
	{
		arg = (ARM_VanillaSwaptionArgSFRM*) &*itsCurrentArg;
	}
    


	double swapVol = ((ARM_ModelParamsSFRM*)GetModelParams() )->LocalVolatity(evalTime,swapResetTime,*arg,arg->GetMu() );
	swapVol *= sqrt((swapResetTime-evalTime)/K_YEAR_LEN);

	/// get all the data
	ARM_VectorPtr swapfwd	;
	ARM_VectorPtr fixAnnuity;
	if (isConstantNotional)
	{
		swapfwd	= arg->GetSwapFwd();
		fixAnnuity= arg->GetFixAnnuity();
	}
	else
	{
		swapfwd	= arg->GetSwapFwd();
		fixAnnuity= arg->GetFixAnnuityWithNominal();
		swapNotional = 1;
	}

	size_t nbStates			= states->size();
	ARM_VectorPtr swaption(new std::vector<double>(nbStates));
	size_t i;

	if( swapVol <= K_NEW_DOUBLE_TOL )
	{
		for(i=0;i<nbStates;i++)
			(*swaption)[i]=swapNotional*(*fixAnnuity)[i]*
			CC_Max<double>(callPut*((*swapfwd)[i]- strikesPerState(i,0)),0.0);
	}
	else
	{
		for(i=0;i<nbStates;i++)
			(*swaption)[i]=swapNotional*BlackSholes_Formula((*swapfwd)[i]+arg->GetAverageShift(),
				swapVol,(*fixAnnuity)[i],strikesPerState(i,0)+arg->GetAverageShift(),callPut);
	}

	if (isArgLocal)
		delete arg;

	return swaption;    
}


////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: GetVanillaSwaptionArg
///	Returns: ARM_VanillaSwaptionArgSFRM
///	Action : computes all the required data for the pricing of a swaption
///				except the vol
////////////////////////////////////////////////////
ARM_VanillaSwaptionArgSFRM* ARM_SFRM::GetVanillaSwaptionArg( 
	const string& curveName,
	double evalTime,
	double swapResetTime,
	double swapNotional,
    const std::vector<double>& fixNotional,
	const std::vector<double>& floatNotional,
	double floatStartTime,
	double floatEndTime,
	const std::vector<double>& floatResetTimes,
	const std::vector<double>& floatStartTimes,
	const std::vector<double>& floatEndTimes,
	const std::vector<double>& floatIntTerms,
	const std::vector<double>& fixPayTimes,
	const std::vector<double>& fixPayPeriods,
    const ARM_GP_Matrix& strikesPerState,
    int callPut,
	const ARM_PricingStatesPtr& states,
	bool isConstantNotional	) const
{
	double asOfDate	= GetAsOfDate().GetJulian();
	ARM_Date startDate(asOfDate+floatStartTime);
	ARM_Date endDate(asOfDate+floatEndTime);

	CC_NS(std,auto_ptr)<ARM_DateStrip> pFloatDateStrip( GetFloatDateStrip(startDate,endDate) );

	double resetTime = swapResetTime;
	std::vector<double>& pFixPayTimes	 = static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(fixPayTimes).Clone());
	std::vector<double>& pFixPayPeriods	 = static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(fixPayPeriods).Clone());
	std::vector<double>& pFloatResetTimes = static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(floatResetTimes).Clone());
	std::vector<double>& pFloatStartTimes = static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(floatStartTimes).Clone());
	std::vector<double>& pFloatEndTimes	 = static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(floatEndTimes).Clone());
	std::vector<double>& pFloatIntTerms	 = static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(floatIntTerms).Clone());

	std::vector<double>& pfixNotional		 = static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(fixNotional).Clone());
	std::vector<double>& pfloatNotional	 = static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(floatNotional).Clone());

    ///Compute the fixed frequency
    int i;
    double periodmax = (*pFixPayPeriods)[0];
    for(i=1; i<pFixPayTimes->size(); ++i)
        periodmax = CC_Max(periodmax,(*pFixPayPeriods)[i]);

    int fixFreq = ROUND(1.0/periodmax);
	/// period of more than every two months, than it should be one month
	if( fixFreq > 8 )
		fixFreq = 12;

	size_t nbFlow = pFloatEndTimes->size();

	/// computes the averageShift
	double averageShift		 = ((ARM_ModelParamsSFRM*) GetModelParams())->AverageShift(*pFloatResetTimes);
    double startTime		 = (*pFloatStartTimes)[0];
	double endTime			 = (*pFloatEndTimes)[nbFlow-1];


	ARM_Currency* pCcy	= const_cast<ARM_SFRM* const>(this)->GetZeroCurve()->GetCurrencyUnit();
	int floatDayCount	= pCcy->GetLiborIndexDayCount();
	std::vector<double> margin = std::vector<double>(1,0.0);
	bool isDbleNotional =false;
	
	std::vector<double>& fwdStartTimes	= static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(floatStartTimes).Clone());
	std::vector<double>& fwdEndTimes		= static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(floatEndTimes).Clone());
	std::vector<double>& floatPayTimes	= static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(floatEndTimes).Clone());


	std::vector<double> fwdPayPeriods = std::vector<double>(fwdStartTimes->size());

	for(i=0; i<fwdStartTimes->size(); ++i )
	{
		fwdPayPeriods.Elt(i) = CountYearsWithoutException( floatDayCount, ARM_Date((fwdStartTimes->Elt(i)+asOfDate)), ARM_Date((fwdEndTimes->Elt(i)+asOfDate)) );
	}

    
	ARM_VectorPtr fixAnnuity (NULL);
	ARM_VectorPtr swapFwd (NULL); 
	ARM_VectorPtr fixAnnuityWithNominal (NULL);

	if (isConstantNotional)
	{
		fixAnnuity = Annuity(curveName, evalTime, fixPayTimes, fixPayPeriods, states);
		ARM_VectorPtr clonedfixAnnuity = ARM_CountedPtr<std::vector<double>>( (std::vector<double>&) fixAnnuity->Clone() );
		///Compute a stub to add to fixAnnuity
		ARM_Date tmpDate(asOfDate + pFixPayTimes->Elt(0));
		char* CurrencyName = pCcy->GetCcyName();  
		tmpDate.AddPeriod(-fixFreq,CurrencyName);
		if(fabs(tmpDate-startDate)>FRMVOL_LAG_THRESHOLD)
		{  
			int fixDayCount = pCcy->GetFixedDayCount();
			double delta = CountYears(fixDayCount,tmpDate,startDate);
			ARM_VectorPtr ZcFwdStart    = GetFixingFunctor()->DiscountFactor(curveName,evalTime,floatStartTime,states);
			ARM_VectorPtr ZcFwdPay    = GetFixingFunctor()->DiscountFactor(curveName,evalTime,pFixPayTimes->Elt(0),states);
			
			for(int i=0; i<states->size(); ++i)
			{
				(*clonedfixAnnuity)[i] = (*clonedfixAnnuity)[i] + delta*((*ZcFwdPay)[i]-(*ZcFwdStart)[i]);
			}
		} 
		
		swapFwd = SwapRateInPlaceWithComputedAnnuity(	curveName, evalTime,
			startTime, endTime, *pFixPayTimes, *pFixPayPeriods, *fwdStartTimes,
			*fwdEndTimes, fwdPayPeriods, *floatPayTimes, *pFloatIntTerms,
			margin, isDbleNotional, clonedfixAnnuity, //// make sure you clone it if neeeded to keep the value of the annuity
			states);
	}
	else
	{
		fixAnnuityWithNominal = AnnuityWithNominal(curveName, evalTime, fixPayTimes, fixPayPeriods, fixNotional,states)		;
		
		ARM_VectorPtr clonedfixAnnuity = ARM_CountedPtr<std::vector<double>>( (std::vector<double>&) fixAnnuityWithNominal->Clone() );
		///Compute a stub to add to fixAnnuity
		ARM_Date tmpDate(asOfDate + pFixPayTimes->Elt(0));
		char* CurrencyName = pCcy->GetCcyName();  
		tmpDate.AddPeriod(-fixFreq,CurrencyName);
		double Notional0 = fixNotional[0];
		if(fabs(tmpDate-startDate)>FRMVOL_LAG_THRESHOLD)
		{  
			int fixDayCount = pCcy->GetFixedDayCount();
			double delta = CountYears(fixDayCount,tmpDate,startDate);
			ARM_VectorPtr ZcFwdStart    = GetFixingFunctor()->DiscountFactor(curveName,evalTime,floatStartTime,states);
			ARM_VectorPtr ZcFwdPay    = GetFixingFunctor()->DiscountFactor(curveName,evalTime,pFixPayTimes->Elt(0),states);
			
			for(int i=0; i<states->size(); ++i)
			{
				(*clonedfixAnnuity)[i] = (*clonedfixAnnuity)[i] + delta*Notional0*((*ZcFwdPay)[i]-(*ZcFwdStart)[i]);
			}
		} 

		swapFwd = SwapRateInPlaceWithComputedAnnuityAndNominal(	curveName, evalTime,
			startTime, endTime, *pFixPayTimes, *pFixPayPeriods, *fwdStartTimes,
			*fwdEndTimes, fwdPayPeriods, *floatPayTimes, *pFloatIntTerms,
			margin, clonedfixAnnuity, //// make sure you clone it if neeeded to keep the value of the annuity
			floatNotional,states);

	}
    /// creates the swaption arg
	ARM_VectorPtr nullstrike(new std::vector<double>(1));
	ARM_VanillaSwaptionArgSFRM* arg = new ARM_VanillaSwaptionArgSFRM(
		resetTime, startTime, endTime, averageShift, pFixPayTimes, pFixPayPeriods, pFloatResetTimes, 
		pFloatStartTimes, pFloatEndTimes, pFloatIntTerms, ARM_VectorPtr(NULL), fixAnnuity, swapFwd,fixFreq, fixAnnuityWithNominal,swapFwd,pfixNotional,pfloatNotional);	

	/// computes the mu
	ARM_VectorPtr mu = ((ARM_ModelParamsSFRM*) GetModelParams())->ComputeVolSwapvolFRA(*arg,*this,isConstantNotional);

	/// AAAAAAAAAAAARRRRRRRRRGGGGGGGGG the ComputeVolSwapvolFRA change the swap fwd!
	arg->SetMu( mu );
	arg->SetSwapFwd( swapFwd );

	return arg;
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionIR
///	Routine: VanillaSumOption
///	Returns: a vector of sum option values
///	Action : Closed form formula for standard
///          Sum Option cap/floor 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SFRM::VanillaSumOption(
		const string& curveName,
		double evalTime,
		int capFloor,
		const std::vector<double>& coeffs,
		const std::vector<double>& fwdResetTimes,
		const std::vector<double>& fwdStartTimes,
		const std::vector<double>& fwdEndTimes,
		double payTime,
		const std::vector<double>& fwdPeriods,
		const std::vector<double>& strikesPerState,
		double volatilityRatio,
		double* sumFwd,
		double* sumVol,
		const ARM_PricingStatesPtr& states) const
{
	size_t statesSize = states->size();

	size_t nbPeriods = coeffs.size();

#if defined( __GP_STRICT_VALIDATION )
	if( fwdResetTimes.size() != nbPeriods )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " fwd reset times size != coeffs size" );

	if( fwdStartTimes.size() != nbPeriods )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " fwd start times size != coeffs size" );

	if( fwdEndTimes.size() != nbPeriods )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " fwd end times size != coeffs size" );

	if( coeffs.size() != nbPeriods )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " fwd end times size != coeffs size" );

	if( strikesPerState.size() != statesSize )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " strike per states size != states size" );

	if (evalTime > fwdResetTimes(0))
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " eval time > first fwd reset time." );

#endif

	size_t i,j,k;

	ARM_Currency* pCcy	= const_cast<ARM_SFRM* const>(this)->GetZeroCurve()->GetCurrencyUnit();
	int floatDayCount = pCcy->GetLiborIndexDayCount();

	ARM_VectorPtrVector libors(nbPeriods);

	ARM_Vector SumLiborWithShift(statesSize), SumLibor(statesSize), shifts(nbPeriods);

	for (i = 0; i < nbPeriods; ++i)
	{
		shifts[i] = ((ARM_ModelParamsSFRM*) GetModelParams() )->ShiftValue(fwdResetTimes[i]);
		libors[i] = Libor(curveName,evalTime,fwdStartTimes[i],fwdEndTimes[i],fwdPeriods[i],fwdResetTimes[i],payTime,states);

		for (j = 0; j < statesSize; ++j)
		{
			SumLibor[j] += (*(libors[i]))[j]*coeffs[i];
			SumLiborWithShift[j] += ((*(libors[i]))[j]+shifts[i])*coeffs[i];
			(*(libors[i]))[j]+=shifts[i];
			(*(libors[i]))[j]*=coeffs[i];
		}
	}

	double stdDevAverage;

	ARM_VectorPtr ZcPay = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,payTime,states);

	ARM_VectorPtr SumOption(new std::vector<double>(statesSize));

	const double Beta = ( (ARM_ModelParamsSFRM*) GetModelParams() )->BetaValue(fwdResetTimes[nbPeriods-1]);

	for (i = 0; i < statesSize; ++i)
	{
		stdDevAverage = 0;

		if ((fabs(SumLibor[i]) > K_DOUBLE_TOL) && (fabs(strikesPerState[i]) > K_DOUBLE_TOL))
		{
			for (j = 0; j < nbPeriods; ++j)
			{
				stdDevAverage += (*(libors[j]))[i]*(*(libors[j]))[i]*((ARM_ModelParamsSFRM*)GetModelParams())->IntegratedVariance(evalTime,fwdResetTimes[j],fwdResetTimes[j]);
				for (k = j+1; k < nbPeriods; ++k)
				{
					stdDevAverage += 2*(*(libors[j]))[i]*(*(libors[k]))[i]*((ARM_ModelParamsSFRM*)GetModelParams())->IntegratedCovariance(evalTime,fwdResetTimes[j],fwdResetTimes[j],fwdResetTimes[k]);
				}
			}

			if (fabs(Beta) > K_DOUBLE_TOL)
			{
				double shift = SumLibor[i]*(1.0/Beta-1.0);

				stdDevAverage /= (SumLibor[i]+shift);
				stdDevAverage /= (SumLibor[i]+shift);
				
				stdDevAverage=sqrt(stdDevAverage);

				(*SumOption)[i] = BlackSholes_Formula(SumLibor[i]/Beta,volatilityRatio*stdDevAverage,(*ZcPay)[i],strikesPerState[i] + SumLibor[i]*(1.0/Beta-1.0),capFloor);
				if(sumFwd)
					*sumFwd = SumLibor[i];
				if(sumVol)
					*sumVol=volatilityRatio*stdDevAverage/sqrt(fwdResetTimes[nbPeriods-1]/K_YEAR_LEN);
			}
			else
			{
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " Beta is null!" );
			}
		}
		else
		{
			(*SumOption)[i] = CC_Max((SumLibor[i]-strikesPerState[i])*capFloor,0.0)*(*ZcPay)[i];
		}
	}

	return SumOption;
}



////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: Clone
///	Returns: 
///	Action : clones the object (virtual constructor)
////////////////////////////////////////////////////
ARM_Object* ARM_SFRM::Clone() const
{ return new ARM_SFRM(*this); }



////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: DiscountFactor Version Richard Modifiee
///	Returns: a vector of Zc(t,T)
///	Action : Closed form formula for DF.
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SFRM::DiscountFactor( const string& curveName, double evalTime,
	double maturityTime, const ARM_PricingStatesPtr& states ) const
{	
	if (!GetNumMethod().IsNull())
	{
		if ((states != ARM_PricingStatesPtr(NULL)) && (states->size() > 1) && (GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDBCKWDLOOKING) && (GetNumMethod()->GetPricingDirCurrLoop() == ARM_NumMethod::GP_BCKWDLOOKING))
		{
			throw Exception(__LINE__,__FILE__,ERR_INVALID_ARGUMENT, ARM_USERNAME + ": The DF function should never be called in the backward loop of the AMC.");
		}
	}

	//ResetDFMap();
	ARM_DFMap::iterator found = itsDFMap.find(maturityTime);
    if(found == itsDFMap.end())
        /// Insert the new function in the map
	{
		/// Waiting for the access to the yield curve with curveName
		ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
		double zc;
		
		/// other simple case with no stochasticity
		if(	evalTime<=K_NEW_DOUBLE_TOL )
		{
			zc = ZcCurve->DiscountPrice(maturityTime/K_YEAR_LEN);
			size_t payoffSize = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
			return ARM_VectorPtr( new std::vector<double>(payoffSize,zc) );
			
		}
		else if( states == ARM_PricingStatesPtr(NULL) )
		{
			zc = ZcCurve->DiscountPrice(maturityTime/K_YEAR_LEN);
			double zct = ZcCurve->DiscountPrice(evalTime/K_YEAR_LEN);
			return ARM_VectorPtr( new std::vector<double>(1,zc/zct) );
		}
		else
		{
			size_t nbStates = states->size();
			
			if ((nbStates == 0) ||
				(itsFwdResetTimes == ARM_VectorPtr(NULL)) ||
				(itsFwdStartTimes == ARM_VectorPtr(NULL)) ||
				(itsFwdEndTimes == ARM_VectorPtr(NULL)) ||
				(itsFwdEndTimes == ARM_VectorPtr(NULL)) ||
				(itsFwdValues == ARM_VectorPtr(NULL)) ||
				(itsFwdIntTerms == ARM_VectorPtr(NULL)) ||
				(itsFwdShifts == ARM_VectorPtr(NULL)))
			{
				return ARM_VectorPtr( new std::vector<double>(0) );
			}
			
#if defined(__GP_STRICT_VALIDATION)
			int FwdNbs	= itsFwdStartTimes->size();
			if(maturityTime > itsFwdEndTimes->Elt(FwdNbs-1) + K_NEW_DOUBLE_TOL)
				throw Exception(__LINE__,__FILE__,ERR_INVALID_ARGUMENT, ARM_USERNAME + ": payment date > terminal date!");
#endif
			
			/// very trivial case?
			if( maturityTime-evalTime < K_NEW_DOUBLE_TOL )
				return ARM_VectorPtr( new std::vector<double>(nbStates,1.0) );
			
			
			//// General Case
			int fwdMinIndex = states->GetPricingStatesContext( GetModelRank() )->ToSFRMPricingStatesContext()->GetFwdMinIndex();			
			int maxIdx;
			
			// If Maturity Time is before the first alive fwd start date
			// we use the first alive fwd to do interpolation
			if (maturityTime < (*itsFwdEndTimes)[fwdMinIndex])
				maxIdx = fwdMinIndex;
			else
			{
				/// find the index such as (*itsEndStartTimes)[maxIdx] < maturityTime 
				maxIdx = upper_boundPosWithPrecision(*itsFwdEndTimes,maturityTime );
			}
			// Look for the standard forward which contains maturity date
			int lastStdFwd = maxIdx;
			
			while (itsFwdStatus[lastStdFwd] == ARM_SFRM::ARM_SFRM_AUX_FWD)
				++lastStdFwd;
			
			ARM_VectorPtr values(new std::vector<double>(nbStates));
			ARM_Currency* pCcy	= const_cast<ARM_SFRM* const>(this)->GetZeroCurve()->GetCurrencyUnit();
			int floatDayCount	= pCcy->GetLiborIndexDayCount();
			
			double dfEvalTime = GetZeroCurve()->DiscountPrice(evalTime/K_YEAR_LEN);
			double dfMaturityTime = GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN);
			size_t i;
			size_t modelNb = GetModelNb();
			double zc;
			
			/// computes basically the product of 1+delta(i)*Lib(i) from evaltime to maturityTime
			
			
#if defined(__GP_STRICT_VALIDATION)
			if(fwdMinIndex > maxIdx)
				throw Exception(__LINE__,__FILE__,ERR_INVALID_ARGUMENT, ARM_USERNAME + ": All the fwds are dead! Please advise.");
#endif
			int firstStdFwd = 0, nextStdFwdAfterLast = 0;		
			
			// Look for the first alive standard forward
			firstStdFwd = fwdMinIndex;
			
			while (itsFwdStatus[firstStdFwd] == ARM_SFRM::ARM_SFRM_AUX_FWD)
				++firstStdFwd;
			
			// Compute the standard forward just after maxIdx (it will be used for 
			// calibration )
			if (lastStdFwd < itsFwdStartTimes->size()-1)
			{
				nextStdFwdAfterLast = lastStdFwd+1;
				
				while (itsFwdStatus[nextStdFwdAfterLast] == ARM_SFRM::ARM_SFRM_AUX_FWD)
					++nextStdFwdAfterLast;
				
				if (itsFwdStatus[nextStdFwdAfterLast] == ARM_SFRM::ARM_SFRM_AUX_FWD)
					nextStdFwdAfterLast=lastStdFwd;
			}
			else
			{
				nextStdFwdAfterLast = lastStdFwd;
			}
			
			/// We first compute the first standard fwd after the eval date and the last
			// standard fwd before the maturity date
			double firstStartTime	= (*itsFwdStartTimes)[firstStdFwd];
			double firstEndTime		= (*itsFwdEndTimes)[firstStdFwd];
			double lastStartTime	= (*itsFwdStartTimes)[lastStdFwd];
			double lastEndTime		= (*itsFwdEndTimes)[lastStdFwd];
			
			if ((firstStartTime-maturityTime)>K_NEW_DOUBLE_TOL)
			{
				firstStartTime= maturityTime;
				lastEndTime = firstStartTime;
				lastStdFwd=firstStdFwd-1;
			}
				
			/// computes the start stub
			double fwdStartStub, deltaStartStub, shiftStartStub, fwdStartStubWithAlea;
			
			/// do we have a startStub?
			if( evalTime<firstStartTime)
			{
				deltaStartStub	= CountYears(floatDayCount,evalTime,firstStartTime);
				fwdStartStub	= (dfEvalTime/ZcCurve->DiscountPrice(firstStartTime/K_YEAR_LEN)-1.0)/deltaStartStub;
				shiftStartStub	= (*itsFwdShifts)[fwdMinIndex];
			}
			
			double deltai,fwd;
			//// linear interpolation
			
			if(maxIdx==fwdMinIndex && ((*itsFwdEndTimes)[fwdMinIndex]-maturityTime)>-FRMVOL_LAG_THRESHOLD)
			{
				
				for(i=0; i<nbStates; ++i )
				{
					zc = 1.0;
					double firsttime = (*itsFwdStartTimes)[fwdMinIndex];
					if( evalTime < firsttime )
					{
						double deltaStub	= CountYears(floatDayCount,evalTime,firsttime);
						double fwdStub	= (dfEvalTime/ZcCurve->DiscountPrice(firsttime/K_YEAR_LEN)-1.0)/deltaStub;
						double shiftStub	= (*itsFwdShifts)[fwdMinIndex];
						
						fwdStartStubWithAlea = (states->GetModelState(i,modelNb+fwdMinIndex)+shiftStub)
							/((*itsFwdValues)[fwdMinIndex]+shiftStub);
						zc	/= 1.0+deltaStub*((fwdStub+shiftStartStub)*fwdStartStubWithAlea-shiftStub);
					}

					deltai	= itsFwdIntTerms->Elt(fwdMinIndex);
					fwd		= states->GetModelState(i,modelNb+fwdMinIndex);
					zc	   /= 1.0+deltai*fwd;

					lastEndTime = (*itsFwdEndTimes)[fwdMinIndex];
					lastStdFwd=fwdMinIndex;
					(*values)[i] = zc;
				}
			}
			
			else
			{
				///// Look at the Map for last standard forward		
				ARM_DFMap::iterator found2 = itsDFMap.find(lastEndTime);
				if(found2 == itsDFMap.end())
				{	
					
					for(i=0; i<nbStates; ++i )
					{
						zc = 1.0;
						/// We compute the stub start using the stochastic part of the first auxiliary or a standard forward.
						if( evalTime < firstStartTime )
						{
							fwdStartStubWithAlea = (states->GetModelState(i,modelNb+fwdMinIndex)+shiftStartStub)
								/((*itsFwdValues)[fwdMinIndex]+shiftStartStub);
							zc	/= 1.0+deltaStartStub*((fwdStartStub+shiftStartStub)*fwdStartStubWithAlea-shiftStartStub);
						}
						
						/// We first use the simulated fwd except the last one which can be the end stub
						for( int j=firstStdFwd; j<=lastStdFwd; ++j)
						{
							if (j == firstStdFwd || itsFwdStatus[j] != ARM_SFRM::ARM_SFRM_AUX_FWD)
							{
								deltai	= itsFwdIntTerms->Elt(j);
								fwd		= states->GetModelState(i,modelNb+j);
								zc	   /= 1.0+deltai*fwd;
							}
						}    				
						(*values)[i] = zc;
					}
					/////// Store in the map the last forward
					ARM_VectorPtr values2(new std::vector<double>(*values));
					pair< double,ARM_VectorPtr > value2(lastEndTime,values2);
					pair< ARM_DFMapIter,bool > result = itsDFMap.insert(value2);
					if(!result.second)
						ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": can't insert a DF Vector in the map");								
					
				}
				else
					values = ARM_VectorPtr(new std::vector<double>(*(found2->second)));
				
			}			
			double fwdEndStub	= 0.0, deltaEndStub = 0.0, shiftEndStub = 0.0, fwdEndStubWithAlea;			
			
			/// do we have a endStub?
			if( maturityTime < lastEndTime )
			{
				deltaEndStub	= CountYears(floatDayCount,maturityTime,lastEndTime);
				fwdEndStub		= (dfMaturityTime/ZcCurve->DiscountPrice(lastEndTime/K_YEAR_LEN)-1.0)/deltaEndStub;
				shiftEndStub	= (*itsFwdShifts)[maxIdx];
				
				for(i=0; i<nbStates; ++i )
				{
					
					deltaEndStub	= CountYears(floatDayCount,maturityTime,lastEndTime);
					fwdEndStub		= (dfMaturityTime/ZcCurve->DiscountPrice(lastEndTime/K_YEAR_LEN)-1.0)/deltaEndStub;
					shiftEndStub	= (*itsFwdShifts)[lastStdFwd];
					
					fwdEndStubWithAlea = (states->GetModelState(i,modelNb+lastStdFwd)+shiftEndStub)
						/((*itsFwdValues)[lastStdFwd]+shiftEndStub);
					(*values)[i]	*= 1.0+deltaEndStub*((fwdEndStub+shiftEndStub)*fwdEndStubWithAlea-shiftEndStub);
				}
			}
			found = itsDFMap.find(maturityTime);
			if(found == itsDFMap.end())	
			{
				pair< double,ARM_VectorPtr > value(maturityTime,values);
				pair< ARM_DFMapIter,bool > result = itsDFMap.insert(value);
				if(!result.second)
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": can't insert a DF Vector in the map");
				return ARM_VectorPtr(new std::vector<double>(*(result.first->second)));
			}
			else
				return ARM_VectorPtr(new std::vector<double>(*(found->second)));
		}
	}
	else //// the maturity is found in the map
		return ARM_VectorPtr(new std::vector<double>(*(found->second)));
}


////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: LocalDrifts
///	Returns: A vector saving the local drift of
///          the state variable
///	Action : Compute local relative drifts of the
///          state variable between each time step
////////////////////////////////////////////////////
void ARM_SFRM::IntegratedLocalDrifts(
	const std::vector<double>& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts ) const
{
	/// nothing as we need the full path to compute localDrift
	/// and compute drift on the fly
	size_t factorsNb=FactorCount();
	size_t nbSteps	= timeSteps.size();
	size_t modelNb	= GetModelNb();

    /// Initialise a Relative Drift matrix to 1 (no drift)
    relativeDrifts	= ARM_GP_MatrixPtr( new ARM_GP_Matrix( (nbSteps-1)*(modelNb+1), factorsNb, 1.0 ) );
    absoluteDrifts	= ARM_GP_MatrixPtr( NULL);
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: PostInit
///	Returns: nothing
///	Action : Function to init after the numerical method
///			 has been initialised.. non const
////////////////////////////////////////////////////
void ARM_SFRM::PostInit()
{
	/// however initialise the initial forward
	size_t totalFwds=itsFwdStartTimes->size();
	itsFwdValues=ARM_VectorPtr(new std::vector<double>(totalFwds));
	itsFwdShifts=ARM_VectorPtr(new std::vector<double>(totalFwds));
	double dfStart,dfEnd;

	for(size_t i=0; i<totalFwds; ++i)
    {	
		dfStart	  = GetZeroCurve()->DiscountPrice((*itsFwdStartTimes)[i]/K_YEAR_LEN);
		dfEnd	  = GetZeroCurve()->DiscountPrice((*itsFwdEndTimes)[i]/K_YEAR_LEN);
		(*itsFwdValues)[i]= (dfStart/dfEnd-1.0)/(*itsFwdIntTerms)[i];
		(*itsFwdShifts)[i]= ((ARM_ModelParamsSFRM*) GetModelParams())->ShiftValue((*itsFwdResetTimes)[i]);
	}

	ComputeNumeraireTimeIndex();
	if(GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_BCKWDLOOKING )
	{
		ComputeFwdsCorrelVarCoeffs();
		ComputeFwdDrifts();
	}

	if ((GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDLOOKING ) ||
		(GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDBCKWDLOOKING ))
		SigneItsModelStateLocalStdDevs();		
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: NumMethodStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local and global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_SFRM::NumMethodStateLocalVariances(
    const std::vector<double>& timeSteps,
    ARM_MatrixVector& localVariances ) const
{
	size_t totalFwds=itsFwdStartTimes->size(),
		factorsNb=FactorCount();
	size_t nbSteps	= timeSteps.size();
	size_t modelNb	= GetModelNb();
    double step		= timeSteps[0],
		nextStep;
	size_t offsetIndex	= (nbSteps-1)*modelNb;

#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "localDrifts.size() != offsetIndex" );
#endif
	localVariances.resize((nbSteps-1)*(modelNb+1));
	size_t startTimePos=0;

#if defined(__GP_STRICT_VALIDATION)
	if( step != 0.0 )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": first time step != 0" );
#endif

	for(size_t i=0;i<nbSteps-1 ;++i)
	{
		/// initalize the toTime
		nextStep=timeSteps[i+1];

		/// initialize everything
		localVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(factorsNb);

		/// get the first bigger time
		ARM_VectorPtr variancePerFactor = ((ARM_ModelParamsSFRM*)GetModelParams())->IntegratedTreeStatesVariancePerFactor(step,nextStep);
		for(size_t k=0;k<factorsNb;++k)
			(*localVariances[offsetIndex+i])(k,k) =(*variancePerFactor)[k];
		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: NumMethodStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local and global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_SFRM::ModelStateLocalVariances(
    const std::vector<double>& timeSteps,
    ARM_MatrixVector& localVariances ) const
{
	size_t totalFwds=itsFwdStartTimes->size(),
		factorsNb=FactorCount();
	/// resize localVariances and std deviation
	size_t timeStepsSize = timeSteps.size();
	size_t modelNb		= GetModelNb();
	size_t offsetIndex	= (timeStepsSize-1)*modelNb;
	size_t startTimePos=0;
	double	fromTime= timeSteps[0],toTime;

#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "localDrifts.size() != offsetIndex" );
#endif
    localVariances.resize((timeStepsSize-1)*(modelNb+1));
	size_t i;


#if defined(__GP_STRICT_VALIDATION)
	if( fromTime != 0.0 )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": first time step != 0" );
#endif

	for(i=0;i<timeStepsSize-1 ;++i)
	{
		/// initalize the toTime
		toTime  = timeSteps[i+1];

		/// initialize everything
		localVariances[i] = new ARM_GP_Matrix(totalFwds,factorsNb);

		/// get the first bigger time
		while(startTimePos< itsFwdResetTimes->size()
			&& (*itsFwdResetTimes)[startTimePos] < toTime )
			++startTimePos;

		for(size_t j=startTimePos; j<totalFwds; ++j)
		{
			ARM_VectorPtr variancePerFactor = ((ARM_ModelParamsSFRM*)GetModelParams())->FwdsVarCorrelCoeffs(
				(*itsFwdResetTimes)[j],
				(*itsFwdResetTimes)[j]);

			for(size_t k=0;k<factorsNb;++k)
				(*localVariances[i])(j,k) =(*variancePerFactor)[k];
		}
		fromTime= toTime;
	}	
}



////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: ComputeFwdsCorrelVarCoeffs
///	Returns: A vector of 1D matrix
///	Action : Compute in the Fwds variances  per factor 
///          the coefficients depending on the FWD reset dates 
////////////////////////////////////////////////////
void ARM_SFRM::ComputeFwdsCorrelVarCoeffs( void ) const
{
	size_t totalFwds=itsFwdStartTimes->size();
	size_t coeffSize=const_cast<ARM_SFRM *>(this)->itsFwdsCorrelVarCoeffs.size();
	for(size_t i=0; i<coeffSize; ++i)
	{
		delete (const_cast<ARM_SFRM *>(this)->itsFwdsCorrelVarCoeffs[i]);
		const_cast<ARM_SFRM *>(this)->itsFwdsCorrelVarCoeffs[i] = NULL;
	}

	size_t factorsNb=FactorCount();
	const_cast<ARM_SFRM *>(this)->itsFwdsCorrelVarCoeffs.resize(factorsNb); 
	int signeCorrel =1;

	for(size_t j=0; j<factorsNb; ++j)
	{
		const_cast<ARM_SFRM *>(this)->itsFwdsCorrelVarCoeffs[j] = new ARM_GP_TriangularMatrix(totalFwds,0);
	}

	ARM_GP_TriangularMatrix coeffs(totalFwds);
	std::vector<double> tmpFwdsCoeffs(totalFwds);

	for(size_t k=0; k<totalFwds; ++k)
	{
		for(size_t l=0; l<=k; ++l)
		{
			ARM_VectorPtr coeffperfactor = ((ARM_ModelParamsSFRM*)GetModelParams())->FwdsVarCorrelCoeffs(
				(*itsFwdResetTimes)[k],(*itsFwdResetTimes)[l]);

			for(size_t j=0;j<factorsNb;++j)
				(*itsFwdsCorrelVarCoeffs[j])(k,l) =(*coeffperfactor)[j];
		}

		CC_NS(std,auto_ptr)<std::vector<double>> correl_k(((ARM_ModelParamsSFRM*)GetModelParams())->InterpolateCorrelation((*itsFwdResetTimes)[k]));

		for(size_t j=0; j<factorsNb; ++j)
		{
			tmpFwdsCoeffs.Elt(k) += (*itsFwdsCorrelVarCoeffs[j])(k,k);
			
			signeCorrel = ( ((*correl_k)[j] > 0) ? 1 : (((*correl_k)[j] < 0) ? -1 : 0) );

			(*itsFwdsCorrelVarCoeffs[j])(k,k) = sqrt((*itsFwdsCorrelVarCoeffs[j])(k,k))*signeCorrel;
		}		
	}

	const_cast<ARM_SFRM *>(this)->itsFwdsVarCoeffs  = ARM_VectorPtr(new std::vector<double>(tmpFwdsCoeffs));
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: ComputeFwdDrifts
///	Returns: 
///	Action : Compute a deterministic drift
/// for each forward
////////////////////////////////////////////////////
void ARM_SFRM::ComputeFwdDrifts(void)
{
	ARM_Currency* pCcy	= const_cast<ARM_SFRM* const>(this)->GetZeroCurve()->GetCurrencyUnit();
	int floatDayCount	= pCcy->GetLiborIndexDayCount();

	size_t nbFwds=itsFwdValues->size();
	size_t factorsNb = FactorCount();

	itsFwdDrifts = ARM_VectorPtr(new std::vector<double>(nbFwds));

	double drift, driftInc, isFirstDriftInc, firstDriftInc, firstDriftIdx, deltaDrift;

	// Forward Before the numeraire forward
	for (int i = itsNumeraireTimeIndex-1; i >= 0; --i)
	{
		drift = 0.0;
		isFirstDriftInc = true;

		if (i < itsNumeraireTimeIndex-1)
		{
			for (size_t j = i+1; j <= itsNumeraireTimeIndex-1; ++j)
			{
				if (itsFwdStatus[j] == ARM_SFRM_STD_FWD)
				{
					driftInc = 0.0;
					for (size_t k = 0; k < factorsNb; ++k)
						driftInc -= (*itsFwdIntTerms)[j]*((*itsFwdValues)[j]+(*itsFwdShifts)[j])/(1.0+(*itsFwdIntTerms)[j]*(*itsFwdValues)[j])*(*itsFwdsCorrelVarCoeffs[k])(j,i);
					if (isFirstDriftInc)
					{
						firstDriftInc = driftInc;
						firstDriftIdx = j;
						isFirstDriftInc = false;
					}
					drift += driftInc;
				}
			}
		}
		(*itsFwdDrifts)[i] = drift;
		if (itsFwdStatus[i] == ARM_SFRM_AUX_FWD)
		{
			deltaDrift = CountYears(floatDayCount,(*itsFwdStartTimes)[firstDriftIdx],(*itsFwdEndTimes)[i])
						/CountYears(floatDayCount,(*itsFwdStartTimes)[firstDriftIdx],(*itsFwdEndTimes)[firstDriftIdx]);
			(*itsFwdDrifts)[i] = deltaDrift*firstDriftInc;

		}

	}

	// Forward after the numeraire forward
	if (nbFwds > itsNumeraireTimeIndex)
	{
		for (int i = itsNumeraireTimeIndex; i < nbFwds; ++i)
		{
			drift = 0.0;
			isFirstDriftInc = true;
			for (int j = i; j <= itsNumeraireTimeIndex-1; --j)
			{
				if (itsFwdStatus[j] == ARM_SFRM_STD_FWD)
				{
					driftInc = 0.0;
					for (size_t k = 0; k <= factorsNb; ++k)
						driftInc += (*itsFwdIntTerms)[j]*((*itsFwdValues)[j]+(*itsFwdShifts)[j])/(1.0+(*itsFwdIntTerms)[j]*(*itsFwdValues)[j])*(*itsFwdsCorrelVarCoeffs[k])(i,j);
					if (isFirstDriftInc)
					{
						firstDriftInc = driftInc;
						firstDriftIdx = j;
						isFirstDriftInc = false;
					}
					drift += driftInc;
				}
			}
			(*itsFwdDrifts)[i] = drift;
			if (itsFwdStatus[i] == ARM_SFRM_AUX_FWD)
			{
				deltaDrift = CountYears(floatDayCount,(*itsFwdStartTimes)[i],(*itsFwdEndTimes)[firstDriftIdx])
							/CountYears(floatDayCount,(*itsFwdStartTimes)[firstDriftIdx],(*itsFwdEndTimes)[firstDriftIdx]);
				(*itsFwdDrifts)[i] += deltaDrift*firstDriftInc;

			}
		}
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: ComputeFwdsCorrelVarCoeffs
///	Returns: A vector of 1D matrix
///	Action : Compute in the Fwds variances  per factor 
///          the coefficients depending on the FWD reset dates 
////////////////////////////////////////////////////
void ARM_SFRM::SigneItsModelStateLocalStdDevs( void ) const
{
	size_t totalFwds=itsFwdStartTimes->size();
	size_t factorsNb=FactorCount();
	const ARM_MatrixVector& localStdDev = GetModelStateLocalStdDevs();

	int timeStepSize = (localStdDev).size();

	size_t timeIndex = 0;
	std::vector<double>& correl_k;
	while( timeIndex<timeStepSize)
	{
		for(size_t k=0; k<totalFwds; ++k )
		{
			correl_k = ((ARM_ModelParamsSFRM*)GetModelParams())->InterpolateCorrelation((*itsFwdResetTimes)[k]);
			for(size_t j=0; j<factorsNb; ++j )
			{
				int signeCorrel = ( ((*correl_k)[j] > 0) ? 1 : (((*correl_k)[j] < 0) ? -1 : 0) );
				(*localStdDev[timeIndex])(k,j) *= signeCorrel;
			}
			delete correl_k;
		}
		++timeIndex;
	}
}




////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: GetFloatDateStrip
///	Returns: a datestrip according to the floating convention
///	Action : computes the datestrip corresponding to the index and the 
///			default convention!
////////////////////////////////////////////////////
ARM_DateStrip* ARM_SFRM::GetFloatDateStrip( const ARM_Date& startDate, const ARM_Date& endDate ) const
{
	ARM_IRIndex* pIndex	= ((ARM_ModelParamsSFRM*)GetModelParams())->GetIRIndex();
	
	//// AAAAAAAAAAAAARRRRRRRRRGGGGGGGGGGG ugly const_cast
	ARM_Currency* pCcy	= const_cast<ARM_SFRM* const>(this)->GetZeroCurve()->GetCurrencyUnit();
	int floatResetFreq	= pIndex->GetResetFrequency();
	int floatDayCount	= pCcy->GetLiborIndexDayCount();
	char* floatResetCal = pCcy->GetResetCalName(pIndex->GetIndexType());
	int fwdRule			= pCcy->GetFwdRule();
	int intRule			= pIndex->GetIntRule();
	int stubRule		= K_SHORTSTART;
	int resetGap		= pIndex->GetResetGap();
	int floatPayFreq	= pIndex->GetPayFrequency();
	int floatPayGap		= pIndex->GetPayGap();
	char* floatPayCal	= pCcy->GetPayCalName(pIndex->GetIndexType());

#if defined(__GP_STRICT_VALIDATION)
	/// checkt that the index and 
	if(strcmp(GetZeroCurve()->GetCurrencyUnit()->GetCcyName(),
			pIndex->GetCurrencyUnit()->GetCcyName() ) != 0 )
	{
		CC_Ostringstream os;
		os	<< ARM_USERNAME << ": curve with " << GetZeroCurve()->GetCurrencyUnit()->GetCcyName()
			<< " ccy while index with " << pIndex->GetCurrencyUnit()->GetCcyName() << " ccy";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str());
	}
#endif

	/// create datestrip
	ARM_DateStrip* pFloatDateStrip = new ARM_DateStrip( startDate, endDate, floatResetFreq, floatDayCount, floatResetCal,
		fwdRule, intRule, stubRule, resetGap, floatPayFreq, floatPayGap, floatPayCal );

	delete floatResetCal;
	delete floatPayCal;

	return pFloatDateStrip;
}
/*////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: PricingTimeSteps
///	Returns :
///	Action  : returns PricingTimeSteps
////////////////////////////////////////////////////
std::vector<double>& ARM_SFRM::PricingTimeSteps(const ARM_TimeInfoPtrVector& timeInfos)
{
	ARM_DiscretisationScheme& discretisationScheme = ARM_EventAndModelTime();
	std::vector<double>& ptimeSteps = discretisationScheme.ModelTimesFromTimesInfo(timeInfos, *this );
	return ptimeSteps;
}*/


////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: ComputeModelTimes
///	Returns :
///	Action  : 
///		This method computes the SFRM forward scheduler based on the time infos of the generic security and
///		the natural datestrip of the model.
///		In this scheduler there are two kind of forwards:
///			_ the standard forward used in the diffusion to update the drift
///			_ the auxiliary forward just used for the pricing (all the reset dates inside a natural period)
///		1) Compute the end date of the generic security 
///		2) Get the market conventions for the date strip based on the index and the currency
///		3) We adjust the end date to make sure the last forward of the scheduler is standard
///		4) We compute the datestrip of the model
///		5) We merge reset dates and date strip with the following rules:
///				_ when the  reset date is closed to a datestrip date we keep the last one
///				_ in the other case we create an auxiliary forward
///		6) We compute start dates and end dates of the forward based on the reset dates. For standard 
///		forward the start date match with the end date of the previous one
///
////////////////////////////////////////////////////
std::vector<double>& ARM_SFRM::ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos )
{
	/// find the max pay time
	double endTime = 0.0;
	int i;
	
	/// 1) Compute the end date of the generic security 
	for(i=0; i<timeInfos.size(); ++i)
	{
		const std::vector<double>& payTimes = timeInfos[i]->GetPayTimes();
		/// uses the fact that the pay times are sorted hence the last one is  the biggest one
		if( endTime< payTimes[payTimes.size()-1])
			endTime = payTimes[payTimes.size()-1];
	}
	
	///	2) Get the market conventions for the date strip based on the index and the currency
	double startTime = timeInfos[0]->GetEventTime();
	double asOfDate  = GetAsOfDate().GetJulian();
	ARM_Date startDate;

	ARM_IRIndex* pIndex	= ((ARM_ModelParamsSFRM*)GetModelParams())->GetIRIndex();
	ARM_Currency* pCcy	= const_cast<ARM_SFRM* const>(this)->GetZeroCurve()->GetCurrencyUnit();
	char* floatResetCal = pCcy->GetResetCalName(pIndex->GetIndexType());
	int resetGap		= pIndex->GetResetGap();
	int floatResetFreq	= pIndex->GetResetFrequency();
	int floatDayCount	= pCcy->GetLiborIndexDayCount();

	if (itsIsFixStartDate)
	{
		startDate = itsFixStartDate;
	}
	else
	{
		startDate = ARM_Date(asOfDate+startTime);
		startDate.GapBusinessDay(-resetGap, floatResetCal );
	}
	

	/// 3) We adjust the end date to make sure the last forward of the scheduler is standard
	ARM_Date endDate;
	if (itsIsFixEndDate)
	{
		endDate = itsFixEndDate;
	}
	else 
	{
		endDate = asOfDate+endTime;

		double lastResetTime = timeInfos[timeInfos.size()-1]->GetEventTime();
		
		/// if a period is distant from more than FRMVOL_LAG_THRESHOLD 
		/// then you should add one more period.
		int resetTimePeriods = (lastResetTime-startTime + FRMVOL_LAG_THRESHOLD )/K_YEAR_LEN*floatResetFreq;
		if( lastResetTime - startTime -resetTimePeriods*K_YEAR_LEN/floatResetFreq > FRMVOL_LAG_THRESHOLD )
			++resetTimePeriods;
		ARM_Date minEndDate(startDate);
		minEndDate.AddPeriodMult(floatResetFreq,resetTimePeriods+1,pCcy->GetCcyName());
		while( minEndDate.GetJulian() < endDate.GetJulian()-FRMVOL_LAG_THRESHOLD)
		{
			minEndDate.AddPeriodMult(floatResetFreq,1,pCcy->GetCcyName());
			
		}

		if ((minEndDate-endDate)>FRMVOL_LAG_THRESHOLD)
		{
			endDate = minEndDate;
		}
		
		//We calculate the true number of period: we round because of integers.  
		int NbPeriod = (int)((endDate-startDate)/K_YEAR_LEN*floatResetFreq+0.5);
		ARM_Date nonAdjEndDate(startDate);
		nonAdjEndDate.AddPeriodMult(floatResetFreq,NbPeriod,pCcy->GetCcyName());

		endDate = nonAdjEndDate;
	}
	CC_NS(std,auto_ptr)<ARM_DateStrip> pFloatDateStrip( GetFloatDateStrip(startDate,endDate) );

	/// get Start, End Dates and interest terms
	std::vector<double>& tmpDateStripResetDates = pFloatDateStrip->GetResetDates();
	*tmpDateStripResetDates-= asOfDate;
	
	vector<double> tmpResetDates;
	vector<double> tmpStartDates;
	vector<double> tmpEndDates;
	vector<double> tmpIntTerms;
	
	tmpResetDates.reserve(tmpDateStripResetDates->size()+timeInfos.size());
	itsFwdStatus.resize(0);
	
	//FIXMEFRED: mig.vc8: no reserve in deque
	//itsFwdStatus.reserve(tmpDateStripResetDates->size()+timeInfos.size());
	
	///	5) We merge reset dates and date strip with the following rules:
	///				_ when the  reset date is closed to a datestrip date we keep the first one!
	///				_ in the other case we create an auxiliary forward
	int resetDateIndex = 0, datestripIndex;

	for (datestripIndex = 0; datestripIndex < tmpDateStripResetDates->size(); ++datestripIndex)
	{
		while ( resetDateIndex < timeInfos.size() &&
			timeInfos[resetDateIndex]->GetEventTime() < (*tmpDateStripResetDates)[datestripIndex]-FRMVOL_LAG_THRESHOLD )
		{
			tmpResetDates.push_back(timeInfos[resetDateIndex]->GetEventTime());
			itsFwdStatus.push_back(ARM_SFRM::ARM_SFRM_AUX_FWD);
			resetDateIndex++;
		}
		
		bool ContainsNoNearByResetDate=true;
		
		while(	resetDateIndex < timeInfos.size() &&
			timeInfos[resetDateIndex]->GetEventTime() < (*tmpDateStripResetDates)[datestripIndex]+FRMVOL_LAG_THRESHOLD )
		{
			if( ContainsNoNearByResetDate == false )
			{
				/// since therewas already a fwd, it means that the previous one was auxiliary!
				itsFwdStatus[itsFwdStatus.size()-1] =ARM_SFRM::ARM_SFRM_AUX_FWD;
				tmpResetDates.push_back(timeInfos[resetDateIndex]->GetEventTime());
				itsFwdStatus.push_back(ARM_SFRM::ARM_SFRM_STD_FWD);
			}
			else
			{
				tmpResetDates.push_back(timeInfos[resetDateIndex]->GetEventTime());
				itsFwdStatus.push_back(ARM_SFRM::ARM_SFRM_STD_FWD);
				ContainsNoNearByResetDate = false;
			}
			++resetDateIndex;
		}
		
		if(ContainsNoNearByResetDate)
		{
			/// we insert only if we are away from a reset date by the FRMVOL_LAG_THRESHOLD
			if(tmpResetDates.size() && (*tmpDateStripResetDates)[datestripIndex] > tmpResetDates[tmpResetDates.size()-1] +  FRMVOL_LAG_THRESHOLD )
			{
				tmpResetDates.push_back( (*tmpDateStripResetDates)[datestripIndex]);
				itsFwdStatus.push_back(ARM_SFRM::ARM_SFRM_STD_FWD);
			}
		}

	}
	
	while(  resetDateIndex < timeInfos.size() )
	{
		tmpResetDates.push_back(timeInfos[resetDateIndex]->GetEventTime());
		itsFwdStatus.push_back(ARM_SFRM::ARM_SFRM_AUX_FWD);
		++resetDateIndex;
	}

	/// 5.1) we check that the last forward is std
#if defined(__GP_STRICT_VALIDATION)
	if (itsFwdStatus[tmpResetDates.size()-1] != ARM_SFRM::ARM_SFRM_STD_FWD)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": the last fwd has to be STD!");

#endif
	/// 6) We compute start dates and end dates
	tmpStartDates.resize(tmpResetDates.size());
	for (i = 0; i < tmpResetDates.size(); ++i)
	{
		ARM_Date startDate(asOfDate+tmpResetDates[i]);
		startDate.GapBusinessDay(-resetGap, floatResetCal);
		tmpStartDates[i] = startDate.GetJulian()-asOfDate;
	}
	
	tmpEndDates.resize(tmpResetDates.size());
	tmpIntTerms.resize(tmpResetDates.size());
	
	double lastEndDate = tmpStartDates[tmpResetDates.size()-1];
	tmpEndDates[tmpResetDates.size()-1] = pFloatDateStrip->GetFlowEndDates()->Elt(pFloatDateStrip->GetFlowEndDates()->size()-1) - asOfDate;
	

	/*if( fabs(tmpEndDates[tmpResetDates.size()-1] - tmpStartDates[tmpResetDates.size()-1] - K_YEAR_LEN/floatResetFreq ) > FRMVOL_LASTLAG_THRESHOLD)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": start date and end date does not seem to be std on the last fwd!");*/

	tmpIntTerms[tmpResetDates.size()-1] = CountYears(floatDayCount, asOfDate+tmpStartDates[tmpResetDates.size()-1], asOfDate+tmpEndDates[tmpResetDates.size()-1]);


	/// we always check that end date and computedEndDate for std forward are not distant from 
	/// more than FRMVOL_LAG_THRESHOLD. If so, takes the computedEndDate!
	for (i = tmpResetDates.size()-2; i >= 0; --i)
	{
		if (itsFwdStatus[i] == ARM_SFRM::ARM_SFRM_STD_FWD)
		{
			ARM_Date computedEndDate(tmpStartDates[i]+asOfDate);
			computedEndDate.AddPeriod(floatResetFreq,pCcy->GetCcyName());
			if( fabs( computedEndDate.GetJulian()-asOfDate - lastEndDate )>FRMVOL_LAG_THRESHOLD)
				tmpEndDates[i] = computedEndDate.GetJulian()-asOfDate;
			else
				tmpEndDates[i] = lastEndDate;

			lastEndDate = tmpStartDates[i];
		}
		else if (itsFwdStatus[i] == ARM_SFRM::ARM_SFRM_AUX_FWD)
		{
			ARM_Date endDate(tmpStartDates[i]+asOfDate);
			endDate.AddPeriod(floatResetFreq,pCcy->GetCcyName());
			tmpEndDates[i] = endDate.GetJulian()-asOfDate;
			
		}
		tmpIntTerms[i] = CountYears(floatDayCount, asOfDate+tmpStartDates[i], asOfDate+tmpEndDates[i]);
	}

	// We want to make sure that the end date of the last forward is at least end date of the deal
	// description !!!
	tmpEndDates[tmpEndDates.size()-1] = CC_Max<double>(tmpEndDates[tmpEndDates.size()-1],endTime);
	tmpIntTerms[tmpResetDates.size()-1] = CountYears(floatDayCount, asOfDate+tmpStartDates[tmpResetDates.size()-1], asOfDate+tmpEndDates[tmpResetDates.size()-1]);
	
	itsFwdResetTimes  = ARM_VectorPtr(new std::vector<double>(tmpResetDates));
	itsFwdStartTimes  = ARM_VectorPtr(new std::vector<double>(tmpStartDates));
	itsFwdEndTimes	  = ARM_VectorPtr(new std::vector<double>(tmpEndDates));
	itsFwdIntTerms	  = ARM_VectorPtr(new std::vector<double>(tmpIntTerms));
	
	delete floatResetCal;
	
	return static_cast<ARM_GP_Vector*>(itsFwdResetTimes->Clone());
}


/// bad_alloc handler
int bad_alloc_handler(size_t) {
  throw std::bad_alloc();
  return 0;
}

////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: FirstPricingStates
///	Returns : 
///	Action  : create the first pricing state
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_SFRM::FirstPricingStates(  size_t bucketSize ) const
{
	/// Change the allocation handler
	_PNH _old_new_handler;
	_old_new_handler = _set_new_handler(bad_alloc_handler);

	ARM_PricingStatesPtr initStates;
	ResetDFMap();
#ifdef __GP_SFRM_STOCH_TERMS
	ARM_GP_Matrix* stochasticTerms = NULL;
#endif

	try
	{
		const size_t nbPayoffs=0;
		size_t nbModelStates=itsFwdValues->size();
		size_t factorsNb		= FactorCount();
		initStates = ARM_PricingStatesPtr( new ARM_PricingStates(bucketSize,nbModelStates,nbPayoffs,factorsNb) );
		initStates->EnqueuePricingStatesContext( ARM_PricingStatesContextPtr( new ARM_SFRMPricingStatesContext () ) );

#ifdef __GP_SFRM_STOCH_TERMS
		stochasticTerms = new ARM_GP_Matrix(bucketSize,nbModelStates);
#endif

		for(size_t i=0; i<bucketSize; ++i )
			for(size_t j=0;j<nbModelStates; ++j)
			{
				initStates->SetModelState(i,j,(*itsFwdValues)[j]);
#ifdef __GP_SFRM_STOCH_TERMS
				(*stochasticTerms)(i,j) = 0.0;
#endif
			}
	}
	catch (bad_alloc e)
	{
		// Reuse the default allocation handler
		_set_new_handler(_old_new_handler);

#ifdef __GP_SFRM_STOCH_TERMS
		delete stochasticTerms;
		stochasticTerms = NULL;
#endif
		throw Exception(__LINE__, __FILE__,ERR_INVALID_DATA, ARM_USERNAME +   ": its was impossible to allocate the pricing state.");
	}

#ifdef __GP_SFRM_STOCH_TERMS
///	delete CC_MUTABLE( ARM_SFRM, itsStochasticTerms );
	delete itsStochasticTerms;
	CC_MUTABLE( ARM_SFRM, itsStochasticTerms ) = stochasticTerms;
#endif

	_set_new_handler(_old_new_handler);

	// SFRM states context creation
	initStates->EnqueuePricingStatesContext( ARM_PricingStatesContextPtr( new ARM_SFRMPricingStatesContext() ));

	// Fwd Min Index Init
	if(GetNumMethod()->GetPricingDirection() ==  ARM_NumMethod::GP_FWDLOOKING || 
		GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDBCKWDLOOKING )
    {
		initStates->GetPricingStatesContext( GetModelRank() )->ToSFRMPricingStatesContext()->SetFwdMinIndex( 0 );
	}
	else if(GetNumMethod()->GetPricingDirection() ==  ARM_NumMethod::GP_BCKWDLOOKING)
    {			
		initStates->GetPricingStatesContext( GetModelRank() )->ToSFRMPricingStatesContext()->SetFwdMinIndex( itsFwdEndTimes->size() -1 );
    }


	return initStates;
}

////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: ComputeNumeraireTimeIndex
///	Returns : int
///	Action  : Compute the numeraire time index based on the numeraire
/// current date. 
////////////////////////////////////////////////////
void ARM_SFRM::ComputeNumeraireTimeIndex() const
{
	double numeraireTime = GetNumeraire()->GetMaturity();	
	itsNumeraireTimeIndex = lower_boundPosWithPrecision(*itsFwdEndTimes,numeraireTime);
	while((itsNumeraireTimeIndex <= itsFwdStatus.size()) && 
		itsFwdStatus[itsNumeraireTimeIndex-1] != ARM_SFRM::ARM_SFRM_STD_FWD)
	{
		CC_MUTABLE( ARM_SFRM, itsNumeraireTimeIndex )++;
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: MCModelStatesFromToNextTime
///	Returns : void
///	Action  : No default implementation since inappropriate for SFRM  model
////////////////////////////////////////////////////
void ARM_SFRM::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
		throw Exception(__LINE__, __FILE__,ERR_INVALID_ARGUMENT,"time index is negative!");
	if( timeIndex >= GetNumMethod()->GetTimeSteps()->size()-1)
		throw Exception(__LINE__, __FILE__,ERR_INVALID_ARGUMENT,"time index bigger than max size!");
	
	/// test pricing method direction
	if(ARM_NumMethod::GP_FWDLOOKING != GetNumMethod()->GetPricingDirCurrLoop() )
		throw Exception(__LINE__, __FILE__,ERR_INVALID_ARGUMENT,"implemented only for forward looking method..More precisely MC method!");


#endif
	ResetDFMap();
	double fromTime = GetNumMethod()->GetTimeStep(timeIndex);
	double nextTime = GetNumMethod()->GetTimeStep(timeIndex+1);
	
	const ARM_MatrixVector& numMethLocalVar	= GetNumMethodStateLocalVars();
	const ARM_MatrixVector& numMethLocalStdDev	= GetNumMethodStateLocalStdDevs();
	const ARM_MatrixVector& modelLocalVar	= GetModelStateLocalVars();
	const ARM_MatrixVector& modelLocalStdDev = GetModelStateLocalStdDevs();

	/// manip to reset the fwdMinIndex
	/// otherwise because this function is called successively
	/// it is performancewise efficient to keep track of the previous
	///	fwdMinIndex. Its value is stored in the SFRM pricing states.
	int fwdMinIndex = states->GetPricingStatesContext( GetModelRank() )->ToSFRMPricingStatesContext()->GetFwdMinIndex();

	while( fwdMinIndex < itsFwdResetTimes->size() 
		&& (*itsFwdResetTimes)[fwdMinIndex]<nextTime)
	{
		++fwdMinIndex;
	}

	states->GetPricingStatesContext( GetModelRank() )->ToSFRMPricingStatesContext()->SetFwdMinIndex( fwdMinIndex );

	if ((GetNumeraire()->GetType() == ARM_Numeraire::RollingEvent) ||
	   (GetNumeraire()->GetType() == ARM_Numeraire::RollingPayment))
	   ComputeNumeraireTimeIndex();

#if defined(__GP_STRICT_VALIDATION)
	/// test fwd min index
	if(  itsNumeraireTimeIndex < fwdMinIndex )
		throw Exception(__LINE__, __FILE__,ERR_INVALID_ARGUMENT, ARM_USERNAME +  ": the numeraire Fwd is dead!");
#endif
	
	int fwdMaxIndex =  itsFwdResetTimes->size();

	size_t factorsNb		= FactorCount();
	
	double shift,delta,fwd,cst,drift,martPart,modelStdDev,stdDev,var;
#ifdef __GP_SFRM_STOCH_TERMS
	double stochasticTerm;
#endif
	size_t nbStates		= states->size();

	/// keeps track of the last standard fwd for auxiliary fwds
	const int UNITIALISED = -1;
	int lastStandardIndex = UNITIALISED;

	std::vector<double> lastDriftPerFactor(factorsNb);
	double deltaDrift = 0.0;
	ARM_Currency* pCcy	= const_cast<ARM_SFRM* const>(this)->GetZeroCurve()->GetCurrencyUnit();
	int floatDayCount	= pCcy->GetLiborIndexDayCount();
	
	/// \latexonly
	/// the diffusion of the forwards in the SFRM model under the probability measure $Q^{N+1}$ 
	/// are given by:
	/// \begin{equation}
	/// \frac{dF_{k}\left( t\right) }{F_{k}\left( t\right) +m_{k}}=\mu _{k}\left(
	/// t\right) dt+\overrightarrow{\sigma }_{k}\left( t\right) d\overrightarrow{W}%
	/// _{t}^{N+1}  \label{eqtn7}
	/// \end{equation}
	/// 
	/// where the drift is given by
	/// 
	/// \begin{equation}
	/// \mu _{k}\left( t\right) =-\overrightarrow{\sigma }_{k}\left( t\right) 
	/// \underset{k+1\leq i\leq N}{\sum }\frac{\delta _{i}\left( F_{i}\left(
	/// t\right) +m_{i}\right) }{1+\delta _{i}F_{i}\left( t\right) }\ 
	/// \overrightarrow{\sigma }_{i}\left( t\right)  \label{eqtn8}
	/// \end{equation}
	/// 
	/// \endlatexonly

	const ARM_NumMethodPtr numMethod= GetNumMethod();
	//ARM_GP_MatrixPtr absoluteDrift	= numMethod->GetNumMethodAbsoluteDrifts();
	//bool useAbsoluteDrift= absoluteDrift!= ARM_GP_MatrixPtr(NULL);

	size_t modelNb = GetModelNb();
	for(size_t k=0; k<nbStates; ++k )
	{
		size_t j; /// reinit the sum
		for(j=0;j<factorsNb; ++j)
			(*itsFactors)[j] = 0;

		int i;
		for(i=itsNumeraireTimeIndex-1; i>=fwdMinIndex; --i)
		{
			shift = (*itsFwdShifts)[i];
			delta = (*itsFwdIntTerms)[i];
		
			/// get previous model states
			fwd		= states->GetModelState(k,i+modelNb);
			cst		= (delta*(fwd+shift))/(1.0+delta*fwd);
			
#ifdef __GP_SFRM_STOCH_TERMS
			stochasticTerm = GetStochasticTerm(k,i);
#endif
			
			/// initialise
			drift	= 0.0;
			martPart= 0.0;

			if (itsFwdStatus[i] == ARM_SFRM::ARM_SFRM_AUX_FWD)
			{
				if( lastStandardIndex == UNITIALISED )
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +
						": error in the schedule... it should have added a std fwd and it did not!" );
				if ((*itsFwdStartTimes)[lastStandardIndex] <= (*itsFwdEndTimes)[i])
					deltaDrift = CountYears(floatDayCount,(*itsFwdStartTimes)[lastStandardIndex],(*itsFwdEndTimes)[i])
							/CountYears(floatDayCount,(*itsFwdStartTimes)[lastStandardIndex],(*itsFwdEndTimes)[lastStandardIndex]);
				else
					deltaDrift = 0.0;
			}

				/// computation of the martingale part!
			for(j=0;j<factorsNb; ++j)
			{
				var	      = (*modelLocalVar[timeIndex])(i,j)*(*numMethLocalVar[timeIndex])(j,j);
				modelStdDev = (*modelLocalStdDev[timeIndex])(i,j);
				stdDev    = (*modelLocalStdDev[timeIndex])(i,j)*(*numMethLocalStdDev[timeIndex])(j,j);
				drift    += (*itsFactors)[j]*stdDev;				

				if (itsFwdStatus[i] == ARM_SFRM::ARM_SFRM_AUX_FWD)
				{
					drift -= lastDriftPerFactor[j]*stdDev*deltaDrift;
				}

				double gaussian = states->GetNumMethodState(k,j+modelNb);
				martPart += modelStdDev*gaussian-0.5*var;//gaussian[k*factorsNb+j]
			}
		
			/*if( useAbsoluteDrift )
				drift += (*absoluteDrift)(timeIndex,modelNb+i);*/

#ifdef __GP_SFRM_STOCH_TERMS
			stochasticTerm +=-drift+martPart;
#endif
			fwd=(fwd+shift)*exp(-drift+martPart)-shift;
			states->SetModelState(k,i+modelNb,fwd);
#ifdef __GP_SFRM_STOCH_TERMS
			const_cast<ARM_SFRM*>(this)->SetStochasticTerm(k,i+modelNb,stochasticTerm);
#endif

			if (itsFwdStatus[i] == ARM_SFRM::ARM_SFRM_STD_FWD)
			{
				for(j=0;j<factorsNb; ++j)
				{
					stdDev    = (*modelLocalStdDev[timeIndex])(i,j)*(*numMethLocalStdDev[timeIndex])(j,j);
					(*itsFactors)[j]   += cst*stdDev;
					lastDriftPerFactor[j] = cst*stdDev;
				}
				lastStandardIndex = i;
			}
		}

		if(fwdMaxIndex>itsNumeraireTimeIndex)
		{
			/// case of a numeraire with a date smaller than the fwd current index
			/// reinit the factors
			for(j=0;j<factorsNb; ++j)
				(*itsFactors)[j] = 0.0;

			/// then does the diffusion for the fwd after the fwd max index
			for(i=itsNumeraireTimeIndex; i<fwdMaxIndex; ++i)
			{

				shift = (*itsFwdShifts)[i];
				delta = (*itsFwdIntTerms)[i];
			
				/// get previous model states
				fwd		= states->GetModelState(k,i+modelNb);
#ifdef __GP_SFRM_STOCH_TERMS
				stochasticTerm = GetStochasticTerm(k,i);
#endif
				cst		= (delta*(fwd+shift))/(1.0+delta*fwd);

				if (itsFwdStatus[i] == ARM_SFRM::ARM_SFRM_STD_FWD)
				{
					for(j=0;j<factorsNb; ++j)
					{
						var	      = (*modelLocalVar[timeIndex])(i,j)*(*numMethLocalVar[timeIndex])(j,j);
						modelStdDev = (*modelLocalStdDev[timeIndex])(i,j);
						stdDev    = (*modelLocalStdDev[timeIndex])(i,j)*(*numMethLocalStdDev[timeIndex])(j,j);
						(*itsFactors)[j]   += cst*stdDev;
						lastDriftPerFactor[j] = cst*stdDev;
			
					}
					lastStandardIndex = i;
				}
				
				/// initialise
				drift	= 0.0;
				martPart= 0.0;

				if (itsFwdStatus[i] == ARM_SFRM::ARM_SFRM_AUX_FWD)
				{
					if( lastStandardIndex == UNITIALISED )
						throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +
							": error in the schedule... it should have added a std fwd and it did not!" );

					if ((*itsFwdStartTimes)[i] >= (*itsFwdEndTimes)[lastStandardIndex])
						deltaDrift = CountYears(floatDayCount,(*itsFwdStartTimes)[i],
															  (*itsFwdEndTimes)[lastStandardIndex])/CountYears(floatDayCount,
															  (*itsFwdStartTimes)[lastStandardIndex],
															  (*itsFwdEndTimes)[lastStandardIndex]);
				}
				
				for(j=0;j<factorsNb; ++j)
				{
					/// computation of the martingale part!
					var	      = (*modelLocalVar[timeIndex])(i,j)*(*numMethLocalVar[timeIndex])(j,j);
					modelStdDev = (*modelLocalStdDev[timeIndex])(i,j);
					stdDev    = (*modelLocalStdDev[timeIndex])(i,j)*(*numMethLocalStdDev[timeIndex])(j,j);
									
					drift    += (*itsFactors)[j]*stdDev;

					if (itsFwdStatus[i] == ARM_SFRM::ARM_SFRM_AUX_FWD)
					{
						drift -= lastDriftPerFactor[j]*stdDev*deltaDrift;
					}

					double gaussian = states->GetNumMethodState(k,j+modelNb);
					martPart += modelStdDev*gaussian-0.5*var;

				}
#ifdef __GP_SFRM_STOCH_TERMS			
				stochasticTerm +=+drift+martPart;
#endif
				fwd=(fwd+shift)*exp(+drift+martPart)-shift;
				states->SetModelState(k,i+modelNb,fwd);
#ifdef __GP_SFRM_STOCH_TERMS
				const_cast<ARM_SFRM*>(this)->SetStochasticTerm(k,i+modelNb,stochasticTerm);
#endif
			}
		}
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: PrepareForMarkovDriftCalculation
///	Returns:  
///	Action : Prepare to calculate the markovian drifts 
////////////////////////////////////////////////////
void ARM_SFRM::PrepareForMarkovDriftCalculation(int timeIndex, double var_At_timeIndex,
	int itsFwdMinIndex,	 int NbFwds, int nbIntegrationSteps, std::vector<double>&& VarCoeff,
	std::vector<double>&& SpotVarCoeff, ARM_GP_Matrix*& FwdVolCoeff, 
	double*& x, double*& w) const
{
    GaussLegendre_Coefficients GL(nbIntegrationSteps);
	x = new double [nbIntegrationSteps];
	w = new double [nbIntegrationSteps];
	double fromTime = GetNumMethod()->GetTimeStep(timeIndex);
	double a=0.5*fromTime;
    double b=a/K_YEAR_LEN;
    for(int i=0;i<nbIntegrationSteps;++i)
    {
        x[i] = GL.get_point(i);
		x[i] = (x[i]+1.0)*a;
		w[i] = GL.get_weight(i);
        w[i] *= b;
    }

	VarCoeff=new std::vector<double>(nbIntegrationSteps);
    SpotVarCoeff=new std::vector<double>(nbIntegrationSteps);
	FwdVolCoeff=new ARM_GP_Matrix(nbIntegrationSteps,NbFwds);

	for(i=0;i<nbIntegrationSteps;++i)
    {
		//calculus of VarCoeff 
		if(K_NEW_DOUBLE_TOL < x[i] &&
				x[i] < fromTime-K_NEW_DOUBLE_TOL)
		{
			VarCoeff->Elt(i)= ((ARM_ModelParamsSFRM*)GetModelParams())->IntegratedLocalVariance(0.0,x[i]);
		}
		else if(x[i] <= K_NEW_DOUBLE_TOL)
			VarCoeff->Elt(i)=0.0;
		else
			VarCoeff->Elt(i)=var_At_timeIndex;

		//calculus of 
		SpotVarCoeff->Elt(i)= ((ARM_ModelParamsSFRM*)GetModelParams())->VolatilitySpotSquared(x[i]);

		double curVarCoeff=VarCoeff->Elt(i);

        if(curVarCoeff > K_NEW_DOUBLE_TOL)
        {
			for(int fwdIdx=0;fwdIdx<NbFwds;++fwdIdx)
            {
				double correl = (*(const_cast<ARM_SFRM *>(this)->itsFwdsVarCoeffs))[fwdIdx+itsFwdMinIndex];
				FwdVolCoeff->Elt(i,fwdIdx)=exp(curVarCoeff/var_At_timeIndex*
								0.5*correl*(var_At_timeIndex-curVarCoeff));
			}
		}
		else
		{
			for(int fwdIdx=0;fwdIdx<NbFwds;++fwdIdx)
				FwdVolCoeff->Elt(i,fwdIdx)=1.0;
		}
	}   
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: TreeStatesToModelStatesDiffusion
///	Returns: void
///	Action : This function transform the tree states
/// into the model states, it means the fwds
////////////////////////////////////////////////////

void ARM_SFRM::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{
	(this->*itsTreeStatesToModelStatesFunc)(states, timeIndex);
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: NeedModelTime
///	Returns: void
///	Action : Do we need to add the model time in 
/// discretisation scheme of the numerical method ?
////////////////////////////////////////////////////
bool ARM_SFRM::NeedModelTime() const
{
// FIXMEFRED: mig.vc8 (30/05/2007 16:03:20):need pointer
	return itsTreeStatesToModelStatesFunc==&ARM::ARM_SFRM::TreeStatesToModelStatesNonParamDrift;
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: NeedModelTime
///	Returns: void
///	Action : Do we need to evaluate the states at 
/// this timeIndex? 
////////////////////////////////////////////////////
// 
bool ARM_SFRM::NeedStatesEval(int timeIndex) const
{
	if (NeedModelTime())
	{
		double timeStep = GetNumMethod()->GetTimeStep(timeIndex);
		for (int i = 0; i < itsFwdResetTimes->size(); i++)
		{
			if (fabs((*itsFwdResetTimes)[i]-timeStep) < K_DOUBLE_TOL)
				return true;
		}
	}

	return false;
}


////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: CreateProbaChange
///	Returns: void
///	Action : This function calculate the integrale of 
/// the probability change with a Gauss Legendre method
/// of the drift
////////////////////////////////////////////////////

double ARM_SFRM::IntegrateProbaChange(
	int n,
	double stdDev,
	const std::vector<double>& states,
	const std::vector<double>& probaChanges,
	bool fromMinMax,
	double MinMax,
	double x0
	) const
{
	double min, max;

	if (fromMinMax)
	{
		min = MinMax;
		max = x0;
	}
	else
	{
		min = x0;
		max = MinMax;
	}

	if (itsLegendreCoeffs.IsNull() || (itsLegendreCoeffs->get_order() != n))
	{
		itsLegendreCoeffs = ARM_CountedPtr<GaussLegendre_Coefficients>(new  GaussLegendre_Coefficients(n,-1,1));
	}

	int nbStates = states.size();

	double x, w, value;
	int stateIndex = 0;

	double delta = (max-min)/2;
	double mid = (max+min)/2;

	double integrale = 0.0;

	for (int i = 0; i < n; i++)
	{
		x = mid + delta*itsLegendreCoeffs->get_point(i);
		w = delta*itsLegendreCoeffs->get_weight(i);

		while ((stateIndex < nbStates-1) && (states[stateIndex]<=x))
			stateIndex++;

		if (stateIndex > 0)
			value = probaChanges[stateIndex-1] + (x-states[stateIndex-1])/(states[stateIndex]-states[stateIndex-1])*(probaChanges[stateIndex]-probaChanges[stateIndex-1]);
		else
			value = probaChanges[stateIndex] + (x-states[stateIndex])/(states[stateIndex+1]-states[stateIndex])*(probaChanges[stateIndex+1]-probaChanges[stateIndex]);


		integrale += 1.0/sqrt(2*PI)*value*exp(-x*x/2.0/stdDev/stdDev)/stdDev*w;
	}

	return integrale;
}


////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: CreateProbaChange
///	Returns: void
///	Action : Create the new probability change
/// based on the previous
/// of the drift
////////////////////////////////////////////////////

void ARM_SFRM::CreateProbaChange(
		ARM_PricingStatesPtr& states,
		const ARM_GP_VectorPtr& proba,
		int fwdIdx,
		double var) const
{
	double shift = (*itsFwdShifts)[fwdIdx];
	double delta = (*itsFwdIntTerms)[fwdIdx];
	double fwd0 = (*itsFwdValues)[fwdIdx];
	double forwardVolCoeff = (*(const_cast<ARM_SFRM *>(this)->itsFwdsVarCoeffs))[fwdIdx];

	double proba7 = ARM_GaussianAnalytics::cdfNormal2(-7.0);

	states->resizeProbaChanges(states->GetProbaChangesSize()+1);

	int prevProbaChangePos = states->GetProbaChangesSize()-2;
	int probaChangePos = states->GetProbaChangesSize()-1;

	size_t i,j;

	double norm = 0.0, probaVal = 0.0, x = 0.0;
	size_t modelNb = GetModelNb();
	size_t factorsNb = FactorCount();
	size_t nbStates = states->size();
	std::vector<double> FwdsCorrelVarCoeffs(factorsNb);

	double stdDev = sqrt(var), volCoeff = 0.0;

	for(j=0;j<factorsNb; ++j)
	{
		FwdsCorrelVarCoeffs.Elt(j) = (*(const_cast<ARM_SFRM *>(this)->itsFwdsCorrelVarCoeffs[j]))(fwdIdx,fwdIdx);
		volCoeff += FwdsCorrelVarCoeffs.Elt(j)*FwdsCorrelVarCoeffs.Elt(j);
	}

	std::vector<double> expTab(nbStates), stochasticTerm(nbStates), probaChanges(nbStates), decStates(nbStates), decStates1(nbStates), decStates2(nbStates);

	for (i = 0; i < nbStates; ++i)
	{
		stochasticTerm[i]=0;
		for(j=0;j<factorsNb; ++j)
		{
			double state_j = states->GetNumMethodState(i,j+modelNb);
			stochasticTerm[i] += (state_j*FwdsCorrelVarCoeffs[j]);
		}

		if (prevProbaChangePos > -1)
			probaChanges[i] = states->GetProbaChange(i,prevProbaChangePos);
		else
			probaChanges[i] = 1.0;
	}

	double dx;

	const double nbStdDev = 10;
	const int nbLegendrePoints = 8;

	std::vector<double> stochTermVec2(2);
	std::vector<double> probaChangesVec2(2);

	double sumPrevProbaChange = 0.0;
	sumPrevProbaChange += ARM_GaussianAnalytics::cdfNormal2(stochasticTerm[0]/stdDev)*probaChanges[0];
	for (i = 1; i < nbStates; ++i)
	{
		stochTermVec2[0] = stochasticTerm[i-1];
		stochTermVec2[1] = stochasticTerm[i];
		probaChangesVec2[0] = probaChanges[i-1];
		probaChangesVec2[1] = probaChanges[i];
		sumPrevProbaChange += IntegrateProbaChange(
									8,
									stdDev,
									stochTermVec2,
									probaChangesVec2,
									true,
									stochasticTerm[i-1],
									stochasticTerm[i]);
	}
	sumPrevProbaChange += ARM_GaussianAnalytics::cdfNormal2(-stochasticTerm[nbStates-1]/stdDev)*probaChanges[nbStates-1];

	ostringstream os;

	os << "SumProbaChange" << prevProbaChangePos << " : " << setprecision(8) << sumPrevProbaChange << endl;

	double cumProbaChange = 0.0;

	double a, b, c, Delta;

	double integrale;

	for (i = 0; i < nbStates; ++i)
	{	
		double norm = 0.0;

		if (prevProbaChangePos > -1)
		{
			if (i > 0)
			{
				stochTermVec2[0] = stochasticTerm[i-1];
				stochTermVec2[1] = stochasticTerm[i];
				probaChangesVec2[0] = probaChanges[i-1];
				probaChangesVec2[1] = probaChanges[i];


				integrale = IntegrateProbaChange(
									nbLegendrePoints,
									stdDev,
									stochTermVec2,
									probaChangesVec2,
									true,
									stochasticTerm[i-1],
									stochasticTerm[i]);

				cumProbaChange += integrale/sumPrevProbaChange;
			}
			else
			{
				integrale = ARM_GaussianAnalytics::cdfNormal2(stochasticTerm[0]/stdDev)*probaChanges[0]/sumPrevProbaChange;
				cumProbaChange = integrale;
			}

			if ((cumProbaChange > proba7) && (cumProbaChange < 1.0-proba7))
				norm = ARM_NormalInvCum::Inverse_erf_Moro(cumProbaChange)*stdDev;
			else if (cumProbaChange  < proba7)
				norm = -7.0*stdDev;
			else if (cumProbaChange > 1.0-proba7)
				norm = 7.0*stdDev;

			dx = stochasticTerm[1]-stochasticTerm[0];

			a = 0.5/var;
			b = -stochasticTerm[i]/var;
			c = log(probaChanges[i]);// log(integrale)+0.5*stochasticTerm[i]/var-log(dx/sqrt(2*PI)/stdDev);

			Delta = sqrt(b*b-4*a*c);

			decStates1[i] = (-b-Delta)/2/a;
			decStates2[i] = (-b+Delta)/2/a;

			decStates[i] = norm - stochasticTerm[i];
		}
		else
		{
			norm = stochasticTerm[i];
		}

		expTab[i] = exp(-0.5*var*volCoeff+volCoeff*norm);
	}
	
	double sumProba = 0.0, sumExp = 0.0;

	std::vector<double> newProbaChanges(nbStates);

	const double maxProbaChange = 1000;

	for (i = 0; i < nbStates; ++i)
	{
		newProbaChanges[i] = (1+delta*((fwd0+shift)*expTab[i]-shift))*probaChanges[i]/(1.0+delta*fwd0);
	}

	double sumNewProbaChange = 0;
	sumNewProbaChange += ARM_GaussianAnalytics::cdfNormal(stochasticTerm[0]/stdDev)*newProbaChanges[0];
	for (i = 1; i < nbStates; ++i)
	{
		stochTermVec2[0] = stochasticTerm[i-1];
		stochTermVec2[1] = stochasticTerm[i];
		probaChangesVec2[0] = newProbaChanges[i-1];
		probaChangesVec2[1] = newProbaChanges[i];
		sumNewProbaChange += IntegrateProbaChange(
									8,
									stdDev,
									stochTermVec2,
									probaChangesVec2,
									true,
									stochasticTerm[i-1],
									stochasticTerm[i]);
	}
	sumNewProbaChange += ARM_GaussianAnalytics::cdfNormal(-stochasticTerm[nbStates-1]/stdDev)*newProbaChanges[nbStates-1];

	os << "SumProbaChange" << probaChangePos << " : " << setprecision(8) << sumNewProbaChange << endl;

	os << "var = " << var << endl;

	for (i = 0; i < nbStates; ++i)
	{
		newProbaChanges[i] /= sumNewProbaChange;
		states->SetProbaChange(i,probaChangePos, newProbaChanges[i]);
	}

	for (i = 0; i < nbStates; ++i)
	{
		os << stochasticTerm[i] << "\t" << newProbaChanges[i] << "\t" <<  decStates[i]  << "\t" << decStates1[i] << "\t" << decStates2[i] << endl;
	}

	ARM_TheEventViewer.Instance()->AddToMessage(os.str());
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: TreeStatesToModelStatesDiffusion
///	Returns: void
///	Action : This function used the diffusion
/// of the drift
////////////////////////////////////////////////////

void ARM_SFRM::TreeStatesToModelStatesNonParamDrift(ARM_PricingStatesPtr& states, int timeIndex) const
{
	ARM_NumerairePtr numeraire = GetNumeraire();
#if defined(__GP_STRICT_VALIDATION)
	if(	numeraire->GetType() != ARM_Numeraire::TerminalZc )
    {		
		CC_Ostringstream os;
		os  << ARM_USERNAME << ": only TerminalZc or EventZC numeraire is "
			<< " supported by the tree method in SFRM at the moment!";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );			
    }
#endif
	if ( numeraire->GetType() == ARM_Numeraire::RollingEvent ||
		 numeraire->GetType() == ARM_Numeraire::RollingPayment )
		ComputeNumeraireTimeIndex();

	int fwdMinIndex = states->GetPricingStatesContext( GetModelRank() )->ToSFRMPricingStatesContext()->GetFwdMinIndex();
	ResetDFMap();
	/// Validation tests
#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
		throw Exception(__LINE__, __FILE__,ERR_INVALID_ARGUMENT,"time index is negative!");
	
	/// test pricing direction
	if(ARM_NumMethod::GP_BCKWDLOOKING != GetNumMethod()->GetPricingDirection() )
		throw Exception(__LINE__, __FILE__,ERR_INVALID_ARGUMENT,"implemented only for backward looking method..More precisely Tree method!");
#endif

	
	//we diffuse first the forward corresponding to the numeraire zero coupon 
	int numeraireFwdIdx = itsNumeraireTimeIndex-1; 
	if (itsFwdStatus[numeraireFwdIdx] == ARM_SFRM::ARM_SFRM_AUX_FWD)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +
			": error in the schedule... the last forward has to be standard!" );
	}
	
	///	If timeIndex=0
	int stateIdx, fwdIdx; 
	size_t NbStates	= states->size();
	size_t modelNb = GetModelNb();
	int fwdMaxIndex =  itsFwdResetTimes->size()-1;
	size_t NbFwds  = fwdMaxIndex+1;
	double correctedDrift;

	std::vector<double> savAcc(NbStates,1.0), prevSavAcc(NbStates,1.0);

	states->resizeModelStates(ModelStatesSize(), NbStates);

	if (timeIndex==0)
	{
		for(stateIdx=0; stateIdx<NbStates; ++stateIdx )
			for(fwdIdx=NbFwds-1; fwdIdx>=fwdMinIndex; fwdIdx--)
				states->SetModelState(stateIdx,modelNb+fwdIdx,0.0);
	}
	
	else
	{
		/// Preliminaries
		/// get data
		double fromTime = GetNumMethod()->GetTimeStep(timeIndex);		
		size_t factorsNb = FactorCount();
		
		/// infered data 
		while( (fwdMinIndex > 0) && (*itsFwdResetTimes)[fwdMinIndex]> fromTime)
			--fwdMinIndex;

		states->GetPricingStatesContext( GetModelRank() )->ToSFRMPricingStatesContext()->SetFwdMinIndex( fwdMinIndex );
		
		size_t NbFwdsToDiffuse  = (NbFwds-fwdMinIndex) < 0.0 ? 0.0 : NbFwds-fwdMinIndex;
		
		/// used variables
		ARM_ZeroCurvePtr zeroCurve = GetZeroCurve();
		ARM_GP_VectorPtr probaSpot;
		const ARM_MatrixVector& globalVars  = GetNumMethod()->GetNumMethodStateGlobalVars();
		double var_At_timeIndex = 0.0;		
		for(int j=0;j<factorsNb; ++j) //variance norm at timeIdx
			var_At_timeIndex += (*globalVars[timeIndex])(j,j);

		double numeraireTime = GetNumeraire()->GetMaturity();
		double numeraireBond=zeroCurve->DiscountPrice(numeraireTime/K_YEAR_LEN);
		
		double shift, delta, FWD_0, ShiftedFWD_0, FWD, martPart, prevMartPart, forwardVolCoeff;


		//states->resizeProbaChanges(0,0);

		int prevProbaChangePos = states->GetProbaChangesSize();

		for (int i = prevProbaChangePos; i < numeraireFwdIdx-fwdMinIndex+1; ++i)
		{
			CreateProbaChange(
				states,
				probaSpot,
				numeraireFwdIdx-i,
				var_At_timeIndex);
		}


		if (states->GetProbaChangesSize() != numeraireFwdIdx-fwdMinIndex+1)
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, "Proba Changes != nbFwds");
		}

			
		for(fwdIdx=numeraireFwdIdx;fwdIdx>=fwdMinIndex;fwdIdx--)
		{
			shift = (*itsFwdShifts)[fwdIdx];
			delta = (*itsFwdIntTerms)[fwdIdx];
			FWD_0 = (*itsFwdValues)[fwdIdx];
			ShiftedFWD_0 = (*itsFwdValues)[fwdIdx]+shift;
			forwardVolCoeff = (*(const_cast<ARM_SFRM *>(this)->itsFwdsVarCoeffs))[fwdIdx];
			
			correctedDrift = 0.0;
			for(stateIdx=0; stateIdx<NbStates; ++stateIdx )		
			{	
				martPart = states->GetProbaChange(stateIdx, numeraireFwdIdx-fwdIdx);
				if (fwdIdx<numeraireFwdIdx)
				{
					prevMartPart = states->GetProbaChange(stateIdx, numeraireFwdIdx-fwdIdx-1);
					if (fabs(prevMartPart) > K_NEW_DOUBLE_TOL)
						martPart /= prevMartPart;
				}

				martPart=(martPart*(1.0+delta*FWD_0)-1.0)/delta;

				martPart=(martPart+shift)/ShiftedFWD_0;
			
				states->SetModelState(stateIdx,fwdIdx+modelNb,martPart);
			}

			for(stateIdx=0; stateIdx<NbStates; ++stateIdx )		
			{
				martPart = states->GetModelState(stateIdx,fwdIdx+modelNb);
				FWD = ShiftedFWD_0*martPart-shift;
				states->SetModelState(stateIdx,fwdIdx+modelNb,FWD);
				savAcc.Elt(stateIdx) = (1.0 + delta*FWD)*savAcc.Elt(stateIdx);
			}
		}
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: TreeStatesToModelStatesMarkovian
///	Returns: void
///	Action : This function used the markovian
/// approximation to evaluate the drift
////////////////////////////////////////////////////

void ARM_SFRM::TreeStatesToModelStatesDriftMarkovian(ARM_PricingStatesPtr& states, int timeIndex) const
{
	ARM_NumerairePtr numeraire = GetNumeraire();
#if defined(__GP_STRICT_VALIDATION)
	if(		numeraire->GetType() != ARM_Numeraire::TerminalZc 
		&&  numeraire->GetType() != ARM_Numeraire::TerminalEventZc
		&&  numeraire->GetType() != ARM_Numeraire::RollingEvent
		&&  numeraire->GetType() != ARM_Numeraire::RollingPayment )
    {		
		CC_Ostringstream os;
		os  << ARM_USERNAME << ": only TerminalZc or EventZC numeraire is "
			<< " supported by the tree method in SFRM at the moment!";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );			
    }
#endif
	if ( numeraire->GetType() == ARM_Numeraire::RollingEvent ||
		 numeraire->GetType() == ARM_Numeraire::RollingPayment )
		ComputeNumeraireTimeIndex();

	int fwdMinIndex = states->GetPricingStatesContext( GetModelRank() )->ToSFRMPricingStatesContext()->GetFwdMinIndex();
	ResetDFMap();
	/// Validation tests
#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
		throw Exception(__LINE__, __FILE__,ERR_INVALID_ARGUMENT,"time index is negative!");
	
	/// test pricing direction
	if(ARM_NumMethod::GP_BCKWDLOOKING != GetNumMethod()->GetPricingDirection() )
		throw Exception(__LINE__, __FILE__,ERR_INVALID_ARGUMENT,"implemented only for backward looking method..More precisely Tree method!");
#endif

	
	//we diffuse first the forward corresponding to the numeraire zero coupon 
	int numeraireFwdIdx = itsNumeraireTimeIndex-1; 
	if (itsFwdStatus[numeraireFwdIdx] == ARM_SFRM::ARM_SFRM_AUX_FWD)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +
			": error in the schedule... the last forward has to be standard!" );
	}
	
	///	If timeIndex=0
	int stateIdx, fwdIdx; 
	size_t NbStates	= states->size();
	size_t modelNb = GetModelNb();
	int fwdMaxIndex =  itsFwdResetTimes->size()-1;
	size_t NbFwds  = fwdMaxIndex+1;

	states->resizeModelStates(ModelStatesSize(), NbStates);

	if (timeIndex==0)
	{
		for(stateIdx=0; stateIdx<NbStates; ++stateIdx )
			for(fwdIdx=NbFwds-1; fwdIdx>=fwdMinIndex; fwdIdx--)
				states->SetModelState(stateIdx,modelNb+fwdIdx,0.0);
	}
	
	else
	{
		/// Preliminaries
		/// get data
		double fromTime = GetNumMethod()->GetTimeStep(timeIndex);		
		size_t factorsNb = FactorCount();				

		ARM_GP_VectorPtr probaSpot          = GetNumMethod()->GetArrowDebreuPrices(timeIndex,*this);
		const ARM_MatrixVector& globalVars  = GetNumMethod()->GetNumMethodStateGlobalVars();
		ARM_ZeroCurvePtr zeroCurve = GetZeroCurve();
		ARM_Currency* pCcy	= zeroCurve->GetCurrencyUnit();
		int floatDayCount	= pCcy->GetLiborIndexDayCount();
		double numeraireTime = GetNumeraire()->GetMaturity();
		
		/// infered data 
		while( (fwdMinIndex > 0) && (*itsFwdResetTimes)[fwdMinIndex]> fromTime)
			--fwdMinIndex;

		states->GetPricingStatesContext( GetModelRank() )->ToSFRMPricingStatesContext()->SetFwdMinIndex( fwdMinIndex );
		
		size_t NbFwdsToDiffuse  = (NbFwds-fwdMinIndex) < 0.0 ? 0.0 : NbFwds-fwdMinIndex;
		
		double numeraireBond=zeroCurve->DiscountPrice(numeraireTime/K_YEAR_LEN);
		
		double var_At_timeIndex = 0.0;		
		for(int j=0;j<factorsNb; ++j) //variance norm at timeIdx
			var_At_timeIndex += (*globalVars[timeIndex])(j,j);
		
		size_t FollowSTDIdx = numeraireFwdIdx; // the last STD FWD diffused
		size_t FollowFollowSTDIdx = FollowSTDIdx;
		
		/// used variables
		ARM_GP_Matrix savAcc(NbStates,NbFwds);
		ARM_GP_Matrix shiftedCalibratedFWD(NbStates,NbFwds);
		std::vector<double> FwdsCorrelVarCoeffs(factorsNb);
		double correctedDrift,correctedDriftShift;
		double forwardVolCoeff, markovDrift;
		
		double shift, delta, ShiftedFWD_0, FWD, martPart;
		double nextShift, nextDelta, NextShiftedFWD_0, shifted_fwd_At_s, NextShiftedFWD, varRatio, driftPart_nextFwdIdx ;
		int nextFwdIdx,i;
		double stochasticTerm;
				
		/// Forward diffusion from fwdMinIndex to fwdMaxIndex
		/// (1)first we diffuse the forward corresponding to the  numeraire numeraireFwdIdx by deterministic calibration
		/// (2)second we diffuse all forwards before numeraireFwdIdx with the markovian drift and deterministic calibration  	    
		/// (3)third we diffuse the deterministic drift correction of all from MAX(fwdMinIndex,numeraireFwdIdx+1) to fwdMaxIndex 
		
		/// (1) states of the forward corresponding to numeraireFwdIdx/////////////////////////
		/// markovian drift =0 
		markovDrift=0.0;
		correctedDrift=0.0;
		correctedDriftShift=0.0;
		shift = (*itsFwdShifts)[numeraireFwdIdx];
		delta = (*itsFwdIntTerms)[numeraireFwdIdx];
		ShiftedFWD_0 = (*itsFwdValues)[numeraireFwdIdx]+shift;
		forwardVolCoeff = (*(const_cast<ARM_SFRM *>(this)->itsFwdsVarCoeffs))[numeraireFwdIdx];
		
		for(j=0;j<factorsNb; ++j)
		{
			FwdsCorrelVarCoeffs.Elt(j) = (*(const_cast<ARM_SFRM *>(this)->itsFwdsCorrelVarCoeffs[j]))(numeraireFwdIdx,numeraireFwdIdx);
		}
		
		for(stateIdx=0; stateIdx<NbStates; ++stateIdx )		
		{
			stochasticTerm =0.0;
			for(j=0;j<factorsNb; ++j)
			{
				double state_j = states->GetNumMethodState(stateIdx,j+modelNb);
				double FwdCorrelcoeff_j = FwdsCorrelVarCoeffs.Elt(j);
				stochasticTerm+= (state_j*FwdCorrelcoeff_j);				
			}
			
			martPart = exp(-markovDrift-0.5*forwardVolCoeff*var_At_timeIndex+stochasticTerm);
			
			shiftedCalibratedFWD.Elt(stateIdx,numeraireFwdIdx) = martPart;
			
			correctedDrift+= (*probaSpot)[stateIdx]*martPart;
		}
		
		/// deterministic drift correction
		double FwdIdxBond = zeroCurve->DiscountPrice((*itsFwdEndTimes)[numeraireFwdIdx]/K_YEAR_LEN);
		
		correctedDrift = FwdIdxBond/numeraireBond/correctedDrift;
		
		for(stateIdx=0; stateIdx<NbStates; ++stateIdx )		
		{
			martPart = shiftedCalibratedFWD.Elt(stateIdx,numeraireFwdIdx);
			FWD = ShiftedFWD_0*martPart*correctedDrift;
			shiftedCalibratedFWD.Elt(stateIdx,numeraireFwdIdx) = FWD;
			FWD = FWD-shift;
			
			states->SetModelState(stateIdx,numeraireFwdIdx+modelNb,FWD);
			
			savAcc.Elt(stateIdx,numeraireFwdIdx) = 1.0 + delta*FWD;
			
			stochasticTerm = log(martPart * correctedDrift);		
		}
		
		/// (2) states of the forwards from fwdMinIndex to numeraireFwdIdx-1
		/// Prepare to calculate the markovian drifts for only the forwards before the numeraireFwdIdx 
		int nbIntegrationSteps; 
		double *x=NULL,*w=NULL;
		std::vector<double> *VarCoeff=NULL,*SpotVarCoeff=NULL;
		ARM_GP_Matrix *FwdVolCoeff=NULL;
		
		if(fwdMinIndex < NbFwds)
		{
			nbIntegrationSteps = (int)(ceil(fromTime/K_YEAR_LEN*INTERMEDIATESTEPS));
			PrepareForMarkovDriftCalculation(timeIndex, 
				var_At_timeIndex,
				fwdMinIndex, 
				NbFwdsToDiffuse, 
				nbIntegrationSteps, 
				VarCoeff,
				SpotVarCoeff,
				FwdVolCoeff,
				x, 
				w);
		}
		
		for(fwdIdx=numeraireFwdIdx-1;fwdIdx>=fwdMinIndex;fwdIdx--)
		{
			for(j=0;j<factorsNb; ++j)
			{
				FwdsCorrelVarCoeffs.Elt(j) = (*(const_cast<ARM_SFRM *>(this)->itsFwdsCorrelVarCoeffs[j]))(fwdIdx,fwdIdx);
			}
			if (itsFwdStatus[fwdIdx] == ARM_SFRM::ARM_SFRM_STD_FWD)
			{
				/// markovian drift
				correctedDrift=0.0;
				correctedDriftShift=0.0;
				shift = (*itsFwdShifts)[fwdIdx];
				delta = (*itsFwdIntTerms)[fwdIdx];
				ShiftedFWD_0 = (*itsFwdValues)[fwdIdx]+shift;
				forwardVolCoeff = (*(const_cast<ARM_SFRM *>(this)->itsFwdsVarCoeffs))[fwdIdx];
				
				for(stateIdx=0; stateIdx<NbStates; ++stateIdx )		
				{	
					markovDrift=0.0;
					for(nextFwdIdx=numeraireFwdIdx; nextFwdIdx>fwdIdx;--nextFwdIdx)
					{
						if (itsFwdStatus[nextFwdIdx] == ARM_SFRM::ARM_SFRM_STD_FWD)
						{
							nextShift = (*itsFwdShifts)[nextFwdIdx];
							NextShiftedFWD_0 =  (*itsFwdValues)[nextFwdIdx]+nextShift;
							nextDelta = (*itsFwdIntTerms)[nextFwdIdx];
							NextShiftedFWD= shiftedCalibratedFWD.Elt(stateIdx,nextFwdIdx);
							
							driftPart_nextFwdIdx=0.0;
							for(i=0; i<nbIntegrationSteps; ++i)
							{
								varRatio = VarCoeff->Elt(i)/var_At_timeIndex;
								shifted_fwd_At_s = NextShiftedFWD_0*pow(NextShiftedFWD/NextShiftedFWD_0,varRatio)
									*FwdVolCoeff->Elt(i,nextFwdIdx-fwdMinIndex);
								
								driftPart_nextFwdIdx+=w[i]*shifted_fwd_At_s/(1+nextDelta*(shifted_fwd_At_s-nextShift))*SpotVarCoeff->Elt(i);
								
							}
							
							double FwdCorrelcoeff=0.0;
							for(j=0;j<factorsNb; ++j)
							{
								FwdCorrelcoeff+=(*(const_cast<ARM_SFRM *>(this)->itsFwdsCorrelVarCoeffs[j]))(nextFwdIdx,fwdIdx);							
							}
							
							markovDrift+=nextDelta*driftPart_nextFwdIdx*FwdCorrelcoeff;
						}
					}/// sum for the drift on only the following STD FWD 
					
					
					stochasticTerm =0.0;
					for(j=0;j<factorsNb; ++j)
					{
						double state_j = states->GetNumMethodState(stateIdx,j+modelNb);
						double FwdCorrelcoeff_j = FwdsCorrelVarCoeffs.Elt(j);
						stochasticTerm+= (state_j*FwdCorrelcoeff_j);				
					}
					
					martPart = exp(-markovDrift-0.5*forwardVolCoeff*var_At_timeIndex+stochasticTerm);
					
					shiftedCalibratedFWD.Elt(stateIdx,fwdIdx) = martPart;
					correctedDrift+= (*probaSpot)[stateIdx]*martPart* savAcc.Elt(stateIdx,FollowSTDIdx);
					
				}/// end stateIdx for markovien drift for STD FWD
				
				/// deterministic drift correction
				double FwdIdxBond = zeroCurve->DiscountPrice((*itsFwdEndTimes)[fwdIdx]/K_YEAR_LEN);
				
				correctedDrift = FwdIdxBond/numeraireBond/correctedDrift;
				
				for(stateIdx=0; stateIdx<NbStates; ++stateIdx )		
				{
					martPart = shiftedCalibratedFWD.Elt(stateIdx,fwdIdx);
					FWD = ShiftedFWD_0*martPart*correctedDrift;
					shiftedCalibratedFWD.Elt(stateIdx,fwdIdx) = FWD;
					FWD = FWD-shift;
					states->SetModelState(stateIdx,fwdIdx+modelNb,FWD);
					
					savAcc.Elt(stateIdx,fwdIdx) = (1.0 + delta*FWD)*savAcc.Elt(stateIdx, FollowSTDIdx);
					stochasticTerm = log(martPart * correctedDrift);				
				} /// end stateIdx for deterministic drift for STD FWD
				
				FollowFollowSTDIdx = FollowSTDIdx;
				FollowSTDIdx = fwdIdx; // the last STD FWD diffused
				
			} /// end STD FWD
			
			
			if (itsFwdStatus[fwdIdx] == ARM_SFRM::ARM_SFRM_AUX_FWD)
			{
				/// markovian drift
				correctedDrift=0.0;
				correctedDriftShift=0.0;
				shift = (*itsFwdShifts)[fwdIdx];
				delta = (*itsFwdIntTerms)[fwdIdx];
				ShiftedFWD_0 = (*itsFwdValues)[fwdIdx]+shift;
				forwardVolCoeff = (*(const_cast<ARM_SFRM *>(this)->itsFwdsVarCoeffs))[fwdIdx];
				
				for(stateIdx=0; stateIdx<NbStates; ++stateIdx )		
				{	
					markovDrift=0.0;
					for(nextFwdIdx=numeraireFwdIdx; nextFwdIdx>FollowSTDIdx;--nextFwdIdx)
					{
						if (itsFwdStatus[nextFwdIdx] == ARM_SFRM::ARM_SFRM_STD_FWD)
						{
							nextShift = (*itsFwdShifts)[nextFwdIdx];
							NextShiftedFWD_0 =  (*itsFwdValues)[nextFwdIdx]+nextShift;
							nextDelta = (*itsFwdIntTerms)[nextFwdIdx];
							NextShiftedFWD= shiftedCalibratedFWD.Elt(stateIdx,nextFwdIdx);
							
							driftPart_nextFwdIdx=0.0;
							for(i=0; i<nbIntegrationSteps; ++i)
							{
								varRatio = VarCoeff->Elt(i)/var_At_timeIndex;
								shifted_fwd_At_s = NextShiftedFWD_0*pow(NextShiftedFWD/NextShiftedFWD_0,varRatio)
									*FwdVolCoeff->Elt(i,nextFwdIdx-fwdMinIndex);
								
								driftPart_nextFwdIdx+=w[i]*shifted_fwd_At_s/(1+nextDelta*(shifted_fwd_At_s-nextShift))*SpotVarCoeff->Elt(i);
								
							}
							
							double FwdCorrelcoeff=0.0;
							for(j=0;j<factorsNb; ++j)
								FwdCorrelcoeff+=(*(const_cast<ARM_SFRM *>(this)->itsFwdsCorrelVarCoeffs[j]))(nextFwdIdx,fwdIdx);							
							
							markovDrift+=nextDelta*driftPart_nextFwdIdx*FwdCorrelcoeff;
						}
					}/// sum for the drift on only the following STD FWD 
					
					/// add of the term deponding on the FollowSTDIdx delta are different then the delta of (*itsFwdIntTerms)[FollowSTDIdx]
					nextShift = (*itsFwdShifts)[FollowSTDIdx];
					NextShiftedFWD_0 =  (*itsFwdValues)[FollowSTDIdx]+nextShift;
					nextDelta = CountYears(floatDayCount,(*itsFwdEndTimes)[fwdIdx],(*itsFwdEndTimes)[FollowSTDIdx]);
					NextShiftedFWD= shiftedCalibratedFWD.Elt(stateIdx,FollowSTDIdx);
					
					driftPart_nextFwdIdx=0.0;
					for(i=0; i<nbIntegrationSteps; ++i)
					{
						varRatio = VarCoeff->Elt(i)/var_At_timeIndex;
						shifted_fwd_At_s = NextShiftedFWD_0*pow(NextShiftedFWD/NextShiftedFWD_0,varRatio)
							*FwdVolCoeff->Elt(i,FollowSTDIdx-fwdMinIndex);
						
						driftPart_nextFwdIdx+=w[i]*shifted_fwd_At_s/(1+nextDelta*(shifted_fwd_At_s-nextShift))*SpotVarCoeff->Elt(i);
						
					}
					
					double FwdCorrelcoeff=0.0;
					for(j=0;j<factorsNb; ++j)
					{
						FwdCorrelcoeff+=(*(const_cast<ARM_SFRM *>(this)->itsFwdsCorrelVarCoeffs[j]))(FollowSTDIdx,fwdIdx);							
					}
					
					markovDrift+=nextDelta*driftPart_nextFwdIdx*FwdCorrelcoeff;
					/// end of the add
					
					stochasticTerm =0.0;
					for(j=0;j<factorsNb; ++j)
					{
						double state_j = states->GetNumMethodState(stateIdx,j+modelNb);
						double FwdCorrelcoeff_j = FwdsCorrelVarCoeffs.Elt(j);
						stochasticTerm+= (state_j*FwdCorrelcoeff_j);				
					}
					
					martPart = exp(-markovDrift-0.5*forwardVolCoeff*var_At_timeIndex+stochasticTerm);
					
					shiftedCalibratedFWD.Elt(stateIdx,fwdIdx) = martPart;					
					double stubDFStart	  = zeroCurve->DiscountPrice((*itsFwdEndTimes)[fwdIdx]/K_YEAR_LEN);
					double stubDFEnd	  = zeroCurve->DiscountPrice((*itsFwdEndTimes)[FollowSTDIdx]/K_YEAR_LEN);
					double stubFwd = 0.0;
					
					if (fabs(nextDelta) > K_NEW_DOUBLE_TOL)
						stubFwd = (stubDFStart/stubDFEnd-1.0)/nextDelta;
					
					double stubSavAcc = (states->GetModelState(stateIdx,FollowSTDIdx+modelNb)+nextShift)/NextShiftedFWD_0;
					stubSavAcc = 1.0+nextDelta*(stubSavAcc*(stubFwd+nextShift)-nextShift);
					stubSavAcc = (FollowSTDIdx==numeraireFwdIdx) ? stubSavAcc : savAcc.Elt(stateIdx,FollowFollowSTDIdx)*stubSavAcc;
					correctedDrift+= (*probaSpot)[stateIdx]*martPart* stubSavAcc;
				}
				
				/// deterministic drift correction
				double FwdIdxBond = zeroCurve->DiscountPrice((*itsFwdEndTimes)[fwdIdx]/K_YEAR_LEN);
				
				correctedDrift = FwdIdxBond/numeraireBond/correctedDrift;
				
				for(stateIdx=0; stateIdx<NbStates; ++stateIdx )		
				{
					martPart = shiftedCalibratedFWD.Elt(stateIdx,fwdIdx);
					FWD = ShiftedFWD_0*martPart*correctedDrift;
					shiftedCalibratedFWD.Elt(stateIdx,fwdIdx) = FWD;
					FWD = FWD-shift;
					states->SetModelState(stateIdx,fwdIdx+modelNb,FWD);
					stochasticTerm = log(martPart * correctedDrift);
				}
				
			}
		}
		
///AFTERNUMINDEX//////////////////////////////////
		FollowSTDIdx = numeraireFwdIdx; // the last STD FWD diffused
		FollowFollowSTDIdx = FollowSTDIdx;
		for(stateIdx=0; stateIdx<NbStates; ++stateIdx )		
		{
			savAcc.Elt(stateIdx,FollowFollowSTDIdx) = 1.0;
		
		}
		for(fwdIdx=numeraireFwdIdx+1;fwdIdx<NbFwds;fwdIdx++)
		{
			for(j=0;j<factorsNb; ++j)
			{
				FwdsCorrelVarCoeffs.Elt(j) = (*(const_cast<ARM_SFRM *>(this)->itsFwdsCorrelVarCoeffs[j]))(fwdIdx,fwdIdx);
			}
			if (itsFwdStatus[fwdIdx] == ARM_SFRM::ARM_SFRM_STD_FWD)
			{
				FollowFollowSTDIdx = FollowSTDIdx;
				FollowSTDIdx = fwdIdx; // the last STD FWD diffused
				/// markovian drift
				correctedDrift=0.0;
				correctedDriftShift=0.0;
				shift = (*itsFwdShifts)[fwdIdx];
				delta = (*itsFwdIntTerms)[fwdIdx];
				ShiftedFWD_0 = (*itsFwdValues)[fwdIdx]+shift;
				forwardVolCoeff = (*(const_cast<ARM_SFRM *>(this)->itsFwdsVarCoeffs))[fwdIdx];
				
				for(stateIdx=0; stateIdx<NbStates; ++stateIdx )		
				{	
					markovDrift=0.0;
					for(nextFwdIdx=numeraireFwdIdx+1; nextFwdIdx<fwdIdx;++nextFwdIdx)
					{
						if (itsFwdStatus[nextFwdIdx] == ARM_SFRM::ARM_SFRM_STD_FWD)
						{
							nextShift = (*itsFwdShifts)[nextFwdIdx];
							NextShiftedFWD_0 =  (*itsFwdValues)[nextFwdIdx]+nextShift;
							nextDelta = (*itsFwdIntTerms)[nextFwdIdx];
							NextShiftedFWD= shiftedCalibratedFWD.Elt(stateIdx,nextFwdIdx);
							
							driftPart_nextFwdIdx=0.0;
							for(i=0; i<nbIntegrationSteps; ++i)
							{
								varRatio = VarCoeff->Elt(i)/var_At_timeIndex;
								shifted_fwd_At_s = NextShiftedFWD_0*pow(NextShiftedFWD/NextShiftedFWD_0,varRatio)
									*FwdVolCoeff->Elt(i,nextFwdIdx-fwdMinIndex);
								
								driftPart_nextFwdIdx+=w[i]*shifted_fwd_At_s/(1+nextDelta*(shifted_fwd_At_s-nextShift))*SpotVarCoeff->Elt(i);
								
							}
							
							double FwdCorrelcoeff=0.0;
							for(j=0;j<factorsNb; ++j)
							{
								FwdCorrelcoeff+=(*(const_cast<ARM_SFRM *>(this)->itsFwdsCorrelVarCoeffs[j]))(nextFwdIdx,fwdIdx);							
							}
							
							markovDrift+=nextDelta*driftPart_nextFwdIdx*FwdCorrelcoeff;
						}
					}/// sum for the drift on only the following STD FWD 

					//plus the coeff of fwdIdx
					nextShift = (*itsFwdShifts)[FollowFollowSTDIdx];
					NextShiftedFWD_0 =  (*itsFwdValues)[FollowFollowSTDIdx]+nextShift;
					NextShiftedFWD = shiftedCalibratedFWD.Elt(stateIdx,FollowFollowSTDIdx);
					double shiftedFWD = (NextShiftedFWD/NextShiftedFWD_0)*ShiftedFWD_0;
					
					driftPart_nextFwdIdx=0.0;
					for(i=0; i<nbIntegrationSteps; ++i)
					{
						varRatio = VarCoeff->Elt(i)/var_At_timeIndex;
						shifted_fwd_At_s = ShiftedFWD_0*pow(shiftedFWD/ShiftedFWD_0,varRatio)
							*FwdVolCoeff->Elt(i,fwdIdx-fwdMinIndex);
						
						driftPart_nextFwdIdx+=w[i]*shifted_fwd_At_s/(1+delta*(shifted_fwd_At_s-shift))*SpotVarCoeff->Elt(i);
						
					}
					
					double FwdCorrelcoeff=0.0;
					for(j=0;j<factorsNb; ++j)
					{
						FwdCorrelcoeff+=(*(const_cast<ARM_SFRM *>(this)->itsFwdsCorrelVarCoeffs[j]))(fwdIdx,fwdIdx);							
					}
					
					markovDrift+=delta*driftPart_nextFwdIdx*FwdCorrelcoeff;///////////
			
					
					stochasticTerm =0.0;
					for(j=0;j<factorsNb; ++j)
					{
						double state_j = states->GetNumMethodState(stateIdx,j+modelNb);
						double FwdCorrelcoeff_j = FwdsCorrelVarCoeffs.Elt(j);
						stochasticTerm+= (state_j*FwdCorrelcoeff_j);				
					}
					
					martPart = exp(markovDrift-0.5*forwardVolCoeff*var_At_timeIndex+stochasticTerm);
					savAcc.Elt(stateIdx,FollowSTDIdx) = savAcc.Elt(stateIdx,FollowFollowSTDIdx)*(1.0 + delta*(ShiftedFWD_0*martPart-shift));
					shiftedCalibratedFWD.Elt(stateIdx,fwdIdx) = martPart;
					correctedDrift+= (*probaSpot)[stateIdx]*martPart/savAcc.Elt(stateIdx,FollowSTDIdx);
					
				}/// end stateIdx for markovien drift for STD FWD
				
				/// deterministic drift correction
				double FwdIdxBond = zeroCurve->DiscountPrice((*itsFwdEndTimes)[fwdIdx]/K_YEAR_LEN);
				
				correctedDrift = FwdIdxBond/numeraireBond/correctedDrift;
				
				for(stateIdx=0; stateIdx<NbStates; ++stateIdx )		
				{
					martPart = shiftedCalibratedFWD.Elt(stateIdx,fwdIdx);
					FWD = ShiftedFWD_0*martPart*correctedDrift;
					shiftedCalibratedFWD.Elt(stateIdx,fwdIdx) = FWD;
					FWD = FWD-shift;
					states->SetModelState(stateIdx,fwdIdx+modelNb,FWD);
					
					savAcc.Elt(stateIdx,fwdIdx) = (1.0 + delta*FWD)*savAcc.Elt(stateIdx, FollowFollowSTDIdx);
					stochasticTerm = log(martPart * correctedDrift);				
				} /// end stateIdx for deterministic drift for STD FWD
				
				/*FollowFollowSTDIdx = FollowSTDIdx;
				FollowSTDIdx = fwdIdx; // the last STD FWD diffused*/
				
			} /// end STD FWD
			
			
			if (itsFwdStatus[fwdIdx] == ARM_SFRM::ARM_SFRM_AUX_FWD)
			{
				/// markovian drift
				correctedDrift=0.0;
				correctedDriftShift=0.0;
				shift = (*itsFwdShifts)[fwdIdx];
				delta = (*itsFwdIntTerms)[fwdIdx];
				ShiftedFWD_0 = (*itsFwdValues)[fwdIdx]+shift;
				forwardVolCoeff = (*(const_cast<ARM_SFRM *>(this)->itsFwdsVarCoeffs))[fwdIdx];
				
				for(stateIdx=0; stateIdx<NbStates; ++stateIdx )		
				{	
					markovDrift=0.0;					
					for(nextFwdIdx=numeraireFwdIdx+1; nextFwdIdx<=FollowSTDIdx;++nextFwdIdx)
					{
						if (itsFwdStatus[nextFwdIdx] == ARM_SFRM::ARM_SFRM_STD_FWD)
						{
							nextShift = (*itsFwdShifts)[nextFwdIdx];
							NextShiftedFWD_0 =  (*itsFwdValues)[nextFwdIdx]+nextShift;
							nextDelta = (*itsFwdIntTerms)[nextFwdIdx];
							NextShiftedFWD= shiftedCalibratedFWD.Elt(stateIdx,nextFwdIdx);
							
							driftPart_nextFwdIdx=0.0;
							for(i=0; i<nbIntegrationSteps; ++i)
							{
								varRatio = VarCoeff->Elt(i)/var_At_timeIndex;
								shifted_fwd_At_s = NextShiftedFWD_0*pow(NextShiftedFWD/NextShiftedFWD_0,varRatio)
									*FwdVolCoeff->Elt(i,nextFwdIdx-fwdMinIndex);
								
								driftPart_nextFwdIdx+=w[i]*shifted_fwd_At_s/(1+nextDelta*(shifted_fwd_At_s-nextShift))*SpotVarCoeff->Elt(i);
								
							}
							
							double FwdCorrelcoeff=0.0;
							for(j=0;j<factorsNb; ++j)
								FwdCorrelcoeff+=(*(const_cast<ARM_SFRM *>(this)->itsFwdsCorrelVarCoeffs[j]))(nextFwdIdx,fwdIdx);							
							
							markovDrift+=nextDelta*driftPart_nextFwdIdx*FwdCorrelcoeff;
						}
					}/// sum for the drift on only the following STD FWD 
					
					/// add of the term deponding on the FollowSTDIdx delta are different then the delta of (*itsFwdIntTerms)[FollowSTDIdx]
					nextShift = (*itsFwdShifts)[FollowSTDIdx];
					NextShiftedFWD_0 =  (*itsFwdValues)[FollowSTDIdx]+nextShift;
					nextDelta = CountYears(floatDayCount,(*itsFwdEndTimes)[fwdIdx],(*itsFwdEndTimes)[FollowSTDIdx]);
					NextShiftedFWD= shiftedCalibratedFWD.Elt(stateIdx,FollowSTDIdx);
					
					driftPart_nextFwdIdx=0.0;
					for(i=0; i<nbIntegrationSteps; ++i)
					{
						varRatio = VarCoeff->Elt(i)/var_At_timeIndex;
						shifted_fwd_At_s = NextShiftedFWD_0*pow(NextShiftedFWD/NextShiftedFWD_0,varRatio)
							*FwdVolCoeff->Elt(i,FollowSTDIdx-fwdMinIndex);
						
						driftPart_nextFwdIdx+=w[i]*shifted_fwd_At_s/(1+nextDelta*(shifted_fwd_At_s-nextShift))*SpotVarCoeff->Elt(i);
						
					}
					
					double FwdCorrelcoeff=0.0;
					for(j=0;j<factorsNb; ++j)
					{
						FwdCorrelcoeff+=(*(const_cast<ARM_SFRM *>(this)->itsFwdsCorrelVarCoeffs[j]))(FollowSTDIdx,fwdIdx);							
					}
					
					markovDrift+=nextDelta*driftPart_nextFwdIdx*FwdCorrelcoeff;
					/// end of the add
					
					stochasticTerm =0.0;
					for(j=0;j<factorsNb; ++j)
					{
						double state_j = states->GetNumMethodState(stateIdx,j+modelNb);
						double FwdCorrelcoeff_j = FwdsCorrelVarCoeffs.Elt(j);
						stochasticTerm+= (state_j*FwdCorrelcoeff_j);				
					}
					
					martPart = exp(markovDrift-0.5*forwardVolCoeff*var_At_timeIndex+stochasticTerm);
					
					shiftedCalibratedFWD.Elt(stateIdx,fwdIdx) = martPart;					
					double stubDFStart	  = zeroCurve->DiscountPrice((*itsFwdEndTimes)[fwdIdx]/K_YEAR_LEN);
					double stubDFEnd	  = zeroCurve->DiscountPrice((*itsFwdEndTimes)[FollowSTDIdx]/K_YEAR_LEN);
					double stubFwd = 0.0;
					
					if (fabs(nextDelta) > K_NEW_DOUBLE_TOL)
						stubFwd = (stubDFStart/stubDFEnd-1.0)/nextDelta;
					
					double stubSavAcc = (states->GetModelState(stateIdx,FollowSTDIdx+modelNb)+nextShift)/NextShiftedFWD_0;
					stubSavAcc = 1.0+nextDelta*(stubSavAcc*(stubFwd+nextShift)-nextShift);
					stubSavAcc = (FollowSTDIdx==numeraireFwdIdx) ? stubSavAcc : savAcc.Elt(stateIdx,FollowFollowSTDIdx)*stubSavAcc;
					correctedDrift+= (*probaSpot)[stateIdx]*martPart/stubSavAcc;
				}
				
				/// deterministic drift correction
				double FwdIdxBond = zeroCurve->DiscountPrice((*itsFwdEndTimes)[fwdIdx]/K_YEAR_LEN);
				
				correctedDrift = FwdIdxBond/numeraireBond/correctedDrift;
				
				for(stateIdx=0; stateIdx<NbStates; ++stateIdx )		
				{
					martPart = shiftedCalibratedFWD.Elt(stateIdx,fwdIdx);
					FWD = ShiftedFWD_0*martPart*correctedDrift;
					shiftedCalibratedFWD.Elt(stateIdx,fwdIdx) = FWD;
					FWD = FWD-shift;
					states->SetModelState(stateIdx,fwdIdx+modelNb,FWD);
					stochasticTerm = log(martPart * correctedDrift);
				}
				
			}
		}
		
		delete [] x;
		delete [] w;
		delete VarCoeff;
		delete SpotVarCoeff;
		delete FwdVolCoeff;
		/// (3) states of the forwards from fwdMinIndex to numeraireFwdIdx-1
	}
}
//::tree


////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: ConvertToBetaOrShiftParam
///	Returns : void
///	Action  : If the model params are in beta space
///			change it to shit space
////////////////////////////////////////////////////
void ARM_SFRM::ConvertToShiftorBetaParam(const ARM_Portfolio& shiftConvPort)
{
    if(((ARM_ModelParamsSFRM*) GetModelParams())->DoesModelParamExist(ARM_ModelParamType::Beta))
        ((ARM_ModelParamsSFRM*) GetModelParams())->ConvertToShiftParam(*this,shiftConvPort);
    else if(((ARM_ModelParamsSFRM*) GetModelParams())->DoesModelParamExist(ARM_ModelParamType::Shift))
        ((ARM_ModelParamsSFRM*) GetModelParams())->ConvertToBetaParam(*this,shiftConvPort);
    else
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME +
                ": an SFRM Model should contain a shift or a beta (but not both)!" );
}

////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: ModelFixTimeStep
///	Returns : modifies if necessary the fixTimeStep!
///	Action  : divides the characteristic period for the model
///				by the fixTimeStep!
///////////////////////////////////////////////////

int ARM_SFRM::ModelFixTimeStep( int fixTimeStep ) const
{
	ARM_IRIndex* pIndex	= ((ARM_ModelParamsSFRM*)GetModelParams())->GetIRIndex();
	int floatResetFreq	= pIndex->GetResetFrequency();
	return K_YEAR_LEN/floatResetFreq/(fixTimeStep+1.0);
}

////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: VarianceToTime
///	Returns : void
///	Action  : No default implementation since inappropriate for SFRM  model
////////////////////////////////////////////////////
double ARM_SFRM::VarianceToTime(double var,double minTime,double maxTime) const
{
    return ((ARM_ModelParamsSFRM*)GetModelParams())->VarianceToTime(var,minTime,maxTime);
}

////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: ComputeUnderlyingVarCovar
///	Returns : void
///	Action  :  compute the correlation between two underlying
////////////////////////////////////////////////////

void ARM_SFRM::ComputeUnderlyingVarCovar( double fromTime, double toTime, 
	double startTime1, double endTime1, double startTime2,
	double endTime2, std::vector<double>& VarCovar, std::vector<double>& swapfwd) const
{
	ARM_ModelParams* modParamsSFRM = const_cast<ARM_SFRM* const>(this)->GetModelParams();
	double asOfDate	= GetAsOfDate().GetJulian();
	ARM_ZeroCurvePtr ZcCurve	= GetZeroCurve();
	ARM_Currency* Ccy			= ZcCurve->GetCurrencyUnit();
	string curveCcy				= Ccy->GetCcyName();
	
	int fixFreq					= Ccy->GetFixedPayFreq();
	long fixDayCount			= Ccy->GetFixedDayCount();
	char fixCalendar[100];
	Ccy->CalcFixPayCal(fixCalendar);
	int fwdDayCount				= Ccy->GetLiborIndexDayCount();

	double evalTime				= 0.0;
	double swapNotional			= 100.0;
	int callPut					= 1;
	std::vector<double> fixNotional (1,0.0);
	std::vector<double> floatNotional (1,0.0);
	ARM_GP_Matrix strikesPerState(1,0.0);
	size_t factorsNb			= FactorCount();
	ARM_PricingStatesPtr states = ARM_PricingStatesPtr( new ARM_PricingStates(1,1,1,factorsNb) );

	//CMS1
	ARM_Date startDate1(asOfDate+startTime1);
	ARM_Date endDate1(asOfDate+endTime1);
	CC_NS(std,auto_ptr)<ARM_DateStrip> dateStrip1( GetFloatDateStrip(startDate1,endDate1) );
	std::vector<double>& floatResetDate1 = dateStrip1->GetResetDates();
	double swapResetTime1 = (*floatResetDate1)[0]-asOfDate;	
	int size1 = floatResetDate1->size();
	std::vector<double>& floatStartTimes1		= dateStrip1->GetFlowStartDates();
	std::vector<double>& floatPayTimes1		= dateStrip1->GetFlowEndDates();
	std::vector<double>& floatIntTerms1		= dateStrip1->GetInterestTerms();


	ARM_DateStrip fixDateStrip1( startDate1, endDate1, fixFreq, fixDayCount, fixCalendar,
		K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,
		fixCalendar );

	std::vector<double>& fixPayTimes1   = fixDateStrip1.GetPaymentDates();
	std::vector<double>& fixPayPeriods1 = fixDateStrip1.GetInterestTerms();
	int i,nbFixFlows1=fixPayTimes1->size();    
	
	for(i=0;i<nbFixFlows1;++i)
		(*fixPayTimes1)[i] = (*fixPayTimes1)[i]-asOfDate;
	
	for(i=0;i<size1;++i)
	{
		(*floatResetDate1)[i] = (*floatResetDate1)[i]-asOfDate;
		(*floatStartTimes1)[i] = (*floatStartTimes1)[i]-asOfDate;
		(*floatPayTimes1)[i] = (*floatPayTimes1)[i]-asOfDate;
	}
	
	ARM_VanillaSwaptionArgSFRM* swaptionArg1;
	swaptionArg1 = GetVanillaSwaptionArg(curveCcy, evalTime, swapResetTime1,
		swapNotional,fixNotional,floatNotional, startTime1, endTime1,*floatResetDate1,*floatStartTimes1,*floatPayTimes1,*floatIntTerms1, *fixPayTimes1, *fixPayPeriods1, strikesPerState,
		callPut, states);

	ARM_GP_VectorPtr mu1 = swaptionArg1->GetMu();
	double var1 = ((ARM_ModelParamsSFRM*)GetModelParams() )->LocalVolatity(fromTime,toTime,*swaptionArg1,mu1);
	var1 = var1 * sqrt((toTime-fromTime)/K_YEAR_LEN);

	ARM_VectorPtr swapfwd1	= swaptionArg1->GetSwapFwd();
	double shift1 = swaptionArg1->GetAverageShift();

	//CMS2		
	ARM_Date startDate2(asOfDate+startTime2);
	ARM_Date endDate2(asOfDate+endTime2);
	CC_NS(std,auto_ptr)<ARM_DateStrip> dateStrip2( GetFloatDateStrip(startDate2,endDate2) );
	std::vector<double>& floatResetDate2 = dateStrip2->GetResetDates();
	double swapResetTime2 = (*floatResetDate2)[0]-asOfDate;
	int size2 = floatResetDate2->size();
	std::vector<double>& floatStartTimes2		= dateStrip2->GetFlowStartDates();
	std::vector<double>& floatPayTimes2		= dateStrip2->GetFlowEndDates();
	std::vector<double>& floatIntTerms2		= dateStrip2->GetInterestTerms();



	ARM_DateStrip fixDateStrip2( startDate2, endDate2, fixFreq, fixDayCount, fixCalendar,
		K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,
		fixCalendar );

	std::vector<double>& fixPayTimes2 = fixDateStrip2.GetPaymentDates();
	std::vector<double>& fixPayPeriods2 = fixDateStrip2.GetInterestTerms();
	int nbFixFlows2=fixPayTimes2->size();    
	
	for(i=0;i<nbFixFlows2;++i)
		(*fixPayTimes2)[i] = (*fixPayTimes2)[i]-asOfDate;
	
	for(i=0;i<size2;++i)
	{
		(*floatResetDate2)[i] = (*floatResetDate2)[i]-asOfDate;
		(*floatStartTimes2)[i] = (*floatStartTimes2)[i]-asOfDate;
		(*floatPayTimes2)[i] = (*floatPayTimes2)[i]-asOfDate;
	}

	ARM_VanillaSwaptionArgSFRM* swaptionArg2;
	swaptionArg2 = GetVanillaSwaptionArg(curveCcy, evalTime, swapResetTime2,swapNotional,fixNotional,floatNotional,
		startTime2,endTime2,*floatResetDate2,*floatStartTimes2,*floatPayTimes2,*floatIntTerms2,*fixPayTimes2,*fixPayPeriods2,strikesPerState,
		callPut, states);

	ARM_GP_VectorPtr mu2 = swaptionArg2->GetMu();
	double var2 = ((ARM_ModelParamsSFRM*)GetModelParams() )->LocalVolatity(fromTime,toTime,*swaptionArg2,mu2);

	var2 = var2 * sqrt((toTime-fromTime)/K_YEAR_LEN);
	
	ARM_VectorPtr swapfwd2	= swaptionArg2->GetSwapFwd();
	double shift2 = swaptionArg2->GetAverageShift();
	
	//covariance_1_2
	double covariance_1_2 =0.0;
	double covariance_1_2prime2 =0.0;

	for(int k=0;k<size1;++k)
	{			
		double mu1_k = (*mu1)[k];
		double resetTime1 = (*floatResetDate1)[k];
		
		for(int j=0;j<size2;++j)
		{
			double mu2_j = (*mu2)[j];
			double resetTime2 = (*floatResetDate2)[j];
			double FWDCovar_k_j = dynamic_cast<ARM_ModelParamsSFRM *>(modParamsSFRM)->IntegratedCovariance(fromTime,toTime,resetTime1,resetTime2);
			double coeff_k_j = (mu1_k*mu2_j);
			covariance_1_2 += coeff_k_j*FWDCovar_k_j;	
		}
	}

	VarCovar[0] = var1;
	VarCovar[1] = var2;
	VarCovar[2] = covariance_1_2;

	swapfwd.resize(4);
	swapfwd[0]=(*swapfwd1)[0];
	swapfwd[1]=(*swapfwd2)[0];
	swapfwd[2]=shift1;
	swapfwd[3]=shift2;

	delete swaptionArg1;
	delete swaptionArg2;
}

////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: ComputeSwaptionVarCovar
///	Returns : void
///	Action  :  compute the correlation between two underlying
////////////////////////////////////////////////////

void ARM_SFRM::ComputeSwaptionVarCovar(double	fromTime,
									   double	toTime,
									   const ARM_VanillaSwaptionArgSFRM& swaptionArg1,
									   const ARM_VanillaSwaptionArgSFRM& swaptionArg2,
									   std::vector<double>& VarCovar,
									   std::vector<double>& swapfwd) const
{
	//long Swaption
	ARM_ModelParams* modParamsSFRM = const_cast<ARM_SFRM* const>(this)->GetModelParams();
	ARM_GP_VectorPtr mu1 = swaptionArg1.GetMu();
	double var1 = ((ARM_ModelParamsSFRM*)modParamsSFRM)->LocalVolatity(fromTime,toTime,swaptionArg1,mu1);
	var1 = var1*sqrt((toTime-fromTime)/K_YEAR_LEN);
	std::vector<double>& floatResetDate1 = const_cast<std::vector<double>& const>(swaptionArg1.GetFloatResetTimes());
	int size1 = floatResetDate1->size();
	ARM_VectorPtr swapfwd1	= swaptionArg1.GetSwapFwd();
	double shift1 = swaptionArg1.GetAverageShift();
	

	//shortSwaption
	ARM_GP_VectorPtr mu2 = swaptionArg2.GetMu();
	double var2 = ((ARM_ModelParamsSFRM*)modParamsSFRM)->LocalVolatity(fromTime,toTime,swaptionArg2,mu2);
	var2 = var2*sqrt((toTime-fromTime)/K_YEAR_LEN);
	std::vector<double>& floatResetDate2 = const_cast<std::vector<double>& const>(swaptionArg2.GetFloatResetTimes());
	int size2 = floatResetDate2->size();
	ARM_VectorPtr swapfwd2	= swaptionArg2.GetSwapFwd();
	double shift2 = swaptionArg2.GetAverageShift();
	
	//covariance_longShort_2
	double covariance_1_2 =0.0;
	
	for(int k=0;k<size1;++k)
	{			
		double mu1_k = (*mu1)[k];
		double resetTime1 = (*floatResetDate1)[k];
		
		for(int j=0;j<size2;++j)
		{
			double mu2_j = (*mu2)[j];
			double resetTime2 = (*floatResetDate2)[j];
			double FWDCovar_k_j = dynamic_cast<ARM_ModelParamsSFRM *>(modParamsSFRM)->IntegratedCovariance(fromTime,toTime,resetTime1,resetTime2);
			double coeff_k_j = (mu1_k*mu2_j);
			covariance_1_2 += coeff_k_j*FWDCovar_k_j;	
		}
	}

	

	VarCovar[0] = var1;
	VarCovar[1] = var2;
	VarCovar[2] = covariance_1_2;

	swapfwd.resize(4);
	swapfwd[0]=(*swapfwd1)[0];
	swapfwd[1]=(*swapfwd2)[0];
	swapfwd[2]=shift1;
	swapfwd[3]=shift2;
}


////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: UnderlyingCorrelation
///	Returns : void
///	Action  :  compute the correlation between two underlying
////////////////////////////////////////////////////

double ARM_SFRM::UnderlyingCorrelation(	string	underlyingType,
	double fromTime, double toTime, double	startTime1, double endTime1,
	double startTime2, double endTime2, double startTime3, double endTime3, double startTime4, double endTime4) const
{

	if ( stringGetUpper(underlyingType)!="FWD" && stringGetUpper(underlyingType)!= "CMS" )
		ARM_THROW( ERR_INVALID_ARGUMENT, "Only FWD or CMS type supported!");

	std::vector<double> VarCovar(3);
	std::vector<double> swapfwd(4);
	
	ComputeUnderlyingVarCovar( fromTime,toTime,startTime1,endTime1,startTime2,endTime2,VarCovar,
		swapfwd);
	double stdev1 = VarCovar[0];
	double stdev2 = VarCovar[1];
	double covariance_1_2 = VarCovar[2]; 
	double correlation_1_2 = covariance_1_2/(stdev1*stdev2);
	return  correlation_1_2 ;	
}


////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: UnderlyingCovariance
///	Returns : void
///	Action  : compute the covariance between two underlying
////////////////////////////////////////////////////

double ARM_SFRM::UnderlyingCovariance( string	underlyingType,
	double fromTime, double toTime, double startTime1, double endTime1,
	double startTime2, double endTime2, double startTime3, double endTime3, double startTime4, double endTime4) const
{
	if ( stringGetUpper(underlyingType)!="FWD" && stringGetUpper(underlyingType)!= "CMS" )
		ARM_THROW( ERR_INVALID_ARGUMENT, "Only FWD or CMS type supported!");

	std::vector<double> VarCovar(3);
	std::vector<double> swapfwd(4);
	ComputeUnderlyingVarCovar( fromTime,toTime,startTime1,endTime1,startTime2,endTime2,VarCovar,
		swapfwd);

	return VarCovar[2];

	/* -------------------------------------------------------------------------
	/// the code below computes the covariance in terms of gaussian dynamics
	/// it takes into accounts the shift parameter on each of the two rates
	/// this can be useful to compute the covariance between a libor index and a spread
	/// leave it as comment not to impact existing spreadsheets...
	double lognormcov = VarCovar[2]; 
	double stdev1 = VarCovar[0];
	double stdev2 = VarCovar[1];
	double correl = lognormcov / (stdev1 * stdev2);

	double fwd1   = swapfwd[0];
	double fwd2   = swapfwd[1];
	double shift1 = swapfwd[2];
	double shift2 = swapfwd[3];

	double covar =		fwd1   * fwd2 
					+	fwd1   * shift2
					+	shift1 * fwd2
					+	shift1 * shift2 ;

	covar *= correl * stdev1 * stdev2 ;
	
	return covar;
	------------------------------------------------------------------------- */
}


////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: VanillaSpreadOption
///	Returns : void
///	Action  : No default implementation
////////////////////////////////////////////////////
ARM_VectorPtr  ARM_SFRM::VanillaSpreadOptionLet(const string& curveName,
												double evalTime,
												int callPut,
												double startTime,
												double endTime,
												double resetTime,
												double payTime,
												double payPeriod,
												double notional,
												double coeffLong,
												double coeffShort,
												const std::vector<double>& strikes,
												double swapLongFloatStartTime,
												double swapLongFloatEndTime,
												const std::vector<double>& swapLongFixPayTimes,
												const std::vector<double>& swapLongFixPayPeriods,
												double swapShortFloatStartTime,
												double swapShortFloatEndTime,
												const std::vector<double>& swapShortFixPayTimes,
												const std::vector<double>& swapShortFixPayPeriods,
												const ARM_PricingStatesPtr& states) const
{
	if( states == ARM_PricingStatesPtr(NULL) )
		return static_cast<ARM_VectorPtr>(new std::vector<double>(1,0.0));
	
    if( !GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
    {
        /// We need to compute the floating leg by forward method and no more by double notional
        /// but we have not at the moment all floating leg datas => throw an error
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "spreadoption pricing not implemented for different discount & fixing curves" );
    }
	
	if(itsUpdateVolSwapvolFRA)
    {
		ARM_VanillaSwaptionArgSFRM* swaptionArg_Long =  ((ARM_VanillaSpreadOptionletArgSFRM*)itsCurrentArg)->GetSwaptionArg_Long();
		ARM_VanillaSwaptionArgSFRM* swaptionArg_Short = ((ARM_VanillaSpreadOptionletArgSFRM*)itsCurrentArg)->GetSwaptionArg_Short();
		ARM_VanillaSwaptionArgSFRM* swaptionArg_ToTimePayTime = ((ARM_VanillaSpreadOptionletArgSFRM*)itsCurrentArg)->GetSwaptionArg_ToTimePayTime();
		ARM_VanillaSwaptionArgSFRM* swaptionArg_StartTimePayTime = ((ARM_VanillaSpreadOptionletArgSFRM*)itsCurrentArg)->GetSwaptionArg_StartTimePayTime();

        /// computes the averageShift 
	    ((ARM_VanillaSwaptionArgSFRM*)swaptionArg_Long)->SetAverageShift( 
			((ARM_ModelParamsSFRM*) GetModelParams())->AverageShift(*((ARM_VanillaSwaptionArgSFRM*)swaptionArg_Long)->GetFloatResetTimes()) );
        
		((ARM_VanillaSwaptionArgSFRM*)swaptionArg_Short)->SetAverageShift( 
			((ARM_ModelParamsSFRM*) GetModelParams())->AverageShift(*((ARM_VanillaSwaptionArgSFRM*)swaptionArg_Short)->GetFloatResetTimes()) );

		((ARM_VanillaSwaptionArgSFRM*)swaptionArg_ToTimePayTime)->SetAverageShift( 
			((ARM_ModelParamsSFRM*) GetModelParams())->AverageShift(*((ARM_VanillaSwaptionArgSFRM*)swaptionArg_ToTimePayTime)->GetFloatResetTimes()) );

		((ARM_VanillaSwaptionArgSFRM*)swaptionArg_StartTimePayTime)->SetAverageShift( 
			((ARM_ModelParamsSFRM*) GetModelParams())->AverageShift(*((ARM_VanillaSwaptionArgSFRM*)swaptionArg_StartTimePayTime)->GetFloatResetTimes()) );
		
		/// Computes  Mu
        ((ARM_VanillaSwaptionArgSFRM*)swaptionArg_Long)->SetMu( 
			((ARM_ModelParamsSFRM*) GetModelParams())->ComputeVolSwapvolFRA(*(ARM_VanillaSwaptionArgSFRM*)swaptionArg_Long,*this) );

		((ARM_VanillaSwaptionArgSFRM*)swaptionArg_Short)->SetMu( 
			((ARM_ModelParamsSFRM*) GetModelParams())->ComputeVolSwapvolFRA(*(ARM_VanillaSwaptionArgSFRM*)swaptionArg_Short,*this) );

		((ARM_VanillaSwaptionArgSFRM*)swaptionArg_ToTimePayTime)->SetMu( 
			((ARM_ModelParamsSFRM*) GetModelParams())->ComputeVolSwapvolFRA(*(ARM_VanillaSwaptionArgSFRM*)swaptionArg_ToTimePayTime,*this) );

		((ARM_VanillaSwaptionArgSFRM*)swaptionArg_StartTimePayTime)->SetMu( 
			((ARM_ModelParamsSFRM*) GetModelParams())->ComputeVolSwapvolFRA(*(ARM_VanillaSwaptionArgSFRM*)swaptionArg_StartTimePayTime,*this) );

			
    }

	ARM_VanillaSpreadOptionletArgSFRM* arg;
	// This flag check if arg is a local variable
	bool isArgLocal = false;
	if( NULL == itsCurrentArg )
	{
		arg =GetVanillaSpreadOptionletArg( 
						curveName,
						evalTime,
						callPut,
						startTime,
						endTime,
						resetTime,
						payTime,
						payPeriod,
						notional,
						coeffLong,
						coeffShort,
						strikes,
						swapLongFloatStartTime,
						swapLongFloatEndTime,
						swapLongFixPayTimes,
						swapLongFixPayPeriods,
						swapShortFloatStartTime,
						swapShortFloatEndTime,
						swapShortFixPayTimes,
						swapShortFixPayPeriods);
		isArgLocal = true;
	}
	else
	{
		arg = (ARM_VanillaSpreadOptionletArgSFRM*) &*itsCurrentArg;
	}

	ARM_VanillaSwaptionArgSFRM* swaptionArg_Long =  ((ARM_VanillaSpreadOptionletArgSFRM*)arg)->GetSwaptionArg_Long();
	ARM_VanillaSwaptionArgSFRM* swaptionArg_Short = ((ARM_VanillaSpreadOptionletArgSFRM*)arg)->GetSwaptionArg_Short();
	ARM_VanillaSwaptionArgSFRM* swaptionArg_ToTimePayTime = ((ARM_VanillaSpreadOptionletArgSFRM*)arg)->GetSwaptionArg_ToTimePayTime();
	ARM_VanillaSwaptionArgSFRM* swaptionArg_StartTimePayTime = ((ARM_VanillaSpreadOptionletArgSFRM*)arg)->GetSwaptionArg_StartTimePayTime();


	//Pricing
	double asOfDate	= GetAsOfDate().GetJulian();
	double fromTime = evalTime;
	double toTime = resetTime;
	double	startTime1 = swapShortFloatStartTime;
	double   endTime1  = swapShortFloatEndTime;
	double	startTime2 = swapLongFloatStartTime;
	double   endTime2  = swapLongFloatEndTime;
	
	double   coeff1  = coeffShort;
	double   coeff2  = coeffLong;

	std::vector<double> VarCovar(3);
	std::vector<double> swapfwd(4);

	/*ComputeSwaptionVarCovar(fromTime,
							toTime,
							*swaptionArg_Short,
							*swaptionArg_Long,
							VarCovar,
							swapfwd);*/
	////////////////////////////////////////////////////////////////////
	ARM_ModelParams* modParamsSFRM = const_cast<ARM_SFRM* const>(this)->GetModelParams();

	//long Swaption	
	ARM_GP_VectorPtr mu_Long = swaptionArg_Long->GetMu();
	double var_Long = ((ARM_ModelParamsSFRM*)modParamsSFRM)->LocalVolatity(fromTime,toTime,*swaptionArg_Long,mu_Long);
	var_Long = var_Long*sqrt((toTime-fromTime)/K_YEAR_LEN);
	std::vector<double>& floatResetDate_Long = const_cast<std::vector<double>& const>(swaptionArg_Long->GetFloatResetTimes());
	int size_Long = floatResetDate_Long->size();
	ARM_VectorPtr swapfwd_Long	= swaptionArg_Long->GetSwapFwd();
	double shift_Long = swaptionArg_Long->GetAverageShift();
	

	//shortSwaption
	ARM_GP_VectorPtr mu_Short = swaptionArg_Short->GetMu();
	double var_Short = ((ARM_ModelParamsSFRM*)modParamsSFRM)->LocalVolatity(fromTime,toTime,*swaptionArg_Short,mu_Short);
	var_Short = var_Short*sqrt((toTime-fromTime)/K_YEAR_LEN);
	std::vector<double>& floatResetDate_Short = const_cast<std::vector<double>& const>(swaptionArg_Short->GetFloatResetTimes());
	int size_Short = floatResetDate_Short->size();
	ARM_VectorPtr swapfwd_Short	= swaptionArg_Short->GetSwapFwd();
	double shift_Short = swaptionArg_Short->GetAverageShift();


	//ToTimePayTime Swaption
	ARM_GP_VectorPtr mu_ToTimePayTime = swaptionArg_ToTimePayTime->GetMu();
	double var_ToTimePayTime = ((ARM_ModelParamsSFRM*)modParamsSFRM)->LocalVolatity(fromTime,toTime,*swaptionArg_ToTimePayTime,mu_ToTimePayTime);
	var_ToTimePayTime = var_ToTimePayTime*sqrt((toTime-fromTime)/K_YEAR_LEN);
	std::vector<double>& floatResetDate_ToTimePayTime = const_cast<std::vector<double>& const>(swaptionArg_ToTimePayTime->GetFloatResetTimes());
	int size_ToTimePayTime = floatResetDate_ToTimePayTime->size();
	ARM_VectorPtr swapfwd_ToTimePayTime	= swaptionArg_ToTimePayTime->GetSwapFwd();
	double shift_ToTimePayTime = swaptionArg_ToTimePayTime->GetAverageShift();

	

	//StartTimePayTime Swaption
	ARM_GP_VectorPtr mu_StartTimePayTime = swaptionArg_StartTimePayTime->GetMu();
	double var_StartTimePayTime = ((ARM_ModelParamsSFRM*)modParamsSFRM)->LocalVolatity(fromTime,toTime,*swaptionArg_StartTimePayTime,mu_StartTimePayTime);
	var_StartTimePayTime = var_StartTimePayTime*sqrt((toTime-fromTime)/K_YEAR_LEN);
	std::vector<double>& floatResetDate_StartTimePayTime = const_cast<std::vector<double>& const>(swaptionArg_StartTimePayTime->GetFloatResetTimes());
	int size_StartTimePayTime = floatResetDate_StartTimePayTime->size();
	ARM_VectorPtr swapfwd_StartTimePayTime	= swaptionArg_StartTimePayTime->GetSwapFwd();
	double shift_StartTimePayTime = swaptionArg_StartTimePayTime->GetAverageShift();
	

	//covariance_long_Short	&& covariance_StartTimePayTime_Long
	double covariance_Long_Short = 0.0;
	double covariance_StartTimePayTime_Long = 0.0;

	for(int k=0;k<size_Long;++k)
	{			
		double mu_Long_k = (*mu_Long)[k];
		double resetTime_Long = (*floatResetDate_Long)[k];
		
		for(int j=0;j<size_Short;++j)
		{
			double mu_Short_j = (*mu_Short)[j];
			double resetTime_Short = (*floatResetDate_Short)[j];
			double FWDCovar_k_j = dynamic_cast<ARM_ModelParamsSFRM *>(modParamsSFRM)->IntegratedCovariance(fromTime,toTime,resetTime_Long,resetTime_Short);
			double coeff_k_j = (mu_Long_k*mu_Short_j);
			covariance_Long_Short += coeff_k_j*FWDCovar_k_j;	
		}

		for(j=0;j<size_StartTimePayTime;++j)
		{
			double mu_StartTimePayTime_j = (*mu_StartTimePayTime)[j];
			double resetTime_StartTimePayTime = (*floatResetDate_StartTimePayTime)[j];
			double FWDCovar_k_j = dynamic_cast<ARM_ModelParamsSFRM *>(modParamsSFRM)->IntegratedCovariance(fromTime,toTime,resetTime_Long,resetTime_StartTimePayTime);
			double coeff_k_j = (mu_Long_k*mu_StartTimePayTime_j);
			covariance_StartTimePayTime_Long += coeff_k_j*FWDCovar_k_j;	
		}
	}


	//covariance_StartTimePayTime_Short
	double covariance_StartTimePayTime_Short =0.0;

	for(k=0;k<size_StartTimePayTime;++k)
	{			
		double mu_StartTimePayTime_k = (*mu_StartTimePayTime)[k];
		double resetTime_StartTimePayTime = (*floatResetDate_StartTimePayTime)[k];
		
		for(int j=0;j<size_Short;++j)
		{
			double mu_Short_j = (*mu_Short)[j];
			double resetTime_Short = (*floatResetDate_Short)[j];
			double FWDCovar_k_j = dynamic_cast<ARM_ModelParamsSFRM *>(modParamsSFRM)->IntegratedCovariance(fromTime,toTime,resetTime_StartTimePayTime,resetTime_Short);
			double coeff_k_j = (mu_StartTimePayTime_k*mu_Short_j);
			covariance_StartTimePayTime_Short += coeff_k_j*FWDCovar_k_j;	
		}
	}

	std::vector<double> VarCovar_P(3);
	std::vector<double> swapfwd_P(4);

	std::vector<double> VarCovar_P_1(3);
	std::vector<double> swapfwd_P_1(4);

	std::vector<double> VarCovar_P_2(3);
	std::vector<double> swapfwd_P_2(4);
	
	//Result_Long_Short
	VarCovar[0] = var_Short;
	VarCovar[1] = var_Long;
	VarCovar[2] = covariance_Long_Short;
	swapfwd.resize(4);
	swapfwd[0]=(*swapfwd_Short)[0];
	swapfwd[1]=(*swapfwd_Long)[0];
	swapfwd[2]=shift_Short;
	swapfwd[3]=shift_Long;

	//Result_ToTimePayTime_ToTimePayTime
	VarCovar_P[0] = var_ToTimePayTime;
	VarCovar_P[1] = var_ToTimePayTime;
	VarCovar_P[2] = 0;
	swapfwd_P.resize(4);
	swapfwd_P[0]=(*swapfwd_ToTimePayTime)[0];
	swapfwd_P[1]=(*swapfwd_ToTimePayTime)[0];
	swapfwd_P[2]=shift_ToTimePayTime;
	swapfwd_P[3]=shift_ToTimePayTime;

	// Result_StartTimePayTime_Short
	VarCovar_P_1[0] = var_StartTimePayTime;
	VarCovar_P_1[1] = var_Short;
	VarCovar_P_1[2] = covariance_StartTimePayTime_Short;
	swapfwd_P_1.resize(4);
	swapfwd_P_1[0]=(*swapfwd_StartTimePayTime)[0];
	swapfwd_P_1[1]=(*swapfwd_Short)[0];
	swapfwd_P_1[2]=shift_StartTimePayTime;
	swapfwd_P_1[3]=shift_Short;

	// Result_StartTimePayTime_Long
	VarCovar_P_2[0] = var_StartTimePayTime;
	VarCovar_P_2[1] = var_Long;
	VarCovar_P_2[2] = covariance_StartTimePayTime_Long;
	swapfwd.resize(4);
	swapfwd_P_2[0]=(*swapfwd_StartTimePayTime)[0];
	swapfwd_P_2[1]=(*swapfwd_Long)[0];
	swapfwd_P_2[2]=shift_StartTimePayTime;
	swapfwd_P_2[3]=shift_Long;
	/////////////////////////////////////////////////////////////////////////////////////////////

    double vol1 = VarCovar[0]/sqrt((toTime-evalTime)/K_YEAR_LEN);
	double vol2 = VarCovar[1]/sqrt((toTime-evalTime)/K_YEAR_LEN);
	double fwd1 = swapfwd[0];
	double fwd2 = swapfwd[1];
	double shift1 = swapfwd[2];
	double shift2 = swapfwd[3];
	
	//BS convexity Adjustment (As ARM&Summit)
	
	
	/*ComputeSwaptionVarCovar(fromTime,
							toTime,
							*swaptionArg_ToTimePayTime,
							*swaptionArg_ToTimePayTime,
							VarCovar_P,
							swapfwd_P);*/

    double theta_R = (toTime-evalTime)/K_YEAR_LEN;
	double theta_P = (payTime-toTime)/360.0;
	double Vol_P = VarCovar_P[0]/sqrt((toTime-evalTime)/K_YEAR_LEN);
	double swap_P =swapfwd_P[0];
	double shift_P = swapfwd_P[2];
	Vol_P=Vol_P*(shift_P+swap_P)/swap_P;
	double payLagConv = sqrt(exp(pow(theta_P*Vol_P*swap_P*sqrt(theta_R)/(1+theta_P*swap_P),2))-1);
	double vol1ForP = vol1*(fwd1+shift1)/fwd1;
	double vol2ForP = vol2*(fwd2+shift2)/fwd2;

	//correlPay1	
	/*ComputeSwaptionVarCovar(fromTime,
							toTime,
							*swaptionArg_StartTimePayTime,
							*swaptionArg_Short,
							VarCovar_P_1,
							swapfwd_P_1);*/

	double Correl_1 = VarCovar_P_1[2]/(VarCovar_P_1[0]*VarCovar_P_1[1]);

	//correlPay2
	/*ComputeSwaptionVarCovar(fromTime,
							toTime,
							*swaptionArg_StartTimePayTime,
							*swaptionArg_Long,
							VarCovar_P_2,
							swapfwd_P_2);*/


	double Correl_2 = VarCovar_P_2[2]/(VarCovar_P_2[0]*VarCovar_P_2[1]);

	double payLagConv_1 =1-Correl_1*payLagConv*sqrt(exp(vol1ForP*vol1ForP*theta_R)-1);
	double payLagConv_2 =1-Correl_2*payLagConv*sqrt(exp(vol2ForP*vol2ForP*theta_R)-1);

	//int decompfreq = Ccy->GetFixedPayFreq(),
	int decompfreq =1;
	int Tenor_1 = static_cast<int>(floor((endTime1-startTime1)/K_YEAR_LEN));
	int Tenor_2 = static_cast<int>(floor((endTime2-startTime2)/K_YEAR_LEN));
	double natLagConv_1 = 1-1/decompfreq*Tenor_1*fwd1/((1+decompfreq*fwd1)*(pow(1+1/decompfreq*fwd1,Tenor_1)-1));
	double natLagConv_2 = 1-1/decompfreq*Tenor_2*fwd2/((1+decompfreq*fwd2)*(pow(1+1/decompfreq*fwd2,Tenor_2)-1));
	
	natLagConv_1 *=(exp(vol1ForP*vol1ForP*theta_R)-1);
	natLagConv_2 *=(exp(vol2ForP*vol2ForP*theta_R)-1);
	
	natLagConv_1 +=1;
	natLagConv_2 +=1;
	double AdjFwd1ARM = (fwd1)*payLagConv_1*natLagConv_1;
	double AdjFwd2ARM = (fwd2)*payLagConv_2*natLagConv_2;		
	
	double AdjFwd1 = AdjFwd1ARM;
	double AdjFwd2 = AdjFwd2ARM;	
	double Correl = VarCovar[2]/(VarCovar[0]*VarCovar[1]);
	double optMat = toTime/K_YEAR_LEN;
	size_t nbStates			= states->size();
	ARM_VectorPtr spreadOption(new std::vector<double>(nbStates));
	size_t i;

	ARM_VectorPtr zcPay	= GetDiscountFunctor()->DiscountFactor(curveName,evalTime,payTime,states);

	for(i=0;i<nbStates;i++)
		(*spreadOption)[i]=notional*(*zcPay)[i]*payPeriod*SpreadOption(coeff1*(AdjFwd1+shift1), coeff2*(AdjFwd2+shift2), 
			vol1, vol2,0.0, 0.0, Correl, 0.0,strikes[i]+coeff2*shift2-coeff1*shift1, optMat, callPut);
	
	if (isArgLocal)
		delete arg;

	return spreadOption;
}


////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: GetVanillaSpreadOptionletArg
///	Returns: ARM_VanillaSwaptionArgSFRM
///	Action : computes all the required data for the pricing of a swaption
///				except the vol
////////////////////////////////////////////////////
ARM_VanillaSpreadOptionletArgSFRM* ARM_SFRM::GetVanillaSpreadOptionletArg( 
	const string& curveName,
		double evalTime,
		int	callPut,
		double startTime,
		double endTime,
		double resetTime,
		double payTime,
		double payPeriod,
		double notional,
		double coeffLong,
		double coeffShort,
		const std::vector<double>& strikes,
		double swapLongFloatStartTime,
		double swapLongFloatEndTime,
		const std::vector<double>& swapLongFixPayTimes,
		const std::vector<double>& swapLongFixPayPeriods,
		double swapShortFloatStartTime,
		double swapShortFloatEndTime,
		const std::vector<double>& swapShortFixPayTimes,
		const std::vector<double>& swapShortFixPayPeriods) const
{
	double asOfDate	= GetAsOfDate().GetJulian();
	ARM_ZeroCurvePtr ZcCurve	= GetZeroCurve();
	ARM_Currency* Ccy			= ZcCurve->GetCurrencyUnit();
	string curveCcy				= Ccy->GetCcyName();
	
	int fixFreq					= Ccy->GetFixedPayFreq();
	long fixDayCount			= Ccy->GetFixedDayCount();
	char fixCalendar[100];
	Ccy->CalcFixPayCal(fixCalendar);
	int fwdDayCount				= Ccy->GetLiborIndexDayCount();

	//create the Fix leg of CMS 

	//CMSLong (Fix Leg)
	ARM_Date startDate1(asOfDate+swapLongFloatStartTime);
	ARM_Date endDate1(asOfDate+swapLongFloatEndTime);
	ARM_DateStrip fixDateStrip1( startDate1, endDate1, fixFreq, fixDayCount, fixCalendar,
		K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,
		fixCalendar );
	std::vector<double>& longFixPayTimes= static_cast<ARM_GP_Vector*>(fixDateStrip1.GetPaymentDates()->Clone());
	std::vector<double>& longFixPayPeriods = static_cast<ARM_GP_Vector*>(fixDateStrip1.GetInterestTerms()->Clone()); 
	int i,nbFixFlowsLong=longFixPayTimes->size();    
	
	for(i=0;i<nbFixFlowsLong;++i)
		(*longFixPayTimes)[i] = (*longFixPayTimes)[i]-asOfDate;

	//CMSShort (Fix leg)
	ARM_Date startDate2(asOfDate+swapShortFloatStartTime);
	ARM_Date endDate2(asOfDate+swapShortFloatEndTime);
	ARM_DateStrip fixDateStrip2( startDate2, endDate2, fixFreq, fixDayCount, fixCalendar,
		K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,
		fixCalendar );
	std::vector<double>& shortFixPayTimes   = static_cast<ARM_GP_Vector*>(fixDateStrip2.GetPaymentDates()->Clone());
	std::vector<double>& shortFixPayPeriods = static_cast<ARM_GP_Vector*>(fixDateStrip2.GetInterestTerms()->Clone());
	int nbFixFlowsShort=shortFixPayTimes->size();    
	
	for(i=0;i<nbFixFlowsShort;++i)
		(*shortFixPayTimes)[i] = (*shortFixPayTimes)[i]-asOfDate;

	//Create ARM_VanillaSpreadOptionArg

	std::vector<double>& NullStrikes = new std::vector<double>(strikes);
	std::vector<double>& resetTimesVector					= new std::vector<double>(1,resetTime);
	std::vector<double>& payTimesVector					= new std::vector<double>(1,payTime);
	std::vector<double>& payPeriodsVector					= new std::vector<double>(1,payPeriod);
	std::vector<double>& notionalVector					= new std::vector<double>(1,notional);
	std::vector<double>& coeffLongVector					= new std::vector<double>(1,coeffLong);
	std::vector<double>& coeffShortVector					= new std::vector<double>(1,coeffShort);
	std::vector<double>& swapLongFloatStartTimeVector		= new std::vector<double>(1,swapLongFloatStartTime);
	std::vector<double>& swapLongFloatEndTimeVector		= new std::vector<double>(1,swapLongFloatEndTime);
	std::vector<double>& swapShortFloatStartTimeVector	= new std::vector<double>(1,swapShortFloatStartTime);
	std::vector<double>& swapShortFloatEndTimeVector		= new std::vector<double>(1,swapShortFloatEndTime);	

	ARM_VectorVector swapLongFixPayTimesVector;
	ARM_VectorVector swapLongFixPayPeriodsVector;
	ARM_VectorVector swapShortFixPayTimesVector;
	ARM_VectorVector swapShortFixPayPeriodsVector;
	
	swapLongFixPayTimesVector.push_back(longFixPayTimes);
	swapLongFixPayPeriodsVector.push_back(longFixPayPeriods);
	swapShortFixPayTimesVector.push_back(shortFixPayTimes);
	swapShortFixPayPeriodsVector.push_back(shortFixPayPeriods);
	
	
	ARM_VanillaSpreadOptionletArgSFRM* arg = new ARM_VanillaSpreadOptionletArgSFRM(
		curveName,
		evalTime,
		callPut,
		resetTime,
		swapLongFloatStartTime,
		endTime,
		resetTimesVector,
		payTimesVector,
		payPeriodsVector,
		notionalVector,
		coeffLongVector,
		coeffShortVector,
		NullStrikes,
		swapLongFloatStartTimeVector,
		swapLongFloatEndTimeVector,
		swapLongFixPayTimesVector,
		swapLongFixPayPeriodsVector,
		swapShortFloatStartTimeVector,
		swapShortFloatEndTimeVector,
		swapShortFixPayTimesVector,
		swapShortFixPayPeriodsVector);
	
	////Create attributes of ARM_VanillaSpreadOptionletArgSFRM :

	//create the swaptions arg	
	double swapNotional			= 100.0;
	std::vector<double> fixNotional (1,0.0);
	std::vector<double> floatNotional (1,0.0);
	ARM_GP_Matrix strikesPerState(1,0.0);
	size_t factorsNb			= FactorCount();
	ARM_PricingStatesPtr states = ARM_PricingStatesPtr( new ARM_PricingStates(1,1,1,factorsNb) );

	//itsSwaptionArg_Long	
	CC_NS(std,auto_ptr)<ARM_DateStrip> dateStrip1( GetFloatDateStrip(startDate1,endDate1) );
	std::vector<double>& floatResetDate1 = dateStrip1->GetResetDates();
	std::vector<double>& floatStartTimes1		= dateStrip1->GetFlowStartDates();
	std::vector<double>& floatPayTimes1		= dateStrip1->GetFlowEndDates();
	std::vector<double>& floatIntTerms1		= dateStrip1->GetInterestTerms();

		
	int size1 = floatResetDate1->size();
	for(i=0;i<size1;++i)
	{
		(*floatResetDate1)[i] = (*floatResetDate1)[i]-asOfDate;
		(*floatStartTimes1)[i] = (*floatStartTimes1)[i]-asOfDate;
		(*floatPayTimes1)[i] = (*floatPayTimes1)[i]-asOfDate;
	}

	double swapResetTime1 = (*floatResetDate1)[0];
	ARM_VanillaSwaptionArgSFRM* swaptionArgLong;
	swaptionArgLong = GetVanillaSwaptionArg(curveCcy, evalTime, swapResetTime1,
		swapNotional,fixNotional,floatNotional, swapLongFloatStartTime, swapLongFloatEndTime,*floatResetDate1,*floatStartTimes1,*floatPayTimes1,*floatIntTerms1, *longFixPayTimes, *longFixPayPeriods, strikesPerState,
		callPut, states);
	arg->SetSwaptionArg_Long(*swaptionArgLong);

	
	//itsSwaptionArg_Short		
	CC_NS(std,auto_ptr)<ARM_DateStrip> dateStrip2( GetFloatDateStrip(startDate2,endDate2) );
	std::vector<double>& floatResetDate2		= dateStrip2->GetResetDates();
	std::vector<double>& floatStartTimes2		= dateStrip2->GetFlowStartDates();
	std::vector<double>& floatPayTimes2		= dateStrip2->GetFlowEndDates();
	std::vector<double>& floatIntTerms2		= dateStrip2->GetInterestTerms();


	int size2 = floatResetDate2->size();
	for(i=0;i<size2;++i)
	{
		(*floatResetDate2)[i] = (*floatResetDate2)[i]-asOfDate;
		(*floatStartTimes2)[i] = (*floatStartTimes2)[i]-asOfDate;
		(*floatPayTimes2)[i] = (*floatPayTimes2)[i]-asOfDate;
	}

	double swapResetTime2 = (*floatResetDate2)[0];	
		
	ARM_VanillaSwaptionArgSFRM* swaptionArgShort;
	swaptionArgShort = GetVanillaSwaptionArg(curveCcy, evalTime, swapResetTime2,
		swapNotional,fixNotional,floatNotional, swapShortFloatStartTime, swapShortFloatEndTime, *floatResetDate2,*floatStartTimes2,*floatPayTimes2,*floatIntTerms2,*shortFixPayTimes, *shortFixPayPeriods, strikesPerState,
		callPut, states);
	arg->SetSwaptionArg_Short(*swaptionArgShort);


	//itsSwaptionArg_ToTimePayTime
	ARM_Date startDateToTime(asOfDate+resetTime);
	ARM_Date endDateToTime(asOfDate+payTime);
	ARM_DateStrip fixDateStripToTime( startDateToTime, endDateToTime, fixFreq, fixDayCount, fixCalendar,
		K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,
		fixCalendar );
	std::vector<double>& fixPayTimesToTime   = fixDateStripToTime.GetPaymentDates();
	std::vector<double>& fixPayPeriodsToTime = fixDateStripToTime.GetInterestTerms();
	int nbFixFlowsToTime=fixPayTimesToTime->size();    
	
	for(i=0;i<nbFixFlowsToTime;++i)
		(*fixPayTimesToTime)[i] = (*fixPayTimesToTime)[i]-asOfDate;
	CC_NS(std,auto_ptr)<ARM_DateStrip> dateStripToTime( GetFloatDateStrip(startDateToTime,endDateToTime) );
	std::vector<double>& floatResetDateToTime = dateStripToTime->GetResetDates();
	std::vector<double>& floatStartDateToTime		= dateStripToTime->GetFlowStartDates();
	std::vector<double>& floatEndDateToTime		= dateStripToTime->GetFlowEndDates();
	std::vector<double>& floatIntTermsToTime		= dateStripToTime->GetInterestTerms();

	int sizeToTime = floatResetDateToTime->size();
	for(i=0;i<sizeToTime;++i)
	{
		(*floatResetDateToTime)[i] = (*floatResetDateToTime)[i]-asOfDate;
		(*floatStartDateToTime)[i] = (*floatStartDateToTime)[i]-asOfDate;
		(*floatEndDateToTime)[i] = (*floatEndDateToTime)[i]-asOfDate;
	}

	double swapResetTimeToTime = (*floatResetDateToTime)[0];	
		
	ARM_VanillaSwaptionArgSFRM* swaptionArgToTime;
	swaptionArgToTime = GetVanillaSwaptionArg(curveCcy, evalTime, swapResetTimeToTime,
		swapNotional,fixNotional,floatNotional, resetTime, payTime, *floatResetDateToTime,*floatStartDateToTime,*floatEndDateToTime,*floatIntTermsToTime,*fixPayTimesToTime, *fixPayPeriodsToTime, strikesPerState,
		callPut, states);

	arg->SetSwaptionArg_ToTimePayTime(*swaptionArgToTime);


	//itsSwaptionArg_StartTimePayTime
	ARM_Date startDateStartTime(asOfDate+swapLongFloatStartTime);
	ARM_Date endDateStartTime(asOfDate+payTime);
	ARM_DateStrip fixDateStripStartTime( startDateStartTime, endDateStartTime, fixFreq, fixDayCount, fixCalendar,
		K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,
		fixCalendar );
	std::vector<double>& fixPayTimesStartTime   = fixDateStripStartTime.GetPaymentDates();
	std::vector<double>& fixPayPeriodsStartTime = fixDateStripStartTime.GetInterestTerms();
	int nbFixFlowsStartTime=fixPayTimesStartTime->size();    
	
	for(i=0;i<nbFixFlowsStartTime;++i)
		(*fixPayTimesStartTime)[i] = (*fixPayTimesStartTime)[i]-asOfDate;
	CC_NS(std,auto_ptr)<ARM_DateStrip> dateStripStartTime( GetFloatDateStrip(startDateStartTime,endDateStartTime) );
	std::vector<double>& floatResetDateStartTime = dateStripStartTime->GetResetDates();
	std::vector<double>& floatStartDateStartTime = dateStripStartTime->GetFlowStartDates();
	std::vector<double>& floatEndDateStartTime = dateStripStartTime->GetFlowEndDates();
	std::vector<double>& floatIntTermsStartToTime		= dateStripStartTime->GetInterestTerms();

	int sizeStartTime = floatResetDateStartTime->size();
	for(i=0;i<sizeStartTime;++i)
	{
		(*floatResetDateStartTime)[i] = (*floatResetDateStartTime)[i]-asOfDate;
		(*floatStartDateStartTime)[i] = (*floatStartDateStartTime)[i]-asOfDate;
		(*floatEndDateStartTime)[i] = (*floatEndDateStartTime)[i]-asOfDate;
	}
	double swapResetTimeStartTime = (*floatResetDateStartTime)[0];	
	
	ARM_VanillaSwaptionArgSFRM* swaptionArgStartTime;
	swaptionArgStartTime = GetVanillaSwaptionArg(curveCcy, evalTime, swapResetTimeStartTime,
		swapNotional,fixNotional,floatNotional, swapLongFloatStartTime, payTime, *floatResetDateStartTime,*floatStartDateStartTime,*floatEndDateStartTime,*floatIntTermsStartToTime,*fixPayTimesStartTime, *fixPayPeriodsStartTime, strikesPerState,
		callPut, states);

	arg->SetSwaptionArg_StartTimePayTime(*swaptionArgStartTime);

	return arg;
}


////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: NumMethodStateGlobalVariances
///	Returns : void
///	Action  : No default implementation since inappropriate for SFRM  model
////////////////////////////////////////////////////
void ARM_SFRM::NumMethodStateGlobalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& globalVariances ) const
{
	size_t factorsNb=FactorCount();
	size_t nbSteps	= timeSteps.size();
	size_t modelNb	= GetModelNb();
    double step		= timeSteps[0],
		   nextStep;
	size_t offsetIndex2	= nbSteps*modelNb;

#if defined(__GP_STRICT_VALIDATION)
	if( globalVariances.size()!= offsetIndex2 ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "localDrifts.size() != offsetIndex" );
#endif

	globalVariances.resize(nbSteps*(modelNb+1));

	/// fills the variance
	size_t i,k;
    globalVariances[offsetIndex2+0]=new ARM_GP_TriangularMatrix(factorsNb,0.0);
    for(i=0;i<nbSteps-1;++i)
    {
        nextStep=timeSteps[i+1];

		/// get the first bigger time
		ARM_VectorPtr variancePerFactor = ((ARM_ModelParamsSFRM*)GetModelParams())->IntegratedTreeStatesVariancePerFactor(step,nextStep);

        /// [i+1] => variance from 0 -> ti+1
		globalVariances[offsetIndex2+i+1] = new ARM_GP_TriangularMatrix(factorsNb);
        variancePerFactor = ((ARM_ModelParamsSFRM*)GetModelParams())->IntegratedTreeStatesVariancePerFactor(0.0,nextStep);
		for(k=0;k<factorsNb;++k)
			(*globalVariances[offsetIndex2+i+1])(k,k) =(*variancePerFactor)[k];
		step=nextStep;
    }	
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: PreProcessing
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_SFRM::PreProcessing(ARM_ModelFitter& modelFitter)
{
	const ARM_StdPortfolioPtr portfolio = modelFitter.GetPortfolio();
    /// To update the shift or beta Curve using the new Volatility Curve size.
    ARM_ModelParamVector::const_iterator foundvol = modelFitter.FindCalibParamWType( ARM_ModelParamType::Volatility );
    if(foundvol!= modelFitter.UnknownCalibParamIterator())
        ConvertToShiftorBetaParam(*portfolio);

    ARM_ModelParamVector::const_iterator foundShift = modelFitter.FindCalibParamWType( ARM_ModelParamType::Shift );
    ARM_ModelParamVector::const_iterator foundBeta = modelFitter.FindCalibParamWType( ARM_ModelParamType::Beta );
    if(foundShift!= modelFitter.UnknownCalibParamIterator() || foundBeta!= modelFitter.UnknownCalibParamIterator())
        itsUpdateVolSwapvolFRA = true;
	
	/// pre-computes data for each instruments in the model fitter
    if ( portfolio == ARM_StdPortfolioPtr(NULL))   
		throw Exception(__LINE__, __FILE__,ERR_INVALID_ARGUMENT,
			"You trying to calculate all coeff of relation Volfra/Volswap without Portfolio!");    

    size_t NbP = portfolio->GetSize();
	ARM_VanillaArgVector  ArgsVec(NbP);
    for( size_t i=0; i<NbP; i++)
    {
        if( ARM_Swaption* swaption = dynamic_cast<ARM_Swaption*>(portfolio->GetAsset(i)))
			ArgsVec[i] = ComputeCachedSwaptionData(swaption);
		else if ( ARM_SpreadOption* spreadoption = dynamic_cast<ARM_SpreadOption*>(portfolio->GetAsset(i)))
        	ArgsVec[i] = ComputeCachedSpreadOptionletData(spreadoption);

		else if ( ARM_CapFloor* capFloor = dynamic_cast<ARM_CapFloor*>(portfolio->GetAsset(i)))
        {
			/// after testing we found out that using ComputeCachedCapDigitalData(capFloor)
			///	was less efficient than before .... hence no caching for this one!
			/// itsArgsVec[i] = ARM_VanillaArgPtr(NULL);
			ArgsVec[i] = ComputeCachedCapData(capFloor);
		}
        else if ( ARM_Digital* digital = dynamic_cast<ARM_Digital*>(portfolio->GetAsset(i)))
        {
			/// after testing we found out that using ComputeCachedCapDigitalData(digital)
			///	was less efficient than before .... hence no caching for this one!
			/// itsArgsVec[i] = ARM_VanillaArgPtr(NULL);
			ArgsVec[i] = ComputeCachedDigitalData(digital);
		}
        else if ( ARM_CorridorLeg* corridorleg = dynamic_cast<ARM_CorridorLeg*>(portfolio->GetAsset(i)))
        {
            ArgsVec[i] = ComputeCachedCorridorLegData(corridorleg);
        }
		else if ( ARM_SumOpt* sumOpt = dynamic_cast<ARM_SumOpt*>(portfolio->GetAsset(i)))
		{
			ArgsVec[i] = NULL;
		}
		else
            throw Exception(__LINE__, __FILE__,ERR_INVALID_ARGUMENT,
				" current calibration supports only cap/floor and swaption!");   
    }
    modelFitter.SetVanillaArgVector(ArgsVec);

    ///After validate modelfitter, we call this function to manage correctly optimisation.
    GetModelParams()->PreProcessing(modelFitter,modelFitter.GetFactorNb());
}

////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routine : ValidateCalibMethod
///	Returns :void
///	Action  : call DefaultValidateWithModel
////////////////////////////////////////////////////
void ARM_SFRM::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	calibMethod.DefaultValidateWithModel(*this);
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: ComputeCachedSwaptionData
///	Returns: 
///	Action : stores all the cached data for the computation of a swaption
///				crucial for fast calibration!
////////////////////////////////////////////////////
ARM_VanillaSwaptionArgSFRM* ARM_SFRM::ComputeCachedSwaptionData( ARM_Swaption* swaption )
{
	double asOfDate	= GetAsOfDate().GetJulian();
	ARM_VanillaSwaptionArg* defaultArg = ARM_ConverterFromKernel::ConvertVanillaSwaption( swaption,  asOfDate );
	ARM_VanillaSwaptionArgSFRM* arg	= new ARM_VanillaSwaptionArgSFRM( *defaultArg );
	delete defaultArg;

	/// computes the additional arguments: averageShift,mu,fixAnnuity!
	arg->SetAverageShift( ((ARM_ModelParamsSFRM*) GetModelParams())->AverageShift(*arg->GetFloatResetTimes()) );
	arg->SetMu( ((ARM_ModelParamsSFRM*) GetModelParams())->ComputeVolSwapvolFRA(*arg,*this, arg->GetIsConstantNotional()) );
	ARM_PricingStatesPtr dumStates;
	if(! arg->GetIsConstantNotional())
		arg->SetFixAnnuityWithNominal (AnnuityWithNominal("Unused", arg->GetEvalTime(), *arg->GetFixPayTimes(), *arg->GetFixPayPeriods(), *arg->GetFixNotional(),dumStates)	)	;
	else
		arg->SetFixAnnuity( Annuity( "Unused", arg->GetEvalTime(), *arg->GetFixPayTimes(), *arg->GetFixPayPeriods(), dumStates) );	

	return arg;
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: ComputeCachedCapData
///	Returns: ARM_VanillaCapArgSFRM
///	Action : stores all the cached data for the computation of 
///          a cap floor for fast calibration!
////////////////////////////////////////////////////
ARM_VanillaCapArgSFRM* ARM_SFRM::ComputeCachedCapData( ARM_CapFloor* capFloor )
{
	double asOfDate					= GetAsOfDate().GetJulian();

	/// get the standard vanilla cap digital arg from the converter
	ARM_VanillaCapArg* defaultArg = ARM_ConverterFromKernel::ConvertVanillaCapFloor( capFloor, asOfDate );
	ARM_VanillaCapArgSFRM* arg = new ARM_VanillaCapArgSFRM( *defaultArg );
	delete defaultArg;

	/// computes the additional cached data libors, zcPays, shifts
	size_t size						= arg->GetResetTimes()->size();
	vector<ARM_VectorPtr> libors(size);
	vector<ARM_VectorPtr> zCPays(size);

	ARM_VectorPtr zcPay, libor, zcFwdEnd;
	ARM_PricingStatesPtr dumStates;

	for( size_t i=0; i<size; ++i )
	{
		zcPay	= GetDiscountFunctor()->DiscountFactor("",arg->GetEvalTime(), (*arg->GetPayTimes())[i], dumStates);
		libor	= GetFixingFunctor()->DiscountFactor("",arg->GetEvalTime(), (*arg->GetStartTimes())[i], dumStates);

		if( (*arg->GetEndTimes())[i] == (*arg->GetPayTimes())[i] && GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
			zcFwdEnd = zcPay;
        else
			zcFwdEnd = GetFixingFunctor()->DiscountFactor("",arg->GetEvalTime(),(*arg->GetEndTimes())[i],dumStates);

		size_t nbStates=zcPay->size();
		for(size_t j=0;j<nbStates; ++j)
			(*libor)[j] = ((*libor)[j]/(*zcFwdEnd)[j]-1.0)/(*arg->GetPayPeriods())[i];

		libors[i]		= libor;
		zCPays[i]		= zcPay;
	}

	arg->SetLibors( libors );
	arg->SetZCPays( zCPays );
	return arg;
}



////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: ComputeCachedDigitalData
///	Returns: ARM_VanillaDigitalArgSFRM
///	Action : stores all the cached data for the computation of 
///          a digital crucial for fast calibration!
////////////////////////////////////////////////////
ARM_VanillaDigitalArgSFRM* ARM_SFRM::ComputeCachedDigitalData( ARM_Digital* digital )
{
	double asOfDate					= GetAsOfDate().GetJulian();
	
	/// get the standard vanilla cap digital arg from the converter
	ARM_VanillaDigitalArg* defaultArg = ARM_ConverterFromKernel::ConvertVanillaDigital( digital,  asOfDate );
	ARM_VanillaDigitalArgSFRM* arg = new ARM_VanillaDigitalArgSFRM( *defaultArg );
	delete defaultArg;

	/// computes the additional cached data libors, zcPays, shifts
	size_t size						= arg->GetResetTimes()->size();
	vector<ARM_VectorPtr> libors(size);
	vector<ARM_VectorPtr> zCPays(size);

	ARM_VectorPtr zcPay, libor, zcFwdEnd;
	ARM_PricingStatesPtr dumStates;

	for( size_t i=0; i<size; ++i )
	{
		zcPay	= GetDiscountFunctor()->DiscountFactor("",arg->GetEvalTime(), (*arg->GetPayTimes())[i], dumStates);
		libor	= GetFixingFunctor()->DiscountFactor("",arg->GetEvalTime(), (*arg->GetStartTimes())[i], dumStates);

		if( (*arg->GetEndTimes())[i] == (*arg->GetPayTimes())[i] && GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
			zcFwdEnd = zcPay;
        else
			zcFwdEnd = GetFixingFunctor()->DiscountFactor("",arg->GetEvalTime(),(*arg->GetEndTimes())[i],dumStates);

		size_t nbStates=zcPay->size();
		for(size_t j=0;j<nbStates; ++j)
			(*libor)[j] = ((*libor)[j]/(*zcFwdEnd)[j]-1.0)/(*arg->GetPayPeriods())[i];

		libors[i]		= libor;
		zCPays[i]		= zcPay;
	}

	arg->SetLibors( libors );
	arg->SetZCPays( zCPays );
	return arg;
}


////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: ComputeCachedCorridorLegData
///	Returns: 
///	Action : stores all the cached data for the computation of a corridor eg
///				crucial for fast calibration!
////////////////////////////////////////////////////
ARM_VanillaCorridorLegArgSFRM* ARM_SFRM::ComputeCachedCorridorLegData(ARM_CorridorLeg* corridorLeg )
{
	double asOfDate					= GetAsOfDate().GetJulian();
	/// get the standard vanilla cap digital arg from the converter
	ARM_VanillaCorridorLegArg* defaultArg = ARM_ConverterFromKernel::ConvertVanillaCorridorLeg( corridorLeg,  asOfDate );
	
	/// creates the corresponding ARM_VanillaCapDigitalArgSFRM*
	ARM_VanillaCorridorLegArgSFRM* arg = new ARM_VanillaCorridorLegArgSFRM( *defaultArg );
	
	delete defaultArg;
	
	/// computes the additional cached data libors, zcPays, shifts
	size_t size	= arg->GetResetTimes()->size();
	vector<ARM_VectorPtr> libors(size);
	vector<ARM_VectorPtr> zCPays(size);
	
    vector<vector<ARM_VectorPtr> > ref_liborss(size);
	vector<vector<ARM_VectorPtr> > ref_DfFwdss(size);
    vector<vector<ARM_VectorPtr> > ref_zcEndss(size);    
    double periodPayEnd;
	
	ARM_VectorPtr zcPay, libor, zcFwdEnd;
    ARM_VectorPtr ref_libor, ref_DfFwd, ref_zcEnd;
	ARM_PricingStatesPtr dumStates;
	
	for( size_t i=0; i<size; ++i )
	{
		zcPay	= GetDiscountFunctor()->DiscountFactor("",arg->GetEvalTime(), (*arg->GetPayTimes())[i], dumStates);
		libor	= GetFixingFunctor()->DiscountFactor("",arg->GetEvalTime(), (*arg->GetStartTimes())[i], dumStates);
		
		if( (*arg->GetEndTimes())[i] == (*arg->GetPayTimes())[i] && GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
			zcFwdEnd = zcPay;
        else
			zcFwdEnd = GetFixingFunctor()->DiscountFactor("",arg->GetEvalTime(),(*arg->GetEndTimes())[i],dumStates);
        
        size_t nbStates = zcPay->size();
        for(size_t j=0;j<nbStates; ++j)
        {
            if(arg->GetIndexPaymentType() != IDXFIXED)
                (*libor)[j] = ((*libor)[j]/(*zcFwdEnd)[j]-1.0)/(*arg->GetFwdPaymentPeriod())[i] + (*arg->GetCouponMargin())[i];
            else
                (*libor)[j] = (*arg->GetCouponMargin())[i];   
        }
		
		libors[i]		= libor;
		zCPays[i]		= zcPay;
		
        size_t size_i	= (arg->GetRefIdxResettimes())[i]->size();
        double payTime = (*arg->GetPayTimes())[i];
        vector<ARM_VectorPtr>  ref_libors(size_i);
		vector<ARM_VectorPtr>  ref_DfFwds(size_i);
        vector<ARM_VectorPtr>  ref_zcEnds(size_i);
        ARM_VectorPtr   ref_DfFwd;
        for( size_t k=0; k<size_i; ++k )
        {
            ref_libor  = GetDiscountFunctor()->DiscountFactor("",arg->GetEvalTime(), (*(arg->GetRefIdxStarttimes())[i])[k], dumStates);
            ref_zcEnd  = GetDiscountFunctor()->DiscountFactor("",arg->GetEvalTime(), (*(arg->GetRefIdxEndtimes())[i])[k], dumStates);
            ref_DfFwd  = ARM_VectorPtr(new std::vector<double>(nbStates));
            
            for(j=0 ; j < nbStates; ++j)
                (*ref_libor)[j] = ((*ref_libor)[j]/(*ref_zcEnd)[j]-1.0)/(*(arg->GetRefFwdPeriods())[i])[k];
			
            double ref_endTime = (*(arg->GetRefIdxEndtimes())[i])[k];
            if (fabs(ref_endTime - payTime) > 5.0)
            {
                if (ref_endTime > payTime)
                {
                    periodPayEnd = CountYears(KACTUAL_360, payTime, ref_endTime);   
                    for(j=0 ; j < nbStates; ++j)
                    {
                        (*ref_DfFwd)[j] = (*ref_zcEnd)[j]/(*zcPay)[j];
                        (*ref_zcEnd)[j] = (1.0/(*ref_DfFwd)[j]-1.0) / periodPayEnd;
                        (*ref_DfFwd)[j]*=(periodPayEnd);
                    }
                }
                else 
                {
                    periodPayEnd = CountYears(KACTUAL_360, ref_endTime, payTime);   
                    for(j=0 ; j < nbStates; ++j)
                    {
                        (*ref_DfFwd)[j] = (*ref_zcEnd)[j] /(*zcPay)[j];
                        (*ref_zcEnd)[j] = ((*ref_DfFwd)[j]-1.0)/periodPayEnd;
                        (*ref_DfFwd)[j]*=(-periodPayEnd);
                    }
                }
            }
            ref_libors[k]= ref_libor ;
            ref_DfFwds[k]= ref_DfFwd;
            ref_zcEnds[k]= ref_zcEnd;
        }
        ref_liborss[i]= ref_libors;
        ref_DfFwdss[i]= ref_DfFwds;
        ref_zcEndss[i]= ref_zcEnds;
    }
	
	arg->SetLibors( libors );
	arg->SetZCPays( zCPays );
	
    arg->SetRef_Libors(ref_liborss);
	arg->SetRef_ZCEnds(ref_zcEndss); 
	arg->SetRef_DfFwds(ref_DfFwdss);
	
    return arg;
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: ComputeCachedSpreadOptionletData
///	Returns: 
///	Action : stores all the cached data for the computation of a swaption
///				crucial for fast calibration!
////////////////////////////////////////////////////
ARM_VanillaSpreadOptionletArgSFRM* ARM_SFRM::ComputeCachedSpreadOptionletData( ARM_SpreadOption* spreadoption)
{
	double asOfDate	= GetAsOfDate().GetJulian();
	ARM_VanillaSpreadOptionArg* defaultArg = ARM_ConverterFromKernel::ConvertVanillaSpreadOption( spreadoption,  asOfDate );
	ARM_VanillaSpreadOptionletArgSFRM* arg	= new ARM_VanillaSpreadOptionletArgSFRM( *defaultArg );
	delete  defaultArg;

	string curveName = arg->GetCurveName();
	double evalTime =arg->GetEvalTime();
	int	callPut=arg->GetCallPut();
	double resetTime= (*(arg->GetResetTimes()))[0];
	double payTime =(*(arg->GetPayTimes()))[0];
	double payPeriod =(*(arg->GetPayPeriods()))[0];
	double notional =(*(arg->GetNotional()))[0];
	double coeffLong = (*(arg->GetCoeffLong()))[0];
	double coeffShort = (*(arg->GetCoeffShort()))[0];
	const std::vector<double>& strikes = arg->GetStrikes()[0];
	double swapLongFloatStartTime = (*(arg->GetSwapLongFloatStartTime()))[0];
	double swapLongFloatEndTime = (*(arg->GetSwapLongFloatEndTime()))[0];
	const std::vector<double>& swapLongFixPayTimes= *arg->GetSwapLongFixPayTimes()[0];
	const std::vector<double>& swapLongFixPayPeriods = *arg->GetSwapLongFixPayPeriods()[0];
	double swapShortFloatStartTime = (*(arg->GetSwapShortFloatStartTime()))[0];
	double swapShortFloatEndTime = (*(arg->GetSwapShortFloatEndTime()))[0];
	const std::vector<double>& swapShortFixPayTimes = *arg->GetSwapShortFixPayTimes()[0];
	const std::vector<double>& swapShortFixPayPeriods = *arg->GetSwapShortFixPayPeriods()[0];

	arg = GetVanillaSpreadOptionletArg( 
						curveName,
						evalTime,
						callPut,
						0,
						0,
						resetTime,
						payTime,
						payPeriod,
						notional,
						coeffLong,
						coeffShort,
						strikes,
						swapLongFloatStartTime,
						swapLongFloatEndTime,
						swapLongFixPayTimes,
						swapLongFixPayPeriods,
						swapShortFloatStartTime,
						swapShortFloatEndTime,
						swapShortFixPayTimes,
						swapShortFixPayPeriods);

	return arg;
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: PostProcessing
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_SFRM::PostProcessing(const ARM_ModelFitter& modelFitter)
{
    itsCurrentArg	= NULL;
    itsUpdateVolSwapvolFRA = false;

	GetModelParams()->PostProcessing(modelFitter,this);
    ARM_ModelParamVector::const_iterator foundBeta = modelFitter.FindCalibParamWType( ARM_ModelParamType::Shift );
    if(foundBeta!= modelFitter.UnknownCalibParamIterator())
         ((ARM_ModelParamsSFRM*) GetModelParams())->ConvertToBetaParam(*this, *modelFitter.GetPortfolio()); 
}

////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: AdviseCurrentCalibSecIndex
///	Returns : void
///	Action  : Advises the model that the current index of the calibration security
///             is the index given ... The model advises just the model params
////////////////////////////////////////////////////
void ARM_SFRM::AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter)
{
    itsCurrentArg = modelFitter.GetVanillaArgVector(index);
}

////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: Re_InitialiseCalibParams
///	Returns : void
///	Action  : Re_Intialise the the calib params
////////////////////////////////////////////////////
void ARM_SFRM::Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter)
{
    ARM_ModelParamVector::const_iterator foundBeta = modelFitter.FindCalibParamWType( ARM_ModelParamType::Beta );
	if( foundBeta!= modelFitter.UnknownCalibParamIterator())
         ((ARM_ModelParamsSFRM*) GetModelParams())->ConvertToShiftParam(*this, *modelFitter.GetPortfolio()); 

	/// To normailize volatility curve
    ARM_ModelParamVector::const_iterator foundVol = modelFitter.FindCalibParamWType( ARM_ModelParamType::Volatility );
   
    if(foundVol != modelFitter.UnknownCalibParamIterator() )
    {
		ARM_CurveModelParam* calibVol = dynamic_cast<ARM_CurveModelParam*>(*foundVol);

		/// are we really a curveCalibParam!
		if( !calibVol )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME 
				+ ": this should be a curve Calib Param!");
        for(int i=0; i < (*foundVol)->size(); ++i)
        {
            double guess_initial= calibVol->GetInitialCurve()->GetOrdinates()[i];
			double timelag		= calibVol->GetInitialCurve()->GetAbscisses()[i];
            double beta = static_cast<ARM_ModelParamsSFRM*>(GetModelParams())->BetaValue(timelag);
            calibVol->SetValueAtPoint(i,guess_initial*beta);
        }
    }
};

////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: AdviseCurrentCalib
///	Returns : void
///	Action  : Advises the model that of the calibration  ... 
///           The model advises just the model params
////////////////////////////////////////////////////
void ARM_SFRM::AdviseCurrentCalib(ARM_ModelFitter& modelFitter)
{
    ARM_ModelParamVector::const_iterator foundBeta = modelFitter.FindCalibParamWType( ARM_ModelParamType::Beta );
    if( foundBeta!= modelFitter.UnknownCalibParamIterator())
         ((ARM_ModelParamsSFRM*) GetModelParams())->ConvertToShiftParam(*this, *modelFitter.GetPortfolio()); 
	
	ARM_ModelParamVector::const_iterator foundMeanRev = modelFitter.FindCalibParamWType( ARM_ModelParamType::MeanReversion );
	ARM_ModelParamVector::const_iterator foundVol = modelFitter.FindCalibParamWType( ARM_ModelParamType::Volatility );
	ARM_ModelParamVector::const_iterator foundCorrel = modelFitter.FindCalibParamWType( ARM_ModelParamType::Correlation);
	ARM_ModelParamVector::const_iterator foundBrownianCorrel = modelFitter.FindCalibParamWType( ARM_ModelParamType::BrownianCorrelation);

	if(    (foundMeanRev!= modelFitter.UnknownCalibParamIterator() ) 
		|| (foundVol!= modelFitter.UnknownCalibParamIterator())
		|| (foundCorrel!= modelFitter.UnknownCalibParamIterator())
		|| (foundBrownianCorrel != modelFitter.UnknownCalibParamIterator()) ) 
			GetModelParams()->PostProcessing(modelFitter,this); 
}


////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: ReInit
///	Returns : ARM_PricingStatesPtr
///	Action  : ReInit is for multi-loop numerical method to init lightwise
///				the model and its numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_SFRM::ReInit()
{
	ResetDFMap();
	return GetNumMethod()->ReInit(*this);
}


////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: AdviseBreakPointTimes
///	Returns : void
///	Action  : sets the corresponding suggested break point times to the model param
///////////////////////////////////////
void ARM_SFRM::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio,
				 ARM_ModelParam* inputModelParam,
				 size_t factorNb)
{
	ARM_CurveModelParam* modelParam = dynamic_cast<ARM_CurveModelParam*>(inputModelParam);
	if( !modelParam )
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "expected an ARM_CurveModelParam!");

    double asOfDate = GetAsOfDate().GetJulian();
    int size1       = portfolio->GetSize();  
    std::vector<double>  tmpdates;
    int i;
    
    switch( modelParam->GetType() )
    {
	case ARM_ModelParamType::Beta:
	case ARM_ModelParamType::Shift:
	case ARM_ModelParamType::BrownianCorrelation:
	case ARM_ModelParamType::Volatility:
        {
            if (portfolio->GetSize())
            {
			    double date = portfolio->GetAsset(0)->GetResetDates()->Elt(0) - asOfDate;
			    tmpdates.push_back(date);
			    for(i=1; i<size1; i++) 
			    {
				    double resetlag = portfolio->GetAsset(i)->GetResetDates()->Elt(0) - asOfDate;
				    if(fabs (date - resetlag) > FRMVOL_LAG_THRESHOLD)
				    {
					    tmpdates.push_back(resetlag);
					    date = resetlag;
				    }
				    else
				    {
					    /// ignore this instrument
					    portfolio->SetWeight(0.0,i);
				    }
			    }
			    modelParam->UpdateValues(&tmpdates);
            }
        }
        break;
		
		/// just return NULL
	case ARM_ModelParamType::Correlation:
	case ARM_ModelParamType::MeanReversion:
		break;
	default:
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
            "Unknown type... an SFRM only supports volatility, meanreversion, beta, shift, correlation, brownian correlation" );
    }
}

////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: SetFixStartDate
///	Returns : 
///	Action  : Set the the fixed start date
////////////////////////////////////////////////////
void ARM_SFRM::SetFixStartDate(
	const ARM_Date& fixStartDate)
{
	itsFixStartDate = fixStartDate;
	itsIsFixStartDate = true;
}


////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: SetFixEndDate
///	Returns : 
///	Action  : Set the the fixed end date
////////////////////////////////////////////////////
void ARM_SFRM::SetFixEndDate(
	const ARM_Date& fixEndDate)
{
	itsFixEndDate = fixEndDate;
	itsIsFixEndDate = true;
}

////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: VolatilitiesAndCorrelationTimesSteps
///	Returns : 
///	Action  : Volatilities and correlation time steps for PDEs
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_SFRM::VolatilitiesAndCorrelationTimesSteps() const
{

	// Time steps where the model params are changing.
	const ARM_ModelParams* modelParams = GetModelParams();
	const ARM_ModelParamsSFRM* modelParamsSFRM = dynamic_cast<const ARM_ModelParamsSFRM*> (modelParams);

	if( !modelParamsSFRM )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"HW Model params Expected !");

	return modelParamsSFRM->ModelParamsTimeSteps();
}

////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: VolatilitiesAndCorrelationTimesSteps
///	Returns : 
///	Action  : function for the generic tree (2nd generation)
////////////////////////////////////////////////////
void ARM_SFRM::VolatilitiesAndCorrelations( const std::vector<double>& timeSteps, 
		ARM_GP_MatrixPtr& vols,
		ARM_GP_MatrixPtr& d1Vols,
		ARM_GP_MatrixPtr& correls,
		bool linearVol) const
{
	if (!linearVol)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "ARM_SFRM::VolatilitiesAndCorrelations : only linearVol case is implemented");

	std::vector<double> times = ((ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility)).GetCurve()->GetAbscisses();
	std::vector<double> values= ((ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility)).GetCurve()->GetOrdinates();

	/// Since volatilies are stepUpRight!!
	for( std::vector<double>::iterator iter = values.begin() ; iter != values.end()-1 ; ++iter)
		(*iter) = *(iter+1);

	std::vector<double> volsVec,d1VolsVec;
	VectorValuesAndDerivativesLinearMidPoints(times,values,timeSteps,volsVec,d1VolsVec);
	for(size_t i=0; i<d1VolsVec.size(); ++i )
		d1VolsVec[i] *= K_YEAR_LEN;

	/// factor by line and zero correl because in one factor!
	vols	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(1,volsVec.size(),&volsVec[0]) );
	d1Vols	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(1,d1VolsVec.size(),&d1VolsVec[0]));
	correls	= ARM_GP_MatrixPtr( NULL );
}

////////////////////////////////////////////////////
///	Class   : ARM_SFRM
///	Routines: VolatilitiesAndCorrelationTimesSteps
///	Returns : 
///	Action  : function for the generic tree (2nd generation)
////////////////////////////////////////////////////
void ARM_SFRM::EulerLocalDrifts(const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const
{
	relativeDrifts= ARM_GP_MatrixPtr( new ARM_GP_Matrix(timeSteps.size(),1, 0.0 ) );
	absoluteDrifts= ARM_GP_MatrixPtr( new ARM_GP_Matrix(timeSteps.size(),1, 0.0 ) );
}


CC_END_NAMESPACE()

