/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file infswopvol.cpp
 *
 *  \brief object to compute the inflation swaption vol
 *		from inflation year on year and zero coupon vol
 *
 *	\author  N. Belgrade, E. Benhamou
 *	\version 1.0
 *	\date September 2004
 *
 *	\version 2.0
 *	\date March 2005
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

#include "gpinflation/infswopvol.h"

// gpinflation
#include "gpinflation/infBSModel.h"
#include "gpinflation/infcurv.h"
#include "gpinflation/infdata.h"
#include "gpinflation/infidx.h"
#include "gpinflation/sparsevolcube.h"

/// gpbase
#include "gpbase/checkarg.h"
#include "gpbase/defaultargs.h"
#include "gpbase/autocleaner.h"
#include "gpbase/gpvector.h"
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/gplinalgconvert.h"

/// kernel
#include <glob/expt.h>
#include <glob/dates.h>
#include <crv/volint.h>
#include <crv/correlmanager.h>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
/// Static Variables
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class      : ARM_InfVolComputation_Producer
/// Static	   : DefaultTenors
///	Description: table for default tenors
////////////////////////////////////////////////////
double ARM_InfVolComputation_Producer::DefaultTenors[] =
{
	1, 2, 3, 4, 5, 6, 7, 10, 15, 20, 25, 30
};

////////////////////////////////////////////////////
///	Class      : ARM_InfVolComputation_Producer
/// Static	   : DefaultTenorsSize
///	Description: size of the default tenors data
////////////////////////////////////////////////////
size_t ARM_InfVolComputation_Producer::DefaultTenorsSize = sizeof(DefaultTenors)/ sizeof(DefaultTenors[0]) ;


////////////////////////////////////////////////////
///	Class      : ARM_InfVolComputation_Producer
/// Static	   : DefaultExpiries
///	Description: table for default expiries
////////////////////////////////////////////////////
double ARM_InfVolComputation_Producer::DefaultExpiries[] = 
{
	1./12., 2./12., 3./12., 6./12.,
	1, 2, 3, 4, 5, 6, 7, 10, 15, 20, 25, 30
};


////////////////////////////////////////////////////
///	Class      : ARM_InfVolComputation_Producer
/// Static	   : DefaultExpiriesSize
///	Description: size of the default expiries data
////////////////////////////////////////////////////
size_t ARM_InfVolComputation_Producer::DefaultExpiriesSize = sizeof(DefaultExpiries)/ sizeof(DefaultExpiries[0]) ;


////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	Class      : ARM_InfVolComputation_Producer
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_Std
///	Routine: Constructor
////////////////////////////////////////////////////

ARM_InfVolComputation_Producer::ARM_InfVolComputation_Producer(ARM_InfBSModel* infIRBSModel)
:	ARM_Object(), itsInfBSModel(infIRBSModel)
{	SetName(ARM_INFSWWOPTVOL_PRODUCER);	}



////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_Std
///	Routine: Copy Constructor
////////////////////////////////////////////////////
ARM_InfVolComputation_Producer::ARM_InfVolComputation_Producer(
	const ARM_InfVolComputation_Producer& rhs )
:	ARM_Object(rhs), itsInfBSModel(rhs.itsInfBSModel )
{}



////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_Std
///	Routine: Assignment operator
////////////////////////////////////////////////////
ARM_InfVolComputation_Producer& ARM_InfVolComputation_Producer::operator=(
	const ARM_InfVolComputation_Producer& rhs )
{
	if( this!= &rhs )
	{
		ARM_Object::operator =( rhs );
		itsInfBSModel = rhs.itsInfBSModel;
	}
	return *this;
}

///////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_Std
///	Routine: Destructor
////////////////////////////////////////////////////
ARM_InfVolComputation_Producer::~ARM_InfVolComputation_Producer()
{}



////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer
///	Routine: View
////////////////////////////////////////////////////
void ARM_InfVolComputation_Producer::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[200];
	
	/// first determine that the file is not already opened
    if ( NULL == ficOut )
    {
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w");
    }
    else
		fOut = ficOut;

	/// use the method to string
	/// and just says the type and what is in it
    fprintf(fOut, "%s", toString().c_str() );

	/// to allow to have nested view
    if ( NULL == ficOut )
       fclose(fOut);
}




///////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_Std
///	Routine: GetDefaultDataIfMissing
///	Returns: double
///	Action : Get DefaultData If Missing
////////////////////////////////////////////////////

void ARM_InfVolComputation_Producer::GetDefaultDataIfMissing(
	ARM_GP_Vector*& tenors,
	ARM_GP_Vector*& expiries )
{
	/// get default data if missing
	/// could be either a null pointor or an empty vector!
	/// both for tenors and expiries
	CopyArgsPointorInPlace<ARM_GP_Vector,double*>(tenors,
		ARM_InfVolComputation_Producer::DefaultTenors,
		ARM_InfVolComputation_Producer::DefaultTenorsSize);

	CopyArgsPointorInPlace<ARM_GP_Vector,double*>(expiries,
		ARM_InfVolComputation_Producer::DefaultExpiries,
		ARM_InfVolComputation_Producer::DefaultExpiriesSize);
}



////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
///	Class      : ARM_InfVolComputation_Producer_Std
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_Std
///	Routine: Constructor
////////////////////////////////////////////////////

ARM_InfVolComputation_Producer_Std::ARM_InfVolComputation_Producer_Std(	ARM_InfBSModel*	infIRBSModel )
:	ARM_InfVolComputation_Producer(infIRBSModel)
{
	SetName(ARM_INFSWWOPTVOL_STD_PRODUCER);
}


////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_Std
///	Routine: Copy Constructor
////////////////////////////////////////////////////
ARM_InfVolComputation_Producer_Std::ARM_InfVolComputation_Producer_Std(
	const ARM_InfVolComputation_Producer_Std& rhs )
:	ARM_InfVolComputation_Producer(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_Std
///	Routine: Assignement operator
////////////////////////////////////////////////////
ARM_InfVolComputation_Producer_Std& ARM_InfVolComputation_Producer_Std::operator=( 
	const ARM_InfVolComputation_Producer_Std& rhs )
{
	if( this != &rhs )
		ARM_InfVolComputation_Producer::operator =( rhs );
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_Std
///	Routine: Destructor
////////////////////////////////////////////////////
ARM_InfVolComputation_Producer_Std::~ARM_InfVolComputation_Producer_Std()
{
	/// the pointor is a shared ressource hence no need to delete it!
}

////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_Std
///	Routine: Clone
////////////////////////////////////////////////////
ARM_Object* ARM_InfVolComputation_Producer_Std::Clone()
{
	return new ARM_InfVolComputation_Producer_Std(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_Std
///	Routine: Clone
////////////////////////////////////////////////////
string ARM_InfVolComputation_Producer_Std::toString( const string& indent )
{
	return "ARM_InfVolComputation_Producer_Std";
}


///////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_Std
///	Routine: VolYoY_to_VolSwp
///	Returns: double
///	Action : general function to compute the volatility of inf swaption from the one 
///				of year on year. We assume a correlation of 1 between the various year on
///				year. This is conservative if we want to compute an upper boundary on the
///				inflation volatility.
////////////////////////////////////////////////////

double ARM_InfVolComputation_Producer_Std::VolYoY_to_VolSwp(
	/// fix leg
	ARM_GP_Vector* pFixLegDF, 
	ARM_GP_Vector* pFixLegIRYearFrac, 
	ARM_GP_Vector* pFixLegDFVol,

	//// inflation
	ARM_GP_Vector* pFwdCPIRatio,		/// this is supposed to be convexity corrected
	ARM_GP_Vector* pVol_YoY,
	ARM_GP_Vector* pInfLegDF,
	ARM_GP_Vector* pInfLegDFVol,
	ARM_Matrix* AvgFloatCPIIRCorrel,
	ARM_Matrix* AvgFloatFixCPIIRCorrel,
	ARM_Matrix* AvgFloatCPICPICorrel,
	ARM_Matrix* AvgFixIRIRCorrel,
	ARM_Matrix* AvgFloatIRIRCorrel,
	ARM_Matrix* AvgFloatFixIRIRCorrel )
{
	/// Validation function that checks that the two vectors have the same size
	/// fix leg part
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pFixLegIRYearFrac,	*pFixLegDF,	"pFixLegIRYearFrac",		"pFixLegDF" );
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pFixLegDFVol,		*pFixLegDF,	"FixLegDFIntegratedVol",	"pFixLegDF" );

	/// inf leg part
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pVol_YoY,				*pFwdCPIRatio,	"Vol_YoY",		"FwdCPIRatio" );
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pInfLegDF,				*pFwdCPIRatio,	"InfLegDF",		"FwdCPIRatio" );
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pInfLegDFVol,			*pFwdCPIRatio,	"pInfLegDFVol",	"FwdCPIRatio" );

	size_t i;
	double FixedLegAnnuity = 0.;
	ARM_GP_Vector fixLegTerms(pFixLegDF->size());
	
	/// Construction of an annuity
	for( i=0; i<pFixLegDF->size(); ++i )
	{
		fixLegTerms[i]	 = pFixLegDF->Elt(i) * pFixLegIRYearFrac->Elt(i);
		FixedLegAnnuity += fixLegTerms[i];
	}


	ARM_GP_Vector iTerms(pFwdCPIRatio->size());
	double FloatLegPV = 0.;
	/// iterms = CPI Ratio * DF(i)
	for( i=0; i<pFwdCPIRatio->size(); ++i )
	{
		iTerms[i] =  pFwdCPIRatio->Elt(i) * pInfLegDF->Elt(i);
		FloatLegPV += iTerms[i] ;
	}


	/// sqrVolFixLeg is the square of the fix leg vol
	/// volFixLeg = Sum(i=0...nFix, DF(i)*VolBond(i))/FixLegAnnuity
	size_t j;
	double sqrVolFixLeg = 0.0;
	for( i=0; i<pFixLegDF->size(); ++i )
		for( j=0; j<pFixLegDF->size(); ++j )
			sqrVolFixLeg += fixLegTerms[i] * fixLegTerms[j]
				* pFixLegDFVol->Elt(i) * pFixLegDFVol->Elt(j)	* AvgFixIRIRCorrel->Elt(i,j);
	sqrVolFixLeg /= FixedLegAnnuity*FixedLegAnnuity;

	/// the correlation used here is the one between the inflation CPI ratio and the interest rate discount factor
	/// sqrVolInfLeg is the square of the inflation leg vol
	double sqrVolInfLeg = 0.0;
	for( i=0; i<pFwdCPIRatio->size(); ++i )
		for( j=0;j<pFwdCPIRatio->size(); ++j )
			sqrVolInfLeg += iTerms[i] * iTerms[j] * (
				  pVol_YoY->Elt(i) * pVol_YoY->Elt(j)			* AvgFloatCPICPICorrel->Elt(i,j)
				+ pVol_YoY->Elt(i) * pInfLegDFVol->Elt(j)		* AvgFloatCPIIRCorrel->Elt(i,j)
				+ pInfLegDFVol->Elt(i) * pVol_YoY->Elt(j)		* AvgFloatCPIIRCorrel->Elt(i,j)
				+ pInfLegDFVol->Elt(i) * pInfLegDFVol->Elt(j)	* AvgFloatIRIRCorrel->Elt(i,j));
	sqrVolInfLeg /= FloatLegPV*FloatLegPV;

	double covarianceTerm = 0.0;
	for( i=0; i<pFwdCPIRatio->size(); ++i )
		for ( j=0; j<pFixLegDF->size(); ++j )
			covarianceTerm += iTerms[i] * fixLegTerms[j] * (
				  pVol_YoY->Elt(i) * pFixLegDFVol->Elt(j)		* AvgFloatFixCPIIRCorrel->Elt(i,j)
				+ pInfLegDFVol->Elt(i) * pFixLegDFVol->Elt(j)	* AvgFloatFixIRIRCorrel ->Elt(i,j) );
	covarianceTerm /= FloatLegPV*FixedLegAnnuity;

	double vol  = sqrVolInfLeg + sqrVolFixLeg - 2.*covarianceTerm;

	/// test before sqrt!
	if( vol < K_NEW_DOUBLE_TOL )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": vol is not positive!");
	vol	 = sqrt(vol);
	return vol;
}

///////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_Std
///	Routine: VolZC_to_VolOATSwp
///	Returns: double
///	Action : general function to compute the volatility of inf OAT swaption from the one 
///				of zero coupon. We assume a correlation of 1 between the various zero coupon
///				year. This is conservative if we want to compute an upper boundary on the
///				inflation volatility.
////////////////////////////////////////////////////
double ARM_InfVolComputation_Producer_Std::VolZC_to_VolOATSwp(										
		ARM_GP_Vector* pFixLegDF,
		ARM_GP_Vector* pFixLegIRYearFrac,
		ARM_GP_Vector* pFixLegDFVol,
		ARM_GP_Vector* pFwdCPI,
		ARM_GP_Vector* pVol_ZC,
		ARM_GP_Vector* pInfLegDF,
		ARM_GP_Vector* pInfLegDFVol,
		ARM_Matrix* AvgFloatCPIIRCorrel,
		ARM_Matrix* AvgFloatFixCPIIRCorrel,
		ARM_Matrix* AvgFloatCPICPICorrel,
		ARM_Matrix* AvgFixIRIRCorrel,
		ARM_Matrix* AvgFloatIRIRCorrel,
		ARM_Matrix* AvgFloatFixIRIRCorrel,
		double	coupon,
		bool	mode)
{
	/// Validation function that checks that the two vectors have the same size
	/// fix leg part
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pFixLegIRYearFrac,	*pFixLegDF,	"pFixLegIRYearFrac",	 "pFixLegDF" );
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pFixLegDFVol,		*pFixLegDF,	"FixLegDFIntegratedVol", "pFixLegDF" );

	/// inf leg part
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pVol_ZC,		*pFwdCPI,	"Vol_ZC",		"FwdCPI"	);
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pInfLegDF,		*pFwdCPI,	"InfLegDF",		"FwdCPI"	);
	CC_NS(ARM_Check,CheckSameArgSize)< ARM_GP_Vector, ARM_GP_Vector >( *pInfLegDFVol,	*pFwdCPI,	"pInfLegDFVol",	"FwdCPI"	);

	

	size_t i=0, j=0;
	double FixedLegAnnuity = 0.;
	ARM_GP_Vector fixLegTerms(pFixLegDF->size());
	ARM_GP_Vector Coupon (pVol_ZC->size(), 1.);

	/// Construction of an annuity
	for( i=0; i<pFixLegDF->size(); ++i )
	{
		fixLegTerms[i]	 = pFixLegDF->Elt(i) * pFixLegIRYearFrac->Elt(i);
		FixedLegAnnuity += fixLegTerms[i];
	}

	/// Store the coupon yield
	if (mode ==true) 
	{
		for( i=0; i<Coupon.size(); ++i )
		{
			{Coupon[i] = coupon;}
		}
		Coupon[Coupon.size()-1] +=1.;
	}

	ARM_GP_Vector iTerms( pFwdCPI->size());
	double FloatLegPV = 0.;

	/// Construction of a periodic cash flows 
	for( i=0; i<pFwdCPI->size(); i++)
	{
		iTerms[i] = Coupon[i]*pInfLegDF->Elt(i)*pFwdCPI->Elt(i);
		FloatLegPV += iTerms[i];
	}

	/// sqrVolFixLeg is the square of the fix leg vol
	/// volFixLeg = Sum(i=0...nFix, DF(i)*VolBond(i))/FixLegAnnuity
	double sqrVolFixLeg = 0.0;

	for( i=0; i<pFixLegDF->size(); ++i )
		for( j=0; j<pFixLegDF->size(); ++j )
			sqrVolFixLeg += fixLegTerms[i] * fixLegTerms[j] * pFixLegDFVol->Elt(i) * pFixLegDFVol->Elt(j) * AvgFixIRIRCorrel->Elt(i,j);
	sqrVolFixLeg /= FixedLegAnnuity*FixedLegAnnuity;


	/// The correlation used here is the one between the inflation CPI and the interest rate discount factor
	/// sqrVolInfLeg is the square of the inflation leg vol
	double sqrVolInfLeg = 0.0;

	for( i=0; i<pFwdCPI->size(); ++i )
		for( j=0;j<pFwdCPI->size(); ++j )
			sqrVolInfLeg += iTerms[i] * iTerms[j] * (
				  pVol_ZC->Elt(i) * pVol_ZC->Elt(j)				* AvgFloatCPICPICorrel->Elt(i,j)
				+ pVol_ZC->Elt(i) * pInfLegDFVol->Elt(j)		* AvgFloatCPIIRCorrel->Elt(i,j)
				+ pInfLegDFVol->Elt(i) * pVol_ZC->Elt(j)		* AvgFloatCPIIRCorrel->Elt(i,j)
				+ pInfLegDFVol->Elt(i) * pInfLegDFVol->Elt(j)	* AvgFloatIRIRCorrel->Elt(i,j));
	sqrVolInfLeg /= FloatLegPV*FloatLegPV;

	/// Computes 
	double covarianceTerm = 0.0;
	for( i=0; i<pFwdCPI->size(); ++i )
		for ( j=0; j<pFixLegDF->size(); ++j )
			covarianceTerm += iTerms[i] * fixLegTerms[j] * (
				  pVol_ZC->Elt(i) * pFixLegDFVol->Elt(j)		* AvgFloatFixCPIIRCorrel->Elt(i,j)
				+ pInfLegDFVol->Elt(i) * pFixLegDFVol->Elt(j)	* AvgFloatFixIRIRCorrel ->Elt(i,j) );
	covarianceTerm /= FloatLegPV*FixedLegAnnuity;

	double vol = sqrVolFixLeg + sqrVolInfLeg - 2. * covarianceTerm;

	/// test before sqrt!
	if( vol < K_NEW_DOUBLE_TOL )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": vol is not positive!");
	vol	 = sqrt(vol);

	return (vol);
}



////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_Std
///	Routine: GenerateSwopVolCurve
///	Returns: ARM_VolCurve
///	Action : general function to compute the volatility of 
///			inf swaption from the one of year on year
////////////////////////////////////////////////////

ARM_VolCurve* ARM_InfVolComputation_Producer_Std::GenerateSwopVolCurve( 
	const ARM_Date&			asOfDate,
	ARM_GP_Vector*			tenors,
	ARM_GP_Vector*			expiries )
{
	/// get the default data if missing
	GetDefaultDataIfMissing( tenors, expiries );
	
	ARM_Matrix* volatilities = new ARM_Matrix( expiries->size(), tenors->size() );
	ARM_AutoCleaner< ARM_Matrix > HoldVolatilities( volatilities );

	/// get default of the currency of the IR market
	/// it is assumed that the convention on the CPI frequency is K_ANNUAL
	double CPIFreq		= K_ANNUAL;
	double CPIPeriod	= 1.0/CPIFreq;
	ARM_Date numDate,denomDate;
	double asOf			= itsInfBSModel->GetZeroCurve()->GetAsOfDate().GetJulian();
	double numYearFrac,denomYearFrac;
	string CPIIndexName = itsInfBSModel->GetInfFwdCurv()->GetInfIdxName();

	ARM_InfIdx* CPIIndex= itsInfBSModel->GetInfFwdCurv()->GetInfIdx();
	ARM_AutoCleaner<ARM_InfIdx> HoldCPIIndex( CPIIndex );
	long dailyInterpType= InfData::GetDailyInterpolation(CPIIndexName.c_str());
	double lookupStrike = 0.0;
	double ATMStrikeDefault = 0.0;
	double GAP = -91.0;

	double fixIRFreq	= CPIFreq;
	double fixIRPeriod	= 1.0/fixIRFreq;
	int fixIRDayCount	= itsInfBSModel->GetZeroCurve()->GetCurrencyUnit()->GetFixedDayCount();
	double IRFactor		= fixIRDayCount == KACTUAL_360? 365.0/360.0 : 1;
	double fixIRperiodDC= fixIRPeriod*IRFactor;

	ARM_CorrelMatrix* CPIIRCorrelMat	= itsInfBSModel->GetInfIRCorrelMatrix(CPIIndex,CPIIndex,  "INF/IR" );
	ARM_CorrelMatrix* CPICPICorrelMat	= itsInfBSModel->GetInfInfCorrelMatrix(CPIIndex,CPIIndex, "INF/INF"  );

	char* calendar= "INF";

	double longestTenor			= tenors->Elt(tenors->size()-1);
	size_t IRDataSize			= longestTenor*fixIRFreq;

	double timeToStart	= ARM_SparseVolCube::spotTime;
	ARM_Date lastKnownDate	= itsInfBSModel->GetVolatility()->GetLastKnownDate();
	int dayCount		= itsInfBSModel->GetInfFwdCurv()->GetMonthlyInterpType();

	ARM_GP_Vector* DF			= new ARM_GP_Vector(IRDataSize+1);
	ARM_AutoCleaner< ARM_GP_Vector > HoldFixLegDF( DF );

	ARM_GP_Vector* ForwardVol		= new ARM_GP_Vector(IRDataSize+1);
	ARM_AutoCleaner< ARM_GP_Vector > HoldFixLegDFVol( ForwardVol );

	ARM_GP_Vector* fixLegIRYearFrac	= new ARM_GP_Vector(IRDataSize,fixIRperiodDC);
		ARM_AutoCleaner< ARM_GP_Vector > HoldFixLegIRYearFrac( fixLegIRYearFrac );

	ARM_GP_Vector* ProgressiveFixLeg = new ARM_GP_Vector(IRDataSize+1);
	ARM_AutoCleaner< ARM_GP_Vector > HoldProgressiveFixLeg( ProgressiveFixLeg );

	ARM_GP_Vector* payTimes = new ARM_GP_Vector(IRDataSize+1);
	ARM_AutoCleaner< ARM_GP_Vector > HoldpayTimes( payTimes );

			/// 1) computes the volatility using the relationship
	for( size_t i=0; i<expiries->size(); ++i )
	{
		/// computes all the IR data for the longest tenors
		double payDate_0 = (*expiries)[i];
		payTimes->Elt(0) = payDate_0;
		DF->Elt(0) = itsInfBSModel->GetDiscountCurve()->DiscountPrice(payDate_0);
		ForwardVol->Elt(0) = itsInfBSModel->GetIRModel()->GetVolatility()->ComputeVolatility(payDate_0, ATMStrikeDefault,fixIRPeriod )/CC_NS( ARM_Constants, volBase );
		ProgressiveFixLeg->Elt(0) = 0.0 ;
		
		int j,k,l;
		for( k=1; k<IRDataSize+1; ++k )
		{
			double payDate_k= payDate_0 + k*fixIRPeriod;
			payTimes->Elt(k) = payDate_k;
			DF->Elt(k)	= itsInfBSModel->GetDiscountCurve()->DiscountPrice(payDate_k);
			ForwardVol->Elt(k) = itsInfBSModel->GetIRModel()->GetVolatility()->ComputeVolatility(payDate_k, ATMStrikeDefault,fixIRPeriod )/CC_NS( ARM_Constants, volBase );
			ProgressiveFixLeg->Elt(k) = ProgressiveFixLeg->Elt(k-1)+fixIRperiodDC*DF->Elt(k);
		}

		/// computes all the data for the CPIs
		size_t CPIDataSize				= longestTenor*CPIFreq;
		
		ARM_GP_Vector* CPIFwdRatio	= new ARM_GP_Vector(CPIDataSize+1);
		ARM_AutoCleaner< ARM_GP_Vector > HoldCPIFwdCPIRatio( CPIFwdRatio);
		
		ARM_GP_Vector* CPIVol_YoY		= new ARM_GP_Vector(CPIDataSize+1);
		ARM_AutoCleaner< ARM_GP_Vector > HoldCPIVol_YoY( CPIVol_YoY );
		
		ARM_GP_Vector* ProgressiveFloatLeg = new ARM_GP_Vector(CPIDataSize+1);
		ARM_AutoCleaner< ARM_GP_Vector > HoldProgressiveFloatLeg( ProgressiveFloatLeg );

		ARM_GP_Vector* publishLagDateVector = new ARM_GP_Vector(CPIDataSize+1);
		ARM_AutoCleaner< ARM_GP_Vector > HoldPublishLagDateVector( publishLagDateVector );
	
		denomYearFrac			= payTimes->Elt(0);
		denomDate				= ARM_Date(asOf+ denomYearFrac* K_YEAR_LEN);
		denomDate.GapBusinessDay(GAP, calendar);
		if( (denomDate.GetJulian() < asOf) )
		{
			denomDate = ARM_Date(asOf);			
		}
		double denomVolLookup	= CountYearsWithoutException( dayCount, lastKnownDate, denomDate );

		ARM_Date publishLagDate = itsInfBSModel->GetModelDateWPublishLag( denomDate, CPIIndex );
		publishLagDateVector->Elt(0)= CountYearsWithoutException( dayCount, asOf, publishLagDate);

		ProgressiveFloatLeg->Elt(0) =0.0;		
		
		for( l=1; l<(CPIDataSize+1); ++l )
		{
			numYearFrac				= payTimes->Elt(l);			
			numDate					= ARM_Date(asOf+  numYearFrac * K_YEAR_LEN);
			numDate.GapBusinessDay(GAP, calendar);
			if( (numDate.GetJulian() < asOf) )
			{
				numDate = ARM_Date(asOf);
			}		

			ARM_Date payDate_l = ARM_Date(asOf+  numYearFrac * K_YEAR_LEN);

			double numVolLookup	= CountYearsWithoutException( dayCount, lastKnownDate, numDate );
			
			CPIFwdRatio->Elt(l)	= itsInfBSModel->FwdCPIRatio( numDate, denomDate, payDate_l, 1.0, 0.0, dailyInterpType, GETDEFAULTVALUE, CPIIndex )[0];
			CPIVol_YoY->Elt(l)	= itsInfBSModel->GetVolatility()->ComputeVolatility( denomVolLookup, lookupStrike, CPIPeriod )/ CC_NS( ARM_Constants, volBase );
			
			
			ProgressiveFloatLeg->Elt(l) = ProgressiveFloatLeg->Elt(l-1)+CPIFwdRatio->Elt(l)*DF->Elt(l);
			denomYearFrac			= numYearFrac;
			denomDate = numDate;
			denomVolLookup = numVolLookup; 

			publishLagDate = itsInfBSModel->GetModelDateWPublishLag( denomDate, CPIIndex );
			publishLagDateVector->Elt(l)= CountYearsWithoutException( dayCount, asOf, publishLagDate);
		}

		ARM_Matrix* AvgCPICPICorrel		= new ARM_Matrix(CPIDataSize,CPIDataSize,1.0);
		ARM_AutoCleaner< ARM_Matrix > HoldAvgCPICPICorrel( AvgCPICPICorrel );

		ARM_Matrix* AvgCPIIRCorrel		= new ARM_Matrix(CPIDataSize,CPIDataSize,1.0);
		ARM_AutoCleaner< ARM_Matrix > HoldAvgCPIIRCorrel( AvgCPIIRCorrel );

		for( l=1; l<(CPIDataSize+1); ++l )
		{
			double T_l			= publishLagDateVector->Elt(l);
			for( k=1; k<IRDataSize+1; ++k )
			{
				double T_k			= publishLagDateVector->Elt(k);
				double T_kmoins1	= publishLagDateVector->Elt(k-1);

				//Correl INF/INF
				AvgCPICPICorrel->Elt(l-1,k-1)	= CPICPICorrelMat->ComputeCorrelData(T_l,T_k)/CC_NS(ARM_Constants,correlBase);
					
				//Correl INF/IR
				AvgCPIIRCorrel->Elt(l-1,k-1)	= CPIIRCorrelMat->ComputeCorrelData(T_l,T_kmoins1)/CC_NS(ARM_Constants,correlBase);	
			}
		}

		/*double maturityFactor = 1/payDate_0;
		double sqrtMaturityFactor = sqrt(maturityFactor);*/
		/// start calling the function to compute volatility with correct data
		for(j=0; j<tenors->size(); ++j )
		{
			/// need to correct bond vol to be really bond vol and not forward volatility
			/// IR data
			double var =0.0;
			double alpha_l, alpha_k, beta_l, beta_k;
			double volCPI_l, volCPI_k, volIR_l, volIR_k;
			double correl_INF_IR, correl_INF_INF;
			
			size_t TenorSize_j = tenors->Elt(j) * CPIFreq+1;
			double fixPonderation_j  =	1.0/(TenorSize_j-1);
			double floatleg_j	= ProgressiveFloatLeg->Elt(TenorSize_j-1);
			double fixleg_j		= ProgressiveFixLeg->Elt(TenorSize_j-1);

			for(l=1;l<TenorSize_j;l++)
			{
				double floatleg_l		=	floatleg_j-ProgressiveFloatLeg->Elt(l-1);
				double fixleg_l		=	fixleg_j-ProgressiveFixLeg->Elt(l-1);
				alpha_l			=	DF->Elt(l)*CPIFwdRatio->Elt(l)/floatleg_j;
				beta_l			=	(DF->Elt(l-1)-DF->Elt(l))/DF->Elt(l-1)*(fixleg_l/fixleg_j-floatleg_l/floatleg_j);
				volCPI_l		=	CPIVol_YoY->Elt(l);
				volIR_l			=	ForwardVol->Elt(l-1);

				double maturityFactor_T_l		= sqrt(1.0/publishLagDateVector->Elt(l));
							
				for(k=1;k<TenorSize_j;k++)
				{
					double floatleg_k		=	floatleg_j-ProgressiveFloatLeg->Elt(k-1);
					double fixleg_k		=	fixleg_j-ProgressiveFixLeg->Elt(k-1);
					alpha_k	=	DF->Elt(k)*CPIFwdRatio->Elt(k)/floatleg_j;
					beta_k	=	(DF->Elt(k-1)-DF->Elt(k))/DF->Elt(k-1)*(fixleg_k/fixleg_j-floatleg_k/floatleg_j);
					volCPI_k	= CPIVol_YoY->Elt(k);;
					volIR_k		= ForwardVol->Elt(k-1);
										
					double maturityFactor_T_k		= sqrt(1.0/publishLagDateVector->Elt(k));
				
					//Correl INF/INF
					correl_INF_INF		= AvgCPICPICorrel->Elt(l-1,k-1);
						
					//Correl INF/IR
					correl_INF_IR		= AvgCPIIRCorrel->Elt(l-1,k-1);
					
					var+=maturityFactor_T_l*maturityFactor_T_k*(alpha_l*alpha_k*correl_INF_INF*volCPI_l*volCPI_k)+
						maturityFactor_T_l*(alpha_l*beta_k*volCPI_l*volIR_k*correl_INF_IR)+
						(beta_l*beta_k*volIR_l*volIR_k);
				}
			}

			if( var < K_NEW_DOUBLE_TOL )
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": vol is not positive!");
			
			var	 = sqrt(var)*CC_NS(ARM_Constants,volBase);
			volatilities->Elt(i,j)= var;
		}
	}

	/// 2) set the volatility value to the matrix
	ARM_Vector* expiriesVec	= To_pARM_Vector( expiries );
	delete expiries;
	ARM_Vector* tenorsVec	= To_pARM_Vector( tenors );
	delete tenors;

	ARM_VolCurve* volCurve = new ARM_VolLInterpol(asOfDate, expiriesVec, tenorsVec, volatilities, K_STK_TYPE_PRICE, K_ATMF_VOL );
	volCurve->SetVolatilities(volatilities);
	HoldVolatilities.Release();
	
	return volCurve;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_Std
///	Routine: GenerateSwopVolCube
///	Returns: ARM_VolCurve
///	Action : general function to compute the volatility of 
///			inf swaption from the one of year on year
////////////////////////////////////////////////////

ARM_VolCube* ARM_InfVolComputation_Producer_Std::GenerateSwopVolCube( 
	const ARM_Date&			asOfDate,
	ARM_GP_Vector*			tenors,
	ARM_GP_Vector*			expiries,
	ARM_GP_Vector*			smiledTenors,
	ARM_GP_Vector*			strikes)
{
	GetDefaultDataIfMissing( tenors, expiries ); 
	ARM_GP_Vector* tenorsCopy = (ARM_GP_Vector*) (tenors->Clone());
	ARM_GP_Vector* expiriesCopy = (ARM_GP_Vector*) (expiries->Clone());
	ARM_VolCurve* ATMVol =  ARM_InfVolComputation_Producer_Std::GenerateSwopVolCurve( asOfDate,tenorsCopy, expiriesCopy);
		
	vector<ARM_VolCurve*> SmiledVols;
	int nbSmiledVols = smiledTenors->size();
	SmiledVols.resize(nbSmiledVols);

	int expiriesSize = expiries->size();
	int strikesSize = strikes->size();

	////////////////////////Generate DAta for coeff relation Vol swap/ Vol FRA
	double CPIFreq		= K_ANNUAL;
	double CPIPeriod	= 1.0/CPIFreq;
	ARM_Date numDate,denomDate;
	double asOf			= itsInfBSModel->GetZeroCurve()->GetAsOfDate().GetJulian();
	double numYearFrac,denomYearFrac;
	string CPIIndexName = itsInfBSModel->GetInfFwdCurv()->GetInfIdxName();

	ARM_InfIdx* CPIIndex= itsInfBSModel->GetInfFwdCurv()->GetInfIdx();
	ARM_AutoCleaner<ARM_InfIdx> HoldCPIIndex( CPIIndex );
	long dailyInterpType= InfData::GetDailyInterpolation(CPIIndexName.c_str());
	double lookupStrike = 0.0;
	double ATMStrikeDefault = 0.0;
	double GAP = -91.0;

	double fixIRFreq	= CPIFreq;
	double fixIRPeriod	= 1.0/fixIRFreq;
	int fixIRDayCount	= itsInfBSModel->GetZeroCurve()->GetCurrencyUnit()->GetFixedDayCount();
	double IRFactor		= fixIRDayCount == KACTUAL_360? 365.0/360.0 : 1;
	double fixIRperiodDC= fixIRPeriod*IRFactor;

	ARM_CorrelMatrix* CPIIRCorrelMat	= itsInfBSModel->GetInfIRCorrelMatrix(CPIIndex,CPIIndex,  "INF/IR" );
	ARM_CorrelMatrix* CPICPICorrelMat	= itsInfBSModel->GetInfInfCorrelMatrix(CPIIndex,CPIIndex, "INF/INF"  );

	char* calendar= "INF";

	CC_NS(std,sort)(smiledTenors->begin(),smiledTenors->end()); 
	double longestTenor			= smiledTenors->Elt(smiledTenors->size()-1);
	size_t IRDataSize			= longestTenor*fixIRFreq;

	double timeToStart	= ARM_SparseVolCube::spotTime;
	ARM_Date lastKnownDate	= itsInfBSModel->GetVolatility()->GetLastKnownDate();
	int dayCount		= itsInfBSModel->GetInfFwdCurv()->GetMonthlyInterpType();

	ARM_Matrix* DF			= new ARM_Matrix(expiriesSize, IRDataSize+1);
	ARM_AutoCleaner< ARM_Matrix > HoldFixLegDF( DF );

	ARM_Matrix* ForwardVol		= new ARM_Matrix(expiriesSize,IRDataSize+1);
	ARM_AutoCleaner< ARM_Matrix > HoldFixLegDFVol( ForwardVol );

	
	ARM_Matrix* ProgressiveFixLeg = new ARM_Matrix(expiriesSize, IRDataSize+1);
	ARM_AutoCleaner< ARM_Matrix > HoldProgressiveFixLeg( ProgressiveFixLeg );

	ARM_Matrix* payTimes = new ARM_Matrix(expiriesSize, IRDataSize+1);
	ARM_AutoCleaner< ARM_Matrix > HoldpayTimes( payTimes );

	size_t CPIDataSize				= longestTenor*CPIFreq;
		
	ARM_Matrix* CPIFwdRatio	= new ARM_Matrix(expiriesSize,CPIDataSize+1);
	ARM_AutoCleaner< ARM_Matrix > HoldCPIFwdCPIRatio( CPIFwdRatio);
	
	ARM_Matrix* CPIVol_YoY		= new ARM_Matrix(expiriesSize,CPIDataSize+1);
	ARM_AutoCleaner< ARM_Matrix > HoldCPIVol_YoY( CPIVol_YoY );
	
	ARM_Matrix* ProgressiveFloatLeg = new ARM_Matrix(expiriesSize,CPIDataSize+1);
	ARM_AutoCleaner< ARM_Matrix > HoldProgressiveFloatLeg( ProgressiveFloatLeg );

	ARM_Matrix* publishLagDateVector = new ARM_Matrix(expiriesSize, CPIDataSize+1);
	ARM_AutoCleaner< ARM_Matrix > HoldPublishLagDateVector( publishLagDateVector );

	ARM_Matrix* dateVolLookup = new ARM_Matrix(expiriesSize, CPIDataSize+1);
	ARM_AutoCleaner< ARM_Matrix > HolddateVolLookup( dateVolLookup );

			/// 1) computes the volatility using the relationship
	for( size_t i=0; i<expiriesSize; ++i )
	{
		/// computes all the IR data for the longest tenors
		double payDate_0 = (*expiries)[i];
		payTimes->Elt(i,0) = payDate_0;
		DF->Elt(i,0) = itsInfBSModel->GetDiscountCurve()->DiscountPrice(payDate_0);
		ForwardVol->Elt(i,0) = itsInfBSModel->GetIRModel()->GetVolatility()->ComputeVolatility(payDate_0, ATMStrikeDefault,fixIRPeriod )/CC_NS( ARM_Constants, volBase );
		ProgressiveFixLeg->Elt(i,0) = 0.0 ;
		
		int k,l;
		for( k=1; k<IRDataSize+1; ++k )
		{
			double payDate_k= payDate_0 + k*fixIRPeriod;
			payTimes->Elt(i,k) = payDate_k;
			DF->Elt(i,k)	= itsInfBSModel->GetDiscountCurve()->DiscountPrice(payDate_k);
			ForwardVol->Elt(i,k) = itsInfBSModel->GetIRModel()->GetVolatility()->ComputeVolatility(payDate_k, ATMStrikeDefault,fixIRPeriod )/CC_NS( ARM_Constants, volBase );
			ProgressiveFixLeg->Elt(i,k) = ProgressiveFixLeg->Elt(i,k-1)+fixIRperiodDC*DF->Elt(i,k);
		}

		/// computes all the data for the CPIs	
		denomYearFrac			= payTimes->Elt(i,0);
		denomDate				= ARM_Date(asOf+ denomYearFrac* K_YEAR_LEN);
		denomDate.GapBusinessDay(GAP, calendar);
		if( (denomDate.GetJulian() < asOf) )
		{
			denomDate = ARM_Date(asOf);			
		}
		double denomVolLookup	= CountYearsWithoutException( dayCount, lastKnownDate, denomDate );

		dateVolLookup->Elt(i,0) = denomVolLookup;

		ARM_Date publishLagDate = itsInfBSModel->GetModelDateWPublishLag( denomDate, CPIIndex );
		publishLagDateVector->Elt(i,0)= CountYearsWithoutException( dayCount, asOf, publishLagDate);

		ProgressiveFloatLeg->Elt(i,0) =0.0;		
		
		for( l=1; l<(CPIDataSize+1); ++l )
		{
			numYearFrac				= payTimes->Elt(i,l);			
			numDate					= ARM_Date(asOf+  numYearFrac * K_YEAR_LEN);
			numDate.GapBusinessDay(GAP, calendar);
			if( (numDate.GetJulian() < asOf) )
			{
				numDate = ARM_Date(asOf);
			}		

			ARM_Date payDate_l = ARM_Date(asOf+  numYearFrac * K_YEAR_LEN);

			double numVolLookup	= CountYearsWithoutException( dayCount, lastKnownDate, numDate );
			
			CPIFwdRatio->Elt(i,l)	= itsInfBSModel->FwdCPIRatio( numDate, denomDate, payDate_l, 1.0, 0.0, dailyInterpType, GETDEFAULTVALUE, CPIIndex )[0];
			CPIVol_YoY->Elt(i,l)	= itsInfBSModel->GetVolatility()->ComputeVolatility( denomVolLookup, lookupStrike, CPIPeriod )/ CC_NS( ARM_Constants, volBase );
			
			
			ProgressiveFloatLeg->Elt(i,l) = ProgressiveFloatLeg->Elt(i,l-1)+CPIFwdRatio->Elt(i,l)*DF->Elt(i,l);
			denomYearFrac			= numYearFrac;
			denomDate = numDate;
			denomVolLookup = numVolLookup; 

			publishLagDate = itsInfBSModel->GetModelDateWPublishLag( denomDate, CPIIndex );
			publishLagDateVector->Elt(i,l)= CountYearsWithoutException( dayCount, asOf, publishLagDate);
			dateVolLookup->Elt(i,l) = denomVolLookup;

		}		
	}
	
	ARM_Vector* expiriesVec = To_pARM_Vector(expiries);
	ARM_Vector* strikesVec = To_pARM_Vector(strikes);
	ARM_Vector* smiledTenorsVec = To_pARM_Vector(smiledTenors);

	
	//ARM_AutoCleaner< ARM_Matrix > HoldVolatilities( volatilities );

	for(int t=0; t<nbSmiledVols; t++)
	{
		ARM_Matrix* volatilities = new ARM_Matrix( expiriesSize, strikesSize, 0.0 );
		ARM_VolCurve* volCurve_t = new ARM_VolLInterpol(asOfDate,  (ARM_Vector*) (expiriesVec->Clone()), (ARM_Vector*) (strikesVec->Clone()), volatilities, K_STK_TYPE_PRICE, K_SMILE_VOL);
		volCurve_t->SetVolatilities(volatilities);
		//HoldVolatilities.Release();

		SmiledVols[t] = volCurve_t;
	}
	

	ARM_Matrix* CPISmiledVol = new ARM_Matrix(CPIDataSize+1, strikesSize);
	ARM_AutoCleaner< ARM_Matrix > HoldCPISmiledVol( CPISmiledVol );
	double ATM_Strike =0.0;

	ARM_GP_Vector* alphaVector	= new ARM_GP_Vector(CPIDataSize+1);
	ARM_AutoCleaner< ARM_GP_Vector > HoldalphaVector( alphaVector);

	ARM_GP_Vector* betaVector	= new ARM_GP_Vector(CPIDataSize+1);
	ARM_AutoCleaner< ARM_GP_Vector > HoldbetaVector( betaVector);

	ARM_Matrix* AvgCPICPICorrel		= new ARM_Matrix(CPIDataSize,CPIDataSize,1.0);
	ARM_AutoCleaner< ARM_Matrix > HoldAvgCPICPICorrel( AvgCPICPICorrel );

	ARM_Matrix* AvgCPIIRCorrel		= new ARM_Matrix(CPIDataSize,CPIDataSize,1.0);
	ARM_AutoCleaner< ARM_Matrix > HoldAvgCPIIRCorrel( AvgCPIIRCorrel );

	for( i=0; i<expiriesSize; ++i )
	{
		double payDate_0 = (*expiries)[i];			
		/// computes all the IR data for the longest tenors	
		
		for( int l=1; l<(CPIDataSize+1); ++l )
		{
			double T_l			= publishLagDateVector->Elt(i,l);
			for( int k=1; k<IRDataSize+1; ++k )
			{
				double T_k			= publishLagDateVector->Elt(i,k);
				double T_kmoins1	= publishLagDateVector->Elt(i,k-1);

				//Correl INF/INF
				AvgCPICPICorrel->Elt(l-1,k-1)	= CPICPICorrelMat->ComputeCorrelData(T_l,T_k)/CC_NS(ARM_Constants,correlBase);
					
				//Correl INF/IR
				AvgCPIIRCorrel->Elt(l-1,k-1)	= CPIIRCorrelMat->ComputeCorrelData(T_l,T_kmoins1)/CC_NS(ARM_Constants,correlBase);	
			}
		}
		
		for(int t=0; t<nbSmiledVols; t++)
		{
			double tenor_t = smiledTenorsVec->Elt(t);
			size_t TenorSize_t		= tenor_t*CPIFreq+1;
			double ATM_Vol_Tenor_t =  ATMVol->ComputeVolatility(payDate_0, ATM_Strike, tenor_t)/ CC_NS( ARM_Constants, volBase );
			double floatleg_t	= ProgressiveFloatLeg->Elt(i,TenorSize_t-1);
			double fixleg_t		= ProgressiveFixLeg->Elt(i,TenorSize_t-1);
			double Swap_Tenor_t = floatleg_t/fixleg_t;

			for( int l=1; l<(TenorSize_t); ++l )
			{
				double T_l			= publishLagDateVector->Elt(i,l);

				double ATM_Vol_l = CPIVol_YoY->Elt(i,l);
				double CPIFwdRatio_l = CPIFwdRatio->Elt(i,l);
				double coeff_l = sqrt(1/T_l)*ATM_Vol_l/ATM_Vol_Tenor_t;
				double denomVolLookup_l = dateVolLookup->Elt(i,l-1);
				for( int p=0; p<strikesSize; ++p )
				{
					double Strike_l_p = coeff_l*(*strikes)[p];
					CPISmiledVol->Elt(l,p)= itsInfBSModel->GetVolatility()->ComputeVolatility( denomVolLookup_l, Strike_l_p, CPIPeriod)/ CC_NS( ARM_Constants, volBase );
				}

				double floatleg_l		=	floatleg_t-ProgressiveFloatLeg->Elt(i,l-1);
				double fixleg_l		=	fixleg_t-ProgressiveFixLeg->Elt(i,l-1);
				double alpha_l			=	DF->Elt(i,l)*CPIFwdRatio->Elt(i,l)/floatleg_t;
				alphaVector->Elt(l) = alpha_l;
				double beta_l	=	(DF->Elt(i,l-1)-DF->Elt(i,l))/DF->Elt(i,l-1)*(fixleg_l/fixleg_t-floatleg_l/floatleg_t);
				betaVector->Elt(l) = beta_l;
			}

			
			
			double volCPI_l, volCPI_k, volIR_l, volIR_k;
			double correl_INF_IR, correl_INF_INF;

			double maturityFactor_T_k, maturityFactor_T_l, alpha_l, alpha_k, beta_l, beta_k;

			double var;
			for( int p=0; p<strikesSize; ++p )
			{
				var =0.0;
				for(l=1;l<TenorSize_t;l++)
				{
					alpha_l			=	alphaVector->Elt(l);
					beta_l			=	betaVector->Elt(l);
					volCPI_l		=	CPISmiledVol->Elt(l,p);
					volIR_l			=	ForwardVol->Elt(i,l-1);

					maturityFactor_T_l		= sqrt(1.0/publishLagDateVector->Elt(i,l));
									
					for(int k=1;k<TenorSize_t;k++)
					{
						alpha_k	=	alphaVector->Elt(k);
						beta_k	=	betaVector->Elt(k);
						volCPI_k	= CPISmiledVol->Elt(k,p);
						volIR_k		= ForwardVol->Elt(i,k-1);

						maturityFactor_T_k		= sqrt(1.0/publishLagDateVector->Elt(i,k));
																	
						//Correl INF/INF
						correl_INF_INF		= AvgCPICPICorrel->Elt(l-1,k-1);
							
						//Correl INF/IR
						correl_INF_IR		= AvgCPIIRCorrel->Elt(l-1,k-1);

						var+=maturityFactor_T_l*maturityFactor_T_k*(alpha_l*alpha_k*correl_INF_INF*volCPI_l*volCPI_k)+
							maturityFactor_T_l*(alpha_l*beta_k*volCPI_l*volIR_k*correl_INF_IR)+
							(beta_l*beta_k*volIR_l*volIR_k);
					}
				}

				if( var < K_NEW_DOUBLE_TOL )
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": vol is not positive!");
			
								

				var	 = sqrt(var)*CC_NS(ARM_Constants,volBase);
				
				
				SmiledVols[t]->GetVolatilities()->Elt(i,p)=var-ATM_Vol_Tenor_t*CC_NS(ARM_Constants,volBase);
			}		
		}				
	}	

	delete expiries;
	delete tenors;
	delete strikes;
	delete smiledTenors;	

	ARM_VolCube* result = new ARM_VolCube(ATMVol,
										  &SmiledVols[0],
										  nbSmiledVols,
										  smiledTenorsVec,
                                          K_ATMF_VOL);//all attribute are cloned we have to delete them

	delete expiriesVec;
	delete strikesVec;
	delete smiledTenorsVec;
	delete ATMVol;

	
	vector<ARM_VolCurve*>::iterator it;
        
   for (it = SmiledVols.begin(); it < SmiledVols.end(); it++)
   {
       delete (*it);
   }
	
	return result;
}
/*
ARM_VolCube* ARM_InfVolComputation_Producer_Std::GenerateSwopVolCube( 
	const ARM_Date&			asOfDate,
	ARM_GP_Vector*			tenors,
	ARM_GP_Vector*			expiries,
	ARM_GP_Vector*			smiledTenors,
	ARM_GP_Vector*			strikes)
{
	GetDefaultDataIfMissing( tenors, expiries ); 
	ARM_GP_Vector* tenorsCopy = (ARM_GP_Vector*) (tenors->Clone());
	ARM_GP_Vector* expiriesCopy = (ARM_GP_Vector*) (expiries->Clone());
	ARM_VolCurve* ATMVol =  ARM_InfVolComputation_Producer_Std::GenerateSwopVolCurve( asOfDate,tenorsCopy, expiriesCopy);
		
	vector<ARM_VolCurve*> SmiledVols;
	int nbSmiledVols = smiledTenors->size();
	SmiledVols.resize(nbSmiledVols);

	int expiriesSize = expiries->size();
	int strikesSize = strikes->size();

	////////////////////////Generate DAta for coeff relation Vol swap/ Vol FRA
	double CPIFreq		= K_ANNUAL;
	double CPIPeriod	= 1.0/CPIFreq;
	ARM_Date numDate,denomDate;
	double asOf			= itsInfBSModel->GetZeroCurve()->GetAsOfDate().GetJulian();
	double numYearFrac,denomYearFrac;
	string CPIIndexName = itsInfBSModel->GetInfFwdCurv()->GetInfIdxName();

	ARM_InfIdx* CPIIndex= itsInfBSModel->GetInfFwdCurv()->GetInfIdx();
	ARM_AutoCleaner<ARM_InfIdx> HoldCPIIndex( CPIIndex );
	long dailyInterpType= InfData::GetDailyInterpolation(CPIIndexName.c_str());
	double lookupStrike = 0.0;
	double ATMStrikeDefault = 0.0;
	double GAP = -91.0;

	double fixIRFreq	= CPIFreq;
	double fixIRPeriod	= 1.0/fixIRFreq;
	int fixIRDayCount	= itsInfBSModel->GetZeroCurve()->GetCurrencyUnit()->GetFixedDayCount();
	double IRFactor		= fixIRDayCount == KACTUAL_360? 365.0/360.0 : 1;
	double fixIRperiodDC= fixIRPeriod*IRFactor;

	ARM_CorrelMatrix* CPIIRCorrelMat	= itsInfBSModel->GetInfIRCorrelMatrix(CPIIndex,CPIIndex,  "INF/IR" );
	ARM_CorrelMatrix* CPICPICorrelMat	= itsInfBSModel->GetInfInfCorrelMatrix(CPIIndex,CPIIndex, "INF/INF"  );

	char* calendar= "INF";

	CC_NS(std,sort)(smiledTenors->begin(),smiledTenors->end()); 
	double longestTenor			= smiledTenors->Elt(smiledTenors->size()-1);
	size_t IRDataSize			= longestTenor*fixIRFreq;

	double timeToStart	= ARM_SparseVolCube::spotTime;
	ARM_Date lastKnownDate	= itsInfBSModel->GetVolatility()->GetLastKnownDate();
	int dayCount		= itsInfBSModel->GetInfFwdCurv()->GetMonthlyInterpType();

	ARM_Matrix* DF			= new ARM_Matrix(expiriesSize, IRDataSize+1);
	ARM_AutoCleaner< ARM_Matrix > HoldFixLegDF( DF );

	ARM_Matrix* ForwardVol		= new ARM_Matrix(expiriesSize,IRDataSize+1);
	ARM_AutoCleaner< ARM_Matrix > HoldFixLegDFVol( ForwardVol );

	
	ARM_Matrix* ProgressiveFixLeg = new ARM_Matrix(expiriesSize, IRDataSize+1);
	ARM_AutoCleaner< ARM_Matrix > HoldProgressiveFixLeg( ProgressiveFixLeg );

	ARM_Matrix* payTimes = new ARM_Matrix(expiriesSize, IRDataSize+1);
	ARM_AutoCleaner< ARM_Matrix > HoldpayTimes( payTimes );

	size_t CPIDataSize				= longestTenor*CPIFreq;
		
	ARM_Matrix* CPIFwdRatio	= new ARM_Matrix(expiriesSize,CPIDataSize+1);
	ARM_AutoCleaner< ARM_Matrix > HoldCPIFwdCPIRatio( CPIFwdRatio);
	
	ARM_Matrix* CPIVol_YoY		= new ARM_Matrix(expiriesSize,CPIDataSize+1);
	ARM_AutoCleaner< ARM_Matrix > HoldCPIVol_YoY( CPIVol_YoY );
	
	ARM_Matrix* ProgressiveFloatLeg = new ARM_Matrix(expiriesSize,CPIDataSize+1);
	ARM_AutoCleaner< ARM_Matrix > HoldProgressiveFloatLeg( ProgressiveFloatLeg );

	ARM_Matrix* publishLagDateVector = new ARM_Matrix(expiriesSize, CPIDataSize+1);
	ARM_AutoCleaner< ARM_Matrix > HoldPublishLagDateVector( publishLagDateVector );

	ARM_Matrix* dateVolLookup = new ARM_Matrix(expiriesSize, CPIDataSize+1);
	ARM_AutoCleaner< ARM_Matrix > HolddateVolLookup( dateVolLookup );

			/// 1) computes the volatility using the relationship
	for( size_t i=0; i<expiriesSize; ++i )
	{
		/// computes all the IR data for the longest tenors
		double payDate_0 = (*expiries)[i];
		payTimes->Elt(i,0) = payDate_0;
		DF->Elt(i,0) = itsInfBSModel->GetDiscountCurve()->DiscountPrice(payDate_0);
		ForwardVol->Elt(i,0) = itsInfBSModel->GetIRModel()->GetVolatility()->ComputeVolatility(payDate_0, ATMStrikeDefault,fixIRPeriod )/CC_NS( ARM_Constants, volBase );
		ProgressiveFixLeg->Elt(i,0) = 0.0 ;
		
		int k,l;
		for( k=1; k<IRDataSize+1; ++k )
		{
			double payDate_k= payDate_0 + k*fixIRPeriod;
			payTimes->Elt(i,k) = payDate_k;
			DF->Elt(i,k)	= itsInfBSModel->GetDiscountCurve()->DiscountPrice(payDate_k);
			ForwardVol->Elt(i,k) = itsInfBSModel->GetIRModel()->GetVolatility()->ComputeVolatility(payDate_k, ATMStrikeDefault,fixIRPeriod )/CC_NS( ARM_Constants, volBase );
			ProgressiveFixLeg->Elt(i,k) = ProgressiveFixLeg->Elt(i,k-1)+fixIRperiodDC*DF->Elt(i,k);
		}

		/// computes all the data for the CPIs	
		denomYearFrac			= payTimes->Elt(i,0);
		denomDate				= ARM_Date(asOf+ denomYearFrac* K_YEAR_LEN);
		denomDate.GapBusinessDay(GAP, calendar);
		if( (denomDate.GetJulian() < asOf) )
		{
			denomDate = ARM_Date(asOf);			
		}
		double denomVolLookup	= CountYearsWithoutException( dayCount, lastKnownDate, denomDate );

		dateVolLookup->Elt(i,0) = denomVolLookup;

		ARM_Date publishLagDate = itsInfBSModel->GetModelDateWPublishLag( denomDate, CPIIndex );
		publishLagDateVector->Elt(i,0)= CountYearsWithoutException( dayCount, asOf, publishLagDate);

		ProgressiveFloatLeg->Elt(i,0) =0.0;		
		
		for( l=1; l<(CPIDataSize+1); ++l )
		{
			numYearFrac				= payTimes->Elt(i,l);			
			numDate					= ARM_Date(asOf+  numYearFrac * K_YEAR_LEN);
			numDate.GapBusinessDay(GAP, calendar);
			if( (numDate.GetJulian() < asOf) )
			{
				numDate = ARM_Date(asOf);
			}		

			ARM_Date payDate_l = ARM_Date(asOf+  numYearFrac * K_YEAR_LEN);

			double numVolLookup	= CountYearsWithoutException( dayCount, lastKnownDate, numDate );
			
			CPIFwdRatio->Elt(i,l)	= itsInfBSModel->FwdCPIRatio( numDate, denomDate, payDate_l, 1.0, 0.0, dailyInterpType, GETDEFAULTVALUE, CPIIndex )[0];
			CPIVol_YoY->Elt(i,l)	= itsInfBSModel->GetVolatility()->ComputeVolatility( denomVolLookup, lookupStrike, CPIPeriod )/ CC_NS( ARM_Constants, volBase );
			
			
			ProgressiveFloatLeg->Elt(i,l) = ProgressiveFloatLeg->Elt(i,l-1)+CPIFwdRatio->Elt(i,l)*DF->Elt(i,l);
			denomYearFrac			= numYearFrac;
			denomDate = numDate;
			denomVolLookup = numVolLookup; 

			publishLagDate = itsInfBSModel->GetModelDateWPublishLag( denomDate, CPIIndex );
			publishLagDateVector->Elt(i,l)= CountYearsWithoutException( dayCount, asOf, publishLagDate);
			dateVolLookup->Elt(i,l) = denomVolLookup;

		}		
	}
	
	ARM_Vector* expiriesVec = To_pARM_Vector(expiries);
	ARM_Vector* strikesVec = To_pARM_Vector(strikes);
	ARM_Vector* smiledTenorsVec = To_pARM_Vector(smiledTenors);

	
	//ARM_AutoCleaner< ARM_Matrix > HoldVolatilities( volatilities );

	for(int t=0; t<nbSmiledVols; t++)
	{
		ARM_Matrix* volatilities = new ARM_Matrix( expiriesSize, strikesSize, 0.0 );
		ARM_VolCurve* volCurve_t = new ARM_VolLInterpol(asOfDate,  (ARM_Vector*) (expiriesVec->Clone()), (ARM_Vector*) (strikesVec->Clone()), volatilities, K_STK_TYPE_PRICE, K_SMILE_VOL);
		volCurve_t->SetVolatilities(volatilities);
		//HoldVolatilities.Release();

		SmiledVols[t] = volCurve_t;
	}
	

	ARM_Matrix* CPISmiledVol = new ARM_Matrix(CPIDataSize+1, strikesSize);
	ARM_AutoCleaner< ARM_Matrix > HoldCPISmiledVol( CPISmiledVol );
	double ATM_Strike =0.0;

	ARM_GP_Vector* alphaVector	= new ARM_GP_Vector(CPIDataSize+1);
	ARM_AutoCleaner< ARM_GP_Vector > HoldalphaVector( alphaVector);

	ARM_GP_Vector* betaVector	= new ARM_GP_Vector(CPIDataSize+1);
	ARM_AutoCleaner< ARM_GP_Vector > HoldbetaVector( betaVector);

	ARM_Matrix* AvgCPICPICorrel		= new ARM_Matrix(CPIDataSize,CPIDataSize,1.0);
	ARM_AutoCleaner< ARM_Matrix > HoldAvgCPICPICorrel( AvgCPICPICorrel );

	ARM_Matrix* AvgCPIIRCorrel		= new ARM_Matrix(CPIDataSize,CPIDataSize,1.0);
	ARM_AutoCleaner< ARM_Matrix > HoldAvgCPIIRCorrel( AvgCPIIRCorrel );


	FILE* f=fopen("c:\\temp\\dumpSWOPTVOL.txt","w");
	fclose(f);
	int conteur =0;
	for( i=0; i<expiriesSize; ++i )
	{
		double payDate_0 = (*expiries)[i];			
		/// computes all the IR data for the longest tenors	
		
		
		for( int l=1; l<(CPIDataSize+1); ++l )
		{
			double T_l			= publishLagDateVector->Elt(i,l);
			for( int k=1; k<IRDataSize+1; ++k )
			{
				double T_k			= publishLagDateVector->Elt(i,k);
				double T_kmoins1	= publishLagDateVector->Elt(i,k-1);

				//Correl INF/INF
				AvgCPICPICorrel->Elt(l-1,k-1)	= CPICPICorrelMat->ComputeCorrelData(T_l,T_k)/CC_NS(ARM_Constants,correlBase);
					
				//Correl INF/IR
				AvgCPIIRCorrel->Elt(l-1,k-1)	= CPIIRCorrelMat->ComputeCorrelData(T_l,T_kmoins1)/CC_NS(ARM_Constants,correlBase);	
			}
		}
		
		for(int t=0; t<nbSmiledVols; t++)
		{
			double tenor_t = smiledTenorsVec->Elt(t);
			size_t TenorSize_t		= tenor_t*CPIFreq+1;
			double ATM_Vol_Tenor_t =  ATMVol->ComputeVolatility(payDate_0, ATM_Strike, tenor_t)/ CC_NS( ARM_Constants, volBase );
			double floatleg_t	= ProgressiveFloatLeg->Elt(i,TenorSize_t-1);
			double fixleg_t		= ProgressiveFixLeg->Elt(i,TenorSize_t-1);
			double Swap_Tenor_t = floatleg_t/fixleg_t;

			for( int l=1; l<(TenorSize_t); ++l )
			{
				double T_l			= publishLagDateVector->Elt(i,l);

				double ATM_Vol_l = CPIVol_YoY->Elt(i,l);
				double CPIFwdRatio_l = CPIFwdRatio->Elt(i,l);
				double coeff_l = sqrt(1/T_l)*ATM_Vol_l/ATM_Vol_Tenor_t;
				double denomVolLookup_l = dateVolLookup->Elt(i,l-1);
				for( int p=0; p<strikesSize; ++p )
				{
					double Strike_l_p = coeff_l*(*strikes)[p];

					ARM_VolCube* volcube = dynamic_cast<ARM_VolCube*> (itsInfBSModel->GetVolatility());
					CPISmiledVol->Elt(l,p)= volcube->VolatilityFunction( denomVolLookup_l, Strike_l_p, CPIPeriod, true )/ CC_NS( ARM_Constants, volBase );
					f=fopen("c:\\temp\\dumpSWOPTVOL.txt","a");
					fprintf(f," CreateVolCube");
					conteur++;
					fprintf(f," conteur=%3d\t \n", conteur);
					fprintf(f," expiry=%7.10lf\t tenor=%7.10lf\t sizeTenor_t=%3d\t sizeTenor_l=%3d\t strikesSize=%3d\t strikesSize_p=%3d\t  CPIPeriod=%7.10lf\t\n", payDate_0, tenor_t, TenorSize_t, l,strikesSize,p, CPIPeriod);
					fprintf(f," strike_p=%7.10lf\t coeff_l=%7.10f\t Strike_l_p=%7.10lf\t volDate=%7.10lf\t  smiledVolYoY=%7.5lf\t  \n", (*strikes)[p], coeff_l,Strike_l_p, denomVolLookup_l,CPISmiledVol->Elt(l,p) );
					fclose(f);

				}

				double floatleg_l		=	floatleg_t-ProgressiveFloatLeg->Elt(i,l-1);
				double fixleg_l		=	fixleg_t-ProgressiveFixLeg->Elt(i,l-1);
				double alpha_l			=	DF->Elt(i,l)*CPIFwdRatio->Elt(i,l)/floatleg_t;
				alphaVector->Elt(l) = alpha_l;
				double beta_l	=	(DF->Elt(i,l-1)-DF->Elt(i,l))/DF->Elt(i,l-1)*(fixleg_l/fixleg_t-floatleg_l/floatleg_t);
				betaVector->Elt(l) = beta_l;
			}

			
			
			double volCPI_l, volCPI_k, volIR_l, volIR_k;
			double correl_INF_IR, correl_INF_INF;

			double maturityFactor_T_k, maturityFactor_T_l, alpha_l, alpha_k, beta_l, beta_k;

			double var;
			for( int p=0; p<strikesSize; ++p )
			{
				var =0.0;
				for(l=1;l<TenorSize_t;l++)
				{
					alpha_l			=	alphaVector->Elt(l);
					beta_l			=	betaVector->Elt(l);
					volCPI_l		=	CPISmiledVol->Elt(l,p);
					volIR_l			=	ForwardVol->Elt(i,l-1);

					maturityFactor_T_l		= sqrt(1.0/publishLagDateVector->Elt(i,l));
									
					for(int k=1;k<TenorSize_t;k++)
					{
						alpha_k	=	alphaVector->Elt(k);
						beta_k	=	betaVector->Elt(k);
						volCPI_k	= CPISmiledVol->Elt(k,p);
						volIR_k		= ForwardVol->Elt(i,k-1);

						maturityFactor_T_k		= sqrt(1.0/publishLagDateVector->Elt(i,k));
																	
						//Correl INF/INF
						correl_INF_INF		= AvgCPICPICorrel->Elt(l-1,k-1);
							
						//Correl INF/IR
						correl_INF_IR		= AvgCPIIRCorrel->Elt(l-1,k-1);

						var+=maturityFactor_T_l*maturityFactor_T_k*(alpha_l*alpha_k*correl_INF_INF*volCPI_l*volCPI_k)+
							maturityFactor_T_l*(alpha_l*beta_k*volCPI_l*volIR_k*correl_INF_IR)+
							(beta_l*beta_k*volIR_l*volIR_k);

						//fprintf(f," expiry=%7.10lf\t tenor=%7.10lf\t sizeTenor_t=%3d\t sizeTenor_l=%3d\t sizeTenor_k=%3d\t \n", payDate_0, tenor_t, TenorSize_t, l,k);
						//fprintf(f," maturityFactor_T_l=%7.10lf\t maturityFactor_T_k=%7.10f\t alpha_l=%7.10lf\t alpha_k=%7.10lf\t  correl_INF_INF=%7.5lf\t  volCPI_l=%7.5lf\t  volCPI_k=%7.5lf\t beta_k=%7.5lf\t beta_l=%7.5lf\t volIR_l=%7.5lf\t volIR_k=%7.5lf\t correl_INF_IR=%7.5lf\t ATM_Vol_Tenor_t=%7.5lf\t var =%7.5lf\t \n",maturityFactor_T_l,maturityFactor_T_k,alpha_l,alpha_k,correl_INF_INF,volCPI_l,volCPI_k,beta_k,beta_l, volIR_l,volIR_k, correl_INF_IR,ATM_Vol_Tenor_t, var );
					}
				}

				if( var < K_NEW_DOUBLE_TOL )
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": vol is not positive!");
			
								

				var	 = sqrt(var)*CC_NS(ARM_Constants,volBase);
				
				
				SmiledVols[t]->GetVolatilities()->Elt(i,p)=var-ATM_Vol_Tenor_t*CC_NS(ARM_Constants,volBase);
			}		
		}				
	}	
	delete expiries;
	delete tenors;
	delete strikes;
	delete smiledTenors;	

	ARM_VolCube* result = new ARM_VolCube(ATMVol,
										  &SmiledVols[0],
										  nbSmiledVols,
										  smiledTenorsVec,
                                          K_ATMF_VOL);

	delete expiriesVec;
	delete strikesVec;
	delete smiledTenorsVec;
	delete ATMVol;

	vector<ARM_VolCurve*>::iterator it;
        
   for (it = SmiledVols.begin(); it < SmiledVols.end(); it++)
   {
       delete (*it);
   }
	
	return result;
}
*/
/*
for(l=1;l<TenorSize_j;l++)
			{
				double floatleg_l		=	floatleg_j-ProgressiveFloatLeg->Elt(l-1);
				double fixleg_l		=	fixleg_j-ProgressiveFixLeg->Elt(l-1);
				alpha_l			=	DF->Elt(l)*CPIFwdRatio->Elt(l)/floatleg_j;
				beta_l			=	(DF->Elt(l-1)-DF->Elt(l))/DF->Elt(l-1)*(fixleg_l/fixleg_j-floatleg_l/floatleg_j);
				volCPI_l		=	CPIVol_YoY->Elt(l);
				volIR_l			=	ForwardVol->Elt(l-1);

				//For correl Calcul
				ARM_GP_Vector* CPIVol_ZC		= new ARM_GP_Vector(CPIDataSize+1);
				ARM_AutoCleaner< ARM_GP_Vector > HoldCPIVol_ZC( CPIVol_ZC );
				CPIVol_ZC->Elt(0)	= itsInfBSModel->GetVolatility()->ComputeVolatility( timeToStart, lookupStrike, denomVolLookup )/ CC_NS( ARM_Constants, volBase );
				CPIVol_ZC->Elt(l)	= itsInfBSModel->GetVolatility()->ComputeVolatility( timeToStart, lookupStrike, numVolLookup )/ CC_NS( ARM_Constants, volBase );
				double Vol_ZC_l = CPIVol_ZC->Elt(l);
				double Vol_ZC_lmoins1 = CPIVol_ZC->Elt(l-1);
				double T_l			= publishLagDateVector->Elt(l);
				double T_lmoins1	= publishLagDateVector->Elt(l-1);

				double rho_l_lmoins1 = CPICPICorrelMat->ComputeCorrelData(T_l,T_lmoins1)/CC_NS(ARM_Constants,correlBase);
				double var_l = Vol_ZC_l*Vol_ZC_l-2*rho_l_lmoins1*Vol_ZC_l*Vol_ZC_lmoins1+Vol_ZC_lmoins1*Vol_ZC_lmoins1;

				for(k=1;k<TenorSize_j;k++)
				{
					double floatleg_k		=	floatleg_j-ProgressiveFloatLeg->Elt(k-1);
					double fixleg_k		=	fixleg_j-ProgressiveFixLeg->Elt(k-1);
					alpha_k	=	DF->Elt(k)*CPIFwdRatio->Elt(k)/floatleg_j;
					beta_k	=	(DF->Elt(k-1)-DF->Elt(k))/DF->Elt(k-1)*(fixleg_k/fixleg_j-floatleg_k/floatleg_j);
					volCPI_k	= CPIVol_YoY->Elt(k);;
					volIR_k		= ForwardVol->Elt(k-1);
										
					//////////////////////For correl Calcul//////////////////////////
					double Vol_ZC_k			= CPIVol_ZC->Elt(k);
					double Vol_ZC_kmoins1	= CPIVol_ZC->Elt(k-1);

					double T_k			= publishLagDateVector->Elt(k);
					double T_kmoins1	= publishLagDateVector->Elt(k-1);

					double rho_k_kmoins1 = CPICPICorrelMat->ComputeCorrelData(T_k,T_kmoins1)/CC_NS(ARM_Constants,correlBase);
					double var_k = Vol_ZC_k*Vol_ZC_k-2*rho_k_kmoins1*Vol_ZC_k*Vol_ZC_kmoins1+Vol_ZC_kmoins1*Vol_ZC_kmoins1;
	

					//Correl INF/INF
					double rho_l_k, rho_l_kmoins1, rho_lmoins1_k, rho_lmoins1_kmoins1;
					rho_l_k				= CPICPICorrelMat->ComputeCorrelData(T_l,T_k)/CC_NS(ARM_Constants,correlBase);
					rho_l_kmoins1		= CPICPICorrelMat->ComputeCorrelData(T_l,T_kmoins1)/CC_NS(ARM_Constants,correlBase);
					rho_lmoins1_k		= CPICPICorrelMat->ComputeCorrelData(T_lmoins1,T_k)/CC_NS(ARM_Constants,correlBase);
					rho_lmoins1_kmoins1	= CPICPICorrelMat->ComputeCorrelData(T_lmoins1,T_kmoins1)/CC_NS(ARM_Constants,correlBase);

					double covar_l_k = rho_l_k*Vol_ZC_l*Vol_ZC_k-rho_l_kmoins1*Vol_ZC_l*Vol_ZC_kmoins1-
						rho_lmoins1_k*Vol_ZC_lmoins1*Vol_ZC_k+rho_lmoins1_kmoins1*Vol_ZC_lmoins1*Vol_ZC_kmoins1;

					correl_INF_INF	=	covar_l_k/sqrt(var_k*var_l);

					//Correl INF/IR
					double mu_l_kmoins1, mu_lmoins1_kmoins1; 
					mu_l_kmoins1		= CPIIRCorrelMat->ComputeCorrelData(T_l,T_kmoins1)/CC_NS(ARM_Constants,correlBase);
					mu_lmoins1_kmoins1	= CPIIRCorrelMat->ComputeCorrelData(T_lmoins1,T_kmoins1)/CC_NS(ARM_Constants,correlBase);
					double mu_covar_l_k = (mu_l_kmoins1*Vol_ZC_l-mu_lmoins1_kmoins1*Vol_ZC_lmoins1)*volIR_k;

					correl_INF_IR	=	mu_covar_l_k/(volIR_k*sqrt(var_l));

					var+=maturityFactor*(alpha_l*alpha_k*correl_INF_INF*volCPI_l*volCPI_k)+
						sqrtMaturityFactor*(alpha_l*beta_k*volCPI_l*volIR_k*correl_INF_IR)+
						(beta_l*beta_k*volIR_l*volIR_k);
				}
			}
*/
////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_Std
///	Routine: GenerateSwopOATVolCurve
///	Returns: ARM_VolCurve* 
///	Action : general function to compute the volatility of 
///			inf swaption from the one of zero coupon
////////////////////////////////////////////////////

ARM_VolCurve* ARM_InfVolComputation_Producer_Std::GenerateOATSwopVolCurve(
		const ARM_Date&		asOfDate,
		ARM_GP_Vector*		tenors,
		ARM_GP_Vector*		expiries,
		double	coupon, 
		bool	mode)
{
	/// get the default data if missing
	GetDefaultDataIfMissing( tenors, expiries );
	
	ARM_Matrix* volatilities = new ARM_Matrix( expiries->size(), tenors->size());
	ARM_AutoCleaner< ARM_Matrix > HoldVolatilities( volatilities );

	/// get default of the currency of the IR market
	double fixIRFreq	=	itsInfBSModel->GetZeroCurve()->GetCurrencyUnit()->GetFixedPayFreq();	
	double fixIRPeriod	= 1.0/fixIRFreq;
	int fixIRDayCount	= itsInfBSModel->GetZeroCurve()->GetCurrencyUnit()->GetFixedDayCount();
	double IRFactor		= fixIRDayCount == KACTUAL_360? 365.0/360.0 : 1;
	double fixIRperiodDC= fixIRPeriod*IRFactor;

	/// it is assumed that the convention on the CPI frequency is K_ANNUAL
	double CPIFreq		= K_ANNUAL;
	double CPIPeriod	= 1.0/CPIFreq;
	double asOf			= itsInfBSModel->GetZeroCurve()->GetAsOfDate().GetJulian();
	int infDayCount		= itsInfBSModel->GetInfFwdCurv()->GetMonthlyInterpType();
	ARM_Date numDate, denomDate, lastKnownDate, refDate;
	double numYearFrac, denomYearFrac, denomYearFrac2;
	string CPIIndexName = itsInfBSModel->GetInfFwdCurv()->GetInfIdxName();

	/// recover the inflation index property
	ARM_InfIdx* CPIIndex= itsInfBSModel->GetInfFwdCurv()->GetInfIdx();
	ARM_AutoCleaner< ARM_InfIdx > HoldCPIIndex( CPIIndex );

	long dailyInterpType= InfData::GetDailyInterpolation( CPIIndexName.c_str());
	double lookupStrike = 0., ATMStrikeDefault = 0., refDateFrac = 0.;


	/// data for the computation of the volatilities
	ARM_GP_Vector* pFixLegDF			= NULL; 
	ARM_GP_Vector* pFixLegIRYearFrac	= NULL;  
	ARM_GP_Vector* pFixLegDFVol			= NULL; 
	ARM_GP_Vector* pFwdCPI				= NULL; 
	ARM_GP_Vector* pVol_ZC				= NULL; 
	ARM_GP_Vector* pInfLegDF			= NULL; 
	ARM_GP_Vector* pInfLegDFVol			= NULL; 
	ARM_GP_Vector* pAvgCPIIRCorrel		= NULL; 

	ARM_CorrelMatrix* CPIIRCorrelMat	= itsInfBSModel->GetInfIRCorrelMatrix(CPIIndex,CPIIndex,  "INF/IR" );
	ARM_CorrelMatrix* CPICPICorrelMat	= itsInfBSModel->GetInfInfCorrelMatrix(CPIIndex,CPIIndex, "INF/INF");

	size_t i, j, k, l;
	/// 1) computes the volatility using the relationship
	for( i=0; i<expiries->size(); ++i )
	{
		/// computes all the IR data for the longest tenors
		double longestTenor = tenors->Elt(tenors->size()-1);
		size_t IRDataSize	= longestTenor*fixIRFreq;
		
		ARM_GP_Vector* fixLegDF			= new ARM_GP_Vector(IRDataSize);
		ARM_AutoCleaner< ARM_GP_Vector > HoldfixLegDF( fixLegDF );

		ARM_GP_Vector* fixLegDFVol			= new ARM_GP_Vector(IRDataSize);
		ARM_AutoCleaner< ARM_GP_Vector > HoldfixLegDFVol( fixLegDFVol );

		ARM_GP_Vector* fixLegIRYearFrac = new ARM_GP_Vector(IRDataSize, fixIRperiodDC);
		ARM_AutoCleaner< ARM_GP_Vector > HoldfixLegIRYearFrac( fixLegIRYearFrac );

		for( j=0; j<IRDataSize; ++j )
		{
			fixLegDF->Elt(j)	= itsInfBSModel->GetDiscountCurve()->DiscountPrice(expiries->Elt(i)+(j+1)*fixIRPeriod);
			fixLegDFVol->Elt(j)	= itsInfBSModel->GetVolatility()->ComputeVolatility(expiries->Elt(i), ATMStrikeDefault, (j+1)*fixIRPeriod )/CC_NS( ARM_Constants, volBase ) * ( 1.0-fixLegDF->Elt(j));
		}

		/// computes all the INF data for the longest tenors
		size_t InfDataSize = longestTenor*CPIFreq;

		ARM_GP_Vector* FwdCPI			= new ARM_GP_Vector(InfDataSize);
		ARM_AutoCleaner< ARM_GP_Vector > HoldFwdCPI( FwdCPI );

		ARM_GP_Vector* CPIFwdVol_ZC		= new ARM_GP_Vector(InfDataSize);
		ARM_AutoCleaner< ARM_GP_Vector > HoldCPIFwdVol_ZC( CPIFwdVol_ZC );

		ARM_GP_Vector* InfLegDF			= new ARM_GP_Vector(InfDataSize);
		ARM_AutoCleaner< ARM_GP_Vector > HoldInfLegDF( InfLegDF );

		ARM_GP_Vector* InfLegDFVol		= new ARM_GP_Vector(InfDataSize);
		ARM_AutoCleaner< ARM_GP_Vector > HoldInfLegDFVol( InfLegDFVol );

		ARM_GP_Vector* AvgCPIIRCorrel	= new ARM_GP_Vector(InfDataSize);
		ARM_AutoCleaner< ARM_GP_Vector > HoldAvgCPIIRCorrel( AvgCPIIRCorrel );

		ARM_Matrix* AvgCPICPICorrel		= new ARM_Matrix(InfDataSize,InfDataSize,1.0);
		ARM_AutoCleaner< ARM_Matrix > HoldAvgCPICPICorrel( AvgCPICPICorrel );

		double timeToStart = asOfDate.GetJulian(), TVolLookup;
		lastKnownDate	= itsInfBSModel->GetVolatility()->GetLastKnownDate();

		for( j=0; j<InfDataSize; ++j )
		{
			InfLegDF->Elt(j)		= itsInfBSModel->GetDiscountCurve()->DiscountPrice(expiries->Elt(i)+(j+1)*CPIFreq);
			InfLegDFVol->Elt(j)		= itsInfBSModel->GetVolatility()->ComputeVolatility( expiries->Elt(i), ATMStrikeDefault, (j+1)*CPIFreq ) / CC_NS( ARM_Constants, volBase ) * ( 1.-InfLegDF->Elt(j) );
			numYearFrac				= expiries->Elt(i)+(j+1)*CPIFreq;
			denomYearFrac			= expiries->Elt(i)+j*CPIFreq;
			numDate					= ARM_Date(asOf +  numYearFrac * K_YEAR_LEN);
			//denomDate				= ARM_Date(asOf +  denomYearFrac * K_YEAR_LEN);
			refDateFrac				= expiries->Elt(i);
			refDate					= ARM_Date(asOf +  denomYearFrac * K_YEAR_LEN);
			ARM_Date numDatePublishLag = itsInfBSModel->GetModelDateWPublishLag( numDate, CPIIndex );

			/// à revoir (s)
			FwdCPI->Elt(j)			= itsInfBSModel->FwdCPIRatio( numDate, refDate, numDatePublishLag, 1., 0., dailyInterpType, GETDEFAULTVALUE, CPIIndex )[0];
			/// compute the fwd ZC vols as a sum of quadratic average of the year on year vols
			for( k=0; k<=j; ++k )
			{
				denomYearFrac2		= expiries->Elt(i)+k*CPIFreq;
				TVolLookup			= itsInfBSModel->GetVolatility()->ComputeVolatility( denomYearFrac2, lookupStrike, CPIFreq)/ CC_NS( ARM_Constants, volBase );
				CPIFwdVol_ZC->Elt(j)	+= denomYearFrac2*TVolLookup*TVolLookup;
			}
			CPIFwdVol_ZC->Elt(j) = sqrt(CPIFwdVol_ZC->Elt(j)/denomYearFrac);

			AvgCPIIRCorrel->Elt(j)	= CPIIRCorrelMat->ComputeCorrelData(numYearFrac,numYearFrac)/CC_NS(ARM_Constants,correlBase);
			/// fill the upper triangle of the matrix
			for( k=0; k<j; ++k )
			{
				denomYearFrac2		= expiries->Elt(i)+k*CPIFreq;
				AvgCPICPICorrel->Elt(k,j) = CPICPICorrelMat->ComputeCorrelData(denomYearFrac2,denomYearFrac)/CC_NS(ARM_Constants,correlBase);
			}
		}

		/// start calling the function to compute volatility with correct data
		for (j=0;j<tenors->size();j++)
		{
			/// need to correct bond vol to be really bond vol and not forward volatility
			/// IR data
			size_t fixLegDataSize = tenors->Elt(j) * fixIRFreq;
			CopyArgsPointorInPlaceNoDelete< ARM_GP_Vector, ARM_GP_Vector >( pFixLegDF, *fixLegDF,fixLegDataSize  );
			ARM_AutoCleaner< ARM_GP_Vector > HoldpFixLegDF( pFixLegDF );
			CopyArgsPointorInPlaceNoDelete< ARM_GP_Vector, ARM_GP_Vector >( pFixLegIRYearFrac, *fixLegIRYearFrac, fixLegDataSize  );
			ARM_AutoCleaner< ARM_GP_Vector > HoldpFixLegIRYearFrac( pFixLegIRYearFrac);
			CopyArgsPointorInPlaceNoDelete< ARM_GP_Vector, ARM_GP_Vector >( pFixLegDFVol, *fixLegDFVol, fixLegDataSize  );
			ARM_AutoCleaner< ARM_GP_Vector > HoldpFixLegDFVol( pFixLegDFVol );

			/// Inf data
			size_t infLegDataSize = tenors->Elt(j) * CPIFreq;
			CopyArgsPointorInPlaceNoDelete< ARM_GP_Vector, ARM_GP_Vector >( pFwdCPI			, *FwdCPI,			infLegDataSize );	
			ARM_AutoCleaner< ARM_GP_Vector > HoldpFwdCPIRatio( pFwdCPI );
			CopyArgsPointorInPlaceNoDelete< ARM_GP_Vector, ARM_GP_Vector >( pVol_ZC			, *CPIFwdVol_ZC,		infLegDataSize );	
			ARM_AutoCleaner< ARM_GP_Vector > HoldpVol_YoY( pVol_ZC );
			CopyArgsPointorInPlaceNoDelete< ARM_GP_Vector, ARM_GP_Vector >( pInfLegDF		, *InfLegDF,		infLegDataSize );	
			ARM_AutoCleaner< ARM_GP_Vector > HoldpInfLegDF( pInfLegDF );
			CopyArgsPointorInPlaceNoDelete< ARM_GP_Vector, ARM_GP_Vector >( pInfLegDFVol	, *InfLegDFVol,		infLegDataSize );	
			ARM_AutoCleaner< ARM_GP_Vector > HoldpInfLegDFVol( pInfLegDFVol );
			CopyArgsPointorInPlaceNoDelete< ARM_GP_Vector, ARM_GP_Vector >( pAvgCPIIRCorrel	, *AvgCPIIRCorrel,	infLegDataSize );	
			ARM_AutoCleaner< ARM_GP_Vector > HoldAvgCPIIRCorrel( pAvgCPIIRCorrel );

			/// Fill the different correlation matrixes
			/// CPI vs. DF float
			ARM_Matrix* AvgFloatCPIIRCorrel		= new ARM_Matrix(pInfLegDF->size(),pInfLegDF->size(), 1. ); 
			ARM_AutoCleaner< ARM_Matrix > HoldAvgFloatCPIIRCorrel( AvgFloatCPIIRCorrel );
			
			for( k=0; k<pInfLegDF->size(); ++k )
				for( l=0; l<pInfLegDF->size(); ++l )
						AvgFloatCPIIRCorrel->Elt(k,l)= AvgCPIIRCorrel->Elt(k);

			/// CPI vs. DF fix
			ARM_Matrix* AvgFloatFixCPIIRCorrel	= new ARM_Matrix(pInfLegDF->size(),pFixLegDF->size(), 1. );
			ARM_AutoCleaner< ARM_Matrix > HoldAvgFloatFixCPIIRCorrel( AvgFloatFixCPIIRCorrel );
			
			for( k=0; k<pInfLegDF->size(); ++k )
				for( l=0; l<pFixLegDF->size(); ++l )
						AvgFloatFixCPIIRCorrel->Elt(k,l)= AvgCPIIRCorrel->Elt(k);

			/// CPI vs. CPI
			ARM_Matrix* AvgFloatCPICPICorrel	= new ARM_Matrix(pInfLegDF->size(),pInfLegDF->size(), 1.0 ); 
			ARM_AutoCleaner< ARM_Matrix > HoldAvgFloatCPICPICorrel( AvgFloatCPICPICorrel );
			
			for( k=0; k<pInfLegDF->size(); ++k )
			{
				for( l=0; l<k; ++l )
					AvgFloatCPICPICorrel->Elt(l,k) = AvgCPICPICorrel->Elt(l,k);
				for( l=k+1; l<pInfLegDF->size(); ++l )
					AvgFloatCPICPICorrel->Elt(l,k) = AvgCPICPICorrel->Elt(k,l);
			}

			/// DF vs. DF: Parfect correlation !
			ARM_Matrix* AvgFixIRIRCorrel		= new ARM_Matrix(pFixLegDF->size(),pFixLegDF->size(), 1.0 ); 
			ARM_AutoCleaner< ARM_Matrix > HoldAvgFixIRIRCorrel( AvgFixIRIRCorrel );

			ARM_Matrix* AvgFloatIRIRCorrel		= new ARM_Matrix(pInfLegDF->size(),pInfLegDF->size(), 1.0 ); 
			ARM_AutoCleaner< ARM_Matrix > HoldAvgFloatIRIRCorrel( AvgFloatIRIRCorrel );
			
			ARM_Matrix* AvgFloatFixIRIRCorrel	= new ARM_Matrix(pInfLegDF->size(),pFixLegDF->size(), 1.0 ); 
			ARM_AutoCleaner< ARM_Matrix > HoldAvgFloatFixIRIRCorrel( AvgFloatFixIRIRCorrel );


			/// Call the fucntion VolZC_to_VolOATSwp to compute the needed vol
			volatilities->Elt(i,j) = VolZC_to_VolOATSwp(
				pFixLegDF, 
				pFixLegIRYearFrac, 
				pFixLegDFVol,
				pFwdCPI,
				pVol_ZC,
				pInfLegDF,
				pInfLegDFVol,
				AvgFloatCPIIRCorrel,
				AvgFloatFixCPIIRCorrel,
				AvgFloatCPICPICorrel,
				AvgFixIRIRCorrel,
				AvgFloatIRIRCorrel,
				AvgFloatFixIRIRCorrel,
				coupon,
				mode)*CC_NS(ARM_Constants,volBase);

		}

	}

	/// 2) set the volatility value to the matrix
	ARM_Vector* expiriesVec	= To_pARM_Vector( expiries );
	delete expiries;
	ARM_Vector* tenorsVec	= To_pARM_Vector( tenors );
	delete tenors;

	ARM_VolCurve* volCurve = new ARM_VolLInterpol(asOfDate, expiriesVec, tenorsVec, volatilities, K_STK_TYPE_PRICE, K_ATMF_VOL );
	volCurve->SetVolatilities(volatilities);
	HoldVolatilities.Release();
	

	return (volCurve);

}

////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_Std
///	Routine: GenerateSwopOATVolCurve
///	Returns: ARM_VolCurve* 
///	Action : general function to compute the volatility of 
///			inf swaption from the one of zero coupon
////////////////////////////////////////////////////

ARM_VolCurve* ARM_InfVolComputation_Producer_Std::CompleteOATSwopVolCurve(
		const ARM_Date&	asOfDate,
		ARM_VolCurve*	vol,
		double			coupon)
{
	ARM_GP_Vector* tenors = new ARM_GP_Vector(ARM_InfVolComputation_Producer::DefaultExpiriesSize);
//	ARM_AutoCleaner< ARM_GP_Vector > Holdtenors( tenors );
	ARM_GP_Vector* expiries = new ARM_GP_Vector(ARM_InfVolComputation_Producer::DefaultExpiriesSize);
//	ARM_AutoCleaner< ARM_GP_Vector > Holdexpiries( expiries );

	ARM_InfVolComputation_Producer::GetDefaultDataIfMissing(
	tenors,
	expiries );

	ARM_VolCurve* newVol = ARM_InfVolComputation_Producer_Std::GenerateOATSwopVolCurve(
		asOfDate,
		tenors,
		expiries,
		coupon, 
		true);

	return (newVol);
}
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
///	Class      : ARM_InfVolComputation_Producer_EqualWeight
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_EqualWeight
///	Routine: Constructor
////////////////////////////////////////////////////

ARM_InfVolComputation_Producer_EqualWeight::ARM_InfVolComputation_Producer_EqualWeight(ARM_InfBSModel* infIRBSModel)
:	ARM_InfVolComputation_Producer(infIRBSModel)
{
	SetName(ARM_INFSWWOPTVOL_EQUALWEIGHT_PRODUCER);
}


////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_EqualWeight
///	Routine: Copy Constructor
////////////////////////////////////////////////////
ARM_InfVolComputation_Producer_EqualWeight::ARM_InfVolComputation_Producer_EqualWeight(
	const ARM_InfVolComputation_Producer_EqualWeight& rhs )
:	ARM_InfVolComputation_Producer(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_EqualWeight
///	Routine: Assignement operator
////////////////////////////////////////////////////
ARM_InfVolComputation_Producer_EqualWeight& ARM_InfVolComputation_Producer_EqualWeight::operator=( 
	const ARM_InfVolComputation_Producer_EqualWeight& rhs )
{
	if( this != &rhs )
		ARM_InfVolComputation_Producer::operator =( rhs );
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_EqualWeight
///	Routine: Destructor
////////////////////////////////////////////////////
ARM_InfVolComputation_Producer_EqualWeight::~ARM_InfVolComputation_Producer_EqualWeight()
{}

////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_EqualWeight
///	Routine: Clone
////////////////////////////////////////////////////
ARM_Object* ARM_InfVolComputation_Producer_EqualWeight::Clone()
{
	return new ARM_InfVolComputation_Producer_EqualWeight(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_EqualWeight
///	Routine: Clone
////////////////////////////////////////////////////
string ARM_InfVolComputation_Producer_EqualWeight::toString( const string& indent )
{
	return "ARM_InfVolComputation_Producer_EqualWeight";
}



////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_EqualWeight
///	Routine: GenerateSwopVolCurve
///	Returns: double
///	Action : general function to compute the volatility of inf swaption from the one 
///				of year on year
////////////////////////////////////////////////////

ARM_VolCurve* ARM_InfVolComputation_Producer_EqualWeight::GenerateSwopVolCurve( 
	const ARM_Date&		asOfDate,
	ARM_GP_Vector*			tenors,
	ARM_GP_Vector*			expiries )
{
	GetDefaultDataIfMissing( tenors, expiries );
	ARM_Matrix* volatilities= new ARM_Matrix( expiries->size(), tenors->size() );
	ARM_AutoCleaner< ARM_Matrix > HoldVolatilities( volatilities );

	/// it is assumed that the convention on the CPI frequency is K_ANNUAL
	double CPIFreq		= K_ANNUAL;
	double CPIPeriod	= 1.0/CPIFreq;

	/// data for the computation of the volatilities
	ARM_GP_Vector* pVol_YoY			= NULL; 
	double lookupStrike				= 0;
	double denomYearFrac;

	/// 1) computes the volatility using the relationship
	for( size_t i=0; i<expiries->size(); ++i )
	{
		/// computes all the IR data for the longest tenors
		double longestTenor			= tenors->Elt(tenors->size()-1);

		/// computes all the data for the CPIs
		size_t CPIDataSize			= longestTenor*CPIFreq;

		ARM_GP_Vector* CPIVol_YoY	= new ARM_GP_Vector(CPIDataSize);
		ARM_AutoCleaner< ARM_GP_Vector > HoldCPIVol_YoY( CPIVol_YoY );

		for( size_t j=0; j<CPIDataSize; ++j )
		{
			denomYearFrac			= expiries->Elt(i)+j*CPIFreq;
			CPIVol_YoY->Elt(j)		= itsInfBSModel->GetVolatility()->ComputeVolatility( denomYearFrac, lookupStrike, CPIPeriod )/ CC_NS( ARM_Constants, volBase );
		}

		double sum = 0;
		/// start calling the function to compute volatility with correct data
		for(j=0; j<tenors->size(); ++j )
		{
			sum += CPIVol_YoY->Elt(j);
			volatilities->Elt(i,j) = sum/(j+1) *CC_NS(ARM_Constants,volBase);;
		}
	}

	/// Set the data to the volatilities matrix
	ARM_Vector* expiriesVec	= To_pARM_Vector( expiries );
	delete expiries;
	ARM_Vector* tenorsVec	= To_pARM_Vector( tenors );
	delete tenors;

	ARM_VolCurve* volCurve = new ARM_VolLInterpol(asOfDate, expiriesVec, tenorsVec, volatilities, K_STK_TYPE_PRICE, K_ATMF_VOL );
	volCurve->SetVolatilities(volatilities);
	HoldVolatilities.Release();
	return volCurve;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_EqualWeight
///	Routine: GenerateSwopVolCube
///	Returns: double
///	Action : general function to compute the volatility of inf swaption from the one 
///				of year on year
////////////////////////////////////////////////////
ARM_VolCube* ARM_InfVolComputation_Producer_EqualWeight::GenerateSwopVolCube( 
	const ARM_Date&			asOfDate,
	ARM_GP_Vector*			tenors,
	ARM_GP_Vector*			expiries,
	ARM_GP_Vector*			smiledTenors,
	ARM_GP_Vector*			strikes)
{
	ARM_Date LastKnownDate;
	ARM_Date AsOf;
	vector< ARM_VolCurve* >* Vols =NULL;
	ARM_Vector* tmpUnderlyings = To_pARM_Vector(strikes);
	ARM_VolCube* result = new ARM_VolCube( Vols, tmpUnderlyings, LastKnownDate );
	delete tmpUnderlyings;
	return result;
}


///////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	Class      : ARM_InfVolComputation_Producer_SimpleWeight
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////



////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_SimpleWeight
///	Routine: Constructor
////////////////////////////////////////////////////

ARM_InfVolComputation_Producer_SimpleWeight::ARM_InfVolComputation_Producer_SimpleWeight(
	ARM_InfBSModel* infIRBSModel, bool usedSquareVersion )
:	ARM_InfVolComputation_Producer(infIRBSModel), itsUsedSquareVersion(usedSquareVersion)
{
	SetName(ARM_INFSWWOPTVOL_SIMPLEWEIGHT_PRODUCER);
}


////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_SimpleWeight
///	Routine: Copy Constructor
////////////////////////////////////////////////////
ARM_InfVolComputation_Producer_SimpleWeight::ARM_InfVolComputation_Producer_SimpleWeight(
	const ARM_InfVolComputation_Producer_SimpleWeight& rhs )
:	ARM_InfVolComputation_Producer(rhs),
	itsUsedSquareVersion(rhs.itsUsedSquareVersion)
{}


////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_SimpleWeight
///	Routine: Assignement operator
////////////////////////////////////////////////////
ARM_InfVolComputation_Producer_SimpleWeight& ARM_InfVolComputation_Producer_SimpleWeight::operator=( 
	const ARM_InfVolComputation_Producer_SimpleWeight& rhs )
{
	if( this != &rhs )
	{
		ARM_InfVolComputation_Producer::operator =( rhs );
		itsUsedSquareVersion = rhs.itsUsedSquareVersion;
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_SimpleWeight
///	Routine: Destructor
////////////////////////////////////////////////////
ARM_InfVolComputation_Producer_SimpleWeight::~ARM_InfVolComputation_Producer_SimpleWeight()
{}

////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_SimpleWeight
///	Routine: Clone
////////////////////////////////////////////////////
ARM_Object* ARM_InfVolComputation_Producer_SimpleWeight::Clone()
{
	return new ARM_InfVolComputation_Producer_SimpleWeight(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_SimpleWeight
///	Routine: Clone
////////////////////////////////////////////////////
string ARM_InfVolComputation_Producer_SimpleWeight::toString( const string& indent )
{
	return "ARM_InfVolComputation_Producer_SimpleWeight";
}


////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_SimpleWeight
///	Routine: GenerateSwopVolCurve
///	Returns: double
///	Action : general function to compute the volatility of inf swaption from the one 
///				of year on year
////////////////////////////////////////////////////
ARM_VolCurve* ARM_InfVolComputation_Producer_SimpleWeight::GenerateSwopVolCurve( 
	const ARM_Date&	asOfDate,
	ARM_GP_Vector*		tenors,
	ARM_GP_Vector*		expiries )
{
	/// get the default data if missing
	GetDefaultDataIfMissing( tenors, expiries );
	ARM_Matrix* volatilities= new ARM_Matrix( expiries->size(), tenors->size() );
	ARM_AutoCleaner< ARM_Matrix > HoldVolatilities( volatilities );

	/// it is assumed that the convention on the CPI frequency is K_ANNUAL
	double CPIFreq		= K_ANNUAL;
	double CPIPeriod	= 1.0/CPIFreq;
	double lookupStrike	= 0.0;
	double denomYearFrac;

	/// 1) computes the volatility using the relationship
	for( size_t i=0; i<expiries->size(); ++i )
	{
		/// computes all the data for the CPIs
		double longestTenor			= tenors->Elt(tenors->size()-1);
		size_t CPIDataSize			= longestTenor*CPIFreq;

		ARM_GP_Vector* CPIVol_YoY	= new ARM_GP_Vector(CPIDataSize);
		ARM_AutoCleaner< ARM_GP_Vector > HoldCPIVol_YoY( CPIVol_YoY );

		ARM_GP_Vector* CPILegDF		= new ARM_GP_Vector(CPIDataSize);
		ARM_AutoCleaner< ARM_GP_Vector > HoldCPILegDF( CPILegDF );

		for( size_t j=0; j<CPIDataSize; ++j )
		{
			CPILegDF->Elt(j)		= itsInfBSModel->GetDiscountCurve()->DiscountPrice(expiries->Elt(i)+(j+1)*CPIFreq);
			denomYearFrac			= expiries->Elt(i)+j*CPIFreq;
			CPIVol_YoY->Elt(j)		= itsInfBSModel->GetVolatility()->ComputeVolatility( denomYearFrac, lookupStrike, CPIPeriod )/ CC_NS( ARM_Constants, volBase );
		}

		/// start calling the function to compute volatility with correct data
		if( itsUsedSquareVersion )
		{
			double sumSquareDFSigma = 0.0, sumDF = 0.0;
			ARM_InfIdx* CPIIndex	= itsInfBSModel->GetInfFwdCurv()->GetInfIdx();
			ARM_CorrelMatrix* correlMat	= itsInfBSModel->GetInfInfCorrelMatrix(CPIIndex,CPIIndex, "INF/INF" );

			/// to compute recursively the square, we add at each step the two symmetric term plus the square term!
			for(j=0; j<tenors->size(); ++j )
			{
				for( size_t k=0; k<j; ++k )
				{
					sumSquareDFSigma += 2.0* CPILegDF->Elt(j)*CPILegDF->Elt(k)*CPIVol_YoY->Elt(j)*CPIVol_YoY->Elt(k)
						* correlMat->ComputeCorrelData(expiries->Elt(i)+j*CPIFreq, expiries->Elt(i)+k*CPIFreq)/ CC_NS(ARM_Constants,correlBase);
				}
				sumSquareDFSigma += CPILegDF->Elt(j)*CPILegDF->Elt(j)*CPIVol_YoY->Elt(j)*CPIVol_YoY->Elt(j);
				sumDF	   += CPILegDF->Elt(j);
				volatilities->Elt(i,j) = sqrt(sumSquareDFSigma)/sumDF * CC_NS(ARM_Constants,volBase);
			}
		}		
		/// simple computation as the weighted average of year on year volatilities with df over annuity!
		else
		{
			double sumDFSigma = 0.0, sumDF = 0.0;

			for(j=0; j<tenors->size(); ++j )
			{
				sumDFSigma += CPILegDF->Elt(j) * CPIVol_YoY->Elt(j);
				sumDF	   += CPILegDF->Elt(j);
				volatilities->Elt(i,j) = sumDFSigma/sumDF * CC_NS(ARM_Constants,volBase);;
			}
		}
	}

	/// Set the data to the volatilities matrix
	ARM_Vector* expiriesVec	= To_pARM_Vector( expiries );
	delete expiries;
	ARM_Vector* tenorsVec	= To_pARM_Vector( tenors );
	delete tenors;

	ARM_VolCurve* volCurve = new ARM_VolLInterpol(asOfDate, expiriesVec, tenorsVec, volatilities, K_STK_TYPE_PRICE, K_ATMF_VOL );
	volCurve->SetVolatilities(volatilities);

	HoldVolatilities.Release();
	return volCurve;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_SimpleWeight
///	Routine: GenerateSwopVolCube
///	Returns: double
///	Action : general function to compute the volatility of inf swaption from the one 
///				of year on year
////////////////////////////////////////////////////

ARM_VolCube* ARM_InfVolComputation_Producer_SimpleWeight::GenerateSwopVolCube( 
	const ARM_Date&			asOfDate,
	ARM_GP_Vector*			tenors,
	ARM_GP_Vector*			expiries,
	ARM_GP_Vector*			smiledTenors,
	ARM_GP_Vector*			strikes)
{
	ARM_Date LastKnownDate;
	ARM_Date AsOf;
	vector< ARM_VolCurve* >* Vols =NULL;
	ARM_Vector* tmpUnderlyings = To_pARM_Vector(strikes);
	ARM_VolCube* result = new ARM_VolCube( Vols, tmpUnderlyings, LastKnownDate );
	delete tmpUnderlyings;
	return result;
}

///////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	Class      : ARM_InfOATVolComputation_Producer_SimpleWeight
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////



////////////////////////////////////////////////////
///	Class  : ARM_InfOATVolComputation_Producer_SimpleWeight
///	Routine: Constructor
////////////////////////////////////////////////////

ARM_InfOATVolComputation_Producer_SimpleWeight::ARM_InfOATVolComputation_Producer_SimpleWeight(
	ARM_InfBSModel* infIRBSModel, bool usedSquareVersion )
:	ARM_InfVolComputation_Producer(infIRBSModel), itsUsedSquareVersion(usedSquareVersion)
{
	SetName(ARM_INFSWWOPTVOL_SIMPLEWEIGHT_PRODUCER);
}


////////////////////////////////////////////////////
///	Class  : ARM_InfOATVolComputation_Producer_SimpleWeight
///	Routine: Copy Constructor
////////////////////////////////////////////////////
ARM_InfOATVolComputation_Producer_SimpleWeight::ARM_InfOATVolComputation_Producer_SimpleWeight(
	const ARM_InfOATVolComputation_Producer_SimpleWeight& rhs )
:	ARM_InfVolComputation_Producer(rhs),
	itsUsedSquareVersion(rhs.itsUsedSquareVersion)
{}


////////////////////////////////////////////////////
///	Class  : ARM_InfOATVolComputation_Producer_SimpleWeight
///	Routine: Assignement operator
////////////////////////////////////////////////////
ARM_InfOATVolComputation_Producer_SimpleWeight& ARM_InfOATVolComputation_Producer_SimpleWeight::operator=( 
	const ARM_InfOATVolComputation_Producer_SimpleWeight& rhs )
{
	if( this != &rhs )
	{
		ARM_InfVolComputation_Producer::operator =( rhs );
		itsUsedSquareVersion = rhs.itsUsedSquareVersion;
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfOATVolComputation_Producer_SimpleWeight
///	Routine: Destructor
////////////////////////////////////////////////////
ARM_InfOATVolComputation_Producer_SimpleWeight::~ARM_InfOATVolComputation_Producer_SimpleWeight()
{}

////////////////////////////////////////////////////
///	Class  : ARM_InfOATVolComputation_Producer_SimpleWeight
///	Routine: Clone
////////////////////////////////////////////////////
ARM_Object* ARM_InfOATVolComputation_Producer_SimpleWeight::Clone()
{
	return new ARM_InfOATVolComputation_Producer_SimpleWeight(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_InfOATVolComputation_Producer_SimpleWeight
///	Routine: Clone
////////////////////////////////////////////////////
string ARM_InfOATVolComputation_Producer_SimpleWeight::toString( const string& indent )
{
	return "ARM_InfOATVolComputation_Producer_SimpleWeight";
}


////////////////////////////////////////////////////
///	Class  : ARM_InfOATVolComputation_Producer_SimpleWeight
///	Routine: GenerateOATSwopVolCurve
///	Returns: double
///	Action : general function to compute the volatility of inf swaption from the one 
///				of year on year
////////////////////////////////////////////////////
ARM_VolCurve* ARM_InfOATVolComputation_Producer_SimpleWeight::GenerateOATSwopVolCurve( 
	const ARM_Date&	asOfDate,
	ARM_GP_Vector*	tenors,
	ARM_GP_Vector*	expiries )
{
	/// get the default data if missing
	GetDefaultDataIfMissing( tenors, expiries );
	ARM_Matrix* volatilities= new ARM_Matrix( expiries->size(), tenors->size() );
	ARM_AutoCleaner< ARM_Matrix > HoldVolatilities( volatilities );

	/// it is assumed that the convention on the CPI frequency is K_ANNUAL
	double	CPIFreq			= K_ANNUAL;
	double	CPIPeriod		= 1.0/CPIFreq;
	double	lookupStrike	= 0.;
	double	denomYearFrac, denomYearFrac2;
	double	longestTenor	= tenors->Elt(tenors->size()-1);
	double	timeToStart		= asOfDate.GetJulian(), TVolLookup;
	int		infDayCount		= itsInfBSModel->GetInfFwdCurv()->GetMonthlyInterpType();
	ARM_Date lastKnownDate	= itsInfBSModel->GetVolatility()->GetLastKnownDate(), denomDate;
	//ARM_Date& asOf = asOfDate;
	double asOf = itsInfBSModel->GetZeroCurve()->GetAsOfDate().GetJulian();


	/// 1) computes the volatility using the relationship
	for( size_t i=0; i<expiries->size(); ++i )
	{
		/// computes all the data for the CPIs
		
		size_t CPIDataSize			= longestTenor*CPIFreq;

		ARM_GP_Vector* CPIFwdVol_ZC	= new ARM_GP_Vector(CPIDataSize);
		ARM_AutoCleaner< ARM_GP_Vector > HoldCPIFwdVol_ZC( CPIFwdVol_ZC );

		ARM_GP_Vector* CPILegDF		= new ARM_GP_Vector(CPIDataSize);
		ARM_AutoCleaner< ARM_GP_Vector > HoldCPILegDF( CPILegDF );

		for( size_t j=0; j<CPIDataSize; ++j )
		{
			CPILegDF->Elt(j)		= itsInfBSModel->GetDiscountCurve()->DiscountPrice(expiries->Elt(i)+(j+1)*CPIFreq);
			denomYearFrac			= expiries->Elt(i)+j*CPIFreq;
			//denomDate				= ARM_Date(asOf + denomYearFrac * K_YEAR_LEN);
			//TVolLookup				= CountYearsWithoutException( infDayCount, lastKnownDate, denomDate );
			//CPIFwdVol_ZC->Elt(j)		= itsInfBSModel->GetVolatility()->ComputeVolatility( asOf, lookupStrike, TVolLookup )/ CC_NS( ARM_Constants, volBase );

			/// compute the fwd ZC vols as a sum of quadratic average of the year on year vols
			for( size_t k=0; k<=j; ++k )
			{
				denomYearFrac2		= expiries->Elt(i)+k*CPIFreq;
				TVolLookup			= itsInfBSModel->GetVolatility()->ComputeVolatility( denomYearFrac2, lookupStrike, CPIFreq)/ CC_NS( ARM_Constants, volBase );
				CPIFwdVol_ZC->Elt(j)	+= denomYearFrac2*TVolLookup*TVolLookup;
			}
			CPIFwdVol_ZC->Elt(j) = sqrt(CPIFwdVol_ZC->Elt(j)/denomYearFrac);
			
		}

		/// start calling the function to compute volatility with correct data
		if( itsUsedSquareVersion )
		{
			double sumSquareDFSigma = 0.0, sumDF = 0.0;
			ARM_InfIdx* CPIIndex	= itsInfBSModel->GetInfFwdCurv()->GetInfIdx();
			ARM_CorrelMatrix* correlMat	= itsInfBSModel->GetInfInfCorrelMatrix(CPIIndex,CPIIndex, "INF/INF" );

			/// to compute recursively the square, we add at each step the two symmetric term plus the square term!
			for(j=0; j<tenors->size(); ++j )
			{
				for( size_t k=0; k<j; ++k )
				{
					sumSquareDFSigma += 2.0* CPILegDF->Elt(j)*CPILegDF->Elt(k)*CPIFwdVol_ZC->Elt(j)*CPIFwdVol_ZC->Elt(k)
						* correlMat->ComputeCorrelData(expiries->Elt(i)+j*CPIFreq, expiries->Elt(i)+k*CPIFreq)/ CC_NS(ARM_Constants,correlBase);
				}
				sumSquareDFSigma += CPILegDF->Elt(j)*CPILegDF->Elt(j)*CPIFwdVol_ZC->Elt(j)*CPIFwdVol_ZC->Elt(j);
				sumDF	   += CPILegDF->Elt(j);
				volatilities->Elt(i,j) = sqrt(sumSquareDFSigma)/sumDF * CC_NS(ARM_Constants,volBase);
			}
		}		
		/// simple computation as the weighted average of year on year volatilities with df over annuity!
		else
		{
			double sumDFSigma = 0.0, sumDF = 0.0;

			for(j=0; j<tenors->size(); ++j )
			{
				sumDFSigma += CPILegDF->Elt(j) * CPIFwdVol_ZC->Elt(j);
				sumDF	   += CPILegDF->Elt(j);
				volatilities->Elt(i,j) = sumDFSigma/sumDF * CC_NS(ARM_Constants,volBase);;
			}
		}
	}

	/// Set the data to the volatilities matrix
	ARM_Vector* expiriesVec	= To_pARM_Vector( expiries );
	delete expiries;
	ARM_Vector* tenorsVec	= To_pARM_Vector( tenors );
	delete tenors;

	ARM_VolCurve* volCurve = new ARM_VolLInterpol(asOfDate, expiriesVec, tenorsVec, volatilities, K_STK_TYPE_PRICE, K_ATMF_VOL );
	volCurve->SetVolatilities(volatilities);

	HoldVolatilities.Release();
	return volCurve;
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
///	Class      : ARM_InfOATVolComputation_Producer_EqualWeight
////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_InfOATVolComputation_Producer_EqualWeight
///	Routine: Constructor
////////////////////////////////////////////////////

ARM_InfOATVolComputation_Producer_EqualWeight::ARM_InfOATVolComputation_Producer_EqualWeight(ARM_InfBSModel* infIRBSModel)
:	ARM_InfVolComputation_Producer(infIRBSModel)
{
	SetName(ARM_INFSWWOPTVOL_EQUALWEIGHT_PRODUCER);
}


////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_EqualWeight
///	Routine: Copy Constructor
////////////////////////////////////////////////////
ARM_InfOATVolComputation_Producer_EqualWeight::ARM_InfOATVolComputation_Producer_EqualWeight(
	const ARM_InfOATVolComputation_Producer_EqualWeight& rhs )
:	ARM_InfVolComputation_Producer(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_InfOATVolComputation_Producer_EqualWeight
///	Routine: Assignement operator
////////////////////////////////////////////////////
ARM_InfOATVolComputation_Producer_EqualWeight& ARM_InfOATVolComputation_Producer_EqualWeight::operator=( 
	const ARM_InfOATVolComputation_Producer_EqualWeight& rhs )
{
	if( this != &rhs )
		ARM_InfVolComputation_Producer::operator =( rhs );
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfOATVolComputation_Producer_EqualWeight
///	Routine: Destructor
////////////////////////////////////////////////////
ARM_InfOATVolComputation_Producer_EqualWeight::~ARM_InfOATVolComputation_Producer_EqualWeight()
{}

////////////////////////////////////////////////////
///	Class  : ARM_InfOATVolComputation_Producer_EqualWeight
///	Routine: Clone
////////////////////////////////////////////////////
ARM_Object* ARM_InfOATVolComputation_Producer_EqualWeight::Clone()
{
	return new ARM_InfOATVolComputation_Producer_EqualWeight(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_InfVolComputation_Producer_EqualWeight
///	Routine: Clone
////////////////////////////////////////////////////
string ARM_InfOATVolComputation_Producer_EqualWeight::toString( const string& indent )
{
	return "ARM_InfOATVolComputation_Producer_EqualWeight";
}



////////////////////////////////////////////////////
///	Class  : ARM_InfOATVolComputation_Producer_EqualWeight
///	Routine: VolZC_to_VolOATSwp
///	Returns: double
///	Action : general function to compute the volatility of inf swaption from the one 
///				of zero coupon 
////////////////////////////////////////////////////

ARM_VolCurve* ARM_InfOATVolComputation_Producer_EqualWeight::GenerateOATSwopVolCurve( 
	const ARM_Date&		asOfDate,
	ARM_GP_Vector*		tenors,
	ARM_GP_Vector*		expiries )
{
	GetDefaultDataIfMissing( tenors, expiries );
	ARM_Matrix* volatilities= new ARM_Matrix( expiries->size(), tenors->size() );
	ARM_AutoCleaner< ARM_Matrix > HoldVolatilities( volatilities );

	/// it is assumed that the convention on the CPI frequency is K_ANNUAL
	double CPIFreq		= K_ANNUAL;
	double CPIPeriod	= 1.0/CPIFreq;

	/// data for the computation of the volatilities
	ARM_GP_Vector* pVol_ZC			= NULL; 
	double lookupStrike				= 0;
	double denomYearFrac, denomYearFrac2;
	double asOf = itsInfBSModel->GetZeroCurve()->GetAsOfDate().GetJulian(), TVolLookup;
	ARM_Date lastKnownDate	= itsInfBSModel->GetVolatility()->GetLastKnownDate(), denomDate;
	int infDayCount		= itsInfBSModel->GetInfFwdCurv()->GetMonthlyInterpType();

	/// 1) computes the volatility using the relationship
	for( size_t i=0; i<expiries->size(); ++i )
	{
		/// computes all the IR data for the longest tenors
		double longestTenor			= tenors->Elt(tenors->size()-1);

		/// computes all the data for the CPIs
		size_t CPIDataSize			= longestTenor*CPIFreq;

		ARM_GP_Vector* CPIFwdVol_ZC	= new ARM_GP_Vector(CPIDataSize);
		ARM_AutoCleaner< ARM_GP_Vector > HoldCPIFwdVol_ZC( CPIFwdVol_ZC );

		for( size_t j=0; j<CPIDataSize; ++j )
		{
			denomYearFrac			= expiries->Elt(i)+j*CPIFreq;
			//denomDate				= ARM_Date(asOf + denomYearFrac * K_YEAR_LEN);
			//TVolLookup				= CountYearsWithoutException( infDayCount, lastKnownDate, denomDate );
			//CPIFwdVol_ZC->Elt(j)		= itsInfBSModel->GetVolatility()->ComputeVolatility( asOf, lookupStrike, TVolLookup )/ CC_NS( ARM_Constants, volBase );

			/// compute the fwd ZC vols as a sum of quadratic average of the year on year vols
			for( size_t k=0; k<=j; ++k )
			{
				denomYearFrac2		= expiries->Elt(i)+k*CPIFreq;
				TVolLookup			= itsInfBSModel->GetVolatility()->ComputeVolatility( denomYearFrac2, lookupStrike, CPIFreq)/ CC_NS( ARM_Constants, volBase );
				CPIFwdVol_ZC->Elt(j)	+= denomYearFrac2*TVolLookup*TVolLookup;
			}
			CPIFwdVol_ZC->Elt(j) = sqrt(CPIFwdVol_ZC->Elt(j)/denomYearFrac);
			
		}

		double sum = 0;
		/// start calling the function to compute volatility with correct data
		for(j=0; j<tenors->size(); ++j )
		{
			sum += CPIFwdVol_ZC->Elt(j);
			volatilities->Elt(i,j) = sum/(j+1) *CC_NS(ARM_Constants,volBase);;
		}
	}

	/// Set the data to the volatilities matrix
	ARM_Vector* expiriesVec	= To_pARM_Vector( expiries );
	delete expiries;
	ARM_Vector* tenorsVec	= To_pARM_Vector( tenors );
	delete tenors;

	ARM_VolCurve* volCurve = new ARM_VolLInterpol(asOfDate, expiriesVec, tenorsVec, volatilities, K_STK_TYPE_PRICE, K_ATMF_VOL );
	volCurve->SetVolatilities(volatilities);
	HoldVolatilities.Release();
	return volCurve;
}
CC_END_NAMESPACE()

/*---------------------------------------------------------------*/
/*---- End Of File ----*/
