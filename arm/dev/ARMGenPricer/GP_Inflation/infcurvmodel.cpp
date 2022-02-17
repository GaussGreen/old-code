/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file infbsmodel.h
 *  \brief Black Scholes model for the inflation
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date September 2003
 */

/// gpinflation
#include "gpinflation/infcurvmodel.h"
#include "gpinflation/infcurv.h"

/// gpbase
#include "gpbase/rateconversion.h"
#include "gpbase/gpvector.h"


CC_BEGIN_NAMESPACE( ARM )

///////////////////////////////////////////
/// Helper class to avoid code duplication
/// between models
///////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : InfFwdModel
///	Routine: FwdCPI
///	Returns: double
///	Action : computation of the forward CPI
////////////////////////////////////////////////////
double InfFwdModel::FwdCPI(
	const ARM_Date& resetDate, 
	long dailyInterpType,
	const string& CPILag,
	const string& DCFLag,
	ARM_InfIdx* infIdx )
{
	return itsInfFwdCurve->CPIInterpolate( resetDate, DCFLag, dailyInterpType, CPILag );
}

 
////////////////////////////////////////////////////
///	Class  : InfFwdModel
///	Routine: FwdCPI
///	Returns: double
///	Action : computation of the forward CPI
////////////////////////////////////////////////////
double InfFwdModel::FwdCPI(
	const ARM_Date& resetDate, 
	long dailyInterpType,
	ARM_InfIdx* infIdx	)
{
	return itsInfFwdCurve->CPIInterpolate( resetDate, resetDate, dailyInterpType );
}


////////////////////////////////////////////////////
///	Class  : InfFwdModel
///	Routine: GetCPIIndexValue
///	Returns: double
///	Action : computation of the CPI Index Value
////////////////////////////////////////////////////
double InfFwdModel::GetCPIIndexValue() const
{
	return itsInfFwdCurve->GetCPIIndexValue();
}


////////////////////////////////////////////////////
///	Class  : InfFwdModel
///	Routine: Constructor
///	Returns: 
///	Action : the inflation curve is clone for safetyness
////////////////////////////////////////////////////
InfFwdModel::InfFwdModel( ARM_InfCurv*	infFwdCurv )
:	itsInfFwdCurve( infFwdCurv == NULL? NULL: (ARM_InfCurv*) infFwdCurv->Clone() )
{}

////////////////////////////////////////////////////
///	Class  : InfFwdModel
///	Routine: Copy Constructor
////////////////////////////////////////////////////
InfFwdModel::InfFwdModel( const InfFwdModel& rhs )
:	itsInfFwdCurve( rhs.itsInfFwdCurve == NULL? NULL: (ARM_InfCurv*) rhs.itsInfFwdCurve->Clone() )
{}



////////////////////////////////////////////////////
///	Class  : InfFwdModel
///	Routine: Assignment operator
////////////////////////////////////////////////////
InfFwdModel& InfFwdModel::operator=( const InfFwdModel& rhs )
{
	if( this != &rhs )
	{
		CleanUp();
		itsInfFwdCurve = rhs.itsInfFwdCurve == NULL? NULL: (ARM_InfCurv*) rhs.itsInfFwdCurve->Clone();
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : InfFwdModel
///	Routine: CleanUp
///	Returns: 
///	Action : delete appropriate pointor
////////////////////////////////////////////////////
void InfFwdModel::CleanUp()
{
	delete itsInfFwdCurve;
	itsInfFwdCurve = NULL;
}



////////////////////////////////////////////////////
///	Class  : InfFwdModel
///	Routine: Destructor
////////////////////////////////////////////////////
InfFwdModel::~InfFwdModel()
{
	CleanUp();
}



////////////////////////////////////////////////////
///	Class  : InfFwdModel
///	Routine: Set
///	Returns: 
///	Action : the inflation curve is clone for safetyness
////////////////////////////////////////////////////
void InfFwdModel::SetInfFwdCurv( ARM_InfCurv*	infFwdCurv )
{
	delete itsInfFwdCurve;
	itsInfFwdCurve = infFwdCurv;
}

////////////////////////////////////////////////////
///	Class  : InfFwdModel
///	Routine: GetInfFwdCurv
///	Returns: ARM_InfCurv*
///	Action : accessor method
////////////////////////////////////////////////////
ARM_InfCurv* InfFwdModel::GetInfFwdCurv( ) const
{
	return itsInfFwdCurve;
}



///////////////////////////////////////////
/// the real infCurvModel
///////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_InfCurvModel
///	Routine: Constructor
////////////////////////////////////////////////////
ARM_InfCurvModel::ARM_InfCurvModel(ARM_ZeroCurve* discountCurve, ARM_InfCurv* infFwdCurv )
:	ARM_Model( discountCurve ), 
	InfFwdModel( infFwdCurv )
{
	SetName(ARM_INFCURVMODEL);	
}


////////////////////////////////////////////////////
///	Class  : ARM_InfCurvModel
///	Routine: Copy Constructor
////////////////////////////////////////////////////
ARM_InfCurvModel::ARM_InfCurvModel(const ARM_InfCurvModel& rhs)
:	ARM_Model( rhs ), InfFwdModel( rhs )
{}

////////////////////////////////////////////////////
///	Class  : ARM_InfCurvModel
///	Routine: Operator=
////////////////////////////////////////////////////
ARM_InfCurvModel& ARM_InfCurvModel::operator = (const ARM_InfCurvModel &rhs )
{
	if( this !=	 &rhs )
	{
	    ARM_Model::operator = ( rhs );
	    InfFwdModel::operator = ( rhs );
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfCurvModel
///	Routine: Destructor
////////////////////////////////////////////////////
ARM_InfCurvModel::~ARM_InfCurvModel()
{}


////////////////////////////////////////////////////
///	Class  : InfFwdModel
///	Routine: Clone
///	Returns: 
///	Action : the inflation curve is clone for safetyness
////////////////////////////////////////////////////
ARM_Object* ARM_InfCurvModel::Clone()
{
	return new ARM_InfCurvModel(*this);
}



////////////////////////////////////////////////////
///	Class  : InfFwdModel
///	Routine: DiscountedCPI
///	Returns: 
///	Action : computes the discounted CPI
////////////////////////////////////////////////////
double ARM_InfCurvModel::DiscountedCPI(
	const ARM_Date& resetDate, 
	const ARM_Date& paymentDate, 
	long dailyInterpType,
	const string& CPILag,
	const string& DCFLag,
	ARM_InfIdx* infIdx )
{
	double fwdCPI		= FwdCPI( resetDate, dailyInterpType, CPILag, DCFLag );
    ARM_ZeroCurve* zc	= GetZeroCurve();
	double df			= zc->DiscountPrice( const_cast<ARM_Date&>(paymentDate) );
	return df*fwdCPI;
}




////////////////////////////////////////////////////
///	Class  : InfFwdModel
///	Routine: View
///	Returns: 
///	Action : view details about the object
////////////////////////////////////////////////////
void ARM_InfCurvModel::View(char* id, FILE* ficOut)
{
    FILE* fOut;
    char fOutName[200];
	
	/// first determine that the file is not already opened
    if ( ficOut == NULL )
    {
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w");
    }
    else
		fOut = ficOut;

	/// since the view method is common
	/// for all inflationswap
	/// get the type of the swap
    fprintf(fOut, "\n\n =======> INFLATION YC MODEL<====== \n" );
    
	fprintf(fOut, "Corresponding interest rate zero curve :\n" );
	GetZeroCurve()->View( id, fOut );

    fprintf(fOut, "Corresponding inflation forward CPI curve :\n" );
	GetInfFwdCurv()->View( id, fOut );

	/// to allow to have nested view
    if ( ficOut == NULL )
       fclose(fOut);
}



////////////////////////////////////////////////////
///	Class  : InfFwdModel
///	Routine: FwdCPIRatio
///	Returns: 
///	Action : Computes the fwd CPI Ratio
////////////////////////////////////////////////////
ARM_GP_Vector ARM_InfCurvModel::FwdCPIRatio( 
	const ARM_Date& numDate,
	const ARM_Date& denomDate,
	const ARM_Date& paymentDate,
	double multiple,
	double spread,
	long dailyInterpType,
	double denomFixing,
	ARM_InfIdx* infIdx )
{
	ARM_GP_Vector result(3);
	
	ARM_Date cutoffDate	= GetStartDate();
	double CPIdenom;
	/// we allow overwritting of fixing only in the past!
	if( denomDate < cutoffDate )
		CPIdenom = denomFixing == GETDEFAULTVALUE? FwdCPI( denomDate, dailyInterpType ) : denomFixing;
	else
		CPIdenom = FwdCPI( denomDate, dailyInterpType );

	double CPInum   = FwdCPI(numDate, dailyInterpType );
	double CPIRatio	= multiple * ( CPInum / CPIdenom ) + spread;

	result[0] = CPIRatio;
	result[1] = CPInum;
	result[2] = CPIdenom;
	return result;
}



////////////////////////////////////////////////////
///	Class  : ARM_InfCurvModel
///	Routine: ExpectedFwdYield
///	Returns: 
///	Action : Computes the fwd yield on the yield curve
////////////////////////////////////////////////////

double ARM_InfCurvModel::ExpectedFwdYield(
	ARM_Date& resetDate,
	ARM_Date& maturity,
	ARM_Date& payDate,
	int compMeth,
	int dayCount,
	int DomOrFrgRate,
	int discYC,
	int YieldDecomp,
	double Margin,
	StoreFwdRateInfo* StoreInfo,
	int IndexFreq)
{
    double reset	= (resetDate.GetJulian()-GetStartDate().GetJulian())/K_YEAR_LEN;
    double mat		= (maturity.GetJulian()-GetStartDate().GetJulian())/K_YEAR_LEN;
    double yieldMat = CountYears(dayCount, resetDate, maturity);
    double fwdYield = ForwardYield(0.0, reset, mat, yieldMat, compMeth, dayCount, DomOrFrgRate);
    double decapSpreadFwdYield = ConvertRateToRate((fwdYield+Margin)/CC_NS(ARM_Constants,rateBase), 1.0,K_COMP_ANNUAL, YieldDecomp) * CC_NS(ARM_Constants,rateBase);

    if (StoreInfo)
       StoreInfo->StoreFwdRate(fwdYield, 0.0);

    return(decapSpreadFwdYield);
}



////////////////////////////////////////////////////
///	Class  : ARM_InfCurvModel
///	Routine: ForwardYield
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

double ARM_InfCurvModel::ForwardYield(
	double calcDate, 
	double resetDate,
    double maturityDate, 
	double yieldMaturity, 
    int compMeth, 
	int dayCount, 
	int DomOrFrg )
{
    if (resetDate + K_NEW_DOUBLE_TOL < calcDate || maturityDate < resetDate )
        throw Exception(__LINE__, __FILE__, ERR_YEAR_TERMS_PAIR,
           "Settlement should be less than ForwardDate and ZeroMaturity");

    if( yieldMaturity <= 0.0 )
       return 0.0;

	/// special treatment of the livretA Curve
    if ( GetZeroCurve()->GetName() == ARM_LIVRETACURVE )
       return GetZeroCurve()->ForwardYield(calcDate, resetDate, 0);

	double zp1 = resetDate<=calcDate? 1.0 : ZeroPrice(calcDate, resetDate, DomOrFrg);
    double zp2 = ZeroPrice(calcDate, maturityDate, DomOrFrg);

	/// standard yield
	double fwdYield = CC_NS(ARM_Constants,rateBase)*log(zp1/zp2)/yieldMaturity;

	/// decompounding treatment
	switch( compMeth )
	{
	case K_COMP_CONT:
		break;
	case K_COMP_PROP:
		fwdYield = CC_NS(ARM_Constants,rateBase)*(exp(0.01*fwdYield*yieldMaturity)-1.0)/yieldMaturity;
		break;
	default:
		fwdYield = CC_NS(ARM_Constants,rateBase)*compMeth*(exp(0.01*fwdYield/compMeth)-1.0);
		break;
	}
    
	return (fwdYield);
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------*/
/*---- End Of File ----*/

