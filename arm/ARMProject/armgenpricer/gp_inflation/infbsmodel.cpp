/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *	\file infbsmodel.h
 *
 *  \brief Black Scholes model for the inflation
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date September 2003
 */

#include "gpinflation/infbsmodel.h"

/// gpbase
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/stringconvert.h"
#include "gpbase/gpvector.h"
#include "gpbase/gplinalgconvert.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/gpmatrixlinalg.h"
#include <gpbase/cloneutilityfunc.h>

/// gpinflation
#include "gpinflation/infcurv.h"
#include "gpinflation/infcapfloor.h"
#include "gpinflation/infswaption.h"
#include "gpinflation/infswopvol.h"
#include "gpinflation/infdata.h"
#include "gpinflation/sparsevolcube.h"
#include "gpinflation/infidx.h"
#include "gpinflation/assetinfo.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/spreadoption_lognormal_interface.h"
#include "gpclosedforms/spreadoption_lognormal_formula.h"
#include "gpclosedforms/spreadoption_shiftedlognormal_interface.h"
#include "gpclosedforms/vanilla_shifted_lognormal.h"
#include "gpclosedforms/vanille_normal_interface.h"
#include "gpclosedforms/vanilla_normal.h"


/// kernel
#include <mod/bsconvadjust.h>
#include <crv/correlmanager.h>
#include <crv/volcube.h>

/// STL
#include <functional>



CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_InfBSModel::ARM_InfBSModel(	
	const ARM_Date&		asOfDate,
	ARM_ZeroCurve*		discountCurv, 
	ARM_InfCurv*		infFwdCurv,
	ARM_VolCurve*		infCapVolCurv,
	ARM_CorrelManager*	correlManager,
	ARM_BSModel*		IRModel,
	ARM_VolCurve*		infSwoptCurve,
	ARM_VolCurve*		IRSwoptCurve )
:
	ARM_BSModel( discountCurv, infCapVolCurv, K_PRICE ),
	InfOptionModel( infFwdCurv ),
	itsCorrelManager( correlManager? (ARM_CorrelManager*) correlManager->Clone() : NULL ), 
	itsIRModel(	IRModel? (ARM_BSModel*) IRModel->Clone() : NULL ),
	itsInfSwoptVolCurve( infSwoptCurve? (ARM_VolCurve*) infSwoptCurve->Clone() : NULL ),
	itsIRSwoptVolCurve( IRSwoptCurve? (ARM_VolCurve*) IRSwoptCurve->Clone() : NULL ),
	itsIRSwoptVolCube ( NULL ),
	isCorrelMatrixValidated(false),
	isGoodCorrelMatrix(true)
{
	/// the asOfDate given in the constructor is set in the model
	/// GetStartDate should therefore provide the model AsOfDate
    SetStartDate( const_cast<ARM_Date&>(asOfDate) );

	// By default we use a BS conv adjust manager, it will be destroyed
	// with the model
    SetConvAdjustManager(new ARM_BSConvAdjust, false);

	SetName(ARM_INFBSMODEL);
};


////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_InfBSModel::ARM_InfBSModel(
	const ARM_Date&		asOfDate,
	ARM_ZeroCurve*		discountCurv, 
	ARM_InfCurv*		infFwdCurv,
	ARM_VolCurve*		infCapVolCurv,
	ARM_CorrelManager*	correlManager,
	ARM_BSModel*		IRModel,
	ARM_VolCurve*		infSwoptCurve,
	ARM_VolCube*		IRSwoptCube)
	:
ARM_BSModel( discountCurv, infCapVolCurv, K_PRICE ),
	InfOptionModel( infFwdCurv ),
	itsCorrelManager(		correlManager?	(ARM_CorrelManager*)	correlManager->Clone()	: NULL ), 
	itsIRModel(				IRModel?		(ARM_BSModel*)			IRModel->Clone()		: NULL ),
	itsInfSwoptVolCurve(	infSwoptCurve?	(ARM_VolCurve*)			infSwoptCurve->Clone()	: NULL ),
	itsIRSwoptVolCube(		IRSwoptCube?	(ARM_VolCube*)			IRSwoptCube->Clone()	: NULL ),
	itsIRSwoptVolCurve( NULL ),
	isCorrelMatrixValidated(false),
	isGoodCorrelMatrix(true)
{
	/// the asOfDate given in the constructor is set in the model
	/// GetStartDate should therefore provide the model AsOfDate
    SetStartDate( const_cast<ARM_Date&>(asOfDate) );

	// By default we use a BS conv adjust manager, it will be destroyed
	// with the model
    SetConvAdjustManager(new ARM_BSConvAdjust, false);

	SetName(ARM_INFBSMODEL);
};
	
////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

void ARM_InfBSModel::ValidateCorrelManagerForVarNotionalSwaption( ARM_InfIdx* infIdx, ARM_IRIndex* otherIndex ) 
{
	ARM_CorrelMatrix* IRSwpCorrelMatrix = GetIRSwpCorrelMatrix( otherIndex, "IR/IR_SWOPT") ;
	ARM_CorrelMatrix* InfIRCorrelMatrix = GetInfIRCorrelMatrix(infIdx,otherIndex, "INF/IR_SWOPT"  );
	ARM_CorrelMatrix* InfSwpCorrelMatrix = GetInfSwpCorrelMatrix( infIdx,  "INF/INF_SWOPT" ) ;


	ARM_GP_Matrix  IRSwpMatrix = To_ARM_GP_Matrix(IRSwpCorrelMatrix->GetVolatilities());
	ARM_GP_Matrix  InfIRMatrix = To_ARM_GP_Matrix(InfIRCorrelMatrix->GetVolatilities());
	ARM_GP_Matrix  InfSwpMatrix = To_ARM_GP_Matrix(InfSwpCorrelMatrix->GetVolatilities());

	ARM_GP_Matrix TestedMatrix = *((ARM_GP_Matrix *) InfSwpMatrix.Clone()); 
	ARM_GP_Matrix InfIRMatrixTranspose = *((ARM_GP_Matrix *) InfIRMatrix.Clone()); 
	InfIRMatrixTranspose = InfIRMatrixTranspose.transpose();

	int infircolsize = InfIRMatrix.GetColsNb();
	for(int j = 0; j<infircolsize; j++)
	{
		ARM_GP_Vector colVect = (*(InfIRMatrix.GetColumn(j)));
		TestedMatrix.push_backColumn(const_cast<ARM_GP_Vector&> (colVect));
	}

	int infcolsize = IRSwpMatrix.GetColsNb();
	for( j = 0; j<infcolsize; j++)
	{
		ARM_GP_Vector colVect = (*(IRSwpMatrix.GetColumn(j)));
		InfIRMatrixTranspose.push_backColumn(const_cast<ARM_GP_Vector&> (colVect));
	}

	int infircolsizenew = InfIRMatrixTranspose.GetRowsNb();
	for(int i = 0; i<infircolsizenew; i++)
	{
		ARM_GP_Vector rowVect = (*(InfIRMatrixTranspose.GetRow(i)));
		TestedMatrix.push_backRow(const_cast<ARM_GP_Vector&> (rowVect));
	}

	int colNb = TestedMatrix.GetColsNb();
	ARM_GP_Vector eigenValues(colNb);
	ARM_GP_Matrix* outMatrix = NULL;

	outMatrix = ACPTransformation( &TestedMatrix, eigenValues );

	int eigenValuesSize = eigenValues.size();
	int k =0;
	while((k<eigenValuesSize) && (eigenValues[k]> (-K_NEW_DOUBLE_TOL)) )
		k++;

	if (k<eigenValuesSize)
		isGoodCorrelMatrix = false;

	isCorrelMatrixValidated = true;

	delete outMatrix;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: Copy Constructor
///	Returns: 
///	Action : Copy Constructor
////////////////////////////////////////////////////
ARM_InfBSModel::ARM_InfBSModel( const ARM_InfBSModel& rhs)
:	
	ARM_BSModel( rhs ), 
	InfOptionModel( rhs ),
	itsCorrelManager( rhs.itsCorrelManager == NULL? NULL : (ARM_CorrelManager*) rhs.itsCorrelManager->Clone() ), 
	itsIRModel( rhs.itsIRModel == NULL? NULL : (ARM_BSModel*) rhs.itsIRModel->Clone()),
	itsInfSwoptVolCurve( rhs.itsInfSwoptVolCurve == NULL? NULL : (ARM_VolCurve*) rhs.itsInfSwoptVolCurve->Clone() ),
	itsIRSwoptVolCurve( rhs.itsIRSwoptVolCurve== NULL ? NULL : (ARM_VolCurve*) rhs.itsIRSwoptVolCurve->Clone() ),
	itsIRSwoptVolCube( rhs.itsIRSwoptVolCube== NULL ? NULL : (ARM_VolCube*) rhs.itsIRSwoptVolCube->Clone() ),
	isCorrelMatrixValidated(rhs.isCorrelMatrixValidated),
	isGoodCorrelMatrix(rhs.isGoodCorrelMatrix)
{}


////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: Operator = 
///	Returns: 
///	Action : Operator = 
////////////////////////////////////////////////////
ARM_InfBSModel& ARM_InfBSModel::operator = (const ARM_InfBSModel &rhs )
{
	if( this !=	 &rhs )
	{
		CleanUp();
	    ARM_BSModel::operator=( rhs );
		InfOptionModel::operator=( rhs  );
		itsCorrelManager	= rhs.itsCorrelManager == NULL ? NULL : (ARM_CorrelManager*) rhs.itsCorrelManager->Clone();
		itsIRModel			= rhs.itsIRModel == NULL ? NULL : (ARM_BSModel*) rhs.itsIRModel->Clone() ;
		itsInfSwoptVolCurve = rhs.itsInfSwoptVolCurve == NULL ? NULL : (ARM_VolCurve*) rhs.itsInfSwoptVolCurve->Clone();
		itsIRSwoptVolCurve	= rhs.itsIRSwoptVolCurve == NULL ? NULL : (ARM_VolCurve*) rhs.itsIRSwoptVolCurve->Clone();
		itsIRSwoptVolCube	= rhs.itsIRSwoptVolCube == NULL ? NULL : (ARM_VolCube*) rhs.itsIRSwoptVolCube->Clone();
		isCorrelMatrixValidated = rhs.isCorrelMatrixValidated;
		isGoodCorrelMatrix = rhs.isGoodCorrelMatrix;
	}
	return *this;
}



////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: CleanUp
///	Returns: 
///	Action : routine to delete pointor
////////////////////////////////////////////////////
void ARM_InfBSModel::CleanUp() 
{
	delete itsCorrelManager;
	delete itsIRModel;
	delete itsInfSwoptVolCurve;
	delete itsIRSwoptVolCurve;
	delete itsIRSwoptVolCube;
}



////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: Dstructor 
///	Returns: 
///	Action : Destructor: itsInfFwdCurve is not deleted because it is shared
////////////////////////////////////////////////////
ARM_InfBSModel::~ARM_InfBSModel()
{
	CleanUp();
}



////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: Clone 
///	Returns: 
///	Action : Clone
////////////////////////////////////////////////////
ARM_Object* ARM_InfBSModel::Clone()
{
	return new ARM_InfBSModel( *this );
}


////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: View method 
///	Returns: 
///	Action : View method
////////////////////////////////////////////////////

void ARM_InfBSModel::View(char* id, FILE* ficOut)
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
    fprintf(fOut, "\n\n =======> INFLATION BS MODEL<====== \n" );

	fprintf(fOut, "Corresponding interest rate zero curve :\n" );
	GetZeroCurve()->View( id, fOut );

    fprintf(fOut, "Corresponding inflation forward CPI curve :\n" );
	GetInfFwdCurv()->View( id, fOut );


	fprintf(fOut, "Corresponding inflation Vol sparse Cube :\n" );
	GetVolatility()->View( id, fOut );

	

	if( itsCorrelManager )
    {
		fprintf(fOut, "Corresponding correl manager :\n" );
		itsCorrelManager->View( id, fOut );
	}

	if( itsIRModel )
    {
		fprintf(fOut, "Corresponding interest rate stochastic model:\n" );
		itsIRModel->View( id, fOut );
	}

	if( itsInfSwoptVolCurve )
    {
		fprintf(fOut, "Corresponding inflation swaption vol curve:\n" );
		itsInfSwoptVolCurve ->View( id, fOut );
	}
	
	if( itsIRSwoptVolCurve )
	{
		fprintf(fOut, "Corresponding Interest rate vol curve:\n" );
		itsIRSwoptVolCurve->View( id, fOut );
	}

	if( itsIRSwoptVolCube )
	{
		fprintf(fOut, "Corresponding Interest rate vol cube:\n" );
		itsIRSwoptVolCube->View( id, fOut );
	}
	
	/// to allow to have nested view
    if ( ficOut == NULL )
       fclose(fOut);
}


////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: DiscountedCPI
///	Returns: double
///	Action : DiscountedCPI
////////////////////////////////////////////////////

double ARM_InfBSModel::DiscountedCPI(
	const ARM_Date& resetDate, 
	const ARM_Date& paymentDate, 
	long dailyInterpType,
	const string& CPILag,
	const string& DCFLag,
	ARM_InfIdx* infIdx	)
{
	double fwdCPI		= FwdCPI( resetDate, dailyInterpType, CPILag, DCFLag );
    ARM_ZeroCurve* zc	= GetZeroCurve();
	double df			= zc->DiscountPrice( const_cast<ARM_Date&>(paymentDate) );
	return df*fwdCPI;
}


////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: GetModelTimeWPublishLag
///	Returns: ARM_Date
///	Action : GetModelDateWPublishLag adds on top of the date the publish lag corresponding to the index
////////////////////////////////////////////////////
ARM_Date ARM_InfBSModel::GetModelDateWPublishLag( const ARM_Date& date, ARM_InfIdx* infIdx )
{
	string indexName	= infIdx->GetIndexName();
	string publishLag	= InfData::GetPublishLag( indexName.c_str() );
	ARM_Date tmpDate( date );
	tmpDate.AddPeriod( publishLag );
	return tmpDate;
}


////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: GetModelTimeWPublishLag
///	Returns: double
///	Action : Computes the model time as the date + Publish Lag of the index - model As Of date
////////////////////////////////////////////////////
double ARM_InfBSModel::GetModelTimeWPublishLag( const ARM_Date& date, ARM_InfIdx* infIdx )
{
	int dayCount				= GetInfFwdCurv()->GetMonthlyInterpType();
	ARM_Date tmpDateWPublishLag	= GetModelDateWPublishLag( date, infIdx );
	return CountYearsWithoutException( dayCount, GetStartDate(), tmpDateWPublishLag	); 
}


////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: FwdCPIRatio
///	Returns: 
///	Action : Computes inflation interest rates correlation between two dates!
///				FwdCPIRatio method to compute the convexity adjustment 
///				\brief
///				\latexonly 
///				
///				\subsection{Intuition}
///				\The intuition is the following: 
///				
///				The forward $CPI\left( t,T_{i}\right) $ fixing at time $T_{i}$ is obviously
///				a martingale under its payment probability measure $Q_{T_{i}}$. Similarly,
///				for the forward $CPI\left( t,T_{j}\right)$ fixing at time $T_{j}$ under the
///				probability measure $Q_{T_{j}},$ but not $Q_{T_{i}}$.
///				
///				Consequently, the expected value under the probability measure $Q_{T_{i}}$
///				of the ratio of the two CPIs (with time $T_{i}>T_{j}$) has to take into
///				account various convexity adjustments\footnote{%
///				The forward of the ratio of CPI is not equal to the ratio of the forward
///				CPIs. One calls more or less improperly this adjustment a convexity
///				adjustment by extension from the one used in interest rates for various
///				change of measures like CMS and in-arrears.}:
///				
///				\begin{itemize}
///				
///				\item  $CPI\left( t,T_{j}\right) $ is not a martingale under the $Q_{T_{i}}$
///				measure. Hence it has to be adjusted to account for the change of measure
///				between $Q_{T_{j}}$ and $Q_{T_{i}}.$ This adjustment should intuitively
///				depend on the covariance between the forward interest bond volatility
///				(between $T_{j}$ and $T_{i}$) and the forward inflation rate in the
///				denominator $CPI\left( t,T_{j}\right) $. This change of measure is similar
///				to the CMS adjustment.
///				
///				\item  In addition, we pay $CPI\left( t,T_{i}\right) /CPI\left(
///				t,T_{j}\right) $.\ Because of the correlation between these two inflation
///				forward rates, we need to account for the their joint move.\ The adjustment
///				should intuitively be depending on the covariance between these two CPI\
///				rates.\ This is similar in a sense to a quanto adjustment.
///				\end{itemize}
///				
///				The adjustment is therefore computed in two steps:
///				
///				\begin{itemize}
///				\item  change of measure between $Q_{T_{j}}$ and $Q_{T_{i}}$.
///				
///				\item  computation of the expected ratio.
///				\end{itemize}
///				
///				\subsection{Assumptions}
/// 
///				The forward $CPI\left( t,T_{i}\right) $ fixing at time $T_{i}$ being a
///				martingale under its payment probability measure $Q_{T_{i}}$, its dynamic
///				is as follows: 
///				\begin{equation}
///				\frac{dCPI\left( t,T_{i}\right) }{CPI\left( t,T_{i}\right) }=\sigma \left(
///				t,T_{i}\right) dW_{T_{i}}\left( t\right) ,  \label{Mg_handwriting}
///				\end{equation}
///				
///				where $\left( W_{T_{i}}\left( t\right) \right) _{1\leq i\leq n}$\ is a $n$
///				-dimensional Brownian motion under $Q_{T_{i}}$. 
///				
///				The computation of the expectation of the forward CPI is given
///				by: 
///				\begin{equation*}
///				\mathbb{E}^{Q_{T_{i}}}\left[ \tfrac{CPI\left( T_{i},T_{i}\right) }{CPI\left(
///				T_{j},T_{j}\right) }\right] /\tfrac{CPI\left( 0,T_{j}\right) }{CPI\left(
///				0,T_{i}\right) }=\exp \left\{ 
///				\begin{array}{c}
///				\int_{0}^{T_{j}}\sigma \left( s,T_{j}\right) \left\{ \sigma \left(
///				s,T_{j}\right) -\rho _{j,i}^{Inf}\sigma \left( s,T_{i}\right) \right\} \\ 
///				+\rho _{i}^{B,I}\sigma \left( s,T_{j}\right) \left\{ \Gamma \left(
///				s,T_{i}\right) -\Gamma \left( s,T_{j}\right) \right\} ds
///				\end{array}
///				\right\} 
///				\end{equation*}
///				
///				\subsection{Assumption about the forward bond volatility}
///				
///				Because of the uncertainty on the estimation of the instantaneous
///				correlation between the forward bond and the CPI, we take as an input the
///				new integrated correlation $\zeta _{T_{j},T_{i}}$. In the case of constant
///				CPI\ volatility $\sigma \left( s,T_{j}\right) $ and forward bond volatility $%
///				\left( \Gamma \left( s,T_{i}\right) -\Gamma \left( s,T_{j}\right) \right) ,$
///				notice that the two correlations: the instantaneous $\rho _{i}^{B,I}$ and
///				integrated one $\zeta _{T_{j},T_{i}}$ are the same.
///				
///				The volatility of the forward bond can be read directly from the volatility
///				structure of the caplets.\ This comes from the fact that 
///				\begin{equation*}
///				B\left( t,T_{j},T_{i}\right) =\prod_{k=j..i-1}\frac{1}{1+\delta _{k}F\left(
///				t,T_{k},T_{k}+\delta _{k}\right) }, 
///				\end{equation*}
///				where $F\left( t,T_{k},T_{k}+\delta _{k}\right) $ is the forward libor of
///				period $\delta _{k}$ fixing at time $T_{k}$ and paid at time $T_{k}+\delta
///				_{k}$ and observed at time $t$.
///				
///				Applying Ito (and looking only at the stochastic part) leads immediately to 
///				\begin{equation*}
///				dB\left( t,T_{j},T_{i}\right) =\sum_{k=j..i-1}\frac{B\left(
///				t,T_{j},T_{i}\right) }{1+\delta _{k}F\left( t,T_{k},T_{k}+\delta _{k}\right) 
///				}\delta _{k}\sigma ^{F}\left( t,T_{k}\right) F\left( t,T_{k},T_{k}+\delta
///				_{k}\right) dB_{Q_{T_{i}}}^{k}\left( t\right) , 
///				\end{equation*}
///				where $\sigma ^{F}\left( t,T_{k}\right) $ is the lognormal volatility of the
///				forward libor $F\left( t,T_{k},T_{k}+\delta _{k}\right) $ and where the
///				diffusion is taken under the $Q_{T_{i}}$ probability measure. This means
///				that the forward bond volatility is approximately given by 
///				\begin{equation*}
///				\Gamma _{BS}\left( 0,T_{j},T_{i}\right) =\sum_{k=j..i-1}\frac{\delta
///				_{k}F\left( 0,T_{k},T_{k}+\delta _{k}\right) \sigma ^{F}\left(
///				0,T_{k}\right) }{1+\delta _{k}F\left( 0,T_{k},T_{k}+\delta _{k}\right) }, 
///				\end{equation*}
///				where we have approximated the forward by its current value.
///				
///				\endlatexonly
////////////////////////////////////////////////////

ARM_GP_Vector ARM_InfBSModel::FwdCPIRatio( 
	const ARM_Date& numDate,
 	const ARM_Date& denomDate,
	const ARM_Date& paymentDate,
	double multiple,
	double spread,
 	long dailyInterpType,
 	double denomFixing,
	ARM_InfIdx* infIdx )
{
	/// initialise as a vector of three elements
	ARM_GP_Vector result(3,0.0);
	ARM_Date lastKnownDate	= GetVolatility()->GetLastKnownDate();
	ARM_Date modelAsOfDate	= GetStartDate();
	string indexName		= infIdx->GetIndexName();
	
	/// we allow overwritting of fixing only in the past!
	/// are we in the past?
	if( denomDate < modelAsOfDate )
	{
		/// in the past, no convexity correction!
		double rawCPInum	= FwdCPI(numDate, dailyInterpType );
		double rawCPIdenom	= denomFixing == GETDEFAULTVALUE? FwdCPI( denomDate, dailyInterpType ) : denomFixing;
		double rawCPIRatio	= multiple * ( rawCPInum / rawCPIdenom ) + spread;
		result[0] = rawCPIRatio ;
		result[1] = rawCPInum;
		result[2] = rawCPIdenom;
		return result;
	}
	else
	{
		double rawCPInum	= FwdCPI( numDate, dailyInterpType );
		double rawCPIdenom	= FwdCPI( denomDate, dailyInterpType );
		double rawCPIRatio	= multiple * ( rawCPInum / rawCPIdenom ) + spread;


		/// right now use the strike for 0.0!
		double lookUpStrike = 0.0;
		int dayCount		= GetInfFwdCurv()->GetMonthlyInterpType();
		double TjVolLookup	= CountYearsWithoutException( dayCount, lastKnownDate, denomDate );
		double TjFromAsOf	= CountYearsWithoutException( dayCount, modelAsOfDate, denomDate ); 
		ARM_Date denomDatewPublishLag = GetModelDateWPublishLag( denomDate, infIdx );
		double TjForMaturity= CountYearsWithoutException( dayCount, modelAsOfDate, denomDatewPublishLag );

		double TiVolLookup	= CountYearsWithoutException( dayCount, lastKnownDate, numDate ); 
		double TiFromAsOf	= CountYearsWithoutException( dayCount, modelAsOfDate, numDate ); 
		double PaymentFromAsOf = CountYearsWithoutException( dayCount, modelAsOfDate, paymentDate ); 
		ARM_Date numDatewPublishLag = GetModelDateWPublishLag( numDate, infIdx );
		double TiForMaturity= CountYearsWithoutException( dayCount, modelAsOfDate, numDatewPublishLag );
		double tenor		= TiVolLookup-TjVolLookup;
		double tenorPayment = PaymentFromAsOf-TiForMaturity;
		
		/// the vol used is the one of the zero coupon market
		/// defined with a time to start from the sparse volcube!
		double timeToStart	= ARM_SparseVolCube::spotTime;
		double volTj		= GetVolatility()->ComputeVolatility( timeToStart, lookUpStrike , TjVolLookup )
			/ CC_NS( ARM_Constants, volBase ) ;
		double volTi		= GetVolatility()->ComputeVolatility( timeToStart, lookUpStrike , TiVolLookup )
			/ CC_NS( ARM_Constants, volBase ) ;
		double volYtY		= GetVolatility()->ComputeVolatility( TjVolLookup, lookUpStrike , tenor )
			/ CC_NS( ARM_Constants, volBase ) ;
		
		double RhoInfIR			= GetInfIRCorrel(TjForMaturity, TjForMaturity, infIdx, infIdx, "INF/IR" );
		double RhoInfIRPaymenti	= GetInfIRCorrel(TiForMaturity,TiForMaturity, infIdx, infIdx, "INF/IR" );
		double RhoInfIRPaymentj	= GetInfIRCorrel(TjForMaturity,TiForMaturity, infIdx, infIdx, "INF/IR" );

		double forward = 0.0, forwardPayment = 0.0;
		double volIR = 0.0, volIRPayment = 0.0;

		/// lookup of the volBond
		if( itsIRModel )
		{
			///// FORCED To const cast!
			forward	= itsIRModel->ExpectedFwdYield( denomDatewPublishLag, numDatewPublishLag, numDatewPublishLag )
				/ CC_NS( ARM_Constants, rateBase );
			/// use ATM vol for interest rates!
			volIR	= itsIRModel->GetVolatility()->ComputeVolatility( TjForMaturity, tenor) 
				/ CC_NS( ARM_Constants, volBase );

			if(paymentDate.GetJulian() < numDatewPublishLag.GetJulian())
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					ARM_USERNAME + ": resetNumgap have to be greater than 2 months");
			}

			/// for payment lag shorter than 6 month, we do not compute it!
			const double tenorNoPaymentLagThreshold = 7.0;
			if( (paymentDate.GetJulian()-numDatewPublishLag.GetJulian()) >= tenorNoPaymentLagThreshold )
			{				/// use ATM vol for interest rates!
				volIRPayment = itsIRModel->GetVolatility()->ComputeVolatility( TiForMaturity , tenorPayment ) 
					/ CC_NS( ARM_Constants, volBase );
				forwardPayment = itsIRModel->ExpectedFwdYield( numDatewPublishLag, const_cast<ARM_Date&>( paymentDate ), const_cast<ARM_Date&>( paymentDate ) )
				/ CC_NS( ARM_Constants, rateBase );

			}
			else
			{
				tenorPayment = 0;
				volIRPayment = 0;
			}

		}
		
		double adj = 1.0;
		ARM_GP_Vector Input(14);
		if(fabs(numDatewPublishLag.GetJulian()-denomDatewPublishLag.GetJulian()-K_YEAR_LEN) < 7.0)
		{			
			Input[0] = TjForMaturity;
			Input[1] = TiForMaturity;
			Input[2] = tenor;
			Input[3] = volTj;
			Input[4] = volTi;
			Input[5] = volYtY;
			Input[6] = forward;
			Input[7] = volIR;
			Input[8] = RhoInfIR;
			
			/// part for the payment lag
			Input[9]  = tenorPayment;
			Input[10] = forwardPayment;
			Input[11] = volIRPayment;
			Input[12] = RhoInfIRPaymentj;
			Input[13] = RhoInfIRPaymenti;			
		}
		else
		{
			//New Calcul For volYtY from the CPI Vols
			string ccy			= infIdx->GetCurrencyUnit()->GetCcyName();
			string indexName	= GetInfFwdCurv()->GetInfIdxName();
			string intraMktTag	= ccy < indexName ? ccy + "_" + indexName : indexName + "_"+ ccy;				
			double RhoInfInf = 1.0;
			if(itsCorrelManager)
			{
				if(itsCorrelManager->GetCorrelData("INF/INF",intraMktTag) )
					RhoInfInf = GetInfInfCorrel(TjForMaturity, TjForMaturity, infIdx, infIdx, "INF/INF" );
			}
			volYtY	= TiForMaturity*volTi*volTi-2*RhoInfInf*sqrt(TiForMaturity*TjForMaturity)*volTi*volTj+TjForMaturity*volTj*volTj;

			//new cacul for VolIR as in Vol Swap Vol FRA
			if(itsIRModel)
			{
				double volBond	= 0.0;
				double YF_TERM =1.0;
				int k=0;
				double T_kmoins1 = TjForMaturity;
				ARM_Date date_kmoins1 = denomDatewPublishLag;
				double T_k = T_kmoins1;
				ARM_Date date_k = date_kmoins1;
				date_k.AddYears(YF_TERM);
				while (date_k.GetJulian() <numDatewPublishLag.GetJulian()+7.0)
				{			
					double delta_k = YF_TERM;
					double forward_kmoins1	= itsIRModel->ExpectedFwdYield( date_kmoins1, date_k, date_k )/ CC_NS( ARM_Constants, rateBase );
					double yf_kmoins1 = (date_kmoins1.GetJulian()-modelAsOfDate.GetJulian())/K_YEAR_LEN;
					double volIR_kmoins1	= itsIRModel->GetVolatility()->ComputeVolatility( yf_kmoins1 , delta_k)/ CC_NS( ARM_Constants, volBase );
					double RhoInfIR_kmoins1 = GetInfIRCorrel(TjForMaturity, yf_kmoins1, infIdx, infIdx, "INF/IR" );

					volBond+= delta_k*forward_kmoins1/(1.0+delta_k*forward_kmoins1)*volIR_kmoins1*RhoInfIR_kmoins1;

					date_kmoins1 = date_k;
					k+=1;
					date_k.AddYears(YF_TERM);
				}

				if( ( numDatewPublishLag.GetJulian()-date_k.GetJulian() ) > 7.0 )
				{
					date_kmoins1 = date_k;
					k+=1;
					date_k=numDatewPublishLag;
					double delta_k = (date_k.GetJulian()-date_kmoins1.GetJulian())/K_YEAR_LEN;
					double forward_kmoins1	= itsIRModel->ExpectedFwdYield( date_kmoins1, date_k, date_k )/ CC_NS( ARM_Constants, rateBase );
					double yf_kmoins1 = (date_kmoins1.GetJulian()-modelAsOfDate.GetJulian())/K_YEAR_LEN;
					double volIR_kmoins1	= itsIRModel->GetVolatility()->ComputeVolatility( yf_kmoins1 , delta_k)/ CC_NS( ARM_Constants, volBase );
					double RhoInfIR_kmoins1 = GetInfIRCorrel(TjForMaturity, yf_kmoins1, infIdx, infIdx, "INF/IR" );

					volBond+= delta_k*forward_kmoins1/(1.0+delta_k*forward_kmoins1)*volIR_kmoins1*RhoInfIR_kmoins1;
				}			
				RhoInfIR	= 1.0;
				volIR = ( 1.0 + tenor * forward )* volBond /(tenor * forward );
			}
		
			Input[0] = TjForMaturity;
			Input[1] = TiForMaturity;
			Input[2] = tenor;
			Input[3] = volTj;
			Input[4] = volTi;			
			Input[5] = sqrt(volYtY/tenor);

			Input[6] = forward;			
			Input[7] = volIR;
			Input[8] = RhoInfIR;
			
			/// part for the payment lag
			Input[9]  = tenorPayment;
			Input[10] = forwardPayment;
			Input[11] = volIRPayment;
			Input[12] = RhoInfIRPaymentj;
			Input[13] = RhoInfIRPaymenti;	
		}

		ARM_Vector tmpInput = To_ARM_Vector(Input);
		adj = GetConvAdjustManager()->FwdCPIRatioAdjust(this, &tmpInput);
		double CPIRatio		= multiple * ( rawCPInum / rawCPIdenom * adj ) + spread;

		/// stores result and returns
		result[0] = CPIRatio ;
		result[1] = rawCPInum;
		result[2] = rawCPIdenom;
		return result;
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: ComputeVolAndStrikeForCap
///	Returns: void
///	Action : Computes the total vol, the pricing strike and the tenor
////////////////////////////////////////////////////
void ARM_InfBSModel::ComputeVolAndStrikeForCap( 
	/// standard part
	double CPIForward,
	double strike,
	int callput,
	ARM_InfIdx* infIdx,			
	StoreInfoObj& storeInfo,
	/// part specific to cap
	const ARM_Date& numDate,	
	const ARM_Date& denomDate,	
	int optionType,				
	double renormalisationFactor,
	/// part to overwrite
	double& totalVol,
	double& pricingStrike,
	double& tenor )
{
	/// 1) computes the tenor
	int dayCount			= GetInfFwdCurv()->GetMonthlyInterpType();
	ARM_Date lastKnownDate	= GetVolatility()->GetLastKnownDate();
	ARM_Date modelAsOfDate	= GetStartDate();
	string indexName		= infIdx->GetIndexName();

	double TjVolLookup		= CountYearsWithoutException( dayCount, lastKnownDate, denomDate ); 
	ARM_Date denomDatewPublishLag = GetModelDateWPublishLag( denomDate, infIdx );
	double TjForMaturity	= CountYearsWithoutException( dayCount, modelAsOfDate, denomDatewPublishLag );

	double TiVolLookup		= CountYearsWithoutException( dayCount, lastKnownDate, numDate ); 
	double TiFromAsOf		= CountYearsWithoutException( dayCount, modelAsOfDate, numDate ); 
	ARM_Date numDatewPublishLag = GetModelDateWPublishLag( numDate, infIdx );
	double TiForMaturity	= CountYearsWithoutException( dayCount, modelAsOfDate, numDatewPublishLag );
	tenor					= TiVolLookup-TjVolLookup;

	/// test past cash flows or negative tenor
	if( TiForMaturity < 0.0 || tenor < 0.0 )
	{
		double data[2] = { -1.0, -1.0 };
		storeInfo.Store( data );
		tenor = -1;
		return ;
	}
	double leverage = 1.0;

	/// 2) get the appropriate pricingStrike
	switch( optionType )
	{
	case K_ZEROCOUPON_LEG:
		pricingStrike = pow( 1 + strike/CC_NS( ARM_Constants, rateBase ), tenor );
		break;
	case K_YEARTOYEAR_LEG:
		{
			StoreVolInfo* storeInfoVol = (StoreVolInfo*) (&storeInfo);

			ARM_InfCapFloor* infCapFloor = storeInfoVol->GetInfCapFloor();
			ARM_InfLeg* infLeg = (ARM_InfLeg*) (infCapFloor->GetSwapLeg());
			double constant = infLeg->GetConstant();
			leverage = infLeg->GetMultiple();
			pricingStrike = strike/CC_NS( ARM_Constants, rateBase ) - constant;
		}
		
		break;
	case K_OATTYPE_LEG: 
		pricingStrike = strike/CC_NS( ARM_Constants, rateBase );
		break;
	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			ARM_USERNAME + ": unknown type of option... can either be a zeroCoupon or a Year to year.. Please advise");
	}



	/// 3) compute vol
	/// if the strike is negative, there is no option, hence
	/// no need to compute the volatility
	if( pricingStrike<0)
	{
		double data[2] = { -1.0, -1.0 };
		storeInfo.Store( data );
		return ;
	}

	/// every time the time to start is less than 1d
	/// change it to one day
	/// hence for a year to year optionlet that has
	/// already fixed its denominator
	/// we use the volatility of the zero coupon!
	
	if( TjForMaturity <= 0 )
	{
		tenor		   += TjForMaturity;
		TjForMaturity	= StringMaturityToYearTerm( "1d" );
		TjVolLookup		= StringMaturityToYearTerm( "1d" );
	}

	if( TjVolLookup <= 0 )
	{
		TjVolLookup		= StringMaturityToYearTerm( "1d" );
	}

	/// see the comments at the top of the function
	double lookUpStrike = pow( renormalisationFactor * pricingStrike/leverage, 1.0/tenor) -1.0;
	double vol			= GetVolatility()->ComputeVolatility( TjVolLookup, lookUpStrike , tenor )/CC_NS( ARM_Constants, volBase );

	/// in order to track pricing information
	double data[2];
	data[0] = vol*CC_NS( ARM_Constants, volBase );
	data[1] = lookUpStrike;
	storeInfo.Store( data );
	totalVol	= vol * sqrt( tenor );
}



////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: ComputeVolForSwaption
///	Returns: double
///	Action : Computes the vol of a swaption given a curve
////////////////////////////////////////////////////
double ARM_InfBSModel::ComputeVolForSwaption(
	double strike,
	double optionMaturity,
	double swaptionMaturity,
	ARM_VolCurve* volCurve,
	const string& volName )
{
	if( !volCurve )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			ARM_USERNAME + ": model does not have any swaption volatility for " + volName + ", please advise!");
	
	if(optionMaturity>0)
		return volCurve->ComputeVolatility( optionMaturity,strike,swaptionMaturity )/CC_NS( ARM_Constants, volBase );
	else 
		return 0.0;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: ComputeVolForSwaptionWithSmile
///	Returns: double
///	Action : Computes the vol of a swaption given a curve
////////////////////////////////////////////////////
double ARM_InfBSModel::ComputeVolForSwaptionWithSmile(
		double strike,
		double optionMaturity,
		double swaptionMaturity,
		ARM_VolCube* volCube,
		const string& volName )
{
	if( !volCube )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			ARM_USERNAME + ": model does not have any swaption volatility for " + volName + ", please advise!");
	
	if(optionMaturity>0)
		return volCube->ComputeVolatility(optionMaturity,strike*CC_NS( ARM_Constants, volBase ),swaptionMaturity)/CC_NS( ARM_Constants, volBase );
		//VolatilityFunctionByStrike(optionMaturity, strike, swaptionMaturity)/CC_NS( ARM_Constants, volBase );		
	else 
		return 0.0;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: FwdOptionPrice
///	Returns: double
///	Action : this routines prices a caplet
///			in the case of an option which has already fixed	
/// 		we need to get the accurate strike to account for the 	
/// 		smile. The strikeRenormalisation factor is used here	
/// 
/// 		Let us explain how we compute it:	
/// 		in the case of a zero coupon, we assume that the convention 	
/// 		is to use a strike given by lookupStrike = pricingStrike^(1/tenor)-1	
/// 		hence for a 6 month zero coupon option, tenor = 0.5	
/// 		and the strike is precisely the strike since the pricing strike is given by	
/// 		pricingStrike = (1.0+strike)^tenor, hence lookupStrike = strike	
/// 
/// 		in the same situation for a year to year option which 	
/// 		has already fixed six months ago	
/// 		pricingStrike = (1.0+ strike) however, we need to account for the difference of fixing	
/// 		between the current fixing and the one previously	
/// 		therefore on top of the tenor effect, we have a renormalisation due to the fixing	
/// 		hence a renormalised strike of YtYFixingCPI/CurrentCPI*(1+strike),	
/// 		hence a lookupStrike = (YtYFixingCPI/CurrentCPI*(1+strike))^(1/Tenor)-1	
/// 		or in terms of the pricing strike = (renormalisationFactor *(pricingStrike)^(1/Tenor)-1	
/// 		slightly different from the original strike	
////////////////////////////////////////////////////

double ARM_InfBSModel::SingleAssetOptionPrice(
	double CPIForward,
	double strike,
	int callput,
	double discounting,
	ARM_InfIdx* infIdx,
	ARM_Object* optionContext,
	StoreInfoObj& storeInfo	)
{
	double totalVol = 0.0;
	double pricingStrike = 0.0;
	double tenor = 0.0;
	double vol = 0.0;

	ARM_InfCapFloorContext* infCapFloorContext = dynamic_cast<ARM_InfCapFloorContext*>(optionContext);
	ARM_InfSwaptionContext* infSwaptionContext = dynamic_cast<ARM_InfSwaptionContext*>(optionContext);

	if(infCapFloorContext)
	{
		ComputeVolAndStrikeForCap( 
			/// part standard to the option
			CPIForward,
			strike,
			callput,
			infIdx,
			storeInfo,
			
			/// part specific to the cap
			infCapFloorContext->GetNumDate(),
			infCapFloorContext->GetDenomDate(),
			infCapFloorContext->GetOptionType(),
			infCapFloorContext->GetRenormalisationFactor(),
			
			/// part to be overwritten
			totalVol,
			pricingStrike, 
			tenor );
	}
	else 
		if(infSwaptionContext)
		{
			tenor			= infSwaptionContext->GetOptionMaturity();
			pricingStrike	= strike;
			/// parcing: variable or constant notional 
			if (infSwaptionContext->HasBeenComputed1())
			{
				vol = infSwaptionContext->GetVolAsset1();
			}
			else
			{
				double strikeForVol = infSwaptionContext->GetStrikeForVol1()/CC_NS( ARM_Constants, rateBase );
				vol = ComputeVolForSwaption(strikeForVol,
										infSwaptionContext->GetExpiryAsset1(),
										infSwaptionContext->GetTenorAsset1(),
										itsInfSwoptVolCurve,
										"inflation" );
			}
			
			totalVol	= vol * sqrt(tenor);
			
			/// the information to store in the asset
			double data[10];
			data[0]= infSwaptionContext->GetFwdAsset1()*CC_NS( ARM_Constants, rateBase );
			data[1]= vol * CC_NS( ARM_Constants, volBase );
			data[2]= infSwaptionContext->GetTenorAsset1();
			data[3]= -1.0;
			data[4]= -1.0;
			data[5]= -1.0;
			data[6]= -1.0;
			data[7]= discounting;
			data[8]= infSwaptionContext->GetOptionMaturity();
			data[9]= infSwaptionContext->GetPricingStrike();
			storeInfo.Store( data );
		}
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"unknown type of option... can either be an inflation cap or an inflation swaption!");

	/// test for past deal
	if( tenor < 0 )
		return 0.0;

	/// test for zero strike
	if( pricingStrike == 0 )
		pricingStrike = K_NEW_DOUBLE_TOL;

	/// test for negative strike that have no sense
	if( pricingStrike < 0 )
	{
		/// if the strike is negative... try to return the intrinsic value
		/// if positive!
		double intrinsic = (CPIForward-pricingStrike) * callput;
		if( intrinsic >= 0 )
			return intrinsic * CC_NS( ARM_Constants, rateBase )*discounting;

		/// otherwise return an exception
		else
		{
			char msg[255];
			sprintf( msg, "%s: strike is equal to %f should be at least positive with the model assumption Probably you want to use a swap!.. Please advise",
				ARM_USERNAME.c_str(), pricingStrike );
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg );
		}
	}
	
	/// compute the BS value
	double result;

	/// parcing: variable or constant notional 
	if (infSwaptionContext && infSwaptionContext->HasBeenComputed1())
	{
		//result = discounting * Export_normal_VanillaOption(CPIForward, pricingStrike, vol, tenor, callput);
		result = BlackSholes_Formula( CPIForward, totalVol, discounting, pricingStrike, callput );
	}
	else
	{
		result = BlackSholes_Formula( CPIForward, totalVol, discounting, pricingStrike, callput );
	}

	return CC_NS( ARM_Constants, rateBase ) * result;	
	
}

///////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: Recovers Inf Flows from Inf Swp Rates
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
void ARM_InfBSModel::ComputeInfFlowsFromSwapRates(
		const ARM_GP_Vector&  dfTerms,
		const ARM_GP_Vector&  fixPayPeriod,
		const ARM_GP_Vector&  swpRates,
		ARM_GP_Vector& infFlows)
{
	infFlows[0] = swpRates[0]*fixPayPeriod[0];
	int i(0), size(swpRates.size());

	/// Size test
	if (size != dfTerms.size()) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_InfBSModel / Swap Rates and DFs must have the same size" );

	/*double level= dfTerms[0]*fixPayPeriod[0];
	double levelRatio = level;*/

	double level_imoins1= dfTerms[0]*fixPayPeriod[0];
	double level_i = level_imoins1;
	
	for (i=1; i<size; i++)
	{
		level_i = level_imoins1+dfTerms[i]*fixPayPeriod[i];
		infFlows[i] = (level_i*swpRates[i]-level_imoins1*swpRates[i-1])/dfTerms[i];
		level_imoins1 = level_i;
	}

}

///////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: Computes the sensitivities of the 
///			Inf notional swp rate 
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
void ARM_InfBSModel::ComputeInfSensiCoefs(
		const ARM_GP_Vector&  dfTerms,
		const ARM_GP_Vector&  fixPayPeriods,
		const ARM_GP_Vector&  infCoefs,
		ARM_GP_Vector&  swpRates,
		const double& shift,
		const double& initVlue,
		ARM_GP_Vector& sensiCoefs)
{
	int				size = swpRates.size();
	double			X(0.), numeraire(0.);
	ARM_GP_Vector	infFlows(size);

	for (int i=0; i<size; i++)
	{
		/// Bumps the slipping inf swap rates 
		swpRates[i] += shift;

		/// Recomputes the inf swap rates without notional
		ComputeInfFlowsFromSwapRates( dfTerms, fixPayPeriods, swpRates, infFlows);
		
		/// Reconstitutes the inf swp rate with notionel
		for (int j=0; j<size; j++)
		{
			X			+= infCoefs[j]		* infFlows[j];
			numeraire	+= fixPayPeriods[j] * infCoefs[j];
		}
		X /= numeraire;

		/// Reconstitutes the sensivity as a derivative 
		sensiCoefs[i] = (X - initVlue) / shift;

		/// Renormalizes
		swpRates[i] -= shift;
		X = numeraire = 0.;
	}
}

///////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: Computes Inf Swopt Vol in the case of 
///				Variable Notional
///	Returns: double
///	Action : Same principle as the IR basket swpot
////////////////////////////////////////////////////
double ARM_InfBSModel::ComputeInfVolWithVarNotional(
		const ARM_GP_Vector& communNotional,
		double floatStartTime,		
		const ARM_GP_Vector& floatEndTimes,
		const ARM_GP_Vector& floatPayTimes,
		const ARM_GP_Vector* dfTerms,	
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector& fwdRates,
		double strike,
		ARM_InfIdx* infIdx,
		ARM_Object* optionContext,
		ARM_GP_Vector& coeffs,
		ARM_GP_Vector& tenors,
		ARM_GP_Vector& vols)
{

	/// Variables Declaration
	int				size(communNotional.size()), i(0), j(0);
	coeffs.resize(size);
	tenors.resize(size);
	vols.resize(size);
	bool success;

	double			level(0.), infSum(0.), X(0.), X0(0.), numeraire0(0.), shift(0.), expiry(0.), basketVol(0.),
					nStDev(0.), strike_i(0.), basketFwd(0.), sqrTenor(0.);

	ARM_GP_Vector	infCoefs(size, 0.), fixCoefs(size, 0.), infFlows(size, 0.), df(size, 0.), swpRates(size, 0.),
					sensiCoefs(size, 0.), swpAtmNormVols(size, 0.), swpNormVols(size, 0.);

	ARM_GP_Matrix	infSwpCors(size, size, 0.);

	ARM_Date		numDate, denomDate, payDate, startDate(floatStartTime);
	
	/// Recovers the inf yoy flows, their coeffs, the slipping swap rates
	for (i=0; i<size; i++)
	{
		fixCoefs[i]	 = (*dfTerms)		[i] * fixPayPeriods[i];
		infCoefs[i]  = communNotional[i] * (*dfTerms)[i];
		infFlows[i]  = fabs(fwdRates[i]);

	}

	/// Stores slipping inf swp rates' correlations 
	for (i=0; i<size; i++)
	{
		double tenor_i = (floatEndTimes[i] - startDate.GetJulian()) / K_YEAR_LEN; //YK
		(tenors)[i] = tenor_i;

		for (j=0; j<size; j++)
		{
			double tenor_j = (floatEndTimes[j] - startDate.GetJulian()) / K_YEAR_LEN; //YK
			
			infSwpCors(i, j) = ARM_InfBSModel::GetInfSwpCorrel(tenor_i, tenor_j, infIdx, "INF/INF_SWOPT");
			if (infSwpCors(i, j)> 1.) infSwpCors(i, j) =  1.;
			if (infSwpCors(i, j)<-1.) infSwpCors(i, j) = -1.;
		}
	}

	/// Some valifation tests
	for (i=0; i<size; i++)
	{
		if (communNotional[i]<0.0)
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_InfBSModel / Variable Notional Swaption: Fixed leg is required to have only positive notionals");
			
		//numeraire0 += (*fixPayPeriods)[i] * communNotional[i] * (*dfTerms)[i];
		numeraire0 += fixPayPeriods[i] * infCoefs[i];
	}

	if (numeraire0 < 1e-12)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_InfBSModel /  Variable Notional Swaption: At least one fixed leg notionals should be > 0");

	
	/// Reconstitutes the slipping swap rates without notional 
	for (i=0; i<size; i++)
	{
		infSum		+= (*dfTerms)[i] * infFlows[i];
		level		+= fixCoefs[i];
		swpRates[i]	=  infSum / level;
	} 
	
	/// Recovers the inf yoy flows from slipping swap rates and dfs
	ComputeInfFlowsFromSwapRates( *dfTerms, fixPayPeriods, swpRates, infFlows);
	
	/// Constitutes the inf leg with notionel
	for (i=0; i<size; i++)
		X0 += infCoefs[i] * infFlows[i];

	/// Swap Rate := Inf Leg / Numeraire !
	X0		  /= numeraire0;

	/// Fwd of variable notional swap rate
	basketFwd  = X0;

	/// Computes the sensitivities of the notional swp rate 
	shift = 0.00001; /// 0.1 bp
	ComputeInfSensiCoefs(*dfTerms, fixPayPeriods, infCoefs, swpRates, shift, X0, sensiCoefs);
	

	/// Computes normal volatilities the slipping swap rates
	// First approximatoion
	ARM_InfSwaptionContext* infSwaptionContext = dynamic_cast<ARM_InfSwaptionContext*>(optionContext);
	expiry = infSwaptionContext->GetExpiryAsset1();
	sqrTenor = sqrt(expiry);

	for (i=0; i<size; i++)
	{
		if (fabs(sensiCoefs[i])>1e-15)
		{
			//double tenor = (floatEndTimes[i] - GetStartDate().DMYToJulian()) / K_YEAR_LEN;
			double tenor = (floatEndTimes[i] - startDate.DMYToJulian()) / K_YEAR_LEN;
			double vol	 = ARM_InfBSModel::ComputeVolForSwaption(strike, expiry, tenor, itsInfSwoptVolCurve, "inflation");
			swpAtmNormVols[i]	 = vol * (swpRates[i]+1.0);
		}	
	} 
	
	/// Computes variance of basket var
	for (i=0; i<size; i++)
	{
		if (fabs(sensiCoefs[i])>1e-15)
		{
			/// Diagonal terms
			double tmp	 = sensiCoefs[i] * swpAtmNormVols[i];
			basketVol	+= tmp * tmp;
			
			/// Crossed terms
			for (j=0; j<i; j++)
			{
				if (fabs(sensiCoefs[j])>1e-15)
				{
					basketVol += 2. * infSwpCors(i,j) 
									* sensiCoefs[i]	* swpAtmNormVols[i] 
									* sensiCoefs[j] * swpAtmNormVols[j];
				}
			}
			
		}	
	} 

	/// Deduces the basket vol 
	basketVol = sqrt(basketVol);
////////////////////////////////////////////////////////////////////////////////////////////
	// Computes strikes & build swpNormVols
	for (i=0; i<size; i++)
	{
		if(fabs(sensiCoefs[i])>1e-15)
		{
			nStDev = strike / (basketVol * sqrTenor);
			if (nStDev > 6.)  nStDev = 6.;
			if (nStDev < -6.) nStDev = -6.;
			
			if(sensiCoefs[i]<0) nStDev *= -1.;

			/// Computes strike
			strike_i = nStDev * swpAtmNormVols[i] * sqrTenor;
			if (strike_i < 1e-3) strike_i = 1e-3;

			double moyness_i = 100. * (strike_i - swpRates[i]);
			//double tenor	 = (floatEndTimes[i] - GetStartDate().DMYToJulian()) / K_YEAR_LEN;
			double tenor	 = (floatEndTimes[i] - startDate.DMYToJulian()) / K_YEAR_LEN;
			double vol(0.);
			
			// Computes lognormal vol @ strike
			vol = ARM_InfBSModel::ComputeVolForSwaption(
				strike_i,
				expiry,
				tenor,
				itsInfSwoptVolCurve,
				"inflation" );
						
			vol *= sqrTenor;

			/// Conversts into normal vol
			int callput		= (swpRates[i]>strike_i) ? -1 : 1;
			double bs		= BlackSholes_Formula( swpRates[i]+1.0, vol, 1., strike_i+1.0, callput);
			success = true;
			swpNormVols[i]	= VanillaImpliedVol_N( swpRates[i]+1.0, bs, strike_i+1.0, expiry, callput, &swpAtmNormVols[i], &success);
		}
	}
	
	// Re-initialize
	basketVol = 0.;

	// Second approximatoion: recomputes variance of basket var with exact ATM vols
	for (i=0; i<size; i++)
	{
		if (fabs(sensiCoefs[i])>1e-15)
		{
			/// Diagonal terms
			double tmp	 = sensiCoefs[i] * swpNormVols[i];
			basketVol	+= tmp * tmp;
			
			/// Crossed terms
			for (j=0; j<i; j++)
			{
				if (fabs(sensiCoefs[j])>1e-15)
				{
					basketVol += 2. * infSwpCors(i, j)
									* sensiCoefs[i] * swpNormVols[i]
									* sensiCoefs[j] * swpNormVols[j];
				}
			}
			
		}	
	} 

	/// Deduces the basket vol 
	basketVol = sqrt(basketVol);
	coeffs = sensiCoefs;	
	vols =swpNormVols;
///////////////////////////////////////////////////////////////////////
	/*coeffs = sensiCoefs;	
	vols =swpAtmNormVols;*/

	/// Test on vol level
	if (basketVol<1e-15) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_InfBSModel / Variable Notional Swaption: unresolved problem");

	return basketVol;
}

///////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: Recovers DF Flows from IR Swp Rates
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
void ARM_InfBSModel::ComputeDiscountFactorsFromSwapRates (
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector& swapRates,
		double startDf, 
		ARM_GP_Vector& dfs)
{
	int size(swapRates.size());
	if (size != dfs.size()) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_InfBSModel / fixPayPeriods and swapRates must have the same size" );

	double swpRatio(0.), dfResult(0.);
	
	for (size_t i(0); i<size; i++)
	{
		if (i==0)
			dfs[i] = startDf / ( 1. + fixPayPeriods[i] * swapRates[i] ) ;
		else
		{
			swpRatio  = swapRates[i] / swapRates[i-1];
			dfResult  = swpRatio * dfs[i-1] + (1. - swpRatio ) * startDf;
			dfResult  /= 1. + fixPayPeriods[i] * swapRates[i];
			dfs[i]	  = dfResult;
		}
	}
}

///////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: Computes the Covariance IR/INF for amort Notional  
///			swap IR/INF
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

double ARM_InfBSModel::ComputeInfIrCovarWithVarNotional (
			ARM_GP_Vector& infcoeffs,
			ARM_GP_Vector& inftenors,
			ARM_GP_Vector& infvols,
			ARM_GP_Vector& ircoeffs,
			ARM_GP_Vector& irtenors,
			ARM_GP_Vector& irvols,
			ARM_InfIdx* infIdx,
			ARM_IRIndex* otherIndex)
{
	int infsize(infcoeffs.size());
	int irsize(ircoeffs.size());
	if ((infsize != inftenors.size()) || (infsize != infvols.size()) || (irsize != irtenors.size()) || (irsize != irvols.size()) ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_InfBSModel::ComputeInfIrCovarWithVarNotional : coeefs, tenors and vols must have the same size" );

	double covar = 0.0;

	// Second approximatoion: recomputes variance of basket var with exact ATM vols
	for (int i=0; i<infsize; i++)
	{
		if (fabs(infcoeffs[i])>1e-15)
		{
			/// Diagonal terms
			double inftmp	 = infcoeffs[i] * infvols[i];
			double inftenor_i = inftenors[i];			
			/// Crossed terms
			for (int j=0; j<irsize; j++)
			{
				double irtenor_j = irtenors[j];
				double inf_ir_correl_ij = GetInfIRCorrel(inftenor_i,irtenor_j, infIdx, otherIndex, "INF/IR_SWOPT" );
				if (fabs(ircoeffs[j])>1e-15)
				{
					double irtmp	 = ircoeffs[j] * irvols[j];
					
					covar += inf_ir_correl_ij*inftmp*irtmp;
				}
			}
			
		}	
	} 
	return covar;
}

///////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: Computes the sensitivities of the 
///			IR notional swp rate 
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
void ARM_InfBSModel::ComputeIRSensiCoefs(
		const double dfCoefStart,
		const double dfStart,
		const ARM_GP_Vector& fixCoefs,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector&  dfCoefs,
		ARM_GP_Vector&  swpRates,
		const double& shift,
		const double& initVlue,
		ARM_GP_Vector& sensiCoefs)
{
	int				size = swpRates.size(), i(0), j(0);
	double			X(0.), numeraire(0.);
	ARM_GP_Vector	dfs(size);
	
	for (i=0; i<size; i++)
	{
		// shift forward swap rate
		swpRates[i] += shift;
		
		ComputeDiscountFactorsFromSwapRates(fixPayPeriods, swpRates, dfStart, dfs);
			
		// compute X
		X += dfCoefStart * dfStart;
		for (j=0; j<size; j++)
			X += dfCoefs[j] * dfs[j];

		numeraire = 0.0;
		
		for (j=0; j<size; j++)
			numeraire += fixCoefs[j] * dfs[j];
			
		X /= numeraire;

		// compute VolCoeff
		sensiCoefs[i] =  (X - initVlue) / shift;
		
		// Renormalizes
		swpRates[i] -= shift;

		X = 0.;
	}

	
}

///////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: Computes IR Swopt Vol in the case of 
///				Variable Notional
///	Returns: double
///	Action : Same principle as the IR basket swpot
////////////////////////////////////////////////////
double ARM_InfBSModel::ComputeIRVolWithVarNotional(
		const ARM_GP_Vector& irNotional,
		double floatStartTime,
		double floatEndTime,
		const ARM_GP_Vector& floatStartTimes,
		const ARM_GP_Vector& floatEndTimes,
		const ARM_GP_Vector& IrFixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		double strike,
		ARM_IRIndex* irIdx,
		ARM_Object* optionContext,
		ARM_GP_Vector& coeffs,
		ARM_GP_Vector& tenors,
		ARM_GP_Vector& vols)
{
	/// Variables Declaration
	int				size(IrFixPayTimes.size()), i(0), j(0);

	double			level(0.), infSum(0.), X(0.), X0(0.), numeraire0(0.), shift(0.), expiry(0.), basketVol(0.), 
					dfCoefStart(0.), maturity(0.), dfStart(0.), nStDev(0.), basketFwd(0.), sqrTenor(0.), 
					strike_i(0.);
	bool success;

	ARM_GP_Vector	dfCoefs(size, 0.), fixCoefs(size, 0.), dfs(size, 0.), swpRates(size, 0.),
					sensiCoefs(size, 0.), swpAtmNormVols(size, 0.), swpNormVols(size, 0.);

	ARM_GP_Matrix	irSwpCors(size, size, 0.);

	ARM_Date		numDate, denomDate, payDate, startDate(floatStartTimes[0]);

	ARM_ZeroCurve*	zc  = GetZeroCurve();

	//ARM_IRIndex*  irIdx = zc->get>GetCurrencyUnit();
	
	/// Recovers the the coeffs of the DF terms
	dfCoefStart = irNotional[0];
	for (i=0; i<size-1; i++)
	{
		dfCoefs[i]  = irNotional[i+1] - irNotional[i];
	}
	dfCoefs[size-1] = - irNotional[size-1];
	

	/// Recovers the slipping swap rates
	maturity = floatStartTime - GetStartDate().DMYToJulian();
	dfStart = zc->DiscountPrice(maturity / K_YEAR_LEN);
	
	for (i=0; i<size; i++)
	{	
		maturity	= (IrFixPayTimes[i] - GetStartDate().DMYToJulian()) / K_YEAR_LEN;
		dfs[i]		=  zc->DiscountPrice(maturity);
		level		+= dfs[i] * fixPayPeriods[i];
		swpRates[i] = (dfStart - dfs[i]) / level;
		fixCoefs[i] = irNotional[i] * fixPayPeriods[i];
	}

	
	/// Stores slipping ir swp rates' correlations 
	for (i=0; i<size; i++)
	{
		//GetStartDate().DMYToJulian()
		double tenor_i = (floatEndTimes[i] - floatStartTime) / K_YEAR_LEN; //YK
		tenors[i] = tenor_i;

		for (j=0; j<size; j++)
		{
			double tenor_j = (floatEndTimes[j] - floatStartTime) / K_YEAR_LEN; //YK
						
			irSwpCors(i, j) = ARM_InfBSModel::GetIRSwpCorrel(tenor_i, tenor_j, irIdx, "IR/IR_SWOPT");
			if (irSwpCors(i, j)> 1.) irSwpCors(i, j) =  1.;
			if (irSwpCors(i, j)<-1.) irSwpCors(i, j) = -1.;
		}
	}

	/// Recovers the df terms form the slipping swap rates
	ComputeDiscountFactorsFromSwapRates( fixPayPeriods, swpRates, dfStart, dfs);

	/// Some valifation tests
	for (i=0; i<size; i++)
	{
		if (irNotional[i]<0.0)
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_InfBSModel / Variable Notional Swaption: Fixed leg is required to have only positive notionals");
		
		numeraire0 += fixPayPeriods[i] * irNotional[i] * dfs[i];
	}

	if (numeraire0 < 1e-12)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_InfBSModel /  Variable Notional Swaption: At least one fixed leg notionals should be > 0");

	/// Constitutes the float leg with notionel
	X0 += dfCoefStart * dfStart;
	for (j=0; j<size; j++)
		X0 += dfCoefs[j] * dfs[j];

	/// Swap Rate := Folat Leg / Numeraire !
	X0 /= numeraire0;

	/// Fwd of variable notional swap rate
	basketFwd = X0;

	/// Computes the sensitivities of the notional swp rate 
	shift = 0.00001; /// 0.1 bp
	ComputeIRSensiCoefs( dfCoefStart, dfStart, fixCoefs, fixPayPeriods, dfCoefs, swpRates, shift, X0, sensiCoefs);

	/// Computes normal volatilities the slipping swap rates
	// First approximatoion
	ARM_InfSwaptionContext* infSwaptionContext = dynamic_cast<ARM_InfSwaptionContext*>(optionContext);
	expiry	  = infSwaptionContext->GetExpiryAsset2();
	sqrTenor = sqrt(expiry);

	double ATMMoyenness =0.0;

	for (i=0; i<size; i++)
	{
		if (fabs(sensiCoefs[i])>1e-15)
		{
			double tenor = (floatEndTimes[i] - startDate.DMYToJulian()) / K_YEAR_LEN;
			double vol (0.);
			if (itsIRSwoptVolCube!=NULL)
			{
				vol	 = ARM_InfBSModel::ComputeVolForSwaptionWithSmile(ATMMoyenness, expiry, tenor, itsIRSwoptVolCube,"interest rates" );
			}
			else
			{
				vol	 = ARM_InfBSModel::ComputeVolForSwaption(ATMMoyenness, expiry, tenor, itsIRSwoptVolCurve,"interest rates" );
			}

			swpAtmNormVols[i]	 = vol * swpRates[i];
		}	
	}

	/// Computes variance of basket var
	for (i=0; i<size; i++)
	{
		if (fabs(sensiCoefs[i])>1e-15)
		{
			/// Diagonal terms
			double tmp	 = sensiCoefs[i] * swpAtmNormVols[i];
			basketVol	+= tmp * tmp;
			
			/// Crossed terms
			for (j=0; j<i; j++)
			{
				if (fabs(sensiCoefs[j])>1e-15)
				{
					basketVol += 2. * irSwpCors(i, j)
									* sensiCoefs[i] * swpAtmNormVols[i] 
									* sensiCoefs[j] * swpAtmNormVols[j];
				}
			}
			
		}	
	} 

	// Deduces the basket vol 
	basketVol = sqrt(basketVol);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Computes strikes & build swpNormVols
	for (i=0; i<size; i++)
	{
		if(fabs(sensiCoefs[i])>1e-15)
		{
			nStDev = strike/ (basketVol * sqrTenor);
			if (nStDev > 6.)  nStDev = 6.;
			if (nStDev < -6.) nStDev = -6.;
			
			if(sensiCoefs[i]<0) nStDev *= -1.;

			/// Computes strike
			strike_i = swpRates[i] + nStDev * swpAtmNormVols[i] * sqrTenor;
			if (strike_i < 1e-3) strike_i = 1e-3;

			double moyness_i = (strike_i - swpRates[i]);
			double tenor	 = (floatEndTimes[i] - startDate.DMYToJulian()) / K_YEAR_LEN;
			double vol(0.);
			
			// Computes lognormal vol @ strike
			if (itsIRSwoptVolCube)
			{
				vol = ARM_InfBSModel::ComputeVolForSwaptionWithSmile(
					moyness_i,
					expiry,
					tenor,
					itsIRSwoptVolCube,
					"interest rates" );
			}
			else
			{
				vol = ARM_InfBSModel::ComputeVolForSwaption(
				strike,
				expiry,
				tenor,
				itsIRSwoptVolCurve,
				"interest rates" );
			}
			
			vol *= sqrTenor;			

			/// Conversts into normal vol
			int callput		= (swpRates[i]>strike_i) ? -1 : 1;
			double bs		= BlackSholes_Formula( swpRates[i], vol, 1., strike_i, callput);
			
			success = true;
			swpNormVols[i]	= VanillaImpliedVol_N( swpRates[i], bs, strike_i, expiry, callput, &swpAtmNormVols[i],&success);
		}
	}

	// Re-initialize
	basketVol = 0.;

	// Second approximatoion: recomputes variance of basket var with exact ATM vols
	for (i=0; i<size; i++)
	{
		if (fabs(sensiCoefs[i])>1e-15)
		{
			/// Diagonal terms
			double tmp	 = sensiCoefs[i] * swpNormVols[i];
			basketVol	+= tmp * tmp;
			
			/// Crossed terms
			for (j=0; j<i; j++)
			{
				if (fabs(sensiCoefs[j])>1e-15)
				{
					basketVol += 2. * irSwpCors(i, j)
									* sensiCoefs[i] * swpNormVols[i]
									* sensiCoefs[j] * swpNormVols[j];
				}
			}
			
		}	
	} 
	basketVol = sqrt(basketVol);

	coeffs = sensiCoefs;	
	vols =swpNormVols;

	/*coeffs = sensiCoefs;	
	vols =swpAtmNormVols;*/

	/// Deduces the basket vol 
	

////////////////////////////////////////////////////////////////////////////////////////
	/// Test on vol level
	if (basketVol<1e-15) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_InfBSModel / Variable Notional Swaption: unresolved problem");

	return basketVol;
}

////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: TwoAssetsOptionPrice
///	Returns: 
///	Action : Computes the option price for two assets
////////////////////////////////////////////////////
double ARM_InfBSModel::TwoAssetsOptionPrice(
	double CPIForward,
	double secondAssetFwd, 
	double strike,
	int callput, 
	double discounting, 
	ARM_InfIdx* infIdx,
	ARM_IRIndex* secondIndex,
	ARM_Object* optionContext, 
	StoreInfoObj& storeInfo	)
{
	if( ARM_InfSwaptionContext* infSwaptionContext = dynamic_cast<ARM_InfSwaptionContext*>(optionContext) )
	{
		if( dynamic_cast<ARM_InfIdx*>(secondIndex) )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				ARM_USERNAME + ": inflation vs inflation can only be priced with a multi bs model!");
		
		double tenor	= infSwaptionContext->GetOptionMaturity();
			
		/// Compute the volatility for the inf leg coupon in the case of the OAT
		double coupon, timeToStart	= ARM_SparseVolCube::spotTime, volAsset1, volAsset2, strike1, strike2;

		/// Renormalize the spread moyness of the inf part accordingly to the Libor's one
		
		//strike1 = CPIForward+strike-secondAssetFwd;
		strike1 = infSwaptionContext->GetStrikeForVol1();

		if (infSwaptionContext->GetFirstAssetType() == K_OATTYPE_LEG)
		{
			coupon = infSwaptionContext->GetFirstAssetCoupon();

			ARM_InfVolComputation_Producer_Std* infVolCompProdStd = new ARM_InfVolComputation_Producer_Std( this );
			
			ARM_VolCurve* InfSwoptVolCurveWithCoupon = infVolCompProdStd->CompleteOATSwopVolCurve(
			itsAsOfDate,
			itsInfSwoptVolCurve,
			coupon);

			volAsset1 = ARM_InfBSModel::ComputeVolForSwaption(
			strike1, 
			infSwaptionContext->GetExpiryAsset1(),
			infSwaptionContext->GetTenorAsset1(),
			InfSwoptVolCurveWithCoupon,
			"inflation" );

		}
		else
		{
			volAsset1 = ARM_InfBSModel::ComputeVolForSwaption(
			strike1,
			infSwaptionContext->GetExpiryAsset1(),
			infSwaptionContext->GetTenorAsset1(),
			itsInfSwoptVolCurve,
			"inflation" );
		};

		//strike2 = -(secondAssetFwd-CPIForward+strike);
		strike2 = infSwaptionContext->GetStrikeForVol2()-secondAssetFwd;
			
		if (itsIRSwoptVolCube!=NULL)
		{
			volAsset2 = ARM_InfBSModel::ComputeVolForSwaptionWithSmile(
			strike2,
			infSwaptionContext->GetExpiryAsset2(),
			infSwaptionContext->GetTenorAsset2(),
			itsIRSwoptVolCube,
			"interest rates" );
		}
		else
		{
			volAsset2= ARM_InfBSModel::ComputeVolForSwaption(
			strike2,
			infSwaptionContext->GetExpiryAsset2(),
			infSwaptionContext->GetTenorAsset2(),
			itsIRSwoptVolCurve,
			"interest rates" );
		}

		double value(0.); 
		double correlMarket = GetInfIRCorrel( infSwaptionContext->GetTenorAsset1(), 
			infSwaptionContext->GetTenorAsset2(), infIdx, secondIndex, "INF/IR_SWOPT" );
		
		double correl = correlMarket;
		if (infSwaptionContext->HasBeenComputedTwice())
		{	
			volAsset1		= infSwaptionContext->GetVolAsset1();
			volAsset2		= infSwaptionContext->GetVolAsset2();
			correl			= infSwaptionContext->GetCorrelation();

			double optionType = ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION;
			/// Uses the closed formula library on the spread option valuation		
			int n			= 120;
			value			= Export_LogNormal_SpreadOption(CPIForward, secondAssetFwd, volAsset1, volAsset2, 
															correl, strike, tenor, callput, optionType, n);
		}

		else
		{		
			double optionType = ARM_CF_SpreadDigitalOption_Formula::SPREADOPTION;

			/// Uses the closed formula library on the spread option valuation		
			int n			= 120;
			value			= Export_LogNormal_SpreadOption(CPIForward, secondAssetFwd, volAsset1, volAsset2, 
															correl, strike, tenor, callput, optionType, n);
		}



		/// store the appropriate information
		double data[10];
		data[0]= infSwaptionContext->GetFwdAsset1();
		data[1]= volAsset1  * CC_NS( ARM_Constants, volBase );
		data[2]= infSwaptionContext->GetTenorAsset1();
		data[3]= infSwaptionContext->GetFwdAsset2();
		data[4]= volAsset2 * CC_NS( ARM_Constants, volBase );
		data[5]= infSwaptionContext->GetTenorAsset2();
		data[6]= correl * CC_NS( ARM_Constants, correlBase );
		data[7]= discounting;
		data[8]= infSwaptionContext->GetOptionMaturity();
		data[9]= strike * CC_NS( ARM_Constants, rateBase );
		storeInfo.Store( data );

		return ( discounting * value * CC_NS( ARM_Constants, rateBase ) );
	}
	else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			ARM_USERNAME + ": only swaption requires two asset pricing!");
}


////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: GetInfIRCorrel
///	Returns: 
///	Action : Computes inflation interest rates correlation between two dates!
////////////////////////////////////////////////////
double ARM_InfBSModel::GetInfIRCorrel( 
	double TjFromAsOf, 
	double TiFromAsOf,
	ARM_InfIdx* infIdx,
	ARM_IRIndex* otherIndex,
	const string& mktTag ) const
{
	/// when there is no correlmanager, return a correlation of zero!
	if( itsCorrelManager )
		return GetInfIRCorrelMatrix(infIdx,otherIndex,mktTag)->ComputeCorrelData( TjFromAsOf, TiFromAsOf)
			/ CC_NS(ARM_Constants,correlBase);
	else
		return 0.0;
}


////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: GetInfIRCorrelMatrix
///	Returns: 
///	Action : gets the inflation interest rates correlation
////////////////////////////////////////////////////
ARM_CorrelMatrix* ARM_InfBSModel::GetInfIRCorrelMatrix( ARM_InfIdx* infIdx, ARM_IRIndex* otherIndex, const string& mktTag  ) const
{
	string ccy			= otherIndex->GetCurrencyUnit()->GetCcyName();
	string indexName	= GetInfFwdCurv()->GetInfIdxName();
	string intraMktTag	= ccy < indexName ? ccy + "_" + indexName : indexName + "_"+ ccy;

	if( itsCorrelManager )
		return itsCorrelManager->ComputeCorrelData( mktTag,  intraMktTag );
	else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			ARM_USERNAME + ": could not find a correlmanager!");
}

////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: GetInfSwpCorrelMatrix
///	Returns: 
///	Action : gets the inflation swaps correlation
////////////////////////////////////////////////////
ARM_CorrelMatrix* ARM_InfBSModel::GetInfSwpCorrelMatrix( ARM_InfIdx* infIdx, const string& mktTag ) const
{
	string indexName	= GetInfFwdCurv()->GetInfIdxName();
	string intraMktTag	= indexName + "_"+ indexName;

	if( itsCorrelManager )
		return itsCorrelManager->ComputeCorrelData( mktTag,  intraMktTag );
	else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			ARM_USERNAME + ": could not find a correlmanager!");
}

////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: GetInfSwpCorrelMatrix
///	Returns: 
///	Action : gets the interest rates swaps correlation
////////////////////////////////////////////////////
ARM_CorrelMatrix* ARM_InfBSModel::GetIRSwpCorrelMatrix( ARM_IRIndex* otherIndex, const string& mktTag  ) const
{

	string ccy			= otherIndex->GetCurrencyUnit()->GetCcyName();
	string intraMktTag	= ccy + "_"+ ccy;

	if( itsCorrelManager )
		return itsCorrelManager->ComputeCorrelData( mktTag,  intraMktTag );
	else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			ARM_USERNAME + ": could not find a correlmanager!");

}
////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: GetInfInfCorrelMatrix
///	Returns: double
///	Action : Get the inflatin interest rate correlation
////////////////////////////////////////////////////
ARM_CorrelMatrix* ARM_InfBSModel::GetInfInfCorrelMatrix( 
	ARM_InfIdx* infIdx1, ARM_InfIdx* infIdx2, const string& mktTag ) const
{
	if( infIdx1->GetIndexName() != infIdx2->GetIndexName() )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": inflation index should be the same otherwise use a multibsmodel! " );

	string indexName	= infIdx1->GetIndexName();
	string intraMktTag	= indexName + "_" + indexName;

	if( GetCorrelManager() )
		return GetCorrelManager()->ComputeCorrelData( mktTag,  intraMktTag ) ;
	else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": correlmanager not found! " );
}



////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: GetInfInfCorrel
///	Returns: double
///	Action : Get the inflatin interest rate correlation
////////////////////////////////////////////////////
double ARM_InfBSModel::GetInfInfCorrel( double TjFromAsOf, double TiFromAsOf, 
	ARM_InfIdx* infIdx1, ARM_InfIdx* infIdx2, const string& mktTag) const
{
	return GetInfInfCorrelMatrix(infIdx1, infIdx2, mktTag )->ComputeCorrelData( TjFromAsOf, TiFromAsOf ) 
		/ CC_NS(ARM_Constants,correlBase);
}


////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: GetInfSwpCorrel
///	Returns: double
///	Action : Get the inflation swap correlation
////////////////////////////////////////////////////
double ARM_InfBSModel::GetInfSwpCorrel(
				double tenor1,
				double tenor2,
				ARM_InfIdx* infIdx,
				const string& mktTag) const
{
	
	return GetInfSwpCorrelMatrix(infIdx, mktTag )->ComputeCorrelData(tenor1, tenor2)
		/ CC_NS(ARM_Constants,correlBase);
}


////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: GetIRSwpCorrel
///	Returns: double
///	Action : Get the interest rate swap correlation
////////////////////////////////////////////////////
double ARM_InfBSModel::GetIRSwpCorrel(
				double tenor1,
				double tenor2,
				ARM_IRIndex* irIdx,
				const string& mktTag) const
{
	return GetIRSwpCorrelMatrix(irIdx, mktTag )->ComputeCorrelData( tenor1, tenor2)
		/ CC_NS(ARM_Constants,correlBase);
}


////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: GetVolatility
///	Returns: ARM_VolCurve*
///	Action : returns the correct volatility curve

////////////////////////////////////////////////////

ARM_VolCurve* ARM_InfBSModel::GetVolatility(int mode ) const
{
	switch( mode )
	{
	/// IR Swaption
	case K_SWOPT:
		{
			if( !itsIRSwoptVolCurve )
			{
				ARM_THROW( ERR_INVALID_ARGUMENT, "the inflation bs model has no interest rate swpation volatility data!" );
			}
			else
				return itsIRSwoptVolCurve;
		}
		break;
		
	/// standard case
	case K_INF_YOY: case -1:
		return ARM_BSModel::GetVolatility( mode );
		break;
	case K_INF_OAT:
		return ARM_BSModel::GetVolatility( mode );
		break;
	/// Cannot do anything else!
	case K_IRG: default:
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, "inflation bs model only support inflation and interest rate volatility data!");
		}
		break;
	}
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------*/
/*---- End Of File ----*/

