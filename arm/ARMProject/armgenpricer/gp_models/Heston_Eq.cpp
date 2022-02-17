/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Heston_Eq.cpp
 *
 *  \brief Heston Equity 
 *
 *	\author  E. Ezzine
 *	\version 1.0
 *	\date May 2006
 */

/// to remove identified warning
#include "gpbase/removeidentifiedwarning.h"

/// gpinfra
#include "gpinfra/pricingmodeltype.h"
#include "gpinfra/discretisationscheme.h"
#include "gpinfra/pricingstates.h"

#include "gpclosedforms/normal.h"
#include "gpclosedforms/heston_pricer.h"


/// gpmodel
#include "gpmodels/Heston_Eq.h"
/*
#include <iostream>
#include <fstream>
#include <iomanip>
*/



CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Eq
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_HestonModel_Eq::ARM_HestonModel_Eq(const ARM_ZeroCurvePtr& zc, 
									   ARM_ModelParamsHeston_Eq* modelParam)//ARM_Heston_ModelParams* modelParam)
:	ARM_EqFxBase(zc,modelParam)
{
	if(	!dynamic_cast<const ARM_ModelParamsHeston_Eq*>( modelParam ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": accepted model params are only Heston Equity model param" );
}


////////////////////////////////////////////////////
///	Class   : ARM_HestonModel_Eq
///	Routines: GetSettlementCalendar
///	Returns : string
///	Action  : get the settlement calendar
////////////////////////////////////////////////////
string ARM_HestonModel_Eq::GetSettlementCalendar(const string& modelName) const
{
#if defined(__GP_STRICT_VALIDATION)
	if( GetZeroCurve() == ARM_ZeroCurvePtr(NULL) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": no discount curve set!" );
	if( !GetZeroCurve()->GetCurrencyUnit() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": no ccy unit set on the discount curve set!" );
#endif

	return GetZeroCurve()->GetCurrencyUnit()->GetCcyName();
}


////////////////////////////////////////////////////
///	Class   : ARM_HestonModel_Eq
///	Routines: GetSettlementGap
///	Returns : double
///	Action  : get the settlement gap
////////////////////////////////////////////////////
double ARM_HestonModel_Eq::GetSettlementGap(const string& modelName) const
{
#if defined(__GP_STRICT_VALIDATION)
	if( GetZeroCurve() == ARM_ZeroCurvePtr(NULL) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": no discount curve set!" );
	if( !GetZeroCurve()->GetCurrencyUnit() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": no ccy unit set on the discount curve set!" );
#endif

	return GetZeroCurve()->GetCurrencyUnit()->GetSpotDays();
}


////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Eq
///	Routine: GetType
///	Returns: int
///	Action : tells the type of the model
////////////////////////////////////////////////////
int ARM_HestonModel_Eq::GetType() const
{
	return MT_EQUITY_MODEL;
}

////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Eq
///	Routine: TreeStatesToModelStates
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_HestonModel_Eq::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{   
	ARM_THROW( ERR_INVALID_ARGUMENT,"ARM_HestonModel_Eq,TreeStatesToModelStates :  not implemented ARM_HestonModel_Eq Model!");
}

////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Eq
///	Routine: VarianceToTime
///	Returns: a time
///	Action : Compute the time such that
///          var(t)=var
////////////////////////////////////////////////////
double ARM_HestonModel_Eq::VarianceToTime(double var,double minTime,double maxTime) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT,"ARM_HestonModel_Eq,VarianceToTime not implemented ARM_HestonModel_Eq Model!");

}

////////////////////////////////////////////////////
///	Class   : ARM_HestonModel_Eq
///	Routines: FirstPricingStates
///	Returns : void
///	Action  : FirstPricingStates
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_HestonModel_Eq::FirstPricingStates( size_t bucketSize ) const
{
	const size_t nbPayoffs			= 0;
	size_t factorsNb				= FactorCount();
	size_t nbModelStates			= ModelStatesSize();
	
	const ARM_ModelParams* params		= GetModelParams();
	const ARM_ModelParams_Eq* eqFunc	= dynamic_cast<const ARM_ModelParams_Eq*>( params );
	if(  !eqFunc )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": eqFunc is a null pointor" );
	double spot						= eqFunc->GetSpot();
	double initvol					= ((ARM_Heston_ModelParams*) GetModelParams())->InitialVol();
	double scaling					= ((ARM_Heston_ModelParams*) GetModelParams())->Scaling(0);

	ARM_PricingStatesPtr initStates = ARM_PricingStatesPtr( new ARM_PricingStates(bucketSize,nbModelStates,nbPayoffs,factorsNb) );

	//Each model state is initialized to 0.0
	for(size_t i=0; i<bucketSize; ++i )
	{
		initStates->SetModelState(i,0,spot);
		initStates->SetModelState(i,1,initvol*initvol*scaling*scaling);
	}
	return initStates;
}

////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Eq
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Default initialisation of the model and the
///          associated numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_HestonModel_Eq::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
	 /// test the numeraire and its type!
	ARM_NumerairePtr numeraire=GetNumeraire();
    if( numeraire == ARM_NumerairePtr(NULL) )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numeraire not set in the Equity model!");

	//// Initialise the numeraire
	numeraire->Init(*(GetDiscountFunctor()),payModelName,timeInfos);
	size_t modelNb = GetModelNb();

    int nbEvents=timeInfos.size();
    bool isSpotUse = nbEvents == 0 || (nbEvents==1 && timeInfos[0]->GetEventTime() <= K_NEW_DOUBLE_TOL);
    if(!isSpotUse)
    {
        // Numerical method and numeraire are needed to go on
        ARM_NumMethodPtr numMethod=GetNumMethod();
		if( numMethod == ARM_NumMethodPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numerical method not set in FxEQuity model!");

		/// creates the model schedule (smart pointor for exception safety!)
		ARM_DiscretisationScheme& discretisationScheme = ARM_EventTime();
		CC_NS(std,auto_ptr)<std::vector<double>> ptimeSteps( discretisationScheme.ModelTimesFromTimesInfo(timeInfos, *this) );
    
		/// Set the basic schedule in the numerical method and...
		numMethod->SetTimeSteps(*ptimeSteps);

		double firstInductTime = timeInfos[0]->GetEventTime();

		/// ...initialise it
		return numMethod->Init(*this,firstInductTime);
    }
    else
    {
        // Compute a single model states set to (0.0,...,0.0)
        int nbDir = FactorCount();
        ARM_PricingStatesPtr initStates(new ARM_PricingStates(1,nbDir,0));
        for(int i=0;i<nbDir;++i)
            initStates->SetModelState(0,i,0.0);
		
        return initStates;
    }
}

////////////////////////////////////////////////////
///	Class   : ARM_HestonModel_Eq
///	Routines: NumMethodStateLocalVariances
///	Returns : void
///	Action  : NumMethodStateLocalVariances
////////////////////////////////////////////////////
void ARM_HestonModel_Eq::NumMethodStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const
{
	int factorsNb		= FactorCount();
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
	
    localVariances.resize((timeStepsSize-1)*(modelNb+1));
	for(int i = 0; i < timeStepsSize-1 ; ++i)
	{
		/// initialize everything
		localVariances[i] = new ARM_GP_Matrix(factorsNb,factorsNb,0.);

		for (int k = 0; k < factorsNb; k++)
			(*localVariances[i])(k,k) = 1.;
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_HestonModel_Eq
///	Routines: NumMethodStateGlobalVariances
///	Returns : void
///	Action  : NumMethodStateGlobalVariances
////////////////////////////////////////////////////
void ARM_HestonModel_Eq::NumMethodStateGlobalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& globalVariances ) const
{
	// useless?
	int factorsNb		= FactorCount();
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
	int	offsetIndex		= (timeStepsSize-1)*modelNb;
    globalVariances.resize(timeStepsSize*(modelNb+1));
	for (size_t i=0;i<timeStepsSize;i++)
		globalVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(factorsNb,1.0);
}

////////////////////////////////////////////////////
///	Class   : ARM_HestonModel_Eq
///	Routines: NumMethodStateLocalGlobalVariances
///	Returns : void
///	Action  : NumMethodStateLocalGlobalVariances
////////////////////////////////////////////////////
void ARM_HestonModel_Eq::NumMethodStateLocalGlobalVariances(
        const std::vector<double>& timeSteps,
        ARM_MatrixVector& localVariances,
        ARM_MatrixVector& variances) const
{
	NumMethodStateLocalVariances (	timeSteps, localVariances	);
	NumMethodStateGlobalVariances(	timeSteps, variances		);
}

////////////////////////////////////////////////////
///	Class   : ARM_HestonModel_Eq
///	Routines: ModelStateLocalVariancesAndStdDev
///	Returns : void
///	Action  : computes the local variance and std dev 
////////////////////////////////////////////////////
void ARM_HestonModel_Eq::ModelStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps)
{
	ARM_PricingModel::ModelStateLocalVariancesAndStdDev(timeSteps);	
}

////////////////////////////////////////////////////
///	Class   : ARM_HestonModel_Eq
///	Routines: NumMethodStateLocalVariancesAndStdDev
///	Returns : void
///	Action  : computes the local variance and std dev 
////////////////////////////////////////////////////
void ARM_HestonModel_Eq::NumMethodStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps )
{
	ARM_MatrixVector numMethodStateLocalVars;
	
	/// computes the local variance
	NumMethodStateLocalVariances( timeSteps, numMethodStateLocalVars);
	
	/// set the result
	SetNumMethodStateLocalVars(numMethodStateLocalVars);
	
}

////////////////////////////////////////////////////
///	Class   : ARM_HestonModel_Eq
///	Routines: ModelStateLocalVariances
///	Returns : void
///	Action  : ModelStateLocalVariancess
////////////////////////////////////////////////////
void ARM_HestonModel_Eq::ModelStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const
{
	int factorsNb		= FactorCount();
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
	
    localVariances.resize((timeStepsSize-1)*(modelNb+1));
	for(int i = 0; i < timeStepsSize-1 ; ++i)
	{
		/// initialize everything
		localVariances[i] = new ARM_GP_Matrix(factorsNb,factorsNb,0.);

		for (int k = 0; k < factorsNb; k++)
			(*localVariances[i])(k,k) = 1.;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Eq
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
void ARM_HestonModel_Eq::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
	double t			= GetNumMethod()->GetTimeStep(timeIndex);
	double T			= GetNumMethod()->GetTimeStep(timeIndex+1);
	size_t factorsNb	= FactorCount();
	size_t statesNb		= states->size();
	size_t modelNb		= GetModelNb();

	const ARM_ModelParams* params		= GetModelParams();
	const ARM_ModelParams_Eq* eqFunc	= dynamic_cast<const ARM_ModelParams_Eq*>( params );
	if(  !eqFunc )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": eqFunc is a null pointor" );

	double spot		= eqFunc->GetSpot();
	
	double q		= ((ARM_Heston_ModelParams*) GetModelParams())->MC_scheme();
	
	double rho		= ((ARM_Heston_ModelParams*) GetModelParams())->Correlation(t);
	double eps		= ((ARM_Heston_ModelParams*) GetModelParams())->VolOfVol(t);
	double lambda	= ((ARM_Heston_ModelParams*) GetModelParams())->VolMeanReversion(t);
	double ltvol	= ((ARM_Heston_ModelParams*) GetModelParams())->LongTermVol(t);
	double theta	= ltvol*ltvol;

	double beta		= ((ARM_Heston_ModelParams*) GetModelParams())->Beta(t);
	double scale	= ((ARM_Heston_ModelParams*) GetModelParams())->Scaling(t);
	theta *= scale*scale;
	eps *= scale;

	int EQ = modelNb;
	int VAR = 1+modelNb;

	double oldEq,oldVar,newEq,newVar,oldX,newX;

	double dt = (T-t)/K_YEAR_LEN;
	double sqrt_dt = sqrt(dt);
	double coeff0 = exp(-lambda*dt);
	double coeff3 = fabs(lambda)<K_NEW_DOUBLE_TOL?dt:(1-coeff0)/lambda;
	double coeff1 = eps*eps*coeff0*coeff3;
	double coeff2 = 0.5*theta*eps*eps*lambda*coeff3*coeff3;

	double dW,dB,E_var,V_var,S_var;
	double G1,G2;

	if (fabs(q)<K_NEW_DOUBLE_TOL)
	{
		for( size_t i=0;i<statesNb; ++i )
		{
			oldVar	= states->GetModelState(i,VAR);
			oldEq	= states->GetModelState(i,EQ);

			oldX	= oldEq + (1.-beta)/beta*spot;
				
			dB = sqrt_dt*states->GetNumMethodState(i,VAR);
			dW = sqrt_dt*states->GetNumMethodState(i,EQ);
				
			E_var = theta+(oldVar-theta)*coeff0;
			V_var = oldVar*coeff1+coeff2;
			S_var = sqrt(log(1.+V_var/E_var/E_var)/dt);

			newVar	= E_var * exp( -0.5*S_var*S_var*dt + S_var*dB );
			newX	= oldX * exp( -0.5*oldVar*dt + sqrt(oldVar)*(rho*dB+sqrt(1.-rho*rho)*dW ) );
			newEq	= newX - (1.-beta)/beta*spot;

			states->SetModelState(i,VAR,newVar);
			states->SetModelState(i,EQ,newEq);
		}
	}
	else if (fabs(q-1.)<K_NEW_DOUBLE_TOL)
	{
		for( size_t i=0;i<statesNb; ++i )
		{
			oldVar	= states->GetModelState(i,VAR);
			oldEq	= states->GetModelState(i,EQ);
			oldX	= oldEq + (1.-beta)/beta*spot;
				
			G1 = states->GetNumMethodState(i,VAR);
			G2 = states->GetNumMethodState(i,EQ);
				
			E_var = theta+(oldVar-theta)*coeff0;
			V_var = oldVar*coeff1+coeff2;
			
			double psi=V_var/E_var/E_var;
			if(psi<1.5)
			{
				double tipmo=(2./psi)-1.;
				double b=sqrt(tipmo+sqrt(tipmo+1.)*sqrt(tipmo));
				double a=E_var/(1. +b*b);
				newVar = a*(b + G1)*(b+ G1);
			}
			else
			{
				double invbeta=0.5*((psi+1)*E_var);
				double p=(psi-1.)/(psi+1.);
				double u=NormalCDF(G1);
				newVar=u<p?0.0:invbeta*(log((1.-p)/(1.-u)));
			}
			double K0=-rho*lambda*dt*theta/eps;
			double K1=0.5*dt*((rho*lambda)/(eps) -0.5)-(rho/eps);
			double K2=K1+2*(rho/eps);
			double K3=0.5*dt*(1.-rho*rho);
			double K4=K3;
			double tmp=K0+K1*oldVar+K2*newVar+sqrt(K3*(oldVar+newVar))*G2;

			newX	= oldX*exp(tmp);
			newEq	= newX - (1.-beta)/beta*spot;

			states->SetModelState(i,VAR,newVar);
			states->SetModelState(i,EQ,newEq);
		}
	}
	else
	{
		double flag;
		if (fabs(q-2.)<K_NEW_DOUBLE_TOL)
			flag=0.;
		else 
			flag=1.;

		for( size_t i=0;i<statesNb; ++i )
		{
			oldVar	= states->GetModelState(i,VAR);
			oldEq	= states->GetModelState(i,EQ);

			oldX	= oldEq + (1.-beta)/beta*spot;
				
			dB = sqrt_dt*states->GetNumMethodState(i,VAR);
			dW = sqrt_dt*states->GetNumMethodState(i,EQ);
				
			newVar	= oldVar + lambda*(theta-oldVar)*dt + eps*sqrt(oldVar)*dB;
			if (newVar<0.)
				newVar = -flag*newVar;

			newX	= oldX * exp( -0.5*oldVar*dt + sqrt(oldVar)*(rho*dB+sqrt(1.-rho*rho)*dW ) );
			newEq	= newX - (1.-beta)/beta*spot;

			states->SetModelState(i,VAR,newVar);
			states->SetModelState(i,EQ,newEq);
		}
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Eq
///	Routine: ComputeFwdAtTime
///	Returns: double
///	Action : 
////////////////////////////////////////////////////
double ARM_HestonModel_Eq::ComputeFwdAtTime( double evalTime  ) const
{
	const ARM_ModelParams* params = GetModelParams();
	const ARM_ModelParams_Eq* eqFunc = dynamic_cast<const ARM_ModelParams_Eq*>( params );

	if(  !eqFunc )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": eqFunc is a null pointor" );

	return eqFunc->Forward( evalTime );
}


////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Eq
///	Routine: Forward
///	Returns: a vector of forward (t,T)
///	Action : computes the forward of the Heston mode
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HestonModel_Eq::Forward(
		const string& modelName, 
        double evalTime,
	    double expiryTime,
	    double settlementTime,
	    double payTime,
        const ARM_PricingStatesPtr& states) const
{
	if( evalTime<=K_NEW_DOUBLE_TOL || states == ARM_PricingStatesPtr(NULL) )
	{
		size_t stateSize = states->size();
		double forwardValue = ComputeFwdAtTime(settlementTime);
		return ARM_VectorPtr( new std::vector<double>(stateSize,forwardValue) );
	}
	else
	{
		size_t modelNb		= GetModelNb();
		size_t stateSize	= states->size();
		int EQ				= modelNb;
		int VAR				= 1+modelNb;
		double forwardEnd	= ComputeFwdAtTime( settlementTime );
		double forwardStart = ComputeFwdAtTime( evalTime );
		double spot;
		
		ARM_GP_VectorPtr fwdVector = ARM_GP_VectorPtr( new std::vector<double>(stateSize, 0.0 ));
		for( size_t i=0; i<stateSize; ++i )
		{
			spot = states->GetModelState(i,EQ);
			(*fwdVector )[i] = spot* forwardEnd/forwardStart;
		}
		return fwdVector;
	}
}

ARM_VectorPtr ARM_HestonModel_Eq::DiscountFactor( 
		const string& curveName,
        double evalTime, 
		const std::vector<double>&  maturityTimes,
        const ARM_PricingStatesPtr& states) const
{
	size_t i, 
		nbStates = ((states!=ARM_PricingStatesPtr(NULL))? states->size(): 1);

    /// Waiting for the access to the yield curve with curveName
    ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();

	double zct = ZcCurve->DiscountPrice(evalTime/K_YEAR_LEN),
		zcT = 1;
	
	ARM_VectorPtr result(new std::vector<double>(nbStates));
	
	for(i = 0; i < nbStates; i++)
	{
		zcT = ZcCurve->DiscountPrice(maturityTimes[i]/K_YEAR_LEN);
		(*result)[i] = zcT/zct;
	}

	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Eq
///	Routine: CallVectorial
///	Returns: a vector of call vectorial
///	Action : computes equity or fx call
////////////////////////////////////////////////////
//NDC
ARM_VectorPtr ARM_HestonModel_Eq::CallVectorial(
	const string& modelName,
    double evalTime,
	const std::vector<double>& expiryTimes,
	const std::vector<double>& settlementTimes,
	const std::vector<double>& strikePerState,
	int callPut,
	const std::vector<double>& payTimes,
    const ARM_PricingStatesPtr& states,
    ARM_PricingContext* context) const
{
	size_t i, size = states->size();
	ARM_VectorPtr result(new std::vector<double>(size,0.));

	if ( 0 != size) {
		/*
		std::string fileName = "c:/Debug/StatAverageVector_operator";
#if defined( _DEBUG )
		fileName += ".dbg";
#else
		fileName += ".release";
#endif
		fileName += ".txt";

		std::ofstream os;
		os.open(fileName.c_str(),ios::out | ios::app );

		if (!os) {
			std::string msg = "Failed to open file : " + fileName;
			os.close();
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg.c_str() );
		}

		os << std::endl << "ARM_HestonModel_Eq::CallVectorial(...)" << std::endl;
		*/

		//os << "states size = " << size << std::endl;

		// Check sizes
		size_t expiryTimesSize	= expiryTimes.size(),
			strikePerStateSize	= strikePerState.size(),
			settlementTimesSize = settlementTimes.size(),
			payTimesSize		= payTimes.size();


		if ( expiryTimesSize > 1 ) {
			if ( expiryTimesSize != size ) {
				CC_Ostringstream oss;
				oss << "expiryTimes size (" << expiryTimesSize << ") != states size (" << size << ")\n"; 
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + oss.str() );
			}
		}
		else if ( expiryTimesSize == 0 ) {
			CC_Ostringstream oss;
			oss << "expiryTimes size must be greater than zero)\n"; 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + oss.str() );
		}
		if ( strikePerStateSize > 1 ) {
			if ( strikePerStateSize != size ) {
				CC_Ostringstream oss;
				oss << "strikePerState size (" << strikePerStateSize << ") != states size (" << size << ")\n"; 
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + oss.str() );
			}
		}
		else if ( strikePerStateSize == 0 ) {
				CC_Ostringstream oss;
				oss << "strikePerSate size must be greater than zero)\n"; 
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + oss.str() );
		}
		if ( settlementTimesSize > 1 ) {
			if ( settlementTimesSize != size ) {
				CC_Ostringstream oss;
				oss << "settlementTimes size (" << settlementTimesSize << ") != states size (" << size << ")\n"; 
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + oss.str() );
			}
		}
		else if ( settlementTimesSize == 0 ) {
			CC_Ostringstream oss;
			oss << "settlementTimes size must be greater than zero)\n"; 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + oss.str() );
		}
		if ( payTimesSize > 1 ) {
			if ( payTimesSize != size ) {
				CC_Ostringstream oss;
				oss << "payTimes size (" << payTimesSize << ") != states size (" << size << ")\n"; 
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + oss.str() );
			}
		}
		else if ( payTimesSize == 0 ) {
			CC_Ostringstream oss;
			oss << "payTimes size must be greater than zero)\n"; 
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + oss.str() );
		}


		size_t modelNb		= GetModelNb();
		int EQ				= modelNb;
		int VAR				= 1+modelNb;
		
		double opt, K;
		double spot, var, expiryTime, tau;

		for(i = 0; i < size; i++)
		{

			expiryTime		= ((expiryTimesSize == 1) ? expiryTimes[0] : expiryTimes[i]);

			K				= ((strikePerStateSize == 1) ? strikePerState[0] : strikePerState[i]);
			spot			= states->GetModelState(i,EQ);
			var				= states->GetModelState(i, VAR);

			tau				= (expiryTime - evalTime) / 365.;

			const ARM_Heston_ModelParams* const params = dynamic_cast<const ARM_Heston_ModelParams * const >(this->GetModelParams());

			double volMeanReversion = params->VolMeanReversion(expiryTime);
			double longTermVol		= params->LongTermVol(expiryTime);
			double correlation		= params->Correlation(expiryTime);
			double volOfVol			= params->VolOfVol(expiryTime);
			double beta				= params->Beta(expiryTime);
			double scaling			= params->Scaling(expiryTime);

			double longTermVariance = longTermVol*longTermVol;

			var						/= scaling*scaling;

			/*
			os << i << "expiryTime = " << expiryTime << std::endl
					<< "K = " << K << std::endl
					<< "spot = " << spot << std::endl
					<< "var = " << var << std::endl
					<< "tau = " << tau << std::endl
					<< "volMeanReversion = " << volMeanReversion << std::endl
					<< "longTermVol = " << longTermVol << std::endl
					<< "correlation = " << correlation << std::endl
					<< "volOfVol = " << volOfVol << std::endl
					<< "beta = " << beta << std::endl
					<< "scaling = " << scaling << std::endl
					<< "longTermVariance = " << longTermVariance << std::endl;
			*/

			ARM_HestonOptionPricer pricer(
										tau, spot, K, callPut, var,
										volMeanReversion,
										longTermVariance,
										correlation,
										volOfVol,
										beta,
										scaling
										 );

			opt = pricer.price();
			opt = (opt < 0.) ? callPut * (spot - K) > 0. ? callPut * (spot - K) : 0. : opt;

			//os << "Heston price = " << opt << std::endl;

			(*result)[i] = opt;

		}

		//os.close();
	}

	return result;
}



////////////////////////////////////////////////////
///	Class  : ARM_HestonModel_Eq
///	Routine: CallVectorial
///	Returns: a vector of call vectorial
///	Action : computes equity or fx call
////////////////////////////////////////////////////
//NDC
ARM_VectorPtr ARM_HestonModel_Eq::GreekVectorial(
									const string& modelName,
									double evalTime,
									const std::vector<double>& expiryTimes,
									const std::vector<double>& settlementTimes,
									const std::vector<double>& strikePerState,
									int callPut,
									const std::vector<double>& payTimes,
									const string& greekType,
									const ARM_PricingStatesPtr& states,
									ARM_PricingContext* context) const
{
	if ( greekType != "Delta")
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "Only DELTA is available for greeks calculation" );

	size_t i, size = states->size();
	ARM_VectorPtr result(new std::vector<double>(size));

	if ( 0 != size) {
		// Check sizes
		size_t expiryTimesSize	= expiryTimes.size(),
			strikePerStateSize	= strikePerState.size(),
			settlementTimesSize = settlementTimes.size(),
			payTimesSize		= payTimes.size();


		if ( expiryTimesSize > 1 )
			if ( expiryTimesSize != size ) {
				CC_Ostringstream oss;
				oss << "expiryTimes size (" << expiryTimesSize << ") != states size (" << size << ")\n"; 
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + oss.str() );
			}
		if ( strikePerStateSize > 1 )
			if ( strikePerStateSize != size ) {
				CC_Ostringstream oss;
				oss << "strikePerState size (" << strikePerStateSize << ") != states size (" << size << ")\n"; 
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + oss.str() );
			}
		if ( settlementTimesSize > 1 )
			if ( settlementTimesSize != size ) {
				CC_Ostringstream oss;
				oss << "settlementTimes size (" << settlementTimesSize << ") != states size (" << size << ")\n"; 
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + oss.str() );
			}
		if ( payTimesSize > 1 )
			if ( payTimesSize != size ) {
				CC_Ostringstream oss;
				oss << "payTimes size (" << payTimesSize << ") != states size (" << size << ")\n"; 
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + oss.str() );
			}


		size_t modelNb		= GetModelNb();
		int EQ				= modelNb;
		int VAR				= 1+modelNb;
		
		double opt, K, delta;
		double spot, var, expiryTime, tau;
		double deltaShift = 0.0001; // 1 bp

		for(i = 0; i < size; i++)
		{
			expiryTime		= ((expiryTimesSize == 1) ? expiryTimes[0] : expiryTimes[i]);

			K				= ((strikePerStateSize == 1) ? strikePerState[0] : strikePerState[i]);
			spot			= states->GetModelState(i,EQ);
			var				= states->GetModelState(i, VAR);

			tau				= (expiryTime - evalTime) / 365.;

			const ARM_Heston_ModelParams* const params = dynamic_cast<const ARM_Heston_ModelParams * const >(this->GetModelParams());

			double volMeanReversion = params->VolMeanReversion(expiryTime);
			double longTermVol		= params->LongTermVol(expiryTime);
			double correlation		= params->Correlation(expiryTime);
			double volOfVol			= params->VolOfVol(expiryTime);
			double beta				= params->Beta(expiryTime);
			double scaling			= params->Scaling(expiryTime);

			double longTermVariance = longTermVol*longTermVol;

			var						/= scaling*scaling;

			spot += 0.5*deltaShift;

			ARM_HestonOptionPricer pricer(
										tau, spot, K, callPut, var,
										volMeanReversion,
										longTermVariance,
										correlation,
										volOfVol,
										beta,
										scaling
										 );

			opt = pricer.price();

			delta = opt < 0. ? callPut * (spot - K) > 0. ? callPut * (spot - K) : 0. : opt;

			spot -= deltaShift;

			ARM_HestonOptionPricer pricer2(
										tau, spot, K, callPut, var,
										volMeanReversion,
										longTermVariance,
										correlation,
										volOfVol,
										beta,
										scaling
										 );

			opt = pricer2.price();

			opt = opt < 0. ? callPut * (spot - K) > 0. ? callPut * (spot - K) : 0. : opt;

			delta = (delta - opt)/deltaShift;

			(*result)[i] = delta;

		}
	}

	return result;
}



////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHeston_Eq
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsHeston_Eq::ARM_ModelParamsHeston_Eq( const ARM_ModelParamVector& params, ARM_ZeroCurvePtr domCurve, double spot  )
:	ARM_ModelParams_Eq(domCurve, spot),
	ARM_Heston_ModelParams(params)
{
	if( params.size() <7 )
		ARM_THROW( ERR_INVALID_ARGUMENT, " expected at least 7 model parameters!" );

	if( !DoesModelParamExist(ARM_ModelParamType::QParameter) )
	{
		std::vector<double> abs(1,0.0);
		std::vector<double> ord(1,0.0);

		ARM_CurveModelParam modeParam(ARM_ModelParamType::QParameter, &ord, &abs);

		SetModelParam(&modeParam,ARM_ModelParamType::QParameter);
	}

	ValidateModelParams();
}

////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHeston_Eq
///	Routines: Validate
///	Returns :
///	Action  : validate the model params to check that this is compatible with the Heston model
////////////////////////////////////////////////////
void ARM_ModelParamsHeston_Eq::ValidateModelParams() const
{
	if( !DoesModelParamExist(ARM_ModelParamType::Volatility) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": HestonModelParam: requires volatility parameter!");

	if( !DoesModelParamExist(ARM_ModelParamType::InitialVol) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": HestonModelParam: requires initial vol parameter!");

	if( !DoesModelParamExist(ARM_ModelParamType::LongTermVol) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": HestonModelParam: requires a long term vol!");

	if( !DoesModelParamExist(ARM_ModelParamType::VolOfVol) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": HestonModelParam: requires a vol of vol!");

	if( !DoesModelParamExist(ARM_ModelParamType::VolMeanReversion) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": HestonModelParam: requires a vol MR!");

	if( !DoesModelParamExist(ARM_ModelParamType::Correlation) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": HestonModelParam: requires a correlation!");

	if( !DoesModelParamExist(ARM_ModelParamType::Beta) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": HestonModelParam: requires beta parameter!");
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHeston_Eq
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsHeston_Eq::ARM_ModelParamsHeston_Eq( const ARM_ModelParamsHeston_Eq& rhs )
:	ARM_ModelParams_Eq(rhs),
	ARM_Heston_ModelParams(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHeston_Eq
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsHeston_Eq::~ARM_ModelParamsHeston_Eq()
{}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHeston_Eq
///	Routine: toString
///	Returns: string
///	Action : Display the contents
////////////////////////////////////////////////////
string ARM_ModelParamsHeston_Eq::toString(const string& indent,const string& nextIndent) const
{
	string str;

	str += ARM_ModelParams_Eq::toString(indent,nextIndent);
	str += ARM_Heston_ModelParams::toString(indent,nextIndent);

	return str;
}
CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

