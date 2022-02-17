/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
*/

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

/// gpbase
#include "gpbase/comparisonfunctor.h"
#include "gpbase/cloneutilityfunc.h"

/// gpcalib
#include "gpcalib/vanillasecuritydensity.h"

/// gpmodels
#include "gpmodels/ProcessBuilderSmiledFRM.h"
#include "gpmodels/ModelParamsSmiledFRM.h"

/// gpclosedforms
#include "gpclosedforms/normal.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/vanilla_normal.h"
#include "gpclosedforms/gaussian_integrals.h"

/// gpnumlib
#include "gpnumlib/solver.h"
#include "gpnumlib/newtonraphson.h"
#include "gpnumlib/numfunction.h"
#include "gpnumlib/dichotomy.h"
#include "gpnumlib/brent.h"


//#include <crv/zerocurv.h>
/// standard libraries
#include <cmath>


const double DEFAULT_PRECISION				= 1.e-6;

#define GLPOINTSNB		16
#define INVSQRT2PI		0.398942280401433


CC_BEGIN_NAMESPACE( ARM )



////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilder
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ProcessBuilder::ARM_ProcessBuilder(double mrs,const ARM_Curve& volCurve)
:	itsMeanReversion(mrs),
	itsVolCurve(volCurve),
	itsVarCurve(ARM_Curve()),
	itsHorizonDate(0),
	itsFwd(0),
	itsNumFwd(0),
	itsDelta(0)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilder
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ProcessBuilder::ARM_ProcessBuilder(const ARM_ProcessBuilder& rhs)
:	itsMeanReversion(rhs.itsMeanReversion),
	itsVolCurve(rhs.itsVolCurve),
	itsVarCurve(rhs.itsVarCurve),
	itsHorizonDate(rhs.itsHorizonDate),
	itsFwd(rhs.itsFwd),
	itsNumFwd(rhs.itsNumFwd),
	itsDelta(rhs.itsDelta)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilder
///	Routine: PrecomputeIntegratedVariance
///	Returns: 
///	Action : Computes variance of underlying gaussian factor at breakpoint dates
////////////////////////////////////////////////////

void ARM_ProcessBuilder::ComputeIntegratedVariance()
{
	std::vector<double> timeSteps(itsVolCurve.GetAbscisses());
	std::vector<double> volFunc(itsVolCurve.GetOrdinates());
	size_t nbSteps = timeSteps.size();
	std::vector<double> varFunc(nbSteps);
	double prec,time,dt;
	if (fabs(MeanReversion()) < K_NEW_DOUBLE_TOL ){
		varFunc[0]=0.;
		prec=varFunc[0];
		time=timeSteps[0];
		for (size_t index=0;index<nbSteps;index++){
			dt=(timeSteps[index]-time)/K_YEAR_LEN;
			varFunc[index]=prec+volFunc[index]*volFunc[index]*dt;
			prec=varFunc[index];
			time=timeSteps[index];
		}
	}else{
		varFunc[0]=0.;
		prec=varFunc[0];
		time=timeSteps[0];
		for (size_t index=0;index<nbSteps;index++){
			dt=(timeSteps[index]-time)/K_YEAR_LEN;
			varFunc[index]=exp(-2*MeanReversion()*dt)*prec+volFunc[index]*volFunc[index]*(1.-exp(-2.*MeanReversion()*dt))/(2.*MeanReversion());
			prec=varFunc[index];
			time=timeSteps[index];
		}
	}
	itsVarCurve.SetAbscisses(timeSteps);
	itsVarCurve.SetOrdinates(varFunc);

}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilder
///	Routine: IntegratedVarianceFromZero
///	Returns: 
///	Action : Variance of underlying gaussian factor
////////////////////////////////////////////////////

double ARM_ProcessBuilder::IntegratedVarianceFromZero(double t) const{
	if (t<K_NEW_DOUBLE_TOL)
		return  0.0;

	double prec,dt,vol;
	
	std::vector<double> timeSteps(itsVolCurve.GetAbscisses());
	int index = lower_boundPosWithPrecision(timeSteps,t,K_FRM_TOL)-1;
	
	if(index == -1)
	{
	    prec=0.;
		dt=t/K_YEAR_LEN;
		vol=(itsVolCurve.GetOrdinates())[0];
	}
	else
	{
        prec=(itsVarCurve.GetOrdinates())[index];
		dt=(t-timeSteps[index])/K_YEAR_LEN;
		vol=Volatility(t);

	}

	if (fabs(MeanReversion()) < K_NEW_DOUBLE_TOL ){
		return prec+vol*vol*dt;
	}
	return exp(-2*MeanReversion()*dt)*prec+vol*vol*(1.-exp(-2.*MeanReversion()*dt))/(2.*MeanReversion());
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilder
///	Routine: Calibrate
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_ProcessBuilder::Calibrate(	double horizonDate,
												const ARM_VanillaSecurityDensity& terminalDensity,
												const std::vector<double>& evalDates,
												const std::vector<double>& numDates,
												bool doPDE)
{
		itsHorizonDate			= horizonDate;
		ComputeIntegratedVariance();

		std::vector<double> delta = terminalDensity.getInterestTerms();

		itsFwd					= terminalDensity.ComputeForwardRate();//(fwd-1.)/delta[0];
		itsNumFwd				= itsFwd;
		if (delta.size()==1)
			itsDelta			= delta[0];
		
}

////////////////////////////////////////////////////
///	Class   : ARM_ProcessBuilder
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_ProcessBuilder::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << indent << "( ResetDate : \t" << CC_NS(std,setprecision)(0) <<itsHorizonDate << " ) \t"<<CC_NS(std,setprecision)(6);
    os << indent << "ThFwd : \t" << CC_NS(std,fixed) << itsFwd << " \t";
	os << indent << "NumFwd : \t" << CC_NS(std,fixed) << itsNumFwd << " \t";
	os << indent << "ShiftEq : \t" << CC_NS(std,fixed) << GetEqShift() << " \t";
	os << indent << "VolEq : \t" << CC_NS(std,fixed) << GetEqVol() << " \n";
	return os.str();
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ProcessBuilderPDE::ARM_ProcessBuilderPDE(double mrs,const ARM_Curve& volCurve, bool withRescalling)
:	ARM_ProcessBuilder(mrs,volCurve),
	itsWithRescalling(withRescalling),
	itsDates(0),
	itsGlobalVariance(0),
	itsGridSize(0),
	itsNbStdDev(0),
	itsTheta(0),
	itsGrid(0),
	itsL1(0),
	itsL2(0),
	itsL3(0),
	itsR1(0),
	itsR2(0),
	itsR3(0),
	itsR(0),
	itsGam(0),
	itsXG(0),
	itsDFRatio(0),
	itsDxLogDFRatio(0),
	itsRate(0),
	itsDxRate(0),
	itsEqShift(0),
	itsEqVol(0),
	itsEqScaling(1),
	itsDiffusionTerm(0)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ProcessBuilderPDE::ARM_ProcessBuilderPDE(const ARM_ProcessBuilderPDE& rhs)
:	ARM_ProcessBuilder(rhs),
	itsWithRescalling(rhs.itsWithRescalling),
	itsDates(rhs.itsDates),
	itsGlobalVariance(rhs.itsGlobalVariance),
	itsGridSize(rhs.itsGridSize),
	itsNbStdDev(rhs.itsNbStdDev),
	itsTheta(rhs.itsTheta),
	itsGrid(rhs.itsGrid),
	itsL1(0),
	itsL2(0),
	itsL3(0),
	itsR1(0),
	itsR2(0),
	itsR3(0),
	itsR(0),
	itsGam(0),
	itsXG(0),
	itsEqShift(rhs.itsEqShift),
	itsEqVol(rhs.itsEqVol),
	itsEqScaling(rhs.itsEqScaling),
	itsDiffusionTerm(rhs.itsDiffusionTerm)
{
	DuplicateCloneablePtrVectorInPlace<std::vector<double>> (rhs.itsDFRatio, itsDFRatio);
	DuplicateCloneablePtrVectorInPlace<std::vector<double>> (rhs.itsDxLogDFRatio, itsDxLogDFRatio);
	DuplicateCloneablePtrVectorInPlace<std::vector<double>> (rhs.itsRate, itsRate);
	DuplicateCloneablePtrVectorInPlace<std::vector<double>> (rhs.itsDxRate, itsDxRate);
}
	
////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: GetProxyShift
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
double ARM_ProcessBuilderPDE::GetProxyShift(double strike1, double strike2,const ARM_VanillaSecurityDensity& density)
{
	bool success (true);
	double tex			= Horizon()/K_YEAR_LEN;
	double mkt1			= density.getDensityFunctor()->Call_Option(strike1, itsFwd, tex);
	double mkt2			= density.getDensityFunctor()->Call_Option(strike2, itsFwd, tex);
	
	SLN_Approx func(itsFwd,tex,strike1,mkt1,strike2,mkt2,density);

	UnaryFuncWithNumDerivative<double> funcWithDev(func);

	T_DichotomySolver< UnaryFuncWithNumDerivative<double> > solver(funcWithDev,0,DEFAULT_PRECISION,DEFAULT_PRECISION);
	solver.setInitialGuess(0.,0.,4.);

	double result	= solver.Solve();
	itsEqShift		= result;
	itsEqVol		= func.GetVol(itsEqShift);
	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: GetProxyShift
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
double ARM_ProcessBuilderPDE::SetProxyVol(double strike,const ARM_VanillaSecurityDensity& density)
{
	bool success (true);
	double tex			= Horizon()/K_YEAR_LEN;
	double mkt			= density.getDensityFunctor()->Call_Option(strike, itsFwd, tex);
	
	SLN_Approx func(itsFwd,tex,strike,mkt,strike,mkt,density);

	itsEqShift		= 0.;
	itsEqVol		= func.GetVol(itsEqShift);
	return itsEqShift;
}




////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: GetIthMoment
///	Returns: double
///	Action : cf name
////////////////////////////////////////////////////
double ARM_ProcessBuilderPDE::GetIthMoment(size_t n) const
{
	int nbPoints = 16;
	double globalVol = sqrt(itsGlobalVariance);
	
	double numStateMin = -itsNbStdDev*globalVol;
	double numStateMax = itsNbStdDev*globalVol;
		
	GaussLegendre_Coefficients glc( GLPOINTSNB, numStateMin, numStateMax);
		
	double result = 0;
	double resultByState=0;
	double state = 0;
	double liborByState = 0;

	double OneOnVar = 1./itsGlobalVariance;
	double OneOnStdDev = sqrt( OneOnVar );
		
	std::vector<double> states (GLPOINTSNB);
	for (size_t i(0); i<GLPOINTSNB; i++)
		states[i] = glc.get_point(i);

	ARM_GP_VectorPtr dfRatio = DFRatioAtResetDate(itsDates.size()-1,states);
		
	for( i=0 ; i<GLPOINTSNB; ++i )
	{
		state = glc.get_point(i);
		liborByState  = ((*dfRatio)[i]-1.)/itsDelta;
		resultByState = 1.;
		for (size_t k=0;k<n;k++)
			resultByState *= ( liborByState - itsFwd );
		result += glc.get_weight(i) * resultByState * exp(-0.5*state*state*OneOnVar);
	}
	result *= OneOnStdDev*INVSQRT2PI;

	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderSmiledPDE
///	Routine: Calibrate
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
double ARM_ProcessBuilderPDE::RescaleVolatilityCalib(ARM_Curve* volatilityCalib)
{
	std::vector<double> old = volatilityCalib->GetOrdinates();
	std::vector<double> date = volatilityCalib->GetAbscisses();

	size_t size=old.size();
	for(size_t index=size-1;index>0;index--)
		old[index] *= GetDiffusionTerm(date[index-1]);

	double sum=0.;
	for(index=1;index<size;index++)
		sum += old[index] * old[index] * (date[index]-date[index-1])/K_YEAR_LEN;

	volatilityCalib->SetOrdinates(old);

	return sum;

	
	
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: ComputeProxyForCalib
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_ProcessBuilderPDE::ComputeProxyForCalib(size_t calibProxy,const ARM_VanillaSecurityDensity& density,ARM_Curve* volatilityCalib,double strike)
{
	switch (calibProxy)
	{
		case ARM_ModelParamsSmiled::ATM:
			ComputeProxyAtm(density,false);
			break;
		case ARM_ModelParamsSmiled::AtmBlack:
			ComputeProxyAtm(density,true);
			break;
		case ARM_ModelParamsSmiled::MomentMatching:
			ComputeProxyMoment(density);
			break;
		case ARM_ModelParamsSmiled::LocalVolatility:
			ComputeProxyLocal(density,volatilityCalib,false,false);
			break;
		case ARM_ModelParamsSmiled::LocalVolatilityBlack:
			ComputeProxyLocal(density,volatilityCalib,false,true);
			break;
		case ARM_ModelParamsSmiled::LocalVolatilityWithRescaling:
			ComputeProxyLocal(density,volatilityCalib,true,false);
			break;
		case ARM_ModelParamsSmiled::LocalVolatilityBlackWithRescaling:
			ComputeProxyLocal(density,volatilityCalib,true,true);
			break;
		case ARM_ModelParamsSmiled::EffectiveSkew:
			ComputeProxyEffective(density,volatilityCalib,false);
			break;
		case ARM_ModelParamsSmiled::EffectiveSkewWithRescaling:
			ComputeProxyEffective(density,volatilityCalib,true);
			break;
		case ARM_ModelParamsSmiled::GaussBasketAtm:
			ComputeProxyGauss(density,true,strike);
			break;
		case ARM_ModelParamsSmiled::GaussBasketMoneyness:
			ComputeProxyGauss(density,false,strike);
			break;
	}
}

double SLN_Approx::operator() (double shift) const 
{
	double price2	= BS(itsFwd + shift, itsStrike2 + shift, itsTime, GetVol(shift));
	return price2-itsTarget2;
}

double SLN_Approx::GetVol(double shift) const
{
	bool success(true);
	double sigma = VanillaImpliedVol_BS(itsFwd + shift, itsStrike1 + shift, itsTime, itsTarget1, 1, NULL, &success);
	return sigma;
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: IdxFromDate
///	Returns: size_t
///	Action : returns the position of date into dateVector
///  throw exception if not found
////////////////////////////////////////////////////
size_t ARM_ProcessBuilderPDE::IdxFromDate( const std::vector<double>& dateVector, double date ) const
{
	size_t k=0;
	while( k<dateVector.size() && dateVector.Elt(k) < date )
		++k;

	if( k >= dateVector.size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "ARM_ProcessBuilderPDE::IdxFromDate: date not found in the map!");

	return k;
}


////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: DFRatio
///	Returns: 
///	Action : Returns martingale process at t
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_ProcessBuilderPDE::DFRatio(double mat,const ARM_PricingStatesPtr& states, size_t modelNb, size_t k, bool extrapol) const
{
	double t;
	if(mat<0)
	    ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ProcessBuilderPDE,DFRatio: Time index is negative!");
	if(mat>Horizon())
	{
		if (extrapol)
			t = Horizon();
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ProcessBuilderPDE,DFRatio: Time index bigger than horizon!");
	}
	else
	{
		t = mat;
	}
	if(IsDateInSchedule(t))
	{
		size_t index = IdxFromValueWithTol(itsDates,t,1);
		return DFRatioAtResetDate(index,states,modelNb,k);
	}
	else
	{
		size_t iRight=0;
		while( t > itsDates[iRight] ) {iRight++;};
		size_t iLeft = iRight - 1;
		double ratio = ( t - itsDates[iLeft] ) / ( itsDates[iRight] - itsDates[iLeft] );

		ARM_VectorPtr dxLeft	= DFRatioAtResetDate(iLeft,states,modelNb,k);
		ARM_VectorPtr dxRight	= DFRatioAtResetDate(iRight,states,modelNb,k);
		ARM_VectorPtr interpol( new std::vector<double>(states->size()) );

		std::vector<double>::iterator iter2 = dxRight->begin();
		std::vector<double>::iterator iter3 = interpol->begin();

		for(std::vector<double>::iterator iter1 = dxLeft->begin(); iter1!=dxLeft->end() ; ++iter1, ++iter2, ++iter3)
		{
			(*iter3) = (*iter1) + ratio * ((*iter2)-(*iter1) );
		}
		return interpol;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: Rate
///	Returns: 
///	Action : Returns martingale process at t
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_ProcessBuilderPDE::Rate(double mat,const ARM_PricingStatesPtr& states, size_t modelNb,size_t k, bool extrapol) const
{
	double t;
	if(mat<0)
	    ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ProcessBuilderPDE,DFRatio: Time index is negative!");
	if(mat>Horizon())
	{
		if (extrapol)
			t = Horizon();
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ProcessBuilderPDE,DFRatio: Time index bigger than horizon!");
	}
	else
	{
		t = mat;
	}
	if(IsDateInSchedule(t))
	{
		size_t index = IdxFromValueWithTol(itsDates,t,1);
		std::vector<double>& result = RateAtResetDate(index,states,modelNb,k);
		return ARM_GP_VectorPtr(result);
	}
	else
	{
		size_t iRight=0;
		while( t > itsDates[iRight] ) {iRight++;};
		size_t iLeft = iRight - 1;
		double ratio = ( t - itsDates[iLeft] ) / ( itsDates[iRight] - itsDates[iLeft] );
		std::vector<double>& interpol = RateBetweenResetDates(ratio,iLeft,iRight,states,modelNb,k);

		return ARM_GP_VectorPtr(interpol);
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: DxRate
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_ProcessBuilderPDE::DxRate(double mat,const ARM_PricingStatesPtr& states, size_t modelNb, size_t k, bool extrapol) const
{
	double t;
	if(mat<0)
	    ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ProcessBuilderPDE,DFRatio: Time index is negative!");
	if(mat>Horizon())
	{
		if (extrapol)
			t = Horizon();
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ProcessBuilderPDE,DFRatio: Time index bigger than horizon!");
	}
	else
	{
		t = mat;
	}
	if(IsDateInSchedule(t))
	{
		size_t index = IdxFromValueWithTol(itsDates,t,1);
		std::vector<double>& result = DxRateAtResetDate(index,states,modelNb,k);
		return ARM_GP_VectorPtr(result);
	}
	else
	{
		size_t iRight=0;
		while( t > itsDates[iRight] ) {iRight++;};
		size_t iLeft = iRight - 1;
		double ratio = ( t - itsDates[iLeft] ) / ( itsDates[iRight] - itsDates[iLeft] );
		std::vector<double>& interpol = DxRateBetweenResetDates(ratio,iLeft,iRight,states,modelNb,k);

		return ARM_GP_VectorPtr(interpol);
	}
}
////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: RateAtResetDate
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
std::vector<double>& ARM_ProcessBuilderPDE::RateAtResetDate(size_t index,const ARM_PricingStatesPtr& states, size_t modelNb, size_t k) const
{
	size_t i, j;
	size_t newSize							= states->size();
	ARM_GP_MatrixPtr newNumMethodStates		= states->GetModelStates();

	double curNumMethState, prevNumMethState, nextNumMethState;
	double storageXmin = itsGrid[0];
	double storageXmax = itsGrid[itsGridSize - 1];
	double storageDx = (storageXmax - storageXmin) / (itsGridSize - 1);

	std::vector<double>& result = new std::vector<double>( newSize );

	for (i = 0; i<newSize; i++)
	{	
		curNumMethState = newNumMethodStates->Elt(modelNb+k,i);
	
		/// flat extrapol if out of range
		if (curNumMethState<=storageXmin)
			result->Elt(i) =  itsRate[index]->Elt(0);
		else if (curNumMethState>=storageXmax)
			result->Elt(i) =  itsRate[index]->Elt(itsGridSize-1);
		else
		{
			//standard interpol
			j = (size_t)ceil((curNumMethState - storageXmin) / storageDx) ;
			nextNumMethState = itsGrid[j];
			prevNumMethState = itsGrid[j-1];
			
			result->Elt(i) = (  (curNumMethState - prevNumMethState) * itsRate[index]->Elt(j)
							  + (nextNumMethState - curNumMethState) * itsRate[index]->Elt(j-1) ) / (nextNumMethState - prevNumMethState);
		
		}
	}
	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: RateAtResetDate
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
std::vector<double>& ARM_ProcessBuilderPDE::RateBetweenResetDates(double ratio,size_t iLeft,size_t iRight,const ARM_PricingStatesPtr& states, size_t modelNb, size_t k) const
{
	size_t i, j;
	size_t newSize							= states->size();
	ARM_GP_MatrixPtr newNumMethodStates		= states->GetModelStates();

	double curNumMethState, prevNumMethState, nextNumMethState;
	double storageXmin = itsGrid[0];
	double storageXmax = itsGrid[itsGridSize - 1];
	double storageDx = (storageXmax - storageXmin) / (itsGridSize - 1);

	std::vector<double>& result = new std::vector<double>( newSize );

	for (i = 0; i<newSize; i++)
	{	
		curNumMethState = newNumMethodStates->Elt(modelNb+k,i);
	
		/// flat extrapol if out of range
		if (curNumMethState<=storageXmin)
			result->Elt(i) =  (1 - ratio) * itsRate[iLeft]->Elt(0) 
								+ ratio * itsRate[iRight]->Elt(0);
		else if (curNumMethState>=storageXmax)
			result->Elt(i) =  (1 - ratio) * itsRate[iLeft]->Elt(itsGridSize-1)
								+ ratio * itsRate[iRight]->Elt(itsGridSize-1);
		else
		{
			//standard interpol
			j = (size_t)ceil((curNumMethState - storageXmin) / storageDx) ;
			nextNumMethState = itsGrid[j];
			prevNumMethState = itsGrid[j-1];
			
			result->Elt(i) = (  (curNumMethState - prevNumMethState) * ((1 - ratio) * itsRate[iLeft]->Elt(j) + ratio * itsRate[iRight]->Elt(j))
							  + (nextNumMethState - curNumMethState) * ((1 - ratio) * itsRate[iLeft]->Elt(j-1) + ratio * itsRate[iRight]->Elt(j-1))
							  ) / (nextNumMethState - prevNumMethState);
		
		}
	}
	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: RateAtResetDate
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
std::vector<double>& ARM_ProcessBuilderPDE::DxRateBetweenResetDates(double ratio,size_t iLeft,size_t iRight,const ARM_PricingStatesPtr& states, size_t modelNb, size_t k) const
{
	size_t i, j;
	size_t newSize							= states->size();
	ARM_GP_MatrixPtr newNumMethodStates		= states->GetModelStates();

	double curNumMethState, prevNumMethState, nextNumMethState;
	double storageXmin = itsGrid[0];
	double storageXmax = itsGrid[itsGridSize - 1];
	double storageDx = (storageXmax - storageXmin) / (itsGridSize - 1);

	std::vector<double>& result = new std::vector<double>( newSize );

	for (i = 0; i<newSize; i++)
	{	
		curNumMethState = newNumMethodStates->Elt(modelNb+k,i);
	
		/// flat extrapol if out of range
		if (curNumMethState<=storageXmin)
			result->Elt(i) =  (1 - ratio) * itsDxRate[iLeft]->Elt(0) 
								+ ratio * itsDxRate[iRight]->Elt(0);
		else if (curNumMethState>=storageXmax)
			result->Elt(i) =  (1 - ratio) * itsDxRate[iLeft]->Elt(itsGridSize-1)
								+ ratio * itsDxRate[iRight]->Elt(itsGridSize-1);
		else
		{
			//standard interpol
			j = (size_t)ceil((curNumMethState - storageXmin) / storageDx) ;
			nextNumMethState = itsGrid[j];
			prevNumMethState = itsGrid[j-1];
			
			result->Elt(i) = (  (curNumMethState - prevNumMethState) * ((1 - ratio) * itsDxRate[iLeft]->Elt(j) + ratio * itsDxRate[iRight]->Elt(j))
							  + (nextNumMethState - curNumMethState) * ((1 - ratio) * itsDxRate[iLeft]->Elt(j-1) + ratio * itsDxRate[iRight]->Elt(j-1))
							  ) / (nextNumMethState - prevNumMethState);
		
		}
	}
	return result;
}
////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: DxRateAtResetDate
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
std::vector<double>& ARM_ProcessBuilderPDE::DxRateAtResetDate(size_t index,const ARM_PricingStatesPtr& states, size_t modelNb, size_t k) const
{
	size_t i, j;
	size_t newSize							= states->size();
	ARM_GP_MatrixPtr newNumMethodStates		= states->GetModelStates();

	double curNumMethState, prevNumMethState, nextNumMethState;
	double storageXmin = itsGrid[0];
	double storageXmax = itsGrid[itsGridSize - 1];
	double storageDx = (storageXmax - storageXmin) / (itsGridSize - 1);

	std::vector<double>& result = new std::vector<double>( newSize );

	for (i = 0; i<newSize; i++)
	{	
		curNumMethState = newNumMethodStates->Elt(modelNb+k,i);
	
		/// flat extrapol if out of range
		if (curNumMethState<=storageXmin)
			result->Elt(i) =  itsDxRate[index]->Elt(0);
		else if (curNumMethState>=storageXmax)
			result->Elt(i) =  itsDxRate[index]->Elt(itsGridSize-1);
		else
		{
			//standard interpol
			j = (size_t)ceil((curNumMethState - storageXmin) / storageDx) ;
			nextNumMethState = itsGrid[j];
			prevNumMethState = itsGrid[j-1];
			
			result->Elt(i) = (  (curNumMethState - prevNumMethState) * itsDxRate[index]->Elt(j)
							  + (nextNumMethState - curNumMethState) * itsDxRate[index]->Elt(j-1) ) / (nextNumMethState - prevNumMethState);
		
		}
	}
	return result;
}
////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: DFRatioAtResetDate
///	Returns: 
///	Action : returns drift adjustment term for current martingale process (when date is stored)
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_ProcessBuilderPDE::DFRatioAtResetDate(size_t index,const std::vector<double>& states) const
{
	size_t i, j;
	size_t newSize							= states.size();
	
	double curNumMethState, prevNumMethState, nextNumMethState;
	double storageXmin = itsGrid[0];
	double storageXmax = itsGrid[itsGridSize - 1];
	double storageDx = (storageXmax - storageXmin) / (itsGridSize - 1);

	std::vector<double>& result = new std::vector<double>( newSize );

	for (i = 0; i<newSize; i++)
	{	
		curNumMethState = states[i];
	
		/// flat extrapol if out of range
		if (curNumMethState<=storageXmin)
			result->Elt(i) =  itsDFRatio[index]->Elt(0);
		else if (curNumMethState>=storageXmax)
			result->Elt(i) =  itsDFRatio[index]->Elt(itsGridSize-1);
		else
		{
			//standard interpol
			j = (size_t)ceil((curNumMethState - storageXmin) / storageDx) ;
			nextNumMethState = itsGrid[j];
			prevNumMethState = itsGrid[j-1];
			
			result->Elt(i) = (  (curNumMethState - prevNumMethState) * itsDFRatio[index]->Elt(j)
							  + (nextNumMethState - curNumMethState) * itsDFRatio[index]->Elt(j-1) ) / (nextNumMethState - prevNumMethState);
		
		}
	}
	return ARM_GP_VectorPtr(result);
}


////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: DxLogDFRatioAtResetDate
///	Returns: 
///	Action : returns drift adjustment term for current martingale process (when date is stored)
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_ProcessBuilderPDE::DFRatioAtResetDate(size_t index,const ARM_PricingStatesPtr& states, size_t modelNb, size_t k) const
{
	size_t i, j;
	size_t newSize							= states->size();
	ARM_GP_MatrixPtr newNumMethodStates		= states->GetModelStates();

	double curNumMethState, prevNumMethState, nextNumMethState;
	double storageXmin = itsGrid[0];
	double storageXmax = itsGrid[itsGridSize - 1];
	double storageDx = (storageXmax - storageXmin) / (itsGridSize - 1);

	std::vector<double>& result = new std::vector<double>( newSize );

	for (i = 0; i<newSize; i++)
	{	
		curNumMethState = newNumMethodStates->Elt(modelNb+k,i);
	
		/// flat extrapol if out of range
		if (curNumMethState<=storageXmin)
			result->Elt(i) =  itsDFRatio[index]->Elt(0);
		else if (curNumMethState>=storageXmax)
			result->Elt(i) =  itsDFRatio[index]->Elt(itsGridSize-1);
		else
		{
			//standard interpol
			j = (size_t)ceil((curNumMethState - storageXmin) / storageDx) ;
			nextNumMethState = itsGrid[j];
			prevNumMethState = itsGrid[j-1];
			
			result->Elt(i) = (  (curNumMethState - prevNumMethState) * itsDFRatio[index]->Elt(j)
							  + (nextNumMethState - curNumMethState) * itsDFRatio[index]->Elt(j-1) ) / (nextNumMethState - prevNumMethState);
		
		}
	}
	return ARM_GP_VectorPtr(result);
}


////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: DxLogDFRatio
///	Returns: 
///	Action : returns drift adjustment term for current martingale process
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_ProcessBuilderPDE::DxLogDFRatio(double mat,const ARM_PricingStatesPtr& states, size_t modelNb, size_t k, bool extrapol) const
{
	double t;
	if(mat<0)
	    ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ProcessBuilderPDE,DFRatio: Time index is negative!");
	if(mat>Horizon())
	{
		if (extrapol)
			t = Horizon();
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ProcessBuilderPDE,DFRatio: Time index bigger than horizon!");
	}
	else
	{
		t = mat;
	}
	if(IsDateInSchedule(t))
	{
		size_t index = IdxFromValueWithTol(itsDates,t,1);
		return DxLogDFRatioAtResetDate(index,states,modelNb,k);
	}
	else
	{
		size_t iRight=0;
		while( t > itsDates[iRight] ) {iRight++;};
		size_t iLeft = iRight - 1;
		double ratio = ( t - itsDates[iLeft] ) / ( itsDates[iRight] - itsDates[iLeft] );

		ARM_VectorPtr dxLeft	= DxLogDFRatioAtResetDate(iLeft,states,modelNb,k);
		ARM_VectorPtr dxRight	= DxLogDFRatioAtResetDate(iRight,states,modelNb,k);
		ARM_VectorPtr interpol( new std::vector<double>(states->size()) );

		std::vector<double>::iterator iter2 = dxRight->begin();
		std::vector<double>::iterator iter3 = interpol->begin();

		for(std::vector<double>::iterator iter1 = dxLeft->begin(); iter1!=dxLeft->end() ; ++iter1, ++iter2, ++iter3)
		{
			(*iter3) = (*iter1) + ratio * ((*iter2)-(*iter1) );
		}
		return interpol;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: DxLogDFRatioAtResetDate
///	Returns: 
///	Action : returns drift adjustment term for current martingale process (when date is stored)
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_ProcessBuilderPDE::DxLogDFRatioAtResetDate(size_t index,const ARM_PricingStatesPtr& states, size_t modelNb, size_t k) const
{
	size_t i, j;
	size_t newSize							= states->size();
	ARM_GP_MatrixPtr newNumMethodStates		= states->GetModelStates();

	double curNumMethState, prevNumMethState, nextNumMethState;
	double storageXmin = itsGrid[0];
	double storageXmax = itsGrid[itsGridSize - 1];
	double storageDx = (storageXmax - storageXmin) / (itsGridSize - 1);

	std::vector<double>& result = new std::vector<double>( newSize );

	for (i = 0; i<newSize; i++)
	{
		curNumMethState = newNumMethodStates->Elt(modelNb+k,i);
	
		/// flat extrapol if out of range
		if (curNumMethState<=storageXmin)
			result->Elt(i) =  itsDxLogDFRatio[index]->Elt(0);
		else if (curNumMethState>=storageXmax)
			result->Elt(i) =  itsDxLogDFRatio[index]->Elt(itsGridSize-1);
		else
		{
			//standard interpol
			j = (size_t)ceil((curNumMethState - storageXmin) / storageDx) ;
			nextNumMethState = itsGrid[j];
			prevNumMethState = itsGrid[j-1];
			result->Elt(i) = (  (curNumMethState - prevNumMethState) * itsDxLogDFRatio[index]->Elt(j)
							  + (nextNumMethState - curNumMethState) * itsDxLogDFRatio[index]->Elt(j-1) ) / (nextNumMethState - prevNumMethState);
		}
	}
	return ARM_GP_VectorPtr(result);
}



////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: setNumericalModelFitter
///	Returns: void
///	Action : updates PDE params
////////////////////////////////////////////////////

void ARM_ProcessBuilderPDE::SetPDEParams(size_t gridSize,double nbStdDev,double theta)
{
	itsGridSize=gridSize;
	itsNbStdDev=nbStdDev;
	itsTheta=theta;
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: Calibrate
///	Returns: 
///	Action : computes and stores E[f(XT)|Ft] (implies solving a 1F PDE)
////////////////////////////////////////////////////
void ARM_ProcessBuilderPDE::Calibrate(	double horizonDate,
												const ARM_VanillaSecurityDensity& density,
												const std::vector<double>& evalDates,
												const std::vector<double>& numDates,
												bool doPDE)
{
	ARM_ProcessBuilder::Calibrate( horizonDate, density, evalDates, numDates);

	itsGlobalVariance		= IntegratedVarianceFromZero(Horizon());

	if (doPDE)
	{
		if (evalDates[0]>0.) itsDates.push_back(0.);
		for (size_t t=0;t<evalDates.size();t++){
			if (evalDates[t]<horizonDate)
				itsDates.push_back(evalDates[t]);
		}
	}
	itsDates.push_back(horizonDate);

	size_t size = itsGridSize;
	size_t n    = itsDates.size();
	size_t mid  = size/2;

	itsRate.resize(n);
	
	BuildGrid();
	
	ARM_GP_VectorPtr func = GetTerminalFunc(density);
	for(size_t index=n;index>0;index--)
	{
		if (index<n)
			Backward(itsDates[index-1],itsDates[index],func,numDates);
		StoreRate(index-1,func);
	}
	
	if (doPDE)
	{
		itsDxRate.resize(n);
		itsDFRatio.resize(n);
		itsDxLogDFRatio.resize(n);

		itsNumFwd = (*func)[size/2];
		for(index=n;index>0;index--)
		{
			if (itsWithRescalling)
				RescaleRate(index-1,itsFwd/itsNumFwd);
			StoreDxRate(index-1,itsRate[index-1]);
			StoreDFRatio(index-1,itsRate[index-1]);
			StoreDxLogDFRatio(index-1,itsDFRatio[index-1]);
		}
	}

	FreeMatrix();
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: ComputeProxyAtm
///	Returns: 
///	Action :
////////////////////////////////////////////////////
void ARM_ProcessBuilderPDE::ComputeProxyAtm(const ARM_VanillaSecurityDensity& density,bool black)
{
	if (black)
		itsEqShift			= SetProxyVol(itsFwd,density);
	else
		itsEqShift			= GetProxyShift(itsFwd-1e-4,itsFwd+1e-4,density);
	
	itsEqScaling		= itsEqVol*sqrt(Horizon()/K_YEAR_LEN/itsGlobalVariance);
	itsDiffusionTerm	= ( itsEqShift + itsFwd ) * itsEqScaling;
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: ComputeProxyGauss
///	Returns: 
///	Action :
////////////////////////////////////////////////////
void ARM_ProcessBuilderPDE::ComputeProxyGauss(const ARM_VanillaSecurityDensity& density,bool atm,double moneyness)
{
	bool success (true);
	double tex			= Horizon()/K_YEAR_LEN;
	double mkt			= density.getDensityFunctor()->Call_Option(itsFwd, itsFwd, tex);
	double volatm		= VanillaImpliedVol_N(itsFwd,mkt,itsFwd,tex,1,NULL, &success);
	if (atm || fabs(moneyness) < K_DOUBLE_TOL)
	{
		itsEqVol	= volatm;
	}
	else
	{
		double K	= CC_Max<double>(itsFwd + moneyness * volatm*sqrt(tex),0.01);
		mkt			= density.getDensityFunctor()->Call_Option(K, itsFwd, tex);
		double vol	= VanillaImpliedVol_N(itsFwd,mkt,K,tex,1,NULL, &success);
		itsEqVol	= vol;
	}
	itsEqShift			= 0;
	itsEqScaling		= itsEqVol*sqrt(Horizon()/K_YEAR_LEN/itsGlobalVariance);
	itsDiffusionTerm	= itsEqScaling;
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: Calibrate
///	Returns: 
///	Action : computes and stores E[f(XT)|Ft] (implies solving a 1F PDE)
////////////////////////////////////////////////////
void ARM_ProcessBuilderPDE::ComputeProxyMoment(const ARM_VanillaSecurityDensity& density)
{
	double m2			= GetIthMoment(2);
	double tex			= Horizon()/K_YEAR_LEN;
	
	itsEqShift			= 0.;
	itsEqVol			= sqrt(log( 1 + m2 / (itsFwd*itsFwd))/tex);
	itsEqScaling		= itsEqVol*sqrt(tex)/sqrt(itsGlobalVariance);
	itsDiffusionTerm	= ( itsEqShift + itsFwd ) * itsEqScaling;
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: ComputeProxyLocal
///	Returns: 
///	Action :
////////////////////////////////////////////////////
void ARM_ProcessBuilderPDE::ComputeProxyLocal(const ARM_VanillaSecurityDensity& density,ARM_Curve* volatilityCalib,bool rescaling,bool black)
{
	double cumVar = RescaleVolatilityCalib(volatilityCalib);
	if (black)
		itsEqShift			= SetProxyVol(itsFwd,density);
	else
		itsEqShift			= GetProxyShift(itsFwd-1e-4,itsFwd+1e-4,density);
	if (rescaling)
	{
		bool success (true);
		double tex = Horizon()/K_YEAR_LEN;
		double mkt = density.getDensityFunctor()->Call_Option(itsFwd, itsFwd, tex);
		double vol = VanillaImpliedVol_BS(itsEqShift+itsFwd,itsEqShift+itsFwd,tex, mkt, 1, NULL, &success);
		double factor = (itsEqShift+itsFwd)*vol*sqrt(tex/cumVar);
		itsEqScaling = factor;
	}
	itsDiffusionTerm = itsEqScaling;
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: ComputeProxyEffective
///	Returns: 
///	Action :
////////////////////////////////////////////////////
void ARM_ProcessBuilderPDE::ComputeProxyEffective(const ARM_VanillaSecurityDensity& density,ARM_Curve* volatilityCalib,bool rescaling)
{
	double cumVar		= RescaleVolatilityCalib(volatilityCalib);
	std::vector<double> vol	= volatilityCalib->GetOrdinates();
	std::vector<double> date	= volatilityCalib->GetAbscisses();
	size_t	size		= IdxFromValueWithTol(date,Horizon(),1) + 1;

	std::vector<double> v(size,0.);

	double sum1 = 0.;
	for(size_t index=1;index<size;index++)
	{
		sum1		= sum1
					+ v[index-1] * vol[index] * vol[index] *(date[index]-date[index-1])/K_YEAR_LEN
					+ vol[index] * vol[index] * vol[index] *vol[index] * (date[index]*date[index]-date[index-1]*date[index-1])/K_YEAR_LEN/K_YEAR_LEN/2;

		v[index]	= v[index-1] + vol[index] * vol[index] *(date[index]-date[index-1])/K_YEAR_LEN;
	}

	double z	=	1.;
	double sum2 = 0.;
	for(index=1;index<size;index++)
		sum2 = sum2 + (z * GetDLocalVolDRate(date[index-1]) + (1-z) * GetDLocalVolDRate(date[index])) * (
						v[index-1] * vol[index] * vol[index] *(date[index]-date[index-1])/K_YEAR_LEN
						+	vol[index] * vol[index] * vol[index] *vol[index] * (date[index]*date[index]-date[index-1]*date[index-1])/K_YEAR_LEN/K_YEAR_LEN/2);

	sum2 = sum2 / sum1 * itsFwd ;		// = b

	itsEqShift = (1-sum2)/sum2*itsFwd ;

	if (rescaling)
	{
		bool success (true);
		double tex = Horizon()/K_YEAR_LEN;
		double mkt = density.getDensityFunctor()->Call_Option(itsFwd, itsFwd, tex);
		double vol = VanillaImpliedVol_BS(itsEqShift+itsFwd,itsEqShift+itsFwd,tex, mkt, 1, NULL, &success);
		double factor = (itsEqShift+itsFwd)*vol*sqrt(tex/cumVar);
		itsEqScaling = factor;
	}
	else
		itsEqScaling = 1.;
	
	itsDiffusionTerm = itsEqScaling;
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: InitTerminalLaw
///	Returns: 
///	Action : computes functional at reset date
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ProcessBuilderPDE::GetTerminalFunc(const ARM_VanillaSecurityDensity& terminalDensity)
{
	double globalStdDev		= sqrt(GlobalVariance());

	std::vector<double> delta		= terminalDensity.getInterestTerms();
	
	ARM_GP_VectorPtr aux	= ARM_GP_VectorPtr(new ARM_GP_Vector(itsGridSize));

	for (size_t i = 0 ; i < itsGridSize ; i++)
		(*aux)[i] = NormalCDF( itsGrid[i] / globalStdDev );

	ARM_VectorPtr result	= terminalDensity.getDensityFunctor()->Quantile(aux,itsFwd,Horizon()/K_YEAR_LEN);

	if (delta.size()==1)
		for (i = 0 ; i < itsGridSize ; i++)
			(*result)[i] = CC_Max((*result)[i],(1e-6-1.)/delta[0]);

	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: Backward
///	Returns: 
///	Action : computes and stores E[f(XT)|Ft] (implies solving a 1F PDE)
////////////////////////////////////////////////////
void ARM_ProcessBuilderPDE::Backward(double first,double last,ARM_GP_VectorPtr& func,const std::vector<double>& numDates)
{
	size_t iFirst	= IdxFromDate(numDates,first/K_YEAR_LEN);
	size_t iLast	= IdxFromDate(numDates,last/K_YEAR_LEN);
	size_t size		= GridSize();
	double dx		= 2*NbStdDev()*sqrt(GlobalVariance())/GridSize();
	double dt,sigma,a,
		   lambda	= MeanReversion();
	
//BE CAREFUL : volatility must be constant between first and last
	sigma  = Volatility((first+last)/2.);
	a      = 0.5*sigma*sigma/dx/dx;

	for(size_t i = 0 ; i < size-2 ; i++)
	{
		itsL1[i] = Theta()*(a-0.5*(-lambda*itsGrid[i+1]/dx));
		itsL2[i] = -Theta()*(2.*a);
		itsL3[i] = Theta()*(a+0.5*(-lambda*itsGrid[i+1]/dx));
		itsR1[i] = (1-Theta())*(a-0.5*(-lambda*itsGrid[i+1]/dx));
		itsR2[i] = -(1-Theta())*(2*a);	
		itsR3[i] = (1-Theta())*(a+0.5*(-lambda*itsGrid[i+1]/dx));
	}
	itsL1[0]		  =0;
	itsL3[size-2-1]  =0;
	itsL2[0]        += Theta()*(a-0.5*(-lambda*itsGrid[1]/dx));
	itsL2[size-2-1] += Theta()*(a+0.5*(-lambda*itsGrid[size-2]/dx));
		

	size_t t      = iLast;
	while(--t>=iFirst && t!=0)
	{
		dt	   = (numDates[t+1]-numDates[t]);
		for (i=0;i<size-2;i++)
			itsR[i]=-( dt*(itsR1[i]*(*func)[i]+itsR2[i]*(*func)[i+1]+itsR3[i]*(*func)[i+2])+(*func)[i+1]);

		itsR[0]		+= dt*Theta()*(a - 0.5*(-lambda*itsGrid[1]/dx)) * ((*func)[1]-(*func)[0]);
		itsR[size-2-1]	-= dt*Theta()*(a + 0.5*(-lambda*itsGrid[size-2]/dx)) * ((*func)[size-2+1]-(*func)[size-2]);

		//inversion
		double bet;

		itsXG[0]=itsR[0]/(bet=(dt*itsL2[0]-1.));
		for (size_t j=1;j<size-2;j++)
		{
			itsGam[j]=dt*itsL3[j-1]/bet;
			bet=(dt*itsL2[j]-1.)-dt*itsL1[j]*itsGam[j];
			itsXG[j]=(itsR[j]-dt*itsL1[j]*itsXG[j-1])/bet;
		}
		for (j=size-2-1;j>0;j--) {
			itsXG[j-1] -= itsGam[j]*itsXG[j]; 
		}

		(*func)[0]		= itsXG[0]-((*func)[1]-(*func)[0]);
		(*func)[size-1]	= itsXG[size-2-1]+((*func)[size-1]-(*func)[size-2]);

		for (i=1;i<size-1;i++) (*func)[i]=itsXG[i-1];
	}	

}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: StoreDFRatio
///	Returns: 
///	Action : stores DF ratios at predefined dates
////////////////////////////////////////////////////

void ARM_ProcessBuilderPDE::StoreDFRatio(size_t index_storage,const ARM_GP_VectorPtr& func)
{
	ARM_GP_VectorPtr ratio	= ARM_GP_VectorPtr(new ARM_GP_Vector(itsGridSize));

	for (size_t k=0;k<itsGridSize;k++) 
		(*ratio)[k]			= 1.+itsDelta*(*func)[k];

	itsDFRatio[index_storage]=ratio;

}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: StoreRate
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

void ARM_ProcessBuilderPDE::StoreRate(size_t index_storage,const ARM_GP_VectorPtr& func){

	ARM_GP_VectorPtr ratio	= ARM_GP_VectorPtr(new ARM_GP_Vector(itsGridSize));

	for (size_t k=0;k<itsGridSize;k++) 
		(*ratio)[k]			= (*func)[k];

	itsRate[index_storage]=ratio;

}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: StoreDxRate
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

void ARM_ProcessBuilderPDE::StoreDxRate(size_t index_storage,const ARM_GP_VectorPtr& func){

	ARM_GP_VectorPtr dx	= ARM_GP_VectorPtr(new ARM_GP_Vector(itsGridSize));

	for (size_t k=1;k<itsGridSize-1;k++)
		(*dx)[k]			= ( (*func)[k+1]-(*func)[k-1]) / (itsGrid[k+1]-itsGrid[k-1]);

	(*dx)[0]				= ( (*func)[1]-(*func)[0]) / (itsGrid[1]-itsGrid[0]);
	(*dx)[itsGridSize-1]	= ( (*func)[itsGridSize-1]-(*func)[itsGridSize-2] ) / (itsGrid[itsGridSize-1]-itsGrid[itsGridSize-2]);

	itsDxRate[index_storage]=dx;
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: StoreDFRatio
///	Returns: 
///	Action : stores DF ratios at predefined dates
////////////////////////////////////////////////////

void ARM_ProcessBuilderPDE::StoreDxLogDFRatio(size_t index_storage,const ARM_GP_VectorPtr& func){

	ARM_GP_VectorPtr dxlog	= ARM_GP_VectorPtr(new ARM_GP_Vector(itsGridSize));

	for (size_t k=1;k<itsGridSize-1;k++)
		(*dxlog)[k]			= (log((*func)[k+1])-log((*func)[k-1])) / (itsGrid[k+1]-itsGrid[k-1]);

	(*dxlog)[0]				= (log((*func)[1])-log((*func)[0])) / (itsGrid[1]-itsGrid[0]);
	(*dxlog)[itsGridSize-1] = (log((*func)[itsGridSize-1])-log((*func)[itsGridSize-2])) / (itsGrid[itsGridSize-1]-itsGrid[itsGridSize-2]);

	itsDxLogDFRatio[index_storage]=dxlog;

}


////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: Rescale
///	Returns: 
///	Action : Resecale the forward value
////////////////////////////////////////////////////

void ARM_ProcessBuilderPDE::RescaleRate(size_t index_storage, double scale){


	for (size_t k=0;k<itsGridSize;k++) 
		(*itsRate[index_storage])[k] *= scale;

}


////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: GetIndex
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
size_t ARM_ProcessBuilderPDE::GetIndex(double x, const std::vector<double>& v) const
{
	size_t size = v.size();

	if (v.size()>0)
	{
		size_t index=0;
		if (x > v[size-1])
			index = size-1;
		else
		{
			while (v[index]<x)
				index++;
		}
		return index;
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ProcessBuilderPDE::GetIndex : vector of size 0");
}
	

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderSLN
///	Routine: GetDiffusionTerm
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
double ARM_ProcessBuilderPDE::GetDiffusionTerm(double t) const
{
	double result;
	if(IsDateInSchedule(t))
	{
		size_t index = IdxFromValueWithTol(itsDates,t,1);
		size_t mid	 = GetIndex(itsFwd,(*itsRate[index]));
		result		 = (*itsDxRate[index])[mid];
	}
	else
		result=1.;
	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderSLN
///	Routine: GetDLocalVolDRate
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
double ARM_ProcessBuilderPDE::GetDLocalVolDRate(double t) const
{
	double result;
	if(IsDateInSchedule(t))
	{
		size_t index = IdxFromValueWithTol(itsDates,t,1);
		size_t mid	 = GetIndex(itsFwd,(*itsRate[index]));
		if (mid==0) mid++;
		
		result = ((*itsDxRate[index])[mid]-(*itsDxRate[index])[mid-1])
				/((*itsRate[index])[mid]-(*itsRate[index])[mid-1])
				/(*itsDxRate[index])[mid];
	}
	else
		result=1.;
	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: buildGrid
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_ProcessBuilderPDE::BuildGrid()
{
	size_t size = GridSize();
	double step = 2*NbStdDev()*sqrt(GlobalVariance())/(size-1.0);
	
	itsGrid.resize(size);

	for( size_t i=0 ; i  <size ; ++i )
		itsGrid[i] = step*(i-(size-1)/2.);

	itsL1.resize(size-2);
	itsL2.resize(size-2);
	itsL3.resize(size-2);
	//Right Tridiag
	itsR1.resize(size-2);
	itsR2.resize(size-2);
	itsR3.resize(size-2);
	//Right Vector
	itsR.resize(size-2);
	//For inversion
	itsGam.resize(size-2);
	itsXG.resize(size-2);
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderPDE
///	Routine: FreeMatrix
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_ProcessBuilderPDE::FreeMatrix()
{
	itsL1.resize(0);
	itsL2.resize(0);
	itsL3.resize(0);
	itsR1.resize(0);
	itsR2.resize(0);
	itsR3.resize(0);
	itsR.resize(0);
	itsGam.resize(0);
	itsXG.resize(0);
}

////////////////////////////////////////////////////
///	Class   : ARM_ProcessBuilderPDE
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_ProcessBuilderPDE::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
	os << ARM_ProcessBuilder::toString(indent,nextIndent);
	size_t size1 = itsGrid.size();
	size_t size2 = itsRate.size();
	size_t size3 = itsRate[0]->size();
	size_t size4 = itsRate[size2-1]->size();
    os << indent << "Rate Functionals\n";
	size_t n=20;
	
	if (size1>0 && size2>0 && (size3==size1) && (size4==size1) && size1>n)
	{
		for (size_t k=0;k<n;k++)
		{
			os << indent << "\t" << CC_NS(std,fixed) << itsGrid[k*(size1-1)/(n-1)]/sqrt(itsGlobalVariance);
			os << indent << "\t\t" << (*itsRate[0])[k*(size1-1)/(n-1)];
			os << indent << "\t\t" <<(*itsRate[size2-1])[k*(size1-1)/(n-1)] <<"\n";
		}
	}
	os <<"\n";
	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderSLN
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ProcessBuilderSLN::ARM_ProcessBuilderSLN(double mrs,const ARM_Curve& volCurve)
:	ARM_ProcessBuilder(mrs,volCurve),
	itsVol(0),
	itsShift(0),
	itsScaling(1)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderSLN
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ProcessBuilderSLN::ARM_ProcessBuilderSLN(const ARM_ProcessBuilderSLN& rhs)
:	ARM_ProcessBuilder(rhs),
	itsVol(rhs.itsVol),
	itsShift(rhs.itsShift),
	itsScaling(rhs.itsScaling)
{
}
	

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderSLN
///	Routine: DFRatio
///	Returns: 
///	Action : Returns martingale process at t
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_ProcessBuilderSLN::DFRatio(double t,const ARM_PricingStatesPtr& states, size_t modelNb, size_t k,bool extrapol) const
{
	size_t i;

	size_t newSize							= states->size();
	ARM_GP_MatrixPtr newNumMethodStates		= states->GetModelStates();

	std::vector<double>& result = new std::vector<double>( newSize );
	double curNumMethState;
	double var = itsScaling * itsScaling * IntegratedVarianceFromZero(t);
	for (i = 0; i<newSize; i++)
	{	
		curNumMethState = newNumMethodStates->Elt(modelNb+k,i);
		result->Elt(i) = 1. + itsDelta * ( (itsShift + itsFwd) * exp( itsScaling*curNumMethState - 0.5 * var) - itsShift);
	}
	return ARM_GP_VectorPtr(result);
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderSLN
///	Routine: DFRatio
///	Returns: 
///	Action : Returns martingale process at t
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_ProcessBuilderSLN::Rate(double t,const ARM_PricingStatesPtr& states, size_t modelNb, size_t k,bool extrapol) const
{
	size_t i;

	size_t newSize							= states->size();
	ARM_GP_MatrixPtr newNumMethodStates		= states->GetModelStates();

	std::vector<double>& result = new std::vector<double>( newSize );
	double curNumMethState;
	double var = itsScaling * itsScaling * IntegratedVarianceFromZero(t);
	for (i = 0; i<newSize; i++)
	{	
		curNumMethState = newNumMethodStates->Elt(modelNb+k,i);
		result->Elt(i) =  (itsShift + itsFwd) * exp( itsScaling*curNumMethState - 0.5 * var)-itsShift;
	}
	return ARM_GP_VectorPtr(result);
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderSLN
///	Routine: DFRatio
///	Returns: 
///	Action : Returns martingale process at t
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_ProcessBuilderSLN::DxRate(double t,const ARM_PricingStatesPtr& states, size_t modelNb, size_t k,bool extrapol) const
{
	size_t i;
	size_t newSize = states->size();
	
	ARM_GP_VectorPtr result	= Rate(t,states,modelNb,k,extrapol);
	for (i = 0; i<newSize; i++)
		result->Elt(i) =  itsScaling * (result->Elt(i) + itsShift);

	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderSLN
///	Routine: DFRatio
///	Returns: 
///	Action : Returns martingale process at t
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_ProcessBuilderSLN::DxLogDFRatio(double t,const ARM_PricingStatesPtr& states, size_t modelNb, size_t k, bool extrapol) const
{
	size_t i;
	size_t newSize	= states->size();
	
	ARM_GP_VectorPtr result	= DFRatio(t,states,modelNb,k,extrapol);
	for (i = 0; i<newSize; i++)
		result->Elt(i) = itsScaling * ( 1. - (1. -  itsShift * itsDelta) / result->Elt(i) );
	
	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_ProcessBuilderSLN
///	Routine: Calibrate
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_ProcessBuilderSLN::Calibrate(	double horizonDate,
												const ARM_VanillaSecurityDensity& terminalDensity,
												const std::vector<double>& evalDates,
												const std::vector<double>& numDates,
												bool doPDE)
{
	ARM_ProcessBuilder::Calibrate( horizonDate, terminalDensity, evalDates, numDates);

	ARM_ShiftedLNDensityFunctor* pdensity;
	pdensity = dynamic_cast <ARM_ShiftedLNDensityFunctor*> (&*terminalDensity.getDensityFunctor());

	if (!pdensity){
		ARM_THROW( ERR_INVALID_ARGUMENT,"ARM_ProcessBuilderSLN,Calibrate : BIG PROBLEM not shifted lognormal");
	}
	else
	{
		itsVol					= pdensity->getVol();
		itsShift				= pdensity->getShift();;
		double time = Horizon();
		if (Horizon() >= K_DOUBLE_TOL)
			itsScaling				= itsVol*sqrt(Horizon()/K_YEAR_LEN/IntegratedVarianceFromZero(Horizon()));
		else
			itsScaling				= 1.0;
	}	
}




CC_END_NAMESPACE()

