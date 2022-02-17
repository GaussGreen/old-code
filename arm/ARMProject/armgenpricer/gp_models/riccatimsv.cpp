/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  function solution of the Riccati Equation with Piece-Wise functions A(t), B(t) and C(t)
 *	df/dt = A(t) * f(t)� + B(t) * f(t) + C(t)
 *  
 *	\file riccati.h
 *
 *  \brief
 *
 *	\author  A. Triki
 *	\version 1.0
 *	\date October 2005
 */
 
#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include <cmath>
#include "gpmodels/riccatimsv.h"
#include <glob/expt.h>   // for the exceptions
#include <algorithm>


/// gpclosedforms
#include "gpclosedforms/gaussian_integrals.h"

using namespace std;


CC_BEGIN_NAMESPACE(ARM)


const double DELTA_LIMIT    = 1.0e-10;

////////////////////////////////////////////////////
///	Class  : ARM_Riccati::ARM_Funct0
///	Routine: ComputeA, ComputeCst
///	Returns: 
///	Action : No root case : 
///          compute B(t,Cst) Given the cst already computed 
////////////////////////////////////////////////////
double ARM_RiccatiMSV::ARM_Funct0MSV::ComputeB(double t1 , double t2) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Function Not Implemented!");
}


////////////////////////////////////////////////////
///	Class  : ARM_Riccati::ARM_Funct1
///	Routine: ComputeA, ComputeCst
///	Returns: 
///	Action : 1 root case :
///          compute B(t,Cst) Given the cst already computed 
////////////////////////////////////////////////////
double ARM_RiccatiMSV::ARM_Funct1MSV::ComputeB(double t1 , double t2) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Function Not Implemented!");
}

////////////////////////////////////////////////////
///	Class  : ARM_Riccati::ARM_Funct2Sup
///	Routine: ComputeA, ComputeCst
///	Returns: 
///	Action : 2 roots case and a constant solution
///          equals to the greater root :
///          compute B(t,Cst) Given the cst already computed 
////////////////////////////////////////////////////
double ARM_RiccatiMSV::ARM_Funct2SupMSV::ComputeB(double t1 , double t2) const
{
	return 0.;
}

////////////////////////////////////////////////////
///	Class  : ARM_Riccati::ARM_Funct2Inf
///	Routine: ComputeA, ComputeCst
///	Returns: 
///	Action : 2 roots case and a constant solution
///          equals to the lower root :
///          compute B(t,Cst) Given the cst already computed 
////////////////////////////////////////////////////
double ARM_RiccatiMSV::ARM_Funct2InfMSV::ComputeB(double t1 , double t2) const
{
	return 0.;
}


////////////////////////////////////////////////////
///	Class  : ARM_Riccati::ARM_Funct2Out
///	Routine: ComputeA, ComputeCst
///	Returns: 
///	Action : 2 roots case and a solution outside roots :
///          compute B(t,Cst) Given the cst already computed 
////////////////////////////////////////////////////
double ARM_RiccatiMSV::ARM_Funct2OutMSV::ComputeB(double t1 , double t2) const
{
	return (1/itsData->itsAt * ( log(fabs(sinh(itsData->itsAt*itsData->itsDiffRoots*(t1-itsCst))))  - log(fabs(sinh(itsData->itsAt*itsData->itsDiffRoots*(t2-itsCst))))));
}
////////////////////////////////////////////////////
///	Class  : ARM_Riccati::ARM_Funct2In
///	Routine: ComputeA, ComputeCst
///	Returns: 
///	Action : 2 roots case and a solution outside roots :
///          compute B(t,Cst) Given the cst already computed 
////////////////////////////////////////////////////
double ARM_RiccatiMSV::ARM_Funct2InMSV::ComputeB(double t1, double t2) const
{
	return (1/itsData->itsAt * ( log(fabs(cosh(itsData->itsAt*itsData->itsDiffRoots*(t1-itsCst))))  - log(fabs(cosh(itsData->itsAt*itsData->itsDiffRoots*(t2-itsCst))))));
}

////////////////////////////////////////////////////
///	Class  : ARM_RiccatiMSV
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_RiccatiMSV::CopyNoCleanUp(const ARM_RiccatiMSV& rhs)
{
	itsDtCurve			= (ARM_Curve*) ((rhs.itsDtCurve)->Clone());
	
	itsModel			= rhs.itsModel;
	itsCurveName		= rhs.itsCurveName;
	itsFloatStartTime	= rhs.itsFloatStartTime;
	itsFloatEndTime		= rhs.itsFloatEndTime;
	itsFixPayTimes		= rhs.itsFixPayTimes;
	itsFixPayPeriods	= rhs.itsFixPayPeriods;
	itsStates			= rhs.itsStates;
	itsCstShiftValue	= rhs.itsCstShiftValue;
	itsTempDt			= rhs.itsTempDt;
	itsMu				= rhs.itsMu;
}

/////////////////////////////////////////////////////
///	Class  : ARM_RiccatiMSV
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_RiccatiMSV::ARM_RiccatiMSV(const ARM_Curve& At,const ARM_Curve& Bt,const ARM_Curve& Ct,const ARM_Curve& Dt)
: ARM_Riccati(At,Bt,Ct)
{
	itsDtCurve		= (ARM_Curve*) (Dt.Clone());
	itsTempAt		= At.GetOrdinates().Elt(0);
	itsTempBt		= Bt.GetOrdinates().Elt(0);
	itsTempCt		= Ct.GetOrdinates().Elt(0);
	itsTempDt		= Dt.GetOrdinates().Elt(0);
	itsMu			= 1.0;
	itsFunctionDatasMSV.resize (itsFunctionDatas.size());
	//// Initialisation
	for(int l=0;l<itsFunctionDatas.size();++l)
		itsFunctionDatasMSV[l].itsPartIntegral = 0.0;

}

/////////////////////////////////////////////////////
///	Class  : ARM_RiccatiMSV
///	Routine: Constructor for the time-dependent Ricati equation
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_RiccatiMSV::ARM_RiccatiMSV(const ARM_Curve& At,const ARM_Curve& Bt,const ARM_Curve& Ct,const ARM_Curve& Dt,
		ARM_SVModels* model,
		string curveName,
		double floatStartTime,
		double floatEndTime,
		std::vector<double> fixPayTimes,
		std::vector<double> fixPayPeriods,
		ARM_PricingStatesPtr states,
		double cstShiftValue,
		double mu)

:	ARM_Riccati(At,Bt,Ct,false),
	itsFloatStartTime(floatStartTime),
	itsFloatEndTime(floatEndTime),
	itsFixPayTimes(fixPayTimes),
	itsFixPayPeriods(fixPayPeriods),
	itsCstShiftValue(cstShiftValue),
	itsMu(mu)
{
	itsDtCurve		= (ARM_Curve*) (Dt.Clone());
	itsTempAt		= At.GetOrdinates().Elt(0);
	itsTempBt		= Bt.GetOrdinates().Elt(0);
	itsTempCt		= Ct.GetOrdinates().Elt(0);
	itsTempDt		= Dt.GetOrdinates().Elt(0);

	itsTempAt		= 0.5 * itsTempAt * itsTempAt 	;
	itsTempCt		= - itsTempCt * itsTempCt * itsMu;


	itsStates		= states;
	itsCurveName	= curveName;
	itsModel		= model;	

	size_t schedSize = itsSchedule.size();
    itsFunctionDatas.resize(schedSize);
	itsFunctionDatasMSV.resize (schedSize);
		//// Initialisation
	for(int l=0;l<itsFunctionDatas.size();++l)
		itsFunctionDatasMSV[l].itsPartIntegral = 0.0;

}

////////////////////////////////////////////////////
///	Class  : ARM_RiccatiMSV
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_RiccatiMSV::ARM_RiccatiMSV(const ARM_RiccatiMSV& rhs)
: ARM_Riccati(rhs)
{
	itsDtCurve		= (ARM_Curve*) ((rhs.itsDtCurve)->Clone());
		
	itsTempDt			= rhs.itsTempDt;
	itsModel			= rhs.itsModel;
	itsCurveName		= rhs.itsCurveName;
	itsFloatStartTime	= rhs.itsFloatStartTime;
	itsFloatEndTime		= rhs.itsFloatEndTime;
	itsFixPayTimes		= rhs.itsFixPayTimes;
	itsFixPayPeriods	= rhs.itsFixPayPeriods;
	itsStates			= rhs.itsStates;
	itsCstShiftValue	= rhs.itsCstShiftValue;
	itsMu				= rhs.itsMu;
}


////////////////////////////////////////////////////
///	Class  : ARM_RiccatiMSV
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_RiccatiMSV::~ARM_RiccatiMSV()
{
	delete itsDtCurve;
}

////////////////////////////////////////////////////
///	Class  : ARM_RiccatiMSV
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_RiccatiMSV& ARM_RiccatiMSV::operator=(const ARM_RiccatiMSV& rhs)
{
    if( this != &rhs )
	{
		ARM_Riccati::operator=( rhs );
		/// Copy in the new ones.
		CopyNoCleanUp( rhs );
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_RiccatiMSV
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_RiccatiMSV::Clone() const
{
	return new ARM_RiccatiMSV(*this);
}
////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_RiccatiMSV
///	Routine: 
///	Returns: 
///	Action : ComputeRoots
///////////////////////////////////////////////////////////////////////////////////////
void ARM_RiccatiMSV::ComputeRoots()
{
    for(size_t i=0;i<itsSchedule.size();++i)
    {
		if(itsFunctionDatas[i].itsDelta > DELTA_LIMIT)
        {
            itsFunctionDatas[i].itsDelta=sqrt(itsFunctionDatas[i].itsDelta);
            itsFunctionDatas[i].itsNbRoots=2;
			if (itsFunctionDatas[i].itsAt>0)
			{
				itsFunctionDatas[i].itsRootInf=(-itsFunctionDatas[i].itsBt-itsFunctionDatas[i].itsDelta)/(2*itsFunctionDatas[i].itsAt);
				itsFunctionDatas[i].itsRootSup=(-itsFunctionDatas[i].itsBt+itsFunctionDatas[i].itsDelta)/(2*itsFunctionDatas[i].itsAt);
			}
			else 
			{
				itsFunctionDatas[i].itsRootSup=(-itsFunctionDatas[i].itsBt-itsFunctionDatas[i].itsDelta)/(2*itsFunctionDatas[i].itsAt);
				itsFunctionDatas[i].itsRootInf=(-itsFunctionDatas[i].itsBt+itsFunctionDatas[i].itsDelta)/(2*itsFunctionDatas[i].itsAt);
			}
			itsFunctionDatas[i].itsSumRoots =  (itsFunctionDatas[i].itsRootSup + itsFunctionDatas[i].itsRootInf)/2;
			itsFunctionDatas[i].itsDiffRoots =  (itsFunctionDatas[i].itsRootSup - itsFunctionDatas[i].itsRootInf)/2;

        }
        else if(itsFunctionDatas[i].itsDelta < -DELTA_LIMIT)
        {
            itsFunctionDatas[i].itsDelta=sqrt(-itsFunctionDatas[i].itsDelta);
            itsFunctionDatas[i].itsNbRoots=0;
            itsFunctionDatas[i].itsRootInf=0.0;
            itsFunctionDatas[i].itsRootSup=0.0;
			itsFunctionDatas[i].itsSumRoots =  -itsFunctionDatas[i].itsBt/(2*itsFunctionDatas[i].itsAt);
			itsFunctionDatas[i].itsDiffRoots =  itsFunctionDatas[i].itsDelta/(2*itsFunctionDatas[i].itsAt);

        }
        else
        {
            itsFunctionDatas[i].itsDelta=0.0;
            itsFunctionDatas[i].itsNbRoots=1;
            itsFunctionDatas[i].itsRootInf=-itsFunctionDatas[i].itsBt/(2*itsFunctionDatas[i].itsAt);
            itsFunctionDatas[i].itsRootSup=itsFunctionDatas[i].itsRootInf;
			itsFunctionDatas[i].itsSumRoots =  (itsFunctionDatas[i].itsRootSup + itsFunctionDatas[i].itsRootInf)/2;
			itsFunctionDatas[i].itsDiffRoots =  (itsFunctionDatas[i].itsRootSup - itsFunctionDatas[i].itsRootInf)/2;

        }
	}
}
////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_RiccatiMSV
///	Routine: 
///	Returns: 
///	Action : Inits function Parameters
///			WARNING: At, Bt and Ct aren't the coefficient of the Riccati Equation
///			At : Vol of Vol, Bt : Vol Mean Reversion,	Ct : Volatility
///			df/dt = 0.5 * volOfvol� * f� + VolMrs * f - Vol�* Mu
/////////////////////////////////////////////////////////////////////////////////////////
void ARM_RiccatiMSV::InitFunctionDatas()
{
    /// Store parameters for each interval ]ti, ti+1]
    /// They are stepwise right constant

    size_t schedSize = itsSchedule.size();
    itsFunctionDatas.resize(schedSize);
    double t,at,bt,ct,delta,at2;
    for(size_t i=0;i<schedSize;++i)
    {
        t=itsSchedule[i];
        itsFunctionDatas[i].itsTime=t;
        
		at = itsAtCurve->Interpolate(t); // VolOfVol
		at = at * at * 0.5;// 0.5 VolOfVol�
		itsFunctionDatas[i].itsAt=at;
		at2 = 2 * at;

		bt = itsBtCurve->Interpolate(t);
		itsFunctionDatas[i].itsBt=bt;

		ct = itsCtCurve->Interpolate(t);
		ct = - ct * ct * itsMu;					// - Mu * sigma� 
		itsFunctionDatas[i].itsCt=ct;

	   	delta = bt*bt - 2 * at2 *ct;
		itsFunctionDatas[i].itsDelta=delta;   
	}

	ComputeRoots();
    
    /// Erase saved functions if any
	ARM_FunctionsMapIter firstMapElem = itsFunctions.begin();
    ARM_FunctionsMapIter lastMapElemn = itsFunctions.end();
    itsFunctions.erase(firstMapElem,lastMapElemn);
}

////////////////////////////////////////////////////
///	Class   : ARM_RiccatiMSV
///	Routine : UpdateCstVolatility
///	Returns : string
///	Action  : Update FunctionsDatas with the new cst volatility for NR Bootstrapping
////////////////////////////////////////////////////
void ARM_RiccatiMSV::UpdateCstVolatility(double vol)
{
    for(size_t i=0;i<itsSchedule.size();++i)
    {
		itsFunctionDatas[i].
			itsCt=- vol * vol * itsMu;
		itsFunctionDatas[i].itsDelta = itsFunctionDatas[i].itsBt * itsFunctionDatas[i].itsBt - 4 * itsFunctionDatas[i].itsAt * itsFunctionDatas[i].itsCt ;
	}
	ComputeRoots();
	
	ARM_FunctionsMapIter firstMapElem = itsFunctions.begin();
    ARM_FunctionsMapIter lastMapElemn = itsFunctions.end();
    itsFunctions.erase(firstMapElem,lastMapElemn);
}

////////////////////////////////////////////////////
///	Class   : ARM_RiccatiMSV
///	Routine : derivs
///	Returns : string
///	Action  : Function to Implement for Runge Kutta Resolution
///			We have to Set At , Bt and Ct at each time Step
////////////////////////////////////////////////////
void ARM_RiccatiMSV::derivs(double x, std::vector<double>& yt, std::vector<double>& dyt) const
{
		UpdateTempParams(x);
		(*dyt)[0] = itsTempAt*(*yt)[0]*(*yt)[0] + itsTempBt*(*yt)[0] + itsTempCt;
		(*dyt)[1] = itsTempDt*(*yt)[0];
}


////////////////////////////////////////////////////
///	Class  : ARM_RiccatiMSV
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_RiccatiMSV::UpdateTempParams(double x) const
{
	/// Only Ct is time-dependent and has to be updated
	//itsTempAt = itsModel->UpdateRiccatiAt(x);
	itsTempCt = itsModel->UpdateRiccatiCt(itsCurveName,x,itsFloatStartTime,itsFloatEndTime,itsFixPayTimes,itsFixPayPeriods);
	itsTempCt = - itsMu * itsTempCt *itsTempCt ;
}

////////////////////////////////////////////////////
///	Class  : ARM_RiccatiMSV
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_RiccatiMSV::ComputeAdditionalDatas()
{
}

////////////////////////////////////////////////////
///	Class  : ARM_RiccatiMSV
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_Riccati::ARM_FunctionsMapIter ARM_RiccatiMSV::GenerateFunction(double maturity) 
{
    /// Localise time just before maturity
    size_t nbFunc=itsFunctionDatas.size();
    int iFirst;
    for(iFirst=0;iFirst<nbFunc && itsFunctionDatas[iFirst].itsTime + K_NEW_DOUBLE_TOL < maturity*K_YEAR_LEN;++iFirst);
    if(iFirst>=nbFunc) iFirst=nbFunc-1;

    double initCond;
    vector< ARM_FunctTypePtr > function(iFirst+1);
    double T0=maturity;
    for(int i=iFirst;i>=0;--i)
    {
        /// Compute the initial condition
        if(i==iFirst)
            initCond = 0.0;
        else
        {
            T0 = itsFunctionDatas[i].itsTime / K_YEAR_LEN;
            initCond = ((ARM_FunctTypeRiccati*) (&*function[i+1]))->ComputeA(T0);
        }

        if(itsFunctionDatas[i].itsNbRoots == 2)
        {
            /// Affect the right function depending of the initial condition
            if(itsFunctionDatas[i].itsRootInf < initCond && initCond < itsFunctionDatas[i].itsRootSup)
			{
				ARM_Funct2InMSV* func = new ARM_Funct2InMSV(&(itsFunctionDatas[i]));
				/// Compute integration constant
		        func->InitCst(T0,initCond);

				if(i==0)
					itsFunctionDatasMSV[i].itsPartIntegral= (func->ComputeB(0,T0) + itsFunctionDatas[i].itsSumRoots * T0) * itsFunctionDatas[i].itsBt;
                else
					itsFunctionDatasMSV[i].itsPartIntegral= (func->ComputeB(itsFunctionDatas[i-1].itsTime/ K_YEAR_LEN,T0)+ itsFunctionDatas[i].itsSumRoots * (T0 - itsFunctionDatas[i-1].itsTime/ K_YEAR_LEN))*itsFunctionDatas[i].itsBt;

				function[i] = ARM_FunctTypePtr(func);
			}
			else if(initCond < itsFunctionDatas[i].itsRootInf || itsFunctionDatas[i].itsRootSup < initCond)
			{
				ARM_Funct2OutMSV* func = new ARM_Funct2OutMSV(&(itsFunctionDatas[i]));
				/// Compute integration constant
				func->InitCst(T0,initCond);
				if(i==0)
					itsFunctionDatasMSV[i].itsPartIntegral= (func->ComputeB(0,T0)+ itsFunctionDatas[i].itsSumRoots * T0) * itsFunctionDatas[i].itsBt;
                else
					itsFunctionDatasMSV[i].itsPartIntegral= (func->ComputeB(itsFunctionDatas[i-1].itsTime/ K_YEAR_LEN,T0)+ itsFunctionDatas[i].itsSumRoots * (T0 - itsFunctionDatas[i-1].itsTime/ K_YEAR_LEN)) * itsFunctionDatas[i].itsBt;

				function[i] = ARM_FunctTypePtr(func);
			}
            else if(itsFunctionDatas[i].itsRootInf == initCond)
			{
				ARM_Funct2InfMSV* func = new ARM_Funct2InfMSV(&(itsFunctionDatas[i]));
			    /// Compute integration constant
				func->InitCst(T0,initCond);
				function[i] = ARM_FunctTypePtr(func);
			}
            else
			{
                ARM_Funct2SupMSV* func = new ARM_Funct2SupMSV(&(itsFunctionDatas[i]));
				/// Compute integration constant
				func->InitCst(T0,initCond);
				function[i] = ARM_FunctTypePtr(func);
			}
        }
        else if(itsFunctionDatas[i].itsNbRoots == 1)
        {
			ARM_Funct1MSV* func = new ARM_Funct1MSV(&(itsFunctionDatas[i]));
			func->InitCst(T0,initCond);
			function[i] = ARM_FunctTypePtr(func);
		}
        else
		{
            ARM_Funct0MSV* func = new ARM_Funct0MSV(&(itsFunctionDatas[i]));
			func->InitCst(T0,initCond);
			function[i] = ARM_FunctTypePtr(func);
		}

    }

    pair< double,vector< ARM_FunctTypePtr > > value(maturity,function);
    pair< ARM_FunctionsMapIter,bool > result = itsFunctions.insert(value);
    if(!result.second)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": can't insert a new Riccati solution in the map");

    return result.first;
}

////////////////////////////////////////////////////
///	Class  : ARM_RiccatiMSV
///	Routine: 
///	Returns: 
///	Action : computes B(0,T)
/// WARNING We suppose that we have already generated the function
////////////////////////////////////////////////////
double ARM_RiccatiMSV::B0T(double T)
{
	ARM_FunctionsMap::iterator found = itsFunctions.find(T);
    if(found == itsFunctions.end())
        /// Insert the new function in the map
        found = GenerateFunction(T);


    double B0T=0.;
	for (int i=0;i<itsFunctionDatas.size();++i)
		B0T += itsFunctionDatasMSV[i].itsPartIntegral;

#ifdef COMPUTATION_TIME_TEST
    itsFunctions.erase(found);
#endif

	return (-B0T);
}




CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/