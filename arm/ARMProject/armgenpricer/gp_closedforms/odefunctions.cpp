/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  function solution of the Riccati Equation with Piece-Wise functions A(t), B(t) and C(t)
 *	df/dt = A(t) * f(t)² + B(t) * f(t) + C(t)
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
#include "gpclosedforms/riccati.h"
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
///          solve Cst so that A(t,Cst)=a,
///          compute A(t,Cst),
////////////////////////////////////////////////////
void ARM_Riccati::ARM_Funct0::InitCst(double t,double a)
{
	double temp = 1/(2 * itsData->itsAt * itsData->itsDiffRoots);
    itsCst = t-temp*atan((a-itsData->itsSumRoots)/(itsData->itsDiffRoots));
}

double ARM_Riccati::ARM_Funct0::ComputeA(double t) const
{
    return (itsData->itsSumRoots + itsData->itsDiffRoots*tan(2*itsData->itsDiffRoots*itsData->itsAt*(t-itsCst)));
}


////////////////////////////////////////////////////
///	Class  : ARM_Riccati::ARM_Funct1
///	Routine: ComputeA, ComputeCst
///	Returns: 
///	Action : 1 root case :
///          solve Cst so that A(t,Cst)=a,
///          compute A(t,Cst),
////////////////////////////////////////////////////
void ARM_Riccati::ARM_Funct1::InitCst(double t,double a)
{
	itsCst=t - 1.0/(itsData->itsAt*(itsData->itsRootInf-a));
}

double ARM_Riccati::ARM_Funct1::ComputeA(double t) const
{
    t /= K_YEAR_LEN;
    return (itsData->itsRootInf - 1/(itsData->itsAt*(t-itsCst)));
}

////////////////////////////////////////////////////
///	Class  : ARM_Riccati::ARM_Funct2Sup
///	Routine: ComputeA, ComputeCst
///	Returns: 
///	Action : 2 roots case and a constant solution
///          equals to the greater root :
///          solve Cst so that A(t,Cst)=a,
///          compute A(t,Cst),
////////////////////////////////////////////////////
void ARM_Riccati::ARM_Funct2Sup::InitCst(double t,double a)
{
	itsCst = itsData->itsRootSup;
}

double ARM_Riccati::ARM_Funct2Sup::ComputeA(double t) const
{
    return itsData->itsRootSup;
}

////////////////////////////////////////////////////
///	Class  : ARM_Riccati::ARM_Funct2Inf
///	Routine: ComputeA, ComputeCst
///	Returns: 
///	Action : 2 roots case and a constant solution
///          equals to the lower root :
///          solve Cst so that A(t,Cst)=a,
///          compute A(t,Cst),
////////////////////////////////////////////////////
void ARM_Riccati::ARM_Funct2Inf::InitCst(double t,double a)
{
    itsCst = itsData->itsRootInf;
}

double ARM_Riccati::ARM_Funct2Inf::ComputeA(double t) const
{
    return itsData->itsRootInf;
}


////////////////////////////////////////////////////
///	Class  : ARM_Riccati::ARM_Funct2Out
///	Routine: ComputeA, ComputeCst
///	Returns: 
///	Action : 2 roots case and a solution outside roots :
///          solve Cst so that A(t,Cst)=a,
///          compute A(t,Cst),
////////////////////////////////////////////////////
void ARM_Riccati::ARM_Funct2Out::InitCst(double t,double a)
{
	double temp = 1/(2 * itsData->itsAt * itsData->itsDiffRoots);
    itsCst = t-temp*acoth((itsData->itsSumRoots - a)/(itsData->itsDiffRoots));
}

double ARM_Riccati::ARM_Funct2Out::ComputeA(double t) const
{
    return (itsData->itsSumRoots - itsData->itsDiffRoots / tanh(2 * itsData->itsAt * itsData->itsDiffRoots * (t - itsCst))); 
}
////////////////////////////////////////////////////
///	Class  : ARM_Riccati::ARM_Funct2In
///	Routine: ComputeA, ComputeCst
///	Returns: 
///	Action : 2 roots case and a solution outside roots :
///          solve Cst so that A(t,Cst)=a,
///          compute A(t,Cst),
////////////////////////////////////////////////////
void ARM_Riccati::ARM_Funct2In::InitCst(double t,double a)
{
	double temp = 1/(2 * itsData->itsAt * itsData->itsDiffRoots);
    itsCst = t-temp*atanh(-(itsData->itsSumRoots - a)/(itsData->itsDiffRoots));
}

double ARM_Riccati::ARM_Funct2In::ComputeA(double t) const
{
    return (itsData->itsSumRoots + itsData->itsDiffRoots * tanh(2 * itsData->itsAt * itsData->itsDiffRoots * (t - itsCst))); 
}

////////////////////////////////////////////////////
///	Class  : ARM_Riccati
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_Riccati::CopyNoCleanUp(const ARM_Riccati& rhs)
{
    itsFunctionDatas    = rhs.itsFunctionDatas;
    itsFunctions        = rhs.itsFunctions;

	itsAtCurve		= (ARM_Curve*) ((rhs.itsAtCurve)->Clone());
	itsBtCurve		= (ARM_Curve*) ((rhs.itsBtCurve)->Clone());
	itsCtCurve		= (ARM_Curve*) ((rhs.itsCtCurve)->Clone());
}

////////////////////////////////////////////////////
///	Class  : ARM_Riccati
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_Riccati::ARM_Riccati(const ARM_Curve& At,const ARM_Curve& Bt,const ARM_Curve& Ct)
{
	itsAtCurve		= (ARM_Curve*) (At.Clone());
	itsBtCurve		= (ARM_Curve*) (Bt.Clone());
	itsCtCurve		= (ARM_Curve*) (Ct.Clone());

	GenerateSchedule();
    InitFunctionDatas(itsSchedule);
}


////////////////////////////////////////////////////
///	Class  : ARM_Riccati
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_Riccati::ARM_Riccati(const ARM_Riccati& rhs)
: ARM_RootObject(rhs)
{
    CopyNoCleanUp(rhs);
}


////////////////////////////////////////////////////
///	Class  : ARM_Riccati
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_Riccati::~ARM_Riccati()
{
	delete itsAtCurve;
	delete itsBtCurve;
	delete itsCtCurve;
}

////////////////////////////////////////////////////
///	Class  : ARM_Riccati
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_Riccati& ARM_Riccati::operator=(const ARM_Riccati& rhs)
{
    if( this != &rhs )
	{
		ARM_RootObject::operator=( rhs );		
		this->~ARM_Riccati();
		new (this) ARM_Riccati(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_Riccati
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_Riccati::Clone() const
{
	return new ARM_Riccati(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_Riccati
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_Riccati::GenerateSchedule()
{
	/// Compute a merged schedule of model parameters
    const ARM_GP_Vector& schedAt  = itsAtCurve->GetAbscisses();
    const ARM_GP_Vector& schedBt  = itsBtCurve->GetAbscisses();
    const ARM_GP_Vector& schedCt = itsCtCurve->GetAbscisses();
    ARM_GP_Vector tmpSched(schedAt.size()+schedBt.size());
	itsSchedule = ARM_GP_Vector(tmpSched.size()+schedCt.size());
    CC_NS(std,merge)(schedAt.begin(),schedAt.end(),schedBt.begin(),schedBt.end(),tmpSched.begin());
    CC_NS(std,merge)(tmpSched.begin(),tmpSched.end(),schedCt.begin(),schedCt.end(),itsSchedule.begin());
    ARM_GP_Vector::iterator last=CC_NS(std,unique)(itsSchedule.begin(),itsSchedule.end());
	itsSchedule.resize( last-itsSchedule.begin() );
    if(itsSchedule[0] < K_NEW_DOUBLE_TOL && itsSchedule.size() > 1)
        itsSchedule.erase(itsSchedule.begin());
}

////////////////////////////////////////////////////
///	Class  : ARM_Riccati
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_Riccati::InitFunctionDatas(const ARM_GP_Vector& sched)
{
    /// Store parameters for each interval ]ti, ti+1]
    /// They are stepwise right constant
    size_t schedSize = sched.size();
    itsFunctionDatas.resize(schedSize);
    double t,at,bt,ct,delta,at2;
    for(size_t i=0;i<schedSize;++i)
    {
        t=sched[i];
        itsFunctionDatas[i].itsTime=t;
        
		at = itsAtCurve->Interpolate(t);
		itsFunctionDatas[i].itsAt=at;

		at2 = 2 * at;

		bt = itsBtCurve->Interpolate(t);
		itsFunctionDatas[i].itsBt=bt;

		ct = itsCtCurve->Interpolate(t);
		itsFunctionDatas[i].itsCt=ct;

	   	delta = bt*bt - 2 * at2 *ct;
        if(delta > DELTA_LIMIT)
        {
            delta=sqrt(delta);
            itsFunctionDatas[i].itsNbRoots=2;
            itsFunctionDatas[i].itsRootInf=(-bt-delta)/(at2);
            itsFunctionDatas[i].itsRootSup=(-bt+delta)/(at2);
        }
        else if(delta < -DELTA_LIMIT)
        {
            delta=sqrt(-delta);
            itsFunctionDatas[i].itsNbRoots=0;
            itsFunctionDatas[i].itsRootInf=0.0;
            itsFunctionDatas[i].itsRootInf=0.0;
        }
        else
        {
            delta=0.0;
            itsFunctionDatas[i].itsNbRoots=1;
            itsFunctionDatas[i].itsRootInf=-bt/(at2);
            itsFunctionDatas[i].itsRootSup=itsFunctionDatas[i].itsRootInf;
        }
        itsFunctionDatas[i].itsDelta=delta;   
		itsFunctionDatas[i].itsSumRoots =  (itsFunctionDatas[i].itsRootSup + itsFunctionDatas[i].itsRootInf)/2;
		itsFunctionDatas[i].itsDiffRoots =  (itsFunctionDatas[i].itsRootSup - itsFunctionDatas[i].itsRootInf)/2;

    }

    /// Erase saved functions if any
    ARM_FunctionsMapIter firstMapElem = itsFunctions.begin();
    ARM_FunctionsMapIter lastMapElemn = itsFunctions.end();
    itsFunctions.erase(firstMapElem,lastMapElemn);
}


////////////////////////////////////////////////////
///	Class  : ARM_Riccati
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_Riccati::ARM_FunctionsMapIter ARM_Riccati::GenerateFunction(double maturity) const
{
    /// Localise time just before maturity
    size_t nbFunc=itsFunctionDatas.size();
    int iFirst;
    for(iFirst=0;iFirst<nbFunc && itsFunctionDatas[iFirst].itsTime + K_NEW_DOUBLE_TOL < maturity;++iFirst);
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
            T0 = itsFunctionDatas[i].itsTime;
            initCond = ((ARM_FunctTypeRiccati*) (&*function[i+1]))->ComputeA(T0);
        }

        if(itsFunctionDatas[i].itsNbRoots == 2)
        {
            /// Affect the right function depending of the initial condition
            if(itsFunctionDatas[i].itsRootInf < initCond && initCond < itsFunctionDatas[i].itsRootSup)
                function[i] = ARM_FunctTypePtr((ARM_FunctType*) new ARM_Funct2In(&(itsFunctionDatas[i])));
            else if(initCond < itsFunctionDatas[i].itsRootInf || itsFunctionDatas[i].itsRootSup < initCond)
                function[i] = ARM_FunctTypePtr((ARM_FunctType*) new ARM_Funct2Out(&(itsFunctionDatas[i])));
            else if(itsFunctionDatas[i].itsRootInf == initCond)
                function[i] = ARM_FunctTypePtr((ARM_FunctType*) new ARM_Funct2Inf(&(itsFunctionDatas[i])));
            else
                function[i] = ARM_FunctTypePtr((ARM_FunctType*) new ARM_Funct2Sup(&(itsFunctionDatas[i])));
        }
        else if(itsFunctionDatas[i].itsNbRoots == 1)
            function[i] = ARM_FunctTypePtr((ARM_FunctType*) new ARM_Funct1(&(itsFunctionDatas[i])));
        else
            function[i] = ARM_FunctTypePtr((ARM_FunctType*) new ARM_Funct0(&(itsFunctionDatas[i])));

        /// Compute integration constant
        ((ARM_FunctTypeRiccati*) (&*function[i]))->InitCst(T0,initCond);
    }

    pair< double,vector< ARM_FunctTypePtr > > value(maturity,function);
    pair< ARM_FunctionsMapIter,bool > result = itsFunctions.insert(value);
    if(!result.second)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": can't insert a new Riccati solution in the map");

    return result.first;
}

////////////////////////////////////////////////////
///	Class  : ARM_Riccati
///	Routine: ComputeA
///	Returns: double
///	Action : Compute A(t,T) using the input function
///          A(.,T).
////////////////////////////////////////////////////
double ARM_Riccati::ComputeA(double time,const vector< ARM_FunctTypePtr >& function) const
{
    /// Localise index w.r.t. schedule
    size_t nbFunc=itsFunctionDatas.size();
    int iFirst;
    for(iFirst=0;iFirst<nbFunc && itsFunctionDatas[iFirst].itsTime + K_NEW_DOUBLE_TOL < time;++iFirst);
    if(iFirst>=nbFunc) iFirst=nbFunc-1;

    /// Call the right function
    ARM_FunctTypePtr funct=function[iFirst];
    if(funct == ARM_FunctTypePtr(NULL))
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": no function for this interval");

    return ((ARM_FunctTypeRiccati*) (&*funct))->ComputeA(time);
}



double ARM_Riccati::A(double t, double T) const
{
	ARM_FunctionsMap::iterator found = itsFunctions.find(T);
    if(found == itsFunctions.end())
        /// Insert the new function in the map
        found = GenerateFunction(T);

    double AtT=ComputeA(t,found->second);

#ifdef COMPUTATION_TIME_TEST
    itsFunctions.erase(found);
#endif

    return AtT;
}


////////////////////////////////////////////////////
///	Class   : ARM_Riccati
///	Routine : toString
///	Returns : string
///	Action  : object dump into a string
////////////////////////////////////////////////////
string ARM_Riccati::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << indent << "Riccati Equation\n";
    os << indent << "-----------\n";

    os << indent << "Function Datas\n";
    os << indent << setw(5)     << "At";
    os << indent << setw(10)    << "Bt";
    os << indent << setw(10)    << "Ct";
    os << indent << setw(10)    << "Delta";
    os << indent << setw(10)    << "NbR";
    os << indent << setw(10)    << "RInf";
    os << indent << setw(10)    << "Rsup";
    os << endl;
    for(int i=0; i<itsFunctionDatas.size();++i)
    {
        os << indent << fixed << setw(5) << setprecision(0) << itsFunctionDatas[i].itsAt;
        os << indent << fixed << setw(10) << setprecision(4) << itsFunctionDatas[i].itsBt;
        os << indent << fixed << setw(10) << setprecision(6) << itsFunctionDatas[i].itsCt;
        os << indent << fixed << setw(10) << setprecision(6) << itsFunctionDatas[i].itsDelta;
        os << indent << dec << setw(10) << itsFunctionDatas[i].itsNbRoots;
        os << indent << fixed << setw(10) << setprecision(3) << itsFunctionDatas[i].itsRootInf;
        os << indent << fixed << setw(10) << setprecision(3) << itsFunctionDatas[i].itsRootSup;
        os << endl;
    }

    return os.str();
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/