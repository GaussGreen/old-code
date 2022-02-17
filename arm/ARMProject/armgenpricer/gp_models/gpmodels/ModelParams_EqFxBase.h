
#ifndef _INGPMODELS_MODELPARAMS_EQFXBASE_H
#define _INGPMODELS_MODELPARAMS_EQFXBASE_H

/// gpbase
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "gpbase/curve.h"
#include "gpbase/curveutils.h"
#include "gpbase/ostringstream.h"
#include "gpbase/utilityport.h"
#include "gpbase/cloneutilityfunc.h"

/// gpinfra
#include "gpinfra/curvemodelparam.h"
#include "typedef.h"

/// forward declaration
class ARM_ZeroCurve;

CC_BEGIN_NAMESPACE( ARM )


//-----------------------------------------------------------------------------
// \class ARM_ModelParams_EqFxBase
//-----------------------------------------------------------------------------

class ARM_ModelParams_EqFxBase
{
private:
	double itsSpot;
public:
	ARM_ModelParams_EqFxBase(double spot)
	:	itsSpot(spot) {};
	
	ARM_ModelParams_EqFxBase( const ARM_ModelParams_EqFxBase& rhs );

    ASSIGN_OPERATOR(ARM_ModelParams_EqFxBase)

	inline double GetSpot() const { return itsSpot; }
	inline void SetSpot(double s)  {itsSpot=s;}

	virtual string toString(const string& indent="",const string& nextIndent="") const;
};


///-----------------------------------------------------------------------------
/// \class ARM_ModelParams_Eq
/// \brief Class where the fwd are specified by a curve of dividend
///-----------------------------------------------------------------------------


class ARM_ModelParams_Eq :public ARM_ModelParams_EqFxBase
{
private:
	ARM_ZeroCurvePtr itsZeroCurve; /// shared (not very robust but inherited from ARM Kernel!)

public:
	ARM_ModelParams_Eq( const ARM_ZeroCurvePtr& zeroCurve, double spot );

	ARM_ModelParams_Eq( const ARM_ModelParams_Eq& rhs )
	:	ARM_ModelParams_EqFxBase( rhs), itsZeroCurve(rhs.itsZeroCurve) {};

    ASSIGN_OPERATOR(ARM_ModelParams_Eq)
	virtual ~ARM_ModelParams_Eq(){}

	virtual double Forward( double MaturityTime ) const;

	virtual string toString(const string& indent="",const string& nextIndent="") const;
};

///-----------------------------------------------------------------------------
/// \class ARM_ModelParams_Fx
/// \brief Class where the fwd are specified by the two interest rates curve
///-----------------------------------------------------------------------------

class ARM_ModelParams_Fx :	public ARM_ModelParams_EqFxBase
{
private:
	ARM_ZeroCurvePtr itsDomCurve; /// shared (not very robust but inherited from ARM Kernel!)
	ARM_ZeroCurvePtr itsForCurve; /// shared (not very robust but inherited from ARM Kernel!)
	string itsFXCal;
	double itsSpotDays;
	double itsDfRatioAtSpotDate;

public:
	ARM_ModelParams_Fx( const ARM_ZeroCurvePtr& domCurve = ARM_ZeroCurvePtr(NULL), 
		const ARM_ZeroCurvePtr& forCurve = ARM_ZeroCurvePtr(NULL), 
		double spot= 1.0 );

	ARM_ModelParams_Fx( const ARM_ModelParams_Fx& rhs )
	:	ARM_ModelParams_EqFxBase(rhs), 
//		itsDomCurve(CreateClonedPtr(&*rhs.itsDomCurve)), 
	//	itsForCurve(CreateClonedPtr(&*rhs.itsForCurve)), 
		itsFXCal( rhs.itsFXCal ),
		itsSpotDays( rhs.itsSpotDays ),
		itsDfRatioAtSpotDate( rhs.itsDfRatioAtSpotDate )
	{}
    
	ASSIGN_OPERATOR(ARM_ModelParams_Fx)
	
	virtual ~ARM_ModelParams_Fx(){}

	virtual double Forward( double MaturityTime ) const;

	/// static version of forward() to avoid building an object
	static double Forward(ARM_ZeroCurve* domCurve,ARM_ZeroCurve* forCurve,double spot,double MaturityTime);

	static ARM_Date ComputeSettlementDate(ARM_ZeroCurve* domCurve,ARM_ZeroCurve* forCurve,ARM_Date& date);

	/// accessors
	inline ARM_ZeroCurvePtr GetDomCurve() const { return itsDomCurve; }
	inline ARM_ZeroCurvePtr GetForCurve() const { return itsForCurve; }
	
	inline void SetDomCurve( const ARM_ZeroCurvePtr& domCrv ) { itsDomCurve = domCrv; }
	inline void SetForCurve( const ARM_ZeroCurvePtr& forCrv ) { itsForCurve = forCrv; }

	virtual string toString(const string& indent="",const string& nextIndent="") const;
};

template<typename T>
class ARM_ModelParams_Eq_T : public ARM_ModelParams_Eq, public T
{
public:
	ARM_ModelParams_Eq_T( const ARM_ModelParamVector& params=ARM_ModelParamVector(), 
		const ARM_ZeroCurvePtr& zeroCurve = ARM_ZeroCurvePtr(NULL), 
		double spot=1.0 ) 
	
	:ARM_ModelParams_Eq(zeroCurve,spot),
		T(params)
	{
	}

	ARM_ModelParams_Eq_T( const ARM_ModelParams_Eq_T& rhs ) 
	:
		ARM_ModelParams_Eq(rhs),
		T(rhs)
	{
	}
	ASSIGN_OPERATOR(ARM_ModelParams_Eq_T)
	virtual ~ARM_ModelParams_Eq_T()
	{
	}
    
	/// How many factors?
    virtual size_t FactorCount() const { return 1; };
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {}
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) {};

	/// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_ModelParams_Eq_T(*this); };
	virtual string toString(const string& indent="",const string& nextIndent="") const {
		string str;
		str += ARM_ModelParams_Eq::toString(indent, nextIndent); 
		str += T::toString(indent, nextIndent);
		return str; 
	}
};

template<typename T>
class ARM_ModelParams_Fx_T : public ARM_ModelParams_Fx, public T
{
public:
	ARM_ModelParams_Fx_T( const ARM_ModelParamVector& params=ARM_ModelParamVector(), 
		ARM_ZeroCurvePtr domCurve=ARM_ZeroCurvePtr(NULL), 
		ARM_ZeroCurvePtr fgnCurve=ARM_ZeroCurvePtr(NULL), 
		double spot=1.0 ) 
		
	:ARM_ModelParams_Fx(domCurve,fgnCurve,spot),
		T(params)
	{
	}
	ARM_ModelParams_Fx_T( const ARM_ModelParams_Fx_T& rhs ) :
		ARM_ModelParams_Fx(rhs),
		T(rhs)
	{
	}
	ASSIGN_OPERATOR(ARM_ModelParams_Fx_T)
	virtual ~ARM_ModelParams_Fx_T()
	{
	}

	/// How many factors?
    virtual size_t FactorCount() const { return 1; };
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {}
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) {};

	/// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_ModelParams_Fx_T(*this); };
    virtual string toString(const string& indent="",const string& nextIndent="") const {
		string str;
		str += ARM_ModelParams_Fx::toString(indent, nextIndent); 
		str += T::toString(indent, nextIndent);
		return str; 
	};
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

