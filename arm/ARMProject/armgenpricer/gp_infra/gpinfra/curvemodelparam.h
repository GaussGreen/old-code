#ifndef _INGPINFRA_CURVEMODELPARAM_H
#define _INGPINFRA_CURVEMODELPARAM_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/numericconstant.h"
#include <string>
#include "typedef.h"
#include "modelparam.h"
#include <cmath>				/// for fabs


CC_BEGIN_NAMESPACE( ARM )


///////////////////////////////////////////////////
/// \class ARM_CurveModelParam
/// \brief
/// ARM_CurveModelParam is a model param with a curve
/////////////////////////////////////////////////////
class ARM_CurveModelParam : public ARM_ModelParam
{
private:
    /// Curve to manage correctly model parmeters
    ARM_Curve* itsCurve;
	string itsCurrency;
	string itsName;

	ARM_GP_Vector* itsLowerBound;
    ARM_GP_Vector* itsUpperBound;
	/// inital curve to store the value  
    ARM_Curve* itsInitialCurve; 
	
	void CopyNoCleanUp( const ARM_CurveModelParam& rhs );
    void CleanUp();

public:
	/// dimension
	int ARM_ModelParamDimension() const
	{ 
		return 1;
	}
	/// constructor (default and with argument)
	ARM_CurveModelParam();
	ARM_CurveModelParam( ParamType type,
		ARM_GP_Vector* values,
		ARM_GP_Vector* breakPointTimes	, 
		const string& paramName			= "",
        const string& interpolatorName	= "STEPUPRIGHT",
		ARM_GP_Vector* lowerBound		= NULL,
		ARM_GP_Vector* upperBound		= NULL,
		bool adviseBreakPointTimes		= false,
		const string & currency			= "");

	/// constructor mono-dim
	ARM_CurveModelParam( ParamType type,
		double value,
		const string& paramName		= "",
		double lowerBound			= -ARM_NumericConstants::ARM_BIGGEST_POSITIVE_NUMBER,
		double upperBound			= ARM_NumericConstants::ARM_BIGGEST_POSITIVE_NUMBER,
		bool adviseBreakPointTimes	= false,
		const string &	currency	= "");


	/// copy constructor
    ARM_CurveModelParam( const ARM_CurveModelParam& rhs );
	ARM_CurveModelParam& operator=( const ARM_CurveModelParam& rhs );
	virtual ~ARM_CurveModelParam();

	/// downcast
	virtual ARM_CurveModelParam& ToCurveModelParam();
	virtual const ARM_CurveModelParam& ToCurveModelParam() const;

	/// std virtual function of a model param
	virtual double GetValue( double x1 ) const;
	virtual double GetValue( double x1, double x2 ) const;

	virtual double GetTimeAtPoint( size_t i ) const	;
	virtual double GetTenorAtPoint( size_t i ) const;
	virtual double GetTimeAtPoint( size_t i, size_t j ) const;
	virtual double GetTenorAtPoint( size_t i, size_t j ) const;

	virtual double GetValueAtPoint( size_t i ) const;
	virtual double GetValueAtPoint( size_t i, size_t j ) const;

	virtual void SetValue( double x1, double value );
	virtual void SetValue( double x1, double x2, double value  );
	virtual void SetValueAtPoint( size_t i, double value );
	virtual void SetValueAtPoint( size_t i, size_t j, double value );
    virtual void MergeModelParam( ARM_ModelParam* NewValue );
	virtual size_t size() const;
	virtual size_t rows() const;
	virtual size_t cols() const;

	virtual double GetValue( double x1, double x2, double x3 ) const;
	virtual double GetValueAtPoint( size_t i, size_t j, size_t k) const;
	virtual void SetValue(double x1, double x2, double x3, double value);
	virtual void SetValueAtPoint( size_t i, size_t j, size_t k, double value);

	///	curve specific functions!
	ARM_GP_Vector* GetData( ARM_DataType type ) const;
    inline ARM_Curve* GetCurve() const {return itsCurve;}; 
    void SetCurve(ARM_Curve* curve);

    /// Replace old values using interpolate if necessairy
    virtual void SetAndUpDate(ARM_GP_Vector* values);    
    virtual void SetValuesAndTimes(ARM_GP_Vector* breakPointTimes,ARM_GP_Vector* values);
	/// method to update the its values using interpolate method
    virtual void UpdateValues(ARM_GP_Vector* BreakPointTimes);

	/// standard ARM Object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="", const string& nextIndent="") const;


	ARM_GP_Vector* GetLowerBound(){return itsLowerBound;};
    void SetLowerBound(ARM_GP_Vector* lowerbound){itsLowerBound = lowerbound;};
    
    ARM_GP_Vector* GetUpperBound(){return itsUpperBound;}; 
    void SetUpperBound(ARM_GP_Vector* upperbound){itsUpperBound = upperbound;}; 

	/// initial curve and value
    inline ARM_Curve* GetInitialCurve() const {return itsInitialCurve;} 
	void SetToInitialCurve();

	ARM_CurveModelParam* GetCalibParam( size_t begin, size_t end);

	virtual void BumpModelParam(size_t row = 0, size_t column = 0,
		double shift = 0.0,ARM_BumpParamType bumpType = ARM_ModelParamBump::isNothing); 

	/// Get/Set currency
	
	inline void SetCurrency( const string & currency){ itsCurrency= currency;}
	inline string GetCurrency() const { return itsCurrency; }


};


///////////////////////////////////////////////////
/// \class ARM_RangedCurveModelParam
/// \brief
/// ARM_RangedCurveModelParam prevents from using certain values
/////////////////////////////////////////////////////

/// positive curve model param
class ARM_RangedCurveModelParam : public ARM_CurveModelParam
{
public:
	ARM_RangedCurveModelParam( ParamType type,
		ARM_GP_Vector* values,
		ARM_GP_Vector* breakPointTimes	, 
		const string& paramName			= "",
        const string& interpolatorName	= "STEPUPRIGHT" );

	/// copy constructor
    ARM_RangedCurveModelParam( const ARM_RangedCurveModelParam& rhs )
	:	ARM_CurveModelParam(rhs) {}
	ARM_RangedCurveModelParam& operator=( const ARM_RangedCurveModelParam& rhs )
	{
		if( this != &rhs )
			ARM_CurveModelParam::operator =(rhs);
		return *this;
	}
	virtual ~ARM_RangedCurveModelParam() {};

	virtual double ModifiedValue( double value ) const { return value; }

	/// set function to prevent a negative value
	virtual void SetValue( double x1, double value )
	{ ARM_CurveModelParam::SetValue( x1, ModifiedValue(value) ); }
	virtual void SetValue( double x1, double x2, double value  )
	{ ARM_CurveModelParam::SetValue( x1, x2, ModifiedValue(value) ); }
	virtual void SetValueAtPoint( size_t i, double value )
	{ ARM_CurveModelParam::SetValueAtPoint( i, ModifiedValue(value) ); }
	virtual void SetValueAtPoint( size_t i, size_t j, double value )
	{ ARM_CurveModelParam::SetValueAtPoint( i, j, ModifiedValue(value) ); }

    virtual void SetValuesAndTimes(ARM_GP_Vector* breakPointTimes,ARM_GP_Vector* values);
};


///////////////////////////////////////////////////
/// \class ARM_PositiveCurveModelParam
/// \brief positive curve model param
/////////////////////////////////////////////////////

class ARM_PositiveCurveModelParam: public ARM_RangedCurveModelParam
{
public:
	/// constructor, copy constructor, operator=,destructor
	ARM_PositiveCurveModelParam( ParamType type,
		ARM_GP_Vector* values,
		ARM_GP_Vector* breakPointTimes	, 
		const string& paramName			= "",
        const string& interpolatorName	= "STEPUPRIGHT" )
	:	ARM_RangedCurveModelParam(type,values,breakPointTimes,paramName,interpolatorName ) {}

	ARM_PositiveCurveModelParam( const ARM_PositiveCurveModelParam& rhs )
	:	ARM_RangedCurveModelParam(rhs) {}
	ARM_PositiveCurveModelParam& operator=( const ARM_PositiveCurveModelParam& rhs )
	{
		if( this != &rhs )
			ARM_RangedCurveModelParam::operator =(rhs);
		return *this;
	}
	virtual ~ARM_PositiveCurveModelParam(){};

	/// function to transform model param
	virtual double ModifiedValue( double value ) const
	{ return fabs(value); }
};


///////////////////////////////////////////////////
/// \class ARM_PositiveTruncatedCurveModelParam
/// \brief positive truncated curve model param
/////////////////////////////////////////////////////

class ARM_PositiveTruncatedCurveModelParam: public ARM_RangedCurveModelParam
{
public:
	/// constructor, copy constructor, operator=,destructor
	ARM_PositiveTruncatedCurveModelParam( ParamType type,
		ARM_GP_Vector* values,
		ARM_GP_Vector* breakPointTimes	, 
		const string& paramName			= "",
        const string& interpolatorName	= "STEPUPRIGHT" )
	:	ARM_RangedCurveModelParam(type,values,breakPointTimes,paramName,interpolatorName ) {}

    ARM_PositiveTruncatedCurveModelParam( const ARM_PositiveTruncatedCurveModelParam& rhs )
	:	ARM_RangedCurveModelParam(rhs) {}
	ARM_PositiveTruncatedCurveModelParam& operator=( const ARM_PositiveTruncatedCurveModelParam& rhs )
	{
		if( this != &rhs )
			ARM_RangedCurveModelParam::operator =(rhs);
		return *this;
	}
	virtual ~ARM_PositiveTruncatedCurveModelParam(){};

	/// function to transform model param
	virtual double ModifiedValue( double value ) const
	{ return value<K_NEW_DOUBLE_TOL? K_NEW_DOUBLE_TOL: value; }
};




///////////////////////////////////////////////////
/// \class ARM_BoundedCurveModelParam
/// \brief bounded between two values of model param
/////////////////////////////////////////////////////

class ARM_BoundedCurveModelParam : public ARM_RangedCurveModelParam
{
private:
	double itsMinPotentialValue;
	double itsMaxPotentialValue;

public:
	/// constructor, copy constructor, operator=,destructor
	ARM_BoundedCurveModelParam( ParamType type,
		ARM_GP_Vector* values,
		ARM_GP_Vector* breakPointTimes	, 
		const string& paramName			= "",
        const string& interpolatorName	= "STEPUPRIGHT",
		double minPotentialValue = -1.0,
		double maxPotentialValue =  1.0)
	:	
		ARM_RangedCurveModelParam(type,values,breakPointTimes,paramName,interpolatorName ),
		itsMinPotentialValue(minPotentialValue),
		itsMaxPotentialValue(maxPotentialValue)
	{}

    ARM_BoundedCurveModelParam( const ARM_BoundedCurveModelParam& rhs )
	:	
		ARM_RangedCurveModelParam(rhs),
		itsMinPotentialValue(rhs.itsMinPotentialValue),
		itsMaxPotentialValue(rhs.itsMaxPotentialValue)
	{}

	ARM_BoundedCurveModelParam& operator=( const ARM_BoundedCurveModelParam& rhs )
	{
		if( this != &rhs )
		{
			ARM_RangedCurveModelParam::operator =(rhs);
			itsMinPotentialValue = rhs.itsMinPotentialValue;
			itsMaxPotentialValue = rhs.itsMaxPotentialValue;
		}
		return *this;
	}

	virtual ~ARM_BoundedCurveModelParam(){};

	/// function to transform model param
	virtual double ModifiedValue( double value ) const
	{ return value<itsMinPotentialValue? itsMinPotentialValue: value>itsMaxPotentialValue? itsMaxPotentialValue : value; }
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

