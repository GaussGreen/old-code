/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: modelparamsurface.h,v $
 * Revision 1.1  2003/10/08 16:45:26  ebenhamou
 * Initial revision
 *
 */


/*! \file modelparamsurface.h
 *
 *  \brief file for surface model param
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2004
 */


#ifndef _INGPINFRA_MODELPARAMSURFACE_H
#define _INGPINFRA_MODELPARAMSURFACE_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "typedef.h"
#include "gpbase/surfacetypedef.h"
#include "modelparam.h"
#include "gpinfra/calibparamcst.h"

CC_BEGIN_NAMESPACE( ARM )



///////////////////////////////////////////////////
/// \class ARM_SurfaceModelParam
/// \brief
/// ARM_SurfaceModelParam is a model param with a surface
/////////////////////////////////////////////////////

class ARM_SurfaceModelParam : public ARM_ModelParam
{
private:
    /// Curve to manage correctly model parmeters
    ARM_Surface* itsSurface;
	string itsName;

	double itsLowerBound;
	double itsUpperBound;
	/// inital surface to store the value  
	ARM_Surface* itsInitialSurface;

public:
	/// dimension
	int ARM_ModelParamDimension() const
	{ 
		return 2;
	}
	/// constructor (default and with argument)
	ARM_SurfaceModelParam( ParamType type, 
							ARM_Surface* surface, 
							const string& paramName		= "", 
							double lowerBound			= ARM_CalibParamCst::CalibLowerBound,
							double upperBound			= ARM_CalibParamCst::CalibUpperBound,
							bool adviseBreakPointTimes	= false );
	/// copy constructor
    ARM_SurfaceModelParam( const ARM_SurfaceModelParam& rhs );
	ARM_SurfaceModelParam& operator=( const ARM_SurfaceModelParam& rhs );
	virtual ~ARM_SurfaceModelParam();

	/// downcast to surface model param
	virtual ARM_SurfaceModelParam& ToSurfaceModelParam();
	virtual const ARM_SurfaceModelParam& ToSurfaceModelParam() const;

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

	/// Surface specific
	ARM_Surface* GetSurface() const { return itsSurface; }
	std::vector<double>* GetData( ARM_DataType type, long& rows, long& cols ) const;
	virtual void SetSurface( ARM_Surface* surface );

	double GetLowerBound() const { return itsLowerBound; }
	double GetUpperBound() const { return itsUpperBound; }

	 /// Replace old values using interpolate if necessairy
    void SetAndUpDate(std::vector<double>* values);    

	/// standard ARM Object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="", const string& nextIndent="") const;
};

CC_END_NAMESPACE()


#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

