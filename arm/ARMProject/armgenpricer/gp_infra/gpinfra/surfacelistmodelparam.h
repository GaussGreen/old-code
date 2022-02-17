/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file surfacelistmodelparam.h
 *
 *  \brief 
 *
 *	\author  a. Rajaona
 *	\version 1.0
 *	\date December 2004
 */

#ifndef _INGPINFRA_SURFACELISTMODELPARAM_H
#define _INGPINFRA_SURFACELISTMODELPARAM_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/surfacetypedef.h"
#include "typedef.h"
#include "gpbase/gpvector.h"
#include "modelparam.h"
#include "surfacemodelparam.h"
#include "gpinfra/calibparamcst.h"



CC_BEGIN_NAMESPACE( ARM )


/////////////////////////////////////////////////////////////////////
/// \class ARM_SurfaceListModelParam
/// \brief
/// ARM_SurfaceListModelParam is a model param with a List of surface
/////////////////////////////////////////////////////////////////////

class ARM_SurfaceListModelParam : public ARM_ModelParam
{
private:
    /// Curve to manage correctly model parmeters
	std::vector<double> itsIndex;
    ARM_SurfacePtrVector itsSurfaceList;
	string itsName;
	double itsLowerBound;
	double itsUpperBound;

public:
	/// dimension
	int ARM_ModelParamDimension() const
	{ 
		return 2;
	}
	/// constructor (default and with argument)
	ARM_SurfaceListModelParam() : ARM_ModelParam( ARM_ModelParamType::Unknown, true )
	{};

	ARM_SurfaceListModelParam( 
		ParamType type, 
		const std::vector<double>& index, 
		const ARM_SurfacePtrVector& surfacelist, 
		const string& paramName		= "",
		double lowerBound			= ARM_CalibParamCst::CalibLowerBound,
		double upperBound			= ARM_CalibParamCst::CalibUpperBound,
		bool adviseBreakPointTimes	= false);

	/// copy constructor
    ARM_SurfaceListModelParam( const ARM_SurfaceListModelParam& rhs );
	ARM_SurfaceListModelParam& operator=( const ARM_SurfaceListModelParam& rhs );
	virtual ~ARM_SurfaceListModelParam();

	/// downcast to surface list model param
    virtual ARM_SurfaceListModelParam& ToSurfaceListModelParam() { return *this; }
    virtual const ARM_SurfaceListModelParam& ToSurfaceListModelParam() const { return *this; }

	// Interpole method
	double Interpolate( double x1, double x2, double x3 ) const;
	
	/// std virtual function of a model param
	virtual double GetValue( double x1 ) const;
	virtual double GetValue( double x1, double x2 ) const;
	virtual double GetValue( double x1, double x2, double x3 ) const;
	virtual double GetValue( size_t i,  double x1, double x2 ) const; // GetValue for the original concept of surface list !

	virtual double GetTimeAtPoint( size_t i ) const	;
	virtual double GetTenorAtPoint( size_t i ) const;
	virtual double GetTimeAtPoint( size_t i, size_t j ) const;
	virtual double GetTenorAtPoint( size_t i, size_t j ) const;

	virtual double GetValueAtPoint( size_t i ) const;
	virtual double GetValueAtPoint( size_t i, size_t j ) const;
	virtual double GetValueAtPoint( size_t i, size_t j, size_t k ) const;

	virtual void SetValue( double x1, double value );
	virtual void SetValue( double x1, double x2, double value  );
	virtual void SetValue( double x1, double x2, double x3, double value  );
	virtual void SetValue( size_t i,  double x1, double x2, double value  ); // SetValue for the original concept of surface list !
	
	virtual void SetValueAtPoint( size_t i, double value );
	virtual void SetValueAtPoint( size_t i, size_t j, double value );
	virtual void SetValueAtPoint( size_t i, size_t j, size_t k, double value );

	double GetLowerBound() const { return itsLowerBound; }
	double GetUpperBound() const { return itsUpperBound; }

	double InterpolateAtPoint( size_t i, size_t j, size_t k ) const;

    virtual void MergeModelParam( ARM_ModelParam* NewValue );
	virtual size_t size() const;
	virtual size_t rows() const;
	virtual size_t cols() const;
	size_t rows(size_t surfIdx) const;
	size_t cols(size_t surfIdx) const;

	/// SurfaceList specific
	virtual size_t GetSize() const { return itsSurfaceList.size(); }
	ARM_SurfacePtrVector& GetSurfaceList ()  { return itsSurfaceList; };
	const ARM_SurfacePtr GetSurface (size_t i) const { return itsSurfaceList[i]; }
	std::vector<double>*GetData( ARM_DataType type, double index, long& rows, long& cols ) const;
	virtual void SetSurfaceList( const ARM_SurfacePtrVector& surfacelist );
	void SetSurfaceListIndex( const std::vector<double>& index );

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
