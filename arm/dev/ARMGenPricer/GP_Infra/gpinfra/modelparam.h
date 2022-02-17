/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: modelparam.h,v $
 * Revision 1.1  2003/10/08 16:45:26  ebenhamou
 * Initial revision
 *
 */


/*! \file modelparam.h
 *
 *  \brief model param is a general object for model parameter
 *	\author  E. Benhamou JM Prie
 *	\version 1.0
 *	\date October 2003
 */

#ifndef _INGPINFRA_MODELPARAM_H
#define _INGPINFRA_MODELPARAM_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include "modelparamtype.h"
#include "typedef.h"
#include <string>

CC_BEGIN_NAMESPACE( ARM )

///////////////////////////////////////////////////
/// \class ARM_ModelParam
/// \brief
/// object to represent all the Model Param
/////////////////////////////////////////////////////


/// forward declaration for downcast
class ARM_CurveModelParam;
class ARM_SurfaceModelParam;
class ARM_SurfaceListModelParam;

class ARM_ModelParam : public ARM_RootObject

{
public:
// FIXMEFRED: mig.vc8 (22/05/2007 15:56:04):needs this type to be public
	typedef ARM_ModelParamType::ParamNb ParamType;


private:
	ParamType itsType;
	bool itsAdviseBreakPointTimes;

public:
	virtual int ARM_ModelParamDimension() const=0;
	/// constructor, copy assignment and destructor
	ARM_ModelParam( ParamType type, bool adviseBreakPointTimes );
    ARM_ModelParam( const ARM_ModelParam& rhs );
	ARM_ModelParam& operator=( const ARM_ModelParam& rhs );
	virtual ~ARM_ModelParam();

	/// type methods
	string GetTypeString() const{ return ARM_ModelParamType::GetTypeString( itsType ); }
    inline ParamType GetType() const { return itsType; }
	inline void SetType( ParamType type ) { itsType = type; }
	inline bool GetAdviseBreakPointTimes() const { return itsAdviseBreakPointTimes; }

	/// downcast
	virtual ARM_CurveModelParam& ToCurveModelParam();
	virtual const ARM_CurveModelParam& ToCurveModelParam() const;
	virtual ARM_SurfaceModelParam& ToSurfaceModelParam();
	virtual const ARM_SurfaceModelParam& ToSurfaceModelParam() const;
	virtual ARM_SurfaceListModelParam& ToSurfaceListModelParam();
	virtual const ARM_SurfaceListModelParam& ToSurfaceListModelParam() const;
	
	/// model param value management
	virtual double GetValue( double x1 ) const							= 0;
	virtual double GetValue( double x1, double x2 ) const				= 0;
	
	virtual double GetTimeAtPoint( size_t i ) const					    = 0;
	virtual double GetTenorAtPoint( size_t i ) const					= 0;
	virtual double GetTimeAtPoint( size_t i, size_t j ) const			= 0;
	virtual double GetTenorAtPoint( size_t i, size_t j ) const			= 0;

	virtual double GetValueAtPoint( size_t i ) const					= 0;
	virtual double GetValueAtPoint( size_t i, size_t j ) const			= 0;	

	virtual void SetValue( double x1, double value )					= 0;
	virtual void SetValue( double x1, double x2, double value  )		= 0;
	
	virtual void SetValueAtPoint( size_t i, double value )				= 0;
	virtual void SetValueAtPoint( size_t i, size_t j, double value )	= 0;
	
    virtual void MergeModelParam( ARM_ModelParam* NewValue )			= 0;
	virtual size_t size() const											= 0;
	virtual size_t rows() const											= 0;
	virtual size_t cols() const											= 0;

	virtual double GetValue( double x1, double x2, double x3 ) const	= 0;
	virtual double GetValueAtPoint( size_t i, size_t j, size_t k) const	= 0;
	virtual void SetValue(double x1, double x2, double x3, double value)= 0;
	virtual void SetValueAtPoint( size_t i, size_t j, size_t k, double value)	= 0;

	/// to bump a model Param
	virtual void BumpModelParam(size_t row = 0, size_t column = 0, double shift = 0.0,
		ARM_BumpParamType bumpType = ARM_ModelParamBump::isNothing)  {   
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "BumpModelParam() not available in this context" ); 
	}
	//
	// purpose = surface list 
	// for consistency with existing code, this method should be = 0, and reimplemented in specialized classes
	// but it's much easier like that....
	virtual size_t GetSize() const
	{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetSize() not available in this context" ); }
	virtual double GetValue( size_t i,  double x1, double x2 ) const	
	{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "GetValue(i, x1, x2) not available in this context" ); }
	virtual void SetValue(size_t i, double x1, double x2, double value)
	{	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "SetValue(i, x1, x2) not available in this context" ); }

	virtual string ExportShortName() const { return "LMDLP";}

};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

