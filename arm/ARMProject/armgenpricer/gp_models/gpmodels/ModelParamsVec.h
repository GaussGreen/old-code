/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file modelparamsvec.h
 *
 *  \brief file for the model params
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date February 2004
 */


#ifndef _INGPMODELS_MODELPARAMSVEC_H
#define _INGPMODELS_MODELPARAMSVEC_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpinfra/modelparams.h"
#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
/// \class ARM_ModelParamsVec
/// \brief Model Params Vec is a subsection of a model that
/// allows to communicate easily with the general
/// calibrator
/////////////////////////////////////////////////////
class ARM_ModelParamsVec : public ARM_ModelParams
{
private:
	vector<ARM_ModelParams*> itsParamsVector; 
	void CheckRange( size_t factorNb ) const;

public:
	ARM_ModelParamsVec( const vector<ARM_ModelParams*>& paramsVec )
	:	ARM_ModelParams(), itsParamsVector(paramsVec) {}

	ARM_ModelParamsVec( const ARM_ModelParamsVec& rhs )
	:	ARM_ModelParams(), itsParamsVector(rhs.itsParamsVector){}

	ARM_ModelParamsVec& operator=( const ARM_ModelParamsVec& rhs )
	{
		if( this != &rhs )
		{
			ARM_ModelParams::operator =(rhs);
			itsParamsVector = rhs.itsParamsVector;
		}
		return *this;
	}

	virtual ~ARM_ModelParamsVec() {}

	/// returns the corresponding model parameter vector
    virtual ARM_ModelParamVector GetModelParams(size_t factorNb=0) const;

	/// Parameter names are gettable, but not settable.
	virtual CC_NS( std, string ) GetModelParamName( int type, size_t factorNb = 0 ) const;

	/// ------------------- Get/set a parameter.
    /// Test if a model param with a specific type exists
    virtual bool DoesModelParamExist( int paramType, size_t factorNb=0 ) const;

	/// const version
	virtual const ARM_ModelParam& GetModelParam( int paramType, size_t factorNb=0 ) const;

	/// non const version
	virtual ARM_ModelParam& GetModelParam( int paramType, size_t factorNb=0 );

    virtual void SetModelParam(ARM_ModelParam* param, size_t factorNb=0 );

	virtual void DeleteModelParam( int paramType, size_t factorNb=0 );
	
	/// Iterators for going through the individual parameters.  
	///	Note that these are actually wrappers around the vector iterator
	typedef ARM_ModelParamVector::iterator	iterator;
	typedef ARM_ModelParamVector::const_iterator const_iterator;

	/// compared to the SetValue member function of modelParam, this function
    /// can do some preprocessing (useful for SFRM)
    /// WARNING this is a function on paramType and not paramNum
    /// I repeat this is a function on paramType and not paramNum
    virtual iterator SetModelParamValue( int paramType, size_t i, double value, double time, double tenor = 0.0, size_t factorNb=0 );

    /// merge model param with checking that the merged model
    /// param is not null
    virtual void MergeModelParam(ARM_ModelParam* NewValue, size_t factorNb=0 );

    /// initialise the news parameters to calibrate (not pure virtual to allow the use of 
    virtual void PreProcessing(ARM_ModelFitter& modelFitter);
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model );

	/// standard ARM Object support
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual ARM_Object* Clone() const { return new ARM_ModelParamsVec(*this);}
};

CC_END_NAMESPACE()


#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

