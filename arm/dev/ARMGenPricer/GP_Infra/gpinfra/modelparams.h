/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: modelparams.h,v $
 * Revision 1.1  2003/10/08 16:45:26  ebenhamou
 * Initial revision
 *
 */

/*! \file modelparams.h
 *
 *  \brief file for the model params
 *	\author  E. Benhamou JM Prie
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPINFRA_MODELPARAMS_H
#define _INGPINFRA_MODELPARAMS_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include <string>
#include "typedef.h"


CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_ModelFitter;

////////////////////////////////////////////////////
/// \class ARM_ModelParams
/// \brief Model Params is a subsection of a model that
/// allows to communicate easily with the general
/// calibrator
/////////////////////////////////////////////////////
class ARM_ModelParams : public ARM_RootObject
{
private:
	ARM_ModelParamVector itsParams;
	void CleanUp();
	void CopyNoCleanUp( const ARM_ModelParams& rhs );
	void GetModelParamValidate( int paramType ) const;

public:	
	ARM_ModelParams( const ARM_ModelParamVector& params = ARM_ModelParamVector() );
	ARM_ModelParams( const ARM_ModelParams&	ModelParams );
	ARM_ModelParams& operator=( const ARM_ModelParams& rhs );
	virtual ~ARM_ModelParams();

	/// How many factors?
	virtual size_t FactorCount() const { return 0; }
	
	/// How many parameters
	virtual int size() const { return itsParams.size(); }

	/// Parameter names are gettable, but not settable.
	virtual CC_NS( std, string ) GetModelParamName( int, size_t factorNb = 0 ) const;

	/// The names of all the parameters (gettable but not setable)
	virtual ARM_StringVector GetModelParamNames() const;

	/// ------------------- Get/set a parameter.
    /// Test if a model param with a specific type exists
    virtual bool DoesModelParamExist( int paramType, size_t factorNb=0 ) const;

	/// const version
	virtual const ARM_ModelParam& GetModelParam( int paramType, size_t factorNb=0 ) const;

	/// non const version
	virtual ARM_ModelParam& GetModelParam( int paramType, size_t factorNb=0 );

    virtual ARM_ModelParamVector GetModelParams(size_t factorNb=0) const { return itsParams; }
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
    virtual void PreProcessing(ARM_ModelFitter& modelFitter, int factorNb = 0)=0;
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model, int factorNb = 0 )=0 ;

	/// standard ARM Object support
	virtual string toString(const string& indent="", const string& nextIndent="") const;

protected:
    /// STL iterator support
	inline iterator begin(){ return itsParams.begin(); }
	inline iterator end()	{ return itsParams.end(); }
	inline const_iterator begin() const { return itsParams.begin(); }
	inline const_iterator end() const { return itsParams.end(); }
};

CC_END_NAMESPACE()


#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

