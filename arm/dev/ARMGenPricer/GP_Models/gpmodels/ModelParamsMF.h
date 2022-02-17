/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ModelParamsHW.h,v $
 * Revision 1.1  2003/10/13 07:51:48  jmprie
 * Initial revision
 *
 *
 */



/*! \file ModelParamsHW.h
 *
 *  \brief 
 *
 *	\author  A Schauly
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPMODELS_MODELPARAMSMF_H
#define _INGPMODELS_MODELPARAMSMF_H

/// gpmodels
#include "gpmodels/modelparamshw1f.h"

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class ARM_ModelParamsMF
// \brief Interface class for model parameters of Hull & White models
//-----------------------------------------------------------------------------
class ARM_ModelParamsMF : public ARM_ModelParamsHW1FStd 
{
public:
	ARM_ModelParamsMF( const ARM_ModelParamsMF& rhs );
	ARM_ModelParamsMF( const ARM_ModelParamVector& params = ARM_ModelParamVector() );
	virtual ~ARM_ModelParamsMF();

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_ModelParamsMF(*this); };
	virtual string ExportShortName() const { return "LMF1P";}
	virtual string toString(const string& indent="",const string& nextIndent="") const;

	/// Markov Functional specific
	virtual void SetVolUpToT1AndFreezeGlobVarUpToT2(double T1, double T2, double volUpToT1, bool& varSqueeze);
	virtual void SetVolFromT1toT2AndFreezeGlobVarUpToT2(double T1, double T2, double volFromT1toT2, bool& varSqueeze);
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/