/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *


/*! \file ModelParamsMSV1F.h
 *
 *  \brief 
 *
 *	\author  A Triki
 *	\version 1.0
 *	\date October 2005
 */


#ifndef _INGPMODELS_MODELPARAMSMSV1F_H
#define _INGPMODELS_MODELPARAMSMSV1F_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "modelparamsmsv.h"
#include "gpbase/port.h"
#include "gpcalib/typedef.h"			/// for calibParamVector
#include "gpinfra/modelparamtype.h"		/// for ARM_ModelParamType


/// forward declaration
class ARM_IRIndex;

CC_BEGIN_NAMESPACE( ARM )

class ARM_CurveModelParam;

//-----------------------------------------------------------------------------
// \class ARM_ModelParamsMSV1F
// \brief Interface class for model parameters of the MSV 1F model
//-----------------------------------------------------------------------------
class ARM_ModelParamsMSV1F : public ARM_ModelParamsMSV
{
	ARM_IRIndex* itsIRIndex;							/// SFRM ir index (cloned by the constructor)
public:
	ARM_ModelParamsMSV1F( const ARM_ModelParamsMSV1F& rhs );
	ARM_ModelParamsMSV1F( const ARM_ModelParamVector& params=ARM_ModelParamVector() ,ARM_IRIndex* index= NULL);
	virtual ~ARM_ModelParamsMSV1F();
    ARM_ModelParamsMSV1F& operator = (const ARM_ModelParamsMSV1F& rhs);

	/// Dates at which modelparams values change
	virtual ARM_GP_VectorPtr ModelParamsTimeSteps() const;

	/// How many factors?
	virtual size_t FactorCount() const { return 2; }

	/// Standard ARM object support
	virtual ARM_Object* Clone() const;

	inline const ARM_IRIndex* GetIRIndex() const                    { return itsIRIndex; }
	inline ARM_IRIndex* GetIRIndex()                                { return itsIRIndex; }

    /// Coefficient of the state variable (for Zc closed form formula)
    virtual double BetatT(double t,double T) const;

	double Deriv_BetatT(double t,double T) const;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

