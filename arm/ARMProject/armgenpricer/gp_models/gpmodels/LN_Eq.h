/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file LN_Eq.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


#ifndef _INGPMODELS_LN_EQ_H
#define _INGPMODELS_LN_EQ_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix

//gpbase
#include "gpbase/env.h"
#include "gpbase/port.h"

// gpmodels
#include "gpmodels/EqFxBase.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/BS_ModelParams.h"

/// forward declaration
class ARM_ZeroCurve;

CC_BEGIN_NAMESPACE( ARM )


class ARM_LN_Eq :  public ARM_EqFxBase
{
public:
	typedef ARM_ModelParams_Eq_T<ARM_BS_ModelParams> ARM_ModelParamsBS_Eq;

	ARM_LN_Eq(const ARM_ZeroCurvePtr& zc, ARM_ModelParamsBS_Eq* modelParam);

	ARM_LN_Eq( const ARM_LN_Eq& rhs )
	:	ARM_EqFxBase(rhs) {}
	
	ASSIGN_OPERATOR(ARM_LN_Eq)
	virtual ~ARM_LN_Eq(){};

	/// convention support : calendar + gap support
	virtual string GetSettlementCalendar(const string& modelName="") const;
	virtual double GetSettlementGap(const string& modelName="") const;
	virtual int GetType() const;

	/// ARM Support
	virtual ARM_Object* Clone() const { return new ARM_LN_Eq( *this ); }
	virtual string ExportShortName() const { return "LLNEQ";}
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
