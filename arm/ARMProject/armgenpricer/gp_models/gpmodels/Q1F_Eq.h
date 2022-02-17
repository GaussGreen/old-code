/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Q1F_Eq.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


#ifndef _INGPMODELS_Q1F_EQ_H
#define _INGPMODELS_Q1F_EQ_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix

//gpbase
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/assignop.h"

//gpmodels
#include "gpmodels/EqFxBase.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/ModelParamsQ1F.h"


/// forward declaration
class ARM_ZeroCurve;

CC_BEGIN_NAMESPACE( ARM )


class ARM_QModel1F_Eq :  public ARM_EqFxBase
{
public:
	typedef ARM_ModelParams_Eq_T<ARM_ModelParamsQ1F> ARM_ModelParamsQ1F_Eq;

	ARM_QModel1F_Eq(const ARM_ZeroCurvePtr& zc, ARM_ModelParamsQ1F_Eq* modelParam);

	ARM_QModel1F_Eq( const ARM_QModel1F_Eq& rhs )
	:	ARM_EqFxBase(rhs) {}
	
	ASSIGN_OPERATOR(ARM_QModel1F_Eq)
	virtual ~ARM_QModel1F_Eq(){};

	/// convention support : calendar + gap support
	virtual string ComputeSettlementCalendar(const string& modelName="") const;
	virtual double ComputeSettlementGap(const string& modelName="") const;

	virtual int GetType() const;

	/// ARM Support
	virtual ARM_Object* Clone() const { return new ARM_QModel1F_Eq( *this ); }
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual string ExportShortName() const { return "LQ1EQ";}
};



CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
