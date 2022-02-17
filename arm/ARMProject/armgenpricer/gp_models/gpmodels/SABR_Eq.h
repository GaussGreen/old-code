/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file SABR_Eq.h
 *
 *  \brief prototype model for the generic pricer
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date April 2005
 */

#ifndef _INGPMODELS_SABR_EQMOD_H
#define _INGPMODELS_SABR_EQMOD_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix

// gpbase
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/assignop.h"

// gpmodels
#include "gpmodels/EqFxBase.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/SABR_ModelParams.h"

/// forward declaration
class ARM_ZeroCurve;

CC_BEGIN_NAMESPACE( ARM )


/// \class ARM_SABR_Eq


class ARM_SABR_Eq :  public ARM_EqFxBase
{
public:
	typedef ARM_ModelParams_Eq_T<ARM_SABR_ModelParams> ARM_ModelParamsSABR_Eq;

	ARM_SABR_Eq(const ARM_ZeroCurvePtr& zc, ARM_ModelParamsSABR_Eq* modelParam);

	ARM_SABR_Eq( const ARM_SABR_Eq& rhs )
	:	ARM_EqFxBase(rhs) {}
	
	ASSIGN_OPERATOR(ARM_SABR_Eq)
	virtual ~ARM_SABR_Eq(){};

	/// convention support : calendar + gap support
	virtual string GetSettlementCalendar(const string& modelName="") const;
	virtual double GetSettlementGap(const string& modelName="") const;
	virtual int GetType() const;

	/// Forward function
	virtual ARM_VectorPtr Forward(
		const string& curveName, 
        double evalTime,
	    double expiryTime,
	    double settlementTime,
	    double payTime,
        const ARM_PricingStatesPtr& states) const;


	/// ARM Support
	virtual ARM_Object* Clone() const { return new ARM_SABR_Eq( *this ); }
	virtual string ExportShortName() const { return "LSABE";}
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

