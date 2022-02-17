/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: zccrvfunctor.h,v $
 * Revision 1.1  2004/04/25 07:51:48  ebenhamou, jmprie
 * Initial revision
 *
 *
 */

/*! \file zccrvfunctor.h
 *
 *  \brief 
 *	\author  E Benhamou, JM Prie
 *	\version 1.0
 *	\date April 2004
 */


#ifndef _INGPINFRA_ZCCRVTFUNCTOR_H
#define _INGPINFRA_ZCCRVTFUNCTOR_H

#include "gpbase/port.h"
#include "pricingmodel.h"
#include <string>
CC_USING_NS(std,string)
#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )

struct ARM_ZeroCurveFunctor
{
    ARM_ZeroCurveFunctor(ARM_PricingModel* model):   
        itsModel(model)
        {}

	ARM_ZeroCurveFunctor* Duplicate( const ARM_PricingModel& RhsModel, ARM_PricingModel* LhsModel ) const
	{
		/// two models are equal if pointing on the
		if( itsModel == &RhsModel )
			return new ARM_ZeroCurveFunctor( LhsModel);
		else
			return new ARM_ZeroCurveFunctor( *this );
	}

    ARM_VectorPtr DiscountFactor( const string& curveName, double evalTime, double maturityTime, const ARM_PricingStatesPtr& states) const
    {
#if defined( __GP_STRICT_VALIDATION)
		if( !itsModel )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": model is null!" );
#endif
        return itsModel->DiscountFactor( curveName , evalTime, maturityTime, states );
    }

    bool IsSameModel(const ARM_ZeroCurveFunctor& rhs) const {return itsModel == rhs.itsModel;}



private:
    ARM_PricingModel* itsModel;
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

