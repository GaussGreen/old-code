/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * \file gramfunctorbase.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */



#ifndef _INGPINFRA_GRAMFUNCTORBASE_H
#define _INGPINFRA_GRAMFUNCTORBASE_H

/// this header comes first as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

/// env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "typedef.h"
#include <string>
CC_USING_NS(std,string)

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_PricingModel;

/// base clase for all Grammar functor
struct ARM_GramFctor
{
    ARM_GramFctor() : itsPayModelName(""), itsAlreadyComputed(false) {}
    ARM_GramFctor(const ARM_GramFctor& rhs) :   itsPayModelName(rhs.itsPayModelName),
                                                itsAlreadyComputed(rhs.itsAlreadyComputed) {}

	/// pure virtual to force redefinition
	virtual ARM_GramFctorArg operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes ) = 0;

    virtual ARM_NodeInfo GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
		double evalDate, vector<ARM_ExpNodePtr>& nodes );

	/// function to clone the gram functor if it is using some mutable arguments pure virtual to force redefinition
	virtual ARM_GramFctorPtr Clone() const = 0;

	/// virtual to avoid slicing
	virtual ~ARM_GramFctor() {};

    inline void SetPayModelName(const string& payModelName) { itsPayModelName=payModelName; }
    inline const string& GetPayModelName() { return itsPayModelName; }

    inline void SetAlreadyComputed(bool alreadyComputed) { itsAlreadyComputed=alreadyComputed; }
    inline bool GetAlreadyComputed() { return itsAlreadyComputed; }

	virtual const string& GetFuncName() { return itsFuncName; }

private:
	/// name for error writting
    static string itsFuncName;

	bool itsAlreadyComputed;

    string itsPayModelName;
};


CC_END_NAMESPACE()

#endif


/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

