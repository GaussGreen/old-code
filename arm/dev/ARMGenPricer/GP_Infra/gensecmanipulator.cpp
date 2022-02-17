/*!
 *
 * Copyright (c) IXIS CIB December 2004 Paris
 *
 *
 *	\file gensecurity.cpp
 *
 *  \brief Generic security
 *
 *	\author  A. Schauly
 *	\version 1.0
 *	\date December 2004
 */

//// gpinfra
#include "gpinfra/typedef.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/gramnode.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/gensecmanipulator.h"

/// STL
#include <cstdlib>
#include <algorithm>
#include <functional>

CC_BEGIN_NAMESPACE( ARM )


void ARM_GenSecManipulator::ChangeAmericanIntoTrigger( ARM_GenSecurity& gensec, ARM_DealDescriptionPtr& dealdesc )
{
	ARM_GP_NodeMatrixPtr ParseTree = gensec.getParseTree();

	ARM_GP_NodeMatrix::iterator endIter = ParseTree->end();

	for (size_t i = 0; i < ParseTree->GetRowsNb(); ++i)
		for (size_t j = 0; j< ParseTree->GetColsNb(); ++j)
			(*ParseTree)(i,j)->ChangeExerciseIntoTrigger( *dealdesc );
}

CC_END_NAMESPACE()