/*!
 *
 * Copyright (c) IXIS CIB January 2005 Paris
 *
 *	\file dealdescriptionpreprocessor.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou, A. Schauly
 *	\version 1.0
 *	\date January 2005
 */

#ifndef _INGPINFRA_DEALDESCRIPTIONPREPROCESSOR_H
#define _INGPINFRA_DEALDESCRIPTIONPREPROCESSOR_H


#include "gpbase/port.h"
#include "typedef.h"
#include <string>
CC_USING_NS(std,string)


CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_DealDescription;
class ARM_GenPricerLexer;


class ARM_DescriptionPreprocessor
{
private:
	static string ProcessExerciseRightArg( ARM_GenPricerLexer& Lexer, int requiredEndToken  );
	static string ProcessExerciseLeftArg(  ARM_GenPricerLexer& Lexer, int requiredEndToken  );

public:
	static string ProcessMaxPV( const string& text );
	static ARM_DealDescriptionPtr ChangeMaxPVPreprocessing( const ARM_DealDescription& dealDes );
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

