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


#include "gpinfra/dealdescription.h"
#include "gpinfra/dealdescriptionpreprocessor.h"
#include "gpinfra/dealdescriptionpreprocessor.h"
#include "gpinfra/lexerdec.h"

#include "gpbase/ostringstream.h"
#include "gpbase/stringmanip.h"

/// ARM Kernel
#include <glob/expt.h>

CC_BEGIN_NAMESPACE( ARM )



////////////////////////////////////////////////////
///	Class  : ARM_DescriptionPreprocessor 
///	Routine: ChangeMaxPVPreprocessing 
///	Returns: ARM_DealDescriptionPtr
///	Action : process the text and change the max pv
////////////////////////////////////////////////////

ARM_DealDescriptionPtr ARM_DescriptionPreprocessor::ChangeMaxPVPreprocessing( const ARM_DealDescription& dealDes )
{
	ARM_StringVector text				= dealDes.GetText();
	size_t rowsNb						= dealDes.GetRowsNb();
	size_t colsNb						= dealDes.GetColsNb();

	for( size_t i=1; i<rowsNb; ++i )
		for( size_t j=1; j<colsNb; ++j )
			text[i*colsNb+j] = ARM_DescriptionPreprocessor::ProcessMaxPV(text[i*colsNb+j] );

	ARM_DealDescription* newDealDescription = new ARM_DealDescription( text, dealDes.GetFormat(), rowsNb, colsNb );
	return static_cast<ARM_DealDescriptionPtr>(newDealDescription);
}

////////////////////////////////////////////////////
///	Class  : ARM_DescriptionPreprocessor 
///	Routine: ProcessMaxPV 
///	Returns: string
///	Action : process the text and change the max pv
////////////////////////////////////////////////////

string ARM_DescriptionPreprocessor::ProcessMaxPV( const string& text )
{
	CC_Ostringstream os;
	ARM_GenPricerLexer	Lexer;
	Lexer.setInput( text );
	
	int	nextToken = Lexer.nextToken(),
		peekToken;


	while( ARM_GenPricerLexer::EOT != nextToken )
	{
		if( nextToken != ARM_GenPricerLexer::SYMBOL )
		{
			os << Lexer.text();
			nextToken=Lexer.nextToken();
		}
		else
		{
			peekToken = Lexer.peekToken();
			switch( peekToken )
			{
				/// this is reference
			case ARM_GenPricerLexer::LBRACKET:
				{
					os << Lexer.text();
					nextToken = Lexer.nextToken();
					
					while( ARM_GenPricerLexer::RBRACKET != nextToken )
					{
						os << Lexer.text();
						nextToken=Lexer.nextToken();
					}
					
					os << Lexer.text();
					nextToken = Lexer.nextToken();
				}
				break;
				
				/// this is a function
			case ARM_GenPricerLexer::LPAREN:
				{
					string functionName = Lexer.text();
					bool isAnExercisePV = false;
					
					/// Max cases
					if( stringGetUpper(functionName ) == "MAX" )
					{
						/// peek one token away and get the argument after
						string argName = Lexer.peekText();
						string exerciseLeftArg, exerciseRightArg, stdMaxLefArg, stdMaxRightArg;
						
						/// type Max(PV(...), )
						if( stringGetUpper( argName ) == "PV" )
						{
							/// skip the Max(PV( (hence 4 next tokens!)
							Lexer.nextToken();
							Lexer.nextToken();
							Lexer.nextToken();
							Lexer.nextToken();
							exerciseRightArg= ARM_DescriptionPreprocessor::ProcessExerciseRightArg( Lexer, ARM_GenPricerLexer::COMMA );

							/// skip the comma
							Lexer.nextToken();

							exerciseLeftArg	= ARM_DescriptionPreprocessor::ProcessExerciseLeftArg( Lexer, ARM_GenPricerLexer::RPAREN );
							isAnExercisePV	= true;
						}
						/// type Max(...,PV( ) )
						else
						{
							/// skip Max(
							Lexer.nextToken();
							Lexer.nextToken();
							stdMaxLefArg = ARM_DescriptionPreprocessor::ProcessExerciseLeftArg( Lexer, ARM_GenPricerLexer::COMMA );

							/// skip the comma
							Lexer.nextToken();

							/// current text
							argName = Lexer.text();
							
							if( stringGetUpper( argName ) == "PV" )
							{
								exerciseLeftArg	= stdMaxLefArg;
								/// skip PV(
								Lexer.nextToken();
								Lexer.nextToken();
								exerciseRightArg= ARM_DescriptionPreprocessor::ProcessExerciseRightArg( Lexer, ARM_GenPricerLexer::RPAREN );
								isAnExercisePV	= true;
							}
							else
							{
								stdMaxRightArg = ARM_DescriptionPreprocessor::ProcessExerciseLeftArg( Lexer, ARM_GenPricerLexer::RPAREN );
							}
						}
						
						if( isAnExercisePV )
						{
							os << "Exercise(" << exerciseLeftArg<< "," <<  exerciseRightArg << ")";
						}
						else
						{
							os << "Max(" << stdMaxLefArg << "," << stdMaxRightArg << ")";
						}
						nextToken = Lexer.nextToken();
					}
					
					/// PV cases
					else if( stringGetUpper(functionName ) == "PV" )
					{
						ARM_THROW( ERR_INVALID_ARGUMENT, " PV without MAX not allowed!" );
					}
					
					/// other cases
					else
					{
						/// write the function Name
						os << functionName;
						Lexer.nextToken();
						/// write the "("
						os << Lexer.text();
						nextToken = Lexer.nextToken();
						
						int parenthesisNb = 1;
						
						CC_Ostringstream os2;
						while( parenthesisNb )
						{
							switch( nextToken )
							{
							case ARM_GenPricerLexer::RPAREN:
								--parenthesisNb;
								break;
							case ARM_GenPricerLexer::LPAREN:
								++parenthesisNb;
								break;
							default: /// nothing
								break;
							}
							if( parenthesisNb )
							{
								os2 << Lexer.text();
								nextToken=Lexer.nextToken();
							}
						}
						
						os << ProcessMaxPV( os2.str() );
						/// add the last parenthesis
						os << ")";
						nextToken =Lexer.nextToken();
					}
				}
				break;
			default:
				ARM_THROW( ERR_INVALID_ARGUMENT, " a symbol can only be followed by a LPAREN or a LBRACKET!" );
				break;
			}
		}
	}

	return os.str();
}



////////////////////////////////////////////////////
///	Class  : ARM_DescriptionPreprocessor 
///	Routine: ProcessExerciseRightArg 
///	Returns: string
///	Action : process a text of the type PV( ), and check that this is the case!
///				
////////////////////////////////////////////////////

string ARM_DescriptionPreprocessor::ProcessExerciseRightArg( ARM_GenPricerLexer& Lexer, int requiredEndToken  )
{
	CC_Ostringstream os2;
	os2 << Lexer.text();
	int nextToken = Lexer.nextToken();

	while( ARM_GenPricerLexer::RPAREN != nextToken )
	{
		os2 << Lexer.text();
		nextToken = Lexer.nextToken();
	}

	/// skip the RPAREN 
	nextToken = Lexer.nextToken();

	if( nextToken != requiredEndToken )
	{
		CC_Ostringstream os3;
		os3 << " expected after PV the following end token " << Lexer.tokenText( requiredEndToken )
			<< " but found " << Lexer.tokenText( nextToken );
		ARM_THROW(ERR_INVALID_ARGUMENT, os3.str() );
	}

	return os2.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_DescriptionPreprocessor 
///	Routine: ProcessExerciseLeftArg 
///	Returns: string
///	Action : process a text of the type ..., without any PV!
////////////////////////////////////////////////////

string ARM_DescriptionPreprocessor::ProcessExerciseLeftArg( ARM_GenPricerLexer& Lexer, int requiredEndToken  )
{
	CC_Ostringstream os2;
	os2 << Lexer.text();
	int nextToken = Lexer.nextToken();

	while( requiredEndToken != nextToken )
	{
		string elemTxt = Lexer.text();
		if(  stringGetUpper(elemTxt ) == "PV" )
			ARM_THROW(ERR_INVALID_ARGUMENT, "did not expect a PV but found one!" );

		os2 << Lexer.text();
		nextToken = Lexer.nextToken();
	}

	return os2.str();
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

