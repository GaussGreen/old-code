/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file lexerdec.h
 *
 * \brief This files is the skeleton for the automatically generated
 * lexer...
 * to generate the file  use flex -Cfa -f -t -8 -L genpricerlexer.l > genpricerlexer.cpp
 *
 * \author  E. Benhamou
 * \version 1.0
 * \date October 2003
 */

#ifndef _INGPINFRA_LEXERDEC_H
#define _INGPINFRA_LEXERDEC_H

/// this headers has to come first
/// as it contains pragma for stupid warning
#include "gpbase/removeidentifiedwarning.h"

#include "gpbase/port.h"
#include <string>
#include <vector>

/// header for read declaration
#ifdef unix
	#include <unistd.h>		/// defined here in unix
#else
	#include <io.h>			/// defined here in NT
#endif

/// the standard namespace in gen pricer in ARM
CC_BEGIN_NAMESPACE( ARM )


class ARM_GenPricerLexer
{
public:

	enum 
	{
		EOT = 0,
		ERR = 1,
		INTEGER = 2,
		
		DOUBLE,
		SYMBOL,
		MATURITY,
		DATE,
		STRING,

		PLUS,				
		MINUS,				
		MULTIPLIES,
		DIVIDES,
		
		COMMA,
		
		LPAREN,
		RPAREN,
		
		LBRACKET,
		RBRACKET,

		MODULUS,
		EQUAL_TO,
		NOT_EQUAL_TO,
		LESS,
		GREATER,
		LESS_EQUAL,
		GREATER_EQUAL,

		LOGICAL_NOT, 
		LOGICAL_AND,
		LOGICAL_OR,  
		ASSIGN,
		
		LAST_TOKEN_NOT_USED
	};
	

	ARM_GenPricerLexer()
		: itsTokens(),
		itsNextTok(itsTokens.end()),
		itsPeekTok(itsTokens.end())
	{}
	
	~ARM_GenPricerLexer()
	{}
	
	CC_NS(std, string) peekText() const
	{
		return itsPeekTok != itsTokens.end() ? itsPeekTok->second.second : "";
	}

	int peekPos() const
	{
		return itsPeekTok != itsTokens.end() ? itsPeekTok->second.first : -1;
	}
	
	int peekToken() 
	{ 
		return itsPeekTok != itsTokens.end() ? (itsPeekTok++)->first : EOT;          
	}
	
	virtual void unPeek()
	{ 
		itsPeekTok = itsNextTok;                                           
	}
	
	int nextToken()
	{
		int tok;
		if( itsNextTok != itsTokens.end() )
		{
			itsText = itsNextTok->second.second;
			itsPos = itsNextTok->second.first;
			tok = itsNextTok->first;
			++itsNextTok;
		}
		else
			itsText = itsTokenNames[ tok = EOT ];
		itsPeekTok = itsNextTok;
		return tok;
	}
	
	
	static CC_NS(std, string) tokenText( int token )
	{    
		if( token < 0 || token >= LAST_TOKEN_NOT_USED )
			return "<unknown_token>";
		
		return itsTokenNames[token];
	}
	
	CC_NS(std, string) text() const
	{
		return itsText;
	}
	

	void setInput( CC_NS(std, string) const& inp );

private:

	typedef CC_NS( std, pair )< int, CC_NS(std, string) > TokenTextAndPosT;
	typedef CC_NS( std, pair )< int, TokenTextAndPosT > TokenTextPairT;
	typedef CC_STL_VECTOR( TokenTextPairT ) TokenVecType;
	
	TokenVecType itsTokens;
	TokenVecType::iterator itsNextTok;
	TokenVecType::iterator itsPeekTok;
	CC_NS(std, string) itsText;
	int itsPos;
	
	static const char* itsTokenNames[];
};


CC_END_NAMESPACE()


#endif 
