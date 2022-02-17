/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *
 *	\file gensecurity.cpp
 *
 *  \brief Generic security
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#include "gpinfra/gensecurity.h"

//// gpbase
#include "gpbase/env.h"
#include "gpbase/ostringstream.h"
#include "gpbase/pair.h"
#include "gpbase/datestamp.h"
#include "gpbase/stringmanip.h"
#include "gpbase/valuetype.h"

/// gpinfra
#include "gpinfra/gramnode.h"
#include "gpinfra/gramnodeop.h"
#include "gpinfra/retcppcode.h"
#include "gpinfra/modelnrefcall.h"
#include "gpinfra/refnodemanager.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/cstmanager.h"
#include "gpinfra/lexerdec.h"


/// STL
#include <cstdlib>
#include <algorithm>
#include <functional>

CC_BEGIN_NAMESPACE( ARM )


const string ARM_GenSecurity::ParseCellRefText		= "Parse cell reference";
const string ARM_GenSecurity::ParseCellText			= "Parse cell";
const string ARM_GenSecurity::ParseTermText			= "Parse term";
const string ARM_GenSecurity::ParseFunctionText		= "Parse function";
const string ARM_GenSecurity::ParseMaturityText		= "Parse maturity";
const string ARM_GenSecurity::ParseArgText			= "Parse argument";
const string ARM_GenSecurity::UNPAYSpecialWord		= "UNPAY";
const string ARM_GenSecurity::ExerciseSpecialWord	= "Exercise";
const string ARM_GenSecurity::TriggerSpecialWord	= "Trigger";


////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity
///	Static : itsLook4CircularRef says whether we look
///			for circular reference or not! Activated
///			only in debug mode
////////////////////////////////////////////////////
#if defined(_DEBUG) 
	bool ARM_GenSecurity::itsLook4CircularRef = true; 
#else 
	bool ARM_GenSecurity::itsLook4CircularRef = false; 
#endif


////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: IsVersionOutOfDate
///	Returns: bool
///	Action : Tells whether the version is obsolete or not
////////////////////////////////////////////////////

bool ARM_GenSecurity::IsVersionOutOfDate() const
{
	double monthDiff= Compute_MonthDiff(ARM_CurrentDate,ARM_CompileDate);
	const double maxPeriod = 6.0; /// version is obsolete after four months!
	return monthDiff>=maxPeriod;
}


////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ARM_GenSecurity
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////

ARM_GenSecurity::ARM_GenSecurity(
	ARM_CountedPtr< ARM_DealDescription >	DealDescription,
	const string&							payModelName,
	ARM_CstManagerPtr						cstManager,
	bool									ExercBoundaryResetFlag,
    bool									otherPayoffsFlag,
	bool									ivFlag)
:	
	itsDealDescription(DealDescription), 
	itsPayModelName(payModelName),
	itsPrintParseTree(false), 
	itsPricingAdviser(new ARM_PricingAdviser),
	itsGlobalCsts(), 
	itsCstManager(cstManager), 
	itsExercBoundaryResetFlag(ExercBoundaryResetFlag),
    itsOtherPayoffsFlag(otherPayoffsFlag),
	itsIVFlag(ivFlag),
	itsRowDates(ARM_GP_VectorPtr(new ARM_GP_Vector())),
	itsParseTree(ARM_GP_NodeMatrixPtr(new ARM_GP_NodeMatrix()))
{
/*	if( IsVersionOutOfDate() )
	{
		CC_Ostringstream os;
		os << ARM_USERNAME << ": obsolete version: get a new one! Compiled= ";
		os << ARM_CompileDate;

#if defined(_DEBUG)
		os << " (Debug)";
#else
		os << " (Release)";
#endif
		os << " while Today= " << ARM_CurrentDate;
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
*/
	size_t	NumRows = itsDealDescription->GetRowsNb(),
		NumCols = itsDealDescription->GetColsNb(),
		Row, Col;

	/// Does some size validation on the tableau.
	if( NumRows < 2 )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": tableau must have at least 2 rows, 1 for column names and 1 for contents" );
	if( NumCols < 2 )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": tableau must have at least 2 columns, 1 for dates and 1 for cashflows" );

	/// Grab the column names.
	/// checking for unique names! (case insensitive)
	itsColNames = itsDealDescription->GetRow(0)->GetText();
	for( Col = 0; Col < NumCols; ++Col )
	{
		if( ARM_STRING_TYPE != itsDealDescription->GetElemFormat( 0, Col  ) )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": column header should be string: not the case for " + CellName( 0, Col ) + " with text " + CellReducedText( 0, Col )  );

		if( "" == itsColNames[Col] )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": empty string in the column header " + CellName( 0, Col ) + " with text " + CellReducedText( 0, Col )  );

		bool InsertedSucces = itsColNamesMap.insert( pair< const string, int >( StrUpper( itsColNames[Col] ), Col ) ).second;
		if( !InsertedSucces )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + string( ": column header should be different names, case insensitively, not the case for " )
				+ CellName( 0, Col ) + " with text " + CellReducedText( 0, Col ) );
	}

	/// Ensure First Col contains dates
	/// sorted in increasing order
	double dateDouble, previousDateDouble = -1;
	itsRowDates->reserve(NumRows - 1);
	for( Row = 1; Row<NumRows; ++Row )
	{
		if( ARM_DATE_TYPE != itsDealDescription->GetElemFormat( Row, 0 ) )
		{
			CC_Ostringstream os;
			os << "In " << CellName( Row, 0 ) << " with text: " << CellReducedText( Row, 0 )
				<< ": leftmost column must contain known dates";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}
		dateDouble = atof( itsDealDescription->GetElem( Row, 0 ).c_str() );

		/// test dates are sorted in increasing order!
		if( dateDouble <= previousDateDouble )
		{
			CC_Ostringstream os;
			os << CellName( Row, 0 ) << "'s Date =" << ARM_Date( dateDouble ).toString() 
				<< " before " << CellName( Row-1, 0 ) << "'s Date = " << ARM_Date( previousDateDouble ).toString() 
				<< " " << ARM_USERNAME << ": leftmost column should have sorted increasing dates";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}

		itsRowDates->push_back(	dateDouble );
		previousDateDouble = dateDouble;
	}
	

	size_t NumPricedColumns = itsDealDescription->GetNbPricedColumns();

	/// Size the vetor of ARM_ExpNodePtr, and parse trees
	itsParseTree->resize(NumRows-1, NumPricedColumns);
	ARM_RefNNode_Manager refNNode_Manager;
	ARM_ModelCallMap modelCallMap;

	/// First pass parse all and build parse tree
	for( Row=1; Row<NumRows; ++Row )
	{
		for ( Col=0; Col< NumPricedColumns; ++Col )
		{
			/// parse the cell
			(*itsParseTree)(Row-1,Col) = InsertFakeNodeRef(ParseCell( Row, itsDealDescription->GetPricedColumn(Col), refNNode_Manager, modelCallMap ),Row,itsDealDescription->GetPricedColumn(Col));
			/// reset circular references
			refNNode_Manager.ResetCircularReferences();
		}
	}

	/// now that we are done ... set various things for the pricing Adviser!
	/// 1) set all the collected reference nodes to the pricing adviser
	ARM_RefNNode_Manager::ARM_SharedNodeMap::iterator 
		sharedNode		= refNNode_Manager.GetSharedNodeMapBegin(),
		sharedNodeEnd	= refNNode_Manager.GetSharedNodeMapEnd();
	while( sharedNode != sharedNodeEnd )
	{
		itsPricingAdviser->AddCellRefCall( (*sharedNode).second );
		++sharedNode;
	}

	/// 2) set the parse tree pointor to the pricing adviser
	itsPricingAdviser->SetPParseTree( itsParseTree );
	
	/// 3) collect all the dates

	itsPricingAdviser->SetPDates(itsRowDates);

	/// 4) set all the model calls
	ARM_ModelCallMap::iterator 
		modelCall		= modelCallMap.begin(),
		modelCallEnd	= modelCallMap.end();
	while( modelCall != modelCallEnd )
	{
		/// remove duplicate function names
		(*modelCall).second->SortAndRemoveDuplicates();
		itsPricingAdviser->AddModelCall((*modelCall).second );
		++modelCall;
	}

	/// finally, we are done and can set the name!
	CC_ARM_SETNAME(ARM_GENSECURITY);
}

////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: InsertFakeNodeRef
///	Returns: ARM_ExpNodePtr
///	Action : Insert a fake ref node
////////////////////////////////////////////////////
ARM_ExpNodePtr ARM_GenSecurity::InsertFakeNodeRef( const ARM_ExpNodePtr& childNode, int Row, int Col ) const
{

// This is a fake ref node to prevent any calculation in the backward loop of the AMC
#if defined(__GP_SHOW_SHARED_NODE_COORDINATES)
			ARM_RowColCoords ParentCoords( Row, Col );
			return ARM_ExpNodePtr( new ARM_ExpNodeRef( 
				childNode, 
				RowDate(Row), 
				RowDate(Row), 
				CellName( Row, Col ), 
				ParentCoords, 
				false, 
				itsPayModelName ) );
#else
			return ARM_ExpNodePtr( new ARM_ExpNodeRef( 
				childNode, 
				RowDate(Row), 
				RowDate(Row), 
				CellName( Row, Col ),
				false, 
				itsPayModelName ) );
#endif
}


////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ARM_GenSecurity
///	Returns: nothing
///	Action : Constructor
////////////////////////////////////////////////////
ARM_GenSecurity::ARM_GenSecurity( const ARM_GenSecurity& rhs )
:
	ARM_RootObject(rhs),
	itsColNames(rhs.itsColNames),
	itsColNamesMap(rhs.itsColNamesMap),
	//for shared pointers
	itsPrintParseTree(rhs.itsPrintParseTree),
	itsPayModelName(rhs.itsPayModelName),
	itsGlobalCsts(rhs.itsGlobalCsts),
	itsExercBoundaryResetFlag(rhs.itsExercBoundaryResetFlag),
    itsOtherPayoffsFlag(rhs.itsOtherPayoffsFlag)
{
 	itsDealDescription	= rhs.itsDealDescription;	
	itsCstManager		= rhs.itsCstManager;		
	itsPricingAdviser	= rhs.itsPricingAdviser;	
	itsParseTree		= rhs.itsParseTree;			
	itsRowDates			= rhs.itsRowDates;				
	itsPricingAdviser->SetPDates(itsRowDates);
	itsPricingAdviser->SetPParseTree(itsParseTree);
}



////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ~ARM_GenSecurity
///	Returns: nothing
///	Action : Destructor
////////////////////////////////////////////////////

ARM_GenSecurity::~ARM_GenSecurity()
{
	/// nothing as we use smart pointor!
}


////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: operator=
///	Returns: *this with new values
///	Action : 
////////////////////////////////////////////////////

ARM_GenSecurity &ARM_GenSecurity::operator = (const ARM_GenSecurity &rhs)
{
	if( this != &rhs )
	{
		ARM_RootObject::operator=( rhs );
		itsColNames			        = rhs.itsColNames;
		itsColNamesMap		        = rhs.itsColNamesMap;
		itsRowDates			        = rhs.itsRowDates;
		itsParseTree		        = rhs.itsParseTree;
		//for shared pointers
		itsCstManager				= rhs.itsCstManager;		
		itsDealDescription			= rhs.itsDealDescription;	
		itsPricingAdviser			= rhs.itsPricingAdviser;	
		itsPricingAdviser->SetPDates(itsRowDates);
		itsPricingAdviser->SetPParseTree(itsParseTree);
		itsPrintParseTree	        = rhs.itsPrintParseTree;
		itsPayModelName		        = rhs.itsPayModelName;
		itsGlobalCsts		        = rhs.itsGlobalCsts;
	    itsExercBoundaryResetFlag   = rhs.itsExercBoundaryResetFlag;
        itsOtherPayoffsFlag         = rhs.itsOtherPayoffsFlag;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: CellName 
///	Returns: Formatted name of the cell, mainly for error messages.
///	Action : 
////////////////////////////////////////////////////

string ARM_GenSecurity::CellName(
	size_t Row, 
	size_t Col	) const
{
	CC_Ostringstream os;
	os << itsColNames[Col] << "[ " << Row << " ]";
	return os.str(); 	
}


////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: CellText 
///	Returns: string of the cell
///	Action : 
////////////////////////////////////////////////////

string ARM_GenSecurity::CellReducedText(
	size_t Row, 
	size_t Col ) const 
{
	const int length= 150;
	string totalText = itsDealDescription->GetElem(Row,Col);
	if( totalText.size() > length )
		return totalText.substr(0,length);
	else
		return totalText;
}


////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: RequireToken 
///	Returns: int
///	Action : 
////////////////////////////////////////////////////

int ARM_GenSecurity::RequireToken(
	size_t Row, 
	size_t Col, 
	const string& NodeName, 
	int Token, 
	ARM_GenPricerLexer& Lexer ) const
{
	int	LexerToken = Lexer.nextToken();

	if( LexerToken != Token )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + string( ": in " ) + CellName( Row, Col ) + " with text: " + CellReducedText( Row, Col ) + ":  " 
			+ "Parsing a " + NodeName + ". Expected a " + ARM_GenPricerLexer::tokenText( Token ) + " but saw a " 
			+ ARM_GenPricerLexer::tokenText( LexerToken ) );
	return LexerToken;
}


////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: UnexpectedToken 
///	Returns: ARM_ExpNodePtr
///	Action : 
////////////////////////////////////////////////////

ARM_ExpNodePtr ARM_GenSecurity::UnexpectedToken(
	size_t Row, 
	size_t Col,
	int Token, 
	ARM_GenPricerLexer &Lexer ) const
{	
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
		ARM_USERNAME + string( ": in " ) + CellName( Row, Col ) + " with text: " + CellReducedText( Row, Col ) + ": " 
		+ "Unexpected token " + ARM_GenPricerLexer::tokenText( Token ) );
}


////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: RowDate 
///	Returns: Date associated with a row of the tableau body (from col 0)
///	Action : 
////////////////////////////////////////////////////

ARM_Date ARM_GenSecurity::RowDate( size_t Row ) const
{
/// only at debug time
/// for no overhead efficiency!
#ifdef __GP_STRICT_VALIDATION
	if( Row < 1 || Row >= itsDealDescription->GetRowsNb() )
	{
		CC_Ostringstream os;
		os << ARM_USERNAME << ": you're trying to access a nonsense Date for row "
			<< Row << " with boundary between 0 and " << itsDealDescription->GetRowsNb()
			<< ", advise!";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
#endif

	return ((ARM_Date) ((*itsRowDates)[Row-1]));
}



////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ParseCell 
///	Returns: ARM_ExpNodePtr
///	Action : parse the cell at the specified location ( row 0 = first body row )
//			 itsDealDescription must be set up before calling this.
////////////////////////////////////////////////////

ARM_ExpNodePtr ARM_GenSecurity::ParseCell(
	size_t Row,
	size_t Col,
	ARM_RefNNode_Manager& refNNodeManager,
	ARM_ModelCallMap& modelCallMap ) const
{
	/// grammar:
	///		date
	///	|	double
	/// |	string expression

	string Elem = itsDealDescription->GetElem( Row, Col );
	ARM_GP_VALUE_TYPE Type = itsDealDescription->GetElemFormat( Row, Col );
	
	switch( Type )
	{
	case ARM_DATE_TYPE:
		return ARM_ExpNodePtr( new ARM_ExpNodeDate( ARM_Date( atof( Elem.c_str() ) ) ) );

	case ARM_DOUBLE_TYPE:
		return ARM_ExpNodePtr( new ARM_ExpNodeDateOrDouble( atof( Elem.c_str() ) ) );
		
	case ARM_STRING_TYPE:
		/// test for the Terminal Key operator!
		return ParseTerminalKeyAndOtherType(Row, Col, Elem, refNNodeManager, modelCallMap);
		
		/// other cases are errors
	default:
		{
			CC_Ostringstream os;
			os << ARM_USERNAME << ": in "  << CellName( Row, Col ) << " with text: " << CellReducedText( Row, Col ) << "Unknown type.... " << ARM_TYPE::Name( Type ) << " Missing types: accepted are date, double and string";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ParseTerminalKeyAndOtherType
///	Returns: ARM_ExpNodePtr
///	Action : Lex the expression and parse the result.
////////////////////////////////////////////////////

ARM_ExpNodePtr ARM_GenSecurity::ParseTerminalKeyAndOtherType(
	size_t Row, 
	size_t Col, 
	const string& ExpressionString,
	ARM_RefNNode_Manager& refNNodeManager,
	ARM_ModelCallMap& modelCallMap ) const
{
	/// string expr =
	///		UnPay( ) + EOT
	///	|	Exercise( ) + EOT
	///	|	Trigger( ) + EOT
	///	|	Cell
	
	/// We have to create a local copy of the lexer so that we can recursively
	/// lex things when we follow references to other cells we have to
	/// parse.

	ARM_GenPricerLexer	Lexer;
	Lexer.setInput( ExpressionString );
	Lexer.nextToken();
	string  tokenString	= Lexer.text();

	ARM_ExpNodePtr node;

	if(( stringGetUpper(tokenString) == ARM_GenSecurity::UNPAYSpecialWord ) ||
		( stringGetUpper(tokenString) == ARM_GenSecurity::ExerciseSpecialWord ) ||
		( stringGetUpper(tokenString) == ARM_GenSecurity::TriggerSpecialWord ))
	{
		/// reset the lexer!
		Lexer.setInput( ExpressionString );
		node = ParseFunction( Row, Col, Lexer, refNNodeManager, modelCallMap, true );

		/// check that this is the end!
		int LexerToken = Lexer.nextToken();
		
		/// test that this is a final operation!
		if( ARM_GenPricerLexer::EOT != LexerToken )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": in " + CellName( Row, Col ) + " with text: " + CellReducedText( Row, Col ) 
				+ " " + ARM_GenSecurity::UNPAYSpecialWord + " should be a final operation!" );
		
		return node;
	}
	else
	{
		/// other case
		node = ParseCell( Row, Col, ExpressionString, refNNodeManager, modelCallMap );
	}
	return node;
}


////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ParseCell
///	Returns: ARM_ExpNodePtr
///	Action : Lex the expression and parse the result.
////////////////////////////////////////////////////

ARM_ExpNodePtr ARM_GenSecurity::ParseCell(
	size_t Row, 
	size_t Col, 
	const string& ExpressionString,
	ARM_RefNNode_Manager& refNNodeManager,
	ARM_ModelCallMap& modelCallMap ) const
{
	/// string expr =
	///		date string
	///	|	symbolic expr
	
	/// We have to create a local copy of the lexer so that we can recursively
	/// lex things when we follow references to other cells we have to
	/// parse.
	ARM_GenPricerLexer	Lexer;

	Lexer.setInput( ExpressionString );
	
	int		PeekToken	= Lexer.peekToken();

	/// use reference counted pointor
	/// for exception safety

	ARM_ExpNodePtr node;
	if( ARM_GenPricerLexer::DATE == PeekToken )
	{
		Lexer.unPeek();
		node = ParseDateExpr( Row, Col, Lexer );
	}
	else
	{
 		Lexer.unPeek();
		node = ParseExpr( Row, Col, Lexer, refNNodeManager, modelCallMap );
	}

	/// check that this is the end!
	RequireToken( Row, Col, ARM_GenSecurity::ParseCellText, ARM_GenPricerLexer::EOT, Lexer );	
	
	return node;
}


////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ParseExpr
///	Returns: ARM_ExpNodePtr
///	Action : Parse an expression out of the lexer.
////////////////////////////////////////////////////

ARM_ExpNodePtr ARM_GenSecurity::ParseExpr(
	size_t Row, 
	size_t Col, 
	ARM_GenPricerLexer &Lexer,
	ARM_RefNNode_Manager& refNNodeManager,
	ARM_ModelCallMap& modelCallMap ) const
{
	/// symbolic expr =
	///		symb sum expr
	///	|	symb sum expr  Comparaison Operator   symb sum expr 
	
	ARM_ExpNodePtr Sum = ParseSum( Row, Col, Lexer, refNNodeManager, modelCallMap );
	
	int		NextToken = Lexer.peekToken();

	switch( NextToken )
	{
		case ARM_GenPricerLexer::MODULUS:
			{
				Lexer.nextToken();
				Sum = ARM_ExpNodePtr( new ARM_ExpNodeBiOp< CC_NS( ARM_GP, modulus )<double> >( BinaryOpModulus, Sum, 
					ParseSum( Row, Col, Lexer, refNNodeManager, modelCallMap ) ) );
				break;
			}		
		case ARM_GenPricerLexer::EQUAL_TO:	
			{
				Lexer.nextToken();
				Sum = ARM_ExpNodePtr( new ARM_ExpNodeWithDatesBiOp< CC_NS( ARM_GP, equal_to )<double> >( BinaryOpEqualTo, Sum, 
					ParseSum( Row, Col, Lexer, refNNodeManager, modelCallMap ) ) );
				break;
			}		
		case ARM_GenPricerLexer::NOT_EQUAL_TO:
			{
				Lexer.nextToken();
				Sum = ARM_ExpNodePtr( new ARM_ExpNodeWithDatesBiOp< CC_NS( ARM_GP, not_equal_to )<double> >( BinaryOpNotEqualTo, Sum, 
					ParseSum( Row, Col, Lexer, refNNodeManager, modelCallMap ) ) );
				break;
			}		
		case ARM_GenPricerLexer::LESS:		
			{
				Lexer.nextToken();
				Sum = ARM_ExpNodePtr( new ARM_ExpNodeWithDatesBiOp< CC_NS( ARM_GP, less )<double> >( BinaryOpLess, Sum, 
					ParseSum( Row, Col, Lexer, refNNodeManager, modelCallMap ) ) );
				break;
			}		
		case ARM_GenPricerLexer::GREATER:	
			{
				Lexer.nextToken();
				Sum = ARM_ExpNodePtr( new ARM_ExpNodeWithDatesBiOp< CC_NS( ARM_GP, greater )<double> >( BinaryOpGreater , Sum, 
					ParseSum( Row, Col, Lexer, refNNodeManager, modelCallMap ) ) );
				break;
			}		
		case ARM_GenPricerLexer::LESS_EQUAL:
			{
				Lexer.nextToken();
				Sum = ARM_ExpNodePtr( new ARM_ExpNodeWithDatesBiOp< CC_NS( ARM_GP, less_equal )<double> >( BinaryOpLessEqual, Sum, 
					ParseSum( Row, Col, Lexer, refNNodeManager, modelCallMap ) ) );
				break;
			}		
		case ARM_GenPricerLexer::GREATER_EQUAL:
			{
				Lexer.nextToken();
				Sum = ARM_ExpNodePtr( new ARM_ExpNodeWithDatesBiOp< CC_NS( ARM_GP, greater_equal )<double> >( BinaryOpGreaterEqual, Sum, 
					ParseSum( Row, Col, Lexer, refNNodeManager, modelCallMap ) ) );
				break;
			}		
		case ARM_GenPricerLexer::LOGICAL_AND:
			{
				Lexer.nextToken();
				Sum = ARM_ExpNodePtr( new ARM_ExpNodeBiOp< CC_NS( ARM_GP, logical_and )<double> >( BinaryOpLogicalAnd, Sum, 
					ParseSum( Row, Col, Lexer, refNNodeManager, modelCallMap ) ) );
				break;
			}		
		case ARM_GenPricerLexer::LOGICAL_OR:
			{
				Lexer.nextToken();
				Sum = ARM_ExpNodePtr( new ARM_ExpNodeBiOp< CC_NS( ARM_GP, logical_or )<double> >( BinaryOpLogicalOr, Sum, 
					ParseSum( Row, Col, Lexer, refNNodeManager, modelCallMap ) ) );
				break;
			}		
		default: //// do nothing
			break;
	}
	
	Lexer.unPeek();
	return Sum;
}

////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ParseSum
///	Returns: ARM_ExpNodePtr
///	Action : Parse an expression out of the lexer.
////////////////////////////////////////////////////

ARM_ExpNodePtr ARM_GenSecurity::ParseSum(
	size_t Row, 
	size_t Col, 
	ARM_GenPricerLexer &Lexer,
	ARM_RefNNode_Manager& refNNodeManager,
	ARM_ModelCallMap& modelCallMap ) const
{
	/// symbolic sum expr =
	///		symb product
	///	|	symb product (+|- symb product) + (syntax similar to lex)
	
	/// priority is to ParseProduct hence the first call
	ARM_ExpNodePtr Sum = ParseProduct( Row, Col, Lexer, refNNodeManager, modelCallMap );
	int	NextToken = Lexer.peekToken();

	while(	( ARM_GenPricerLexer::PLUS	== NextToken ) || 
			( ARM_GenPricerLexer::MINUS	== NextToken ) 
			)
	{
		/// get next token
		Lexer.nextToken();

		/// parse the node according to its type
		if( ARM_GenPricerLexer::PLUS == NextToken )
			Sum = ARM_ExpNodePtr( new ARM_ExpNodeBiOp< CC_NS( ARM_GP, plus )<double> >( BinaryOpPlus, Sum, 
				ParseProduct( Row, Col, Lexer, refNNodeManager, modelCallMap ) ) );
		else
			Sum = ARM_ExpNodePtr( new ARM_ExpNodeBiOp< CC_NS( ARM_GP, minus )<double> >( BinaryOpMinus, Sum, 
				ParseProduct( Row, Col, Lexer, refNNodeManager, modelCallMap ) ) );

		/// take the next token
		NextToken = Lexer.peekToken();
	}

	Lexer.unPeek();
	return Sum;
}



////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ParseProduct
///	Returns: ARM_ExpNodePtr
///	Action : Parse a product out of the lexer.
////////////////////////////////////////////////////

ARM_ExpNodePtr ARM_GenSecurity::ParseProduct(
	size_t Row, 
	size_t Col, 
	ARM_GenPricerLexer &Lexer,
	ARM_RefNNode_Manager& refNNodeManager,
	ARM_ModelCallMap& modelCallMap ) const
{
	/// symb product =
	///		symb term
	///	|	symb term (*|/ symb term)+ (syntax similar to lex)
	
	ARM_ExpNodePtr Term = ParseTerm( Row, Col, Lexer, refNNodeManager, modelCallMap );
	int	NextToken = Lexer.peekToken();
 
	while(	( ARM_GenPricerLexer::MULTIPLIES== NextToken ) || 
			( ARM_GenPricerLexer::DIVIDES	== NextToken ) 
			)
	{
		/// get next token
		Lexer.nextToken();

		/// parse the node according to its type
		if( ARM_GenPricerLexer::MULTIPLIES == NextToken )
			Term = ARM_ExpNodePtr( new ARM_ExpNodeBiOp< CC_NS( ARM_GP, multiplies )<double> >(BinaryOpMultiplies, Term, 
				ParseTerm( Row, Col, Lexer, refNNodeManager, modelCallMap ) ) );
		else
			Term = ARM_ExpNodePtr( new ARM_ExpNodeBiOp< CC_NS( ARM_GP, divides )<double> >(BinaryOpDivides, Term, 
				ParseTerm( Row, Col, Lexer, refNNodeManager, modelCallMap ) ) );

		/// get the next one!
		NextToken = Lexer.peekToken();
	}
	
	Lexer.unPeek();
	
	return Term;
}



////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ParseTerm
///	Returns: ARM_ExpNodePtr
///	Action : Parse an expression out of the lexer.
////////////////////////////////////////////////////

ARM_ExpNodePtr ARM_GenSecurity::ParseTerm(
	size_t Row, 
	size_t Col, 
	ARM_GenPricerLexer &Lexer,
	ARM_RefNNode_Manager& refNNodeManager,
	ARM_ModelCallMap& modelCallMap ) const
{
	/// symb term =
	///		integer
	///	|	double
	///	|	function( ... ) = symbol + ( .. )
	///	|	reference( ... )= symbol + [ .. ]
	/// |   cst				= symbol from the cst manager
	///	|	model-id
	///	|	( symb sum )
	///	|	- symb sum		= unary minus
	///	|   ! symb sum		= negate!
	
	
	ARM_ExpNodePtr Result;
	
	int Token = Lexer.peekToken();
	
	switch(Token)
	{
	/// number
	case ARM_GenPricerLexer::INTEGER:
	case ARM_GenPricerLexer::DOUBLE:
		{			
			Lexer.nextToken();
			Result = ARM_ExpNodePtr( new ARM_ExpNodeDateOrDouble( atof( Lexer.text().c_str() ) ) );
			break;	
		}
	
	/// symbol
	case ARM_GenPricerLexer::SYMBOL:
		{
			int NextToken = Lexer.peekToken();

			Lexer.unPeek();
			
			switch( NextToken )
			{
			case ARM_GenPricerLexer::LPAREN:
				Result = ParseFunction( Row, Col, Lexer, refNNodeManager, modelCallMap, false );
				break;
				
			case ARM_GenPricerLexer::LBRACKET:
				Result = ParseCellRef( Row, Col, Lexer, refNNodeManager, modelCallMap );
				break;
				
			default:
				Result = ParseCst( Row, Col, Lexer );
				break;
			}
		break;
		}
		
	case ARM_GenPricerLexer::LPAREN:
		{
			Lexer.nextToken();
			Result = ParseExpr( Row, Col, Lexer, refNNodeManager, modelCallMap );

			/// check that we close the parenthesis!
			RequireToken( Row, Col, ARM_GenSecurity::ParseTermText, ARM_GenPricerLexer::RPAREN, Lexer );

			Result = ARM_ExpNodePtr( new ARM_ExpNodeFactor( Result ) );
			break;
		}
		
	/// unary minus	and logical not
	case ARM_GenPricerLexer::MINUS: 
		{
			Lexer.nextToken();
			Result = ParseTerm( Row, Col, Lexer, refNNodeManager, modelCallMap );

			Result = ARM_ExpNodePtr( new ARM_ExpNodeUnaryOp< CC_NS( ARM_GP, negate )<double> >( UnaryOpNegate, Result ) );
			break;

		}
	case ARM_GenPricerLexer::LOGICAL_NOT: 
		{
			Lexer.nextToken();
			Result = ParseTerm( Row, Col, Lexer, refNNodeManager, modelCallMap );

			Result = ARM_ExpNodePtr( new ARM_ExpNodeUnaryOp< CC_NS( ARM_GP, logical_not )<double> >( UnaryOpLogicalNot, Result ) );
			break;
		}
		
	case ARM_GenPricerLexer::STRING:
		{
			Lexer.unPeek();
			
			Result = ARM_ExpNodePtr( ParseString( Lexer ) );
			break;			
		}
		
	default:
		return UnexpectedToken( Row, Col, Token, Lexer );
	}

	return Result;
}


	
////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ParseCellRef
///	Returns: ARM_RowColCoords
///	Action : Parse an CellRef out of the lexer.
////////////////////////////////////////////////////

ARM_RowColCoords ARM_GenSecurity::ParseCellRefCoords(
	size_t Row, 
	size_t Col, 
	ARM_GenPricerLexer &Lexer,
	set< ARM_RowColCoords >& CircularRefCoords ) const
{
	int Token = RequireToken( Row, Col, ARM_GenSecurity::ParseCellRefText, ARM_GenPricerLexer::SYMBOL, Lexer ); 
	string ColName = Lexer.text();
	
	/// string insensitive lookup!
	CC_NS(std,map)< string, int, less< string > >::const_iterator
		Found = itsColNamesMap.find( StrUpper( ColName ) );
	
	if( Found == itsColNamesMap.end() )
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			string( "In " ) + CellName( Row, Col ) + " with text: " + CellReducedText( Row, Col )
			+ " Unknown column name " + ColName );	
	}
	
	RequireToken( Row, Col, ARM_GenSecurity::ParseCellRefText, ARM_GenPricerLexer::LBRACKET, Lexer );
	int	NextToken = Lexer.nextToken();
	size_t	ColRef = (*Found).second,	RowRef = 0;
	
	switch( NextToken )
	{
	case ARM_GenPricerLexer::INTEGER:
		{	
			RowRef = atoi( Lexer.text().c_str() );
			break;
		}
		
	case ARM_GenPricerLexer::SYMBOL:
		{
			/// check that this is exactly i case insensitive
			if(	("i" != Lexer.text() ) &&
				("I" != Lexer.text() ) ) 
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				string( "In " ) + CellName( Row, Col ) + " with text: " + CellReducedText( Row, Col ) + ": " 
				+ "expected i or I, and found " + Lexer.text() );	
			
			int		OpToken = Lexer.peekToken();
			
			switch( OpToken )
			{
			case ARM_GenPricerLexer::PLUS:
			case ARM_GenPricerLexer::MINUS:
				{
					Lexer.nextToken();
					
					RequireToken( Row, Col, ARM_GenSecurity::ParseCellText, ARM_GenPricerLexer::INTEGER, Lexer );
					
					size_t	Offset = atoi( Lexer.text().c_str() );
					
					if( ARM_GenPricerLexer::PLUS == OpToken )
					{
						RowRef = Row + Offset; 
						
						if( RowRef >= itsDealDescription->GetRowsNb() )
							throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
							string( "In " ) + CellName( Row, Col ) + " with text: " + CellReducedText( Row, Col ) + " "
							+ "Reference to " + CellName( RowRef, ColRef ) +  " beyond end of tableau: "
							+ CellName( itsDealDescription->GetRowsNb() -1, ColRef ) ); 
					}
					else
					{
						RowRef = Row - Offset;
						
						if( RowRef <= 0 )
							throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
							string( "In " ) + CellName( Row, Col ) + " with text: " + CellReducedText( Row, Col ) + " "
							+ "Reference to " + CellName( RowRef, ColRef ) +" before start of tableau: "
							+ CellName( 1, ColRef ) ); 
					}
					
					break;
				}
				
			case ARM_GenPricerLexer::RBRACKET:
				{
					RowRef = Row;
					break;
				}
			
			default:
				UnexpectedToken( Row, Col, NextToken, Lexer );					
				break;
			}
			break;
		}
		
	default:
		UnexpectedToken( Row, Col, NextToken, Lexer );
		break;
	}
	
	RequireToken( Row, Col, ARM_GenSecurity::ParseCellRefText, ARM_GenPricerLexer::RBRACKET, Lexer );
	
	/// result
	ARM_RowColCoords RefCoord(RowRef, ColRef);

	/// test circular reference!
	if( ARM_GenSecurity::itsLook4CircularRef && CircularRefCoords.find( RefCoord ) != CircularRefCoords.end() )
	{
		CC_Ostringstream os;
		os << ARM_USERNAME + ": in " << CellName( Row, Col ) << " with text: " 
			<< CellReducedText( Row, Col ) << " circular reference to " << CellName( RowRef, ColRef );
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}

	return RefCoord;
}


////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ParseCellRef
///	Returns: ARM_ExpNodePtr
///	Action : Parse an CellRef out of the lexer.
////////////////////////////////////////////////////

ARM_ExpNodePtr ARM_GenSecurity::ParseCellRef(
	size_t Row, 
	size_t Col, 
	ARM_GenPricerLexer &Lexer,
	ARM_RefNNode_Manager& refNNodeManager,
	ARM_ModelCallMap& modelCallMap ) const
{
	ARM_RowColCoords ParentCoords( Row, Col );

	/// a reference grammar is of the following type
	///		Symbol [ integer ]
	/// |	Symbol [ i ]
	/// |	Symbol [ i +/- integer ]
	
	/// the caller can not be called
	/// hence it belongs to the CircularRefCoords
	if( ARM_GenSecurity::itsLook4CircularRef )
		refNNodeManager.InsertCircularRefCoords( ParentCoords );

	ARM_RowColCoords ChildCoords = ParseCellRefCoords( Row, Col, Lexer, refNNodeManager.GetCircularRefCoords() );
	size_t RowRef = ChildCoords.itsElem1;
	size_t ColRef = ChildCoords.itsElem2;
	ARM_ExpNodePtr ParentNode;

	/// see if this node has not been already inserted and could be shared!
	ARM_SharedNodeInfo sharedNodeInfo( ChildCoords, RowDate(Row) );

	/// did we find a shared node?
	if( !GetSharedNode(sharedNodeInfo,refNNodeManager,ParentNode,ParentCoords) )
	{
		/// could not find a shared node..
		/// we need to create the corresponding node!
		/// Do not allow any reference to cell in the future
		/// as only ONLY ONLY PV can do that currently
		/// and PV is specifically designed for that!
		if( RowRef > Row )
		{

			CC_Ostringstream os;
			os << ARM_USERNAME + ": in " << CellName( Row, Col ) << " with text: " 
				<< CellReducedText( Row, Col ) << " try to reference to " << CellName( RowRef, ColRef ) 
				<< " which is in the future but do not use PV!";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}

		/// we make a copy of the refCoords to do only 
		/// circular reference lookup on the current node!
		set< ARM_RowColCoords > oldCircularRefCoords = refNNodeManager.GetCircularRefCoords();

		/// parse the cell
		ARM_ExpNodePtr ChildNode = ParseCell( RowRef, ColRef, refNNodeManager, modelCallMap );
		
		/// restore previous circular references
		refNNodeManager.SetCircularRefCoords( oldCircularRefCoords );
		
		bool isAPVNode = false;

		/// create the parent node
#if defined(__GP_SHOW_SHARED_NODE_COORDINATES)
		ParentNode = ARM_ExpNodePtr( new ARM_ExpNodeRef( ChildNode, RowDate(Row), RowDate(RowRef), CellName( RowRef, ColRef ),  ParentCoords, isAPVNode, itsPayModelName ) );
#else
		ParentNode = ARM_ExpNodePtr( new ARM_ExpNodeRef( ChildNode, RowDate(Row), RowDate(RowRef), CellName( RowRef, ColRef ), isAPVNode, itsPayModelName ) );
#endif
		/// register the shared node
		RegisterSharedNode( sharedNodeInfo, ParentNode, RowDate(Row), ChildNode, RowDate(RowRef), ParentCoords, ChildCoords, refNNodeManager );

	}

	return ParentNode;
}



////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ParseDateOrMatu
///	Returns: ARM_ExpNodePtr
///	Action : Parse a date expression, i.e. either a date or a 
//			 reference to one or a maturity!
////////////////////////////////////////////////////

ARM_ExpNodePtr ARM_GenSecurity::ParseDateOrMatu( size_t Row, size_t Col, ARM_GenPricerLexer &Lexer ) const
{
	/// Date Or Matu= 
	///		Maturity 
	/// |	Date

	int	NextToken = Lexer.peekToken();
	Lexer.unPeek();

	if( NextToken == ARM_GenPricerLexer::MATURITY )
	{
		Lexer.nextToken();
		return ARM_ExpNodePtr( new ARM_ExpNodeString( Lexer.text()) );
	}
	else
		return ParseDate( Row, Col, Lexer );
}

////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ParseStringOrCurve
///	Returns: ARM_ExpNodePtr
///	Action : 
////////////////////////////////////////////////////

ARM_ExpNodePtr ARM_GenSecurity::ParseStringOrCurve( size_t Row, size_t Col, 
	ARM_GenPricerLexer& Lexer, 
	ARM_RefNNode_Manager& refNNodeManager,
	ARM_ModelCallMap& modelCallMap ) const
{
	/// Date Or Matu= 
	///		Maturity 
	/// |	Date

	int Token = Lexer.peekToken();
	Lexer.unPeek();
	
	if( Token == ARM_GenPricerLexer::SYMBOL )
	{
		return ParseDateOrExpr( Row, Col, Lexer, refNNodeManager, modelCallMap );
		//return ARM_ExpNodePtr( ParseString( Lexer ) );
	}
	else
		return ARM_ExpNodePtr( ParseString( Lexer ) );
		
}



////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ParseDateOrExpr
///	Returns: ARM_ExpNodePtr
///	Action : Parse a date expression, i.e. either a date or a 
//			 reference to one or parse an expression!
////////////////////////////////////////////////////

ARM_ExpNodePtr ARM_GenSecurity::ParseDateOrExpr( size_t Row, size_t Col, 
	ARM_GenPricerLexer& Lexer, 
	ARM_RefNNode_Manager& refNNodeManager,
	ARM_ModelCallMap& modelCallMap ) const
{
	string NextTokenText = Lexer.peekText();
	int	NextToken = Lexer.peekToken();
	Lexer.unPeek();
	
	if(	   NextToken == ARM_GenPricerLexer::DOUBLE 
		|| NextToken == ARM_GenPricerLexer::INTEGER )
	{
		double JulianDouble = atof( NextTokenText.c_str() );
			
		/// corresponds to 1-Jan-90 in julian date format
		const double MININUMDATE = 2447893.0;
		
		/// test that the date is making sense!
		if( JulianDouble > MININUMDATE )
		{
			Lexer.nextToken();
			return ARM_ExpNodePtr( new ARM_ExpNodeDate( ARM_Date( JulianDouble ) ) );
		}
	};

	/// other cases
	return ParseExpr( Row, Col, Lexer, refNNodeManager, modelCallMap );
}



////////////////////////////////////////////////////
///	Class  : ARM_GenericSecurity 
///	Routine: ParseMaturity
///	Returns: ARM_ExpNodePtr 
///	Action : Parse a maturity!
////////////////////////////////////////////////////

ARM_ExpNodePtr ARM_GenSecurity::ParseMaturity( size_t Row, size_t Col, ARM_GenPricerLexer &Lexer ) const
{
	RequireToken( Row, Col, ARM_GenSecurity::ParseMaturityText, ARM_GenPricerLexer::MATURITY, Lexer );
	return ARM_ExpNodePtr( new ARM_ExpNodeString( Lexer.text()) );
}



////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ParseDate
///	Returns: ARM_ExpNodePtr 
///	Action : Parse a date expression, i.e. either a date or a 
//			 reference to one.
////////////////////////////////////////////////////

ARM_ExpNodePtr ARM_GenSecurity::ParseDate( size_t Row, size_t Col, ARM_GenPricerLexer &Lexer ) const
{
	/// Date= 
	///		DateExpr
	/// |	Double interpreted as a date

	int	NextToken = Lexer.peekToken();
	switch( NextToken )
	{
	
	/// is it a date expr?
	case ARM_GenPricerLexer::DATE :
		{
			Lexer.unPeek();
			return ParseDateExpr( Row, Col, Lexer );
		}
	
	/// or a simple integer or double
	case ARM_GenPricerLexer::DOUBLE : case ARM_GenPricerLexer::INTEGER :
		{

			Lexer.nextToken();
			string Elem( Lexer.text() );
			
			double JulianDouble = atof( Elem.c_str() );
			
			/// corresponds to 1-Jan-90 in julian date format
			const double MININUMDATE = 2447893.0;
		
			/// test that the date is making sense!
			if( JulianDouble< MININUMDATE )
			{
				CC_Ostringstream os;
				os << ARM_USERNAME + ": in " << CellName( Row, Col ) << " with text: " 
					<< CellReducedText( Row, Col ) << " try to use a date that is before 1-Jan-1990"
					<< " which is not permitted!";
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
			}
			else
				return ARM_ExpNodePtr( new ARM_ExpNodeDate( ARM_Date( JulianDouble ) ) );
		}

	default :
		return UnexpectedToken( Row, Col, NextToken, Lexer );
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ParseDateExpr
///	Returns: ARM_ExpNodePtr
///	Action : Parse a date expression, i.e. either a date or a 
//			 reference to one.
////////////////////////////////////////////////////

ARM_ExpNodePtr ARM_GenSecurity::ParseDateExpr(
	size_t Row, 
	size_t Col, 
	ARM_GenPricerLexer &Lexer ) const
{
	/// DateExpr = Date
	/// of the type ddMONTyear like 12JAN2003

	int	NextToken = Lexer.peekToken();

	switch( NextToken )
	{

		/// is it a date?
		case ARM_GenPricerLexer::DATE:
		{
			Lexer.nextToken();
			string DateString = Lexer.text();

			int day, year;
			char monthChar[3];
			sscanf( DateString.c_str(), "%d%3s%d", &day, &monthChar, &year );
			int month = GetNumMonthFromStr( monthChar );
			
			/// Y2K issue
			year = year < 100? year + 2000 : year;

			return ARM_ExpNodePtr( new ARM_ExpNodeDate( ARM_Date( day, month, year ) ) );
		}

		/// other means an unexpectedToken
		default:
			return UnexpectedToken( Row, Col, NextToken, Lexer );
	}	
}


////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: GetSharedNode
///	Returns: bool
///	Action : For a given row col coordinate, and a parent call
///				date, tries to see if it can find a shared node!
///				the boolean return result says whether there already
///				existed a shared node or not!
////////////////////////////////////////////////////

bool ARM_GenSecurity::GetSharedNode( 
	const ARM_SharedNodeInfo& sharedNodeInfo,
	ARM_RefNNode_Manager& refNNodeManager,
	ARM_ExpNodePtr& ParentNode,
	const ARM_RowColCoords& pCoords ) const
{
	/// typedef for shorter type
	typedef ARM_RefNNode_Manager::ARM_SharedNodeMap::iterator ARM_SharedNodeMapIter;

	/// is there already a node corresponding to it?
	ARM_SharedNodeMapIter sharedNode = refNNodeManager.findSharedNode( sharedNodeInfo );

	/// did we find a shared node?
	if( sharedNode != refNNodeManager.GetSharedNodeMapEnd() )
	{
		/// ok the node is shared
		/// get the corresponding refCall
		ParentNode = (*sharedNode).second->AddRerenceNGetRefNode( pCoords );
		return true;
	}
	else
		return false;
}


////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: RegisterSharedNode
///	Returns: void
///	Action : Add to the reference and node manager the shared Node
////////////////////////////////////////////////////
void ARM_GenSecurity::RegisterSharedNode(
	const ARM_SharedNodeInfo& sharedNodeInfo,
	const ARM_ExpNodePtr& ParentNode,
	const ARM_Date& ParentDate,
	const ARM_ExpNodePtr& ChildNode,
	const ARM_Date& ChildDate,
	const ARM_RowColCoords& parentRowCoords,
	const ARM_RowColCoords& childRowCoords,
	ARM_RefNNode_Manager& refNNodeManager ) const
{
	/// we need to add to the current list of shared node this one
	/// create the corresponding cellRefCallPtr
	const int UNASSIGNED = -1;
	ARM_CellRefCallPtr cellRefCall( new ARM_CellRefCall( ParentNode, UNASSIGNED ) );

#ifdef __GP_STRICT_VALIDATION
	/// check for successful insertion
	bool success = refNNodeManager.InsertSharedNode( sharedNodeInfo, cellRefCall );
	if( !success )
	{
		CC_Ostringstream os;
		os << "Could not insert shared node with coord = " << sharedNodeInfo.first.toString()
			<< " with parent node'date = " << sharedNodeInfo.second.toString();
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
#else

	/// no check of succesful insertion
	refNNodeManager.InsertSharedNode( sharedNodeInfo, cellRefCall );
#endif
}



////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ParseArg
///	Returns: ARM_ExpNodePtr
///	Action : Parse an argument of a function given its type
////////////////////////////////////////////////////

ARM_ExpNodePtr ARM_GenSecurity::ParseArg( 
	size_t Row, 
	size_t Col, 
	ARM_GenPricerLexer &Lexer, 
	GramFuncArgType argType,
	ARM_RefNNode_Manager& refNNodeManager,
	ARM_ModelCallMap& modelCallMap,
	bool AddRefNode ) const
{
	/// look 2 tokens ahead to determine whether there is a reference
	Lexer.peekToken();
	int TokenLBracket	= Lexer.peekToken();

	/// is it a reference?
	/// a reference is of type
	///		Symbol [ integer ]
	/// |	Symbol [ i ]
	/// |	Symbol [ i +/- integer ]

	/*try
	{*/
	if( ARM_GenPricerLexer::LBRACKET == TokenLBracket )
	{
		int TokenRBracket = Lexer.peekToken();
		
		// goes up to the end of the 
		while( TokenRBracket != ARM_GenPricerLexer::EOT 
			&& TokenRBracket != ARM_GenPricerLexer::RBRACKET )
			TokenRBracket = Lexer.peekToken();

		if( ARM_GenPricerLexer::EOT == TokenRBracket )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + string(": in ") + CellName( Row, Col ) + " with text: " + CellReducedText( Row, Col ) 
			+ " No matching closing RBracket before end of text in"
			+ CellReducedText( Row, Col ) ); 

		/// a real reference should be of type reference followed by comma or RBRACKET
		int peekToken = Lexer.peekToken();
		if(		ARM_GenPricerLexer::COMMA == peekToken 
			||  ARM_GenPricerLexer::RPAREN == peekToken )
		{
			Lexer.unPeek();

			/// parent coordinates
			ARM_RowColCoords ParentCoords( Row, Col );

			/// the caller can not be called
			/// hence it belongs to the CircularRefCoords
			if( ARM_GenSecurity::itsLook4CircularRef )
				refNNodeManager.InsertCircularRefCoords( ParentCoords );

			/// get the coordinate of the child node
			ARM_RowColCoords ChildCoords = ParseCellRefCoords( Row, Col, Lexer, refNNodeManager.GetCircularRefCoords() );
			size_t RowRef = ChildCoords.itsElem1;
			size_t ColRef = ChildCoords.itsElem2;

			ARM_ExpNodePtr ParentNode;

			/// see if this node has not been already inserted and could be shared!
			ARM_SharedNodeInfo sharedNodeInfo( ChildCoords, RowDate(Row) );

			if( !GetSharedNode(sharedNodeInfo,refNNodeManager,ParentNode,ParentCoords) )
			{
				string ElemRef = itsDealDescription->GetElem( RowRef, ColRef );

				/// We have to make a copy of the lexer so that we can recursively
				/// lex things when we follow references to other cells we have to parse.
				ARM_GenPricerLexer LexerRef;
				LexerRef.setInput( ElemRef );

				/// we make a copy of the ChildCoordss to do only 
				/// circular reference lookup on the current node!
				set< ARM_RowColCoords > oldCircularRefCoords = refNNodeManager.GetCircularRefCoords();
				
				/// use reference counted pointor for exception safety!
				ARM_ExpNodePtr ChildNode;

				/// to flag PV node, we need to detect them
				bool isAPVNode = false;

				/// handle specifically the case of PV
				if( GFAT_FUTUREREF_TYPE == argType )
				{
					isAPVNode = true;

					if( RowRef <= Row )
					{
						CC_Ostringstream os;
						os << ARM_USERNAME << ": in " << CellName( Row, Col )
							<< " with text: " << CellReducedText( Row, Col ) 
							<< " try to reference to " << CellName( RowRef, ColRef ) 
							<<  " in the past or present with pv that is only meaningful for future cashflows!";
						throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
					}

					ChildNode = ParseExpr( RowRef, ColRef, LexerRef, refNNodeManager, modelCallMap );
				}
				else if( GFAT_VECTOR_OR_FUTUREREF_TYPE == argType )
				{

					if( RowRef > Row )
                    {
					    isAPVNode = true;
					    ChildNode = ParseExpr( RowRef, ColRef, LexerRef, refNNodeManager, modelCallMap );
                    }
					else
					    ChildNode = ParseArg( RowRef, ColRef, LexerRef, argType, refNNodeManager, modelCallMap, AddRefNode );
				}
				else
				{
					/// Do not allow any reference to cell in the future
					/// as only ONLY ONLY PV can do that currently
					/// And PV is specifically designed for that!
					if( RowRef > Row )
						throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
						ARM_USERNAME + string( ": in " ) + CellName( Row, Col ) + " with text: " + CellReducedText( Row, Col ) 
						+ " try to reference to " + CellName( RowRef, ColRef ) 
						+ " which is in the future but do not use the type Future states (defined in pv)."
						+ " Please advise" ); 

					ChildNode = ParseArg( RowRef, ColRef, LexerRef, argType, refNNodeManager, modelCallMap, AddRefNode );
				}

				/// restore previous circular references
				refNNodeManager.SetCircularRefCoords( oldCircularRefCoords );

				/// check that this is the end!
				RequireToken( RowRef, ColRef, ARM_GenSecurity::ParseArgText, ARM_GenPricerLexer::EOT, LexerRef );	


				/// create the parent node
#if defined(__GP_SHOW_SHARED_NODE_COORDINATES)
				ParentNode = ARM_ExpNodePtr( new ARM_ExpNodeRef( ChildNode, RowDate(Row), RowDate(RowRef), CellName( RowRef, ColRef ), ParentCoords, isAPVNode, itsPayModelName ) );
#else
				ParentNode = ARM_ExpNodePtr( new ARM_ExpNodeRef( ChildNode, RowDate(Row), RowDate(RowRef), CellName( RowRef, ColRef ), isAPVNode, itsPayModelName ) );
#endif

				/// register the shared node
				RegisterSharedNode( sharedNodeInfo, ParentNode, RowDate(Row), ChildNode, RowDate(RowRef), ParentCoords, ChildCoords, refNNodeManager );
			}

			return ParentNode;
		}
	}

	/// other cases
	Lexer.unPeek();
	return ParseArgNoRef( Row, Col, Lexer, argType, refNNodeManager, modelCallMap, AddRefNode );
}




////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ParseArgNoRef
///	Returns: ARM_ExpNodePtr
///	Action : Parse an argument of a function given its type
///			assuming in most cases that there is no reference!
////////////////////////////////////////////////////

ARM_ExpNodePtr ARM_GenSecurity::ParseArgNoRef( 
	size_t Row, 
	size_t Col, 
	ARM_GenPricerLexer &Lexer, 
	GramFuncArgType argType,
	ARM_RefNNode_Manager& refNNodeManager,
	ARM_ModelCallMap& modelCallMap,
	bool addRefNode) const
{
	ARM_ExpNodePtr childNode, parentNode;	

	switch( argType )
	{
    case GFAT_DOUBLE_TYPE: case GFAT_VECTOR_TYPE: case GFAT_VECTOR_OR_FUTUREREF_TYPE:
		childNode = ParseExpr( Row, Col, Lexer, refNNodeManager, modelCallMap );
		break;
	case GFAT_STRING_TYPE:
		childNode = ParseString( Lexer );
		break;
	case GFAT_STRING_OR_CURVE_TYPE:
		childNode = ParseStringOrCurve( Row, Col, Lexer, refNNodeManager, modelCallMap );
		break;
	case GFAT_DATEORMATU_TYPE:
		childNode = ParseDateOrMatu( Row, Col, Lexer );
		break;
	case GFAT_MODEL_TYPE: case GFAT_MULTITOKENSTRING_TYPE:
		childNode = ParseMultiTokenString( Row, Col, Lexer );
		break;
	case GFAT_MATURITY_TYPE:
		childNode = ParseMaturity( Row, Col, Lexer );
		break;
	case GFAT_DATE_TYPE: 
		childNode = ParseDate( Row, Col, Lexer );
		break;
	case GFAT_VECTOR_OR_CURVE_TYPE:
	case GFAT_DATE_OR_VECTOR_TYPE:
	case GFAT_DATE_OR_DOUBLE_TYPE:
	case GFAT_DATESTRIP_TYPE:
		childNode = ParseDateOrExpr( Row, Col, Lexer, refNNodeManager, modelCallMap );
		break;
	/// other cases are error
	default:
		{
			Lexer.nextToken();
			string Elem( Lexer.text() );

			CC_Ostringstream os;
			os  << ARM_USERNAME + ": in " << CellName( Row, Col ) << " with text " << CellReducedText( Row, Col )
				<< " Type for argument " <<  Elem << " invalid in " << CellReducedText( Row, Col )
				<< " You probably specified incorrectly the table function defintion";

			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}
	}

	if (addRefNode)
	{
		parentNode = InsertFakeNodeRef(childNode,Row,Col);
	}
	else
	{
		parentNode = childNode;
	}

	return parentNode;
}


////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ParseCst
///	Returns: ARM_ExpNodePtr
///	Action : Parse a cst
////////////////////////////////////////////////////

ARM_ExpNodePtr ARM_GenSecurity::ParseCst( size_t Row, size_t Col, ARM_GenPricerLexer& Lexer ) const
{
	Lexer.nextToken();
	string CstName( Lexer.text() );

	/// some validation!
	if(itsCstManager == ARM_CstManagerPtr(NULL) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + string(": in") +
		CellName( Row, Col )  + " with text " + CellReducedText( Row, Col )
		+ " could not find a cst mananger!" );

	if( !itsCstManager->DoesCstNameExist( CstName ) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + string(": in") +
		CellName( Row, Col )  + " with text " + CellReducedText( Row, Col )
		+ "could not find cst " + CstName );

	stringNodePtrMap::const_iterator iter = itsGlobalCsts.find(CstName);

	/// does it already exist?
	if( iter == itsGlobalCsts.end() )
	{
		ARM_ExpNodePtr node = ARM_ExpNodePtr( new ARM_ExpNodeDoubleCst( CstName, itsCstManager ) );
		CC_MUTABLE( ARM_GenSecurity, itsGlobalCsts ).insert( pair<const string, ARM_ExpNodePtr>(CstName,node) );
		return node;
	}
	else
	{
		//// increase the countor! (ugly cast because of the smart pointor framework
		((ARM_ExpNodeDoubleCst&)(*(*iter).second)).IncCountor();
		return (*iter).second;
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: ParseFunction
///	Returns: ARM_ExpNodePtr
///	Action : Parse an function out of the lexer.
////////////////////////////////////////////////////

ARM_ExpNodePtr ARM_GenSecurity::ParseFunction(
	size_t Row, 
	size_t Col, 
	ARM_GenPricerLexer &Lexer,
	ARM_RefNNode_Manager& refNNodeManager,
	ARM_ModelCallMap& modelCallMap,
	bool isAnUnPayFunctionCall ) const
{
	/// symb function =
	///		function name ( series of Arg seperated by commas )
	///
	/// where the function name is checked in the global function table!

	Lexer.nextToken();
	string FunctionName( Lexer.text() );
	int posInCell( Lexer.peekPos() );

	/// unless explicitly said, unpay is not allowed
	if( !isAnUnPayFunctionCall && stringGetUpper(FunctionName) == ARM_GenSecurity::UNPAYSpecialWord )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + string(": in") +
		CellName( Row, Col )  + " with text " + CellReducedText( Row, Col )
		+ string( " " ) + ARM_GenSecurity::UNPAYSpecialWord + " only allowed in the last column as a final operation!" );

	ARM_GramFunction Function( FunctionName, itsPayModelName );
	RequireToken( Row, Col, ARM_GenSecurity::ParseFunctionText, ARM_GenPricerLexer::LPAREN, Lexer );	

	ARM_GramFunction::ArgVector::const_iterator
			Arg = Function.Args().begin(),
			End = Function.Args().end();

	int size = Function.Args().size();

#ifdef __GP_STRICT_VALIDATION
	/// Checks if there is no vaarg unless it is the last the last argument
	while ( Arg != End )
	{
		if( Arg->IsVaArg() && ( (Arg+1) != End )) 
		{
			CC_Ostringstream os;
			os << ARM_USERNAME << ": parsing function in" << CellName( Row, Col ) << " with text: " << CellReducedText( Row, Col ) 
				<< ": no vaarg expected at a place different from the last one. "
				<< "Please advise";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}
		Arg++;
	}
	Arg = Function.Args().begin();
#endif

	/// hold pointor of arguments
	vector< ARM_ExpNodePtr > Args(0);

	/// and init to NULL first
	size_t i = 0;

	/// to detect model call
	/// check whether any of the argument is of type model
	bool ModelCall = false;
	string ModelName;

	while( Arg != End )
	{
		int Token;
		vector< ARM_ExpNodePtr > CurrentArgs(0);
		ARM_ExpNodePtr currentArg(NULL);

		GramFuncArgType argType = Arg->Type();

		/// grab params only if we are not at the end!
		/// and we are not at the end if the next token is 
		/// not the RPAREN!

		Token = Lexer.peekToken();
		Lexer.unPeek();

		if ( !Arg->IsVaArg() )
		{
		/// If current arg is not a VaArg, we have currentArgs = currentArg
			if( ARM_GenPricerLexer::RPAREN != Token && ARM_GenPricerLexer::COMMA != Token  )
				currentArg = ParseArg( Row, Col, Lexer, argType, refNNodeManager, modelCallMap, Arg->AddRefNode() );

			CurrentArgs.push_back( currentArg );
		}
		else
		{
		/// If current arg is a VaArg, we have to loop so as to find all currentArg and fill
		/// CurrentArgs with them. 
			if( ARM_GenPricerLexer::RPAREN != Token )
				do
				{
					Token = Lexer.peekToken();
					Lexer.unPeek();
					if( ARM_GenPricerLexer::RPAREN != Token && ARM_GenPricerLexer::COMMA != Token  )
					{
						currentArg = ParseArg( Row, Col, Lexer, argType, refNNodeManager, modelCallMap, Arg->AddRefNode() );
						CurrentArgs.push_back( currentArg );

						Token = Lexer.peekToken();
						Lexer.unPeek();
						if( Token == ARM_GenPricerLexer::COMMA )
							Token = Lexer.nextToken();
					}
					else 
					{
						CC_Ostringstream os;
						os << ARM_USERNAME << ": parsing function in" << CellName( Row, Col ) << " with text: " << CellReducedText( Row, Col ) 
							<< ": did not expect ,, or )) in va arg parameter."
							<< "Please advise";
						throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
					}
				}
				while ( Token == ARM_GenPricerLexer::COMMA );

		/// The loop is supposed to end with a right parenthesis (vaargs are last params only)
			if (Token != ARM_GenPricerLexer::RPAREN )
			{
				CC_Ostringstream os;
				os << ARM_USERNAME << ": parsing function in" << CellName( Row, Col ) << " with text: " << CellReducedText( Row, Col ) 
					<< ": expected a right parenthesis after the last va arg parameter."
					<< "Please advise";
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
			}
			else
				Lexer.unPeek();
		}

		/// check that required arguments are here!
		// CurrentArgs.size() == 0 does this check in case of required vaarg
		// CurrentArgs[0] == ARM_ExpNodePtr(NULL) stands for the novaarg case. 
		if( Arg->Required() && ( CurrentArgs.size() == 0 || CurrentArgs[0] == ARM_ExpNodePtr(NULL) ) )
		{
			CC_Ostringstream os;
			os << ARM_USERNAME << ": parsing function in" << CellName( Row, Col ) << " with text: " << CellReducedText( Row, Col ) 
				<< ": expected arg " << i << " in function " + FunctionName + " and saw nothing"
				<< "Please advise";

			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}

		/// for non required argument that is missing
		/// take the default which can only be string 
		/// for the time being
		/// Note: does not apply for vaarg arguments
		if( !Arg->Required() && !Arg->IsVaArg() && CurrentArgs[0] == ARM_ExpNodePtr(NULL) )
		{
			switch( Arg->Type() )
			{
			case GFAT_STRING_TYPE:
			case GFAT_STRING_OR_CURVE_TYPE:
            case GFAT_MULTITOKENSTRING_TYPE:
			case GFAT_MATURITY_TYPE:
				CurrentArgs[0] = ARM_ExpNodePtr( new ARM_ExpNodeString( Arg->DefaultValue() ) );
				break;

			case GFAT_DOUBLE_TYPE:
			case GFAT_VECTOR_TYPE:
			case GFAT_VECTOR_OR_CURVE_TYPE: // include the case of double type
			case GFAT_DATE_OR_DOUBLE_TYPE:
			case GFAT_DATE_OR_VECTOR_TYPE:
			case GFAT_DATESTRIP_TYPE:
				CurrentArgs[0] = ARM_ExpNodePtr( new ARM_ExpNodeDateOrDouble( atof( Arg->DefaultValue().c_str() ) ) );
				break;

			default:
				{
					CC_Ostringstream os;
					os  << ARM_USERNAME << ": default arguments can only be string, maturity or double! Pb in " << CellName( Row, Col ) << " with text: " << CellReducedText( Row, Col ) 
						<< "in " <<  FunctionName  << " " << ARM_USERNAME <<": please advise";
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
				}
			}

#ifdef __GP_STRICT_VALIDATION
			if ( CurrentArgs.size() != 1 )
			{
				CC_Ostringstream os;
 				os  << ARM_USERNAME << ": CurrentArgs.size() !=1 "  
					<< "in " <<  FunctionName  << " " << ARM_USERNAME <<": please advise";
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
			}
#endif
		}

		/// loop to push back all current args
		for (size_t j = 0; j < CurrentArgs.size(); ++j)
		{
			Args.push_back( CurrentArgs[j] );

			/// handle model call for the pricingAdviser
			if( GFAT_MODEL_TYPE == argType)
			{
				ModelCall = true;
				/// AAAAAAAAARRRRRRGGGGGG because of smartpointor wrap up
				/// we need to use the devil with &* ... not the best
				/// coding type but currently the only one that works!
				
				/// a model call is either a string or
				/// a reference to a string
				if( !dynamic_cast<ARM_ExpNodeRef*>( &*Args[i] ) )
					ModelName = ((ARM_ExpNodeString*) &*Args[i])->GetString();
				else
					ModelName = ((ARM_ExpNodeString*) &*((ARM_ExpNodeRef*) &*Args[i] )->GetChildNode() )->GetString();
			}
		i++;
		}

		/// can be either a RPAREN or a COMMA!
		/// but if it is not a RPAREN, it SHOULD BE a COMMA!
		Token = Lexer.peekToken();
		Lexer.unPeek();

		if( ARM_GenPricerLexer::RPAREN != Token )
			RequireToken( Row, Col, ARM_GenSecurity::ParseFunctionText, ARM_GenPricerLexer::COMMA, Lexer );

		/// increase
		++Arg; 
	}

#ifdef __GP_STRICT_VALIDATION
	/// sanity check: 
	/// checks if Args has the right size. 
	int k=0;
	Arg = Function.Args().begin();

	while ( Arg != End )
	{
		if( Arg->IsVaArg() ) k++;
		Arg++;
	}
	if ( k==0 )
	{
		// if there is no vaarg, then we must have Args.size() == Function.Args().size()
		if ( Args.size() != Function.Args().size() )
		{
			CC_Ostringstream os;
			os  << ARM_USERNAME << ": Args.size != Fucntion.Args().size() " 
				<< "in " <<  FunctionName  << " " << ARM_USERNAME <<": please advise";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}
	}
	else
	{
		// if there is at least one vaarg, then we must have Args.size() >= Function.Args().size()
		if ( Args.size() < Function.Args().size()-k )
		{
			CC_Ostringstream os;
			os  << ARM_USERNAME << ": Args.size < Fucntion.Args().size() " 
				<< "in " <<  FunctionName  << " " << ARM_USERNAME <<": please advise";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}
	}

#endif

	/// should check for RPAREN
	RequireToken( Row, Col, ARM_GenSecurity::ParseFunctionText, ARM_GenPricerLexer::RPAREN, Lexer );

    ARM_ExpNodeFunc* expNodeFunc = (Function.GetGramNodeBuilderFunc())( Function, Args, RowDate(Row).GetJulian());
	expNodeFunc->setRowInDealDesc( Row );
	expNodeFunc->setColInDealDesc( Col );
	expNodeFunc->setPosInCell( posInCell );

    ARM_ExpNodePtr resultNode( expNodeFunc );

	if(ModelCall)
	{
		ARM_ModelCallMap::iterator modelCall = modelCallMap.find( RowDate(Row) );
		if( modelCall != modelCallMap.end() )
		{
			(*modelCall).second->AddModelCall( RowDate(Row),ModelName, 
				FunctionName, resultNode );
		}
		else
		{
			ARM_ModelCallPtr modelCall( new ARM_ModelCall(  RowDate(Row),ModelName, 
				FunctionName, resultNode ) );
#ifdef __GP_STRICT_VALIDATION

			bool success = modelCallMap.insert( CC_NS(std,pair)< const ARM_Date, ARM_ModelCallPtr >( RowDate(Row), modelCall ) ).second;
			if( !success )
			{
				CC_Ostringstream os;
				os  << ARM_USERNAME << ": could not insert model call for = " << RowDate(Row).toString()
					<< " with function = " << FunctionName;
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
			}
#else
			modelCallMap.insert( CC_NS(std,pair)< const ARM_Date, ARM_ModelCallPtr >( RowDate(Row), modelCall ) );
#endif
		}
	}

	if( Function.FuncName() == "PV" )
		AddPVNodeToPricingAdviser( resultNode ); 

	return resultNode;
}


////////////////////////////////////////////////////
///	Class  : ARM_GenericSecurity 
///	Routine: ParseString
///	Returns: ARM_ExpNodePtr
///	Action : Parse an rate out of the lexer.
////////////////////////////////////////////////////

ARM_ExpNodePtr ARM_GenSecurity::ParseString( ARM_GenPricerLexer &Lexer ) const
{
	Lexer.nextToken();
	return ARM_ExpNodePtr( new ARM_ExpNodeString( Lexer.text() ) );					   
}


////////////////////////////////////////////////////
///	Class  : ARM_GenericSecurity 
///	Routine: ParseMultiTokenString
///	Returns: ARM_ExpNodePtr
///	Action : Parse a multi-token string!
////////////////////////////////////////////////////
ARM_ExpNodePtr ARM_GenSecurity::ParseMultiTokenString(
	size_t Row, 
	size_t Col, 
	ARM_GenPricerLexer &Lexer ) const
{
	string result;
	int currentToken = Lexer.peekToken();
	while( ARM_GenPricerLexer::RPAREN!= currentToken
		&& ARM_GenPricerLexer::COMMA != currentToken 
		&& ARM_GenPricerLexer::EOT   != currentToken )
	{
		Lexer.nextToken();
		result += Lexer.text();
		currentToken = Lexer.peekToken();
	}
	Lexer.unPeek();
		
	return ARM_ExpNodePtr( new ARM_ExpNodeString(result) );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity 
///	Routine: Clone
///	Returns: ARM_Object*
///	Action :
/////////////////////////////////////////////////////////////////
ARM_Object* ARM_GenSecurity::Clone() const
{
	return new ARM_GenSecurity(*this);
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_GenSecurity
///	Routine: toString 
///	Returns: string describing the content of the ARM_GenSecurity
///	Action : 
/////////////////////////////////////////////////////////////////

string ARM_GenSecurity::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;

	os << "ARM_Gensecurity :\n\n";

	os << "Exercise Boundary Reset (if any) after pricing : ";
		os << (itsExercBoundaryResetFlag ? "on" : "off");
	os << "\n\n";

    os << "Payment Model Name : " << itsPayModelName << "\n\n";
	os << "Corresponding Deal Description:\n";
	os << itsDealDescription->toString();

	os << "\n\nCircular Reference      : ";
	os << ( ARM_GenSecurity::itsLook4CircularRef ? "on\n" : "off\n" );
	os << "\n\n\nCorresponding Pricing Adviser:\n";
	os << itsPricingAdviser->toString();

	if( itsPrintParseTree )
	{
		os << "\n\n\n<Parse-Tree>\n";
		int i=0, j=0, rowSize = itsParseTree->GetRowsNb(), colSize = itsParseTree->GetColsNb();
		for( i=0; i<rowSize; ++i )
		{
			for( j=0; j<colSize; ++j )
			{
				os << "CashFlow	       : " << i+1 << "\n";
				os << "Date            : " << RowDate(i+1).toString() << "\n";
				os << "Node Description:\n" << (*itsParseTree)(i,j)->toString( "\t");
			}
		}
		os << "</Parse-Tree>\n";
	}

	if(itsCstManager != ARM_CstManagerPtr(NULL) )
	{
		os << "\n\n" <<itsCstManager->toString( "\t" );
	}
	
	return os.str();
}

void ARM_GenSecurity::AddPVNodeToPricingAdviser( const ARM_ExpNodePtr& PVNode ) const
{ 
	itsPricingAdviser->AddPVNode(PVNode); 
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

