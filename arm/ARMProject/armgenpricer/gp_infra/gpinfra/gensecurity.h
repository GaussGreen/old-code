/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gensecurity.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPINFRA_GENSECURITY_H
#define _INGPINFRA_GENSECURITY_H

/// this headers has to come first
/// as it contains pragma for stupid warning
#include "gpbase/removeidentifiedwarning.h"

#include "gpbase/port.h"
#include "gpbase/gpvector.h"
#include "gpbase/rootobject.h"
#include "gpbase/gpmatrix.h"
#include "typedef.h"

/// ARM Kernel
#include <glob/dates.h>

/// STL
#include <vector>
#include <string>
CC_USING_NS( std, string )
#include <map>
CC_USING_NS( std, map )
CC_USING_NS( std, less )
#include <set>
CC_USING_NS( std, set)

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_RefNNode_Manager;
class ARM_GenPricerLexer;
enum GramFuncArgType;

////////////////////////////////////////////////////
/// \class ARM_GenSecurity 
/// \brief
/// This generic security is responsible for getting the
/// description of the deal as a vector of string
/// and parsing it
/////////////////////////////////////////////////////

class ARM_GenSecurity : public ARM_RootObject
{
private:
	/// pointor to the deal description
	/// using reference counting hence auto-cleanup
	ARM_DealDescriptionPtr itsDealDescription;

	/// col names and the coresponding map for fast lookup!
	CC_STL_VECTOR( string ) itsColNames;
	map< string, int, less< string > > itsColNamesMap;

	/// vector of the date in the leftmost column
	ARM_GP_VectorPtr itsRowDates;

	/// parse tree corresponding to the cashflows of the last column
	ARM_GP_NodeMatrixPtr  itsParseTree;

	/// for checking various things and taking care of the pricing part
	ARM_PricingAdviserPtr itsPricingAdviser;

	/// flag to say whether to print the parse tree or not!
	bool itsPrintParseTree;

	/// flag to say whether to reset ExerciseBoundaries after pricing or not!
	bool itsExercBoundaryResetFlag;

    /// flag to say if IntermediatePayoffs and Snapshots must be computed or not
    bool itsOtherPayoffsFlag;
	// flag to say if the Intermediate Values must be computed or not
	bool itsIVFlag;

	/// flag to say whether we look for circular references
	static bool itsLook4CircularRef;


	/// Get the date of a row.
	ARM_Date RowDate( size_t Row ) const;

	/// typedef for modelCall
	typedef CC_NS(std,map)<ARM_Date,ARM_ModelCallPtr,less<ARM_Date> > ARM_ModelCallMap;

	/// Format a cell name.  Mostly for error generation.
	CC_NS( std, string ) CellName( size_t Row, size_t Col ) const;
	/// reduced text is similar to cell name except that it truncates the message
	/// used mostly for excel error message since they are limited to 256 char!
	CC_NS( std, string ) CellReducedText( size_t Row, size_t Col ) const;

	/// Helper subroutines for parsing.
	int RequireToken( size_t Row, size_t Col, const string& NodeName, 
		int Token, ARM_GenPricerLexer& Lexer ) const; 
	ARM_ExpNodePtr UnexpectedToken( size_t Row, size_t Col, int Token, ARM_GenPricerLexer&Lexer ) const;

	/// to get CellRef Coordinates
	ARM_RowColCoords ParseCellRefCoords( size_t Row, size_t Col, ARM_GenPricerLexer &Lexer, set< ARM_RowColCoords >& CircularRefCoords ) const;

	/// Parse a cell of unknown type.  Grab the value from itsParseTimeRows.
	ARM_ExpNodePtr ParseCell( size_t Row, size_t Col, ARM_RefNNode_Manager& refNNodeManager, ARM_ModelCallMap& modelCallMap ) const;

	/// Starting point of the parser, given a string.
	ARM_ExpNodePtr ParseCell( size_t Row, size_t Col, const string &ExpressionString, ARM_RefNNode_Manager& refNNodeManager, ARM_ModelCallMap& modelCallMap ) const;

	/// Parse various components of the grammar.
	ARM_ExpNodePtr ParseTerminalKeyAndOtherType( size_t Row, size_t Col, const string& ExpressionString, ARM_RefNNode_Manager& refNNodeManager, ARM_ModelCallMap& modelCallMap ) const;
	ARM_ExpNodePtr ParseExpr( size_t Row, size_t Col, ARM_GenPricerLexer &Lexer, ARM_RefNNode_Manager& refNNodeManager, ARM_ModelCallMap& modelCallMap ) const;
	ARM_ExpNodePtr ParseProduct( size_t Row, size_t Col, ARM_GenPricerLexer &Lexer, ARM_RefNNode_Manager& refNNodeManager, ARM_ModelCallMap& modelCallMap ) const;
	ARM_ExpNodePtr ParseSum( size_t Row, size_t Col, ARM_GenPricerLexer &Lexer, ARM_RefNNode_Manager& refNNodeManager, ARM_ModelCallMap& modelCallMap ) const;
	ARM_ExpNodePtr ParseTerm(	size_t Row, size_t Col, ARM_GenPricerLexer &Lexer, ARM_RefNNode_Manager& refNNodeManager, ARM_ModelCallMap& modelCallMap ) const;
	ARM_ExpNodePtr ParseCellRef( size_t Row, size_t Col, ARM_GenPricerLexer &Lexer, ARM_RefNNode_Manager& refNNodeManager, ARM_ModelCallMap& modelCallMap ) const;
	ARM_ExpNodePtr ParseFunction( size_t Row, size_t Col, ARM_GenPricerLexer &Lexer, ARM_RefNNode_Manager& refNNodeManager, ARM_ModelCallMap& modelCallMap, bool isAnUnPVFunctionCall ) const;
	ARM_ExpNodePtr ParseArg( size_t Row, size_t Col, ARM_GenPricerLexer &Lexer, GramFuncArgType argType, ARM_RefNNode_Manager& refNNodeManager, ARM_ModelCallMap& modelCallMap, bool addRefNode ) const;
	ARM_ExpNodePtr ParseArgNoRef( size_t Row, size_t Col, ARM_GenPricerLexer &Lexer, GramFuncArgType argType, ARM_RefNNode_Manager& refNNodeManager, ARM_ModelCallMap& modelCallMap, bool addRefNode ) const;

	// Inser a fake ref node
	ARM_ExpNodePtr InsertFakeNodeRef( const ARM_ExpNodePtr& childNode, int Row, int Col ) const;

	/// terminal nodes!
	ARM_ExpNodePtr ParseDate( size_t Row, size_t Col, ARM_GenPricerLexer &Lexer ) const;
	ARM_ExpNodePtr ParseDateOrMatu( size_t Row, size_t Col, ARM_GenPricerLexer &Lexer ) const;
	ARM_ExpNodePtr ParseDateExpr( size_t Row, size_t Col, ARM_GenPricerLexer &Lexer ) const;
	ARM_ExpNodePtr ParseString( ARM_GenPricerLexer &Lexer ) const;  
	ARM_ExpNodePtr ParseMaturity( size_t Row, size_t Col, ARM_GenPricerLexer &Lexer ) const ;
	ARM_ExpNodePtr ParseMultiTokenString( size_t Row, size_t Col, ARM_GenPricerLexer &Lexer ) const;
	ARM_ExpNodePtr ParseDateOrExpr(size_t Row, size_t Col, ARM_GenPricerLexer& Lexer, ARM_RefNNode_Manager& refNNodeManager, ARM_ModelCallMap& modelCallMap ) const;
	ARM_ExpNodePtr ParseCst( size_t Row, size_t Col, ARM_GenPricerLexer& Lexer) const;
	ARM_ExpNodePtr ParseStringOrCurve(size_t Row, size_t Col, ARM_GenPricerLexer& Lexer, ARM_RefNNode_Manager& refNNodeManager, ARM_ModelCallMap& modelCallMap ) const;
	
	/// shared node manangement
	bool GetSharedNode(  const ARM_SharedNodeInfo& sharedNodeInfo, 
		ARM_RefNNode_Manager& refNNodeManager, 
		ARM_ExpNodePtr& ParentNode,
		const ARM_RowColCoords& pCoords ) const;
	void RegisterSharedNode( const ARM_SharedNodeInfo& sharedNodeInfo, const ARM_ExpNodePtr& ParentNode, const ARM_Date& ParentDate, 
		const ARM_ExpNodePtr& ChildNode, const ARM_Date& ChildDate, 
		const ARM_RowColCoords& ChildCoords, const ARM_RowColCoords& ParentCoords,
		ARM_RefNNode_Manager& refNNodeManager ) const;

	/// timestamp to avoid obsolete xlls
	bool IsVersionOutOfDate() const;

	/// name of the discounting currency
    string			itsPayModelName;

	/// the cst manager
	ARM_CstManagerPtr itsCstManager;

	/// vector of all the global csts
	typedef map<string,ARM_ExpNodePtr, less<string> > stringNodePtrMap;
	CC_IS_MUTABLE stringNodePtrMap itsGlobalCsts;
	string toStringCommon(const string& indent, const string& nextIndent, bool detailMode) const;

	static const string ParseCellRefText;
	static const string ParseCellText;
	static const string ParseTermText;	
	static const string ParseFunctionText;
	static const string ParseMaturityText;
	static const string ParseArgText;
	static const string UNPAYSpecialWord;
	static const string ExerciseSpecialWord;
	static const string TriggerSpecialWord;

public:
	/// this is used to get the type of the text from excel
	ARM_GenSecurity( ARM_DealDescriptionPtr DealDescription, const string& payModelName="", ARM_CstManagerPtr cstManager = ARM_CstManagerPtr(NULL), bool ExercBoundaryResetFlag = true, bool otherPayoffsFlag = true, bool ivFlag =true );
	ARM_GenSecurity( const ARM_GenSecurity& rhs );
	virtual ~ARM_GenSecurity();
	ARM_GenSecurity& operator=( const ARM_GenSecurity& rhs );

	/// inline for fast access
	inline ARM_PricingAdviserPtr GetPricingAdviser() const { return itsPricingAdviser; }

	/// standard Root Object support
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="",const string& nextIndent="") const;
	string detailledString(const string& indent, const string& nextIndent ) const;
	virtual ARM_CLASS_NAME GetRootName() { return ARM_GENSECURITY; }

	inline void SetParseTreePrint( bool val ) { itsPrintParseTree=val; }
	static void SetLookForCircularRef( bool val ) { itsLook4CircularRef=val; }

	/// method to get parse tree corresponding to the cashflows of the last column
	inline ARM_GP_NodeMatrixPtr getParseTree( void ) { return itsParseTree; }

	/// method to add a PVNode to PricingAdviser
	void AddPVNodeToPricingAdviser( const ARM_ExpNodePtr& PVNode ) const;

	/// beware that if you want to modify the deal description, you will have to recreate a generic security
	/// as the parse tree is built once and for all.
    inline const ARM_DealDescription& GetDealDescription() const {return *itsDealDescription;}
	inline ARM_DealDescription& GetDealDescription() {return *itsDealDescription;}

	inline const ARM_DealDescriptionPtr& GetDealDescriptionPtr() {return itsDealDescription;}

    inline const string& GetPayModelName() const {return itsPayModelName;}
    inline void SetPayModelName(const string& payModelName) {itsPayModelName=payModelName;}

	inline bool GetExercBoundaryResetFlag() const { return itsExercBoundaryResetFlag; }
	inline void SetExercBoundaryResetFlag( bool ExercBoundaryResetFlag ) { itsExercBoundaryResetFlag = ExercBoundaryResetFlag; }

	inline bool GetOtherPayoffsFlag() const { return itsOtherPayoffsFlag; }
	inline void SetOtherPayoffsFlag( bool otherPayoffsFlag ) { itsOtherPayoffsFlag = otherPayoffsFlag; }

	inline bool GetIVFlag() const { return itsIVFlag; }
	inline void SetIVFlag( bool IVFlag ) { itsIVFlag = IVFlag; }

	inline ARM_CstManagerPtr GetCstManager() const { return itsCstManager; }
	inline void SetCstManager( const ARM_CstManagerPtr& CstManager ) { itsCstManager = CstManager; }

    virtual string ExportShortName() const { return "LGPGS";}

	const ARM_VectorPtr& GetRowDates() const { return itsRowDates; };

};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

