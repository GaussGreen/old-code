/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gramnode.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */

#ifndef _INGPINFRA_GRAMNODE_H
#define _INGPINFRA_GRAMNODE_H

/// this headers has to come first
/// as it contains pragma for stupid warning
#include "gpbase/removeidentifiedwarning.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/countedptr.h"
#include "gpbase/pair.h"

/// kernel
#include <glob/expt.h>
#include <glob/dates.h>

/// current project
#include "typedef.h"
#include "gramfunction.h"
#include "gramfunctorarg.h"
#include "pricingstates.h"
#include "functorop.h"
#include "functordef.h"
#include "gramfunctorargcheck.h"
#include "cstmanager.h"

/// STL
#include <vector>
CC_USING( std::vector )
//CC_USING( std::vector< std::string > )
#include <string>
CC_USING( std::string )
#include <memory>
CC_USING( std::auto_ptr )


CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_CstManager;
class ARM_DealDescription;
class ARM_PricingModel;

////////////////////////////////////////////////////
/// \class ARM_ExpNode
/// \brief
/// this is the abstract class for all parsing nodes used 
/// in the parsing!
/// Real classes can be either:
///		-ARM_ExpNodeRef
///		-ARM_ExpNodeDate
///		-ARM_ExpNodeFactor
///		-ARM_ExpNodeUnaryOp
///		-ARM_ExpNodeString
///		-ARM_ExpNodeBiOp
///		-ARM_ExpNodeFunc
///
///	an expression node gives information about
///		-itsType
///	it can evaluate itsel and returns the information
/////////////////////////////////////////////////////



/// class definition
struct ARM_ExpNode
{
public: 
	enum ResetLevel
	{ ResetLoop, ResetBucket, ResetTotal, ResetCountor };

	/// pure virtual to force redefinition
	virtual ARM_GramFctorArg Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states) = 0;
	
	/// pure virtual to force redefinition
	virtual string toString( const string& indent = "" ) const = 0;
	/// virtual destructor pure virtual to force redefinition
	virtual ~ARM_ExpNode() =0;
	/// function to remove any cached value pure virtual to force redefinition
	virtual bool Reset(ResetLevel partialReset) =0;
	/// function to change Exercise Nodes into Trigger Nodes ; does nothing by default.
	virtual void ChangeExerciseIntoTrigger( ARM_DealDescription& dealDesc ) {};

	/// clone of the node
	virtual ARM_ExpNode* Clone() const = 0;

	/// separator for all ARM_ExpNode
	static const string itsSeparator;
	/// support in place
	static bool ModelDoesNotSupportInPlace( ARM_PricingModel* model );

	/// Accessors to parents and parents' pointers. 
	void AddParent( ARM_ExpNode* expNode ) { itsParents.push_back( expNode ); }
	vector< ARM_ExpNode* > * getParents() { return &itsParents; }
	void AddParentPointer( ARM_ExpNodePtr * expNodePtr ) { itsParentsPointers.push_back( expNodePtr ); }
	vector< ARM_ExpNodePtr* > * getParentsPointers() { return &itsParentsPointers; }

private:
	vector< ARM_ExpNode* > itsParents;
	vector< ARM_ExpNodePtr* > itsParentsPointers;
};

////////////////////////////////////////////////////
///	Class  :ARM_ExpNodePtr
///	Routine: operator<<
///	Returns: void
///	Action : Display the contents of an ExpNode
////////////////////////////////////////////////////

ostream& operator<< (ostream& os, const ARM_ExpNode& expNode );


////////////////////////////////////////////////////
/// \class ARM_ExpNodeRef
/// \brief node representing a reference to another cell
/// of the deal description
/////////////////////////////////////////////////////

class ARM_ExpNodeRef: public ARM_ExpNode 
{
public:
	/// constructor in line for fast access
#if defined(__GP_SHOW_SHARED_NODE_COORDINATES)
	inline ARM_ExpNodeRef( 
		const ARM_ExpNodePtr& childNode, 
		const ARM_Date& parentDate,
		const ARM_Date& childDate,
		const string& childLocationString,
		const ARM_RowColCoords& parentCoords,
		bool isAPVNode,
        const string& payModelName )
		:	itsChildNode(childNode), 
			itsParentDate(parentDate),
			itsValue(NULL), 
			itsIsPreComputed(false),
			itsIsReseted(true), 
			itsHasBeenChangeIntoTrigger(false),
			itsChildLocationString( childLocationString ), 
			itsChildDate(childDate),
			itsCountor(1),
			itsEvalCountor(0),
			itsIsAPVNode(isAPVNode),
			itsParentCoords(1,parentCoords),
            itsPayModelName(payModelName)
	{ itsChildNode->AddParent( this ); itsChildNode->AddParentPointer( &itsChildNode );}
#else
	inline ARM_ExpNodeRef( 
		const ARM_ExpNodePtr& childNode, 
		const ARM_Date& parentDate,
		const ARM_Date& childDate,
		const string& childLocationString,
		bool isAPVNode,
        const string& payModelName )
		:	itsChildNode(childNode), 
			itsParentDate(parentDate),
			itsValue(NULL), 
			itsIsPreComputed(false), 
			itsIsReseted(true),
			itsHasBeenChangeIntoTrigger(false),
			itsChildLocationString( childLocationString ), 
			itsChildDate(childDate),
			itsCountor(1),
			itsEvalCountor(0),
			itsIsAPVNode(isAPVNode),
            itsPayModelName(payModelName)
	{itsChildNode->AddParent( this ); itsChildNode->AddParentPointer( &itsChildNode );}
#endif
	ARM_GramFctorArg EvalAux( ARM_PricingModel* model, const ARM_PricingStatesPtr& states, bool auxiliaryEval=false);
	virtual ARM_GramFctorArg Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states);
	void SetValue( const ARM_GramFctorArg& val, ARM_PricingModel* model, const ARM_PricingStatesPtr& states );
	virtual string toString( const string& indent = "" ) const;

	/// ------ accessors
	inline ARM_Date GetParentDate() const { return itsParentDate; }
	inline ARM_Date GetChildDate() const { return itsChildDate; }
	inline ARM_ExpNodePtr GetChildNode() const { return itsChildNode; }
#if defined(__GP_SHOW_SHARED_NODE_COORDINATES)
	inline vector<ARM_RowColCoords> GetParentCoords() const { return itsParentCoords; }
#endif
	inline string GetChildLocationString() const { return itsChildLocationString; }

	inline bool IsPreComputed() const {return itsIsPreComputed; }
	inline bool IsReseted() const {return itsIsReseted; }
	inline bool IsHasBeenChangeIntoTrigger() const {return itsHasBeenChangeIntoTrigger; }
	inline void IncCountor( const ARM_RowColCoords& pCoords)
	{ 
#if defined(__GP_SHOW_SHARED_NODE_COORDINATES)
		itsParentCoords.push_back(pCoords );
#endif
		++itsCountor; 
	}
	inline int GetCountor() { return itsCountor; }
	virtual void ChangeExerciseIntoTrigger( ARM_DealDescription& dealDesc );
	virtual bool Reset(ResetLevel partialReset);
	virtual ARM_ExpNode* Clone() const { return new ARM_ExpNodeRef(*this); }
	inline ARM_ExpNodePtr& getChildNode() { return itsChildNode; }

private:
	ARM_ExpNodePtr	itsChildNode;
	ARM_Date itsParentDate;
	ARM_GramFctorArg itsValue;
	bool itsIsPreComputed;
	bool itsIsReseted;
	bool itsHasBeenChangeIntoTrigger;
	string itsChildLocationString;
	ARM_Date itsChildDate;
	int itsCountor;
	int itsEvalCountor;
	bool itsIsAPVNode;

    string itsPayModelName;

#if defined(__GP_SHOW_SHARED_NODE_COORDINATES)
	vector<ARM_RowColCoords> itsParentCoords;
#endif
};



////////////////////////////////////////////////////
/// \class	ARM_ExpNodeDateOrDouble
/// \brief	node representing either a date or a double
/////////////////////////////////////////////////////

/// class definition
class ARM_ExpNodeDateOrDouble: public ARM_ExpNode 
{
public:
	ARM_ExpNodeDateOrDouble( double d ) : itsDouble(d) {}
	virtual ARM_GramFctorArg Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states);
	virtual string toString( const string& indent = "" ) const;
	virtual bool Reset(ResetLevel partialReset) { return (partialReset == ResetLoop) ? false : true; }
	virtual ARM_ExpNode* Clone() const { return new ARM_ExpNodeDateOrDouble(*this); }

private:
	double itsDouble;
};



////////////////////////////////////////////////////
/// \class	ARM_ExpNodeDoubleCst
/// \brief	node representing a cst that is given by a cst manager
/////////////////////////////////////////////////////
/// class definition
class ARM_ExpNodeDoubleCst: public ARM_ExpNode 
{
public:
	ARM_ExpNodeDoubleCst( const string& name, ARM_CstManagerPtr cstManager );
	virtual ARM_GramFctorArg Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states);
	virtual string toString( const string& indent = "" ) const;
	virtual bool Reset(ResetLevel partialReset) { return (partialReset == ResetLoop) ? false : true; }
	virtual ARM_ExpNode* Clone() const { return new ARM_ExpNodeDoubleCst(*this); }
	void IncCountor(){ ++itsCountor; }
private:
	ARM_CstManager::iterator itsCstIterator;
	size_t itsCountor;
	size_t itsEvalCountor;
};


////////////////////////////////////////////////////
/// \class	ARM_ExpNodeDate
/// \brief	node for date
/////////////////////////////////////////////////////

class ARM_ExpNodeDate: public ARM_ExpNode 
{
public:
	ARM_ExpNodeDate( const ARM_Date& d ) : itsDate( d ) {};
	virtual ARM_GramFctorArg Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states);
	virtual string toString( const string& indent = "" ) const;
	virtual bool Reset(ResetLevel partialReset) { return (partialReset == ResetLoop) ? false : true; }
	virtual ARM_ExpNode* Clone() const { return new ARM_ExpNodeDate(*this); }

private:
	ARM_Date itsDate;
};



////////////////////////////////////////////////////
/// \class	ARM_ExpNodeFactor
/// \brief	node for factor
/////////////////////////////////////////////////////

class ARM_ExpNodeFactor: public ARM_ExpNode
{
public:
	ARM_ExpNodeFactor( ARM_ExpNodePtr pNode ) : itspNode( pNode ) {itspNode->AddParent( this ); itspNode->AddParentPointer( &itspNode );}
	virtual ARM_GramFctorArg Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states);
	virtual string toString( const string& indent = "" ) const;
	virtual bool Reset(ResetLevel partialReset) { return itspNode->Reset(partialReset); }
	virtual void ChangeExerciseIntoTrigger( ARM_DealDescription& dealDesc );
	virtual ARM_ExpNode* Clone() const { return new ARM_ExpNodeFactor(*this); }

private:
	ARM_ExpNodePtr	itspNode;

};


////////////////////////////////////////////////////
/// \class	ARM_ExpNodeString
/// \brief	node for string
/////////////////////////////////////////////////////

class ARM_ExpNodeString: public ARM_ExpNode
{
public:
	ARM_ExpNodeString( const string& s ) : itsString( s ) {}
	virtual ARM_GramFctorArg Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states);
	virtual string toString( const string& indent = "" ) const;
	inline string GetString() const { return itsString; }
	virtual bool Reset(ResetLevel partialReset) { return (partialReset == ResetLoop) ? false : true; }
	virtual ARM_ExpNode* Clone() const { return new ARM_ExpNodeString(*this); }

private:
	string itsString;
};


////////////////////////////////////////////////////
/// \class	ARM_ExpNodeUnaryOp
/// \brief	node for unary minus and logical not
/////////////////////////////////////////////////////

template <typename UnaryOp>
	class ARM_ExpNodeUnaryOp: public ARM_ExpNode
{
public:
	ARM_ExpNodeUnaryOp(UnaryOp Op, const ARM_ExpNodePtr& pNode ) : itsOp( Op), itspNode( pNode ) {itspNode->AddParent( this ); itspNode->AddParentPointer( &itspNode );} 
	virtual ARM_GramFctorArg Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states);
	virtual string toString( const string& indent = "" ) const;
	virtual bool Reset(ResetLevel partialReset) { return itspNode->Reset(partialReset); }
	virtual void ChangeExerciseIntoTrigger( ARM_DealDescription& dealDesc );
	virtual ARM_ExpNode* Clone() const { return new ARM_ExpNodeUnaryOp<UnaryOp>(*this); }

private:
	ARM_GramFctorArg EvalInPlace( ARM_PricingModel* model, const ARM_PricingStatesPtr& states);
	ARM_GramFctorArg EvalWithClone( ARM_PricingModel* model, const ARM_PricingStatesPtr& states);

	UnaryOp itsOp;
	ARM_ExpNodePtr itspNode;
};


////////////////////////////////////////////////////
/// \class	ARM_ExpNodeBiOp
/// \brief	node for binary operators, either:
///				- Arithmetic ones	: +,-,*,/ mod
/////////////////////////////////////////////////////

template <typename Operator>
	class ARM_ExpNodeBiOp : public ARM_ExpNode
{
public:
	ARM_ExpNodeBiOp( const Operator& op, const ARM_ExpNodePtr& pLeft, const ARM_ExpNodePtr& pRight ) : itsOp(op), itspLeft(pLeft), itspRight(pRight) {itspLeft->AddParent( this );itspRight->AddParent( this ); itspLeft->AddParentPointer( &itspLeft );itspRight->AddParentPointer( &itspRight );}
	virtual ARM_GramFctorArg Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states);
	virtual string toString( const string& indent = "" ) const;
	// With the partial reset flag you are resetting the node
	// asking for.
	virtual bool Reset(ResetLevel partialReset) { bool resetLeft = itspLeft->Reset(partialReset); bool resetRight = itspRight->Reset(partialReset); return resetLeft || resetRight;}
	virtual void ChangeExerciseIntoTrigger( ARM_DealDescription& dealDesc );
	virtual ARM_ExpNode* Clone() const { return new ARM_ExpNodeBiOp<Operator>(*this); }

protected:
	Operator itsOp;
	ARM_ExpNodePtr itspLeft;
	ARM_ExpNodePtr itspRight;
private: 
	ARM_GramFctorArg EvalInPlace( ARM_PricingModel* model, const ARM_PricingStatesPtr& states);
	ARM_GramFctorArg EvalWithClone( ARM_PricingModel* model, const ARM_PricingStatesPtr& states);
};




////////////////////////////////////////////////////
/// \class	ARM_ExpNodeWithDatesBiOp 
/// \brief	same version as ARM_ExpNodeBiOp except that it accepts dates
/////////////////////////////////////////////////////

template <typename Operator>
	class ARM_ExpNodeWithDatesBiOp : public ARM_ExpNodeBiOp<Operator>
{
public:
	ARM_ExpNodeWithDatesBiOp( const Operator& op, const ARM_ExpNodePtr& pLeft, const ARM_ExpNodePtr& pRight ) : ARM_ExpNodeBiOp<Operator>(op, pLeft, pRight ) {itspLeft->AddParent( this );itspRight->AddParent( this );itspLeft->AddParentPointer( &itspLeft );itspRight->AddParentPointer( &itspRight );}
//	virtual ARM_GramFctorArg Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states);
	virtual ARM_ExpNode* Clone() const { return new ARM_ExpNodeWithDatesBiOp<Operator>(*this); }
	
	/// copy and assignement operator
	ARM_ExpNodeWithDatesBiOp( const ARM_ExpNodeWithDatesBiOp<Operator>& rhs ): ARM_ExpNodeBiOp<Operator>(rhs) {}
	ARM_ExpNodeWithDatesBiOp& operator=( const ARM_ExpNodeWithDatesBiOp<Operator>& rhs )
	{
		if( this != &rhs )
			ARM_ExpNodeBiOp::operator=( rhs );
		return *this;
	}
private: 
	ARM_GramFctorArg EvalInPlace( ARM_PricingModel* model, const ARM_PricingStatesPtr& states);
	ARM_GramFctorArg EvalWithClone( ARM_PricingModel* model, const ARM_PricingStatesPtr& states);
};




////////////////////////////////////////////////////
/// \class	ARM_ExpNodeFunc
/// \brief	node representing a function. These functions
///				are defined in the gramfunctiontable
/////////////////////////////////////////////////////

class ARM_ExpNodeFunc: public ARM_ExpNode
{
public:
	ARM_ExpNodeFunc( const ARM_GramFunction& Function, const CC_STL_VECTOR( ARM_ExpNodePtr )& Args, double evalDate);
			
	virtual ARM_GramFctorArg Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states );
	virtual string toString( const string& indent = "" ) const;
	virtual bool Reset(ResetLevel partialReset);
	virtual void ChangeExerciseIntoTrigger( ARM_DealDescription& dealDesc );
	virtual ARM_ExpNode* Clone() const { return new ARM_ExpNodeFunc(*this); }

	inline void setRowInDealDesc( size_t row ) { itsRowInDealDesc=row;}
	inline void setColInDealDesc( size_t col ) { itsColInDealDesc=col;}
	inline size_t getColInDealDesc( ) { return itsColInDealDesc;}
	inline size_t getRowInDealDesc( ) { return itsRowInDealDesc;}
	inline void setPosInCell( size_t pos ) { itsPosInCell = pos; }
	inline size_t getPosInCell( ) { return itsPosInCell; }
	inline double getEvalDate() const { return itsEvalDate; }
	inline ARM_GramFctorPtr& getFctor() { return itsFunction.Fctor(); }
	inline ARM_GramFunction& getFunction() { return itsFunction; }
	inline CC_STL_VECTOR( ARM_ExpNodePtr ) * getArgs() { return &itsArgs; }

	/// specific to Expression Node!
    ARM_NodeInfo GetUsedTimeLags( ARM_PricingModel* mod);

protected:
	ARM_GramFunction itsFunction;				/// function object including the functor
	CC_STL_VECTOR( ARM_ExpNodePtr )  itsArgs;	/// vector on the nodes of the arguments
	double itsEvalDate;
	size_t itsRowInDealDesc, itsColInDealDesc, itsPosInCell;
};
	

////////////////////////////////////////////////////
/// \class	ARM_ExpNodeExerciseFunc
/// \brief	node representing the exercise the first argument (payoff) is
///  evaluated in a forward induction loop and the second argument is 
///  evaluated in a backward induction loop.
/////////////////////////////////////////////////////

class ARM_ExpNodeExerciseFunc: public ARM_ExpNodeFunc
{
public:
	ARM_ExpNodeExerciseFunc( const ARM_GramFunction& Function, const CC_STL_VECTOR( ARM_ExpNodePtr )& Args, double evalDate)
		:	ARM_ExpNodeFunc( Function, Args, evalDate ) {};
	ARM_ExpNodeExerciseFunc( const ARM_ExpNodeExerciseFunc& rhs ): ARM_ExpNodeFunc(rhs ) {}
	ARM_ExpNodeExerciseFunc& operator=(const ARM_ExpNodeExerciseFunc& rhs )
	{
		if( this != &rhs )
			ARM_ExpNodeFunc::operator=( rhs );
		return *this;
	}
	virtual ARM_GramFctorArg Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states );
	virtual bool Reset(ResetLevel resetParent);
	virtual void ChangeExerciseIntoTrigger( ARM_DealDescription& dealDesc );
	virtual ARM_ExpNode* Clone() const { return new ARM_ExpNodeExerciseFunc(*this); }
};



////////////////////////////////////////////////////
/// \class	ARM_ExpNodeLinearPVFunc
/// \brief	node representing the exercise the first argument (payoff) is
///  evaluated in a forward induction loop and the second argument is 
///  evaluated in a backward induction loop.
/////////////////////////////////////////////////////

class ARM_ExpNodeLinearPVFunc: public ARM_ExpNodeFunc
{
public:
	ARM_ExpNodeLinearPVFunc( const ARM_GramFunction& Function, const CC_STL_VECTOR( ARM_ExpNodePtr )& Args, double evalDate)
		:	ARM_ExpNodeFunc( Function, Args, evalDate ) {};
	ARM_ExpNodeLinearPVFunc( const ARM_ExpNodeLinearPVFunc& rhs ): ARM_ExpNodeFunc(rhs ) {}
	ARM_ExpNodeLinearPVFunc& operator=(const ARM_ExpNodeLinearPVFunc& rhs )
	{
		if( this != &rhs )
			ARM_ExpNodeFunc::operator=( rhs );
		return *this;
	}
	virtual ARM_GramFctorArg Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states );
	virtual bool Reset(ResetLevel resetParent);
	virtual ARM_ExpNode* Clone() const { return new ARM_ExpNodeLinearPVFunc(*this); }
};


////////////////////////////////////////////////////
/// \class	ARM_ExpNodeTriggerFunc
/// \brief	node representing the exercise the first argument (payoff) is
///  evaluated in a forward induction loop and the second argument is 
///  evaluated in a backward induction loop.
/////////////////////////////////////////////////////

class ARM_ExpNodeTriggerFunc: public ARM_ExpNodeFunc
{
public:
	ARM_ExpNodeTriggerFunc( const ARM_GramFunction& Function, const CC_STL_VECTOR( ARM_ExpNodePtr )& Args, double evalDate) 
		:	ARM_ExpNodeFunc( Function, Args, evalDate ) {};
	ARM_ExpNodeTriggerFunc( const ARM_ExpNodeTriggerFunc& rhs ): ARM_ExpNodeFunc(rhs ) {}
	ARM_ExpNodeTriggerFunc& operator=(const ARM_ExpNodeTriggerFunc& rhs )
	{
		if( this != &rhs )
			ARM_ExpNodeFunc::operator=( rhs );
		return *this;
	}
	virtual ARM_GramFctorArg Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states );
	virtual bool Reset(ResetLevel resetParent);
	virtual ARM_ExpNode* Clone() const { return new ARM_ExpNodeTriggerFunc(*this); }
};

////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeUnaryOp<T>
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////

template <typename T>
	string ARM_ExpNodeUnaryOp<T>::toString(const string& indent) const
{
	CC_Ostringstream os;
	os << indent << "<Unary-operator-node>\n"
		<< indent << itsSeparator << "Operator: "  << itsOp.toString() << "\n"
		<< indent << itsSeparator << "Node:\n" << itspNode->toString( indent + itsSeparator ) << "\n"
		<< indent << "</Unary-operator-node>\n";
	return os.str();
}



////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeUnaryOp<T>
///	Routine: Eval
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type unary operator
////////////////////////////////////////////////////

template <typename T>
	ARM_GramFctorArg ARM_ExpNodeUnaryOp<T>::Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states )
{
	if( ARM_ExpNode::ModelDoesNotSupportInPlace( model ) )
		return EvalWithClone( model, states );
	else
		return EvalInPlace( model, states );
}

////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeUnaryOp<T>
///	Routine: EvalWithClone
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type unary operator
////////////////////////////////////////////////////

template <typename T>
	ARM_GramFctorArg ARM_ExpNodeUnaryOp<T>::EvalWithClone( ARM_PricingModel* model, const ARM_PricingStatesPtr& states )
{
	ARM_GramFctorArg arg( itspNode->Eval(model, states) );
#ifdef __GP_STRICT_VALIDATION
	string opName = itsOp.toString();
	GPAF_CheckArgType( arg, GFAT_VECTOR_TYPE, opName );
#endif

	if( GFAT_DOUBLE_TYPE == arg.GetType() )
		return ARM_GramFctorArg( itsOp( arg.GetDouble() ) );
	else
	{
		ARM_GP_VectorPtr newvec( static_cast<ARM_GP_Vector*>(arg.GetVector()->Clone()) );
		/// force to use an explicit template argument
		FuncUnaryInPlace( newvec,  itsOp );
		return ARM_GramFctorArg( newvec );
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeUnaryOp<T>
///	Routine: EvalInPlace
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type unary operator
////////////////////////////////////////////////////

template <typename T>
	ARM_GramFctorArg ARM_ExpNodeUnaryOp<T>::EvalInPlace( ARM_PricingModel* model, const ARM_PricingStatesPtr& states )
{
	ARM_GramFctorArg arg( itspNode->Eval(model, states) );
#ifdef __GP_STRICT_VALIDATION
	string opName = itsOp.toString();
	GPAF_CheckArgType( arg, GFAT_VECTOR_TYPE, opName );
#endif

	if( GFAT_DOUBLE_TYPE == arg.GetType() )
		arg.SetDouble( itsOp( arg.GetDouble() ) );
	else
		FuncUnaryInPlace( arg.GetVector(),  itsOp );

	return arg;
}

////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeUnaryOp<T>
///	Routine: ChangeExerciseIntoTrigger
///	Returns: void
///	Action : Changes itsChildNode's Exercises to Trigger
////////////////////////////////////////////////////

template <typename T>
void ARM_ExpNodeUnaryOp<T>::ChangeExerciseIntoTrigger( ARM_DealDescription& dealDesc )
{
	itspNode->ChangeExerciseIntoTrigger( dealDesc );
}


////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeBiOp<T>
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////

template <typename T>
	string ARM_ExpNodeBiOp<T>::toString(const string& indent) const
{
	CC_Ostringstream os;
	os << indent << "<Binary-operator-node>\n"
		<< indent << itsSeparator << "Operator: "  << itsOp.toString() << "\n"
		<< indent << itsSeparator << "Left  Node:\n" << itspLeft->toString( indent + itsSeparator )
		<< indent << itsSeparator << "Right Node:\n" << itspRight->toString( indent + itsSeparator )
		<< indent << "</Binary-operator-node>\n";
	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeBiOp<T>
///	Routine: Eval
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type binary operator
////////////////////////////////////////////////////

template <typename T>
	ARM_GramFctorArg ARM_ExpNodeBiOp<T>::Eval( ARM_PricingModel* model, const ARM_PricingStatesPtr& states )
{
	if( ARM_ExpNode::ModelDoesNotSupportInPlace( model ) )
		return EvalWithClone( model, states );
	else
		return EvalInPlace( model, states );
}


////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeBiOp<T>
///	Routine: EvalInPlace
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type binary operator
////////////////////////////////////////////////////

template <typename T>
	ARM_GramFctorArg ARM_ExpNodeBiOp<T>::EvalInPlace( ARM_PricingModel* model, const ARM_PricingStatesPtr& states )
{
	/// If the nummeth. is fwdbckwd evalNotInPlace ; evalInPlace otherwise
	ARM_GramFctorArg lhs( itspLeft->Eval(model, states) );
	ARM_GramFctorArg rhs( itspRight->Eval(model, states) );

	/// coerce the type
	if( lhs.GetType() == GFAT_DATE_TYPE )
		lhs.SetDouble( lhs.GetDate().GetJulian() );
	if( rhs.GetType() == GFAT_DATE_TYPE )
		rhs.SetDouble( rhs.GetDate().GetJulian() );

	GPAF_CheckArgType( lhs, GFAT_VECTOR_TYPE, "ARM_ExpNodeBiOp<T>" );
	GPAF_CheckArgType( rhs, GFAT_VECTOR_TYPE, "ARM_ExpNodeBiOp<T>" );

	/// handle all the cases
	if(	   lhs.GetType() == GFAT_VECTOR_TYPE 
		&& rhs.GetType() == GFAT_VECTOR_TYPE )
	{	

#ifdef __GP_STRICT_VALIDATION
		string opName = itsOp.toString();
		vector<ARM_GramFctorArg> arg;
		arg.reserve(2);
		arg.push_back(lhs);
		arg.push_back(rhs);
		GPAF_CheckArgVecSameSize( arg, opName );
#endif
		FuncBinaryInPlace( lhs.GetVector(), rhs.GetVector(), itsOp );
	}
	else
		if(    lhs.GetType() == GFAT_DOUBLE_TYPE 
			&& rhs.GetType() == GFAT_VECTOR_TYPE )
		{
			FuncUnaryInPlace( rhs.GetVector(), CC_NS( std, bind1st )( itsOp, lhs.GetDouble() ) );
			lhs = rhs;
		}
		else
			if(	   lhs.GetType() == GFAT_VECTOR_TYPE 
				&& rhs.GetType() == GFAT_DOUBLE_TYPE )
			{
				FuncUnaryInPlace( lhs.GetVector(), CC_NS( std, bind2nd )( itsOp, rhs.GetDouble() ) );
			}
			else
				if( lhs.GetType() == GFAT_DOUBLE_TYPE 
					&& rhs.GetType() == GFAT_DOUBLE_TYPE )
				{
					/// use the functor easily
					lhs.SetDouble( itsOp( lhs.GetDouble(), rhs.GetDouble() ) );
				}
	return lhs;
}

////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeBiOp<T>
///	Routine: EvalWithClone
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type binary operator
////////////////////////////////////////////////////

template <typename T>
	ARM_GramFctorArg ARM_ExpNodeBiOp<T>::EvalWithClone( ARM_PricingModel* model, const ARM_PricingStatesPtr& states )
{
	ARM_GramFctorArg lhs( itspLeft->Eval(model, states) );
	ARM_GramFctorArg rhs( itspRight->Eval(model, states) );

#ifdef __GP_STRICT_VALIDATION
	string opName = itsOp.toString();
	GPAF_CheckArgType( lhs, GFAT_VECTOR_TYPE, opName  );
	GPAF_CheckArgType( rhs, GFAT_VECTOR_TYPE, opName  );
#endif

	/// handle all the cases
	if(	   lhs.GetType() == GFAT_VECTOR_TYPE 
		&& rhs.GetType() == GFAT_VECTOR_TYPE )
	{	

#ifdef __GP_STRICT_VALIDATION
		vector<ARM_GramFctorArg> arg;
		arg.reserve(2);
		arg.push_back(lhs);
		arg.push_back(rhs);
		GPAF_CheckArgVecSameSize( arg, opName );
#endif
		ARM_GP_VectorPtr newvec( ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>( lhs.GetVector()->Clone() ) ) );
		/// force to use explicit template argument
		FuncBinaryInPlace( newvec, rhs.GetVector(), itsOp );
		return ARM_GramFctorArg( newvec );
	}
	else
		if(    lhs.GetType() == GFAT_DOUBLE_TYPE 
			&& rhs.GetType() == GFAT_VECTOR_TYPE )
		{
			/// because the return type is really a pointor
			/// we are forced to use our own version of bind1st that 
			/// operates on pointor hence the name
			ARM_VectorPtr newvec(ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>(rhs.GetVector()->Clone())));
			FuncUnaryInPlace( newvec, CC_NS( std, bind1st )( itsOp, lhs.GetDouble() ) );
			return ARM_GramFctorArg( newvec );
		}
		else
			if(	   lhs.GetType() == GFAT_VECTOR_TYPE 
				&& rhs.GetType() == GFAT_DOUBLE_TYPE )
			{
				/// because the return type is really a pointor
				/// we are forced to use our own version of bind2nd that 
				/// operates on pointor hence the name
				ARM_VectorPtr newvec( ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>(lhs.GetVector()->Clone() ) ) );
				FuncUnaryInPlace( newvec, CC_NS( std, bind2nd )( itsOp, rhs.GetDouble() ) );
				return ARM_GramFctorArg( newvec );
			}
			else
				if( lhs.GetType() == GFAT_DOUBLE_TYPE 
					&& rhs.GetType() == GFAT_DOUBLE_TYPE )
				{
					/// use the functor easily
					return ARM_GramFctorArg( itsOp( lhs.GetDouble(), rhs.GetDouble() ) );
				}
	return lhs;
}


////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeWithDatesBiOp<T>
///	Routine: EvalWithClone
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type binary operator
////////////////////////////////////////////////////

template <typename T>
	ARM_GramFctorArg ARM_ExpNodeWithDatesBiOp<T>::EvalWithClone( ARM_PricingModel* model, const ARM_PricingStatesPtr& states )
{
	ARM_GramFctorArg lhs( itspLeft->Eval(model, states) );
	ARM_GramFctorArg rhs( itspRight->Eval(model, states) );

#ifdef __GP_STRICT_VALIDATION
	string opName = itsOp.toString();
	GPAF_CheckArgType( lhs, GFAT_DATE_VECTOR_OR_DOUBLE_TYPE, opName );
	GPAF_CheckArgType( rhs, GFAT_DATE_VECTOR_OR_DOUBLE_TYPE, opName );
#endif

	/// handle all the cases
	if(	   lhs.GetType() == GFAT_VECTOR_TYPE 
		&& rhs.GetType() == GFAT_VECTOR_TYPE )
	{	

#ifdef __GP_STRICT_VALIDATION
		vector<ARM_GramFctorArg> arg;
		arg.reserve(2);
		arg.push_back(lhs);
		arg.push_back(rhs);
		GPAF_CheckArgVecSameSize( arg, opName );
#endif
		ARM_GP_VectorPtr newvec( ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>(lhs.GetVector()->Clone() ) ) );
		/// force to use explicit template argument
		FuncBinaryInPlace( newvec, rhs.GetVector(), itsOp );
		return ARM_GramFctorArg( newvec );
	}
	else
		if(    lhs.GetType() == GFAT_DOUBLE_TYPE 
			&& rhs.GetType() == GFAT_VECTOR_TYPE )
		{
			/// because the return type is really a pointor
			/// we are forced to use our own version of bind1st that 
			/// operates on pointor hence the name
			ARM_GP_VectorPtr newvec(  ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>(rhs.GetVector()->Clone() ) ) );
			FuncUnaryInPlace( newvec, CC_NS( std, bind1st )( itsOp, lhs.GetDouble() ) );
			return ARM_GramFctorArg( newvec );
		}
		else
			if(	   lhs.GetType() == GFAT_VECTOR_TYPE 
				&& rhs.GetType() == GFAT_DOUBLE_TYPE )
			{
				/// because the return type is really a pointor
				/// we are forced to use our own version of bind2nd that 
				/// operates on pointor hence the name
				ARM_GP_VectorPtr newvec( ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>(lhs.GetVector()->Clone() ) ) );
				FuncUnaryInPlace( newvec, CC_NS( std, bind2nd )( itsOp, rhs.GetDouble() ) );
				return ARM_GramFctorArg( newvec );
			}
			else
				if(    lhs.GetType() == GFAT_DOUBLE_TYPE 
					&& rhs.GetType() == GFAT_DOUBLE_TYPE )
					/// use the functor easily
					return ARM_GramFctorArg( itsOp( lhs.GetDouble(), rhs.GetDouble() ) );
				else 
				/////////////////
				/// date section
				/////////////////
					if(    lhs.GetType() == GFAT_DATE_TYPE 
						&& rhs.GetType() == GFAT_DATE_TYPE )
						/// use the functor easily
						return ARM_GramFctorArg( itsOp( lhs.GetDate().GetJulian(), rhs.GetDate().GetJulian() ) );
					else 
						if(	   lhs.GetType() == GFAT_VECTOR_TYPE 
							&& rhs.GetType() == GFAT_DATE_TYPE )
						{
							/// because the return type is really a pointor
							/// we are forced to use our own version of bind2nd that 
							/// operates on pointor hence the name
							ARM_GP_VectorPtr newvec( ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>(lhs.GetVector()->Clone() ) ) );
							FuncUnaryInPlace( newvec, CC_NS( std, bind2nd )( itsOp, rhs.GetDate().GetJulian() ) );
							return ARM_GramFctorArg( newvec );
						}
						else
							if(    lhs.GetType() == GFAT_DATE_TYPE 
								&& rhs.GetType() == GFAT_VECTOR_TYPE )
							{
								/// because the return type is really a pointor
								/// we are forced to use our own version of bind1st that 
								/// operates on pointor hence the name
								ARM_GP_VectorPtr newvec( ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*>( rhs.GetVector()->Clone() ) ) );
								FuncUnaryInPlace( newvec, CC_NS( std, bind1st )( itsOp, lhs.GetDate().GetJulian() ) );
								return ARM_GramFctorArg( newvec );
							}
							else
								throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "dates can only be compared with dates! not double, not vector!" );
}

////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeWithDatesBiOp<T>
///	Routine: EvalInPlace
///	Returns: ARM_GramFctorArg
///	Action : Evaluate a node of the type binary operator
////////////////////////////////////////////////////

template <typename T>
	ARM_GramFctorArg ARM_ExpNodeWithDatesBiOp<T>::EvalInPlace( ARM_PricingModel* model, const ARM_PricingStatesPtr& states )
{
	ARM_GramFctorArg lhs( itspLeft->Eval(model, states) );
	ARM_GramFctorArg rhs( itspRight->Eval(model, states) );

#ifdef __GP_STRICT_VALIDATION
	string opName = itsOp.toString();
	GPAF_CheckArgType( lhs, GFAT_DATE_VECTOR_OR_DOUBLE_TYPE, opName );
	GPAF_CheckArgType( rhs, GFAT_DATE_VECTOR_OR_DOUBLE_TYPE, opName );
#endif

	/// handle all the cases
	if(	   lhs.GetType() == GFAT_VECTOR_TYPE 
		&& rhs.GetType() == GFAT_VECTOR_TYPE )
	{	

#ifdef __GP_STRICT_VALIDATION
		vector<ARM_GramFctorArg> arg;
		arg.reserve(2);
		arg.push_back(lhs);
		arg.push_back(rhs);
		GPAF_CheckArgVecSameSize( arg, opName );
#endif
		/// force to use explicit template argument
		FuncBinaryInPlace( lhs.GetVector(), rhs.GetVector(), itsOp );
	}
	else
		if(    lhs.GetType() == GFAT_DOUBLE_TYPE 
			&& rhs.GetType() == GFAT_VECTOR_TYPE )
		{
			/// because the return type is really a pointor
			/// we are forced to use our own version of bind1st that 
			/// operates on pointor hence the name
			FuncUnaryInPlace( rhs.GetVector(), CC_NS( std, bind1st )( itsOp, lhs.GetDouble() ) );
			lhs = rhs;
		}
		else
			if(	   lhs.GetType() == GFAT_VECTOR_TYPE 
				&& rhs.GetType() == GFAT_DOUBLE_TYPE )
			{
				/// because the return type is really a pointor
				/// we are forced to use our own version of bind2nd that 
				/// operates on pointor hence the name
				FuncUnaryInPlace( lhs.GetVector(), CC_NS( std, bind2nd )( itsOp, rhs.GetDouble() ) );
			}
			else
				if(    lhs.GetType() == GFAT_DOUBLE_TYPE 
					&& rhs.GetType() == GFAT_DOUBLE_TYPE )
					lhs.SetDouble( itsOp( lhs.GetDouble(), rhs.GetDouble() ) );
				else 
				/////////////////
				/// date section
				/////////////////
					if(    lhs.GetType() == GFAT_DATE_TYPE 
						&& rhs.GetType() == GFAT_DATE_TYPE )
						lhs.SetDate( itsOp( lhs.GetDouble(), rhs.GetDouble() ) );
					else 
						if(	   lhs.GetType() == GFAT_VECTOR_TYPE 
							&& rhs.GetType() == GFAT_DATE_TYPE )
						{
							/// because the return type is really a pointor
							/// we are forced to use our own version of bind2nd that 
							/// operates on pointor hence the name
							FuncUnaryInPlace( lhs.GetVector(), CC_NS( std, bind2nd )( itsOp, rhs.GetDate().GetJulian() ) );
						}
						else
							if(    lhs.GetType() == GFAT_DATE_TYPE 
								&& rhs.GetType() == GFAT_VECTOR_TYPE )
							{
								/// because the return type is really a pointor
								/// we are forced to use our own version of bind1st that 
								/// operates on pointor hence the name
								FuncUnaryInPlace( rhs.GetVector(), CC_NS( std, bind1st )( itsOp, lhs.GetDate().GetJulian() ) );
								lhs = rhs;
							}
							else
								throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "dates can only be compared with dates! not double, not vector!" );
	return lhs;

}

////////////////////////////////////////////////////
///	Class  : ARM_ExpNodeUnaryOp<T>
///	Routine: ChangeExerciseIntoTrigger
///	Returns: void
///	Action : Changes itsChildNode's Exercises to Trigger
////////////////////////////////////////////////////

template <typename T>
void ARM_ExpNodeBiOp<T>::ChangeExerciseIntoTrigger( ARM_DealDescription& dealDesc )
{
	itspLeft->ChangeExerciseIntoTrigger( dealDesc );
	itspRight->ChangeExerciseIntoTrigger( dealDesc );
}


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

