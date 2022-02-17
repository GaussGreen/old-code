/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gramfunction.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPINFRA_GRAMFUNCTION_H
#define _INGPINFRA_GRAMFUNCTION_H 

/// this headers has to come first
/// as it contains pragma for stupid warning
#include "gpbase/removeidentifiedwarning.h"


#include "typedef.h"
#include "gpbase/port.h"
#include "gpbase/rootobject.h"

#include "typedef.h"
#include "gramfunctionarg.h"

#include <map>
CC_USING_NS( std, map )
CC_USING_NS( std, less )

#include <string>
CC_USING_NS( std, string )


CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
struct ARM_GramHelper;
class ARM_ExpNodeFunc;


/////////////////////////////////////////////////////////////////
/// \class	ARM_GramFunction
/// \brief	this is a class that makes the indirection to the 
/// appropriate function
/////////////////////////////////////////////////////////////////

class ARM_GramFunction : public ARM_RootObject
{
public:
	typedef ARM_ExpNodeFunc* (*GramNodeBuilderFunc)( const ARM_GramFunction& Function, const CC_STL_VECTOR( ARM_ExpNodePtr )& Args, double evalDate );

// FIXMEFRED: mig.vc8 (30/05/2007 16:48:14):makes it public
	/// to make declaration shorter
	typedef ARM_GramFunctionArgVector ArgVector;

private:

	/////////////////////////////////////////////////////////////////
	/// \class	ARM_GramFunction::FunctionDeclaration
	/// \brief	A nested class to make things easy !
	/////////////////////////////////////////////////////////////////
	class Declaration
	{
	private:
		/// name of the function
		string itsFuncName;	
		
		/// type of the return
		GramFuncArgType itsReturnType;
		
		/// its description
		string itsDescription;

		/// vector of arguments
		ArgVector itsArgs;

		ARM_GramFctor* itsFctor;

		/// category like excel!
		string itsCategory;

		GramNodeBuilderFunc itsGramNodeBuilderFunc;
		
	public:
		/// default constructor
		Declaration();
		
		/// real constructor
		Declaration( const string& FuncName,
			GramFuncArgType ReturnType, 
			const string& Description, 
			const ArgVector& Args,
			const string& Category,
			ARM_GramFctor* Fctor,
			GramNodeBuilderFunc gramNodeBuilderFunc );
		
		/// copy constructor
		Declaration( const Declaration& rhs );
		
		/// destructor
		~Declaration();
		
		/// assignment operator
		Declaration &operator=( const Declaration& rhs );
		
		/// for easy debugging
		virtual string toString(const string& indent="", const string& nextIndent="") const;
		
		/// accessors inline as almotst nothing
		inline string FuncName() const { return itsFuncName; }
		inline string Description()	const { return itsDescription; }
		inline const ArgVector& Args() const { return itsArgs; }
		inline GramFuncArgType ReturnType()	const { return itsReturnType; }
		inline ARM_GramFctor* Fctor() const { return itsFctor; }
		inline GramNodeBuilderFunc GetGramNodeBuilderFunc() const { return itsGramNodeBuilderFunc; }
	};

	/// to allow ARM_GramHelper to use a Declaration
	friend struct ARM_GramHelper;

	/// typedef for easy naming
	typedef map< string , void*, less< string > >  MapStringToVoid;

	/// pointor on the real function
	const Declaration* itsDeclaration;

 	/// functor (clone of the functor of the table)
 	/// that can get some attributes!
 	ARM_GramFctorPtr itsFctor;

	/// function to validate arguments
	static void ValidateArgStr( const string& t, int argNb, const string& typeName, const string& FuncName );

	/// function to create the function table
	static MapStringToVoid* CreateTheFuncTable();
	
	/// table containing all the functions!
	static MapStringToVoid* TheFuncTable;

public:
	/// constructor from a name
	ARM_GramFunction( const string &FuncName, const string& payCurveName );
	ARM_GramFunction( const ARM_GramFunction& rhs );
	virtual ~ARM_GramFunction();	
	ARM_GramFunction& operator=( const ARM_GramFunction &rhs );

	////for easy debugging
	static bool IsAFunctionName( const string& funcName );

	////accessors
	inline string FuncName() const { return itsDeclaration->FuncName(); }
	inline string Description() const { return itsDeclaration->Description(); }
	inline const ArgVector& Args() const { return itsDeclaration->Args(); }
	inline GramFuncArgType ReturnType() const { return itsDeclaration->ReturnType(); }
	inline GramNodeBuilderFunc GetGramNodeBuilderFunc() const { return itsDeclaration->GetGramNodeBuilderFunc(); }
	inline ARM_GramFctorPtr& Fctor() { return itsFctor; }
    inline const ARM_GramFctorPtr& Fctor() const { return itsFctor; }
	static size_t NumFunctions() { return TheFuncTable->size(); }

	/// iterator support
	typedef MapStringToVoid::const_iterator const_iterator;
	typedef MapStringToVoid::iterator iterator;
	inline const_iterator begin() const { return TheFuncTable->begin(); }
	inline const_iterator end() const { return TheFuncTable->end(); }
	inline iterator begin() { return TheFuncTable->begin(); }
	inline iterator end() { return TheFuncTable->end(); }

	/// for caching operation (performance improvement)
	void Reset( bool totalReset );

	/// standard ARM Object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="", const string& nextIndent="") const { return itsDeclaration->toString(); }

	/// init method
	static void InitTheFuncTable();
	static void ReleaseTheFuncTable();
};



CC_END_NAMESPACE()

#endif 
