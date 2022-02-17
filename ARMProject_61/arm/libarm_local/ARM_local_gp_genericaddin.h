
/*! \file ARM_local_gp_genericaddin.h
 *
 *  \brief file for the generic addin
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date October 2006
 */


#ifndef ARMLOCAL_GP_GENERICADDIN_H
#define ARMLOCAL_GP_GENERICADDIN_H

#include "firstToBeIncluded.h"
#include <gpbase/rootobject.h>
#include <vector>
#include <map>
#include <functional>
using ARM::ARM_RootObject;
using std::vector;
using std::map;
using std::binary_function;

#include <ARM\libarm\ARM_result.h>

///////////////////////////////////////////////////////////////////////////////////
// In this file we define a new infrastructure to build the GENERIC ADDIN:
// With this new fuctionality we can create an EXCEL interface for any ARM
// function with :
// _ a description of the function
// _ a functor implementing its behaviour
//
// The user can use 2 generic addin
// _ ARM_GP_GenericAddin( FunctionName, ParamsNames, ParamsValues )
//
// With the GENERIC ADDIN you can call any function of the library
//
// _ ARM_GP_GenericAddinHelper( FunctionName )
//
// With the GENERIC ADDIN HELPER you can get an help for any function 
// of the library
///////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////
//
// This structure represents the value of any parameter: double, string or object
//
///////////////////////////////////////////////////////////////////////////////////

typedef enum {GA_DOUBLE, GA_STRING, GA_OBJECT, GA_DOUBLE_VECTOR, GA_STRING_VECTOR, GA_TERMINAL} GA_Type;

struct ARM_GenericParamValue
{
	explicit ARM_GenericParamValue(double value);
	explicit ARM_GenericParamValue(const string& value);
	explicit ARM_GenericParamValue(long value);

	void Init();

	void SetType(GA_Type newType) { itsType = newType; }
	GA_Type GetType() const { return itsType; }
	double GetDouble() const;
	const string& GetString() const;
	long GetObjectId() const;

	GA_Type itsType;
	double itsDoubleValue;
	string itsStringValue;
	long itsObjectValue;
};

///////////////////////////////////////////////////////////////////////////////////
//
// Those structures will be used for convinient constructors
//
///////////////////////////////////////////////////////////////////////////////////

struct GenericParamStruct
{
	GA_Type Type;
	const char* ParamName;
	const char* ParamDescription;
	ARM_GenericParamValue* DefaultValue;
};

class ARM_GenericAddinFunctor;

struct GenericAddinStruct
{
	const char* FunctionName;
	const char* FunctionDescription;
	const GenericParamStruct* params;
	ARM_GenericAddinFunctor* Functor;
	bool RetObj;
};

///////////////////////////////////////////////////////////////////////////////////
//
// This class decribes a parameter: name, description, default value
//
///////////////////////////////////////////////////////////////////////////////////

class ARM_GenericParamDesc : public ARM_RootObject
{
public:
	ARM_GenericParamDesc(
		GenericParamStruct paramStruct) 
		: itsType(paramStruct.Type),
		itsParamName(paramStruct.ParamName),
		itsParamDescription(paramStruct.ParamDescription),
		itsDefaultParamValue(paramStruct.DefaultValue) {};

	ARM_GenericParamDesc(
		const ARM_GenericParamDesc& rhs)
		: itsType(rhs.itsType),
		itsParamName(rhs.itsParamName),
		itsParamDescription(rhs.itsParamDescription),
		itsDefaultParamValue(rhs.itsDefaultParamValue) {};

	virtual ~ARM_GenericParamDesc() {}

	GA_Type GetType() const { return itsType; }
	const string& GetParamName() const { return itsParamName; }
	const string& GetParamDescription() const { return itsParamDescription; }
	const ARM_GenericParamValue* GetDefaultParamValue() const { return itsDefaultParamValue; }

	virtual ARM_Object* Clone() const { return new ARM_GenericParamDesc(*this); }
    virtual string toString(const string& indent="", const string& nextIndent="") const { return ""; };

private:
	
	GA_Type itsType;
	string itsParamName;
	string itsParamDescription;
	const ARM_GenericParamValue* itsDefaultParamValue;
};

///////////////////////////////////////////////////////////////////////////////////
//
// This class decribes a set of parameter
//
///////////////////////////////////////////////////////////////////////////////////

class ARM_GenericParams;

class ARM_GenericParamsDesc : ARM_RootObject
{
public:
	ARM_GenericParamsDesc(const string& functionName, const GenericParamStruct* paramsStruct);


	ARM_GenericParamsDesc(const ARM_GenericParamsDesc& rhs)
		: itsFunctionName(rhs.itsFunctionName), 
		itsParams(rhs.itsParams),	
		itsParamNames(rhs.itsParamNames)
	{
	};

	~ARM_GenericParamsDesc() {}

	virtual ARM_Object* Clone() const { return new ARM_GenericParamsDesc(*this); }
    virtual string toString(const string& indent="", const string& nextIndent="") const;

	const ARM_GenericParamDesc& GetParam(const string& name) const;
	vector<string> GetParamNames(bool withDefault) const;

	void fillDefaultParams(ARM_GenericParams* genericParams);

private:
	typedef map<string,ARM_GenericParamDesc> GenericParamMap;
	string itsFunctionName;
	vector<string> itsParamNames;
	GenericParamMap itsParams;
};


///////////////////////////////////////////////////////////////////////////////////
//
// This class defines the map of the parameter : name / value
// It will be build with the parameter of the GENERIC ADDIN
//
///////////////////////////////////////////////////////////////////////////////////

class ARM_GenericParams : ARM_RootObject
{
public:
	ARM_GenericParams(const string& functionName)
		: 
	itsFunctionName(functionName),
	itsMap()
	{
	}

	ARM_GenericParams(const ARM_GenericParams& rhs)
		: 
	itsFunctionName(rhs.itsFunctionName),
	itsMap(rhs.itsMap)
	{
	}

	~ARM_GenericParams() {}

	void SetParamValue(const string& name, const ARM_GenericParamValue& value);
	const ARM_GenericParamValue& GetParamValue(const string& name) const;
	bool IsExists(const string& name) const;

	virtual ARM_Object* Clone() const { return new ARM_GenericParams(*this); }
    virtual string toString(const string& indent="", const string& nextIndent="") const { return "";}

private:
	typedef map<string,ARM_GenericParamValue> GenericParamValueMap;
	string itsFunctionName;
	GenericParamValueMap itsMap;
};

///////////////////////////////////////////////////////////////////////////////////
//
// This class defines the functor associated to the function
//
///////////////////////////////////////////////////////////////////////////////////


/*!
 * General function for all excel interface 
 * intended to avoid the monstruous copy and
 * paste in excel interface
 *
 * the class T has to be a functor that 
 * supports long operator()( ARM_result, long )
 * derived from binary_function to make it easily usable
 * for STL
 *
 */
class ARMResultLong2LongFunc : public binary_function< ARM_result, long, long >
{
public:
	/// pure virtual function to
	/// force redefinition
	virtual	long operator()( ARM_result& result, long objId ) = 0;
};

class ARM_GenericAddinFunctor : public ARMResultLong2LongFunc
{
public:
	typedef enum {DOUBLE, STRING} RetType;
	typedef struct
	{
		RetType type;
		double dblVal;
		string strVal;
	} RetStruct;

	void SetGenericParams(ARM_GenericParams* params) { itsGenericParams = params; }
	ARM_GenericParams* GetGenericParams() const { return itsGenericParams; }

	virtual char* ClassName() const { return NULL; }

	void ResizeRetValues(int nbRows, int nbCols);
	void SetValue(int row, int col, double val);
	void SetValue(int row, int col, const string& val);

	int GetNbRows() const { return itsNbRows; }
	int GetNbCols() const { return itsNbCols; }
	const vector<RetStruct>& GetRetValues() { return itsRetValues; }

	
private:
	ARM_GenericParams* itsGenericParams;

	vector<RetStruct> itsRetValues;
	int itsNbRows;
	int itsNbCols;
};

class ARM_GenericAddinDesc : public ARM_RootObject
{
public:
	ARM_GenericAddinDesc(const GenericAddinStruct& addinStruct);

	ARM_GenericAddinDesc(const ARM_GenericAddinDesc& rhs) 
		: itsFunctionName(rhs.itsFunctionName),
		itsFunctionDescription(rhs.itsFunctionDescription),
		itsParams(rhs.itsParams),
		itsFunctor(rhs.itsFunctor),
		itsRetObj(rhs.itsRetObj) {}

	~ARM_GenericAddinDesc() {};

	const string& GetFunctionName() const { return itsFunctionName; }
	const string& GetFunctionDescription() const { return itsFunctionDescription; }
	const ARM_GenericParamsDesc& GetParams() const { return itsParams; }
	ARM_GenericAddinFunctor* GetFunctor() const { return itsFunctor; }
	bool IsRetObj() const { return itsRetObj; }

	void fillDefaultParams(ARM_GenericParams* genericParams) { itsParams.fillDefaultParams(genericParams); }

	virtual ARM_Object* Clone() const { return new ARM_GenericAddinDesc(*this); }
    virtual string toString(const string& indent="", const string& nextIndent="") const;

private:
	string itsFunctionName;
	string itsFunctionDescription;
	ARM_GenericParamsDesc itsParams;
	ARM_GenericAddinFunctor* itsFunctor;
	bool itsRetObj;
};

class ARM_GenericAddinsDesc : public ARM_RootObject
{
public:
	static void CreateTheAddingTable();
	static void ReleaseTheAddingTable();

	static ARM_GenericAddinsDesc* TheAddinsTable() { return itsAddinsTable; }	

	ARM_GenericAddinsDesc() {};
	ARM_GenericAddinsDesc(const GenericAddinStruct* addinStruct, size_t nbAddins);

	ARM_GenericAddinsDesc(const ARM_GenericAddinsDesc& rhs) 
		: itsMap(rhs.itsMap){}

	~ARM_GenericAddinsDesc() {};

	ARM_GenericAddinDesc GetAddin(const string& name);

	ARM_GenericAddinsDesc ExtractAddinsDesc(const string& name);

	virtual ARM_Object* Clone() const { return new ARM_GenericAddinsDesc(*this); }
    virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual string ExportShortName() const { return "LGAHE";}

private:
	typedef map<string,ARM_GenericAddinDesc> GenericAddinMap;
	GenericAddinMap itsMap;	

	static ARM_GenericAddinsDesc* itsAddinsTable;
};

long ARMLOCAL_GenericAddin_Helper(
		const string& functionName,
		ARM_result&	result, 
        long objId=-1);

long ARMLOCAL_GenericAddin_ParamNames(
		const string& functionName,
		bool withDefault,
		vector<string>& paramNames,
		ARM_result&	result);


#endif