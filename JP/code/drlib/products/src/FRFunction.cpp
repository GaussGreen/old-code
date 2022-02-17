//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRFunction.cpp
//
//   Description : Defines interface to functions that the parser can invoke
//
//   Author      : Mark A Robson
//
//   Date        : 26 July 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FRFunction.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

/** Returns the type of the function (either as far as the parser
    is concerned (if isGeneric() is false) or returns the
    FRIfaces::VarType of the return type.*/
int FRFunction::getGrammarType() const{
    return grammarType;
}

/* returns true if the function is a 'generic' function - in which case
   getGrammarType returns the FRIfaces::VarType
   of the return type. False if the function is integrated into the
   parser */
bool FRFunction::isGeneric() const{
    return generic;
}

//// returns the function's name eg "MAX"
const char* FRFunction::name() const{
    return funcName;
}
//// returns the description of the function
const char* FRFunction::getDescription() const{
    return description;
}

//// how many parameters the function takes
int FRFunction::numArgs() const{
    return numParams;
}
//// returns type of parameter given by index
FRIfaces::VarType FRFunction::argType(int index) const{
    return types[index];
}

//// returns array of types identifying the parameters
const FRIfaces::VarType* FRFunction::argTypes() const{
    return types;
}
/** returns either null or an array  indicating whether each
    of the parameters should be native (precalculated) or left as
    expressions (either for current index or for all indexes) */
const FRFunction::ValCalcType* FRFunction::getValCalcType() const{
    return calcType;
}

//// returns the function that actually does the work
const FRFunction::Func& FRFunction::getFunc() const{
    return func;
}

/** turns enum into string */
const char* FRFunction::argTypeToString(
    FRIfaces::VarType varType){
    switch (varType){
    case FRIfaces::doubleType:
        return "double";
    case FRIfaces::intType:
        return "int";
    case FRIfaces::boolType:
        return "bool";
    case FRIfaces::dateType:
        return "Date";
    case FRIfaces::scheduleType:
        return "Schedule";
    case FRIfaces::tabulatedFuncType:
        return "TabulatedFunc";
    case FRIfaces::doubleArrayType:
        return "DoubleArray";
    case FRIfaces::intArrayType:
        return "IntArray";
    case FRIfaces::boolArrayType:
        return "BoolArray";
    default:
        throw ModelException("FRFunction::argTypeToString", "Unrecognised "
                             "variable type");
    }
}

/** turns enum into string */
const char* FRFunction::argTypeToString(
    FRIfaces::VarType varType, bool isArray){
    if (!isArray){
        return argTypeToString(varType);
    }
    switch (varType){
    case FRIfaces::doubleType:
        return "double[]";
    case FRIfaces::intType:
        return "int[]";
    case FRIfaces::boolType:
        return "bool[]";
    case FRIfaces::dateType:
        return "Date[]";
    case FRIfaces::scheduleType:
        return "Schedule[]";
    case FRIfaces::tabulatedFuncType:
        return "TabulatedFunc[]";
    case FRIfaces::doubleArrayType:
        return "DoubleArray[]";
    default:
        throw ModelException("FRFunction::argTypeToString", "Unrecognised "
                             "variable type");
    }
}


//// map to hold functions against names
static map<string, FRFunction> functionMap;

/** Returns FRFunction in static] table given in its name. Returns null if
    it does not exist */
const FRFunction* FRFunction::find(const string& funcName){
    map<string, FRFunction>::iterator iter(functionMap.find(funcName));
    if (iter == functionMap.end()){
        return 0;
    }
    return &iter->second;
}

/** This is a general version of registerGenericDoubleFunction,
    registerGenericDateFunction etc */
void FRFunction::registerGenericFunction(
    const char*                  funcName,
    const char*                  description,
    int                          numParams,
    const FRIfaces::VarType*     types,
    const char**                 params,
    const ValCalcType*           calcType, // null => native for all params
    FRIfaces::VarType            returnType,
    const Func&                  func){
    for (int i = 0; i < numParams; i++){
        if (!calcType || calcType[i] == native){
            switch (types[i]){
            case FRIfaces::doubleArrayType:
            case FRIfaces::intArrayType:
            case FRIfaces::boolArrayType:
                throw ModelException("FRFunction::registerGenericFunction",
                                     "'Native' for arrays not supported");
            default:
                ; // ok
            }
        }
    }
    functionMap[funcName] = FRFunction(funcName, description, numParams,
                                       types, params, calcType, 
                                       returnType, func);
}

/** Register a function that returns a bool */
void FRFunction::registerGenericBoolFunction(
    const char*              funcName,
    const char*              description,
    int                      numParams,
    const FRIfaces::VarType* types,
    const char**             params,
    const ValCalcType*       calcType, // null => native for all params
    BFunc*                   func){
    Func funcUnion;
    funcUnion.bFunc = func;
    registerGenericFunction(funcName, description, numParams, types, params,
                            calcType, FRIfaces::boolType, funcUnion);
}

/** Register a function that returns an int */
void FRFunction::registerGenericIntFunction(
    const char*              funcName,
    const char*              description,
    int                      numParams,
    const FRIfaces::VarType* types,
    const char**             params,
    const ValCalcType*       calcType, // null => native for all params
    IFunc*                   func){
    Func funcUnion;
    funcUnion.iFunc = func;
    registerGenericFunction(funcName, description, numParams, types, params,
                            calcType, FRIfaces::intType, funcUnion);
}

/** Register a function that returns a double */
void FRFunction::registerGenericDoubleFunction(
    const char*              funcName,
    const char*              description,
    int                      numParams,
    const FRIfaces::VarType* types,
    const char**             params,
    const ValCalcType*       calcType, // null => native for all params
    DFunc*                   func){
    Func funcUnion;
    funcUnion.dFunc = func;
    registerGenericFunction(funcName, description, numParams, types, params,
                            calcType, FRIfaces::doubleType, funcUnion);
}

/** Register a function that returns a double */
void FRFunction::registerGenericDateFunction(
    const char*              funcName,
    const char*              description,
    int                      numParams,
    const FRIfaces::VarType* types,
    const char**             params,
    const ValCalcType*       calcType, // null => native for all params
    DtFunc*                  func){
    Func funcUnion;
    funcUnion.dtFunc = func;
    registerGenericFunction(funcName, description, numParams, types, params,
                            calcType, FRIfaces::doubleType, funcUnion);
}


/** Register a function which is 'built in' to the grammar rules. The
    numParams is just for documentation purposes */
void FRFunction::registerParserSpecificFunction(int         parserType,
                                                const char* name,
                                                const char* description,
                                                int         numParams){
    functionMap[name] = FRFunction(parserType, name, description, numParams);
}
    
 
/** For functions handled explicitly by the grammer */
FRFunction::FRFunction(int             grammarType,
                       const char*     funcName,
                       const char*     description,
                       int             numParams):
    grammarType(grammarType), generic(false), funcName(funcName), 
    description(description), numParams(numParams), types(0){}

/** for general functions */
FRFunction::FRFunction(const char*              funcName,
                       const char*              description,
                       int                      numParams,
                       const FRIfaces::VarType* types,
                       const char**             params,
                       const ValCalcType*       calcType,
                       FRIfaces::VarType        returnType,
                       const Func&              func):
    grammarType(returnType), generic(true), funcName(funcName), 
    description(description),
    numParams(numParams), types(types), params(params), 
    calcType(calcType), func(func){}

FRFunction::FRFunction(){} // needed for silly template
FRFunction::~FRFunction(){}

/** Appends to the supplied arrays, the details about this function */
void FRFunction::getFuncInfo(StringArray&      names,
                             StringArray&      sigs,
                             StringArray&      descs) const{
    string signature;
    if (generic){
        signature += string(argTypeToString((FRIfaces::VarType) grammarType))+
            ' ';
    }
    signature += string(funcName)+'(';
    for (int i = 0; i< numParams; i++){
        if (i > 0){
            signature += ", ";
        }
        if (generic){
            bool isArrayExpression = calcType && calcType[i]== arrayExpression;
            signature += string(argTypeToString(types[i], isArrayExpression))+
                ((params && params[i])? ' '+ string(params[i]): "");
        } else {
            signature += '.';
        }
    }
    signature += ')';
    names.push_back(funcName);
    sigs.push_back(signature);
    descs.push_back(description? description: "Test Function");
}

        

class FRFunction::AddinHelp: public CObject{
public:
    static CClassConstSP const TYPE;
private:
    // optional parameter - function name 
    string functionName;

    static IObjectSP help(AddinHelp* params){
        ObjectArraySP objArray(new ObjectArray(3));
        CStringArraySP funcs(new StringArray());
        CStringArraySP sigs(new StringArray());
        CStringArraySP descs(new StringArray());
        (*objArray)[0] = funcs;
        (*objArray)[1] = sigs;
        (*objArray)[2] = descs;
        if (!params->functionName.empty()){
            const string& name = params->functionName;
            const FRFunction* func = FRFunction::find(name);
            if (!func){
                throw ModelException("FRFunction::AddinHelp",
                                     "No function with name "+name);
            }
            func->getFuncInfo(*funcs, *sigs, *descs);
        } else {
            for (map<string, FRFunction>::iterator iter(functionMap.begin());
                 iter != functionMap.end(); ++iter){
                if (iter->second.getDescription()){
                    // skip those with null description
                    iter->second.getFuncInfo(*funcs, *sigs, *descs);
                }
            }
        }
        return objArray;
    }

    static void load(CClassSP& clazz){
        REGISTER(AddinHelp, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(functionName, "Name of FR function");
        FIELD_MAKE_OPTIONAL(functionName);
        Addin::registerClassObjectMethod("FR_FUNC_HELP",
                                         Addin::FLEX_PAYOFF,
                                         "Provides help for FR functions",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)help);
    }

    static IObject* defaultConstructor(){
        return new AddinHelp();
    }

public: // shut the compiler up
    AddinHelp(): CObject(TYPE){};
};

CClassConstSP const FRFunction::AddinHelp::TYPE = 
CClass::registerClassLoadMethod("FRFunction::AddinHelp",
                                typeid(AddinHelp), load);

DRLIB_END_NAMESPACE
