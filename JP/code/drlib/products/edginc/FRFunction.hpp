//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRFunction.hpp
//
//   Description : Defines interface to functions that the parser can invoke
//
//   Author      : Mark A Robson
//
//   Date        : 26 July 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_FRFUNCTION_HPP
#define EDR_FRFUNCTION_HPP

#include "edginc/FRIfaces.hpp"

DRLIB_BEGIN_NAMESPACE

/** Defines interface to generic functions that the parser can invoke.
    ie to add a function to the FRParser write it to conform to one of 
    DFunc, DtFunc signatures and register using registerGenericXXXXX */
class PRODUCTS_DLL FRFunction{
public:
    /** Defines types of parameters that functions can choose to
        accept.  In particular, a function is passed an array of this
        union.  Note that a function can choose to either to do
        'getValue()' on its parameters or let the infrastructure do
        it. The choice is made via an array of ValCalcType passed in at
        registration time. A null array results in the infrastructure
        calling getValue() for each parameter. Eg for a parameter of
        type doubleType then if 'native' is selected for that parameter,
        the double will appear in field d otherwise the
        expression will be in dExp (if 'expression' is chosen) or the
        rValues field will be populated (if 'arrayExpression' is chosen).
        Finally, note that the array of Val passed to the function has
        one additional parameter at the end of the array - this is the
        current index of the simulation */
    typedef union{
        // note no type - the function specifies when it registers itself
        bool                          b;
        FRIfaces::IRValueBool::RT*    bExp; // bool not yet evaluated
        int                           i;
        FRIfaces::IRValueInt::RT*     iExp; // int not yet evaluated
        double                        d;
        FRIfaces::IRValueDouble::RT*  dExp; // double not yet evaluated
        FRIfaces::IRValueDoubleArray::RT* dExpArray; /* db array not
                                                        yet evaluated */
        FRIfaces::IRValueIntArray::RT*    iExpArray; /* int array not
                                                        yet evaluated */
        FRIfaces::IRValueBoolArray::RT*   bExpArray; /* bool array not
                                                        yet evaluated */
        const DateTime::Date*         dt;
        FRIfaces::IRValueDate::RT*    dtExp; // date not yet evaluated
        const Schedule*               sched;
        FRIfaces::IRValueSchedule*    schedExp; // sched not yet evaluated
        const TabulatedFunc*          tabFunc;
        FRIfaces::IRValueTabulatedFunc* tabFuncExp; //tabFunc not yet evaluated
        const vector<FRIfaces::RValUnion>* rValues; // values for all indexes
    } Val;
    
    typedef enum _ValCalcType{
        native = 0, // ie one of b, i, d etc
        expression, // ie one of bExp, iExp, dExp
        arrayExpression // ie rValues
    } ValCalcType;

    /** the signature of functions that returns a bool */
    typedef bool (BFunc)(Val* values);
    /** the signature of functions that returns an int */
    typedef int (IFunc)(Val* values);
    /** the signature of functions that returns a double */
    typedef double (DFunc)(Val* values);
    /** the signature of functions that returns a bool array. Function
        must ensures dbles is set to correct size. */
    typedef void (BArrayFunc)(Val* values, BoolArray& bools);
    /** the signature of functions that returns an int array. Function
        must ensures ints is set to correct size. */
    typedef void (IArrayFunc)(Val* values, IntArray& ints);
    /** the signature of functions that returns a double array. Function
        must ensures dbles is set to correct size. */
    typedef void (DArrayFunc)(Val* values, DoubleArray& dbles);
    /** the signature of functions that returns a date */
    typedef const DateTime::Date& (DtFunc)(Val* values);

    //// the type of function supported
    typedef union{
        BFunc*       bFunc;
        IFunc*       iFunc;
        DFunc*       dFunc;
        BArrayFunc*  bArrayFunc;
        IArrayFunc*  iArrayFunc;
        DArrayFunc*  dArrayFunc;
        DtFunc*      dtFunc;
    } Func;
        
    /** Defines how the parser passes individual parameters to
     * this class */
    typedef union{
        FRIfaces::IRValueDouble*      dExp;//for r-values evaluating to doubles
        FRIfaces::IRValueBool*        bExp; // for r-values evaluating to bools
        FRIfaces::IRValueInt*         iExp;  // for r-values evaluating to ints
        FRIfaces::IRValueDoubleArray* dExpArray; // DoubleArray
        FRIfaces::IRValueIntArray*    iExpArray; // IntArray
        FRIfaces::IRValueBoolArray*   bExpArray; // BoolArray
        FRIfaces::IRValueDate*        dtExp;// for r-values evaluating to dates
        FRIfaces::IRValueSchedule*    schedExp; // ... evaluating to schedules
        FRIfaces::IRValueTabulatedFunc* tabFuncExp; // etc
        const vector<FRIfaces::RValUnion>* rValues; // values for all indexes
    } Exp;
    
    /** This is a general version of registerGenericDoubleFunction,
        registerGenericDateFunction etc */
    static void registerGenericFunction(
        const char*               funcName,
        const char*               description,
        int                       numParams,
        const FRIfaces::VarType*  types,  // must not go out of scope
        const char**              params, // must not go out of scope
        const ValCalcType*        calcType, // null => native for all params
        FRIfaces::VarType         returnType,
        const Func&               func);

    /** Register a function that returns a bool */
    static void registerGenericBoolFunction(
        const char*              funcName,
        const char*              description,
        int                      numParams,
        const FRIfaces::VarType* types,
        const char**             params,
        const ValCalcType*       calcType, // null => native for all params
        BFunc*                   func);

    /** Register a function that returns an int */
    static void registerGenericIntFunction(
        const char*              funcName,
        const char*              description,
        int                      numParams,
        const FRIfaces::VarType* types,
        const char**             params,
        const ValCalcType*       calcType, // null => native for all params
        IFunc*                   func);

    /** Register a function that returns a double */
    static void registerGenericDoubleFunction(
        const char*              funcName,
        const char*              description,
        int                      numParams,
        const FRIfaces::VarType* types,
        const char**             params,
        const ValCalcType*       calcType, // null => native for all params
        DFunc*                   func);

    /** Register a function that returns a date */
    static void registerGenericDateFunction(
        const char*              funcName,
        const char*              description,
        int                      numParams,
        const FRIfaces::VarType* types,
        const char**             params,
        const ValCalcType*       calcType, // null => native for all params
        DtFunc*                  func);

    /** Register a function which is 'built in' to the grammar rules. The
        numParams is just for documentation purposes */
    static void registerParserSpecificFunction(int         parserType,
                                               const char* name,
                                               const char* description,
                                               int         numParams);
    
    /* returns true if the function is a 'generic' function - in which case
       getGrammarType returns the FRIfaces::VarType
       of the return type. False if the function is integrated into the
       parser */
    bool isGeneric() const;

    /** Returns the type of the function (either as far as the parser
        is concerned (if isGeneric() is false) or returns the
        FRIfaces::VarType of the return type.*/
    int getGrammarType() const;
 
    /** turns enum into string */
    static const char* argTypeToString(
        FRIfaces::VarType varType);

    /** turns enum into string with optional '[]' */
    static const char* argTypeToString(
        FRIfaces::VarType varType, bool isArray);

    //// returns the function's name eg "MAX"
    const char* name() const;

    //// returns the description of the function
    const char* getDescription() const;

    //// how many parameters the function takes
    int numArgs() const;

    //// returns type of parameter given by index
    FRIfaces::VarType argType(int index) const;

    //// returns array of types identifying the parameters
    const FRIfaces::VarType* argTypes() const;

    /** returns either null or an array  indicating whether each
        of the parameters should be native (precalculated) or left as
        expressions (either for current index or for all indexes) */
    const ValCalcType* getValCalcType() const;

    //// returns the function that actually does the work
    const Func& getFunc() const;

    /** Returns FRFunction in static] table given in its name. Returns null if
        it does not exist */
    static const FRFunction* find(const string& funcName);

    /** Appends to the supplied arrays, the details about this function */
    void getFuncInfo(StringArray&      names,
                     StringArray&      sigs,
                     StringArray&      descs) const;

    FRFunction(); // needed for silly template
    ~FRFunction();

private:
    // fields
    int                      grammarType; // the magic value for the grammer
    bool                     generic;
    const char*              funcName;
    const char*              description;
    int                      numParams;
    const FRIfaces::VarType* types; // only if generic
    const char**             params;// only if generic
    const ValCalcType*       calcType; /* only if generic. null => native
                                          for all */
    Func                     func;
    
    FRFunction(int             grammarType,
               const char*     funcName,
               const char*     description,
               int             numParams);
    
    FRFunction(const char*              funcName,
               const char*              description,
               int                      numParams,
               const FRIfaces::VarType* types,
               const char**             params,
               const ValCalcType*       calcType,
               FRIfaces::VarType        returnType,
               const Func&              func);

    class AddinHelp;
};

DRLIB_END_NAMESPACE
#endif
