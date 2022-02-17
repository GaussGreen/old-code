//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRIfaces.hpp
//
//   Description : Defines FlexRules interfaces
//
//   Author      : Mark A Robson
//
//   Date        : 25 March 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_FRIFACES_HPP
#define EDR_FRIFACES_HPP
#include "edginc/YieldCurve.hpp"
#include "edginc/VolRequest.hpp"
#include "edginc/RefLevel.hpp"
#include "edginc/SVGenSpot.hpp"
#include "edginc/SVGenDiscFactor.hpp"

DRLIB_BEGIN_NAMESPACE
class FRController;
class Schedule;
class TabulatedFunc;

class PRODUCTS_DLL FRIfaces{
public:
    /** Used to indicate the type of variable that an
        IRValueExpression represents. Used by the parser to
        identify variable types (see method getType()) */
    typedef enum _VarType{
        doubleType = 0,
        intType,
        boolType,
        dateType,
        scheduleType,
        tabulatedFuncType,
        doubleArrayType,
        intArrayType,
        boolArrayType
    } VarType;

    /** Interface RValue (cf lvalue in C language semantics) holds the
        value of an expression at a specific time point */
    class PRODUCTS_DLL IRValue{
    public:
        /** Returns the value of the variable at this time point. More of
            indicative method as derived types should supply method allowing
            easier extraction of variable (eg as double). Non const as
            method might involve calculation */
        virtual IObjectConstSP get() = 0;

        /** Is the value of this object known before the simulation starts
            eg a constant (This is useful as it allows some optimisations
            to be made) */
        virtual bool isKnown() const = 0;

        virtual ~IRValue();
    };

    typedef attic::oldCountPtr<IRValue>::type IRValueSP; // FIXME

    /** For classes that hold state information that needs to be cleared
        between each simulated path */
    class PRODUCTS_DLL IHoldsState{
    public:
        /** Previously this function was 'void reset()' - it reset any
            variable (so as to prevent reading variables that had not
            been set). However, to avoid churing through vast amounts of
            memory to reset all the variables we've turned it round so now
            the controller specifies where to store the state of the variable.
            Also using char rather than bool since sizeof(bool) = 4 on solaris
        */
        virtual void setReset(char* reset) = 0;

        virtual ~IHoldsState();
    };

     /** Interface LValue (cf lvalue in C language semantics) holds an
        expression representing the value of a variable at a specific
        time point */
    class PRODUCTS_DLL ILValue: virtual public IRValue{
    public:
        /** sets the value of the variable at this time point. More of
            indicative method as derived types should supply method allowing
            easier extraction of variable (eg as double)  */
        virtual void set(const IObjectConstSP& object) = 0;
        /** Previously this function was 'void reset()' - it reset any
            variable (so as to prevent reading variables that had not
            been set). However, to avoid churing through vast amounts of
            memory to reset all the variables we've turned it round so now
            the controller specifies where to store the state of the variable.
            Also using char rather than bool since sizeof(bool) = 4 on solaris
        */
        virtual void setReset(char* reset) = 0;

        virtual ~ILValue();
    };

    /** Interface implemented by LValues whose values are constant ie
        can be determined before the simulation starts */
    class PRODUCTS_DLL ILConstValue: virtual public ILValue{
    public:
        /** Alters the behaviour of isKnown() to return specified values.
            Normally ILConstValue objects will return true for isKnown().
            However, it is useful (for the FR tester) to be able to
            'hide' this functionality */
        virtual void setIsValueKnown(bool isKnown) = 0;

        virtual ~ILConstValue();
    };

    /** more specific interface of RValue for doubles */
    class PRODUCTS_DLL IRValueDouble: virtual public IRValue{
    public:
        typedef double (TGetValue)(void* structure);
        /** Defines smaller 'run time' structure to use. The main purpose
            is to reduce the amount of memory that needs to be processed in
            the simulation */
        struct PRODUCTS_DLL RT{
            TGetValue* func;
        };
        // a->getValue() is equivalent to (a->getRT())->func(a)

        //// get the variable expressed as a double
        virtual double getValue() = 0;
        //// get the run-time object to use ie cut down version of whole class
        virtual RT* getRT() = 0;

        virtual ~IRValueDouble();
    };

    /** more specific interface of RValue for ints */
    class PRODUCTS_DLL IRValueInt: virtual public IRValue{
    public:
        typedef int (TGetValue)(void* structure);
        /** Defines smaller 'run time' structure to use. The main purpose
            is to reduce the amount of memory that needs to be processed in
            the simulation */
        struct PRODUCTS_DLL RT{
            TGetValue* func;
        };
        // a->getValue() is equivalent to (a->getRT())->func(a)

        // get the variable expressed as an int
        virtual int getValue() = 0;

        //// get the run-time object to use ie cut down version of whole class
        virtual RT* getRT() = 0;

        virtual ~IRValueInt();
    };

    /** more specific interface of RValue for bools */
    class PRODUCTS_DLL IRValueBool: virtual public IRValue{
    public:
        typedef bool (TGetValue)(void* structure);
        /** Defines smaller 'run time' structure to use. The main purpose
            is to reduce the amount of memory that needs to be processed in
            the simulation */
        struct PRODUCTS_DLL RT{
            TGetValue* func;
        };
        // a->getValue() is equivalent to (a->getRT())->func(a)

        // get the variable expressed as a double
        virtual bool getValue() = 0;

        //// get the run-time object to use ie cut down version of whole class
        virtual RT* getRT() = 0;

        virtual ~IRValueBool();
    };

    /** more specific interface of RValue for dates */
    class PRODUCTS_DLL IRValueDate: virtual public IRValue{
    public:
        typedef const DateTime::Date& (TGetValue)(void* structure);
        /** Defines smaller 'run time' structure to use. The main purpose
            is to reduce the amount of memory that needs to be processed in
            the simulation */
        struct PRODUCTS_DLL RT{
            TGetValue* func;
        };
        // a->getValue() is equivalent to (a->getRT())->func(a)

        // get the variable expressed as a date
        virtual const DateTime::Date& getValue() = 0;

        //// get the run-time object to use ie cut down version of whole class
        virtual RT* getRT() = 0;

        virtual ~IRValueDate();
    };
    /** more specific interface of RValue for DoubleArrays */
    class PRODUCTS_DLL IRValueDoubleArray: virtual public IRValue{
    public:
        /** provides access to length of array */
        typedef int (TGetSize)(void* structure);
        /** provides access to element given by index parameter. It is
            the responsibility of the caller to ensure that index is within
            bounds */
        typedef double (TGetValue)(void* structure, int index);
        /** Defines smaller 'run time' structure to use. The main purpose
            is to reduce the amount of memory that needs to be processed in
            the simulation */
        struct PRODUCTS_DLL RT{
            TGetSize*  size; // size of the array
            TGetValue* func;
        };
        //// get the run-time object to use ie cut down version of whole class
        virtual RT* getRT() = 0;

        virtual ~IRValueDoubleArray();
    };
    /** more specific interface of RValue for IntArrays */
    class PRODUCTS_DLL IRValueIntArray: virtual public IRValue{
    public:
        /** provides access to length of array */
        typedef int (TGetSize)(void* structure);
        /** provides access to element given by index parameter. It is
            the responsibility of the caller to ensure that index is within
            bounds */
        typedef int (TGetValue)(void* structure, int index);
        /** Defines smaller 'run time' structure to use. The main purpose
            is to reduce the amount of memory that needs to be processed in
            the simulation */
        struct PRODUCTS_DLL RT{
            TGetSize*  size; // size of the array
            TGetValue* func;
        };
        //// get the run-time object to use ie cut down version of whole class
        virtual RT* getRT() = 0;

        virtual ~IRValueIntArray();
    };

    /** more specific interface of RValue for IntArrays */
    class PRODUCTS_DLL IRValueBoolArray: virtual public IRValue{
    public:
        /** provides access to length of array */
        typedef int (TGetSize)(void* structure);
        /** provides access to element given by index parameter. It is
            the responsibility of the caller to ensure that index is within
            bounds */
        typedef bool (TGetValue)(void* structure, int index);
        /** Defines smaller 'run time' structure to use. The main purpose
            is to reduce the amount of memory that needs to be processed in
            the simulation */
        struct PRODUCTS_DLL RT{
            TGetSize*  size; // size of the array
            TGetValue* func;
        };
        //// get the run-time object to use ie cut down version of whole class
        virtual RT* getRT() = 0;

        virtual ~IRValueBoolArray();
    };

    /** more specific interface of RValue for schedules */
    class PRODUCTS_DLL IRValueSchedule: virtual public IRValue{
    public:
        // get the variable expressed as a Schedule
        virtual const Schedule* getValue() = 0;

        virtual ~IRValueSchedule();
    };

    /** more specific interface of RValue for TabulatedFuncs */
    class PRODUCTS_DLL IRValueTabulatedFunc: virtual public IRValue{
    public:
        // get the variable expressed as a TabulatedFunc
        virtual const TabulatedFunc* getValue() = 0;

        virtual ~IRValueTabulatedFunc();
    };

    //// holder for IRValues with dynamic cast already done
    typedef union{
        IRValueBool::RT*        bExp; // bool not yet evaluated
        IRValueInt::RT*         iExp; // int not yet evaluated
        IRValueDouble::RT*      dExp; // double not yet evaluated
        IRValueDate::RT*        dtExp; // date not yet evaluated
        IRValueSchedule*        schedExp; // sched not yet evaluated
        IRValueTabulatedFunc*   tabFuncExp; //tabFunc not yet evaluated
        IRValueDoubleArray::RT* dArrayExp;  // double array
        IRValueIntArray::RT*    iArrayExp;  // double array
        IRValueBoolArray::RT*   bArrayExp;  // double array
    } RValUnion;

    /** more specific interface of LValue for doubles */
    class PRODUCTS_DLL ILValueDouble: virtual public ILValue,
                         virtual public IRValueDouble{
    public:
        typedef void (TSetValue)(void* structure, double value);
        /** 'Derived' from IRValueDouble::RT. This is what
            getRT() returns */
        struct PRODUCTS_DLL RT{
            TGetValue* func;
            TSetValue* setFunc;
        };

        // set the variable using a double
        virtual void setValue(double value) = 0;

        virtual ~ILValueDouble();
    };

    /** more specific interface of LValue for Arrays of doubles */
    class PRODUCTS_DLL ILValueDoubleArray: virtual public ILValue,
                              virtual public IRValueDoubleArray{
    public:
        /** Set the value of the 'index' element in the array using
            the 'index' element of rValue. If index < 0, then assign
            all elements of the array. It is up to the caller to ensure
            that index is valid if a specific index is specified.
            If idx < 0 then this method should validate the lengths
            as appropriate */
        typedef void (TSetValue)(
            void*                        structure,
            int                          index,
            IRValueDoubleArray::RT*      rValue);

        /** 'Derived' from IRValueDoubleArray::RT. This is what
            getRT() returns */
        struct PRODUCTS_DLL RT{
            TGetSize*  size; // size of the array
            TGetValue* func;
            TSetValue* setFunc;
        };

        virtual ~ILValueDoubleArray();
    };

    /** more specific interface of LValue for Arrays of ints */
    class PRODUCTS_DLL ILValueIntArray: virtual public ILValue,
                              virtual public IRValueIntArray{
    public:
        /** Set the value of the 'index' element in the array using
            the 'index' element of rValue. If index < 0, then assign
            all elements of the array. It is up to the caller to ensure
            that index is valid if a specific index is specified.
            If idx < 0 then this method should validate the lengths
            as appropriate */
        typedef void (TSetValue)(
            void*                        structure,
            int                          index,
            IRValueIntArray::RT*         rValue);

        /** 'Derived' from IRValueIntArray::RT. This is what
            getRT() returns */
        struct PRODUCTS_DLL RT{
            TGetSize*  size; // size of the array
            TGetValue* func;
            TSetValue* setFunc;
        };

        virtual ~ILValueIntArray();
    };

    /** more specific interface of LValue for Arrays of bools */
    class PRODUCTS_DLL ILValueBoolArray: virtual public ILValue,
                           virtual public IRValueBoolArray{
    public:
        /** Set the value of the 'index' element in the array using
            the 'index' element of rValue. If index < 0, then assign
            all elements of the array. It is up to the caller to ensure
            that index is valid if a specific index is specified.
            If idx < 0 then this method should validate the lengths
            as appropriate */
        typedef void (TSetValue)(
            void*                        structure,
            int                          index,
            IRValueBoolArray::RT*        rValue);

        /** 'Derived' from IRValueBoolArray::RT. This is what
            getRT() returns */
        struct PRODUCTS_DLL RT{
            TGetSize*  size; // size of the array
            TGetValue* func;
            TSetValue* setFunc;
        };

        virtual ~ILValueBoolArray();
    };

    /** more specific interface of LValue for ints */
    class PRODUCTS_DLL ILValueInt: virtual public ILValue,
                      virtual public IRValueInt{
    public:
        typedef void (TSetValue)(void* structure, int value);
        /** 'Derived' from IRValueDouble::RT. This is what
            getRT() returns */
        struct PRODUCTS_DLL RT{
            TGetValue* func;
            TSetValue* setFunc;
        };
        // set the variable using an int
        virtual void setValue(int value) = 0;

        virtual ~ILValueInt();
    };

    /** more specific interface of LValue for bools */
    class PRODUCTS_DLL ILValueBool: virtual public ILValue,
                       virtual public IRValueBool{
    public:
        typedef void (TSetValue)(void* structure, bool value);
        /** 'Derived' from IRValueDouble::RT. This is what
            getRT() returns */
        struct PRODUCTS_DLL RT{
            TGetValue* func;
            TSetValue* setFunc;
        };
        // set the variable using a bool
        virtual void setValue(bool value) = 0;

        virtual ~ILValueBool();
    };

     /** more specific interface of LValue for dates */
    class PRODUCTS_DLL ILValueDate: virtual public ILValue,
                       virtual public IRValueDate{
    public:
        // set the variable using a DateTime
        virtual void setValue(const DateTime::Date& value) = 0;

        virtual ~ILValueDate();
    };

    /** Interface RValue Expression (cf lvalue in C language
        semantics) holds the rule determing the value of an expression. */
    class PRODUCTS_DLL IRValueExpression: virtual public IObject{
    public:
        static CClassConstSP const TYPE;

        virtual ~IRValueExpression();

        /** gets IRValue from FlexRule which represents the value
            of this expression at the current timepoint. DO NOT FREE. */
        virtual IRValue* getRValue(
            FRController* frCtrl) const = 0;

        /** gets IRValue from FlexRule which represents the value
            of this expression at the specified timepoint. DO NOT FREE. */
        virtual IRValue* getRValue(
            int           index,
            FRController* frCtrl) const = 0;

        /** returns a string that can be used to identify this object */
        virtual const string& getID() const = 0;

        /** Used in conjunction with equals() by FlexRule to keep a hash
            of RValueExpressions. Implementations should ensure hashCode()
            and equals() are implemented so that objects which refer to
            the same IRValue are 'equal' */
        virtual int hashCode() const = 0;

        /** Returns true if the this and the given IRValueExpression refer
            to the IRValue */
        virtual bool equals(const IRValueExpression* name) const = 0;
    };
    typedef smartPtr<IRValueExpression>      IRValueExpressionSP;
    typedef smartConstPtr<IRValueExpression> IRValueExpressionConstSP;
    typedef array<IRValueExpressionSP,
        IRValueExpression> IRValueExpressionArray;
    typedef smartPtr<IRValueExpressionArray> IRValueExpressionArraySP;
    typedef smartConstPtr<IRValueExpressionArray>
    IRValueExpressionArrayConstSP;

    /** Interface LValue (cf lvalue in C language semantics) holds an
        expression representing the value of a variable */
    class PRODUCTS_DLL ILValueExpression:  virtual public IRValueExpression{
    public:
        static CClassConstSP const TYPE;

        virtual ~ILValueExpression();

        /** Returns the type that this variable represents */
        virtual VarType getType() const = 0;

        /** is the variable simulation date independent ie same value at
            each simulation date */
        virtual bool isSDI() const = 0;

        /** creates FlexRule::ILValue from FlexRule which represents
            the value of this variable at the current timepoint. DO
            NOT FREE. */
        virtual ILValue* getLValue(
            FRController* frCtrl) const = 0;

        /** creates FlexRule::ILValue from FlexRule which represents
            the value of this variable at a specific timepoint. DO NOT
            FREE. */
        virtual ILValue* getLValue(
            int           index,
            FRController* frCtrl) const = 0;

        /** returns an array of 3 strings containing the inputs for IMS:
            type, name, initial value */
        virtual StringArray getIMSInput() const = 0;

    };

    typedef smartPtr<ILValueExpression> ILValueExpressionSP;
    typedef array<ILValueExpressionSP,
        ILValueExpression> ILValueExpressionArray;
    typedef smartPtr<ILValueExpressionArray> ILValueExpressionArraySP;

    /** Defines object which has ability to set an lvalue given an rvalue. */
    class PRODUCTS_DLL IAssignment{
    public:
        typedef void (TAssign)(void* structure);
        struct PRODUCTS_DLL RT{
            TAssign* func;
        };
        // a->assign() is equivalent to (a->getRT())->func(a)
        //// assigns rvalue into lvalue
        virtual void assign() = 0;

        //// get the run-time object to use ie cut down version of whole class
        virtual RT* getRT() = 0;
        //// does setReset() on lvalue
        virtual void setReset(char* reset) = 0;

        virtual ~IAssignment();
    };
    typedef refCountPtr<IAssignment> IAssignmentSP;
    typedef vector<IAssignmentSP> IAssignmentArray;
    typedef refCountPtr<IAssignmentArray> IAssignmentArraySP;

    class PRODUCTS_DLL IVarBarrierLevelAssist {
    public:
        virtual string getAssetName() const = 0;
        virtual double getRefLevel(FRController* frCtrl) const = 0;

        virtual ~IVarBarrierLevelAssist() {};
    };

    class PRODUCTS_DLL IVarArrayBarrierLevelAssist {
    public:
        virtual vector<IVarBarrierLevelAssist*> 
            getComponentAssists(FRController* frCtrl) const = 0;

        virtual ~IVarArrayBarrierLevelAssist() {};
    };

    // Describes what debug info is desired
    class PRODUCTS_DLL DebugRequest : public CObject {
    public:
        static CClassConstSP const TYPE;
        friend class DebugRequestHelper;

        void validatePop2Object();

        virtual ILValueExpressionSP getFlagVar() const;
        virtual bool doMinValues() const;
        virtual bool doMaxValues() const;
        virtual bool doAvgValues() const;
        virtual bool doDistns() const;
         // Array of double type vars - each can have distn reported on specific date
        int getDistnSize() const;
        void setDistnSize(int size);
        virtual ILValueExpressionArraySP getDistnVars() const;
        virtual const DateTimeArray& getDistnDates() const; // one per var

    private:
        DebugRequest(): CObject(TYPE), 
                        captureMinValues(false), 
                        captureMaxValues(false),
                        captureAvgValues(false), 
                        captureDistn(false),
                        distnSize(0),
                        captureAllValues(false) {}; // for reflection
        DebugRequest(const DebugRequest& rhs); // not implemented
        DebugRequest& operator=(const DebugRequest& rhs); // not implemented

        // fields
        ILValueExpressionSP      captureDebug;
        bool                     captureMinValues;
        bool                     captureMaxValues;
        bool                     captureAvgValues;
        bool                     captureDistn;
        int                      distnSize; // size of output
        ILValueExpressionArraySP distnVars; // can be empty
        DateTimeArray            distnDates; // parallels distnVars
        // no longer used, but remains to avoid backwards compatability issues
        bool                     captureAllValues; 
    };
    typedef smartPtr<DebugRequest> DebugRequestSP;

    //// provides view onto product
    class PRODUCTS_DLL IProductView{
    public:
        /** Returns the value date ie today */
        virtual const DateTime& getValueDate() const = 0;

        //// create assignment array for current time point. The array is
        //// allowed to contain nulls
        virtual IAssignmentArray* createAssignments(
            FRController* frCtrl) const = 0;

        /** Returns variable which is used for calculating fair value.
            Null is allowed which is useful for testers */
        virtual const ILValueExpression* getPayVariable(
            FRController* frCtrl) const = 0;

        /** Controls debug gathering
            Null is allowed which is useful for testers */
        virtual const DebugRequest* getDebugRequest(
            const FRController* frCtrl) const = 0;

        /** returns the simulation dates */
        virtual DateTimeArrayConstSP getSimDates() const = 0;

        /** returns the notional - used to scale values of each path by */
        virtual double getNotional() const = 0;

         /** returns the payment date for each simulation date */
        virtual DateTimeArrayConstSP getPayDates() const = 0;

        /** populates array of discount factors. Each value being the
            discount factor for a payment at the corresponding date. The return
            value indicates the value of the index which is the first to be
            in the future */
        virtual int getDiscountFactors(
            const DateTimeArrayConstSP& payDates,
            DoubleArray&                discountFactors) const = 0;

        /** Given asset's name and (optional) ccy treatment, returns the
            index for that asset */
        virtual int findAssetIndex(const string& assetName,
                                   const string& ccyTreatment) const = 0;

        /** Given asset index and (optional) vol request, returns the path
            index for this asset/vol request combination. To do: review this.
            Perhaps store vol requests in controller? */
        virtual int findPathIndex(int                   assetIndex,
                                  const CVolRequestSP&  volRequest) const = 0;

        //// Given ccy name find the actual YieldCurve from instrument
        virtual YieldCurveConstSP findYieldCurve(
            const string* ccyName = 0) const = 0; // now optional - default to inst ccy

        /* returns variables used in payoff */
        virtual const ILValueExpressionArray* getVariables() const = 0;

        /** Return a string representing the rule for given simulation date
            and assignment eg a = MAX(b,c)+5. Used for error message */
        virtual string getExpression(
            int             simDateIndex,
            const DateTime& date, // for simDateIndex
            int             assignmentIndex) const = 0;

        /** Returns a reference level sv generator associated with this
            instrument. Returns null if the instrument does not have one
            (typically should only
            be requested in conjunction with spot assets) */
        virtual IRefLevel::IStateVarGenSP getRefLevelGen() const = 0;

        /** Returns a SV generator for instrument's assets' spots. Returns
            null if there are none */
        virtual SVGenSpotSP getSpotGen() const = 0;

        /** Returns a SV generator for computing discount factors on the
            settlement adjusted dates. ie this method must adjust the
            unadjustedDates via the appropriate settlement rules */
        virtual SVGenDiscFactorSP getDiscFactorGen(
            const DateTimeArray& unadjustedDates) const = 0;

        virtual ~IProductView();

    };

    /** Holds the set of rules expressing how each variable is to be
        calculated for a single time point */
    class PRODUCTS_DLL ITimePointRules: virtual public IObject{
    public:
        static CClassConstSP const TYPE;

        /** Creates an array of assignments - one for each rule at this
            timepoint. The array is allowed to contain nulls */
        virtual IAssignmentArray* createAssignments(
            FRController* frCtrl) const = 0;

        /** Returns the number of rules for this given timepoint ie the
            number of assignments */
        virtual int numRules() const = 0;

        /** Returns the lValue for the specified rule as a string ie the
            variable name */
        virtual string lValueID(int ruleNum) const = 0;

        /** Returns the rValue for the specified rule as a string eg X+Y */
        virtual string rValueID(int ruleNum) const = 0;

        virtual ~ITimePointRules();
    };
    typedef smartPtr<ITimePointRules> ITimePointRulesSP;
    typedef array<ITimePointRulesSP, ITimePointRules> ITimePointRulesArray;
    typedef smartPtr<ITimePointRulesArray> ITimePointRulesArraySP;

    /** Holds the set of variables for which a 'no op' should be
        performed - useful when setting up 'sparse' algorithms */
    class PRODUCTS_DLL ITimePointNoOp: virtual public IObject{
    public:
        static CClassConstSP const TYPE;

        /** Returns the number of rules for this given timepoint ie the
            number of assignments */
        virtual int numVars() const = 0;

        /** Returns the lValue for the specified variable ie the
            variables name */
        virtual const string& lValueID(int ruleNum) const = 0;

        virtual ~ITimePointNoOp();
    };
    typedef smartPtr<ITimePointNoOp> ITimePointNoOpSP;
    typedef array<ITimePointNoOpSP, ITimePointNoOp> ITimePointNoOpArray;
    typedef smartPtr<ITimePointNoOpArray> ITimePointNoOpArraySP;

    /** Holds the set of rules expressing how each variable is to be
        calculated for each time point in the simulation */
    class PRODUCTS_DLL IAlgorithm: virtual public IObject{
    public:
        static CClassConstSP const TYPE;

        /** create assignment array for current time point. The array
            is allowed to contain nulls */
        virtual IAssignmentArray* createAssignments(
            FRController* frCtrl) const = 0;

        /** write the set of rules in a human readable form to the supplied
            file. The payoutVariable indicates which is the 'driving'
            variable */
        virtual void writeRules(const string& fileName,
                                const string& payoutVariable) const = 0;

        /** Return a string representing the rule for given simulation date
            and assignment eg a = MAX(b,c)+5. Used for error message */
        virtual string getExpression(
            int             simDateIndex,
            const DateTime& date, // for simDateIndex
            int             assignmentIndex) const = 0;

        /** Report all dates known to the algorithm */
        virtual DateTimeArraySP getAllDates() const = 0;

        /** returns an array of array of strings representing the input for IMS
         the first array contains the date rule IDs
         the second array contains the ID of each rule
         the third array contains the rules */
        virtual StringArrayArray getIMSInput() const = 0;
        virtual void getIMSInput(
            IntArraySP    simDateRuleIds,
            IntArraySP    ruleId,
            StringArraySP ruleDefn) const = 0;

        virtual ~IAlgorithm();
    };
    typedef smartPtr<IAlgorithm> IAlgorithmSP;

};

DRLIB_END_NAMESPACE
#endif
