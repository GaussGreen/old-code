//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FR.hpp
//
//   Description : Defines common concrete FlexRules objects available at the 
//                 spreadsheet side of things
//
//   Author      : Mark A Robson
//
//   Date        : 25 March 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_FR_HPP
#define EDR_FR_HPP
#include "edginc/FRIfaces.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE
#if defined(_MSC_VER) && defined(NDEBUG)
#define EDR_INLINE __forceinline
#else
#define EDR_INLINE inline
#endif

// defining this gives better debug info
#define EDR_FR_VAR_DEBUG_EXTRA

/** Defines common FlexRules objects */
class PRODUCTS_DLL FR{
public:
    /** Common base class for IRValueExpressions. Note that this base class
        overrides the clone method to provide a shallow clone. Therefore
        all methods/fields should be const (unless clone method overridden) */
    class PRODUCTS_DLL RValueBase: public CObject,
                      virtual public FRIfaces::IRValueExpression{
    protected:
        /** creates FlexRule::IRValue which represents the value
            of this expression at the specific timepoint. The getRValue()
            uses this method and then saves it in the FRController. Return
            null if the save has already been done */
        virtual FRIfaces::IRValue* createRValue(
            int           index,
            FRController* frCtrl) const = 0;
        
        /** returns the id for this RValueBase */
        virtual const string& getID() const = 0;

        /* constructor */
        RValueBase(const CClassConstSP& clazz);

    public:
        static CClassConstSP const TYPE;

        /** returns shallow copy */
        virtual IObject* clone() const;

        /** gets FRIfaces::IRValue from FlexRule which represents the value
            of this expression at the current timepoint. DO NOT FREE. */
        virtual FRIfaces::IRValue* getRValue(FRController* frCtrl) const;

        /** gets FRIfaces::IRValue from FlexRule which represents the value
            of this expression at the specified timepoint. DO NOT FREE. */
        virtual FRIfaces::IRValue* getRValue(int           index,
                                             FRController* frCtrl) const;
        
        /** Used in conjunction with equals() by FlexRule to keep a
            hash of RValueExpressions. Here, objects with the pointer value
            are deemed to refer to the same IRValue. */
        virtual int hashCode() const;

        /** Returns true if the this and the given IRValueExpression refer
            to the same IRValue. Done via comparison of pointers */
        virtual bool equals(const FRIfaces::IRValueExpression* name) const;

    private:
        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz);
    };

    /** Common base class for ILValueExpressions. Note that the parent
        of this class overrides the clone method to provide a shallow
        clone. Therefore all methods/fields should be const (unless
        clone method overridden) */
    class PRODUCTS_DLL LValueBase: public RValueBase,
                      virtual public FRIfaces::ILValueExpression{
    protected:
        LValueBase(const CClassConstSP& clazz);

    public:
        static CClassConstSP const TYPE;

        /** is the variable simulation date independent ie same value at
            each simulation date. This implementation returns false.
            Inheriting classes may need to override this */
        virtual bool isSDI() const;
        
        /** Overrides implementation in RValueBase - hashes string
            returned by getID */
        virtual int hashCode() const;

        /** Returns true if both objects have the same id string */
        virtual bool equals(const FRIfaces::IRValueExpression* name) const;

        /** gets FlexRule::ILValue from FlexRule which represents the value
            of this expression at the current timepoint. Implemented via
            cast of return object from getRValue. DO NOT FREE. */
        virtual FRIfaces::ILValue* getLValue(FRController* frCtrl) const;

        /** creates FlexRule::ILValue  from FlexRule which represents the value
            of this variable at a specific timepoint. Implemented via
            cast of return object from getRValue. DO NOT FREE. */
        virtual FRIfaces::ILValue* getLValue(int           index,
                                             FRController* frCtrl) const;
    private:
        /** Invoked when Class is 'loaded' */
        static void load(CClassSP& clazz);
    };

    /** Class to manage memory allocation for Flex Rules. Ensures memory is
        allocated for each rule/expression as close as possible to each 
        other. This should help performance. All memory allocated from here
        is freed in one big block when shutdown is called */
    class PRODUCTS_DLL MemMgr{
    public:
        /** Sets up intitial cache. numRules and numDates need not be precise
            - just an approximation */
        static void intialise(int numRules, int numDates);

        /** frees any memory */
        static void shutdown();

        /** allocates memory of requested size */
        static void* alloc(size_t size);

        /** deallocates memory - currently does nothing */
        static void dealloc(void* ptr);
    private:
        static vector<void*> heaps;    // the memory
        static size_t        pos;      // where we are in heaps.back()
        static size_t        heapSize; // size of heaps
        static const size_t  bytesPerRule; // guess for bytes per rule
    };

    /** abstract base class for classes implementing IRValueDouble */
    class PRODUCTS_DLL RValueDouble: virtual public FRIfaces::IRValueDouble{
    public:
        /** Just wrapper around getValue method */
        virtual IObjectConstSP get();

        /** Is the value of this object known before the simulation starts
            eg a constant. Implementation here returns false */
        virtual bool isKnown() const;
        
        virtual ~RValueDouble();
    };

    /** abstract base class for classes implementing IRValueInt */
    class PRODUCTS_DLL RValueInt: virtual public FRIfaces::IRValueInt{
    public:
        /** Just wrapper around getValue method */
        virtual IObjectConstSP get();

        /** Is the value of this object known before the simulation starts
            eg a constant. Implementation here returns false */
        virtual bool isKnown() const;

        virtual ~RValueInt();
    };
    
    /** abstract base class for classes implementing IRValueBool */
    class PRODUCTS_DLL RValueBool: virtual public FRIfaces::IRValueBool{
    public:
        /** Just wrapper around getValue method */
        virtual IObjectConstSP get();

        /** Is the value of this object known before the simulation starts
            eg a constant. Implementation here returns false */
        virtual bool isKnown() const;

        virtual ~RValueBool();
    };

    /** abstract base class for classes implementing IRValueBool */
    class PRODUCTS_DLL RValueDate: virtual public FRIfaces::IRValueDate{
    public:
        /** Just wrapper around getValue method */
        virtual IObjectConstSP get();

        /** Is the value of this object known before the simulation starts
            eg a constant. Implementation here returns false */
        virtual bool isKnown() const;

        virtual ~RValueDate();
    };
    /** abstract base class for classes implementing IRValueDoubleArray */
    class PRODUCTS_DLL RValueDoubleArray: virtual public FRIfaces::IRValueDoubleArray{
    public:
        /** Just wrapper around getValue method */
        virtual IObjectConstSP get();

        /** Is the value of this object known before the simulation starts
            eg a constant. Implementation here returns false */
        virtual bool isKnown() const;
        
        virtual ~RValueDoubleArray();
    };
    /** abstract base class for classes implementing IRValueIntArray */
    class PRODUCTS_DLL RValueIntArray: virtual public FRIfaces::IRValueIntArray{
    public:
        /** Just wrapper around getValue method */
        virtual IObjectConstSP get();

        /** Is the value of this object known before the simulation starts
            eg a constant. Implementation here returns false */
        virtual bool isKnown() const;
        
        virtual ~RValueIntArray();
    };

    /** abst`ract base class for classes implementing IRValueBoolArray */
    class PRODUCTS_DLL RValueBoolArray: virtual public FRIfaces::IRValueBoolArray{
    public:
        /** Just wrapper around getValue method */
        virtual IObjectConstSP get();

        /** Is the value of this object known before the simulation starts
            eg a constant. Implementation here returns false */
        virtual bool isKnown() const;
        
        virtual ~RValueBoolArray();
    };

    /** Utility Class for a constant double value */
    class PRODUCTS_DLL RConstDouble: public RValueDouble{
    private:
        struct PRODUCTS_DLL MyRT{
            TGetValue* func;
            double     value;

            explicit MyRT(double  value);
            
            static double getValue(void* structure);
        };
        MyRT* rt;
    public:
        // simple constructor
        RConstDouble(double value);

        // get the variable expressed as a double
        virtual double getValue();

        /** Is the value of this object known before the simulation starts
            eg a constant. Here we return true */
        virtual bool isKnown() const;

        //// get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueDouble::RT* getRT();

        virtual ~RConstDouble();
    };

    /** Utility Class for a constant int value */
    class PRODUCTS_DLL RConstInt: public RValueInt{
    private:
        struct PRODUCTS_DLL MyRT{
            TGetValue*    func;
            int           value;

            explicit MyRT(int  value);
            
            static int getValue(void* structure);
        };
        MyRT* rt;
    public:
        // simple constructor
        RConstInt(int value);

        // get the variable expressed as a int
        virtual int getValue();

        /** Is the value of this object known before the simulation starts
            eg a constant. Here we return true */
        virtual bool isKnown() const;

        //// get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueInt::RT* getRT();

        virtual ~RConstInt();
    };

    /** Utility Class for a constant bool value */
    class PRODUCTS_DLL RConstBool: public RValueBool{
    private:
        struct PRODUCTS_DLL MyRT{
            TGetValue*    func;
            bool          value;
            
            explicit MyRT(bool  value);
            
            static bool getValue(void* structure);
        };
        MyRT* rt;
    public:
        // simple constructor
        RConstBool(bool value);

        // get the variable expressed as a bool
        virtual bool getValue();

        /** Is the value of this object known before the simulation starts
            eg a constant. Here we return true */
        virtual bool isKnown() const;

        //// get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueBool::RT* getRT();

        virtual ~RConstBool();
    };

    /** Utility Class for a constant double value */
    class PRODUCTS_DLL RConstDate: public RValueDate{
    private:
        DateTime::Date  value;
    public:
        // simple constructor
        RConstDate(const DateTime::Date& value);

        // get the variable expressed as a double
        virtual const DateTime::Date& getValue();

        /** Is the value of this object known before the simulation starts
            eg a constant. Here we return true */
        virtual bool isKnown() const;

        //// get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueDate::RT* getRT();

        virtual ~RConstDate();
    };
    /** class for classes implementing IRValueDoubleArray - pass in
        array of FRIfaces::IRValueDouble's (or one by one).
        Main use: for arrays constructed on the fly eg {a, b, c} */
    class PRODUCTS_DLL RDoubleArray: public RValueDoubleArray{
        struct MyRT;
        MyRT*                               rt;
        vector<FRIfaces::IRValueDouble*>    vectorRValues;
        bool                                known;
    public:
        /** constructor for initially empty array */
        RDoubleArray();

        /** Add an FRIfaces::IRValueDouble to the object.  Essentially does
            a push_back. Just takes a reference - does not copy or free. */
        void addArg(FRIfaces::IRValueDouble* rValue);

        //// get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueDoubleArray::RT* getRT();
        /** Is the value of this object known before the simulation starts
            eg a constant.  */
        virtual bool isKnown() const;

        virtual ~RDoubleArray();
    };

    /** class for classes implementing IRValueIntArray - pass in
        array of FRIfaces::IRValueInt's (or one by one).
        Main use: for arrays constructed on the fly eg {a, b, c} */
    class PRODUCTS_DLL RIntArray: public RValueIntArray{
        struct MyRT;
        MyRT*                            rt;
        vector<FRIfaces::IRValueInt*>    vectorRValues;
        bool                             known;
    public:
        /** constructor for initially empty array */
        RIntArray();

        /** Add an FRIfaces::IRValueInt to the object.  Essentially does
            a push_back. Just takes a reference - does not copy or free. */
        void addArg(FRIfaces::IRValueInt* rValue);

        //// get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueIntArray::RT* getRT();
        /** Is the value of this object known before the simulation starts
            eg a constant.  */
        virtual bool isKnown() const;

        virtual ~RIntArray();
    };

    /** class for classes implementing IRValueBoolArray - pass in
        array of FRIfaces::IRValueBool's (or one by one).
        Main use: for arrays constructed on the fly eg {a, b, c} */
    class PRODUCTS_DLL RBoolArray: public RValueBoolArray{
        struct MyRT;
        MyRT*                            rt;
        vector<FRIfaces::IRValueBool*>    vectorRValues;
        bool                             known;
    public:
        /** constructor for initially empty array */
        RBoolArray();

        /** Add an FRIfaces::IRValueBool to the object.  Essentially does
            a push_back. Just takes a reference - does not copy or free. */
        void addArg(FRIfaces::IRValueBool* rValue);

        //// get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueBoolArray::RT* getRT();
        /** Is the value of this object known before the simulation starts
            eg a constant.  */
        virtual bool isKnown() const;

        virtual ~RBoolArray();
    };

    /** Class for ordinary double values */
    class PRODUCTS_DLL LValueDouble: public RValueDouble,
                        virtual public FRIfaces::ILValueDouble{
    public:
        class PRODUCTS_DLL RT{
            friend class LValueDouble;
        private:
            TGetValue*  func;
            TSetValue*  setFunc;
            char*       isSet;
            const char* name;
            double      value;
            /** generates exception saying either value not set yet or 
                already set */
            ModelException makeException();

        public:
            // constructor for no initial value
            explicit RT(char* isSet, const char* name);

            //// uses FR::MemMgr
            void* operator new(size_t size);
            void operator delete(void *ptr);

            static EDR_INLINE double getValue(void* structure){
                RT* rt = (RT*)structure;
                if (!*rt->isSet){
                    throw rt->makeException();
                }
                return rt->value;
            }
            
            // set the variable using a double
            static EDR_INLINE void setValue(void* structure, double value) {
                RT* rt = (RT*)structure;
                if (*rt->isSet){
                    throw rt->makeException();
                }
                rt->value = value;
                *rt->isSet = 1;
            }
        };
    private:
        RT*  rt;
        char isSet;
    public:

        ~LValueDouble();

        // constructor for no initial value
        LValueDouble(const char* name);
        //// allows central pooling of which variables are known
        virtual void setReset(char* reset);
        
        ////get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueDouble::RT* getRT();

        /** Is the value of this object known before the simulation starts
            eg a constant (This is useful as it allows some optimisations
            to be made). Here this is driven off whether the variable is set
            or not. This allows other variables to be precalculated if this
            variable can be */
        virtual bool isKnown() const;;

        /** sets the value of the variable at this time point. More of
            indicative method as derived types should supply method allowing
            easier extraction of variable (eg as double)  */
        virtual void set(const IObjectConstSP& object);

        // get the variable expressed as a double
        virtual double getValue();

        // set the variable using a double
        virtual void setValue(double value);

    };
    /** Double Array Variable - where we hold the values directly */
    class PRODUCTS_DLL LValueDoubleArray: public RValueDoubleArray,
                             virtual public FRIfaces::ILValueDoubleArray{
    private:
        struct RT;
        RT*    rt;
        char   isSet;
    public:

        ~LValueDoubleArray();

        // constructor for no initial value
        LValueDoubleArray(const char* name);
        //// allows central pooling of which variables are known
        virtual void setReset(char* reset);
        
        //// get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueDoubleArray::RT* getRT();

        /** Is the value of this object known before the simulation starts
            eg a constant (This is useful as it allows some optimisations
            to be made). Here this is driven off whether the variable is set
            or not. This allows other variables to be precalculated if this
            variable can be */
        virtual bool isKnown() const;;

        /** sets the value of the variable at this time point. More of
            indicative method as derived types should supply method allowing
            easier extraction of variable (eg as double)  */
        virtual void set(const IObjectConstSP& object);

    };
    /** Double Array Variable - where the values are held in an array
        of LValueDoubles */
    class PRODUCTS_DLL LValueDoubleVarArray: public RValueDoubleArray,
                                virtual public FRIfaces::ILValueDoubleArray{
    private:
        struct RT;
        RT*                                 rt;
        vector<FRIfaces::ILValueDouble*>    vectorLValues;
    public:

        ~LValueDoubleVarArray();

        //// constructor (pass in array of variables to use for array)
        LValueDoubleVarArray(
            const vector<FRIfaces::ILValueDouble*>& vectorLValues);

        LValueDoubleVarArray();
        //// allows central pooling of which variables are known
        virtual void setReset(char* reset);
        
        //// get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueDoubleArray::RT* getRT();

        /** Is the value of this object known before the simulation starts
            eg a constant (This is useful as it allows some optimisations
            to be made). Here this is driven off whether the variable is set
            or not. This allows other variables to be precalculated if this
            variable can be */
        virtual bool isKnown() const;;

        /** sets the value of the variable at this time point. More of
            indicative method as derived types should supply method allowing
            easier extraction of variable (eg as double)  */
        virtual void set(const IObjectConstSP& object);

    };

    /** Double Array Variable - where we hold the values directly */
    class PRODUCTS_DLL LValueIntArray: public RValueIntArray,
                             virtual public FRIfaces::ILValueIntArray{
    private:
        struct RT;
        RT*    rt;
        char   isSet;
    public:

        ~LValueIntArray();

        // constructor for no initial value
        LValueIntArray(const char* name);
        //// allows central pooling of which variables are known
        virtual void setReset(char* reset);
        
        //// get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueIntArray::RT* getRT();

        /** Is the value of this object known before the simulation starts
            eg a constant (This is useful as it allows some optimisations
            to be made). Here this is driven off whether the variable is set
            or not. This allows other variables to be precalculated if this
            variable can be */
        virtual bool isKnown() const;;

        /** sets the value of the variable at this time point. More of
            indicative method as derived types should supply method allowing
            easier extraction of variable (eg as int)  */
        virtual void set(const IObjectConstSP& object);

    };

    /** Int Array Variable - where the values are held in an array
        of LValueInts */
    class PRODUCTS_DLL LValueIntVarArray: public RValueIntArray,
                             virtual public FRIfaces::ILValueIntArray{
    private:
        struct RT;
        RT*                                 rt;
        vector<FRIfaces::ILValueInt*>    vectorLValues;
    public:

        ~LValueIntVarArray();

        //// constructor (pass in array of variables to use for array)
        LValueIntVarArray(
            const vector<FRIfaces::ILValueInt*>& vectorLValues);

        LValueIntVarArray();
        //// allows central pooling of which variables are known
        virtual void setReset(char* reset);
        
        //// get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueIntArray::RT* getRT();

        /** Is the value of this object known before the simulation starts
            eg a constant (This is useful as it allows some optimisations
            to be made). Here this is driven off whether the variable is set
            or not. This allows other variables to be precalculated if this
            variable can be */
        virtual bool isKnown() const;;

        /** sets the value of the variable at this time point. More of
            indicative method as derived types should supply method allowing
            easier extraction of variable (eg as int)  */
        virtual void set(const IObjectConstSP& object);

    };

    /** Double Array Variable - where we hold the values directly */
    class PRODUCTS_DLL LValueBoolArray: public RValueBoolArray,
                             virtual public FRIfaces::ILValueBoolArray{
    private:
        struct RT;
        RT*    rt;
        char   isSet;
    public:

        ~LValueBoolArray();

        // constructor for no initial value
        LValueBoolArray(const char* name);
        //// allows central pooling of which variables are known
        virtual void setReset(char* reset);
        
        //// get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueBoolArray::RT* getRT();

        /** Is the value of this object known before the simulation starts
            eg a constant (This is useful as it allows some optimisations
            to be made). Here this is driven off whether the variable is set
            or not. This allows other variables to be precalculated if this
            variable can be */
        virtual bool isKnown() const;;

        /** sets the value of the variable at this time point. More of
            indicative method as derived types should supply method allowing
            easier extraction of variable (eg as int)  */
        virtual void set(const IObjectConstSP& object);

    };

    /** Bool Array Variable - where the values are held in an array
        of LValueBools */
    class PRODUCTS_DLL LValueBoolVarArray: public RValueBoolArray,
                                virtual public FRIfaces::ILValueBoolArray{
    private:
        struct RT;
        RT*                                 rt;
        vector<FRIfaces::ILValueBool*>    vectorLValues;
    public:

        ~LValueBoolVarArray();

        //// constructor (pass in array of variables to use for array)
        LValueBoolVarArray(
            const vector<FRIfaces::ILValueBool*>& vectorLValues);

        LValueBoolVarArray();
        //// allows central pooling of which variables are known
        virtual void setReset(char* reset);
        
        //// get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueBoolArray::RT* getRT();

        /** Is the value of this object known before the simulation starts
            eg a constant (This is useful as it allows some optimisations
            to be made). Here this is driven off whether the variable is set
            or not. This allows other variables to be precalculated if this
            variable can be */
        virtual bool isKnown() const;;

        /** sets the value of the variable at this time point. More of
            indicative method as derived types should supply method allowing
            easier extraction of variable (eg as int)  */
        virtual void set(const IObjectConstSP& object);

    };

    /** Class for ordinary int values */
    class PRODUCTS_DLL LValueInt: public RValueInt,
                        virtual public FRIfaces::ILValueInt{
    public:
        class PRODUCTS_DLL RT{
            friend class LValueInt;
        private:
            TGetValue*  func;
            TSetValue*  setFunc;
            char*       isSet;
            const char* name;
            int         value;
            /** generates exception saying either value not set yet or 
                already set */
            ModelException makeException();

        public:
            // constructor for no initial value
            explicit RT(char* isSet, const char* name);

            //// uses FR::MemMgr
            void* operator new(size_t size);
            void operator delete(void *ptr);

            static EDR_INLINE int getValue(void* structure){
                RT* rt = (RT*)structure;
                if (!*rt->isSet){
                    throw rt->makeException();
                }
                return rt->value;
            }
            
            // set the variable using a int
            static EDR_INLINE void setValue(void* structure, int value) {
                RT* rt = (RT*)structure;
                if (*rt->isSet){
                    throw rt->makeException();
                }
                rt->value = value;
                *rt->isSet = 1;
            }
        };

    private:
        RT* rt;
        char isSet;
    public:

        ~LValueInt();

        // constructor for no initial value
        LValueInt(const char* name);
        //// allows central pooling of which variables are known
        virtual void setReset(char* reset);

        //// get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueInt::RT* getRT();

        /** Is the value of this object known before the simulation starts
            eg a constant (This is useful as it allows some optimisations
            to be made). Here this is driven off whether the variable is set
            or not. This allows other variables to be precalculated if this
            variable can be */
        virtual bool isKnown() const;;

        /** sets the value of the variable at this time point. More of
            indicative method as derived types should supply method allowing
            easier extraction of variable (eg as int)  */
        virtual void set(const IObjectConstSP& object);

        // get the variable expressed as a int
        virtual int getValue();

        // set the variable using a int
        virtual void setValue(int value);

    };

    /** Class for ordinary bool values */
    class PRODUCTS_DLL LValueBool: public RValueBool,
                      virtual public FRIfaces::ILValueBool{
    public:
        class PRODUCTS_DLL RT{
            friend class LValueBool;
        private:
            TGetValue*  func;
            TSetValue*  setFunc;
            char*       isSet;
            const char* name;
            bool        value;
            /** generates exception saying either value not set yet or 
                already set */
            ModelException makeException();

        public:
            // constructor for no initial value
            explicit RT(char* isSet, const char* name);

            //// uses FR::MemMgr
            void* operator new(size_t size);
            void operator delete(void *ptr);

            static EDR_INLINE bool getValue(void* structure){
                RT* rt = (RT*)structure;
                if (!*rt->isSet){
                    throw rt->makeException();
                }
                return rt->value;
            }
            
            // set the variable using a bool
            static EDR_INLINE void setValue(void* structure, bool value) {
                RT* rt = (RT*)structure;
                if (*rt->isSet){
                    throw rt->makeException();
                }
                rt->value = value;
                *rt->isSet = 1;
            }
        };

    private:
        RT* rt;
        char isSet;
    public:

        ~LValueBool();

        // constructor for no initial value
        LValueBool(const char* name);
        //// allows central pooling of which variables are known
        virtual void setReset(char* reset);

        //// get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueBool::RT* getRT();

        /** Is the value of this object known before the simulation starts
            eg a constant (This is useful as it allows some optimisations
            to be made). Here this is driven off whether the variable is set
            or not. This allows other variables to be precalculated if this
            variable can be */
        virtual bool isKnown() const;;

        /** sets the value of the variable at this time point. More of
            indicative method as derived types should supply method allowing
            easier extraction of variable (eg as bool)  */
        virtual void set(const IObjectConstSP& object);

        // get the variable expressed as a bool
        virtual bool getValue();

        // set the variable using a bool
        virtual void setValue(bool value);

    };


    /** Class for ordinary date values */
    class PRODUCTS_DLL LValueDate: public RValueDate,
                      virtual public FRIfaces::ILValueDate{
    public:
        struct PRODUCTS_DLL RT{
            TGetValue*             func;
            char*                  isSet;
            const char*            name;
            const DateTime::Date*  value;

            //// uses FR::MemMgr
            void* operator new(size_t size);
            void operator delete(void *ptr);
            /** generates exception saying either value not set yet or 
                already set */
            ModelException makeException();

            explicit RT(char* isSet, const char* name);

            static EDR_INLINE const DateTime::Date& getValue(void* structure){
                RT* rt = (RT*)structure;
                if (!*rt->isSet){
                    throw rt->makeException();
                }
                return *(rt->value);
            }
            // set the variable using a double
            static EDR_INLINE void setValue(void* structure, 
                                            const DateTime::Date& value) {
                RT* rt = (RT*)structure;
                if (*rt->isSet){
                    throw rt->makeException();
                }
                rt->value = &value;
                *rt->isSet = 1;
            }
        };
    private:
        RT* rt;
        char isSet;
    public:
        LValueDate(const char* name);
        
        ~LValueDate();
        //// allows central pooling of which variables are known
        virtual void setReset(char* reset);

        /** sets the value of the variable at this time point. More of
            indicative method as derived types should supply method allowing
            easier extraction of variable (eg as double)  */
        virtual void set(const IObjectConstSP& object);

        // get the variable expressed as a date
        virtual const DateTime::Date& getValue();

        // set the variable using a double
        virtual void setValue(const DateTime::Date& value);

        //// get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueDate::RT* getRT();
    };

    /** Base class for constant l-values. These need to implement the
        FRIfaces::ILConstValue interface and also fail for set operations */
    class PRODUCTS_DLL LValueConst: virtual public FRIfaces::ILConstValue{
    private:
        bool known;
    public:
        ~LValueConst();

        /** allows default value of isKnown() to be overridden for testing
            purposes */
        virtual void setIsValueKnown(bool isKnown);
        
        /** Is the value of this object known before the simulation starts
            eg a constant (This is useful as it allows some optimisations
            to be made) */
        virtual bool isKnown() const;
        //// allows central pooling of which variables are known
        virtual void setReset(char* reset);

        //// throws exception
        virtual void set(const IObjectConstSP& object);
    };

private:
    friend class FlexAlgorithm;
    friend class FlexInstrument;
    friend class Flex;
    class Instrument;
    typedef smartPtr<Instrument> InstrumentSP;
    class WriteRules;
    class InputIMS;
    class IMSConvert;
    class ProductMC;
    class ProductMCSV;
    class TPNoOp;
    class TPRSimpleArray;
    class TimePointRules;
    class TPRSparseArray;
    class Tester;
    friend class WriteRules;
    friend class InputIMS;
    friend class IMSConvert;
    typedef smartPtr<TimePointRules> TimePointRulesSP;
public:
    typedef array<TimePointRulesSP, TimePointRules> TimePointRulesArray;
    typedef smartPtr<TimePointRulesArray> TimePointRulesArraySP;
};

DRLIB_END_NAMESPACE
#endif
