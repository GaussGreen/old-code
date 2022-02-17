//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRParser.hpp
//
//   Description : Supports creating FR::RValues by parsing string
//                 based expressions
//
//   Author      : Mark A Robson
//
//   Date        : 28 June 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_FRPARSER_HPP
#define EDR_FRPARSER_HPP

#include "edginc/FRController.hpp"  
#include "edginc/FRFunction.hpp"  
#include "edginc/FR.hpp"  
#include "edginc/Maths.hpp"  

DRLIB_BEGIN_NAMESPACE
class PRODUCTS_DLL FRParser: public FR::RValueBase{
private:
    const string expression;
public:
    static CClassConstSP const TYPE;

    /** returns the id for this RValueBase */
    virtual const string& getID() const;

    /** Parses supplied string. Throws FRParseException for parse type errors.
        Object returned needs to be freed */
    static FRIfaces::IRValue* parseRValue(
        const char*   expression,
        FRController* frCtrl);
 
protected:
    /** creates FlexRule::IRValue which represents the value
        of this expression at the specific timepoint. The getRValue()
        uses this method and then saves it in the FRController. Return
        null if the save has already been done */
    virtual FRIfaces::IRValue* createRValue(
        int           index,
        FRController* frCtrl) const;

private:
    static void load(CClassSP& clazz);
    FRParser();
    static IObject* defaultConstructor();

public:
    // now all the template stuff - needs to be public so as to be visile
    // from C-style parser interface. Ideally it goes in this class but
    // I can't get it to compile if I do that...
    class Params;
    class FuncArgs;
};

// this part only needs to be seen by FRParser.cpp 
#ifdef EDR_FRPARSER_SOURCE_FILE

// work around for template bug in msvc 6
#if (defined(_MSC_VER) && (_MSC_VER < 1300))
#define FRPARSER_TEMPLATE_BUG_FIX 1
#endif
/** Template for adding different types together. First template parameter
    is the return type. The next two are the types of the two parameters */
template<class R, class T1 = R, class T2 = T1> struct AddOp{
    inline R operator()(T1 v1, T1 v2) const{
        return (v1+v2);
    }
};

//override for date + int
template <> struct AddOp<const DateTime::Date&, int>{
    mutable DateTime::Date result;
    inline const DateTime::Date& operator()(const DateTime::Date& v1,
                                            int v2)const{
        result = DateTime::Date(v1.getDate()+v2);
        return result;
    }
};

/** Template for MAXing different types. First template parameter
    is the return type. The next two are the types of the two parameters */
template<class R, class T1 = R, class T2 = T1> struct MaxOp{
    inline R operator()(T1 v1, T1 v2) const{
        return (v1 > v2? v1: v2);
    }
};

//override for dates
template <> struct MaxOp<const DateTime::Date&>{
    mutable DateTime::Date result;
    inline const DateTime::Date& operator()(const DateTime::Date& v1, 
                                            const DateTime::Date& v2) const{
        result = v1.getDate() > v2.getDate()? v1: v2;
        return result;
    }
};
    
/** Template for MINing different types. First template parameter
    is the return type. The next two are the types of the two parameters */
template<class R, class T1 = R, class T2 = T1> struct MinOp{
    inline R operator()(T1 v1, T1 v2) const{
        return (v1 < v2? v1: v2);
    }
};

//override for dates
template <> struct MinOp<const DateTime::Date&>{
    mutable DateTime::Date result;
    inline const DateTime::Date& operator()(const DateTime::Date& v1,
                                            const DateTime::Date& v2) const{
        result = v1.getDate() < v2.getDate()? v1: v2;
        return result;
    }
};

/** Template for subtracting different types from each other. First
    template parameter is the return type. The next two are the types
    of the two parameters */
template<class R, class T1 = R, class T2 = T1> struct SubOp{
    inline R operator()(T1 v1, T1 v2) const{
        return (v1-v2);
    }
};

//override for date - int
template <> struct SubOp<const DateTime::Date&, int>{
    mutable DateTime::Date result;
    inline const DateTime::Date& operator()(const DateTime::Date& v1, 
                                            int v2) const{
        result = DateTime::Date(v1.getDate()-v2);
        return result;
    }
};

//override for date - date
template <> struct SubOp<const DateTime::Date&>{
    inline int operator()(const DateTime::Date& v1, 
                          const DateTime::Date& v2) const{
        return (v1.getDate()-v2.getDate());
    }
};

/** Template for multipling types together. First
    template parameter is the return type. The next two are the types
    of the two parameters */
template<class R, class T1 = R, class T2 = T1> struct MultOp{
    inline R operator()(T1 v1, T1 v2) const{
        return (v1*v2);
    }
};

/** Template for multipling arrays together to get a scalar. First
    template parameter is the return type. The next two are the types
    of the two parameters */
template<class R, class T1 = R, class T2 = T1> struct ArrayMultOp{
    inline R operator()(T1 v1, T1 v2) const{
        int size1 = v1->size(v1);
        int size2 = v2->size(v2);
        if (size1 != size2){
            throw ModelException("Scalar product", "Arrays must be of equal"
                                 " length");
        }
        R result = 0;
        for (int i = 0; i < size1; i++){
            result += v1->func(v1, i) * v2->func(v2, i);
        }
        return result;
    }
};


/** Template for dividing types. First
    template parameter is the return type. The next two are the types
    of the two parameters */
template<class R, class T1 = R, class T2 = T1> struct DivOp{
    inline R operator()(T1 v1, T1 v2) const{
        if (v2 == 0){
            throw ModelException("Division by zero!");
        }
        return (v1/v2);
    }
};

/* For doing '^' */
class PowOp{
public:
    inline double operator()(double v1, double v2){
        return pow(v1, v2);
    }
};

/* for doing ! on bool expressions */
class NotOp{
public:
    inline bool operator()(bool p){
        return (!p);
    }
};

/** Template for testing for equality of two types. The two
    template parameter are the types of the two parameters */
template<class T1, class T2 = T1> struct EqOp{
    inline bool operator()(T1 v1, T1 v2) const{
        return (v1 == v2);
    }
};

//override for dates
template <> struct EqOp<const DateTime::Date&>{
    inline bool operator()(const DateTime::Date& v1,
                           const DateTime::Date& v2) const{
        return (v1.equals(v2));
    }
};

/** Template for testing for equality of two types. The two
    template parameter are the types of the two parameters */
template<class T1, class T2 = T1> struct NotEqOp{
    inline bool operator()(T1 v1, T1 v2) const{
        return (v1 != v2);
    }
};

//override for dates
template <> struct NotEqOp<const DateTime::Date&>{
    inline bool operator()(const DateTime::Date& v1,
                           const DateTime::Date& v2) const{
        return (!v1.equals(v2));
    }
};


/** Template for testing for ">" of two types. The two
    template parameter are the types of the two parameters */
template<class T1, class T2 = T1> struct GTOp{
    inline bool operator()(T1 v1, T1 v2) const{
        return (v1 > v2);
    }
};

//override for dates
template <> struct GTOp<const DateTime::Date&>{
    inline bool operator()(const DateTime::Date& v1,
                           const DateTime::Date& v2) const{
        return (v1.isGreater(v2));
    }
};

/** Template for testing for "<" of two types. The two
    template parameter are the types of the two parameters */
template<class T1, class T2 = T1> struct LTOp{
    inline bool operator()(T1 v1, T1 v2) const{
        return (v1 < v2);
    }
};

//override for dates
template <> struct LTOp<const DateTime::Date&>{
    inline bool operator()(const DateTime::Date& v1,
                           const DateTime::Date& v2) const{
        return (!v1.isGreaterOrEqual(v2));
    }
};

/** Template for doing negation ie '-' */
template<class T> struct NegOp{
    inline T operator()(T v) const{
        return (-v);
    }
};

/** Template for doing MAX or MIN on an array, Use GTOp for MAX and
    LTOp for MIN */
template<class R, class RT, class Op> struct MaxMinArrayOp{
    Op op;
    inline R operator()(RT rt) const{
        int size = rt->size(rt);
        if (size == 0){
            throw ModelException("MaxArrayOp", "Array is empty");
        }
        typedef  R (RTGetValue)(void* structure, int index);
        RTGetValue* getValue = rt->func;
        R maxMinVal = getValue(rt, 0);
        for (int i = 1; i < size; i++){
            R value = getValue(rt, i);
            if (op(value, maxMinVal)){
                maxMinVal = value;
            }
        }
        return maxMinVal;
    }
};

/** For pretending that '[]' operation on arrays is a binary operation
    RT: the RT class eg FRIfaces::IRValueDoubleArray::RT
    X: the return type of getValue method on an element in the array */
template<class RT, class X> struct VariableArrayIndexOperator{
    inline X operator()(RT* rtArray, int index) const{
        int size = rtArray->size(rtArray);
        if (index < 0 || index >= size){
            throw ModelException("operator []",
                                 "Attempt to read element with index "+
                                 Format::toString(index)+
                                 " of array of length "+
                                 Format::toString(size));
        }
        return (rtArray->func(rtArray, index));
    }
};

/** template function to get value. First template parameter is the
    variable type (eg a pointer to a class). The second is the return
    type ie what does that variable evaluate to (eg a double) */
template <class V, class T> inline T value(V v){
    return (v->func(v));
}

// specialization for doubles
template<> inline double value<double, double>(double v){
    return v;
}

// specialization for ints
template<> inline int value<int, int>(int v){
    return v;
}

// specialization for bools
template<> inline bool value<bool, bool>(bool v){
    return v;
}

// specialization for dates
template<> inline const DateTime::Date& value<const DateTime::Date&, 
    const DateTime::Date&>(const DateTime::Date& v){
    return v;
}

// specialization for FR::LValueDouble::RT - quicker access
template<> inline double value<FR::LValueDouble::RT*, 
                               double>(FR::LValueDouble::RT* v){
    return FR::LValueDouble::RT::getValue(v);
}

// specialization for FR::LValueInt - quicker access
template<> inline int value<FR::LValueInt::RT*, int>(FR::LValueInt::RT* v){
    return FR::LValueInt::RT::getValue(v);
}

// specialization for FR::LValueBool - quicker access
template<> inline bool value<FR::LValueBool::RT*, bool>(FR::LValueBool::RT* v){
    return FR::LValueBool::RT::getValue(v);
}

// specialization for double arrays - we just return ourselves. Allows binary
// operator template to be used for '[]' operator
template<> inline FRIfaces::IRValueDoubleArray::RT*
value<FRIfaces::IRValueDoubleArray::RT*, 
    FRIfaces::IRValueDoubleArray::RT*>(FRIfaces::IRValueDoubleArray::RT* v){
    return v;
}

// specialization for bool arrays - we just return ourselves. Allows binary
// operator template to be used for '[]' operator
template<> inline FRIfaces::IRValueBoolArray::RT*
value<FRIfaces::IRValueBoolArray::RT*, 
    FRIfaces::IRValueBoolArray::RT*>(FRIfaces::IRValueBoolArray::RT* v){
    return v;
}

// specialization for int arrays - we just return ourselves. Allows binary
// operator template to be used for '[]' operator
template<> inline FRIfaces::IRValueIntArray::RT*
value<FRIfaces::IRValueIntArray::RT*, 
    FRIfaces::IRValueIntArray::RT*>(FRIfaces::IRValueIntArray::RT* v){
    return v;
}

/** Template class for generic unary function. Template parameters:
    I: what the class derives from (eg RValueDouble)
    R: the return type of the evaulation eg double
    V: Type of the parameter (eg RValueInt*)
    OP: the operation to perform. Uses operator (.) on OP
    RT: the return type of getRT()
    R2: the return type to use when doing value<,> on V
*/
template <class I, class R, class V, 
    class OP, class RT, class R2 = R> class Unary: public I{
private:
    struct MyRT{
        typename I::TGetValue* func;
        V                      var;
        OP                     op;
    public:
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }

        explicit MyRT(V var): func(&getValue), var(var){}

        static R getValue(void* structure){
            MyRT* rt = (MyRT*)structure;
            return (rt->op(value<V, R2>(rt->var)));
        }
    };
    MyRT* rt;
 public:
    explicit Unary(V var): rt(new MyRT(var)){}
    ~Unary(){
        delete rt;
    }

    virtual R getValue(){
        return MyRT::getValue(rt);
    }

    virtual RT* getRT(){
        return (RT*)rt;
    }
};

/** Template class for generic 2 factor binary operations.
    Template parameters:
    I: what the class derives from (eg RValueDouble)
    R: the return type of the evaulation eg double
    V1: Type of the first parameter (eg RValueInt::RT*)
    V2: Type of the second parameter (eg RValueInt::RT*)
    OP: the operation to perform. Uses operator (.,.) on OP
    RT: the return type of getRT()
    V1T: What variable1 evaluates to (eg int)
    V2T: What variable2 evaluates to (eg int) 
*/
template <class I, class R, class V1, class V2, 
          class OP, class RT, class V1T = R, class V2T=R>
class Binary: public I{
    struct MyRT{
        typename I::TGetValue* func;
        V1                     var1;
        V2                     var2;
        OP                     op;
    public:
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(V1 var1, V2 var2): 
            func(&getValue), var1(var1), var2(var2){}
    
        static R getValue(void* structure){
            MyRT* rt = (MyRT*)structure;
            return rt->op(value<V1, V1T>(rt->var1), value<V2, V2T>(rt->var2));
        }
    };
    MyRT* rt;
public:
    explicit Binary(V1 var1, V2 var2): rt(new MyRT(var1, var2)){}
    virtual RT* getRT(){
        return (RT*)rt;
    }
    virtual R getValue(){
        return MyRT::getValue(rt);
    }
    ~Binary(){
        delete rt;
    }
};

/** Template class for generic "if then else". Template parameters:
    I: what the class derives from (eg RValueDouble)
    RT: the return type of getRT()
    C: the type of the condition parameter (eg IRValueBool* or bool)
    R: the return type of the evaulation eg double
    V1: Type of the first parameter (eg RValueInt*)
    V2: Type of the first parameter (eg RValueInt*)
*/
template <class I, class RT, class C, class R, class V1, 
          class V2> 
class Condition: public I{
    struct MyRT{
        typename I::TGetValue* func;
        C                      cond;
        V1                     ifTrue;
        V2                     ifFalse;
    public:
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(C cond, V1 ifTrue, V2 ifFalse): 
            func(&getValue), cond(cond), ifTrue(ifTrue), ifFalse(ifFalse){}
        static R getValue(void* structure){
            MyRT* rt = (MyRT*)structure;
            return (value<C, bool>(rt->cond)?
                    value<V1, R>(rt->ifTrue): value<V2, R>(rt->ifFalse));
        }
    };
    MyRT* rt;
public:
    explicit Condition(C cond, V1 ifTrue, V2 ifFalse): 
        rt(new MyRT(cond, ifTrue, ifFalse)){}
    
    virtual RT* getRT(){
        return (RT*)rt;
    }
    virtual R getValue(){
        return MyRT::getValue(rt);
    }
    ~Condition(){
        delete rt;
    }
};

/** Template class for generic "if then else" for arrays. 
    Template parameters:
    I: what the class derives from (eg RValueDouble)
    RT: the return type of getRT()
    R: the return type of the evaluation eg double
    V: Type of the first and second parameter (eg IRValueDoubleArray)
*/
template <class I, class RT, class R, class V> 
class ArrayCondition: public I{
    typedef typename V::RT VRT;
    struct MyRT{
        typename I::TGetSize*           size;
        typename I::TGetValue*          func;
        FRIfaces::IRValueBoolArray::RT* cond;
        VRT*                            ifTrue;
        VRT*                            ifFalse;
    public:
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(FRIfaces::IRValueBoolArray::RT* cond, 
                      VRT* ifTrue, VRT* ifFalse): 
            size(&getSize), func(&getValue), 
            cond(cond), ifTrue(ifTrue), ifFalse(ifFalse){}

        static int getSize(void* structure){
            MyRT* rt = (MyRT*)structure;
            int size1 = rt->ifTrue->size(rt->ifTrue);
            int size2 = rt->ifFalse->size(rt->ifFalse);
            int size3 = rt->cond->size(rt->cond);
            if (size1 != size2 || size1 != size3){
                throw ModelException("Condition Array operator must act on "
                                     "arrays of equal length");
            }
            return size1;
        }

        static R getValue(void* structure, int index){
            MyRT* rt = (MyRT*)structure;
            return (rt->cond->func(rt->cond, index)?
                    rt->ifTrue->func(rt->ifTrue, index):
                    rt->ifFalse->func(rt->ifFalse, index));
        }
    };
    FRIfaces::IRValueBoolArray* cond;
    V*                          ifTrue;
    V*                          ifFalse;
    MyRT*                       rt;
public:
    explicit ArrayCondition(FRIfaces::IRValueBoolArray* cond, 
                            V*                          ifTrue, 
                            V*                          ifFalse):
        cond(cond), ifTrue(ifTrue), ifFalse(ifFalse),
        rt(new MyRT(cond->getRT(), ifTrue->getRT(), ifFalse->getRT())){}
    
    virtual RT* getRT(){
        return (RT*)rt;
    }
    /** Cheap way of implementing this - better performance if we precompute
        the result and store in a separate class. But this can wait */
    virtual bool isKnown() const{
        return (cond->isKnown() && ifTrue->isKnown() && ifFalse->isKnown());
    }
    ~ArrayCondition(){
        delete rt;
    }
};

/** Template class for generic "if then else" for arrays. 
    Template parameters:
    I: what the class derives from (eg RValueDouble)
    RT: the return type of getRT()
    R: the return type of the evaluation eg double
    V1: Type of the first [array] parameter (eg IRValueDoubleArray)
    V2: Type of the second [scalar] parameter (eg IRValueDouble)
    negateCond: Whether condition should be negated
*/
template <class I, class RT, class R, class V1, class V2, bool negateCond> 
class ArrayConditionUsingScalar: public I{
    typedef typename V1::RT V1RT;
    typedef typename V2::RT V2RT;
    struct MyRT{
        typename I::TGetSize*           size;
        typename I::TGetValue*          func;
        FRIfaces::IRValueBoolArray::RT* cond;
        V1RT*                           ifTrue;
        V2RT*                           ifFalse;
    public:
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(FRIfaces::IRValueBoolArray::RT* cond, 
                      V1RT* ifTrue, V2RT* ifFalse): 
            size(&getSize), func(&getValue), 
            cond(cond), ifTrue(ifTrue), ifFalse(ifFalse) {}

        static int getSize(void* structure){
            MyRT* rt = (MyRT*)structure;
            int size1 = rt->ifTrue->size(rt->ifTrue);
            int size2 = rt->cond->size(rt->cond);
            if (size1 != size2){
                throw ModelException("Condition Array operator must act on "
                                     "arrays of equal length");
            }
            return size1;
        }

        static R getValue(void* structure, int index){
            MyRT* rt = (MyRT*)structure;
            // note: logical exclusive 'or' here - allows us to reuse
            // this template for array? scalar: array
            return ((negateCond ^ rt->cond->func(rt->cond, index))?
                    rt->ifTrue->func(rt->ifTrue, index):
                    rt->ifFalse->func(rt->ifFalse));
        }
    };
    FRIfaces::IRValueBoolArray* cond;
    V1*                         ifTrue;
    V2*                         ifFalse;
    MyRT*                       rt;
public:
    /** The negateCond parameter allows use to use the template for
        array? scalar: array (you need to swap ifTrue and ifFalse over) */
    explicit ArrayConditionUsingScalar(
        FRIfaces::IRValueBoolArray* cond, 
        V1*                         ifTrue, 
        V2*                         ifFalse):
        cond(cond), ifTrue(ifTrue), ifFalse(ifFalse),
        rt(new MyRT(cond->getRT(), ifTrue->getRT(), ifFalse->getRT())){}
    
    virtual RT* getRT(){
        return (RT*)rt;
    }
    /** Cheap way of implementing this - better performance if we precompute
        the result and store in a separate class. But this can wait */
    virtual bool isKnown() const{
        return (cond->isKnown() && ifTrue->isKnown() && ifFalse->isKnown());
    }
    ~ArrayConditionUsingScalar(){
        delete rt;
    }
};

/** Template class for generic "if then else" for arrays where the if is on
    a bool. 
    Template parameters:
    I: what the class derives from (eg RValueDouble)
    RT: the return type of getRT()
    R: the return type of the evaluation eg double
    V: Type of the first and second parameter (eg IRValueDoubleArray)
*/
template <class I, class RT, class R, class V> 
class ArrayConditionFromBool: public I{
    typedef typename V::RT VRT;
    struct MyRT{
        typename I::TGetSize*           size;
        typename I::TGetValue*          func;
        FRIfaces::IRValueBool::RT*      cond;
        VRT*                            ifTrue;
        VRT*                            ifFalse;
    public:
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(FRIfaces::IRValueBool::RT* cond, 
                      VRT* ifTrue, VRT* ifFalse): 
            size(&getSize), func(&getValue), 
            cond(cond), ifTrue(ifTrue), ifFalse(ifFalse){}

        static int getSize(void* structure){
            MyRT* rt = (MyRT*)structure;
            int size1 = rt->ifTrue->size(rt->ifTrue);
            int size2 = rt->ifFalse->size(rt->ifFalse);
            if (size1 != size2){
                throw ModelException("Condition Array operator must act on "
                                     "arrays of equal length");
            }
            return size1;
        }

        static R getValue(void* structure, int index){
            MyRT* rt = (MyRT*)structure;
            return (rt->cond->func(rt->cond)?
                    rt->ifTrue->func(rt->ifTrue, index):
                    rt->ifFalse->func(rt->ifFalse, index));
        }
    };
    FRIfaces::IRValueBool*      cond;
    V*                          ifTrue;
    V*                          ifFalse;
    MyRT*                       rt;
public:
    explicit ArrayConditionFromBool(FRIfaces::IRValueBool*      cond, 
                                    V*                          ifTrue, 
                                    V*                          ifFalse):
        cond(cond), ifTrue(ifTrue), ifFalse(ifFalse),
        rt(new MyRT(cond->getRT(), ifTrue->getRT(), ifFalse->getRT())){}
    
    virtual RT* getRT(){
        return (RT*)rt;
    }
    /** Cheap way of implementing this - better performance if we precompute
        the result and store in a separate class. But this can wait */
    virtual bool isKnown() const{
        return (cond->isKnown() && ifTrue->isKnown() && ifFalse->isKnown());
    }
    ~ArrayConditionFromBool(){
        delete rt;
    }
};

/** Template function for building binary 'operator'
    class R1: the type of the return value eg FRIfaces::IRValueBool
    class I: the class that the Binary class should derive from
    class R2: the type that the return of getValue on p1 should be 
    assigned to 
    class R2x: the same as R2 except applies to p2
    class R3: the type that getValue() returns on the returned object
    class C: the type that can be used to make an instance of R1 
    using type op(R2,R2) 
    class Op: the binary operator
    class Alt1: type to try to cast p1 to for performance
    class Alt2: type to try to cast p2 to for performance
    class P1: the type of the first operand (eg FRIfaces::IRValueBool*)
    class P2: the type of the second operand (eg FRIfaces::IRValueBool*)
*/
template <class R1, class I, class R2, class R2x, class R3,
          class C, class Op, class Alt1, class Alt2, class P1, 
          class P2> R1* createBinaryFunc(P1* p1, P2* p2
#ifdef FRPARSER_TEMPLATE_BUG_FIX // work around for template bug
                                         , R1* d1 = 0, I* d2 = 0, 
                                         R2* d3 = 0,
                                         R2x* d7 = 0,
                                         R3* d4 = 0,
                                         C* d5 = 0, Op* d6 = 0,
                                         Alt1* d8 = 0, Alt1* d9 = 0
#endif
              ){
    typedef typename R1::RT R1RT;
    typedef typename P1::RT P1RT;
    typedef typename P2::RT P2RT;
    typedef typename Alt1::RT Alt1RT;
    typedef typename Alt2::RT Alt2RT;
    bool p1Known = p1->isKnown();
    bool p2Known = p2->isKnown();
    Alt1* p1Alt = dynamic_cast<Alt1*>(p1);
    Alt2* p2Alt = dynamic_cast<Alt2*>(p2);
    if (p1Known){
        if (p2Known){
            Op op;
            return new C(op(p1->getValue(), p2->getValue()));
        }
        if (p2Alt){
            return new Binary<I, R3, R2, Alt2RT*, 
                Op, R1RT, R2, R2x>(p1->getValue(), (Alt2RT*)p2Alt->getRT());
        } else {
            return new Binary<I, R3, R2, P2RT*, 
                Op, R1RT, R2, R2x>(p1->getValue(), p2->getRT());
        }
    }
    if (p2Known){
        if (p1Alt){
            return new Binary<I, R3, Alt1RT*, R2x, 
                Op, R1RT, R2, R2x>((Alt1RT*)p1Alt->getRT(), p2->getValue());
        } else {
            return new Binary<I, R3, P1RT*, R2x, Op,
                R1RT, R2, R2x>(p1->getRT(), p2->getValue());
        }
    }
    if (p1Alt){
        if (p2Alt){
            return new Binary<I, R3, Alt1RT*, Alt2RT*,
                Op, R1RT, R2, R2x>((Alt1RT*)p1Alt->getRT(), 
                                   (Alt2RT*)p2Alt->getRT());
        } else {
            return new Binary<I, R3, Alt1RT*, P2RT*, 
                Op, R1RT, R2, R2x>((Alt1RT*)p1Alt->getRT(), p2->getRT());
        }
    } else if (p2Alt){
        return new Binary<I, R3, P1RT*, Alt2RT*,
            Op, R1RT, R2, R2x>(p1->getRT(), (Alt2RT*)p2Alt->getRT());
    }
    // default
    return new Binary<I, R3, P1RT*, P2RT*, Op, R1RT, R2, R2x>(p1->getRT(), 
                                                              p2->getRT());
}


/** Template function for building "?.: ." 'operator'
    class R1: the type of the return value eg FRIfaces::IRValueBool*
    class I: the class that the Condition class should derive from
    class R2: the type that getValue() returns on the object created
    class Alt: type to try to cast p1 or p2 to for performance
    class P1: the type of the first operand (eg FRIfaces::IRValueBool*)
    class P2: the type of the second operand (eg FRIfaces::IRValueBool*)
*/
template <class R1, class I, class R2, class Alt, class P1, class P2> 
static R1* createConditionFunc(FRIfaces::IRValueBool* cond, P1* p1, P2* p2
#ifdef FRPARSER_TEMPLATE_BUG_FIX // work around for template bug
                              ,R1* d1 = 0, I* d2 = 0, R2* d3 = 0,
                               Alt* d4 =0
#endif
    ){
    typedef typename R1::RT R1RT;
    typedef typename P1::RT P1RT;
    typedef typename P2::RT P2RT;
    typedef typename Alt::RT AltRT;
    if (cond->isKnown()){
        return (cond->getValue()? p1: p2);
    }
    FRIfaces::IRValueBool::RT* condRT = cond->getRT();
    bool p1Known = p1->isKnown();
    bool p2Known = p2->isKnown();
    P1RT* p1RT = p1->getRT();
    P2RT* p2RT = p2->getRT();
    FR::LValueBool* condAlt = dynamic_cast<FR::LValueBool*>(cond);
    FR::LValueBool::RT* condAltRT = condAlt? 
        (FR::LValueBool::RT*)condAlt->getRT(): 0;
    Alt* p1Alt = dynamic_cast<Alt*>(p1);
    AltRT* p1AltRT = p1Alt? (AltRT*)p1Alt->getRT(): 0;
    Alt* p2Alt = dynamic_cast<Alt*>(p2);
    AltRT* p2AltRT = p2Alt? (AltRT*)p2Alt->getRT(): 0;
    if (p1Known){
        R2 p1Val = p1->getValue();
        if (p2Known){
            R2 p2Val = p2->getValue();
            if (condAlt){
                return new 
                    Condition<I, R1RT, FR::LValueBool::RT*, R2,
                    R2, R2>(condAltRT, p1Val, p2Val);
            } else {
                return new 
                    Condition<I, R1RT, FRIfaces::IRValueBool::RT*, R2,
                    R2, R2>(condRT, p1Val, p2Val);
            }
        }
        if (condAlt){
            if (p2Alt){
                return new Condition<I, R1RT, FR::LValueBool::RT*, R2, R2, 
                    AltRT*> (condAltRT, p1->getValue(), p2AltRT);
            } else {
                return new Condition<I, R1RT, FR::LValueBool::RT*, R2, R2,
                    P2RT*> (condAltRT, p1Val, p2RT);
            }
        } else {
            if (p2Alt){
                return new Condition<I, R1RT, FRIfaces::IRValueBool::RT*,
                    R2, R2, AltRT*>(condRT, p1Val, p2AltRT);
            } else {
                return new Condition<I, R1RT, FRIfaces::IRValueBool::RT*,
                    R2, R2,P2RT*> (condRT, p1Val, p2RT);
            }
        }
    }
    if (p2Known){
        R2 p2Val = p2->getValue();
        if (condAlt){
            if (p1Alt){
                return new Condition<I, R1RT, FR::LValueBool::RT*, R2, AltRT*,
                    R2> (condAltRT, p1AltRT, p2Val);
            } else {
                return new Condition<I, R1RT, FR::LValueBool::RT*, R2,P1RT*,R2>
                    (condAltRT, p1RT, p2Val);
            }
        } else {
            if (p1Alt){
                return new Condition<I, R1RT, FRIfaces::IRValueBool::RT*, R2, 
                    AltRT*, R2> (condRT, p1AltRT, p2Val);
            } else {
                return new Condition<I, R1RT, FRIfaces::IRValueBool::RT*, R2, 
                    P1RT*, R2> (condRT, p1RT, p2Val);
            }
        }
    }
    if (condAlt){
        if (p1Alt){
            if (p2Alt){
                return new Condition<I, R1RT, FR::LValueBool::RT*, R2, AltRT*, 
                    AltRT*>(condAltRT, p1AltRT, p2AltRT);
            } else {
                return new Condition<I, R1RT, FR::LValueBool::RT*, R2, AltRT*, 
                    P2RT*>(condAltRT, p1AltRT, p2RT);
            }
        } else {
            if (p2Alt){
                return new Condition<I, R1RT, FR::LValueBool::RT*, R2, P1RT*, 
                    AltRT*>(condAltRT, p1RT, p2AltRT);
            } else {
                return new Condition<I, R1RT, FR::LValueBool::RT*, R2, P1RT*, 
                    P2RT*>(condAltRT, p1RT, p2RT);
            }
        }
    } else {
        if (p1Alt){
            if (p2Alt){
                return new Condition<I, R1RT, FRIfaces::IRValueBool::RT*, R2,
                    AltRT*, AltRT*>(condRT, p1AltRT, p2AltRT);
            } else {
                return new Condition<I, R1RT, FRIfaces::IRValueBool::RT*, R2, 
                    AltRT*, P2RT*>(condRT, p1AltRT, p2RT);
            }
        } else {
            if (p2Alt){
                return new Condition<I, R1RT, FRIfaces::IRValueBool::RT*, R2, 
                    P1RT*, AltRT*>(condRT, p1RT, p2AltRT);
            } else {
                return new Condition<I, R1RT, FRIfaces::IRValueBool::RT*, R2,
                    P1RT*, P2RT*>(condRT, p1RT, p2RT);
            }
        }
    }
}

template <class Op> FRIfaces::IRValueBool* createBinaryBoolFunc(
    FRIfaces::IRValueBool* p1, 
    FRIfaces::IRValueBool* p2
#ifdef FRPARSER_TEMPLATE_BUG_FIX // work around for template bug
    , Op*                    d1 = 0
#endif
    ){
    return createBinaryFunc<FRIfaces::IRValueBool, FR::RValueBool, 
        bool, bool, bool, FR::RConstBool, Op, FR::LValueBool, 
        FR::LValueBool>(p1, p2);
}

template <class Op> FRIfaces::IRValueBool* createBinaryInt2BoolFunc(
    FRIfaces::IRValueInt* p1, 
    FRIfaces::IRValueInt* p2
#ifdef FRPARSER_TEMPLATE_BUG_FIX // work around for template bug
    , Op*                    d1 = 0
#endif
    ){
    return createBinaryFunc<FRIfaces::IRValueBool, FR::RValueBool, 
        int, int, bool, FR::RConstBool, Op,
        FR::LValueInt, FR::LValueInt>(p1, p2);
}

template <class Op>
FRIfaces::IRValueBool* createBinaryDouble2BoolFunc(
    FRIfaces::IRValueDouble* p1, 
    FRIfaces::IRValueDouble* p2
#ifdef FRPARSER_TEMPLATE_BUG_FIX // work around for template bug
    , Op*                    d1 = 0
#endif
    ){
    return createBinaryFunc<FRIfaces::IRValueBool, FR::RValueBool, 
        double, double, bool, FR::RConstBool, Op,
        FR::LValueDouble, FR::LValueDouble>(p1, p2);
}
template <class Op>
FRIfaces::IRValueBool* createBinaryDate2BoolFunc(
    FRIfaces::IRValueDate* p1, 
    FRIfaces::IRValueDate* p2
#ifdef FRPARSER_TEMPLATE_BUG_FIX // work around for template bug
    , Op*                    d1 = 0
#endif
    ){
    return createBinaryFunc<FRIfaces::IRValueBool, FR::RValueBool, 
        const DateTime::Date&, const DateTime::Date&, bool, 
        FR::RConstBool, Op, 
        FRIfaces::IRValueDate /* to do */, FRIfaces::IRValueDate>(p1, p2);
}

template <class Op> FRIfaces::IRValueInt* createBinaryIntFunc(
    FRIfaces::IRValueInt* p1, 
    FRIfaces::IRValueInt* p2
#ifdef FRPARSER_TEMPLATE_BUG_FIX // work around for template bug
    , Op*                    d1 = 0
#endif
    ){
    return createBinaryFunc<FRIfaces::IRValueInt, FR::RValueInt, 
        int, int, int, FR::RConstInt, Op, 
        FR::LValueInt, FR::LValueInt>(p1, p2);
}

template <class Op, class P1, class P2>
FRIfaces::IRValueDouble* createBinaryDoubleFunc(P1 p1, P2 p2
#ifdef FRPARSER_TEMPLATE_BUG_FIX // work around for template bug
                                                , Op*   d1 = 0
#endif
    ){
    return createBinaryFunc<FRIfaces::IRValueDouble, FR::RValueDouble, 
        double, double, double, FR::RConstDouble, Op, 
        FR::LValueDouble, FR::LValueDouble>(p1, p2);
}

// sort of specialization of createBinaryDoubleFunc where p1 might be a
// FR::LValueInt,
template <class Op, class P1, class P2>
FRIfaces::IRValueDouble* createBinaryDoubleFuncInt1(P1 p1, P2 p2
#ifdef FRPARSER_TEMPLATE_BUG_FIX // work around for template bug
                                                , Op*   d1 = 0
#endif
    ){
    return createBinaryFunc<FRIfaces::IRValueDouble, FR::RValueDouble, 
        int, double, double, FR::RConstDouble, Op, 
        FR::LValueInt, FR::LValueDouble>(p1, p2);
}

// sort of specialization of createBinaryDoubleFunc where p2 might be a
// FR::LValueInt,
template <class Op, class P1, class P2>
FRIfaces::IRValueDouble* createBinaryDoubleFuncInt2(P1 p1, P2 p2
#ifdef FRPARSER_TEMPLATE_BUG_FIX // work around for template bug
                                                , Op*   d1 = 0
#endif
    ){
    return createBinaryFunc<FRIfaces::IRValueDouble, FR::RValueDouble, 
        double, int, double, FR::RConstDouble, Op, 
        FR::LValueDouble, FR::LValueInt>(p1, p2);
}

/** class P1V and P2V are the types of getValue on P1 and P2 */
template <class Op, class P1V, class P2V, class P1, class P2>
FRIfaces::IRValueDate* createBinaryDateFunc(P1 p1, P2 p2
#ifdef FRPARSER_TEMPLATE_BUG_FIX // work around for template bug
                                            , P1V* d1 = 0
                                            , P2V* d2 = 0
                                            , Op*  d3 = 0
#endif
    ){
    return createBinaryFunc<FRIfaces::IRValueDate, FR::RValueDate, 
        P1V, P2V, const DateTime::Date&, FR::RConstDate, Op,
        FRIfaces::IRValueDate /* to do */, FRIfaces::IRValueDate>(p1, p2);
}

/** class P1V and P2V are the types of getValue on P1 and P2 */
template <class Op, class P1V, class P2V, class P1, class P2>
FRIfaces::IRValueDate* createBinaryDateFuncInt1(P1 p1, P2 p2
#ifdef FRPARSER_TEMPLATE_BUG_FIX // work around for template bug
                                            , P1V* d1 = 0
                                            , P2V* d2 = 0
                                            , Op*  d3 = 0
#endif
    ){
    return createBinaryFunc<FRIfaces::IRValueDate, FR::RValueDate, 
        P1V, P2V, const DateTime::Date&, FR::RConstDate, Op,
        FR::LValueInt, FRIfaces::IRValueDate /* to do */>(p1, p2);
}

/** class P1V and P2V are the types of getValue on P1 and P2 */
template <class Op, class P1V, class P2V, class P1, class P2>
FRIfaces::IRValueDate* createBinaryDateFuncInt2(P1 p1, P2 p2
#ifdef FRPARSER_TEMPLATE_BUG_FIX // work around for template bug
                                            , P1V* d1 = 0
                                            , P2V* d2 = 0
                                            , Op*  d3 = 0
#endif
    ){
    return createBinaryFunc<FRIfaces::IRValueDate, FR::RValueDate, 
        P1V, P2V, const DateTime::Date&, FR::RConstDate, Op,
        FRIfaces::IRValueDate /* to do */, FR::LValueInt>(p1, p2);
}

static FRIfaces::IRValueBool* createBoolCondition(
    FRIfaces::IRValueBool* cond,
    FRIfaces::IRValueBool* p1,
    FRIfaces::IRValueBool* p2){
    return createConditionFunc<FRIfaces::IRValueBool, FR::RValueBool,
        bool, FR::LValueBool>(cond, p1, p2);
}

static FRIfaces::IRValueInt* createIntCondition(
    FRIfaces::IRValueBool* cond,
    FRIfaces::IRValueInt*  p1,
    FRIfaces::IRValueInt*  p2){
    return createConditionFunc<FRIfaces::IRValueInt, FR::RValueInt,
        int, FR::LValueInt>(cond, p1, p2);
}

static FRIfaces::IRValueDouble* createDoubleCondition(
    FRIfaces::IRValueBool*   cond,
    FRIfaces::IRValueDouble* p1,
    FRIfaces::IRValueDouble* p2){
    return createConditionFunc<FRIfaces::IRValueDouble, FR::RValueDouble,
        double, FR::LValueDouble>(cond, p1, p2);
}

/** Template class for '&&' operation. Template parameters:
    V1: Type of the first parameter (eg RValueBool*)
    V2: Type of the first parameter (eg RValueBool*)
    We can't use the Binary class since we want to avoid redundant 
    evaluations
*/
template <class V1, class V2> class And: public FR::RValueBool{
private:
    struct MyRT{
        TGetValue* func;
        V1         var1;
        V2         var2;
    public:
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(V1 var1, V2 var2): 
            func(&getValue), var1(var1), var2(var2){}

        static bool getValue(void* structure){
            MyRT* rt = (MyRT*)structure;
            return (value<V1, bool>(rt->var1) && value<V2, bool>(rt->var2));
        }
    };
    MyRT* rt;
public:
    explicit And(V1 var1, V2 var2): rt(new MyRT(var1, var2)){}
    
    ~And(){
        delete rt;
    }

    virtual bool getValue(){
        return MyRT::getValue(rt);
    }
    virtual FRIfaces::IRValueBool::RT* getRT(){
        return (FRIfaces::IRValueBool::RT*)rt;
    }
};

/** Template class for '||' operation. Template parameters:
    V1: Type of the first parameter (eg RValueBool*)
    V2: Type of the first parameter (eg RValueBool*)
    We can't use the Binary class since we want to avoid redundant 
    evaluations
*/
template <class V1, class V2> class Or: public FR::RValueBool{
private:
    struct MyRT{
        TGetValue* func;
        V1         var1;
        V2         var2;
    public:
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(V1 var1, V2 var2): 
            func(&getValue), var1(var1), var2(var2){}
        static bool getValue(void* structure){
            MyRT* rt = (MyRT*)structure;
            return (value<V1, bool>(rt->var1) || value<V2, bool>(rt->var2));
        }
    };
    MyRT* rt;
public:
    explicit Or(V1 var1, V2 var2): rt(new MyRT(var1, var2)){}
    
    ~Or(){
        delete rt;
    }

    virtual bool getValue(){
        return MyRT::getValue(rt);
    }
    virtual FRIfaces::IRValueBool::RT* getRT(){
        return (FRIfaces::IRValueBool::RT*)rt;
    }
};

/** function for doing && on two bools. */
static FR::RValueBool* createBoolAnd(FRIfaces::IRValueBool* p1, 
                                     FRIfaces::IRValueBool* p2){
    bool p1Known = p1->isKnown();
    bool b1 = p1Known? p1->getValue(): false; // avoid compiler warning
    bool p2Known = p2->isKnown();
    bool b2 = p2Known? p2->getValue(): false; // avoid compiler warning
    if (p1Known){
        if (p2Known){
            return new FR::RConstBool(b1 && b2);
        }
        return new And<bool, FRIfaces::IRValueBool::RT*>(b1, p2->getRT());
    }
    if (p2Known){
        return new And<FRIfaces::IRValueBool::RT*, bool>(p1->getRT(), b2);
    }
    return new And<FRIfaces::IRValueBool::RT*, 
        FRIfaces::IRValueBool::RT*>(p1->getRT(), p2->getRT());
}

/** function for doing || on two bools. */
static FR::RValueBool* createBoolOr(FRIfaces::IRValueBool* p1, 
                                    FRIfaces::IRValueBool* p2){
    bool p1Known = p1->isKnown();
    bool b1 = p1Known? p1->getValue(): false; // avoid compiler warning
    bool p2Known = p2->isKnown();
    bool b2 = p2Known? p2->getValue(): false; // avoid compiler warning
    if (p1Known){
        if (p2Known){
            return new FR::RConstBool(b1 || b2);
        }
        return new Or<bool, FRIfaces::IRValueBool::RT*>(b1, p2->getRT());
    }
    if (p2Known){
        return new Or<FRIfaces::IRValueBool::RT*, bool>(p1->getRT(), b2);
    }
    return new Or<FRIfaces::IRValueBool::RT*, 
        FRIfaces::IRValueBool::RT*>(p1->getRT(), p2->getRT());
}

/** template function to call getValue() on right field of FRIfaces::RValUnion
    V is the type to return eg bool, double 
    R is a bogus class to get around a bug in MSVC templates */
template <class V, class R> inline V invokeGetValue(
    const FRIfaces::RValUnion&  theExp,
    R* dummy1 = 0){
    /* if you're here with a compile error you need to write a
       specialisation */
    return theExp.bExp->func(theExp.bExp);
}

// specialization for doubles
template<> inline double invokeGetValue<double, FR::RValueDouble>(
    const FRIfaces::RValUnion&  theExp,
    FR::RValueDouble*       d1){
    return theExp.dExp->func(theExp.dExp);
}

// specialization for ints
template<> inline int invokeGetValue<int, FR::RValueInt>(
    const FRIfaces::RValUnion&  theExp,
    FR::RValueInt*          d1){
    return theExp.iExp->func(theExp.iExp);
}

// specialization for bools
template<> inline bool invokeGetValue<bool, FR::RValueBool>(
    const FRIfaces::RValUnion&  theExp,
    FR::RValueBool*         d1){
    return theExp.bExp->func(theExp.bExp);
}

// specialization for dates
template<> inline const DateTime::Date& invokeGetValue<const DateTime::Date&,
    FR::RValueDate>(
    const FRIfaces::RValUnion&  theExp,
    FR::RValueDate*         d1){
    return theExp.dtExp->func(theExp.dtExp);
}

/** template function to call getValue() on right field of FRIfaces::RValUnion
    for arrays.
    V is the type to return eg bool, double 
    R is a bogus class to get around a bug in MSVC templates */
template <class V, class R> inline V invokeArrayGetValue(
    const FRIfaces::RValUnion&  theExp,
    int                         index,
    R* dummy1 = 0){
    /* if you're here with a compile error you need to write a
       specialisation */
    return theExp.bArrayExp->func(theExp.bExp, index);
}

// specialization for doubles
template<> inline double invokeArrayGetValue<double, FR::RValueDoubleArray>(
    const FRIfaces::RValUnion&  theExp,
    int                         index,
    FR::RValueDoubleArray*       d1){
    return theExp.dArrayExp->func(theExp.dExp, index);
}

// specialization for ints
template<> inline int invokeArrayGetValue<int, FR::RValueIntArray>(
    const FRIfaces::RValUnion&  theExp,
    int                         index,
    FR::RValueIntArray*         d1){
    return theExp.iArrayExp->func(theExp.iExp, index);
}

// specialization for bools
template<> inline bool invokeArrayGetValue<bool, FR::RValueBoolArray>(
    const FRIfaces::RValUnion&  theExp,
    int                         index,
    FR::RValueBoolArray*        d1){
    return theExp.bArrayExp->func(theExp.bExp, index);
}

/** Template function for creating unary operator
    class I: the class that Unary should derive from eg FR::RValueInt
    class R: the primitive type of P eg int
    class C: the class to use for a constant value of type R
    class Op: the operator to use for the unary function
    class Alt: type to try to cast p to for performance
    class P: the type of the parameter eg FRIfaces::IRValueInt */
template <class I, class R, class C, class Op, class Alt, class P>
P* createUnaryFunc(P* p
#ifdef FRPARSER_TEMPLATE_BUG_FIX // work around for template bug
                  , I* d1 = 0, R* d2 = 0,
                  C* d3 = 0, Op* d4 = 0,
                  Alt* d5 = 0
#endif
    ){
    typedef typename P::RT PRT;
    typedef typename Alt::RT AltRT;
    if (p->isKnown()){
        Op op;
        return new C(op(p->getValue()));
    }
    Alt* pAlt = dynamic_cast<Alt*>(p);
    if (pAlt){
        return new Unary<I, R, AltRT*, Op, PRT>((AltRT*)pAlt->getRT());
    }
    return new Unary<I, R, PRT*, Op, PRT>(p->getRT());
}
/** Template class for doing [] on variables of different types
    Template parameters:
    I: what the class derives from (eg RValueDouble)
    RT: the return type of getRT()
    R: the return type of the evaulation eg double
    V: the type of the expression eg FRIfaces::IRValueDouble
*/
template<class I, class R, class V> class IndexOperator: public I{
    struct MyRT{
        typename I::TGetValue*                 func;
        const vector<FRIfaces::RValUnion>*     rValues;
        const string*                          id; // of variable
        FRIfaces::IRValueInt*                  indexExp;

        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }

        explicit MyRT(const vector<FRIfaces::RValUnion>*  rValues,
                      const string&                       id,
                      FRIfaces::IRValueInt*               indexExp):
            func(&getValue), rValues(rValues), id(&id), indexExp(indexExp){}

        static R getValue(void* structure){
            MyRT* rt = (MyRT*)structure;
            int index = rt->indexExp->getValue();
            if (index < 0 || index >= (int) rt->rValues->size()){
                throw ModelException("FRParser::IndexOperator", 
                                     "Index for '"+*rt->id+"' ("+
                                     Format::toString(index)+
                                     ") is out of bounds");
            }
            const FRIfaces::RValUnion& rValue = (*rt->rValues)[index];
            return invokeGetValue<R, I>(rValue, 0);          
        }
    };
    MyRT*  rt;
public:
    IndexOperator(
        const vector<FRIfaces::RValUnion>*     rValues,
        const string&                      id, // of variable
        FRIfaces::IRValueInt*              indexExp): 
        rt(new MyRT(rValues, id, indexExp)){}
    
    ~IndexOperator(){
        delete rt;
    }

    virtual R getValue(){
        return MyRT::getValue(rt);
    }

    virtual typename V::RT* getRT(){
        return (typename V::RT*)rt;
    }
    
};
        
/** template function for handling  [] operator 
    Template parameters:
    I: what the class created should derive from (eg RValueDouble)
    R: the return type of the evaulation eg double
    V: the type of the expression eg FRIfaces::IRValueDouble
*/
template<class I, class R, class V>  V* createIndexFunc(
    FRIfaces::ILValueExpression* var,
    FRIfaces::IRValueInt*        rValue,
    FRController*                ctrl
#ifdef FRPARSER_TEMPLATE_BUG_FIX // work around for template bug
                                               , I* d1 = 0, R* d2 = 0, 
                                               V* d3 = 0
#endif
){
    if (rValue->isKnown()){
        try{
            int index = rValue->getValue();
            FRIfaces::IRValue* rValue = var->getRValue(index, ctrl);
            // must be a FRIfaces::IRValueBool*
            V& val = dynamic_cast<V&>(*rValue);
            return &val; // see comment next to return statement below
        } catch (exception&){
            // may be this expression is never evaluated so defer until
            // run time
        }
    }
    const vector<FRIfaces::RValUnion>& rValues = 
        ctrl->getRValuesForLValueExpression(var);
    /* Note: caller can use FRController::store on the return value of this
       function since the object obtained from getRValue is stored in the 
       same place so will not be freed twice */
    return new IndexOperator<I, R, V>(&rValues, var->getID(), rValue);
}

/** Template function for doing [] on arrays. Here the [] refers not to the
    time dimension but to the array's dimension at a specific time point 
    I: what the class derives from (eg FR::RValueDouble)
    R: the return type of the evaulation eg double
    V: the RT type of the variable eg FRIfaces::IRValueDoubleArray::RT
    RT: the return type of getRT() eg FRIfaces::IRValueDouble::RT
*/
template<class I, class R, class V, class RT> I* createArrayBracketOperator(
        V*                                rt1,
        FRIfaces::IRValueInt::RT*         rt2
#ifdef FRPARSER_TEMPLATE_BUG_FIX // work around for template bug
        , I* d1 = 0, R* d2 = 0, V* d3 = 0
#endif
    )
{
    return new  Binary<I, R, V*, FRIfaces::IRValueInt::RT*, 
        VariableArrayIndexOperator<V, R>, RT,
        V*, int>(rt1, rt2);
#if 0
    return new Binary<FR::RValueDouble, double,
        FRIfaces::IRValueDoubleArray::RT*, FRIfaces::IRValueInt::RT*, 
        VariableArrayIndexOperator<FRIfaces::IRValueDoubleArray::RT, 
        double>, FRIfaces::IRValueDouble::RT,
        FRIfaces::IRValueDoubleArray::RT*, int>(rt1, rt2);
#endif
}
/** Template class for doing [] on array variables of different types
    (here the [] refers to the time dimension)
    Template parameters:
    I: what the class derives from (eg RValueDoubleArray)
    RT: the return type of getRT()
    R: the return type of the evaulation eg double
    V: the type of the expression eg FRIfaces::IRValueDouble
*/
template<class I, class R, class V> class ArrayIndexOperator: public I{
    struct MyRT{
        typename I::TGetSize*                  size;
        typename I::TGetValue*                 func;
        const vector<FRIfaces::RValUnion>*     rValues;
        const string*                          id; // of variable
        FRIfaces::IRValueInt*                  indexExp;

        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        //// returns the index at which to look up the variable
        int getTimeIndex(){
            int index = indexExp->getValue();
            if (index < 0 || index >= (int) rValues->size()){
                throw ModelException("FRParser::ArrayIndexOperator", 
                                     "Index for '"+*id+"' ("+
                                     Format::toString(index)+
                                     ") is out of bounds");
            }
            return index;
        }
        //// returns the size of the array
        static int getSize(void* structure){
            MyRT* rt = (MyRT*)structure;
            int timeIndex = rt->getTimeIndex();
            // to determine length of array we can pull out any of the
            // array pointers (since they all have the size method as the
            // first element in the struct)
            // (This is a bit hacky - the fix is to define 
            // FRIfaces::IRValueArray class and its RT counterpart)
            FRIfaces::IRValueDoubleArray::RT* rt2 = 
                (FRIfaces::IRValueDoubleArray::RT*) 
                (*rt->rValues)[timeIndex].dArrayExp;
            return (rt2->size(rt2));
        }

        explicit MyRT(const vector<FRIfaces::RValUnion>*  rValues,
                      const string&                       id,
                      FRIfaces::IRValueInt*               indexExp):
            size(&getSize), func(&getValue), 
                 rValues(rValues), id(&id), indexExp(indexExp){}

        static R getValue(void* structure, int index){
            MyRT* rt = (MyRT*)structure;
            int timeIndex = rt->getTimeIndex();
            const FRIfaces::RValUnion& rValue = (*rt->rValues)[timeIndex];
            return (invokeArrayGetValue<R, I>(rValue, index, 0));
        }
    };
    MyRT*  rt;
public:
    ArrayIndexOperator(
        const vector<FRIfaces::RValUnion>*     rValues,
        const string&                      id, // of variable
        FRIfaces::IRValueInt*              indexExp): 
        rt(new MyRT(rValues, id, indexExp)){}
    
    ~ArrayIndexOperator(){
        delete rt;
    }

    virtual typename V::RT* getRT(){
        return (typename V::RT*)rt;
    }
    
};
        
/** template function for handling  [] operator 
    Template parameters:
    I: what the class created should derive from (eg RValueDouble)
    R: the return type of the evaulation eg double
    V: the type of the expression eg FRIfaces::IRValueDouble
*/
template<class I, class R, class V>  V* createArrayIndexFunc(
    FRIfaces::ILValueExpression* var,
    FRIfaces::IRValueInt*        rValue,
    FRController*                ctrl
#ifdef FRPARSER_TEMPLATE_BUG_FIX // work around for template bug
                                               , I* d1 = 0, R* d2 = 0, 
                                               V* d3 = 0
#endif
){
    if (rValue->isKnown()){
        try{
            int index = rValue->getValue();
            FRIfaces::IRValue* rValue = var->getRValue(index, ctrl);
            // must be a FRIfaces::IRValueBool*
            V& val = dynamic_cast<V&>(*rValue);
            return &val; // see comment next to return statement below
        } catch (exception&){
            // may be this expression is never evaluated so defer until
            // run time
        }
    }
    const vector<FRIfaces::RValUnion>& rValues = 
        ctrl->getRValuesForLValueExpression(var);
    /* Note: caller can use FRController::store on the return value of this
       function since the object obtained from getRValue is stored in the 
       same place so will not be freed twice */
    // to do: can we generalise createIndexFunc so as to avoid this
    // template?
    return new ArrayIndexOperator<I, R, V>(&rValues, var->getID(), rValue);
}

/* class to 'cast' ints to doubles */
class FRIntToDouble: public FR::RValueDouble{
    struct MyRT{
        FR::RValueDouble::TGetValue* func;
        FRIfaces::IRValueInt::RT*    intValue;
        
        explicit MyRT(FRIfaces::IRValueInt::RT* intValue):
            func(&getValue), intValue(intValue){}
        static double getValue(void* structure){
            MyRT* rt = (MyRT*)structure;
            return (rt->intValue->func(rt->intValue));
        }
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
    };
    FRIfaces::IRValueInt* intValue;
    MyRT*                 rt;
public:
    FRIntToDouble(FRIfaces::IRValueInt* intValue): 
        intValue(intValue), rt(new MyRT(intValue->getRT())){}

    ~FRIntToDouble(){
        delete rt;
    }
    //// pass through to int
    virtual bool isKnown() const{
        return intValue->isKnown();
    }

    virtual double getValue(){
        return MyRT::getValue(rt);
    }

    virtual FRIfaces::IRValueDouble::RT* getRT(){
        return (FRIfaces::IRValueDouble::RT*)rt;
    }
};

/** Template class for generic unary operations which acts on
    an array and returns an array.
    Template parameters:
    I: what the class derives from (eg RValueDouble)
    R: the return type of the evaluation eg double
    V: Type of the first parameter (eg RValueInt)
    OP: the operation to perform. Uses operator (.) on OP
    RT: the return type of getRT()
*/
template <class I, class R, class V, class OP, class RT>
class ArrayUnary: public I{
    typedef typename V::RT VRT;
    struct MyRT{
        typename I::TGetSize*  size;
        typename I::TGetValue* func;
        VRT*                   var;
        OP                     op;
    public:
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(VRT* var): 
            size(&getSize), func(&getValue), var(var) {}

        static int getSize(void* structure){
            MyRT* rt = (MyRT*)structure;
            return rt->var->size(rt->var);
        }
        static R getValue(void* structure, int index){
            MyRT* rt = (MyRT*)structure;
            return rt->op(rt->var->func(rt->var, index));
        }
    };
    V*    v;
    MyRT* rt;
public:
    explicit ArrayUnary(V* v): v(v), rt(new MyRT(v->getRT())){}
    virtual RT* getRT(){
        return (RT*)rt;
    }

    /** Cheap way of implementing this - better performance if we precompute
        the result and store in a separate class. But this can wait */
    virtual bool isKnown() const{
        return (v->isKnown());
    }

    ~ArrayUnary(){
        delete rt;
    }
};

/** Template class for generic 2 factor binary operations which acts on
    arrays and returns an array.
    Template parameters:
    I: what the class derives from (eg RValueDouble)
    R: the return type of the evaluation eg double
    V1: Type of the first parameter (eg RValueInt*)
    V2: Type of the second parameter (eg RValueInt*)
    OP: the operation to perform. Uses operator (.,.) on OP
    RT: the return type of getRT()
*/
template <class I, class R, class V1, class V2, 
          class OP, class RT>
class ArrayBinary: public I{
    typedef typename V1::RT V1RT;
    typedef typename V2::RT V2RT;
    struct MyRT{
        typename I::TGetSize*  size;
        typename I::TGetValue* func;
        V1RT*                  var1;
        V2RT*                  var2;
        OP                     op;
    public:
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(V1RT* var1, V2RT* var2): 
            size(&getSize), func(&getValue), var1(var1), var2(var2){}

        static int getSize(void* structure){
            MyRT* rt = (MyRT*)structure;
            int size1 = rt->var1->size(rt->var1);
            int size2 = rt->var2->size(rt->var2);
            if (size1 != size2){
                throw ModelException("Binary Array operator must act on "
                                     "arrays of equal length");
            }
            return size1;
        }
        static R getValue(void* structure, int index){
            MyRT* rt = (MyRT*)structure;
            return rt->op(rt->var1->func(rt->var1, index),
                          rt->var2->func(rt->var2, index));
        }
    };
    V1*   v1;
    V2*   v2;
    MyRT* rt;
public:
    explicit ArrayBinary(V1* v1, V2* v2): v1(v1), v2(v2),
                                          rt(new MyRT(v1->getRT(), 
                                                      v2->getRT())){}
    virtual RT* getRT(){
        return (RT*)rt;
    }

    /** Cheap way of implementing this - better performance if we precompute
        the result and store in a separate class. But this can wait */
    virtual bool isKnown() const{
        return (v1->isKnown() && v2->isKnown());
    }

    ~ArrayBinary(){
        delete rt;
    }
};

/** 'OR' Operator for arrays of bools. Cannot use ArrayBinary template
    class because we want to redundant evaluations */
class ArrayOr: public FR::RValueBoolArray{
    struct MyRT{
        FR::RValueBoolArray::TGetSize*     size;
        FR::RValueBoolArray::TGetValue*    func;
        FRIfaces::IRValueBoolArray::RT*    var1;
        FRIfaces::IRValueBoolArray::RT*    var2;
    public:
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(FRIfaces::IRValueBoolArray::RT* var1,
                      FRIfaces::IRValueBoolArray::RT* var2): 
            size(&getSize), func(&getValue), var1(var1), var2(var2){}

        static int getSize(void* structure){
            MyRT* rt = (MyRT*)structure;
            int size1 = rt->var1->size(rt->var1);
            int size2 = rt->var2->size(rt->var2);
            if (size1 != size2){
                throw ModelException("Array 'OR' operator must act on "
                                     "arrays of equal length");
            }
            return size1;
        }
        static bool getValue(void* structure, int index){
            MyRT* rt = (MyRT*)structure;
            return (rt->var1->func(rt->var1, index) ||
                    rt->var2->func(rt->var2, index));
        }
    };
    FRIfaces::IRValueBoolArray*   v1;
    FRIfaces::IRValueBoolArray*   v2;
    MyRT*                         rt;
public:
    explicit ArrayOr(FRIfaces::IRValueBoolArray* v1, 
                     FRIfaces::IRValueBoolArray* v2): 
        v1(v1), v2(v2), rt(new MyRT(v1->getRT(), 
                                    v2->getRT())){}
    virtual RT* getRT(){
        return (RT*)rt;
    }

    /** Cheap way of implementing this - better performance if we precompute
        the result and store in a separate class. But this can wait */
    virtual bool isKnown() const{
        return (v1->isKnown() && v2->isKnown());
    }

    ~ArrayOr(){
        delete rt;
    }
};

/** 'AND' Operator for arrays of bools. Cannot use ArrayBinary template
    class because we want to redundant evaluations */
class ArrayAnd: public FR::RValueBoolArray{
    struct MyRT{
        FR::RValueBoolArray::TGetSize*     size;
        FR::RValueBoolArray::TGetValue*    func;
        FRIfaces::IRValueBoolArray::RT*    var1;
        FRIfaces::IRValueBoolArray::RT*    var2;
    public:
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(FRIfaces::IRValueBoolArray::RT* var1,
                      FRIfaces::IRValueBoolArray::RT* var2): 
            size(&getSize), func(&getValue), var1(var1), var2(var2){}

        static int getSize(void* structure){
            MyRT* rt = (MyRT*)structure;
            int size1 = rt->var1->size(rt->var1);
            int size2 = rt->var2->size(rt->var2);
            if (size1 != size2){
                throw ModelException("Array 'AND' operator must act on "
                                     "arrays of equal length");
            }
            return size1;
        }
        static bool getValue(void* structure, int index){
            MyRT* rt = (MyRT*)structure;
            return (rt->var1->func(rt->var1, index) &&
                    rt->var2->func(rt->var2, index));
        }
    };
    FRIfaces::IRValueBoolArray*   v1;
    FRIfaces::IRValueBoolArray*   v2;
    MyRT*                         rt;
public:
    explicit ArrayAnd(FRIfaces::IRValueBoolArray* v1, 
                     FRIfaces::IRValueBoolArray* v2): 
        v1(v1), v2(v2), rt(new MyRT(v1->getRT(), 
                                    v2->getRT())){}
    virtual RT* getRT(){
        return (RT*)rt;
    }

    /** Cheap way of implementing this - better performance if we precompute
        the result and store in a separate class. But this can wait */
    virtual bool isKnown() const{
        return (v1->isKnown() && v2->isKnown());
    }

    ~ArrayAnd(){
        delete rt;
    }
};

/** Template class for generic 2 factor binary operations which acts on
    an array (1st parameter) and a scalar and returns an array.
    Template parameters:
    I: what the class derives from (eg RValueDouble)
    R: the return type of the evaluation eg double
    V1: Type of the first parameter (eg RValueInt*)
    V2: Type of the second parameter (eg RValueInt*)
    OP: the operation to perform. Uses operator (.,.) on OP
    RT: the return type of getRT()
*/
template <class I, class R, class V1, class V2, 
          class OP, class RT>
class ArrayBinaryOnScalar: public I{
    typedef typename V1::RT V1RT;
    typedef typename V2::RT V2RT;
    struct MyRT{
        typename I::TGetSize*  size;
        typename I::TGetValue* func;
        V1RT*                  var1;
        V2RT*                  var2;
        OP                     op;
    public:
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(V1RT* var1, V2RT* var2): 
            size(&getSize), func(&getValue), var1(var1), var2(var2){}

        static int getSize(void* structure){
            MyRT* rt = (MyRT*)structure;
            return rt->var1->size(rt->var1);
        }
        static R getValue(void* structure, int index){
            MyRT* rt = (MyRT*)structure;
            return rt->op(rt->var1->func(rt->var1, index),
                          rt->var2->func(rt->var2));
        }
    };
    V1*   v1;
    V2*   v2;
    MyRT* rt;
public:
    explicit ArrayBinaryOnScalar(V1* v1, V2* v2): v1(v1), v2(v2),
                                                  rt(new MyRT(v1->getRT(), 
                                                              v2->getRT())){}
    virtual RT* getRT(){
        return (RT*)rt;
    }

    /** Cheap way of implementing this - better performance if we precompute
        the result and store in a separate class. But this can wait */
    virtual bool isKnown() const{
        return (v1->isKnown() && v2->isKnown());
    }

    ~ArrayBinaryOnScalar(){
        delete rt;
    }
};

/** Template class for selecting subsection of an array.
    Template parameters:
    I: what the class derives from (eg FR::RValueDoubleArray)
    R: the return type of the evaluation eg double
    V: Type of the array parameter (eg FRIfaces::IRValueDoubleArray)
    RT: the return type of getRT()
*/
template <class I, class R, class V, class RT>
class SubArray: public I{
    typedef typename V::RT VRT;
    struct MyRT{
        typename I::TGetSize*      size;
        typename I::TGetValue*     func;
        VRT*                       var;
        int                        cachedStart;
        int                        cachedLength;
        FRIfaces::IRValueInt::RT*  start;
        FRIfaces::IRValueInt::RT*  length;
    public:
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(VRT*                      var, 
                      FRIfaces::IRValueInt::RT* start, 
                      FRIfaces::IRValueInt::RT* length): 
            size(&getSize), func(&getValue), var(var), 
            start(start), length(length) {}

        static int getSize(void* structure){
            static const string method("SubArray::getSize");
            // getSize must be invoked before getValue so we can simply 
            // cache the start and length here
            MyRT* rt = (MyRT*)structure;
            rt->cachedStart = rt->start->func(rt->start);
            if (rt->cachedStart < 0){
                throw ModelException(method, "Index to start "
                                     "array must be >= 0");
            }
            rt->cachedLength = rt->length->func(rt->length);
            int actualArrayLength = rt->var->size(rt->var);
            if (rt->cachedStart + rt->cachedLength > actualArrayLength){
                string m("Requesting an array of length "+
                         Format::toString(rt->cachedLength)+ 
                         " starting at index "+
                         Format::toString(rt->cachedStart)+" from "
                         "array of length "+
                         Format::toString(actualArrayLength));
                throw ModelException(method, m);
            }
            return rt->cachedLength;
        }
        static R getValue(void* structure, int index){
            MyRT* rt = (MyRT*)structure;
            return rt->var->func(rt->var, rt->cachedStart+index);
        }
    };
    V*                     v;
    FRIfaces::IRValueInt*  start;
    FRIfaces::IRValueInt*  length;
    MyRT*                  rt;
public:
    explicit SubArray(V*                     v, 
                      FRIfaces::IRValueInt*  start,
                      FRIfaces::IRValueInt*  length):
        v(v), start(start), length(length),
        rt(new MyRT(v->getRT(), start->getRT(), length->getRT())){}
    virtual RT* getRT(){
        return (RT*)rt;
    }

    /** Cheap way of implementing this - better performance if we precompute
        the result and store in a separate class. But this can wait */
    virtual bool isKnown() const{
        return (v->isKnown() && start->isKnown() && length->isKnown());
    }

    ~SubArray(){
        delete rt;
    }
};

/** Internal representation of function - holds incoming variables plus
    any precomputed values. This is the base class  */
class FRInternalFunction: virtual public FRIfaces::IRValue{
public:
    class RTData; // definition below
private:
    FRController*             ctrl;
    const FRFunction*         function;
    vector<FRIfaces::VarType> expTypes; // their types
    BoolArray                 expIsArray; // is the expression an array
    const FRFunction::ValCalcType*    calcType; // how to evaluate args
    RTData*                   rtData;
    bool                      rtDataNeeded; // is RTData needed
    bool                      allArgsKnown;
public:
    /** Data needed at runtime by 'RT' classes */
    class RTData{
    public:
        int                       numParams;
        FRFunction::Val*          values; // array
        vector<FRFunction::Exp>   expressions;
        const FRIfaces::VarType*  reqTypes; // what they're supposed to be

        ~RTData(){
            FR::MemMgr::dealloc(values);
        }

        RTData(const FRFunction*    function): 
            numParams(function->numArgs()),
            values((FRFunction::Val*)FR::MemMgr::alloc(
                (numParams+1)* sizeof(FRFunction::Val))),
            reqTypes(function->argTypes()){}

        //// Populate the values (needs to be invoked at run time) 
        //// unless preCompute says everything has been preprocessed
        void evaluateArgs(){
            for (int i = 0; i < numParams; i++){
                if (expressions[i].dExp){// they're all pointers so any will do
                    // not yet calculated
                    switch (reqTypes[i]){
                    case FRIfaces::doubleType:
                        values[i].d = expressions[i].dExp->getValue();
                        break;
                    case FRIfaces::intType:
                        values[i].i = expressions[i].iExp->getValue();
                        break;
                    case FRIfaces::boolType:
                        values[i].b = expressions[i].bExp->getValue();
                        break;
                    case FRIfaces::doubleArrayType:
                    case FRIfaces::intArrayType:
                    case FRIfaces::boolArrayType:
                        // nothing to do
                        break;
                    case FRIfaces::dateType:
                        values[i].dt = &(expressions[i].dtExp->getValue());
                        break;
                    case FRIfaces::scheduleType:
                        values[i].sched = expressions[i].schedExp->getValue();
                        break;
                    case FRIfaces::tabulatedFuncType:
                        values[i].tabFunc =
                            expressions[i].tabFuncExp->getValue();
                        break;
                    default:
                        throw ModelException("Function::getValue",
                                             "Type not supported");
                    }
                }
            }
        }
    };

    //// little utitlity method - appends to expTypes and expIsArray arrays
    void recordTypeAndArg(FRIfaces::VarType type, 
                          bool              isArray, 
                          FRFunction::Exp&  exp){
        expTypes.push_back(type);
        expIsArray.push_back(isArray);
        rtData->expressions.push_back(exp);
    }
        
public:
    /** verify number and types of arguments - also does any casting of
        parameters needed (eg int to double) */
    void verifyArgs(){
        string wrongNumArgs;
        string wrongType;
        if (rtData->numParams != (int) expTypes.size()){
            wrongNumArgs = Format::toString((int)expTypes.size())+
                " parameter(s) "+
                "supplied to function '"+function->name()+
                "' but function takes "+
                Format::toString(rtData->numParams)+" parameters\n";
        }
        for (int i = 0; i < Maths::min(rtData->numParams, expTypes.size());
             i++){
            bool arrayReq = calcType && 
                (calcType[i] == FRFunction::arrayExpression);
            // check types of param including whether array or not
            if (expIsArray[i] != arrayReq || 
                rtData->reqTypes[i] != expTypes[i]){
                // possible type mismatch - first check for any casts
                if (!expIsArray[i] &&
                    expTypes[i] == FRIfaces::intType &&
                    rtData->reqTypes[i] == FRIfaces::doubleType){
                    rtData->expressions[i].dExp = 
                        new FRIntToDouble(rtData->expressions[i].iExp);
                    // avoid memory leak
                    ctrl->store(rtData->expressions[i].dExp);
                } else {
                    wrongType += "Incorrect type for parameter "+
                        Format::toString(i+1)+" for function "+
                        function->name()+". Function requires type "+
                        FRFunction::argTypeToString(rtData->reqTypes[i],
                                                    arrayReq)+
                        " but parameter is of type "+
                        FRFunction::argTypeToString(expTypes[i]
                                                    ,expIsArray[i])+ "\n";
                }
            }
        }
        wrongNumArgs += wrongType;
        if (!wrongNumArgs.empty()){
            FRParseException e(wrongNumArgs);
            e.addMsg("FRParser::FuncArgs::verifyArgs");
            throw e;
        }
    }
    
    /** Precompute values where possible - returns true if all values
        could be precomputed. (Also sets allKnown if function can be
        evaluated immediately and rtDataNeeded to true if no RTData needed) */
    bool  preCompute(){
        int numPopulated = 0; // how many times expressions is populated
        int numDone = 0; // how many time isKnown is true
        for (int i = 0; i < rtData->numParams; i++){
            bool evalArg = calcType? 
                calcType[i] == FRFunction::native: true;
            if (calcType && calcType[i]  == FRFunction::arrayExpression){
                rtData->values[i].rValues = rtData->expressions[i].rValues;
            } else {
                switch (rtData->reqTypes[i]){
                case FRIfaces::doubleType:
                    rtData->values[i].dExp = 
                        rtData->expressions[i].dExp->getRT();
                    if (rtData->expressions[i].dExp->isKnown()){
                        if (evalArg){
                            rtData->values[i].d =
                                rtData->expressions[i].dExp->getValue();
                            rtData->expressions[i].dExp = 0; // clear
                        }
                        numDone++;
                    }
                    break;
                case FRIfaces::intType:
                    rtData->values[i].iExp = 
                        rtData->expressions[i].iExp->getRT();
                    if (rtData->expressions[i].iExp->isKnown()){
                        if (evalArg){
                            rtData->values[i].i = 
                                rtData->expressions[i].iExp->getValue();
                            rtData->expressions[i].iExp = 0; // clear
                        }
                        numDone++;
                    }
                    break;
                case FRIfaces::boolType:
                    rtData->values[i].bExp = 
                        rtData->expressions[i].bExp->getRT();
                    if (rtData->expressions[i].bExp->isKnown()){
                        if (evalArg){
                            rtData->values[i].b = 
                                rtData->expressions[i].bExp->getValue();
                            rtData->expressions[i].bExp = 0; // clear
                        }
                        numDone++;
                    }
                    break;
                case FRIfaces::doubleArrayType:
                    rtData->values[i].dExpArray = 
                        rtData->expressions[i].dExpArray->getRT();
                    if (rtData->expressions[i].dExpArray->isKnown()){
                        numDone++;
                    }
                    break;
                case FRIfaces::intArrayType:
                    rtData->values[i].iExpArray = 
                        rtData->expressions[i].iExpArray->getRT();
                    if (rtData->expressions[i].iExpArray->isKnown()){
                        numDone++;
                    }
                    break;
                case FRIfaces::boolArrayType:
                    rtData->values[i].bExpArray = 
                        rtData->expressions[i].bExpArray->getRT();
                    if (rtData->expressions[i].bExpArray->isKnown()){
                        numDone++;
                    }
                    break;
                case FRIfaces::dateType:
                    rtData->values[i].dtExp = 
                        rtData->expressions[i].dtExp->getRT();
                    if (rtData->expressions[i].dtExp->isKnown()){
                        if (evalArg){
                            rtData->values[i].dt = 
                                &(rtData->expressions[i].dtExp->getValue());
                            rtData->expressions[i].dtExp = 0; // clear
                        }
                        numDone++;
                    }
                    break;
                case FRIfaces::scheduleType:
                    rtData->values[i].schedExp = 
                        rtData->expressions[i].schedExp;
                    if (rtData->expressions[i].schedExp->isKnown()){
                        if (evalArg){
                            rtData->values[i].sched =
                                rtData->expressions[i].schedExp->getValue();
                            rtData->expressions[i].schedExp = 0; // clear
                        }
                        numDone++;
                    }
                    break;
                case FRIfaces::tabulatedFuncType:
                    rtData->values[i].tabFuncExp = 
                        rtData->expressions[i].tabFuncExp;
                    if (rtData->expressions[i].tabFuncExp->isKnown()){
                        if (evalArg){
                            rtData->values[i].tabFunc = 
                                rtData->expressions[i].tabFuncExp->getValue();
                            rtData->expressions[i].tabFuncExp = 0; // clear
                        }
                        numDone++;
                    }
                    break;
                default:
                    throw ModelException("FRInternalFunction",
                                         "Type not supported");
                }
            }
            if (!evalArg){
                rtData->expressions[i].dExp = 0; // clear any type of pointer
            }
            if (!rtData->expressions[i].dExp){
                numPopulated++;
            }
        }
        /* finally create space for extra 'magic' parameter which is
           the current index */
        FRFunction::Val  indexVal;
        indexVal.i = ctrl->getIndex();
        rtData->values[rtData->numParams] = indexVal;
        allArgsKnown = numDone == rtData->numParams;
        rtDataNeeded = numPopulated == rtData->numParams? false: true;
        return allArgsKnown;
    }

    ~FRInternalFunction(){
        delete rtData;
    }

protected:
    /** Doesn't copy or take ownership of memory. Assumes number and
        types of parameters already checked */
    FRInternalFunction(FRController*        ctrl,
                       const FRFunction*    function):
        ctrl(ctrl), function(function),
        calcType(function->getValCalcType()),
        rtData(new RTData(function)), rtDataNeeded(true), allArgsKnown(false){}

    /** Returns RTData needed by run time classes (can be null) */
    RTData* getRTData(){
        return (rtDataNeeded? rtData: 0);
    }

public:
    virtual bool isKnown() const{
        return allArgsKnown; /* ideally an alternative object is
                                created if known */
    }

    /** gets the return type of the function eg double */
    FRIfaces::VarType getReturnType() const{
        return (FRIfaces::VarType)(function->getGrammarType());
    }

    //// returns the function's name eg "MAX"
    const char* name() const{
        return function->name();
    }

    const FRFunction* getFunction() const{
        return function;
    }

    FRController* getController() const{
        return ctrl;
    }

    //// Add a 'double' argument
    void addArg(FRIfaces::IRValueDouble* rValue){
        FRFunction::Exp exp;
        exp.dExp = rValue;
        recordTypeAndArg(FRIfaces::doubleType, false, exp);
    }

    //// Add an 'int' argument
    void addArg(FRIfaces::IRValueInt* rValue){
        FRFunction::Exp exp;
        exp.iExp = rValue;
        recordTypeAndArg(FRIfaces::intType, false, exp);
    }

    //// Add a 'bool' argument
    void addArg(FRIfaces::IRValueBool* rValue){
        FRFunction::Exp exp;
        exp.bExp = rValue;
        recordTypeAndArg(FRIfaces::boolType, false, exp);
    }

    //// Add a 'double array' argument
    void addArg(FRIfaces::IRValueDoubleArray* dExpArray){
        FRFunction::Exp exp;
        exp.dExpArray = dExpArray;
        recordTypeAndArg(FRIfaces::doubleArrayType, false, exp);
    }

    //// Add an 'int array' argument
    void addArg(FRIfaces::IRValueIntArray* iExpArray){
        FRFunction::Exp exp;
        exp.iExpArray = iExpArray;
        recordTypeAndArg(FRIfaces::intArrayType, false, exp);
    }

    //// Add an 'bool array' argument
    void addArg(FRIfaces::IRValueBoolArray* bExpArray){
        FRFunction::Exp exp;
        exp.bExpArray = bExpArray;
        recordTypeAndArg(FRIfaces::boolArrayType, false, exp);
    }

    //// Add a 'date' argument
    void addArg(FRIfaces::IRValueDate* rValue){
        FRFunction::Exp exp;
        exp.dtExp = rValue;
        recordTypeAndArg(FRIfaces::dateType, false, exp);
    }

    //// Add a 'schedule' argument
    void addArg(FRIfaces::IRValueSchedule* rValue){
        FRFunction::Exp exp;
        exp.schedExp = rValue;
        recordTypeAndArg(FRIfaces::scheduleType, false, exp);
    }

    //// Add a 'tabulated function' argument
    void addArg(FRIfaces::IRValueTabulatedFunc* rValue){
        FRFunction::Exp exp;
        exp.tabFuncExp = rValue;
        recordTypeAndArg(FRIfaces::tabulatedFuncType, false, exp);
    }

    //// Add a variable 'array' argument
    void addArg(FRIfaces::ILValueExpression* var){
        FRFunction::Exp exp;
        exp.rValues = &ctrl->getRValuesForLValueExpression(var);
        recordTypeAndArg(var->getType(), true, exp);
    }

    /** Returns FRFunction::Val* to pass to function  */
    FRFunction::Val*  getValuesBegin(){
        return rtData->values;
    }

};
    
/** class for evaluating generic functions which return a bool */
class FRInternalFunctionBool: public FRInternalFunction,
                              public virtual FRIfaces::IRValueBool{
private:
    struct MyRT{
        TGetValue*                        funcGetValue;
        FRFunction::BFunc*                func;
        FRFunction::Val*                  valuesBegin;
        RTData*                           rtData; // may be null
    public:
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(FRFunction::BFunc* func, 
                      FRFunction::Val*   valuesBegin,
                      RTData*            rtData): 
            funcGetValue(&getValue), func(func), 
            valuesBegin(valuesBegin), rtData(rtData){}

        static bool getValue(void* structure){
            MyRT* rt = (MyRT*)structure;
            if (rt->rtData){
                rt->rtData->evaluateArgs();
            }
            return rt->func(rt->valuesBegin);
        }
    };

    MyRT*  rt;
public:

    /** Creates object ready to accept parameters */
    FRInternalFunctionBool(FRController*        ctrl,
                           const FRFunction*    function):
        FRInternalFunction(ctrl, function), rt(0){}

    ~FRInternalFunctionBool(){
        delete rt;
    }

    //// do the calculation
    virtual bool getValue(){
        return MyRT::getValue(getRT());
    }

    virtual FRIfaces::IRValueBool::RT* getRT(){
        if (!rt){
            rt = new MyRT(getFunction()->getFunc().bFunc, getValuesBegin(),
                          getRTData());
        }
        return (FRIfaces::IRValueBool::RT*)rt;
    }
    
    /** Just wrapper around getValue method */
    IObjectConstSP get(){
        getRT();
        return IObjectConstSP(CBool::create(getValue()));
    }

    /** checks arguments to ensure they are of correct number and type.
        Precomputes any known values. Returns new or existing object -
        which should be used in place of this */
    static FRIfaces::IRValueBool* checkArgs(FRInternalFunction* f){
        try{
            // defer any memory allocations until we've executed code that
            // might fail
            f->verifyArgs(); // check we've got the right args etc
            // precompute values where possible
            bool known = f->preCompute(); // can we call it now?

            // now work on return type
            FRIfaces::VarType retType = f->getReturnType();
            if (retType != FRIfaces::boolType){
                throw FRParseException("Function "+string(f->name())+
                                       " returns type "+
                                       FRFunction::argTypeToString(retType)+
                                       " but type bool is required");
            }
            FRIfaces::IRValueBool* function = 
                &dynamic_cast<FRIfaces::IRValueBool&>(*f);
            if (known){
                // can evaluate now 
                try{
                    return new FR::RConstBool(function->getValue());
                } catch (exception&){
                    // perhaps this function is never evaluated so defer to
                    // run time
                }
            }
            return function;
        } catch (exception& e){
            throw ModelException(e, "FRInternalFunctionBool::create");
        }
    }
};
    
/** class for evaluating generic functions which return a int */
class FRInternalFunctionInt: public FRInternalFunction,
                             public virtual FRIfaces::IRValueInt{
private:
    struct MyRT{
        TGetValue*                       funcGetValue;
        FRFunction::IFunc*               func;
        FRFunction::Val*                 valuesBegin;
        RTData*                          rtData; // may be null
    public:
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(FRFunction::IFunc* func, 
                      FRFunction::Val*   valuesBegin,
                      RTData*            rtData): 
            funcGetValue(&getValue), func(func), 
            valuesBegin(valuesBegin), rtData(rtData){}

        static int getValue(void* structure){
            MyRT* rt = (MyRT*)structure;
            if (rt->rtData){
                rt->rtData->evaluateArgs();
            }
            return rt->func(rt->valuesBegin);
        }
    };

    MyRT*  rt;
public:

    /** Creates object ready to accept parameters */
    FRInternalFunctionInt(FRController*        ctrl,
                           const FRFunction*    function):
        FRInternalFunction(ctrl, function), rt(0){}

    ~FRInternalFunctionInt(){
        delete rt;
    }

    //// do the calculation
    virtual int getValue(){
        return MyRT::getValue(getRT());
    }

    virtual FRIfaces::IRValueInt::RT* getRT(){
        if (!rt){
            rt = new MyRT(getFunction()->getFunc().iFunc, getValuesBegin(),
                          getRTData());
        }
        return (FRIfaces::IRValueInt::RT*)rt;
    }
    
    /** Just wrapper around getValue method */
    IObjectConstSP get(){
        getRT();
        return IObjectConstSP(CInt::create(getValue()));
    }

    /** checks arguments to ensure they are of correct number and type.
        Precomputes any known values. Returns new or existing object -
        which should be used in place of this */
    static FRIfaces::IRValueInt* checkArgs(FRInternalFunction* f){
        try{
            // defer any memory allocations until we've executed code that
            // might fail
            f->verifyArgs(); // check we've got the right args etc
            // precompute values where possible
            bool known = f->preCompute(); // can we call it now?

            // now work on return type
            FRIfaces::VarType retType = f->getReturnType();
            if (retType != FRIfaces::intType){
                throw FRParseException("Function "+string(f->name())+
                                       " returns type "+
                                       FRFunction::argTypeToString(retType)+
                                       " but type int is required");
            }
            FRIfaces::IRValueInt* function = 
                &dynamic_cast<FRIfaces::IRValueInt&>(*f);
            if (known){
                // can evaluate now 
                try{
                    return new FR::RConstInt(function->getValue());
                } catch (exception&){
                    // perhaps this function is never evaluated so defer to
                    // run time
                }
            }
            return function;
        } catch (exception& e){
            throw ModelException(e, "FRInternalFunctionInt::create");
        }
    }
};

/** class for evaluating generic functions which return a double */
class FRInternalFunctionDouble: public FRInternalFunction,
                              public virtual FRIfaces::IRValueDouble{
private:
    struct MyRT{
        TGetValue*                          funcGetValue;
        FRFunction::DFunc*                  func;
        FRFunction::Val*                    valuesBegin;
        RTData*                             rtData; // may be null
    public:
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(FRFunction::DFunc* func, 
                      FRFunction::Val*   valuesBegin,
                      RTData*            rtData): 
            funcGetValue(&getValue), func(func), 
            valuesBegin(valuesBegin), rtData(rtData){}

        static double getValue(void* structure){
            MyRT* rt = (MyRT*)structure;
            if (rt->rtData){
                rt->rtData->evaluateArgs();
            }
            return rt->func(rt->valuesBegin);
        }
    };

    MyRT*  rt;
public:

    /** Creates object ready to accept parameters */
    FRInternalFunctionDouble(FRController*        ctrl,
                             const FRFunction*    function):
        FRInternalFunction(ctrl, function), rt(0){}

    ~FRInternalFunctionDouble(){
        delete rt;
    }

    //// do the calculation
    virtual double getValue(){
        return MyRT::getValue(getRT());
    }

    virtual FRIfaces::IRValueDouble::RT* getRT(){
        if (!rt){
            rt = new MyRT(getFunction()->getFunc().dFunc, getValuesBegin(),
                          getRTData());
        }
        return (FRIfaces::IRValueDouble::RT*)rt;
    }
    
    /** Just wrapper around getValue method */
    IObjectConstSP get(){
        getRT();
        return IObjectConstSP(CDouble::create(getValue()));
    }

    /** checks arguments to ensure they are of correct number and type.
        Precomputes any known values. Returns new or existing object -
        which should be used in place of this */
    static FRIfaces::IRValueDouble* checkArgs(FRInternalFunction* f){
        try{
            // defer any memory allocations until we've executed code that
            // might fail
            f->verifyArgs(); // check we've got the right args etc
            // precompute values where possible
            bool known = f->preCompute(); // can we call it now?

            // now work on return type
            FRIfaces::IRValueDouble* function; 
            FRIfaces::VarType retType = f->getReturnType();
            if (retType == FRIfaces::doubleType){
                function = &dynamic_cast<FRIfaces::IRValueDouble&>(*f);
            } else if (retType == FRIfaces::intType){
                // need to insert cast
                FRIfaces::IRValueInt* intFunc =
                    &dynamic_cast<FRIfaces::IRValueInt&>(*f);
                function = new FRIntToDouble(intFunc);
            } else {
                throw FRParseException("Function "+string(f->name())+
                                       " returns type "+
                                       FRFunction::argTypeToString(retType)+
                                       " but type double is required");
            }
            if (known){
                // can evaluate now 
                try{
                    return new FR::RConstDouble(function->getValue());
                } catch (exception&){
                    // perhaps this function is never evaluated so defer to
                    // run time
                }
            }
            return function;
        } catch (exception& e){
            throw ModelException(e, "FRInternalFunctionDouble::create");
        }
    }
};
    
/** class for evaluating generic functions which return a double */
class FRInternalFunctionDoubleArray:
    public FRInternalFunction,
    public virtual FRIfaces::IRValueDoubleArray,
    public virtual FRIfaces::IHoldsState{
private:
    //// 'derived' from FRIfaces::IRValueDoubleArray::RT
    struct MyRT{
        TGetSize*                       size; // size of the array
        TGetValue*                      getValueFunc;
        FRFunction::DArrayFunc*         func;
        FRFunction::Val*                valuesBegin;
        DoubleArray                     dbles; // to be populated
        char*                           calculated; // allow cache
        RTData*                         rtData; // may be null
    private:
        static int getSize(void* structure){
            MyRT* rt = (MyRT*)structure;
            rt->calculate();
            return (rt->dbles.size());
        }
            
        static double getValue(void* structure, int index){
            MyRT* rt = (MyRT*)structure;
            rt->calculate();
            return (rt->dbles[index]);
        }
    public:
        void calculate(){
            if (!*calculated){
                if (rtData){
                    rtData->evaluateArgs();
                }
                func(valuesBegin, dbles);
                *calculated = 1;
            }
        }
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }

        explicit MyRT(FRFunction::DArrayFunc*         func,
                      FRFunction::Val*                valuesBegin,
                      char*                           calculated,
                      RTData*                         rtData):
            size(&getSize), getValueFunc(&getValue), func(func),
            valuesBegin(valuesBegin), calculated(calculated), rtData(rtData){}
    };
    MyRT*  rt;
    char   calculated; // used for initial determination (ie outside of sim)
 
public:

    /** Creates object ready to accept parameters */
    FRInternalFunctionDoubleArray(FRController*        ctrl,
                                  const FRFunction*    function):
        FRInternalFunction(ctrl, function), rt(0), calculated(0){}

    virtual bool isKnown() const{
        return (calculated? true: false);
    }

    ~FRInternalFunctionDoubleArray(){
        delete rt;
    }

    virtual FRIfaces::IRValueDoubleArray::RT* getRT(){
        if (!rt){
            rt = new MyRT(getFunction()->getFunc().dArrayFunc, 
                          getValuesBegin(),
                          &calculated, getRTData());
        }
        return (FRIfaces::IRValueDoubleArray::RT*) rt;
    }

    /** Return the double array as an object */
    virtual IObjectConstSP get(){
        getRT();
        rt->calculate(); // calculate if not already
        return IObjectConstSP::attachToRef(&rt->dbles);
    }
    /** Pass value given by controller to RT class. Controller resets 
        value pointed to by pointer before each simulation */
    virtual void setReset(char* reset){
        if (!calculated){
            // calculation must be done inside simulation - so must reset
            // after each simulation
            getRT();
            rt->calculated = reset;
        }
    }

    /** checks arguments to ensure they are of correct number and type.
        Precomputes any known values. Returns new or existing object -
        which should be used in place of this */
    static FRIfaces::IRValueDoubleArray* checkArgs(FRInternalFunction* f){
        try{
            // defer any memory allocations until we've executed code that
            // might fail
            f->verifyArgs(); // check we've got the right args etc
            // precompute values where possible
            bool known = f->preCompute(); // can we call it now?

            // now work on return type
            FRIfaces::VarType retType = f->getReturnType();
            if (retType != FRIfaces::doubleArrayType){
                throw FRParseException("Function "+string(f->name())+
                                       " returns type "+
                                       FRFunction::argTypeToString(retType)+
                                       " but type double array is required");
            }
            FRInternalFunctionDoubleArray* function = 
                &dynamic_cast<FRInternalFunctionDoubleArray&>(*f);
            if (known){
                // can evaluate now 
                try{
                    function->getRT(); // build rt
                    function->rt->calculate();
                    function->calculated = 1;
                } catch (exception&){
                    // perhaps this function is never evaluated so defer to
                    // run time
                }
            } else {
                /* get controller to manage state for us (ie the 'calculated'
                   in the RT class) */
                function->getController()->registerStateObject(function);
            }
            return function;
        } catch (exception& e){
            throw ModelException(e, "FRInternalFunctionDouble::create");
        }
    }
};

/** class for evaluating generic functions which return a double */
class FRInternalFunctionIntArray:
    public FRInternalFunction,
    public virtual FRIfaces::IRValueIntArray,
    public virtual FRIfaces::IHoldsState{
private:
    //// 'derived' from FRIfaces::IRValueIntArray::RT
    struct MyRT{
        TGetSize*                       size; // size of the array
        TGetValue*                      getValueFunc;
        FRFunction::IArrayFunc*         func;
        FRFunction::Val*                valuesBegin;
        IntArray                     dbles; // to be populated
        char*                           calculated; // allow cache
        RTData*                         rtData; // may be null
    private:
        static int getSize(void* structure){
            MyRT* rt = (MyRT*)structure;
            rt->calculate();
            return (rt->dbles.size());
        }
            
        static int getValue(void* structure, int index){
            MyRT* rt = (MyRT*)structure;
            rt->calculate();
            return (rt->dbles[index]);
        }
    public:
        void calculate(){
            if (!*calculated){
                if (rtData){
                    rtData->evaluateArgs();
                }
                func(valuesBegin, dbles);
                *calculated = 1;
            }
        }
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }

        explicit MyRT(FRFunction::IArrayFunc*         func,
                      FRFunction::Val*                valuesBegin,
                      char*                           calculated,
                      RTData*                         rtData):
            size(&getSize), getValueFunc(&getValue), func(func),
            valuesBegin(valuesBegin), calculated(calculated), rtData(rtData){}
    };
    MyRT*  rt;
    char   calculated; // used for initial determination (ie outside of sim)
 
public:

    /** Creates object ready to accept parameters */
    FRInternalFunctionIntArray(FRController*        ctrl,
                                  const FRFunction*    function):
        FRInternalFunction(ctrl, function), rt(0), calculated(0){}

    virtual bool isKnown() const{
        return (calculated? true: false);
    }

    ~FRInternalFunctionIntArray(){
        delete rt;
    }

    virtual FRIfaces::IRValueIntArray::RT* getRT(){
        if (!rt){
            rt = new MyRT(getFunction()->getFunc().iArrayFunc, 
                          getValuesBegin(),
                          &calculated, getRTData());
        }
        return (FRIfaces::IRValueIntArray::RT*) rt;
    }

    /** Return the int array as an object */
    virtual IObjectConstSP get(){
        getRT();
        rt->calculate(); // calculate if not already
        return IObjectConstSP::attachToRef(&rt->dbles);
    }
    /** Pass value given by controller to RT class. Controller resets 
        value pointed to by pointer before each simulation */
    virtual void setReset(char* reset){
        if (!calculated){
            // calculation must be done inside simulation - so must reset
            // after each simulation
            getRT();
            rt->calculated = reset;
        }
    }

    /** checks arguments to ensure they are of correct number and type.
        Precomputes any known values. Returns new or existing object -
        which should be used in place of this */
    static FRIfaces::IRValueIntArray* checkArgs(FRInternalFunction* f){
        try{
            // defer any memory allocations until we've executed code that
            // might fail
            f->verifyArgs(); // check we've got the right args etc
            // precompute values where possible
            bool known = f->preCompute(); // can we call it now?

            // now work on return type
            FRIfaces::VarType retType = f->getReturnType();
            if (retType != FRIfaces::intArrayType){
                throw FRParseException("Function "+string(f->name())+
                                       " returns type "+
                                       FRFunction::argTypeToString(retType)+
                                       " but type int array is required");
            }
            FRInternalFunctionIntArray* function = 
                &dynamic_cast<FRInternalFunctionIntArray&>(*f);
            if (known){
                // can evaluate now 
                try{
                    function->getRT(); // build rt
                    function->rt->calculate();
                    function->calculated = 1;
                } catch (exception&){
                    // perhaps this function is never evaluated so defer to
                    // run time
                }
            } else {
                /* get controller to manage state for us (ie the 'calculated'
                   in the RT class) */
                function->getController()->registerStateObject(function);
            }
            return function;
        } catch (exception& e){
            throw ModelException(e, "FRInternalFunctionInt::create");
        }
    }
};
    

/** class for evaluating generic functions which return a bool */
class FRInternalFunctionBoolArray:
    public FRInternalFunction,
    public virtual FRIfaces::IRValueBoolArray,
    public virtual FRIfaces::IHoldsState{
private:
    //// 'derived' from FRIfaces::IRValueBoolArray::RT
    struct MyRT{
        TGetSize*                       size; // size of the array
        TGetValue*                      getValueFunc;
        FRFunction::BArrayFunc*         func;
        FRFunction::Val*                valuesBegin;
        BoolArray                       dbles; // to be populated
        char*                           calculated; // allow cache
        RTData*                         rtData; // may be null
    private:
        static int getSize(void* structure){
            MyRT* rt = (MyRT*)structure;
            rt->calculate();
            return (rt->dbles.size());
        }
            
        static bool getValue(void* structure, int index){
            MyRT* rt = (MyRT*)structure;
            rt->calculate();
            return (rt->dbles[index]);
        }
    public:
        void calculate(){
            if (!*calculated){
                if (rtData){
                    rtData->evaluateArgs();
                }
                func(valuesBegin, dbles);
                *calculated = 1;
            }
        }
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }

        explicit MyRT(FRFunction::BArrayFunc*         func,
                      FRFunction::Val*                valuesBegin,
                      char*                           calculated,
                      RTData*                         rtData):
            size(&getSize), getValueFunc(&getValue), func(func),
            valuesBegin(valuesBegin), calculated(calculated), rtData(rtData){}
    };
    MyRT*  rt;
    char   calculated; // used for initial determination (ie outside of sim)
 
public:

    /** Creates object ready to accept parameters */
    FRInternalFunctionBoolArray(FRController*        ctrl,
                                  const FRFunction*    function):
        FRInternalFunction(ctrl, function), rt(0), calculated(0){}

    virtual bool isKnown() const{
        return (calculated? true: false);
    }

    ~FRInternalFunctionBoolArray(){
        delete rt;
    }

    virtual FRIfaces::IRValueBoolArray::RT* getRT(){
        if (!rt){
            rt = new MyRT(getFunction()->getFunc().bArrayFunc, 
                          getValuesBegin(),
                          &calculated, getRTData());
        }
        return (FRIfaces::IRValueBoolArray::RT*) rt;
    }

    /** Return the bool array as an object */
    virtual IObjectConstSP get(){
        getRT();
        rt->calculate(); // calculate if not already
        return IObjectConstSP::attachToRef(&rt->dbles);
    }
    /** Pass value given by controller to RT class. Controller resets 
        value pointed to by pointer before each simulation */
    virtual void setReset(char* reset){
        if (!calculated){
            // calculation must be done inside simulation - so must reset
            // after each simulation
            getRT();
            rt->calculated = reset;
        }
    }

    /** checks arguments to ensure they are of correct number and type.
        Precomputes any known values. Returns new or existing object -
        which should be used in place of this */
    static FRIfaces::IRValueBoolArray* checkArgs(FRInternalFunction* f){
        try{
            // defer any memory allocations until we've executed code that
            // might fail
            f->verifyArgs(); // check we've got the right args etc
            // precompute values where possible
            bool known = f->preCompute(); // can we call it now?

            // now work on return type
            FRIfaces::VarType retType = f->getReturnType();
            if (retType != FRIfaces::boolArrayType){
                throw FRParseException("Function "+string(f->name())+
                                       " returns type "+
                                       FRFunction::argTypeToString(retType)+
                                       " but type bool array is required");
            }
            FRInternalFunctionBoolArray* function = 
                &dynamic_cast<FRInternalFunctionBoolArray&>(*f);
            if (known){
                // can evaluate now 
                try{
                    function->getRT(); // build rt
                    function->rt->calculate();
                    function->calculated = 1;
                } catch (exception&){
                    // perhaps this function is never evaluated so defer to
                    // run time
                }
            } else {
                /* get controller to manage state for us (ie the 'calculated'
                   in the RT class) */
                function->getController()->registerStateObject(function);
            }
            return function;
        } catch (exception& e){
            throw ModelException(e, "FRInternalFunctionBool::create");
        }
    }
};

/** class for evaluating generic functions which return a date */
class FRInternalFunctionDate: public FRInternalFunction,
                              public virtual FRIfaces::IRValueDate{
private:
    struct MyRT{
        TGetValue*               funcGetValue;
        FRFunction::DtFunc*      func;
        FRFunction::Val*         valuesBegin;
        RTData*                  rtData; // may be null
    public:
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(FRFunction::DtFunc* func, 
                      FRFunction::Val*   valuesBegin,
                      RTData*            rtData): 
            funcGetValue(&getValue), func(func), 
            valuesBegin(valuesBegin), rtData(rtData){}

        static const DateTime::Date& getValue(void* structure){
            MyRT* rt = (MyRT*)structure;
            if (rt->rtData){
                rt->rtData->evaluateArgs();
            }
            return rt->func(rt->valuesBegin);
        }
    };

    MyRT*  rt;
public:

    /** Creates object ready to accept parameters */
    FRInternalFunctionDate(FRController*        ctrl,
                           const FRFunction*    function):
        FRInternalFunction(ctrl, function), rt(0){}

    ~FRInternalFunctionDate(){
        delete rt;
    }

    //// do the calculation
    virtual const DateTime::Date& getValue(){
        return MyRT::getValue(getRT());
    }

    virtual FRIfaces::IRValueDate::RT* getRT(){
        if (!rt){
            rt = new MyRT(getFunction()->getFunc().dtFunc, getValuesBegin(),
                          getRTData());
        }
        return (FRIfaces::IRValueDate::RT*)rt;
    }
    
    /** Just wrapper around getValue method */
    IObjectConstSP get(){
        getRT();
        return IObjectConstSP(new DateTime::Date(getValue()));
    }

    /** checks arguments to ensure they are of correct number and type.
        Precomputes any known values. Returns new or existing object -
        which should be used in place of this */
    static FRIfaces::IRValueDate* checkArgs(FRInternalFunction* f){
        try{
            // defer any memory allocations until we've executed code that
            // might fail
            f->verifyArgs(); // check we've got the right args etc
            // precompute values where possible
            bool known = f->preCompute(); // can we call it now?

            // now work on return type
            FRIfaces::VarType retType = f->getReturnType();
            if (retType != FRIfaces::dateType){
                throw FRParseException("Function "+string(f->name())+
                                       " returns type "+
                                       FRFunction::argTypeToString(retType)+
                                       " but type date is required");
            }
            FRIfaces::IRValueDate* function = 
                &dynamic_cast<FRIfaces::IRValueDate&>(*f);
            if (known){
                // can evaluate now 
                try{
                    return new FR::RConstDate(function->getValue());
                } catch (exception&){
                    // perhaps this function is never evaluated so defer to
                    // run time
                }
            }
            return function;
        } catch (exception& e){
            throw ModelException(e, "FRInternalFunctionDate::create");
        }
    }
};

#endif
DRLIB_END_NAMESPACE
#endif
