//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRController.cpp
//
//   Description : Defines common concrete classes used internally
//
//   Author      : Mark A Robson
//
//   Date        : 25 March 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/FRController.hpp"
#include "edginc/FR.hpp"
#include "edginc/Hashtable.hpp"
#include "edginc/FRParseException.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Malloc.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include "edginc/EventResults.hpp"
#include "edginc/BarrierLevel.hpp"
#include ext_hash_map
#include ext_hash_set

DRLIB_BEGIN_NAMESPACE
FRController::INotifyPayoffCall::~INotifyPayoffCall(){}
FRController::IStateVarClient::~IStateVarClient(){}

/** Uses the object based get and set methods */
class FRController::AssignmentObject: public FRIfaces::IAssignment{
private:
    struct MyRT{
        FRIfaces::IAssignment::TAssign*  func;
        FRIfaces::ILValue*               lvalue;
        FRIfaces::IRValue*               rvalue;

        explicit MyRT(FRIfaces::ILValue* lvalue,
                      FRIfaces::IRValue* rvalue):
            func(&assign), lvalue(lvalue), rvalue(rvalue){}

        static void assign(void* structure){
            MyRT* rt = (MyRT*)structure;
            rt->lvalue->set(rt->rvalue->get());
        }
    };
    MyRT*                    rt;
public:
    AssignmentObject(FRIfaces::ILValue* lvalue, FRIfaces::IRValue* rvalue):
        rt(new MyRT(lvalue, rvalue)){}

    ~AssignmentObject(){
        delete rt;
    }

    //// get the run-time object to use ie cut down version of whole class
    virtual RT* getRT(){
        return (FRIfaces::IAssignment::RT*) rt;
    }

    virtual void assign(){
        MyRT::assign(rt);
    }
    virtual void setReset(char* reset){
        rt->lvalue->setReset(reset);
    }
};
/** Defines object which has ability to set an lvalue given an rvalue
    where the lvalue implements ILValueDouble and the rvalue implements
    IRValueDouble. This class should never be used in the run time loop since
    it is not particularly efficient. One of the specialised ones below should
    be used */
class FRController::AssignmentDouble: public FRIfaces::IAssignment{
private:
    struct MyRT{
        FRIfaces::IAssignment::TAssign*  func;
        FRIfaces::ILValueDouble*         lvalue;
        FRIfaces::IRValueDouble::RT*     rvalue;

        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(FRIfaces::ILValueDouble*     lvalue,
                      FRIfaces::IRValueDouble::RT* rvalue):
            func(&assign), lvalue(lvalue), rvalue(rvalue){}

        static void assign(void* structure){
            MyRT* rt = (MyRT*)structure;
            rt->lvalue->setValue(rt->rvalue->func(rt->rvalue));
        }
    };
    MyRT*    rt;
public:

    AssignmentDouble(FRIfaces::ILValueDouble* lvalue,
                     FRIfaces::IRValueDouble* rvalue):
        rt(new MyRT(lvalue, rvalue->getRT())) {}
    ~AssignmentDouble(){
        delete rt;
    }
    //// get the run-time object to use ie cut down version of whole class
    virtual RT* getRT(){
        return (FRIfaces::IAssignment::RT*) rt;
    }
    // assigns rvalue into lvalue using explicit double based methods
    virtual void assign(){
        MyRT::assign(rt);
    }
    virtual void setReset(char* reset){
        rt->lvalue->setReset(reset);
    }
};

/** Defines object which has ability to set an lvalue given an rvalue
    where the lvalue implements FR::LValueDouble and the rvalue implements
    IRValueDouble */
class FRController::AssignmentDoubleSpecial1: public FRIfaces::IAssignment{
private:
    struct MyRT{
        FRIfaces::IAssignment::TAssign* func;
        FR::LValueDouble::RT*           lvalue;
        FRIfaces::IRValueDouble::RT*    rvalue;

        explicit MyRT(FR::LValueDouble::RT*        lvalue,
                      FRIfaces::IRValueDouble::RT* rvalue):
            func(&assign), lvalue(lvalue), rvalue(rvalue){}

        static void assign(void* structure){
            MyRT* rt = (MyRT*)structure;
            FR::LValueDouble::RT::setValue(rt->lvalue,
                                           rt->rvalue->func(rt->rvalue));
        }
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
    };
    MyRT*                    rt;
    FRIfaces::ILValueDouble* lvalue;
public:
    //// uses FR::MemMgr
    AssignmentDoubleSpecial1(FR::LValueDouble*        lvalue,
                             FRIfaces::IRValueDouble* rvalue):
        rt(new MyRT((FR::LValueDouble::RT*)lvalue->getRT(),
                    rvalue->getRT())),
        lvalue(lvalue){}

    ~AssignmentDoubleSpecial1(){
        delete rt;
    }

    //// get the run-time object to use ie cut down version of whole class
    virtual RT* getRT(){
        return (FRIfaces::IAssignment::RT*) rt;
    }
    // assigns rvalue into lvalue using explicit double based methods
    virtual void assign(){
        MyRT::assign(rt);
    }
    virtual void setReset(char* reset){
        lvalue->setReset(reset);
    }
};

/** Defines object which has ability to set an lvalue given an rvalue
    where the lvalue implements FR::LValueDouble and the rvalue implements
    IRValueDouble */
class FRController::AssignmentDoubleSpecial2: public FRIfaces::IAssignment{
private:
    struct MyRT{
        FRIfaces::IAssignment::TAssign* func;
        FR::LValueDouble::RT*           lvalue;
        FR::RValueDouble::RT*           rvalue;

        explicit MyRT(FR::LValueDouble::RT*       lvalue,
                      FR::RValueDouble::RT*       rvalue):
            func(&assign), lvalue(lvalue), rvalue(rvalue){}

        static void assign(void* structure){
            MyRT* rt = (MyRT*)structure;
            FR::LValueDouble::RT::setValue(
                rt->lvalue, FR::LValueDouble::RT::getValue(rt->rvalue));
        }

        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
    };
    MyRT*                    rt;
    FRIfaces::ILValueDouble* lvalue;
public:
    AssignmentDoubleSpecial2(FR::LValueDouble* lvalue,
                             FR::LValueDouble* rvalue):
        rt(new MyRT((FR::LValueDouble::RT*)lvalue->getRT(),
                    (FR::RValueDouble::RT*)rvalue->getRT())),
        lvalue(lvalue) {}

    ~AssignmentDoubleSpecial2(){
        delete rt;
    }
    //// get the run-time object to use ie cut down version of whole class
    virtual RT* getRT(){
        return (FRIfaces::IAssignment::RT*) rt;
    }
    // assigns rvalue into lvalue using explicit double based methods
    virtual void assign(){
        MyRT::assign(rt);
    }
    virtual void setReset(char* reset){
        lvalue->setReset(reset);
    }
};

/** Defines object which has ability to set an lvalue given an rvalue
    where the lvalue implements ILValueDouble and the rvalue implements
    IRValueInt */
class FRController::AssignmentDoubleFromInt: public FRIfaces::IAssignment{
private:
    struct MyRT{
        FRIfaces::IAssignment::TAssign* func;
        FR::LValueDouble::RT*           lvalue;
        FRIfaces::IRValueInt::RT*       rvalue;

        explicit MyRT(FR::LValueDouble::RT*           lvalue,
                      FRIfaces::IRValueInt::RT*       rvalue):
            func(&assign), lvalue(lvalue), rvalue(rvalue){}

        static void assign(void* structure){
            MyRT* rt = (MyRT*)structure;
            FR::LValueDouble::RT::setValue(rt->lvalue,
                                           rt->rvalue->func(rt->rvalue));
        }
        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
    };
    MyRT*                    rt;
    FRIfaces::ILValueDouble* lvalue;
public:
    AssignmentDoubleFromInt(FR::LValueDouble*        lvalue,
                            FRIfaces::IRValueInt*    rvalue):
        rt(new MyRT((FR::LValueDouble::RT*)lvalue->getRT(),
                    rvalue->getRT())),
        lvalue(lvalue){}

    ~AssignmentDoubleFromInt(){
        delete rt;
    }

    //// get the run-time object to use ie cut down version of whole class
    virtual RT* getRT(){
        return (FRIfaces::IAssignment::RT*) rt;
    }
    // assigns rvalue into lvalue using explicit double based methods
    virtual void assign(){
        MyRT::assign(rt);
    }
    virtual void setReset(char* reset){
        lvalue->setReset(reset);
    }
};
/** Defines object which has ability to set an lvalue given an rvalue
    where the lvalue implements L (eg ILValueDoubleArray) and the
    rvalue implements R (eg IRValueDoubleArray)  */
template <class L, class R> class FRController::AssignmentArray:
    virtual public FRIfaces::IAssignment{
private:
    struct MyRT{
        FRIfaces::IAssignment::TAssign* func;
        typename L::RT*                 lvalue;
        typename R::RT*                 rvalue;

        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(typename L::RT* lvalue, typename R::RT* rvalue):
            func(&assign), lvalue(lvalue), rvalue(rvalue){}

        static void assign(void* structure){
            MyRT* rt = (MyRT*)structure;
            rt->lvalue->setFunc(rt->lvalue, -1, // set all elements
                                rt->rvalue);
        }
    };
    L*    lvalue;
    MyRT* rt;
public:

    AssignmentArray(L* lvalue, R* rvalue):
        lvalue(lvalue),
        rt(new MyRT((typename L::RT*)(lvalue->getRT()), rvalue->getRT())) {}
    ~AssignmentArray(){
        delete rt;
    }
    //// get the run-time object to use ie cut down version of whole class
    virtual RT* getRT(){
        return (FRIfaces::IAssignment::RT*) rt;
    }
    // assigns rvalue into lvalue using explicit double based methods
    virtual void assign(){
        MyRT::assign(rt);
    }
    virtual void setReset(char* reset){
        lvalue->setReset(reset);
    }
};

/** Defines object which has ability to set an lvalue given an rvalue
    where the lvalue implements ILValueInt and the rvalue implements
    IRValueInt. This class should never be used in the run time loop since
    it is not particularly efficient. One of the specialised ones below should
    be used  */
class FRController::AssignmentInt: public FRIfaces::IAssignment{
private:
    struct MyRT{
        FRIfaces::IAssignment::TAssign* func;
        FRIfaces::ILValueInt*           lvalue;
        FRIfaces::IRValueInt::RT*       rvalue;

        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(FRIfaces::ILValueInt*     lvalue,
                      FRIfaces::IRValueInt::RT* rvalue):
            func(&assign), lvalue(lvalue), rvalue(rvalue){}

        static void assign(void* structure){
            MyRT* rt = (MyRT*)structure;
            rt->lvalue->setValue(rt->rvalue->func(rt->rvalue));
        }
    };
    MyRT*                 rt;
public:
    AssignmentInt(FRIfaces::ILValueInt* lvalue,
                  FRIfaces::IRValueInt* rvalue):
        rt(new MyRT(lvalue, rvalue->getRT())) {}

    ~AssignmentInt(){
        delete rt;
    }

    //// get the run-time object to use ie cut down version of whole class
    virtual RT* getRT(){
        return (FRIfaces::IAssignment::RT*) rt;
    }
    // assigns rvalue into lvalue using explicit int based methods
    virtual void assign(){
        MyRT::assign(rt);
    }
    virtual void setReset(char* reset){
        rt->lvalue->setReset(reset);
    }
};

/** Defines object which has ability to set an lvalue given an rvalue
    where the lvalue implements ILValueInt and the rvalue implements
    IRValueInt */
class FRController::AssignmentIntSpecial1: public FRIfaces::IAssignment{
private:
    struct MyRT{
        FRIfaces::IAssignment::TAssign* func;
        FR::LValueInt::RT*              lvalue;
        FRIfaces::IRValueInt::RT*       rvalue;

        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(FR::LValueInt::RT*        lvalue,
                      FRIfaces::IRValueInt::RT* rvalue):
            func(&assign), lvalue(lvalue), rvalue(rvalue){}

        static void assign(void* structure){
            MyRT* rt = (MyRT*)structure;
            FR::LValueInt::RT::setValue(rt->lvalue,
                                        rt->rvalue->func(rt->rvalue));
        }
    };
    MyRT*           rt;
    FR::LValueInt* lvalue;

public:
    AssignmentIntSpecial1(FR::LValueInt*        lvalue,
                          FRIfaces::IRValueInt* rvalue):
        rt(new MyRT((FR::LValueInt::RT*)lvalue->getRT(), rvalue->getRT())),
        lvalue(lvalue){}

    ~AssignmentIntSpecial1(){
        delete rt;
    }

    //// get the run-time object to use ie cut down version of whole class
    virtual RT* getRT(){
        return (FRIfaces::IAssignment::RT*) rt;
    }
    // assigns rvalue into lvalue using explicit int based methods
    virtual void assign(){
        MyRT::assign(rt);
    }
    virtual void setReset(char* reset){
        lvalue->setReset(reset);
    }
};

/** Defines object which has ability to set an lvalue given an rvalue
    where the lvalue implements ILValueBool and the rvalue implements
    IRValueBool. This class should never be used in the run time loop since
    it is not particularly efficient. One of the specialised ones below should
    be used  */
class FRController::AssignmentBool: public FRIfaces::IAssignment{
private:
    struct MyRT{
        FRIfaces::IAssignment::TAssign* func;
        FRIfaces::ILValueBool*          lvalue;
        FRIfaces::IRValueBool::RT*      rvalue;

        //// uses FR::MemMgr
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
        explicit MyRT(FRIfaces::ILValueBool*     lvalue,
                      FRIfaces::IRValueBool::RT* rvalue):
            func(&assign), lvalue(lvalue), rvalue(rvalue){}

        static void assign(void* structure){
            MyRT* rt = (MyRT*)structure;
            rt->lvalue->setValue(rt->rvalue->func(rt->rvalue));
        }
    };
    MyRT*                  rt;
public:
    AssignmentBool(FRIfaces::ILValueBool* lvalue,
                   FRIfaces::IRValueBool* rvalue):
        rt(new MyRT(lvalue, rvalue->getRT())) {}

    ~AssignmentBool(){
        delete rt;
    }
    //// get the run-time object to use ie cut down version of whole class
    virtual RT* getRT(){
        return (FRIfaces::IAssignment::RT*) rt;
    }
    // assigns rvalue into lvalue using explicit bool based methods
    virtual void assign(){
        MyRT::assign(rt);
    }
    virtual void setReset(char* reset){
        rt->lvalue->setReset(reset);
    }
};

/** Defines object which has ability to set an lvalue given an rvalue
    where the lvalue implements ILValueBool and the rvalue implements
    IRValueBool */
class FRController::AssignmentBoolSpecial1: public FRIfaces::IAssignment{
private:
    struct MyRT{
        FRIfaces::IAssignment::TAssign* func;
        FR::LValueBool::RT*             lvalue;
        FRIfaces::IRValueBool::RT*      rvalue;

        explicit MyRT(FR::LValueBool::RT*        lvalue,
                      FRIfaces::IRValueBool::RT* rvalue):
            func(&assign), lvalue(lvalue), rvalue(rvalue){}

        static void assign(void* structure){
            MyRT* rt = (MyRT*)structure;
            FR::LValueBool::RT::setValue(rt->lvalue,
                                         rt->rvalue->func(rt->rvalue));
        }
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
    };
    MyRT*                 rt;
    FR::LValueBool*       lvalue;
public:
    AssignmentBoolSpecial1(FR::LValueBool*        lvalue,
                           FRIfaces::IRValueBool* rvalue):
         rt(new MyRT((FR::LValueBool::RT*)lvalue->getRT(), rvalue->getRT())),
         lvalue(lvalue){}

    ~AssignmentBoolSpecial1(){
        delete rt;
    }
    //// get the run-time object to use ie cut down version of whole class
    virtual RT* getRT(){
        return (FRIfaces::IAssignment::RT*) rt;
    }
    // assigns rvalue into lvalue using explicit bool based methods
    virtual void assign(){
        MyRT::assign(rt);
    }
    virtual void setReset(char* reset){
        lvalue->setReset(reset);
    }
};

/** Defines object which has ability to set an lvalue given an rvalue
    where the lvalue implements ILValueDouble and the rvalue implements
    IRValueDouble. This class should never be used in the run time loop since
    it is not particularly efficient. One of the specialised ones below should
    be used */
class FRController::AssignmentDate: public FRIfaces::IAssignment{
private:
    struct MyRT{
        FRIfaces::IAssignment::TAssign* func;
        FRIfaces::ILValueDate*          lvalue;
        FRIfaces::IRValueDate::RT*      rvalue;

        explicit MyRT(FRIfaces::ILValueDate*     lvalue,
                      FRIfaces::IRValueDate::RT* rvalue):
            func(&assign), lvalue(lvalue), rvalue(rvalue){}

        static void assign(void* structure){
            MyRT* rt = (MyRT*)structure;
            rt->lvalue->setValue(rt->rvalue->func(rt->rvalue));
        }
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
    };
    MyRT*                    rt;
public:
    AssignmentDate(FRIfaces::ILValueDate* lvalue,
                   FRIfaces::IRValueDate* rvalue):
        rt(new MyRT(lvalue, rvalue->getRT())) {}
    ~AssignmentDate(){
        delete rt;
    }
    //// get the run-time object to use ie cut down version of whole class
    virtual RT* getRT(){
        return (FRIfaces::IAssignment::RT*) rt;
    }
    // assigns rvalue into lvalue using explicit Date based methods
    virtual void assign(){
        MyRT::assign(rt);
    }
    virtual void setReset(char* reset){
        rt->lvalue->setReset(reset);
    }
};

/** Defines object which has ability to set an lvalue given an rvalue
    where the lvalue implements ILValueDouble and the rvalue implements
    IRValueDouble */
class FRController::AssignmentDateSpecial1: public FRIfaces::IAssignment{
private:
    struct MyRT{
        FRIfaces::IAssignment::TAssign* func;
        FR::LValueDate::RT*             lvalue;
        FRIfaces::IRValueDate::RT*      rvalue;

        explicit MyRT(FR::LValueDate::RT*        lvalue,
                      FRIfaces::IRValueDate::RT* rvalue):
            func(&assign), lvalue(lvalue), rvalue(rvalue){}

        static void assign(void* structure){
            MyRT* rt = (MyRT*)structure;
            FR::LValueDate::RT::setValue(rt->lvalue,
                                         rt->rvalue->func(rt->rvalue));
        }
        void* operator new(size_t size){
            return FR::MemMgr::alloc(size);
        }
        void operator delete(void *ptr){
            FR::MemMgr::dealloc(ptr);
        }
    };
    MyRT*                    rt;
    FRIfaces::ILValueDate*   lvalue;
public:
    AssignmentDateSpecial1(FRIfaces::ILValueDate* lvalue,
                           FRIfaces::IRValueDate* rvalue):
        rt(new MyRT((FR::LValueDate::RT*)lvalue->getRT(),
                    rvalue->getRT())), lvalue(lvalue) {}
    ~AssignmentDateSpecial1(){
        delete rt;
    }
    //// get the run-time object to use ie cut down version of whole class
    virtual RT* getRT(){
        return (FRIfaces::IAssignment::RT*) rt;
    }
    // assigns rvalue into lvalue using explicit Date based methods
    virtual void assign(){
        MyRT::assign(rt);
    }
    virtual void setReset(char* reset){
        lvalue->setReset(reset);
    }
};

FRIfaces::IAssignment* FRController::createAssignment(
    FRIfaces::ILValue* lvalue, FRIfaces::IRValue* rvalue){
    static const string method("FRController::createAssignment");
    FRIfaces::ILValueDouble* lvalueDouble;
    FRIfaces::ILValueBool*   lvalueBool;
    FRIfaces::ILValueDate*   lvalueDate;
    FRIfaces::ILValueInt*    lvalueInt;
    FRIfaces::ILValueDoubleArray* lvalueDoubleArray;
    FRIfaces::ILValueIntArray*    lvalueIntArray;
    FRIfaces::ILValueBoolArray*   lvalueBoolArray;

    if (rvalue->isKnown()){
        lvalue->set(rvalue->get());
        return 0; // no assignment needed at sim time
    }
    if ((lvalueDouble = dynamic_cast<FRIfaces::ILValueDouble*>(lvalue))){
        FRIfaces::IRValueDouble* rvalueDouble =
            dynamic_cast<FRIfaces::IRValueDouble*>(rvalue);
        FR::LValueDouble* lvalueDoubleSpecial1 =
            dynamic_cast<FR::LValueDouble*>(lvalueDouble);
        if (rvalueDouble){
            if (lvalueDoubleSpecial1){
                FR::LValueDouble* lvalueDoubleSpecial2 =
                    dynamic_cast<FR::LValueDouble*>(rvalueDouble);
                if (lvalueDoubleSpecial2){
                    return new AssignmentDoubleSpecial2(lvalueDoubleSpecial1,
                                                        lvalueDoubleSpecial2);
                }
                return new AssignmentDoubleSpecial1(lvalueDoubleSpecial1,
                                                    rvalueDouble);
            }
            return new AssignmentDouble(lvalueDouble, rvalueDouble);
        } else {
            FRIfaces::IRValueInt* rvalueInt =
                dynamic_cast<FRIfaces::IRValueInt*>(rvalue);
            if (rvalueInt && lvalueDoubleSpecial1){
                return new AssignmentDoubleFromInt(lvalueDoubleSpecial1,
                                                   rvalueInt);
            }
            throw FRParseException("Trying to assign non-double to "
                                   "a double variable");
        }
    } else if ((lvalueInt = dynamic_cast<FRIfaces::ILValueInt*>(lvalue))){
        FRIfaces::IRValueInt* rvalueInt =
            dynamic_cast<FRIfaces::IRValueInt*>(rvalue);
        if (rvalueInt){
            FR::LValueInt* lvalueIntSpecial1 =
                dynamic_cast<FR::LValueInt*>(lvalueInt);
            if (lvalueIntSpecial1){
                return new AssignmentIntSpecial1(lvalueIntSpecial1, rvalueInt);
            }
            return new AssignmentInt(lvalueInt, rvalueInt);
        } else {
            throw FRParseException("Trying to assign non-int to a int "
                                   "variable");
        }
    } else if ((lvalueBool = dynamic_cast<FRIfaces::ILValueBool*>(lvalue))){
        FRIfaces::IRValueBool* rvalueBool =
            dynamic_cast<FRIfaces::IRValueBool*>(rvalue);
        if (rvalueBool){
            FR::LValueBool* lvalueBoolSpecial1 =
                dynamic_cast<FR::LValueBool*>(lvalueBool);
            if (lvalueBoolSpecial1){
                return new AssignmentBoolSpecial1(lvalueBoolSpecial1,
                                                  rvalueBool);
            }
            return new AssignmentBool(lvalueBool, rvalueBool);
        } else {
            throw FRParseException("Trying to assign non-bool to a bool "
                                   "variable");
        }
    } else if ((lvalueDate = dynamic_cast<FRIfaces::ILValueDate*>(lvalue))){
        FRIfaces::IRValueDate* rvalueDate =
            dynamic_cast<FRIfaces::IRValueDate*>(rvalue);
        if (rvalueDate){
            return new AssignmentDate(lvalueDate, rvalueDate);
        } else {
            throw FRParseException("Trying to assign a non-date to a "
                                   "date variable");
        }
    } else if ((lvalueDoubleArray =
              dynamic_cast<FRIfaces::ILValueDoubleArray*>(lvalue))){
        FRIfaces::IRValueDoubleArray* rvalueDoubleArray =
            dynamic_cast<FRIfaces::IRValueDoubleArray*>(rvalue);
        if (!rvalueDoubleArray){
            throw FRParseException("Trying to assign non-double array to "
                                   "a double array variable");
        }
        return new AssignmentArray<FRIfaces::ILValueDoubleArray,
            FRIfaces::IRValueDoubleArray>(lvalueDoubleArray,
                                          rvalueDoubleArray);
    } else if ((lvalueIntArray =
              dynamic_cast<FRIfaces::ILValueIntArray*>(lvalue))){
        FRIfaces::IRValueIntArray* rvalueIntArray =
            dynamic_cast<FRIfaces::IRValueIntArray*>(rvalue);
        if (!rvalueIntArray){
            throw FRParseException("Trying to assign non-int array to "
                                   "a int array variable");
        }
        return new AssignmentArray<FRIfaces::ILValueIntArray,
            FRIfaces::IRValueIntArray>(lvalueIntArray,
                                          rvalueIntArray);
    } else if ((lvalueBoolArray =
              dynamic_cast<FRIfaces::ILValueBoolArray*>(lvalue))){
        FRIfaces::IRValueBoolArray* rvalueBoolArray =
            dynamic_cast<FRIfaces::IRValueBoolArray*>(rvalue);
        if (!rvalueBoolArray){
            throw FRParseException("Trying to assign non-bool array to "
                                   "a bool array variable");
        }
        return new AssignmentArray<FRIfaces::ILValueBoolArray,
            FRIfaces::IRValueBoolArray>(lvalueBoolArray,
                                          rvalueBoolArray);
    }
    return new AssignmentObject(lvalue, rvalue);
}

// hash rules for IRValueExpression
struct FRControllerHash {
    // hash method
    size_t operator()(const FRIfaces::IRValueExpression* p) const{
        return p->hashCode();
    }
    // equality method
    bool operator()(const FRIfaces::IRValueExpression* p1,
                    const FRIfaces::IRValueExpression* p2) const{
        return p1->equals(p2);
    }
};

// hash rules for IRValueSP
struct FRControllerPtrHash {
    // hash method
    size_t operator()(const FRIfaces::IRValueSP& p) const{
        return (size_t)p.get();
    }
    // equality method
    bool operator()(const FRIfaces::IRValueSP& p1,
                    const FRIfaces::IRValueSP& p2) const{
        return (p1.get() == p2.get());
    }
};

typedef hash_map<const FRIfaces::IRValueExpression*, FRIfaces::IRValue*,
    FRControllerHash, FRControllerHash> RValueHash;

// hash map of vectors. One vector for each variable. Vector holds rValue
// for each timepoint ie the object on which getValue() is invoked in order
// to populate the variable. Forced to use a union to avoid casts in simulation
// time expressions (for parse type of functions who want to look at values at
// at different indexes
typedef hash_map<const FRIfaces::ILValueExpression*,
    vector<FRIfaces::RValUnion>, FRControllerHash,
    FRControllerHash> LValuesHash;

typedef hash_map<string, FRIfaces::ILValueExpression*,
                 Hashtable::StringHash> LValueVarsHash;

// set of IRValueSP's - used for memory management by the parser - it just
// sticks any object it makes in here
typedef hash_set<FRIfaces::IRValueSP,
    FRControllerPtrHash, FRControllerPtrHash> HashSetIRValue;

typedef refCountPtr<FRController::INotifyPayoffCall> INotifyPayoffCallSP;
typedef refCountPtr<FRController::IStateVarClient> IStateVarClientSP;

class FRBarrierLevelSupport {
public:
    string   assetName;
    bool     isUp;
    DateTime date;
    double   level;
    const FRIfaces::IVarBarrierLevelAssist* barLevelAssist;

    FRBarrierLevelSupport() {} // for array

    FRBarrierLevelSupport(string          assetName,
                          bool            isUp, 
                          const DateTime& date, 
                          double          level,
                          const FRIfaces::IVarBarrierLevelAssist* barLevelAssist):
        assetName(assetName), isUp(isUp), date(date), 
        level(level), barLevelAssist(barLevelAssist) {}

    BarrierLevel createBarrierLevel(FRController* frCtrl) const {
        return BarrierLevel(isUp, 
                            date, 
                            level * barLevelAssist->getRefLevel(frCtrl), 
                            false); // never continuous
    }
};
typedef hash_map<string, vector<FRBarrierLevelSupport>,
                 Hashtable::StringHash> BarLevelSupportByAssetName;    

class FRController::Imp{
public:
    // fields - NB there is a bit of duplication here - need to review
    // exactly what we need and don't need
    const FRIfaces::IProductView*        productView; // ref
    int                                  numDates;
    DateTimeArrayConstSP                 simDates;
    int                                  numDatesInPast;
    int                                  index;
    LValueVarsHash                       variables;
    vector<RValueHash>                   rValuesByIndex;
    LValuesHash                          lValuesHash;
    vector<FRIfaces::IAssignmentArraySP> algorithm;
    vector< vector<FRIfaces::IAssignment::RT*> > algorithmRT;
    // when non state vars, length = num sim dates, otherwise the number of
    // of future pay dates
    vector<FR::LValueDouble::RT*>        payDouble;
    const MCPathGenerator*     pathGen;   // ref
    double                               notional;
    DateTimeArrayConstSP                 possiblePayDates;
    IntArray                             realPaymentIdxs; // index of dates for which there's a "real" payoff
    // remove discountFactors field when we switch to state vars
    DoubleArray                          discountFactors;
    int                                  firstPayDate; /* retire when we switch
                                                          to state vars */
    CashFlowArray                        knownCashFlows;
    bool                                 useStateVars;
    bool                                 allowResolutionOfConstVars;
    HashSetIRValue                       rValueStore; // for memory management
    vector<INotifyPayoffCallSP>          toBeNotified;
    vector<IStateVarClientSP>            stateVarClients;
    vector<vector<FRIfaces::IHoldsState*> > holdsState;
    char*                                resets;
    int                                  numResets;
    vector<char*>                        resetsByDate;
    double                               fairValueForPast;
    SVDiscFactorSP                    dfSV;  // df state variable
    SVGenDiscFactorSP                       dfGen; // df generator

    // for storing data from events
    vector<IFREventSP>                   events;
    int                                  eventIndex;
    bool                                 triggerEvents;

    // sadly specific for BARRIER_LEVEL request when possible
    BarLevelSupportByAssetName           barLevelsByAssetName;

    // Record of debug info
    vector<FR::LValueBool::RT*>          captureDebug; // similar to payDouble above
    HashtableSP                          conditionedDebugInfo;
    HashtableSP                          minDebugInfo;
    HashtableSP                          maxDebugInfo;
    HashtableSP                          avgDebugInfo;
    HashtableSP                          avgDebugCounts;
    vector<FR::LValueDouble::RT*>        distnVars; // size == number of vars for which distn requested
    vector<int>                          distnDateIdx; // dates on which distn requested
    HashtableSP                          distnDebugInfo; 

    // constructor
    Imp(const FRIfaces::IProductView*  productView,
        const DateTimeArrayConstSP&    simDates,
        bool                           allowResolutionOfConstVars,
        bool                           triggerEvents):
        productView(productView),
        numDates(simDates->size()), simDates(simDates), index(0),
        rValuesByIndex(numDates), algorithm(numDates),
        algorithmRT(numDates),
        pathGen(0), notional(productView->getNotional()),
        possiblePayDates(productView->getPayDates()),
        realPaymentIdxs(0),
        discountFactors(numDates),
        firstPayDate(productView->getDiscountFactors(possiblePayDates,
                                                     discountFactors)),
        knownCashFlows(),
        allowResolutionOfConstVars(allowResolutionOfConstVars),
        holdsState(numDates), resets(0), numResets(0), resetsByDate(numDates),
        fairValueForPast(0.0), events(0), eventIndex(-1),
        triggerEvents(triggerEvents) {
        static const string routine("FRController::FRController");
        try{
            const DateTime& today = productView->getValueDate();
            numDatesInPast = today.numPastDates(*simDates); // defines loops
            // set up hash of variables
            const FRIfaces::ILValueExpressionArray* vars =
                productView->getVariables();
            // set up our memory manager
            FR::MemMgr::intialise(Maths::max(vars->size(),1) /* good guess */,
                                  numDates);
            for (int i = 0; i < vars->size(); i++){
                FRIfaces::ILValueExpression* var = (*vars)[i].get();
                if (!var){
                    throw ModelException(routine, "Variable number "+
                                         Format::toString(i+1)+
                                         " is null (check "
                                         "variables declared in instrument");
                }
                const string& id = var->getID();
                if (variables.find(id) != variables.end()){
                    throw ModelException(routine,
                                         "Duplicate definition of variable "
                                         +id);
                }
                variables[id] = (*vars)[i].get();
                // reserve space for hash of rvalues for lvalues - and set
                // to null. getRValuesForLValueExpression() populates them
                FRIfaces::RValUnion rValUnion;
                rValUnion.bExp = 0; // any ptr will do
                lValuesHash[(*vars)[i].get()] =
                    vector<FRIfaces::RValUnion>(numDates, rValUnion);
            }
        } catch (exception& e){
            FR::MemMgr::shutdown();
            throw ModelException(e, routine);
        }
    }

    /** memory management utility - stores pointer and frees it when
        FRController goes out of scope */
    void store(FRIfaces::IRValue*  rValue){
        // have to be a bit careful with the memory
        FRIfaces::IRValueSP  rValueSP(rValue);
        // insert if and only if there is nothing there already
        pair<HashSetIRValue::iterator, bool> entry =
            rValueStore.insert(rValueSP);
        if (!entry.second){
            // entry already existed
            rValueSP.release(); // ensure we don't free it again
        }
    }


    /** memory management utility - stores pointer and frees it when
        FRController goes out of scope */
    void store(const FRIfaces::IRValueSP&  rValueSP){
        // only store if the SP actually owns the memory
        if (rValueSP.owns()){
            // our hash is really a set of pointers with a "delete all
            // pointers" on free. So we need to break the link with
            // the reference count
            FRIfaces::IRValueSP myRValueSP(rValueSP.release());
            // insert if and only if there is nothing there already -
            // this copes with the pointer having already been stored
            // via the above function
            pair<HashSetIRValue::iterator, bool> entry =
                rValueStore.insert(myRValueSP);
            if (!entry.second){
                // entry already existed
                myRValueSP.release(); // ensure we don't free it
            }
        }
    }

    void validateIndex(int index) const{
        if (index < 0 || index >= numDates){
            throw ModelException("FRController::Imp::validateIndex",
                                 "Index ("+Format::toString(index)+") is out"
                                 " of bounds");
        }
    }
    void notifyToBeNotified(const FRController* ctrl){
        for (unsigned int i = 0; i < toBeNotified.size(); i++){
            toBeNotified[i]->update(ctrl);
        }
    }
    /** reset all variables for timepoints specified by start and end */
    void resetVariables(int start, int end){
        char* resetStart = resetsByDate[start];
        char* resetEnd = end == numDates?
            resets+numResets: resetsByDate[end];
        memset(resetStart, 0, resetEnd - resetStart);
    }
    ~Imp(){
        Malloc::deallocate(resets);
    }
};

/** Returns how many simulation dates there are (ie max value of
    getIndex()) */
int FRController::numDates() const{
    return my->numDates;
}

/** returns the index of the current point */
int FRController::getIndex() const{
    return my->index;
}

/** returns the date corresponding to index */
const DateTime& FRController::getDate(int index) const{
    my->validateIndex(index);
    return (*my->simDates)[index];
}

/** Returns the path generator for the current path. Returns null outside
    of payoff */
const MCPathGenerator*  FRController::getPathGenerator() const{
    if (!my->pathGen){
        throw ModelException("FRController::getPathGenerator", "Path "
                             "Generator not available");
    }
    return my->pathGen;
}

/** Returns the IRValue object for supplied variable at
    specified index. Will return null if no value set */
FRIfaces::IRValue* FRController::getRValue(
    const FRIfaces::IRValueExpression* expression,
    int                                index) const{
    my->validateIndex(index);
    const RValueHash& rValues = my->rValuesByIndex[index];
    RValueHash::const_iterator iter = rValues.find(expression);
    if (iter == rValues.end()){
        return 0;
    }
    return iter->second;
}

/** Returns the IRValue object for supplied variable at
    the current index */
FRIfaces::IRValue* FRController::getRValue(
    const FRIfaces::IRValueExpression* expression) const{
    return getRValue(expression, my->index);
}
/** dynamically casts rValue to type based on lValueExpression and
    populates corresponding value in union */
void FRController::populateRValueUnion(
    const FRIfaces::ILValueExpression* lValue,
    FRIfaces::IRValue*                 rValue,
    FRIfaces::RValUnion&               rValUnion){
    FRIfaces::VarType type = lValue->getType();
    switch (type){
    case FRIfaces::doubleType:
        rValUnion.dExp =
            (dynamic_cast<FRIfaces::IRValueDouble&>(*rValue)).getRT();
        break;
    case FRIfaces::intType:
        rValUnion.iExp =
            (dynamic_cast<FRIfaces::IRValueInt&>(*rValue)).getRT();
        break;
   case FRIfaces::boolType:
        rValUnion.bExp =
            (dynamic_cast<FRIfaces::IRValueBool&>(*rValue)).getRT();
        break;
    case FRIfaces::dateType:
        rValUnion.dtExp =
            (dynamic_cast<FRIfaces::IRValueDate&>(*rValue)).getRT();
        break;
    case FRIfaces::scheduleType:
        rValUnion.schedExp =
            &dynamic_cast<FRIfaces::IRValueSchedule&>(*rValue);
        break;
    case FRIfaces::tabulatedFuncType:
        rValUnion.tabFuncExp =
            &dynamic_cast<FRIfaces::IRValueTabulatedFunc&>(*rValue);
        break;
    case FRIfaces::doubleArrayType:
        rValUnion.dArrayExp =
            (dynamic_cast<FRIfaces::IRValueDoubleArray&>(*rValue)).getRT();
        break;
    case FRIfaces::intArrayType:
        rValUnion.iArrayExp =
            (dynamic_cast<FRIfaces::IRValueIntArray&>(*rValue)).getRT();
        break;
    case FRIfaces::boolArrayType:
        rValUnion.bArrayExp =
            (dynamic_cast<FRIfaces::IRValueBoolArray&>(*rValue)).getRT();
        break;
    default:
        throw ModelException("FRController::populateRValueUnion",
                             "Unrecognised variable type");
    }
}

/** Sets/stores the RValueSP for the specified index for the specified
    expression */
void FRController::setRValue(int                                index,
                             const FRIfaces::IRValueExpression* expression,
                             const FRIfaces::IRValueSP&         rValue){
    static const string routine("FRController::setRValue");
    my->validateIndex(index);
    RValueHash& rValues = my->rValuesByIndex[index];
    RValueHash::const_iterator iter = rValues.find(expression);
    if (iter != rValues.end()){
        throw ModelException(routine, "internal error - "
                             "rvalue for expression has already been set");
    }
    // see if it's an lvalue
    const FRIfaces::ILValueExpression* lValue =
        dynamic_cast<const FRIfaces::ILValueExpression*>(expression);
    if (lValue){
        const string& id = lValue->getID();
        // ensure it was in original list
        if (my->variables.find(id) == my->variables.end()){
            throw ModelException(routine, "Variable "+id+ " was "
                                 "not declared in list of variables");
        }
        // store it in our lValues hash - overwrite is ok
        populateRValueUnion(lValue, rValue.get(),
                            my->lValuesHash[lValue][index]);
        FRIfaces::ILConstValue* constValue =
            dynamic_cast<FRIfaces::ILConstValue*>(rValue.get());
        if (constValue){
            constValue->setIsValueKnown(my->allowResolutionOfConstVars);
        }
    }
    rValues[expression] = rValue.get();
    // and stick in hash set to ensure it gets freed at the end and not before
    store(rValue);
}

/** Sets/stores the RValueSP for current index for the specified
    expression */
void FRController::setRValue(const FRIfaces::IRValueExpression* expression,
                         const FRIfaces::IRValueSP&             rValue) {
    setRValue(my->index, expression, rValue);
}

/** Given variable name returns LValueExpression representing that
 * variable */
FRIfaces::ILValueExpression* FRController::getLValueExpression(
    const string& varName) const{
    // this method is needed by the parser
    LValueVarsHash::const_iterator iter(my->variables.find(varName));
    if (iter == my->variables.end()){
        return 0;
    }
    return iter->second;
}

/** Returns a vector of IRValues which reflect how the given
    variable is calculated at each simulation date */
const vector<FRIfaces::RValUnion>& FRController::getRValuesForLValueExpression(
    const FRIfaces::ILValueExpression* lValueExp) {
    LValuesHash::iterator iter = my->lValuesHash.find(lValueExp);
    if (iter == my->lValuesHash.end()){
        throw ModelException("FRController::getRValuesForLValueExpression",
                             "Internal error");
    }
    // then ensure none are null
    vector<FRIfaces::RValUnion>& valsByIndex = iter->second;
    // identify type
//    for (unsigned int i = 0; i < valsByIndex.size(); i++){
    unsigned int maxIdx = Maths::min(valsByIndex.size(), 1+getIndex());
    for (unsigned int i = 0; i < maxIdx; i++){
        FRIfaces::RValUnion& rValUnion = valsByIndex[i];
        if (!rValUnion.bExp){ // any ptr will do
            // get hold of rValue for this timepoint - this will invoke a
            // setRValue call which populates the rValUnion
            lValueExp->getRValue(i, this);
            if (!rValUnion.bExp){ // any ptr will do
                throw ModelException("FRController::getRValuesFor"
                                     "LValueExpression", "Internal error");
            }
        }
    }
    return iter->second;
}

bool FRController::getTriggerEvents() const {
    return my->triggerEvents;
}

int FRController::getEventIndex() const {
    return my->eventIndex;
}

void FRController::addEvent(IFREvent* event) {
    my->events.push_back(IFREventSP(event));
}

void FRController::retrieveEvents(EventResults* events) {
    // loop through events
    for (unsigned int i = 0; i < my->events.size(); i++) {
        IFREvent* flexEvent = my->events[i].get();
        // did we have an event?
        if (flexEvent->hasEvent()) {
            // report it
            events->addEvent(flexEvent->createEvent());
       }
    }
}

void FRController::addBarrierLevel(const string&   assetName,
                                   bool            isUp,
                                   double          barLevel,
                                   const DateTime& barrierDate,
                                   const FRIfaces::IVarBarrierLevelAssist* barLevelAssist) {
    FRBarrierLevelSupport bls(assetName,
                              isUp,
                              barrierDate,
                              barLevel, 
                              barLevelAssist);
    my->barLevelsByAssetName[assetName].push_back(bls);
}

void FRController::getBarrierLevelReports(vector<string>&              assetNames,
                                          vector<BarrierLevelArraySP>& levels) {
    assetNames.clear();
    levels.clear();
    for (BarLevelSupportByAssetName::const_iterator iter = my->barLevelsByAssetName.begin();
         iter != my->barLevelsByAssetName.end(); ++iter){
        assetNames.push_back(iter->first);
        const vector<FRBarrierLevelSupport>& blSupps = iter->second;
        BarrierLevelArraySP barLevels(new BarrierLevelArray(0));
        for(unsigned int i=0; i<blSupps.size(); i++) {
            (*barLevels).push_back(blSupps[i].createBarrierLevel(this));
        }
        levels.push_back(barLevels);
    }
}

/** Creates FRController object ready for simulation run. */
FRController::FRController(
    const FRIfaces::IProductView* productView,
    bool                          useStateVars,
    bool                          allowResolutionOfConstVars,
    bool                          triggerEvents):
    my(new Imp(productView, productView->getSimDates(),
               allowResolutionOfConstVars, triggerEvents)){
    static const string routine("FRController::FRController");
    my->useStateVars = useStateVars;
    my->eventIndex = triggerEvents ? my->numDatesInPast - 1 : -1;

    // need to build up vector of AssignmentArrays (ie 'algorithm')
    string  parseError; // to build up collection of parse errors
    int numErrors = 0;
    try{
        try{
            for (my->index = 0; my->index < my->numDates; my->index++){
                bool reportDate = false;
                try{
                    FRIfaces::IAssignmentArraySP assigments(
                        productView->createAssignments(this));
                    my->algorithm[my->index] = assigments;
                    my->algorithmRT[my->index] =
                        vector<FRIfaces::IAssignment::RT*>(assigments->size());
                    vector<FRIfaces::IAssignment::RT*>& assignsRT =
                        my->algorithmRT[my->index];
                    for (unsigned int j = 0; j < assigments->size(); j++){
                        if (!((*assigments)[j])){
                            /* we allow null in the array so that we can get
                               the right error message in assignValues() */
                            assignsRT[j] = 0;
                        } else {
                            assignsRT[j] = (*assigments)[j]->getRT();
                            my->numResets++;
                        }
                    }
                    /* add on those objects explicitly registered via
                       registerStateObject */
                    my->numResets += my->holdsState[my->index].size();
                } catch (ModelException& e){
                    ModelException* cause = e.getCause();
                    FRParseException* parseEx =
                        dynamic_cast<FRParseException*>(cause? cause:&e);
                    if (parseEx){
                        // record error and continue
                        parseError += string(e.what())+"\n";
                        reportDate = true;
                        numErrors += parseEx->getNumErrors();
                        if (numErrors > FRParseException::MAX_NUM_ERRORS){
                            throw FRParseException(numErrors,
                                                   parseError+"Max number of "
                                                   "errors exceeded");
                        }
                    } else {
                        throw;
                    }
                }
                if (reportDate && !parseError.empty()){
                    parseError += "When determining assignments "
                        "for date "+(*my->simDates)[my->index].toString()+"\n";
                }
            }
        } catch (exception& e){
            throw ModelException(e, routine,
                                 "Failed when determining assignments "
                                 "for date "+(*my->simDates)[my->index].
                                 toString());
        }
        if (!parseError.empty()){
            throw ModelException(routine, parseError);
        }

        // then look up payment variable at each timepoint
        my->realPaymentIdxs.clear();
        const FRIfaces::ILValueExpression* payVar =
            productView->getPayVariable(this);
        if (payVar){
            my->payDouble.resize(my->numDates);
            for (my->index = 0; my->index < my->numDates;
                 my->index++){
                FRIfaces::IRValue* payValue = getRValue(payVar, my->index);
                if (payValue){
                    FR::LValueDouble* lValueDouble =
                        dynamic_cast<FR::LValueDouble*>(payValue);
                    if (!lValueDouble){
                        throw ModelException(routine,
                                             "Payment variable is "
                                             "not a standard double variable");
                    }
                    // if identically 0 then ignore as a payment date
                    if (!payValue->isKnown() || !Maths::isZero(lValueDouble->getValue())) {
                        my->realPaymentIdxs.push_back(my->index);
                    }
                    if (my->index >= my->firstPayDate) {
                        my->payDouble[my->index] = 
                            (FR::LValueDouble::RT*)lValueDouble->getRT();
                    }
                }
            }

            // now ask for a discount factor state variable generator for the
            // actual dates needed (rather than all) SNN - backed out while issue resolved
            /* Thought: could possibly optimise this by using determinstic rates
               whenever the value of the cashflow is known */
            if (my->useStateVars){
                my->dfGen = productView->
                    getDiscFactorGen(*my->simDates);
            }
        }

        // look up debug var.
        const FRIfaces::DebugRequest* debugReq =
            productView->getDebugRequest(this);
        if (debugReq){
            const FRIfaces::ILValueExpression* debugVar =
                debugReq->getFlagVar().get();
            if (debugVar) {
                my->captureDebug.resize(my->numDates);
                for (my->index = my->firstPayDate; my->index < my->numDates;
                     my->index++){
                    FRIfaces::IRValue* debugValue = getRValue(debugVar, my->index);
                    if (debugValue){
                        FR::LValueBool* lValueBool =
                            dynamic_cast<FR::LValueBool*>(debugValue);
                        if (!lValueBool){
                            throw ModelException(routine,
                                                 "Debug variable is "
                                                 "not a standard boolean variable");
                        }
                        my->captureDebug[my->index] =
                            (FR::LValueBool::RT*)lValueBool->getRT();
                    }
                }
            }
            // If distribution requested then identify which step indexes
            FRIfaces::ILValueExpressionArraySP dVars = debugReq->getDistnVars();
            if (dVars.get()) {
                const DateTimeArray& dDates = debugReq->getDistnDates();
                int lastPastDateIdx;
                my->distnDateIdx = DateTime::createMapping2(*my->simDates.get(),
                                                            dDates,
                                                            productView->getValueDate(),
                                                            lastPastDateIdx);
                my->distnVars.clear();
                for(int j=0; j<dVars->size(); j++) {
                    FRIfaces::IRValue* dValue = getRValue((*dVars)[j].get(), my->distnDateIdx[j]);
                    if (!dValue){
                        throw ModelException(routine,
                                             "Debug distribution is requested for variable " +
                                             (*dVars)[j]->getID() + " at date " 
                                             + Format::toString(my->distnDateIdx[j])  + 
                                             "but it is not defined there.");
                    }
                    FR::LValueDouble* lValueDouble = 
                        dynamic_cast<FR::LValueDouble*>(dValue);
                    if (!lValueDouble){
                        throw ModelException(routine,
                                             "Distribution can only be supplied for standard "
                                             "double type variables, but " + (*dVars)[j]->getID() +
                                             " is not of this type.");
                    }
                    my->distnVars.push_back((FR::LValueDouble::RT*)lValueDouble->getRT());
                }
            }
        }

        // then reserve our room for our resets
        my->resets = (char*) Malloc::allocate(my->numResets);
        memset(my->resets, 0, my->numResets); // clear
        // and set them
        int resetPos = 0;
        for (my->index = 0; my->index < my->numDates; my->index++){
            FRIfaces::IAssignmentArray& assigments =
                *(my->algorithm[my->index]);
            // record where each date starts
            my->resetsByDate[my->index] = my->resets+resetPos;
            // first resets for variables
            for (unsigned int j = 0; j < assigments.size(); j++){
                if (assigments[j].get()){
                    assigments[j]->setReset(my->resets+resetPos);
                    resetPos++;
                }
            }
            // and then resets for IHoldsState
            vector<FRIfaces::IHoldsState*>& holdsState =
                my->holdsState[my->index];
            for (unsigned int i = 0; i < holdsState.size(); i++){
                holdsState[i]->setReset(my->resets+resetPos);
                resetPos++;
            }
        }
    }  catch (exception&){
        this->~FRController(); // work around C++ flaw
        throw;
    }
}

/** add an object to the list of objects that get called when the
    payoff is invoked. If memoryManage is true the supplied object
    will be freed when this object is deleted. */
void FRController::addToNotifyPayoffCall(INotifyPayoffCall* obj,
                                         bool memoryManage){
    INotifyPayoffCallSP objSP;
    if (memoryManage){
        objSP = INotifyPayoffCallSP(obj);
    } else {
        objSP = INotifyPayoffCallSP(obj, NullDeleter());
    }
    my->toBeNotified.push_back(objSP);
}

/** Returns true if use of the state var framework is requested */
bool FRController::stateVarUsed() const{
    return my->useStateVars;
}

/** add an object to the list of flex variables that use state variables.
    If memoryManage is true the supplied object
    will be freed when this object is deleted. */
void FRController::addToStateVars(IStateVarClient* obj, bool memoryManage){
    IStateVarClientSP objSP(memoryManage? IStateVarClientSP(obj):
            IStateVarClientSP(obj, NullDeleter()));
    my->stateVarClients.push_back(objSP);
}

/** Populates the collector with the state variables required by the
    various assets. To be called by the containing IMCProduct */
void FRController::collectStateVars(IStateVariableCollectorSP svCollector) const{
    svCollector->append(my->dfGen.get());
    for (unsigned int i = 0; i < my->stateVarClients.size(); i++){
        my->stateVarClients[i]->collectStateVars(svCollector);
    }
}

/** To be called when the pah generator changes (currently before doing
    the past, and before doing the future).
    To be called by the containing IMCProduct */
void FRController::pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
    // this method gets called for both old and new
    if (my->useStateVars){
        my->dfSV = my->dfGen->getSVDiscFactor(newPathGen);
        for (unsigned int i = 0; i < my->stateVarClients.size(); i++){
            my->stateVarClients[i]->pathGenUpdated(newPathGen);
        }
    }
}

FRController::~FRController() {
    delete my; // must do this before shutdown
    FR::MemMgr::shutdown();
}

/** returns view of product */
const FRIfaces::IProductView* FRController::productView() const{
    return my->productView;
}

/** calculates variable for timepoints between start and end */
void FRController::assignValues(int start, int end){
    unsigned int j = 0;
    try{
        for (my->index = start; my->index < end; my->index++){
             vector<FRIfaces::IAssignment::RT*>& assignments =
                 my->algorithmRT[my->index];
            // loop over rules at each simulation date
             for (j = 0; j < assignments.size(); j++){
                // assign values to variables
                 FRIfaces::IAssignment::RT* assign = assignments[j];
                 if (assign){
                     assign->func(assign);
                 }
             }
        }
    } catch (exception& e){
        string exp(my->productView->getExpression(my->index,
                                                  (*my->simDates)[my->index],
                                                  j));
        throw ModelException(e, "FRController::assignValues", "When "
                             "evaluating "+exp+" on "+
                             (*my->simDates)[my->index].toString());
    }
}

/** tester. Ignores value of productView->getPayVariable() and uses
    supplied value. Can cope with any types of variable for 'pay'
    variable. Returns 2D array of values of 'pay' variable for each date for
    each pay variable. Runs simulation across all dates */
ObjectArraySP FRController::tester(
    const FRIfaces::ILValueExpressionArray* payVars){
    my->pathGen = 0; // switch off
    // loop over all simulation dates
    // clear variables
    my->resetVariables(0, my->numDates);
    // then calculate
    assignValues(0, my->numDates);
    // and then report the values
    ObjectArraySP allResults(new ObjectArray(payVars->size()));
    for (int var = 0; var < payVars->size(); var++){
        ObjectArraySP results(new ObjectArray(my->numDates));
        (*allResults)[var] = results;
        ObjectArray&  values = *results;
        for (int i = 0; i < my->numDates; i++){
            FRIfaces::IRValue* payValue = getRValue((*payVars)[var].get(), i);
            if (payValue){
                values[i] = IObjectSP(payValue->get()->clone());
            }
        }
    }
    return allResults;
}


/** evaluates one path */
void FRController::payoff(const MCPathGenerator*  pathGen,
                          IMCPrices&      prices){
    my->pathGen = pathGen; // record path gen
    my->notifyToBeNotified(this); /* inform anyone who want to know that we're
                                     are doing a new path */
    bool doingPast = pathGen->doingPast();

    // loop over simulation dates
#if 0
    int start = pathGen->begin(0); // all assets use same date
    int end = pathGen->end(0);
#else
    int start = doingPast? 0: my->numDatesInPast;
    int end = doingPast? my->numDatesInPast: my->numDates;
#endif
    if (!doingPast){
        // need to clear variables
        my->resetVariables(start, end);
    }
    // then calculate values
    assignValues(start, end);
    // For KNOWN_CASHFLOWS need purely past
    if (doingPast) {
        my->fairValueForPast = 0.0; // clear
        // Build list of known cash flows. We can strictly treat as "known"
        // only those that have notification date <= today. Of course some
        // payments may be known long in advance but we can't tell
        // (or not easily).
        // The treatment here drops the value
        // of a payment on the pay date, not the notification date.
        const FRIfaces::ILValueExpression* payVar =
            my->productView->getPayVariable(this);
        if (payVar){
            int rpIdx = 0; // index into realPaymentIdxs
            for (my->index = start; my->index < end && rpIdx<my->realPaymentIdxs.size(); my->index++){
                // only consider those with "real payments"
                if (my->index < my->realPaymentIdxs[rpIdx]) continue;
                rpIdx++; 
                
                FRIfaces::IRValue* payValue = getRValue(payVar, my->index);
                if (payValue){
                    FR::LValueDouble* lValueDouble =
                        dynamic_cast<FR::LValueDouble*>(payValue);
                    double payAmount =
                        FR::LValueDouble::RT::getValue(
                            (FR::LValueDouble::RT*)lValueDouble->getRT());
                    CashFlow cf((*(my->possiblePayDates.get()))[my->index],
                                payAmount * my->notional);
                    my->knownCashFlows.push_back(cf);
                }
            }
        }
    }
    // finally do pv of payment amounts
    double value;
    // need to always store a value (even if doing past)
    if (my->useStateVars){
#if 0
        // Would prefer this, but currently issue with mid-life coupons
        // and non-0 settlement
        value = my->fairValueForPast;
        const SVPath& path = my->dfSV->path();
        for (int i = path.begin(); i < path.end(); i++){
            // get payoff value
            double payment = FR::LValueDouble::RT::getValue(my->payDouble[i]);
            value += payment * path[i];
        }
        if (doingPast){
            my->fairValueForPast = value; // save for future
        }
#else
        value = 0.0; // always compute all cashflows
        const SVPath& path = my->dfSV->path();
        for (int i = my->firstPayDate; i < end; i++){
            // get payoff value
            if (my->payDouble[i]){
                double payment = 
                    FR::LValueDouble::RT::getValue(my->payDouble[i]);
                value += payment * path[i];
            }
        }
#endif
    } else {
        value = 0.0; // always compute all cashflows
        for (my->index = my->firstPayDate; my->index < end; my->index++){
            // get payoff value
            if (my->payDouble[my->index]){
                double payment =
                    FR::LValueDouble::RT::getValue(my->payDouble[my->index]);
                value += payment * my->discountFactors[my->index];
            }
        }
    }
    prices.add(value * my->notional); // and record
    // leave pathGen alone as needed in getCurrentValues()
}

/** evaluates one path recording debug info
    Note also that the dfSV is based on all sim dates and not just
    ones for which the payoff var is defined. This means we need to
    check still for a payDouble.
    Also we can be a little less paranoid about performance. */
void FRController::payoffWithDebug(const MCPathGenerator*  pathGen,
                                   IMCPrices&              prices){
    // do the work
    payoff(pathGen, prices);

    // then consider recording it
    const FRIfaces::DebugRequest* debugReq =
        my->productView->getDebugRequest(this);
    if (debugReq){
        // debug recording stops once we have met the condition
        if (!my->conditionedDebugInfo.get()) {
            // derived - has to be done every iter
            // (up to condition - if any - being met)
            updateDebugValues(debugReq->doMinValues(),
                              debugReq->doMaxValues(),
                              debugReq->doAvgValues());

            const FRIfaces::ILValueExpression* debugVar =
                debugReq->getFlagVar().get();
            if (debugVar) {
                // conditioned on a boolean - done only once
                int end = pathGen->doingPast()? my->numDatesInPast: my->numDates;
                for (my->index = my->firstPayDate; my->index < end; my->index++){
                    if (my->captureDebug[my->index]) {
                        bool isDebug = FR::LValueBool::RT::getValue(my->captureDebug[my->index]);
                        if (isDebug) {
                            my->conditionedDebugInfo = HashtableSP(getCurrentValues());
                            // once only
                            break;
                        }
                    }
                }
            }
        }
        if (debugReq->doDistns()) {
            // Distribution records every iter
            FRIfaces::ILValueExpressionArraySP distnLVEs = debugReq->getDistnVars();
            if (!my->distnDebugInfo.get()) {
                // starting this hash
                my->distnDebugInfo = HashtableSP(new Hashtable());
            }
            for(unsigned int i = 0; i<my->distnVars.size(); i++) {
                DoubleArraySP distnSoFar;
                const string& name = (*distnLVEs)[i]->getID();
                if (my->distnDebugInfo->containsKey(name)) {
                    distnSoFar = DoubleArraySP::dynamicCast(my->distnDebugInfo->get(name));
                } else {
                    distnSoFar = DoubleArraySP(new DoubleArray());
                }
                distnSoFar->push_back(FR::LValueDouble::RT::getValue(my->distnVars[i]));
                my->distnDebugInfo->put((*distnLVEs)[i]->getID(),
                                        distnSoFar);
            }
        }
    }
}

/** returns a hashtable holding the value of each variable in the latest
    simulated path */
Hashtable* FRController::getCurrentValues() const{
    if (!my->pathGen){
        throw ModelException("FRController::getCurrentValues",
                             "Path Generator not set");
    }
    // store final values of variables in debug packet
    HashtableSP  finalVals(new Hashtable());
    // set up common strings
    IObjectSP notDefinedObj(CString::create("Not defined"));
    IObjectSP sameAsBefore(CString::create("Unchanged"));
    for (LValueVarsHash::const_iterator iter = my->variables.begin();
         iter != my->variables.end(); ++iter){
        ObjectArraySP valByDate(new ObjectArray(my->numDates));
        // then loop over the simulation dates

        /* lastChanged variable used to decide sameAsBefore.  Note:
           use SP to avoid object going out of scope and then 'new'
           reusing same memory resulting in obj.get() == lastChanged */
        IObjectConstSP lastChanged;

        for (int i = 0; i < my->numDates; i++){
            // locate the corresponding value
            const RValueHash& rValues = my->rValuesByIndex[i];
            RValueHash::const_iterator iterVal = rValues.find(iter->second);
            if (iterVal == rValues.end()){
                lastChanged.reset();
                (*valByDate)[i] = notDefinedObj;
            } else {
                IObjectConstSP obj;
                bool isKnown = iterVal->second->isKnown();
                if (!isKnown){
                    // not sure if this really means it's not been
                    // calculated or not - depends on implementation of
                    // isKnown
                    try {
                        // see if we can calculate it
                        obj = iterVal->second->get();
                    } catch (exception&){
                        // store not defined
                        lastChanged.reset();
                        (*valByDate)[i] = notDefinedObj;
                    }
                } else {
                    obj = iterVal->second->get();
                }
                if (obj.get()){
                    if (i > 0 && obj.get() == lastChanged.get()){
                        // avoid massive files eg with schedules
                        (*valByDate)[i] = sameAsBefore;
                    } else {
                        lastChanged = obj;
                        (*valByDate)[i] = IObjectSP(obj->clone());
                    }
                }
            }
        }
        // insert into hash
        finalVals->put(iter->first, valByDate);
    }
    return finalVals.release();
}

// There are doubtless more elegant solutions to the min/max/avg
void FRController::updateDebugValues(bool doMin,
                                     bool doMax,
                                     bool doAvg) {
    static const string routine("FRController::updateDebugValues");
    if (!doMin && !doMax && !doAvg) {
        return;
    }
    if (!my->pathGen){
        throw ModelException(routine,
                             "Path Generator not set");
    }

    for (LValueVarsHash::const_iterator iter = my->variables.begin();
         iter != my->variables.end(); ++iter){
        // only handle scalar doubles
        FRIfaces::ILValueExpression* lve = getLValueExpression(iter->first);
        if (!lve || lve->getType()!=FRIfaces::doubleType) {
            continue;
        }

        // vals we (may) already know about
        IObjectSP notDefinedObj(CString::create("Not defined"));
        ObjectArraySP minValsSoFar;
        bool minNewHash = false;
        if (doMin) {
            if (!my->minDebugInfo.get()) {
                // starting this hash
                my->minDebugInfo = HashtableSP(new Hashtable());
            }
            if (my->minDebugInfo->containsKey(iter->first)) {
                // take existing one and update
                minValsSoFar = ObjectArraySP::dynamicCast(my->minDebugInfo->get(iter->first));
            } else {
                // starting entry for new variable
                minValsSoFar = ObjectArraySP(new ObjectArray(my->numDates, notDefinedObj));
                minNewHash = true;
            }
        }
        ObjectArraySP maxValsSoFar;
        bool maxNewHash = false;
        if (doMax) {
            if (!my->maxDebugInfo.get()) {
                // starting this hash
                my->maxDebugInfo = HashtableSP(new Hashtable());
            }
            if (my->maxDebugInfo->containsKey(iter->first)) {
                // take existing one and update
                maxValsSoFar = ObjectArraySP::dynamicCast(my->maxDebugInfo->get(iter->first));
            } else {
                // starting entry for new variable
                maxValsSoFar = ObjectArraySP(new ObjectArray(my->numDates, notDefinedObj));
                maxNewHash = true;
            }
        }
        ObjectArraySP avgValsSoFar;
        IntArraySP avgCountsSoFar;
        bool avgNewHash = false;
        if (doAvg) {
            if (!my->avgDebugInfo.get()) {
                // starting this hash
                my->avgDebugInfo = HashtableSP(new Hashtable());
                my->avgDebugCounts = HashtableSP(new Hashtable());
            }
            if (my->avgDebugInfo->containsKey(iter->first)) {
                // take existing one and update
                avgValsSoFar = ObjectArraySP::dynamicCast(my->avgDebugInfo->get(iter->first));
                avgCountsSoFar = IntArraySP::dynamicCast(my->avgDebugCounts->get(iter->first));
            } else {
                // starting entry for new variable
                avgValsSoFar = ObjectArraySP(new ObjectArray(my->numDates, notDefinedObj));
                avgCountsSoFar = IntArraySP(new IntArray(my->numDates, 0));
                avgNewHash = true;
            }
        }

        // then loop over the simulation dates
        int start = my->pathGen->doingPast()? 0: my->numDatesInPast;
        int end = my->pathGen->doingPast()? my->numDatesInPast: my->numDates;
        for (int i = start; i < end; i++){
            // locate the corresponding value
            const RValueHash& rValues = my->rValuesByIndex[i];
            RValueHash::const_iterator iterVal = rValues.find(iter->second);
            bool isDefined = false;
            double val = 0.;
            if (iterVal == rValues.end()){
                // not defined at this step - no contribution to derived debug
            } else {
                // only worth proceeding if we're dealing with a scalar double
                FRIfaces::IRValueDouble* rval = dynamic_cast<FRIfaces::IRValueDouble*>(iterVal->second);
                if (!rval) {
                    // this should have been caught above
                    throw ModelException(routine,
                                         "Non-double type found - internal error!");
                } else {
                    bool isKnown = rval->isKnown();
                    if (!isKnown){
                        // not sure if this really means it's not been
                        // calculated or not - depends on implementation of
                        // isKnown
                        try {
                            // see if we can calculate it
                            val = rval->getValue();
                            isDefined = true;
                        } catch (exception&){
                            // i.e. not defined
                        }
                    } else {
                        val = rval->getValue();
                        isDefined = true;
                    }
                    if (isDefined){
                        // we can use 'val'
                        if (doMin) {
                            CDouble* dbl = dynamic_cast<CDouble*>((*minValsSoFar)[i].get());
                            if (dbl) {
                                // already defined so update
                                double minSoFar = dbl->doubleValue();
                                if (val < minSoFar) {
                                    // rather annoying - I cannot set a value inside 'dbl' so
                                    // need to recreate
                                    (*minValsSoFar)[i] = IObjectSP(CDouble::create(val));
                                }
                            } else {
                                // will be a string ("Not defined") so initialise here
                                (*minValsSoFar)[i] = IObjectSP(CDouble::create(val));
                            }
                        }
                        if (doMax) {
                            CDouble* dbl = dynamic_cast<CDouble*>((*maxValsSoFar)[i].get());
                            if (dbl) {
                                // already defined so update
                                double maxSoFar = dbl->doubleValue();
                                if (val > maxSoFar) {
                                    // rather annoying - I cannot set a value inside 'dbl' so
                                    // need to recreate
                                    (*maxValsSoFar)[i] = IObjectSP(CDouble::create(val));
                                }
                            } else {
                                // will be a string ("Not defined") so initialise here
                                (*maxValsSoFar)[i] = IObjectSP(CDouble::create(val));
                            }
                        }
                        if (doAvg) {
                            CDouble* dbl = dynamic_cast<CDouble*>((*avgValsSoFar)[i].get());
                            double avgSoFar;
                            if (dbl) {
                                int numSoFar = (*avgCountsSoFar)[i];
                                avgSoFar = (dbl->doubleValue()*numSoFar + val)/(double)(numSoFar+1);
                            } else {
                                avgSoFar = val;
                            }
                            (*avgCountsSoFar)[i] += 1; // init to 0 at creation
                            (*avgValsSoFar)[i] = IObjectSP(CDouble::create(avgSoFar));
                        }
                    }
                }
            }
        }
        // update hash
        if (minNewHash) {
            my->minDebugInfo->put(iter->first, minValsSoFar);
        }
        if (maxNewHash) {
            my->maxDebugInfo->put(iter->first, maxValsSoFar);
        }
        if (avgNewHash) {
            my->avgDebugInfo->put(iter->first, avgValsSoFar);
            my->avgDebugCounts->put(iter->first, avgCountsSoFar);
        }
    }
}

Hashtable* FRController::getDebugValues(const string& id) const{
    if (CString::equalsIgnoreCase(id, "Cond")) {
        return my->conditionedDebugInfo.get();
    }
    if (CString::equalsIgnoreCase(id, "Min")) {
        return my->minDebugInfo.get();
    }
    if (CString::equalsIgnoreCase(id, "Max")) {
        return my->maxDebugInfo.get();
    }
    if (CString::equalsIgnoreCase(id, "Avg")) {
        return my->avgDebugInfo.get();
    }
    if (CString::equalsIgnoreCase(id, "Distn")) {
        // A slight complication here - we translate the raw numbers
        // into my->debugRequest->numBuckets
        // 1. For each distn var
        // 2. Get the raw array
        // 3. Identify max and min across all assets
        // 4. Bucket boundaries = min:max divided into numBuckets
        //    - boundaries imply one extra/fewer value of freqs. How resolve?
        //      Freq against a given bucket boundary means count of samples
        //      above adjacent lower boundary, and less than that boundary.
        // 5. Count freq of samples that fall in each bucket
        // 6. Normalise so the freq is a % (sum to 1)
        // 7. Add buckets to output hash table
        // 8. Add each var freq distn to output hash table
        const FRIfaces::DebugRequest* debugReq = 
            my->productView->getDebugRequest(this);
        if (debugReq && debugReq->doDistns()) {
            int numBuckets = debugReq->getDistnSize();
            DoubleArraySP buckets(new DoubleArray(numBuckets, 0.0));
            // Identify limits for the buckets (aka "bins")
            bool   minStarted = false;
            double minBucket;
            bool   maxStarted = false;
            double maxBucket;
            Hashtable::KeyEnumerationSP keys(my->distnDebugInfo->keys());
            while(keys->hasMoreElements()) {
                const string& name = keys->nextElement();
                DoubleArraySP rawDistn = DoubleArraySP::dynamicCast(my->distnDebugInfo->get(name));
                for(int j=0; j<rawDistn->size(); j++) {
                    if (minStarted) {
                        minBucket = Maths::min(minBucket, (*rawDistn)[j]);
                    } else {
                        minBucket = (*rawDistn)[j];
                        minStarted = true;
                    }
                    if (maxStarted) {
                        maxBucket = Maths::max(maxBucket, (*rawDistn)[j]);
                    } else {
                        maxBucket = (*rawDistn)[j];
                        maxStarted = true;
                    }
                }
            }
            double bucketSize = (maxBucket-minBucket)/numBuckets;
            for(int j=0; j<numBuckets; j++) {
                (*buckets)[j] = minBucket + (j+1) * bucketSize;
            }
            // Now process the raw results into a histogram based on these buckets
            HashtableSP distnResults(new Hashtable());
            // store our buckets
            distnResults->put("HistogramBuckets", buckets);
            Hashtable::KeyEnumerationSP kkeys(my->distnDebugInfo->keys()); // how reset 'keys'?
            while(kkeys->hasMoreElements()) {
                const string& name = kkeys->nextElement();
                DoubleArraySP rawDistn = DoubleArraySP::dynamicCast(my->distnDebugInfo->get(name));
                sort(rawDistn->begin(), rawDistn->end());
                DoubleArraySP bucketDistn(new DoubleArray(numBuckets, 0.0));
                for(int i = 0, j=0; i < numBuckets; i++) {
                    for(; j<rawDistn->size(); j++) {
                        if (Maths::isPositive((*rawDistn)[j] - (*buckets)[i])) {
                            break;
                        }
                        // includes the bucket boundary
                        (*bucketDistn)[i] += 1.0;
                    }
                    // normalise
                    (*bucketDistn)[i] /= rawDistn->size();
                }
                distnResults->put(name, bucketDistn);
            }
            return distnResults.release();
        }
    }
    return 0;
}

/** support for these output requests */
// XXX Would like to have this const but getPayVariable messes it up.
void FRController::recordPaymentDates(CControl*     control,
                                      Results*      results) {
    // It's not quite "my->possiblePayDates" since they are the payment dates
    // corresponding to all the sim dates. To be sensible we need to
    // filter according to which dates have rules which include an
    // actual payoff variable. As with the known cash flows we cannot
    // use "my->payDouble" since that is defined for future only (which
    // is fair enough for pricing)
    // The work is now done elsewhere.
    DateTimeArray realPayDates;
    for(int i=0; i<my->realPaymentIdxs.size(); i++) {
        realPayDates.push_back((*(my->possiblePayDates.get()))[my->realPaymentIdxs[i]]);
    }
    OutputRequestUtil::recordPaymentDates(control,
                                          results,
                                          &realPayDates);
}
void FRController::recordKnownCashFlows(CControl*     control,
                                        Results*      results,
                                        const string& ccyName) const {
    OutputRequestUtil::recordKnownCashflows(control,
                                            results,
                                            ccyName,
                                            &my->knownCashFlows);
}

/** memory management utility - stores pointer and frees it when
    FRController goes out of scope */
void FRController::store(FRIfaces::IRValue*  rValue){
    my->store(rValue);
}

/** memory management utility - stores pointer and frees it when
    FRController goes out of scope */
void FRController::store(const FRIfaces::IRValueSP&  rValueSP){
     my->store(rValueSP);
}

/** state management utility for objects which store state which needs to
    be reset between each simulation */
void FRController::registerStateObject(FRIfaces::IHoldsState* obj){
    my->holdsState[my->index].push_back(obj);
}

DRLIB_END_NAMESPACE


