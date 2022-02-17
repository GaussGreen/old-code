//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RootFinder.cpp
//
//   Description : 
//
//   Date        : 23 Nov 2001
//
//
//   $Log: RootFinder.cpp,v $
//   Revision 1.7  2005/04/28 11:40:59  mrobson
//   Added missing template<>
//
//   Revision 1.6  2004/02/18 17:41:40  mrobson
//   Use new macros for interfaces
//
//   Revision 1.5  2002/08/27 15:28:30  mrobson
//    RtSafe_solve can now return status/
//
//   Revision 1.4  2002/01/23 15:11:05  asegger
//   added ZBrent::validatePop2Object()
//
//   Revision 1.3  2001/12/20 10:03:01  rguichar
//   Added ZBrac + removed base classes to the RootFinder from the public interface.
//
//   Revision 1.2  2001/11/30 15:03:38  rguichar
//   Trying to get it to compile. Messaged templated function.
//
//   Revision 1.1  2001/11/28 11:56:38  rguichar
//   First revision.
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_ROOTFINDER_CPP
#include "edginc/RootFinder.hpp"

DRLIB_BEGIN_NAMESPACE

/* Needs explicit instantiation */
template<> CClassConstSP const RootFinder1D::TwoInitVal<Func1D::WtDeriv>::TYPE = CClass::registerInterfaceLoadMethod(
    "RootFinder1D::TwoInitValWtDeriv", typeid(RootFinder1D::TwoInitVal<Func1D::WtDeriv>), RootFinder1D::TwoInitVal<Func1D::WtDeriv>::load);
static CClassConstSP vc7BugFix1=RootFinder1D::TwoInitVal<Func1D::WtDeriv>::TYPE;

template <> void RootFinder1D::TwoInitVal<Func1D::WtDeriv>::load(CClassSP& clazz){
    REGISTER_INTERFACE(RootFinder1D::TwoInitVal<Func1D::WtDeriv>, clazz);
    EXTENDS(IObject);
}

template<> CClassConstSP const RootFinder1D::TwoInitVal<Func1D::NoDeriv>::TYPE = CClass::registerInterfaceLoadMethod(
    "RootFinder1D::TwoInitValNoDeriv", typeid(RootFinder1D::TwoInitVal<Func1D::NoDeriv>), RootFinder1D::TwoInitVal<Func1D::NoDeriv>::load);
static CClassConstSP vc7BugFix2=RootFinder1D::TwoInitVal<Func1D::NoDeriv>::TYPE;

template <> void RootFinder1D::TwoInitVal<Func1D::NoDeriv>::load(CClassSP& clazz){
    REGISTER_INTERFACE(RootFinder1D::TwoInitVal<Func1D::NoDeriv>, clazz);
    EXTENDS(IObject);
}

void Bracketer1D::load(CClassSP& clazz){
    REGISTER_INTERFACE(Bracketer1D, clazz);
    EXTENDS(IObject);
}    

CClassConstSP const Bracketer1D::TYPE = CClass::registerInterfaceLoadMethod(
    "Bracketer1D", typeid(Bracketer1D), Bracketer1D::load);

const double ZBrac::default_FACTOR = 1.6;
const int    ZBrac::default_NUM_TRIES = 50;

void ZBrac::bracket(const Func1D::NoDeriv& func,
                    double&                x1,
                    double&                x2) const{
    ZBrac_bracket(func,
                  x1,
                  x2,
                  FACTOR,
                  NUM_TRIES);
}

ZBrac::ZBrac(double _FACTOR,
                         int    _NUM_TRIES):
CObject(TYPE),
FACTOR(_FACTOR),
NUM_TRIES(_NUM_TRIES){}


class ZBracHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ZBrac, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(Bracketer1D);
        EMPTY_SHELL_METHOD(defaultZBrac);
        FIELD(FACTOR, "FACTOR");
        FIELD_MAKE_OPTIONAL(FACTOR);
        FIELD(NUM_TRIES, "NUM_TRIES");
        FIELD_MAKE_OPTIONAL(NUM_TRIES);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultZBrac(){
        return new ZBrac();
    }
};

CClassConstSP const ZBrac::TYPE = CClass::registerClassLoadMethod(
    "ZBrac", typeid(ZBrac), ZBracHelper::load);

const double ZBracPositive::default_FACTOR = 1.6;
const int    ZBracPositive::default_NUM_TRIES = 50;

void ZBracPositive::bracket(const Func1D::NoDeriv& func,
                            double&                x1,	            // not used on input
                            double&                x2) const{
    ZBracPositive_bracket(func,
                          x1,
                          x2,
                          true,
                          FACTOR,
                          NUM_TRIES);
}

ZBracPositive::ZBracPositive(double FACTOR,
                             int    NUM_TRIES):
CObject(TYPE),
FACTOR(default_FACTOR),
NUM_TRIES(default_NUM_TRIES){}

class ZBracPositiveHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ZBracPositive, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(Bracketer1D);
        EMPTY_SHELL_METHOD(defaultZBracPositive);
        FIELD(FACTOR, "FACTOR");
        FIELD_MAKE_OPTIONAL(FACTOR);
        FIELD(NUM_TRIES, "NUM_TRIES");
        FIELD_MAKE_OPTIONAL(NUM_TRIES);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultZBracPositive(){
        return new ZBracPositive();
    }
};

CClassConstSP const ZBracPositive::TYPE = CClass::registerClassLoadMethod(
    "ZBracPositive", typeid(ZBracPositive), ZBracPositiveHelper::load);

const int RtSafe::default_MAXIT  = 100;  // Maximum allowed number of iterations.

double RtSafe::solve(const Func1D::WtDeriv&  func,
                     double                  x1,
                     double                  x2) const{
    double result;
    RtSafe_solve(func,
                 x1,
                 x2,
                 xacc,
                 true,
                 result,
                 MAXIT);
    return result;
}

RtSafe::RtSafe(const CClassConstSP& clazz):
CObject(clazz),
MAXIT(default_MAXIT){}

RtSafe::RtSafe():
CObject(TYPE),
MAXIT(default_MAXIT){}

RtSafe::RtSafe(double xacc,
               int    MAXIT):
CObject(TYPE),
xacc(xacc),
MAXIT(MAXIT){}

class RtSafeHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(RtSafe, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(RootFinder1D::TwoInitValWtDeriv);
        EMPTY_SHELL_METHOD(defaultRtSafe);
        FIELD(xacc, "xacc");
        FIELD(MAXIT, "MAXIT");
        FIELD_MAKE_OPTIONAL(MAXIT);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultRtSafe(){
        return new RtSafe();
    }
};

CClassConstSP const RtSafe::TYPE = CClass::registerClassLoadMethod(
    "RtSafe", typeid(RtSafe), RtSafeHelper::load);

const int ZBrent::default_ITMAX  = 100;      // Maximum allowed number of iterations.
const double ZBrent::EPS  = DBL_EPSILON;

ZBrent::ZBrent(const CClassConstSP& clazz):
CObject(clazz),
ITMAX(default_ITMAX){}

ZBrent::ZBrent():
CObject(TYPE),
ITMAX(default_ITMAX){}

ZBrent::ZBrent(double tol,
               int    ITMAX):
CObject(TYPE),
tol(tol),
ITMAX(ITMAX){}

void ZBrent::validatePop2Object() {
    static const string method("ZBrent::validatePop2Object");
    
    if (!Maths::isPositive(tol)){
        throw ModelException(method,
                             "The tolerance must be strictly greater than 0.0");
    }

    if (!Maths::isPositive(ITMAX)){
        throw ModelException(method,
                             "The maximum number of iterations must be strictly greater than 0");
    }
}

double ZBrent::solve(const Func1D::NoDeriv&  func,
                     double                  x1,
                     double                  x2) const{
    return ZBrent_solve(func,
                        x1,
                        x2,
                        tol,
                        ITMAX);
}

class ZBrentHelper {
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ZBrent, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(RootFinder1D::TwoInitValNoDeriv);
        EMPTY_SHELL_METHOD(defaultZBrent);
        FIELD(tol, "tol");
        FIELD(ITMAX, "ITMAX");
        FIELD_MAKE_OPTIONAL(ITMAX);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }

    static IObject* defaultZBrent(){
        return new ZBrent();
    }
};

CClassConstSP const ZBrent::TYPE = CClass::registerClassLoadMethod(
    "ZBrent", typeid(ZBrent), ZBrentHelper::load);

DRLIB_END_NAMESPACE
