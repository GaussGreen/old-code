//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Interpolator.cpp
//
//   Description : 
//
//   Date        : 06 June 2002
//
//
//   $Log: Interpolator.cpp,v $
//   Revision 1.9  2004/06/01 16:50:09  mvenardos
//   Redesign to make interpolant either virtual or non
//   virtual
//
//   Revision 1.8  2004/03/08 11:47:17  mrobson
//   Ensure array names are 'consistent' with component type (DRI requirement)
//
//   Revision 1.7  2004/02/18 17:41:40  mrobson
//   Use new macros for interfaces
//
//   Revision 1.6  2004/02/17 19:50:50  mrobson
//   Use correct register method for array
//
//   Revision 1.5  2003/05/21 13:35:23  rguichar
//   Made some methods const, as they should
//
//   Revision 1.4  2002/10/08 14:40:42  evenardo
//   Added constructor to Interpolant
//
//   Revision 1.3  2002/10/08 08:45:41  evenardo
//   Interpolant is a CObject
//
//   Revision 1.2  2002/08/07 17:17:35  evenardo
//   Addin functions: GET_INTERPOLANT and GET_INTERPOLATED_VALUES
//
//   Revision 1.1  2002/06/07 14:07:05  evenardo
//   Interpolator and Interpolant interfaces
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_INTERPOLATOR_CPP
#include "edginc/Interpolator.hpp"
#include "edginc/Addin.hpp"


DRLIB_BEGIN_NAMESPACE
// InterpolantBase
template<> CClassConstSP const InterpolantBase<Interpolator::Interpolant>::TYPE = 
CClass::registerClassLoadMethod(
    "InterpolantBase<Interpolator::Interpolant>", 
    typeid(InterpolantBase<Interpolator::Interpolant>), 
    InterpolantBase<Interpolator::Interpolant>::load);

static CClassConstSP vc7BugFix1 = 
InterpolantBase<Interpolator::Interpolant>::TYPE;

template<> CClassConstSP const InterpolantBase<Interpolator::NullBase>::TYPE = 
CClass::registerClassLoadMethod(
    "InterpolantBase<Interpolator::NullBase>", 
    typeid(InterpolantBase<Interpolator::NullBase>), 
    InterpolantBase<Interpolator::NullBase>::load);
static CClassConstSP vc7BugFix2 = InterpolantBase<Interpolator::NullBase>::TYPE;


// Interpolator
void Interpolator::load(CClassSP& clazz){
    REGISTER_INTERFACE(Interpolator, clazz);
    EXTENDS(IObject);
}    

CClassConstSP const Interpolator::TYPE = CClass::registerInterfaceLoadMethod(
    "Interpolator", typeid(Interpolator), Interpolator::load);

////////////////////////////////////////////////////////////////////////////////

// Interpolant
void Interpolator::Interpolant::load(CClassSP& clazz){
    REGISTER_INTERFACE(Interpolator::Interpolant, clazz);
    EXTENDS(IObject);
}

CClassConstSP const Interpolator::Interpolant::TYPE = CClass::registerInterfaceLoadMethod(
    "Interpolator::Interpolant", typeid(Interpolator::Interpolant), Interpolator::Interpolant::load);

////////////////////////////////////////////////////////////////////////////////

// NullBase
void Interpolator::NullBase::load(CClassSP& clazz){
    REGISTER_INTERFACE(Interpolator::NullBase, clazz);
    EXTENDS(IObject);
}

CClassConstSP const Interpolator::NullBase::TYPE = CClass::registerInterfaceLoadMethod(
    "Interpolator::NullBase", typeid(Interpolator::NullBase), Interpolator::NullBase::load);

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

DEFINE_TEMPLATE_TYPE_WITH_NAME("Interpolator::InterpolantArray", InterpolantArray);

////////////////////////////////////////////////////////////////////////////////

/** ADDIN method for getting an Interpolant */
class GetInterpolantAddin: public CObject{
    static CClassConstSP const TYPE;

    InterpolatorSP interpolator;
    DoubleArray   x;
    DoubleArray   y;

    static IObjectSP getInterpolant(GetInterpolantAddin* params){
        static const string routine = "GetInterpolantAddin::getInterpolant";
        try {
            Interpolator::InterpolantConstSP interpolant(params->interpolator->computeInterp(params->x, params->y));
            
            return Interpolator::InterpolantSP::constCast(interpolant);

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }


    /** for reflection */
    GetInterpolantAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(GetInterpolantAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetInterpolantAddin);
        FIELD(interpolator, "Interpolator");
        FIELD(x, "x");
        FIELD(y, "y");

        Addin::registerClassObjectMethod("GET_INTERPOLANT",
                                         Addin::RISK,
                                         "Get interpolant",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)getInterpolant);

    }

    static IObject* defaultGetInterpolantAddin(){
        return new GetInterpolantAddin();
    }
};

CClassConstSP const GetInterpolantAddin::TYPE = CClass::registerClassLoadMethod(
    "GetInterpolantAddin", typeid(GetInterpolantAddin), load);

////////////////////////////////////////////////////////////////////////////////

/** ADDIN method for interpolating with some interpolant */
class GetInterpolatedValuesAddin: public CObject{
    static CClassConstSP const TYPE;

    Interpolator::InterpolantSP interpolant;
    DoubleArray   x;
    int n;

    static IObjectSP computeValue(GetInterpolatedValuesAddin* params){
        static const string routine = "GetInterpolatedValuesAddin::computeValue";
        try {
            DoubleArray y(params->x.size());
            params->interpolant->value(params->x, params->n, y);
            
            return IObjectSP(y.clone());

        } catch (exception& e){
            throw ModelException(e, routine);
        }
    }


    /** for reflection */
    GetInterpolatedValuesAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(GetInterpolatedValuesAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultGetInterpolatedValuesAddin);
        FIELD(interpolant, "Interpolant");
        FIELD(x, "x");
        FIELD(n, "n");

        Addin::registerClassObjectMethod("GET_INTERPOLATED_VALUES",
                                         Addin::RISK,
                                         "Get interpolated values",
                                         TYPE,
                                         false,
                                         Addin::expandSimple,
                                         (Addin::ObjMethod*)computeValue);

    }

    static IObject* defaultGetInterpolatedValuesAddin(){
        return new GetInterpolatedValuesAddin();
    }
};

CClassConstSP const GetInterpolatedValuesAddin::TYPE = CClass::registerClassLoadMethod(
    "GetInterpolatedValuesAddin", typeid(GetInterpolatedValuesAddin), load);

bool InterpolatorLoad() {
    return Interpolator::TYPE != NULL;
}

DRLIB_END_NAMESPACE
