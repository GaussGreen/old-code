//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : NotApplicable.cpp
//
//   Description : Result object to show that requested result is N/A
//
//   Author      : Andrew J Swain
//
//   Date        : 2 May 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_NOTAPPLICABLE_CPP
#include "edginc/NotApplicable.hpp"
#include "edginc/Array.hpp"

DRLIB_BEGIN_NAMESPACE

NotApplicable::NotApplicable(): CObject(TYPE), cause(modelOrInstrument) {}
NotApplicable::NotApplicable(Cause cause): CObject(TYPE), cause(cause) {}

class NotApplicableHelper{
public:
/** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(NotApplicable, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(CombinableMixedResult);
        EMPTY_SHELL_METHOD(defaultNotApplicable);
    }
    
    static IObject* defaultNotApplicable(){
        return new NotApplicable();
    }
};

/** Overrides CObject implementation */
IObject* NotApplicable::clone() const{
    return new NotApplicable(cause);
}

/** write object out in 'output' format - ie suitable for comparing
    regression files with */
void NotApplicable::outputWrite(const string& linePrefix,
                                const string& prefix,
                                ostream&      stream) const{
    stream << linePrefix << prefix << ": " << "NotApplicable" << endl;
}

/** scale by factor x */
void NotApplicable::scale(double x) {}  // empty

/** add a NotApplicable to this result */
void NotApplicable::add(const CombinableResult& x, double scaleFactor) {
    // nothing to do provided x is of the right type
    dynamic_cast<const NotApplicable&>(static_cast<const IObject&>(x));
}
/** add an object (of general type) to this result. (implementation of 
    CombinableMixedResult) */
IObject* NotApplicable::addResult(const  IObject& x,
                                  double scaleFactor) const{
    // essentially NotApplicable + x = x. Only issue is make sure we can
    // do the scaling
    static const string method("NotApplicable::addResult");
    if (x.getClass()->isArray()){
        IArraySP arrayToAdd(&dynamic_cast<IArray&>(*x.clone()));
        for (int i = 0; i < arrayToAdd->getLength(); i++){
            IObjectSP elt(arrayToAdd->get(i));
            if (!CombinableResult::TYPE->isInstance(elt)){
                throw ModelException(method, "Cannot scale "+
                                     elt->getClass()->getName());
            }
            CombinableResult& cr = dynamic_cast<CombinableResult&>(*elt);
            cr.scale(scaleFactor);
        }
        return arrayToAdd.release();;
    }
    if (!CombinableResult::TYPE->isInstance(&x)){
        throw ModelException("NotApplicable::addResult", "Cannot scale "+
                             x.getClass()->getName());
    }
    CombinableResult& newResult = dynamic_cast<CombinableResult&>(*x.clone());
    newResult.scale(scaleFactor);
    return &newResult;
}

/** Returns cause of NotApplicable. See comments with enum */
NotApplicable::Cause NotApplicable::getCause() const{
    return cause;
}

CClassConstSP const NotApplicable::TYPE = CClass::registerClassLoadMethod(
    "NotApplicable", typeid(NotApplicable), NotApplicableHelper::load);

DRLIB_END_NAMESPACE
