//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Untweakable.hpp
//
//   Description : Result object to show that requested result failed
//
//   Author      : Andrew J Swain
//
//   Date        : 2 May 2001
//
//
//----------------------------------------------------------------------------

#ifndef UNTWEAKABLE_HPP
#define UNTWEAKABLE_HPP
#include "edginc/DECLARE.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/CombinableMixedResult.hpp"
#include <string>

using namespace std;


DRLIB_BEGIN_NAMESPACE

/** Result object to show that requested result failed */
class RISKMGR_DLL Untweakable: public CObject ,
                   public CombinableMixedResult{
public:
    static CClassConstSP const TYPE;

    Untweakable();
    Untweakable(const ModelException& e);
    Untweakable(const string& message);

    /** scale by factor x (implementation of CombinableResult) (does nothing)*/
    virtual void scale(double x);
    /** add Untweakable object to this result (Implementation of
        CombinableResult) (does nothing) */
    virtual void add(const CombinableResult& x, double scaleFactor);
    /** add an object (of general type) to this result. (implementation of 
     CombinableMixedResult) */
    virtual IObject* addResult(const  IObject& x,
                               double scaleFactor) const;

    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const;

    /** the message */
    string getMessage() const;

private:
    friend class UntweakableHelper;
    string message;
};

DECLARE(Untweakable)
#ifndef QLIB_UNTWEAKABLE_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<Untweakable>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<Untweakable>);
#endif

DRLIB_END_NAMESPACE
#endif
