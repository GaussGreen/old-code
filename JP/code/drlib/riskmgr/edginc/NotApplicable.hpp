//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : NotApplicable.hpp
//
//   Description : Result object to show that requested result is N/A
//
//   Author      : Andrew J Swain
//
//   Date        : 2 May 2001
//
//
//----------------------------------------------------------------------------

#ifndef NOTAPPLICABLE_HPP
#define NOTAPPLICABLE_HPP
#include "edginc/Object.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/CombinableMixedResult.hpp"

using namespace std;


DRLIB_BEGIN_NAMESPACE

/** Result object to show that requested result is N/A e.g. asking for
    naked bond rho for a stock */
class RISKMGR_DLL NotApplicable: public CObject,
                     public CombinableMixedResult {
public:
    static CClassConstSP const TYPE;

    /** Gives limited information as to the cause of something being
        not applicable. Please note that this currently is only for 
        internal use (basically to help the results aggregation). Should
        we ever wish to expose this externally we should probably 
        scrap the enum */
    enum Cause{
        modelOrInstrument = 0,
        aggregation
    };

    /** Cause is set to modelOrInstrument. See comments with enum */
    NotApplicable();

    NotApplicable(Cause cause); // See comments with enum 
    /** Overrides CObject implementation */
    virtual IObject* clone() const;

    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const;

    /** scale by factor x */
    virtual void scale(double x);

    /** add NotApplicable object to this result (Implementation of
     CombinableResult) */
    virtual void add(const CombinableResult& x, double scaleFactor);
    /** add an CombinableResult (of general type) to this
     result. (implementation of CombinableMixedResult) */
    virtual IObject* addResult(const  IObject& x, 
                               double scaleFactor) const;
    /** Returns cause of NotApplicable. See comments with enum */
    Cause getCause() const;

private:
    const Cause  cause; // $unregistered
    friend class NotApplicableHelper;
};

typedef smartPtr<NotApplicable> NotApplicableSP;
#ifndef QLIB_NOTAPPLICABLE_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<NotApplicable>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<NotApplicable>);
#endif

DRLIB_END_NAMESPACE
#endif
