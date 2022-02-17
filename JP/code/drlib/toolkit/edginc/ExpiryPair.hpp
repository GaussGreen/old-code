/**
 * @file ExpiryPair.hpp
 */

#ifndef QLIB_ExpiryPair_H
#define QLIB_ExpiryPair_H

#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ExpiryPair)
FORWARD_DECLARE(Expiry)


typedef smartConstPtr<ExpiryPair> ExpiryPairConstSP;
typedef smartPtr<ExpiryPair> ExpiryPairSP;
typedef array<ExpiryPairSP, ExpiryPair> ExpiryPairArray;
typedef smartPtr<ExpiryPairArray> ExpiryPairArraySP;
typedef smartConstPtr<ExpiryPairArray> ExpiryPairArrayConstSP;

#ifndef QLIB_EXPIRYPAIR_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<ExpiryPair>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<ExpiryPair>);
EXTERN_TEMPLATE(class TOOLKIT_DLL array<ExpiryPairSP _COMMA_ ExpiryPair>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<ExpiryPairArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<ExpiryPairArray>);
#else
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<ExpiryPair>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<ExpiryPair>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL array<ExpiryPairSP _COMMA_ ExpiryPair>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<ExpiryPairArray>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<ExpiryPairArray>);
#endif

/**
 * A pair of expiries: one for the option, and another one for the underlying. 
 * Used as a qualifier in computing credit vega sensitivities 
 */

class TOOLKIT_DLL ExpiryPair: public CObject {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

    /**
     * Field which specifies the option expiry
     */

    ExpiryConstSP optExpiry;

    /**
     * Field which specifies the underlying expiry
     */

    ExpiryConstSP ulExpiry;

    /**
     * Constructor.
     *
     * See also SP(), around().
     */

    ExpiryPair(ExpiryConstSP optExpiry, ExpiryConstSP ulExpiry);

    /**
     * Constructor returning a smartPtr.
     */

    static ExpiryPairSP SP(ExpiryConstSP optExpiry, ExpiryConstSP ulExpiry);
    
    /** Returns true if given expiry pair matches this exactly  */
    bool equals(const ExpiryPair* expiry) const;

    string toString() const;

    ~ExpiryPair();
};

DRLIB_END_NAMESPACE

#endif
