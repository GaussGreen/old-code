//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CorrelationSkew.hpp
//
//   Description : Holds implied correlation parameters
//
//   Author      : Oliver Brockhaus
//
//   Date        : 15 Oct 2002
//
//
//----------------------------------------------------------------------------
#ifndef EDR_CORRELATIONSKEW_HPP
#define EDR_CORRELATIONSKEW_HPP

#include "edginc/MarketObject.hpp"
#include "edginc/PhiSkew.hpp"
#include "edginc/PhiSkewPower.hpp"

DRLIB_BEGIN_NAMESPACE

/* 
 * This is just a simple draft implementation of something which 
 * will be great and wonderful one day 
 */
class MARKET_DLL CorrelationSkew: public MarketObject,
                    public virtual PhiSkew::RestorableShift,
                    public virtual PhiSkewPower::RestorableShift 
{
public:
    static CClassConstSP const TYPE;
    friend class CorrelationSkewHelper;

    ~CorrelationSkew();

    /** Validation */
    void validatePop2Object();

    /** returns the correlationSkew as a simple double */
    double getCorrelationSkew() const;
    
    /** returns the correlationSkew as a simple double */
    double getCorrelationSkewPower() const;
    
    /** returns the correlationSkew's name */
    virtual string getName() const;

    // Sensitivity methods
    
    /** Returns the name of the stock/asset - used to determine
        whether to tweak the object */
    virtual string sensName(PhiSkew* shift)const;
    
    /** Shifts the object using given shift */    
    virtual bool sensShift(PhiSkew* shift);
    
    /** Restores the object to its original form */
    virtual void sensRestore(PhiSkew* shift);   

    /** Returns the name of the stock/asset - used to determine
        whether to tweak the object */
    virtual string sensName(PhiSkewPower* shift)const;
    
    /** Shifts the object using given shift */    
    virtual bool sensShift(PhiSkewPower* shift);
    
    /** Restores the object to its original form */
    virtual void sensRestore(PhiSkewPower* shift);   

    /** constructor */
    CorrelationSkew(const string& name,                 // optional
                    double        correlationSkew,
                    double        correlationSkewPower);

public:
    string          name;                   // optional
    double          correlationSkew;
    double          correlationSkewPower;
    CorrelationSkew();

};

typedef smartPtr<CorrelationSkew> CorrelationSkewSP;
typedef smartConstPtr<CorrelationSkew> CorrelationSkewConstSP;
#ifndef QLIB_CORRELATIONSKEW_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<CorrelationSkew>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<CorrelationSkew>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<CorrelationSkew>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<CorrelationSkew>);
#endif

// support for arrays of correlationSkews (note array of smart pointers)
typedef array<CorrelationSkewSP, CorrelationSkew> CorrelationSkewArray;
#ifndef QLIB_CORRELATIONSKEW_CPP
EXTERN_TEMPLATE(class MARKET_DLL array<CorrelationSkewSP _COMMA_ CorrelationSkew>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL array<CorrelationSkewSP _COMMA_ CorrelationSkew>);
#endif
// support for smart pointer to array
typedef smartConstPtr<CorrelationSkewArray> CorrelationSkewArrayConstSP;
typedef smartPtr<CorrelationSkewArray> CorrelationSkewArraySP;
#ifndef QLIB_CORRELATIONSKEW_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<CorrelationSkewArray>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<CorrelationSkewArray>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<CorrelationSkewArray>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<CorrelationSkewArray>);
#endif

// support for wrapper class
typedef MarketWrapper<CorrelationSkew> CorrelationSkewWrapper;
#ifndef QLIB_CORRELATIONSKEW_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<CorrelationSkew>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<CorrelationSkew>);
#endif

/** specialisations of arrayObjectCast */
template <> class MARKET_DLL arrayObjectCast<CorrelationSkewWrapper>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const CorrelationSkewWrapper& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(CorrelationSkewWrapper& value);

    /** Turns the IObjectSP into a DateTime */
    static CorrelationSkewWrapper fromIObject(IObjectSP& value);
};
// arrays of wrappers (note array of structures)
typedef array<CorrelationSkewWrapper, CorrelationSkewWrapper> CorrelationSkewWrapperArray;
#ifndef QLIB_CORRELATIONSKEW_CPP
EXTERN_TEMPLATE(class MARKET_DLL array<CorrelationSkewWrapper _COMMA_ CorrelationSkewWrapper>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL array<CorrelationSkewWrapper _COMMA_ CorrelationSkewWrapper>);
#endif


DRLIB_END_NAMESPACE
#endif // EDR_CORRELATIONSKEW_HPP
