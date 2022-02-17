//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IAggregate.hpp
//
//   Description : Captures collapse of dimension from array to single value in various ways
//
//   Date        : Nov 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_IAGGREGATE_HPP
#define EDR_IAGGREGATE_HPP

#include "edginc/Object.hpp"
#include "edginc/IDoubleArray.hpp"
#include "edginc/IAssetFilter.hpp"

DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL IAggregate {
public:
    virtual double aggregate() = 0;

    virtual ~IAggregate() {};
};

typedef refCountPtr<IAggregate> IAggregateSP;

// Build IAggregate via a Maker class so the IAggregate can
// be identified once and aggregate() needs no parameters. It's a bit like
// building the thing we need in 2 stages.
class PRODUCTS_DLL IAggregateMaker : virtual public IObject {
public:
    static CClassConstSP const TYPE;
    
    virtual IAggregate* getAggregate(IDoubleArray* comps) = 0;

    // Some aggregators will allow assets to be dropped along the way.
    // Many will not.
    virtual IAggregate* getAggregate(IDoubleArray*   comps,
                                     IAssetFilterSP  filter) = 0;

    virtual ~IAggregateMaker() {};

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<IAggregateMaker> IAggregateMakerSP;
#ifndef QLIB_IAGGREGATE_CPP
EXTERN_TEMPLATE(class PRODUCTS_DLL_SP smartPtr<IAggregateMaker>);
#else
INSTANTIATE_TEMPLATE(class PRODUCTS_DLL smartPtr<IAggregateMaker>);
#endif

/*********************************************************************/
// Standard basket, but assumes components are supplied as perfs (for Pct)
class PRODUCTS_DLL BasketAggregateMaker : public CObject,
                             virtual public IAggregateMaker {
public:
    static CClassConstSP const TYPE;
    friend class BasketAggregate;
    friend class BasketEqualAggregate;
    friend class BasketAggregateMakerHelper; // for hiding the registration

    void validatePop2Object();

    BasketAggregateMaker(const DoubleArray& weights);

    virtual IAggregate* getAggregate(IDoubleArray* comps);

    virtual IAggregate* getAggregate(IDoubleArray*   comps,
                                     IAssetFilterSP  filter);

protected:
    DoubleArray weights;

private:
    BasketAggregateMaker(): CObject(TYPE), weights(0) {} // for reflection
    BasketAggregateMaker(const BasketAggregateMaker& rhs); // not implemented
    BasketAggregateMaker& operator=(const BasketAggregateMaker& rhs); // not implemented

};

typedef smartPtr<BasketAggregateMaker> BasketAggregateMakerSP;
typedef array<BasketAggregateMakerSP, BasketAggregateMaker> BasketAggregateMakerArray;
typedef smartPtr<BasketAggregateMakerArray> BasketAggregateMakerArraySP;


/*********************************************************************/
// Rainbowed version of above
class PRODUCTS_DLL RainbowAggregateMaker : public CObject,
                              virtual public IAggregateMaker {
public:
    static CClassConstSP const TYPE;
    friend class RainbowAggregate;
    friend class RainbowAggregateMakerHelper;

    RainbowAggregateMaker(const DoubleArray& weights);

    virtual IAggregate* getAggregate(IDoubleArray* comps);

    virtual IAggregate* getAggregate(IDoubleArray*   comps,
                                     IAssetFilterSP  filter);

protected:
    DoubleArray weights;
    
private:
    RainbowAggregateMaker(): CObject(TYPE) {} // for reflection
    RainbowAggregateMaker(const RainbowAggregateMaker& rhs); // not implemented
    RainbowAggregateMaker& operator=(const RainbowAggregateMaker& rhs); // not implemented
};
typedef smartPtr<RainbowAggregateMaker> RainbowAggregateMakerSP;
typedef array<RainbowAggregateMakerSP, RainbowAggregateMaker> RainbowAggregateMakerArray;
typedef smartPtr<RainbowAggregateMakerArray> RainbowAggregateMakerArraySP;

/*********************************************************************/
// Product (i.e. multiply) version
class PRODUCTS_DLL ProductAggregateMaker : public CObject,
                              virtual public IAggregateMaker {
public:
    static CClassConstSP const TYPE;
    friend class ProductAggregate;
    friend class ProductAggregateMakerHelper;

    ProductAggregateMaker();

    virtual IAggregate* getAggregate(IDoubleArray* comps);

    virtual IAggregate* getAggregate(IDoubleArray*   comps,
                                     IAssetFilterSP  filter);

private:
    ProductAggregateMaker(const ProductAggregateMaker& rhs); // not implemented
    ProductAggregateMaker& operator=(const ProductAggregateMaker& rhs); // not implemented

};
typedef smartPtr<ProductAggregateMaker> ProductAggregateMakerSP;
typedef array<ProductAggregateMakerSP, ProductAggregateMaker> ProductAggregateMakerArray;
typedef smartPtr<ProductAggregateMakerArray> ProductAggregateMakerArraySP;

/*********************************************************************/
// Reinvested version
class PRODUCTS_DLL ReinvestAggregateMaker : public CObject,
                               virtual public IAggregateMaker {
public:
    static CClassConstSP const TYPE;
    friend class ReinvestAggregate;
    friend class ReinvestAggregateMakerHelper;

    ReinvestAggregateMaker(const DoubleArray& weights);

    virtual IAggregate* getAggregate(IDoubleArray* comps);

    virtual IAggregate* getAggregate(IDoubleArray*   comps,
                                     IAssetFilterSP  filter);

protected:
    DoubleArray weights;
    
private:
    ReinvestAggregateMaker(): CObject(TYPE) {} // for reflection
    ReinvestAggregateMaker(const ReinvestAggregateMaker& rhs); // not implemented
    ReinvestAggregateMaker& operator=(const ReinvestAggregateMaker& rhs); // not implemented
};
typedef smartPtr<ReinvestAggregateMaker> ReinvestAggregateMakerSP;
typedef array<ReinvestAggregateMakerSP, ReinvestAggregateMaker> ReinvestAggregateMakerArray;
typedef smartPtr<ReinvestAggregateMakerArray> ReinvestAggregateMakerArraySP;


/*********************************************************************/
// ReinvestRainbowed version
class PRODUCTS_DLL ReinvestRainbowAggregateMaker : public CObject,
                                      virtual public IAggregateMaker {
public:
    static CClassConstSP const TYPE;
    friend class ReinvestRainbowAggregate;
    friend class ReinvestRainbowAggregateMakerHelper;

    ReinvestRainbowAggregateMaker(const DoubleArray& weights);

    virtual IAggregate* getAggregate(IDoubleArray* comps);

    virtual IAggregate* getAggregate(IDoubleArray*   comps,
                                     IAssetFilterSP  filter);

protected:
    DoubleArray weights;
    
private:
    ReinvestRainbowAggregateMaker(): CObject(TYPE) {} // for reflection
    ReinvestRainbowAggregateMaker(const ReinvestRainbowAggregateMaker& rhs); // not implemented
    ReinvestRainbowAggregateMaker& operator=(const ReinvestRainbowAggregateMaker& rhs); // not implemented
};

typedef smartPtr<ReinvestRainbowAggregateMaker> ReinvestRainbowAggregateMakerSP;
typedef array<ReinvestRainbowAggregateMakerSP, ReinvestRainbowAggregateMaker> ReinvestRainbowAggregateMakerArray;
typedef smartPtr<ReinvestRainbowAggregateMakerArray> ReinvestRainbowAggregateMakerArraySP;


/*********************************************************************/
// Returns lowest/highest number above/below some level in a list of double
class PRODUCTS_DLL ConditionalAggregateMaker : public CObject,
                                      virtual public IAggregateMaker {
public:
    static CClassConstSP const TYPE;
    friend class ConditionalAggregate;
    friend class ConditionalAggregateMakerHelper;

    ConditionalAggregateMaker(double level, bool isAbove);

    virtual IAggregate* getAggregate(IDoubleArray* comps);

    virtual IAggregate* getAggregate(IDoubleArray*   comps,
                                     IAssetFilterSP  filter);

protected:
    double  level;
    bool    isAbove;
    
private:
    ConditionalAggregateMaker(); // for reflection
    ConditionalAggregateMaker(const ConditionalAggregateMaker& rhs); // not implemented
    ConditionalAggregateMaker& operator=(const ConditionalAggregateMaker& rhs); // not implemented
};

typedef smartPtr<ConditionalAggregateMaker> ConditionalAggregateMakerSP;
typedef array<ConditionalAggregateMakerSP, ConditionalAggregateMaker> ConditionalAggregateMakerArray;
typedef smartPtr<ConditionalAggregateMakerArray> ConditionalAggregateMakerArraySP;

DRLIB_END_NAMESPACE

#endif

