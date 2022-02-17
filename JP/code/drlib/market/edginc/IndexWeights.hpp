//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IndexWeights.hpp
//
//   Description : A container for (IndexSkewWrapper, Weight) pairs
//
//   Author      : Antoine Gregoire
//
//   Date        : March 2005
//
//
//----------------------------------------------------------------------------

#ifndef EDR_INDEX_WEIGHTS_HPP
#define EDR_INDEX_WEIGHTS_HPP

#include "edginc/IndexSkew.hpp"

DRLIB_BEGIN_NAMESPACE

/** 
 * IndexWeight is just a container for (IndexSkewWrapper, Weight) pairs
 * */
class MARKET_DLL IndexWeights : public CObject, public virtual IGetMarket {
public:
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Destructor */
    ~IndexWeights();
    
    /** Called immediately after object constructed */
    virtual void validatePop2Object();

    /** public constructor */
    IndexWeights(IndexSkewWrapperArraySP indexNames, CDoubleArraySP weights);

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);
    
    /** Gets index names */
    const IndexSkewWrapperArraySP getIndexNames() const;
    
    /** Gets index weights */
    const CDoubleArraySP getWeights() const;
    
private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
    
    /** Only build instances of that class using reflection */
    IndexWeights();
    
    /** Default constructor */
    static IObject* defaultConstructor();

    // ----------------
    // MANDATORY FIELDS
    // ----------------
    
    /** Array of IndexSkew wrappers */
    IndexSkewWrapperArraySP indexNames;
    
    /** Array of IndexSkew weights */
    CDoubleArraySP weights;    
};

DECLARE(IndexWeights);

DRLIB_END_NAMESPACE

#endif
