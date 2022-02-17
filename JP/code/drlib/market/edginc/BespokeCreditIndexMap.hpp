//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : BespokeCreditIndexMap.hpp
//
//   Description : Index mapping 
//
//   Author      : Gordon C Stephens
//
//   Date        : November 2005
//
//----------------------------------------------------------------------------

#ifndef FLEXIBLE_CREDIT_INDEX_MAP_HPP
#define FLEXIBLE_CREDIT_INDEX_MAP_HPP

#include "edginc/CreditIndexMap.hpp"
#include "edginc/CreditIndex.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL BespokeCreditIndexMap : public MarketObject,
                              public virtual ICreditIndexMap
{
    public:
        //------------------------------
        // BespokeCreditIndexMap methods
        //------------------------------
        static CClassConstSP const TYPE;

        //------------------------
        // ICreditIndexMap methods
        //------------------------

        /** Return the credit index associated to the single name
            as defined by this map */
        virtual CreditIndexBaseConstSP getIndex(const string& singleName) const;

        /** Return the single credit index defined by the map */
        virtual CreditIndexBaseConstSP getIndex() const;

        //---------------------
        // MarketObject methods
        //---------------------

        /** get the map name */
        virtual string getName() const;

        /** populate from market cache */
        virtual void getMarket(const IModel* model, const MarketData* market);

        //----------------
        // CObject methods
        //----------------

        /** checks parameters immediately after object is constructed */
        void validatePop2Object();

private:

        //------------------------------
        // BespokeCreditIndexMap methods
        //------------------------------

        BespokeCreditIndexMap();


        //infrastructure support
        static void load(CClassSP& clazz);
        static IObject* defaultConstructor();

        //fields
        string                      name;        //the name of the map
        ICDSParSpreadsWrapperArray  singleNames; //required to be unique
        CreditIndexWrapperArray     indexNames;  //one per singleName
                                                 //(will have repeated entries)

        //transient fields
        //to protect against getting multiple copies of the CreditIndex objects
        //when this object gets cloned (IntArray is much cheaper)
        IntArray                    indexPosition;   //to dereference resolvedIndices
                                                     //one per single name
        CreditIndexWrapperArray     resolvedIndices; //the actual CreditIndex's to use
};

DRLIB_END_NAMESPACE

#endif

