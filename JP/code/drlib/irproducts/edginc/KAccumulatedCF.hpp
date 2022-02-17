//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : KAccumulatedCF.hpp
//
//   Description : Tarn component
//
//----------------------------------------------------------------------------
#ifndef _KAccumulatedCF_HPP
#define _KAccumulatedCF_HPP

#include "edginc/KComponent.hpp"
#include "edginc/Schedule.hpp"

DRLIB_BEGIN_NAMESPACE

// tarn component
class KAccumulatedCF : public KComponent,
                       virtual public FDModel::IIntoProduct
{
public:
    static CClassConstSP const TYPE;    

    void setup(const IModel* model, const MarketData* market);

    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;

    KAccumulatedCF(CClassConstSP const &type = TYPE) 
        : KComponent(type), numState(1), todayValue(0.) {}

private:
    friend class KAccumulatedCFTree;

    /****************** exported fields ************/
    int		        numState;   // number of states to price this product
    double          todayValue;  
    ScheduleSP      boundaryValues; // low / hi values for the cashflow (state variable)

    IProdCreatorSP  underlier; // the component that generates the cashflows 

    /****************** transient fields ************/

    DateTimeArray   cfDates;
    /****************** methods ************/
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new KAccumulatedCF(); }
};
typedef smartPtr<KAccumulatedCF> KAccumulatedCFSP;

DRLIB_END_NAMESPACE

#endif


