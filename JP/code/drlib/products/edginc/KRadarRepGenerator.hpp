#ifndef _TREERADAR_REP_GENERATOR_HPP
#define _TREERADAR_REP_GENERATOR_HPP

#include "edginc/DECLARE.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/FDProduct.hpp"
#include "edginc/FDModel.hpp"
#include "edginc/KComponent.hpp"
#include "edginc/MarketObservable.hpp"
#include "edginc/IFittingVarTransWrapper.hpp"
#include "edginc/IFuncBasisWrapper.hpp"

DRLIB_BEGIN_NAMESPACE

class PRODUCTS_DLL KRadarRepGenerator : public KComponent,
                        virtual public FDModel::IIntoProduct 
{
public:
    // Analogous to KSum except it will do regression in getValue and return null slice.   
    static CClassConstSP const TYPE;
    friend class KRadarRepGeneratorTree;      

    /****************** exported fields ************/
protected:
    CModel::IProdCreatorArray fittingKVec;
    CModel::IProdCreatorSP obsK;
    DateTimeArray obsDates;

    IFuncBasisWrapperSP funcBasis;
    IFittingVarTransWrapperSP fittingTransform;
    
    /****************** transient fields ************/
    /*IndexSpecFXArray fx;
    DateTimeArray valueDates;*/

    /****************** methods ************/
public:
    void add(IProdCreatorSP el);
    /*KRadarRepGenerator() : KComponent(TYPE), constant(0.0) {}*/
    KRadarRepGenerator() : KComponent(TYPE) {}

    virtual void validatePop2Object(void);

    //virtual void addValueDates(const DateTimeArray &valueDate);

protected:
    /* KComponent:: */
    virtual void setup(const IModel* model, const MarketData* market);
    /* FDModel::IIntoProduct:: */
    virtual FDProductSP createProduct(FDModel * model) const;

    /* IProdCreator */
    virtual double getValue(DateTime date, bool &canDo, CashflowInfo &cfi) const;

    /*KRadarRepGenerator(CClassConstSP const &type) : KComponent(type), constant(0.0) {}*/
    KRadarRepGenerator(CClassConstSP const &type) : KComponent(type) {}

    DateTimeArray getResetDates(IProdCreatorSP prod);

    void addDatesIndexSpec(CModel::IProdCreatorArray prodArr, DateTimeArray valueDates);

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new KRadarRepGenerator(TYPE); }
};

DECLARE(KRadarRepGenerator);

DRLIB_END_NAMESPACE
#endif
