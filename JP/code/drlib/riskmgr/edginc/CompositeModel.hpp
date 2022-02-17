//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CompositeModel.hpp
//
//   Description : Class that wraps one or more CModels. To be used when the 
//                 pricing is delegated to the instrument and possibly requires 
//                 more than one model
//
//   Author      : Regis Guichard
//
//   Date        : 31 March 2003
//
//
//----------------------------------------------------------------------------

#ifndef MODEL_COMBO_HPP
#define MODEL_COMBO_HPP
#include "edginc/Model.hpp"

DRLIB_BEGIN_NAMESPACE

/** Implementation of Model algorithm where the algorithm is specific to
    the instrument and yet is based on a log normal methodology */
class RISKMGR_DLL CompositeModel: public CModel{
public:
    static CClassConstSP const TYPE;
    friend class CompositeModelHelper;

    /** the class that the product must be able to create */
    class RISKMGR_DLL IProduct{
    public:
        virtual void price(CompositeModel* model,
                           Control*        control, 
                           CResults*       results) const = 0;
        virtual ~IProduct(){};
    };

    /** interface that the instrument must implement */
    class RISKMGR_DLL IIntoProduct: virtual public CModel::IModelIntoProduct {
    public:
        friend class CompositeModelHelper;
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(CompositeModel* model) const = 0;
    };

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);

    // should not be needed (therefore throws an exception)
    virtual MarketObjectSP GetMarket(const MarketData*    market,
                                     const string&        name,
                                     const CClassConstSP& type) const;
    
    /** throws exception if expected nb of models differ from actual */
    void checkNbModels(int    nbModelsExpected,
                       string className) const;

    /** returns the nber of models */
    int getNbModels() const;

    /** returns the i-th model */
    IModelSP getModel(int i) const;

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingIrrelevant, because (according to the note above) we
     * are using a "log normal methodology".  See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

private:
    CompositeModel(const CompositeModel &rhs);
    CompositeModel& operator=(const CompositeModel& rhs);

    /** for reflection */
    CompositeModel();
    
    // registered type
    CModelArray models;
};

DRLIB_END_NAMESPACE
#endif
