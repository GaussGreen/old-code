//----------------------------------------------------------------------------
//
//   Group       : GCCG DR
//
//   Filename    : ClosedFormEnergy.hpp
//
//   Description : Input info for Energy closed form model
//
//   Author      : Sean Chen
//
//   Date        : August 25, 2005
//
//
//----------------------------------------------------------------------------

#ifndef CLOSEDFORMENERGY_HPP
#define CLOSEDFORMENERGY_HPP
#include "edginc/Model.hpp"

DRLIB_BEGIN_NAMESPACE


class PRODUCTS_DLL ClosedFormEnergy: public CModel{
public:
    static CClassConstSP const TYPE;
    friend class ClosedFormEnergyHelper;

    /** the class that the product must be able to create */
    class PRODUCTS_DLL IProduct{
    public:
        virtual void price(ClosedFormEnergy* model,
                           Control*    control, 
                           CResults*   results) const = 0;
    };

    /** interface that the instrument must implement */
    class PRODUCTS_DLL IIntoProduct{
    public:
        friend class ClosedFormEnergyHelper;
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(ClosedFormEnergy* model) const = 0;
    };

    /** Simple constructor */
    ClosedFormEnergy();

    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);
    
	bool isSmileOff() const { return smileOff; }

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingIrrelevant, because we don't price using a parametric
     * (latent, non-observable) pdf that needs risk mapping.  See
     * IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

private:
    ClosedFormEnergy(const ClosedFormEnergy &rhs);
    ClosedFormEnergy& operator=(const ClosedFormEnergy& rhs);

    bool  smileOff;
};

DRLIB_END_NAMESPACE
#endif
