//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DerivativeAsset.hpp
//
//   Description : interface class for a derivative asset (e.g. AssetCVB)
//
//   Author      : Jay Blumenstein
//
//   Date        : 18 Sep 2002
//
//
//----------------------------------------------------------------------------

#ifndef DERIVATIVE_ASSET_HPP
#define DERIVATIVE_ASSET_HPP

#include "edginc/Control.hpp"
#include "edginc/Model.hpp"
#include "edginc/Instrument.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface objects need to implement in order for them to have a
    theo value and a mtm value. (Don't think interface class is
    restricted to assets despite its name although the objects will
    typically be assets) */
class RISKMGR_DLL DerivativeAsset: virtual public IObject {
public:
    static CClassConstSP const TYPE;

    /** Identifier for MTM results when we are doing ASSET_THEO_MTM. The
        results live in the instrument packet using the supplied identifier */
    static string ASSET_MTM_ID; 

    /** Creates a model that should be used to price a product which has
        a DerivativeAsset in it. Note copy of model is not made */
    static IModel* createModel(const IModelSP& model, 
                               IInstrumentCollectionSP insts);

    /** Returns true if the instrument contains derivatives assets */
    static bool derivativeAssetsExist(IInstrumentCollectionSP insts);

    /** Sets all DerivativeAssets in inst to use supplied value for
        whether to use theoretical price or mtm.   */
    static void setUseTheoAssetPrice(IInstrumentCollectionSP insts,
                                     bool useTheoAssetPrice);

    /** sets whether the asset should use the theoretical price or the
        mtm price. This is called once per block of pricings. At this
        point all market data has been obtained and all other preparatory
        work (eg roll to now) has been performed */
    virtual void setUseTheoAssetPrice(bool useTheoAssetPrice) = 0;

    /** sets the control ptr in the asset to the inputted
        control. This method is invoked on each DerivativeAsset just
        before a pricing call is made to the instrument. The supplied
        control is guaranteed to be valid indefinitely. Note that the
        model is passed here to ensure the right model is being used
        for the right pricings:- some engines store information across
        pricing calls in the model whilst some tweaking code clones
        the model and runs a new pricing call. The DerivativeAsset should
        not clone the supplied inputs - they should be used as is */
    virtual void setControlAndModel(const CControlSP& ctrl,
                                    const IModelSP&   model) = 0;

    /** Returns the name of the derivative asset - this is used to manage
        the asset's model (ie model used to calculate asset's theo price).
        Important when there are multiple DerivativeAssets */
    virtual string getName() const = 0;

    /* Returns the model used to calculate asset's theo price. This is used
       to allow the DerivativeAsset model to manage the model. See 
       setControlAndModel(). */
    virtual IModelSP getModel() const = 0;

    /* Returns the instrument used to calculate asset's theo price. */
    virtual CInstrumentSP getInstrument() = 0;

    /* Returns the MTM price of the derivative asset if it exists */
    virtual double getMTM() const = 0;

    virtual void getMarket(const IModel* model, const MarketData* market) = 0;
        
private:
    class Model;
    class UpdateAssetWithControlAndModel;
    class UpdateAssetWithUseTheoFlag;
    class FindDerivAssets;
    static void load(CClassSP& clazz);
};

typedef smartPtr<DerivativeAsset> DerivativeAssetSP;

DRLIB_END_NAMESPACE

#endif
