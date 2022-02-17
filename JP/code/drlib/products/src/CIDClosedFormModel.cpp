//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : CIDClosedFormModel.cpp
//
//   Description : CID Closed Form Model
//
//   Date        : May 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CIDClosedFormModel.hpp"
#include "edginc/IInstrumentCollection.hpp"

DRLIB_BEGIN_NAMESPACE

///////////////////////////////////////////////////////////////////////////////
void CIDClosedFormModel::Price( CInstrument*  gInstrument, 
                                CControl*     gControl, 
                                CResults*     gResults)
{
    static const string routine =  "CIDClosedFormModel::getMarket";
    try
    {   
        // get an instance of IIntoProduct ------------------------------------
        IIntoProduct* _intoProd = dynamic_cast<IIntoProduct*>(gInstrument);
        QLIB_VERIFY(_intoProd != 0, 
                    "Cannot cast into CIDClosedFormModel::IIntoProduct");
        // create a product (product = instrument + model) --------------------
        IProductSP _product(_intoProd->createProduct(this));
        // finally, price it --------------------------------------------------
        _product->price(this, gControl, gResults);
    }
    catch (exception& e)
    {
        throw ModelException(e, routine);
    }
}
///////////////////////////////////////////////////////////////////////////////
void CIDClosedFormModel::getMarket(const MarketData*  gMarket,
                                   IInstrumentCollectionSP gInstruments)

{
    static const string routine =  "CIDClosedFormModel::getMarket";
    try
    {   
        fieldCIDParametersWrapper.getData(this, gMarket);
    }
    catch (exception& e)
    {
        throw ModelException(e, routine);
    }
}
///////////////////////////////////////////////////////////////////////////////
IModel::WantsRiskMapping CIDClosedFormModel::wantsRiskMapping() const
{
    static const string routine =  "CIDClosedFormModel::wantsRiskMapping";
    try
    {   
        return riskMappingIrrelevant;
    }
    catch (exception& e)
    {
        throw ModelException(e, routine);
    }
}
///////////////////////////////////////////////////////////////////////////////
CIDClosedFormModel::CIDClosedFormModel()
:   CModel(CIDClosedFormModel::TYPE)
,   fieldCIDParametersWrapper(CIDParameters::DEFAULTNAME)
{}
///////////////////////////////////////////////////////////////////////////////
CClassConstSP const CIDClosedFormModel::TYPE 
          = CClass::registerClassLoadMethod("CIDClosedFormModel", 
                                            typeid(CIDClosedFormModel), load);
///////////////////////////////////////////////////////////////////////////////
void CIDClosedFormModel::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(CIDClosedFormModel, clazz);
    SUPERCLASS(CModel);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(       fieldCIDParametersWrapper,
        string("Market Wrapper around CIDParameters. ")
        + "Default name is " + CIDParameters::DEFAULTNAME);
    FIELD_MAKE_OPTIONAL(fieldCIDParametersWrapper);
}
///////////////////////////////////////////////////////////////////////////////
IObject* CIDClosedFormModel::defaultConstructor()
{
    static const string routine = "CIDClosedFormModel::defaultConstructor";
    try
    {
        return new CIDClosedFormModel();
    }
    catch (exception& e)
    {
        throw ModelException(e, routine);
    }
}
///////////////////////////////////////////////////////////////////////////////
const CIDParameters & CIDClosedFormModel::getCIDParameters() const
{
    static const string routine = "CIDClosedFormModel::getCIDParameters";
    try
    {
        return (*fieldCIDParametersWrapper.getSP());
    }
    catch (exception& e)
    {
        throw ModelException(e, routine);
    }
}
///////////////////////////////////////////////////////////////////////////////
CClassConstSP const CIDClosedFormModel::IIntoProduct::TYPE 
      = CClass::registerClassLoadMethod("CIDClosedFormModel::IIntoProduct", 
                               typeid(CIDClosedFormModel::IIntoProduct), load);
///////////////////////////////////////////////////////////////////////////////
void CIDClosedFormModel::IIntoProduct::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(CIDClosedFormModel::IIntoProduct, clazz);
    IMPLEMENTS(CModel::IModelIntoProduct);
}
///////////////////////////////////////////////////////////////////////////////
DRLIB_END_NAMESPACE
