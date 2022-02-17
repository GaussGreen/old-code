//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : CIDClosedFormModel.hpp
//
//   Description : CID Closed Form Model
//
//   Date        : May 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_CIDCLOSEDFORMMODEL_HPP
#define QLIB_CIDCLOSEDFORMMODEL_HPP

#include "edginc/Model.hpp"
#include "edginc/CIDParameters.hpp"
#include "edginc/VirtualDestructorBase.hpp"
//#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/******************************************************************************
*******************************************************************************
** CID Closed Form Model                                                     **
*******************************************************************************
******************************************************************************/

class PRODUCTS_DLL CIDClosedFormModel : public CModel
{
public:
    /************************************************************************** 
     * Instruments that can be priced with CIDClosedFormModel
     * should implement this interface 
    **************************************************************************/
    class IIntoProduct;
    /************************************************************************** 
    * All products = [an instrument + CIDClosedFormModel]
    * should implement this interface 
    **************************************************************************/
    class IProduct;
    typedef smartPtr<IProduct>      IProductSP;
    typedef smartConstPtr<IProduct> IProductConstSP;

    /** 
      * Returns the CIDParameters object inside CIDParametersWrapper ********
    **/
    const CIDParameters & getCIDParameters() const;

///////////////////////////////////////////////////////////////////////////////
//  R E F L E C T I O N   &   D L L                                          //
///////////////////////////////////////////////////////////////////////////////
    static CClassConstSP const TYPE;
    ~CIDClosedFormModel(){;};
private:
    CIDClosedFormModel();
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);
    /** 
      * Clone constructor and operator = are not defined. Do not use them *****
    **/
    CIDClosedFormModel(const CIDClosedFormModel&);
    void operator=(const CIDClosedFormModel&);
    
///////////////////////////////////////////////////////////////////////////////
//  V I R T U A L   A S S E S S O R S                                        //
///////////////////////////////////////////////////////////////////////////////

    /** 
      * Risk mapping **********************************************************
    **/
    WantsRiskMapping wantsRiskMapping() const;

///////////////////////////////////////////////////////////////////////////////
//  V I R T U A L   M O D I F I E R S                                        //
///////////////////////////////////////////////////////////////////////////////

    /** 
      * Invoked after instrument has got its market data. 
      * Allows model to get any extra data required. 
      * Default implementation does nothing ***********************************
    **/
    void getMarket(const MarketData*  gMarket,
                         IInstrumentCollectionSP gInstruments);
    /** 
      * Instantiates a product and prices it **********************************
    **/
    void Price( CInstrument*  gInstrument, 
                CControl*     gControl, 
                CResults*     gResults);

///////////////////////////////////////////////////////////////////////////////
//  F I E L D S                                                              //
///////////////////////////////////////////////////////////////////////////////
    /** 
      * Contains all market data we need to price *****************************
    **/
    CIDParametersWrapper    fieldCIDParametersWrapper;
};
//DECLARE(CIDClosedFormModel);

/******************************************************************************
*******************************************************************************
** Instruments that can be priced with CIDClosedFormModel                    **
** should implement this interface                                           **
*******************************************************************************
******************************************************************************/
class PRODUCTS_DLL CIDClosedFormModel::IIntoProduct
                                     : public virtual CModel::IModelIntoProduct
{
public:
    /** 
      * Main method of this interface. Creates a product **********************
    **/
    virtual IProductSP createProduct(const CIDClosedFormModel * gModel)const=0;
///////////////////////////////////////////////////////////////////////////////
//  R E F L E C T I O N   &   D L L                                          //
///////////////////////////////////////////////////////////////////////////////
    static CClassConstSP const TYPE;
protected:
    IIntoProduct(){;};
    virtual ~IIntoProduct(){;};
private:
    static void load(CClassSP& clazz);
    /** 
      * Clone constructor and operator = are not defined. Do not use them *****
    **/
    IIntoProduct(const IIntoProduct&);
    void operator=(const IIntoProduct&);
};

/****************************************************************************** 
*******************************************************************************
** All products = [an instrument + CIDClosedFormModel]                       **
** should implement this interface                                           **
*******************************************************************************
******************************************************************************/
class PRODUCTS_DLL CIDClosedFormModel::IProduct : public VirtualDestructorBase
{
public:
    /** 
      * Pricing happens inside this method ************************************
    **/
    virtual void price( CIDClosedFormModel* gModel,
                        Control           * gControl, 
                        CResults          * gResults) const = 0;
///////////////////////////////////////////////////////////////////////////////
//  R E F L E C T I O N   &   D L L                                          //
///////////////////////////////////////////////////////////////////////////////
    IProduct(){;};
    virtual ~IProduct(){;};
private:
    /** 
      * Clone constructor and operator = are not defined. Do not use them *****
    **/
    IProduct(const IProduct&);
    void operator=(const IProduct&);
};
///////////////////////////////////////////////////////////////////////////////
DRLIB_END_NAMESPACE
#endif // QLIB_CIDCLOSEDFORMMODEL_HPP
