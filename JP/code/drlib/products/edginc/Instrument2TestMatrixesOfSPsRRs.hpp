//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : Instrument2TestMatrixesOfSPsRRs.hpp
//
//   Description : Instrument we use to test matrices Of SPs and RRs
//
//   Date        : May 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_INSTRUMENT2TESTMATRIXESOFSPSRRS_HPP
#define QLIB_INSTRUMENT2TESTMATRIXESOFSPSRRS_HPP

#include "edginc/Instrument.hpp"
#include "edginc/CIDClosedFormModel.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

/******************************************************************************
*******************************************************************************
** Instrument we use to test matrices Of SPs and RRs                         **
*******************************************************************************
******************************************************************************/
class PRODUCTS_DLL Instrument2TestMatrixesOfSPsRRs 
                              : public CInstrument
                              , public virtual CIDClosedFormModel::IIntoProduct
{
public:
///////////////////////////////////////////////////////////////////////////////
//  A C C E S S O R S                                                        //
///////////////////////////////////////////////////////////////////////////////

    /** 
      * Outputs pricing date. This is an overridden method, but it is used 
      * to check the consistency of the value dates of the instrument and
      * the product. 
      * This method should be 'private' if value dates consistency is checked 
      * automatically *********************************************************
    **/
    DateTime    getValueDate() const;
    /** 
      * TRUE if closed form, FALSE if Monte-Carlo *****************************
    **/
    bool isClosedForm() const;
    /** 
      * Creates the state variable of the right type **************************
    **/ 
    // MCStateVarSPsRRsSP   getStateVar() const;
    /** 
      * Creates the state variable of the right type **************************
    **/ 
    StringArrayConstSP   getNames()    const;
    /** 
      * Creates the state variable of the right type **************************
    **/ 
    DateTimeArrayConstSP getDates()    const;

///////////////////////////////////////////////////////////////////////////////
//  R E F L E C T I O N   &   D L L                                          //
///////////////////////////////////////////////////////////////////////////////
    static CClassConstSP const TYPE;
    virtual ~Instrument2TestMatrixesOfSPsRRs(){;};
private:
    Instrument2TestMatrixesOfSPsRRs();
    static IObject* defaultConstructor();
    static void load(CClassSP& clazz);
    void validatePop2Object();
    /** 
      * Clone constructor and operator = are not defined. Do not use them *****
    **/
    Instrument2TestMatrixesOfSPsRRs(const Instrument2TestMatrixesOfSPsRRs&);
    void operator=(const Instrument2TestMatrixesOfSPsRRs&);

///////////////////////////////////////////////////////////////////////////////
//  O V E R R I D D E N   M E T H O D S                                      //
///////////////////////////////////////////////////////////////////////////////

    /** 
      * Creates a product for CIDClosedFormModel ******************************
    **/
    CIDClosedFormModel::IProductSP createProduct(
                                       const CIDClosedFormModel* gModel) const;
    /** 
      * Populates fields with market data if needed. NOT needed in this case **
    **/
    void GetMarket( const IModel* model,
                    const CMarketDataSP market){};
    /** 
      * Called after GetMarket(). NOT needed in this case *********************
    **/
    void Validate(){};
    /** 
      * YVXXX method of IInstrument. Is it called anywhere? *******************
    **/
    string discountYieldCurveName() const{return "";};

///////////////////////////////////////////////////////////////////////////////
// FIELDS                                                                    //
///////////////////////////////////////////////////////////////////////////////

    /** 
      *All names survival probabilities and recovery rates will be computed for
    **/
    StringArrayConstSP   names;
    /** 
      *All dates survival probabilities and recovery rates will be computed for
    **/
    DateTimeArrayConstSP dates;
    /** 
      * Date the instrument will be valued for ********************************
    **/
    DateTime             valueDate;
};
///////////////////////////////////////////////////////////////////////////////
DECLARE(Instrument2TestMatrixesOfSPsRRs);
///////////////////////////////////////////////////////////////////////////////
DRLIB_END_NAMESPACE
#endif // QLIB_INSTRUMENT2TESTMATRIXESOFSPSRRS_HPP
