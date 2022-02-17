//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : CIDClosedFormProduct2TestMatrixesOfSPsRRs.hpp
//
//   Description : Instrument we use to test matrices Of SPs and RRs
//                                                     + CID Closed Form Model
//
//   Date        : Aug 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_CIDCLOSEDFORMPRODUCT2TESTMATRIXESOFSPSRRS_HPP
#define QLIB_CIDCLOSEDFORMPRODUCT2TESTMATRIXESOFSPSRRS_HPP

#include "edginc/Instrument2TestMatrixesOfSPsRRs.hpp"
#include "edginc/CIDClosedFormModel.hpp"

DRLIB_BEGIN_NAMESPACE

/******************************************************************************
*******************************************************************************
** Instrument we use to test matrices Of SPs and RRs                         **
*******************************************************************************
******************************************************************************/
class PRODUCTS_DLL CIDClosedFormProduct2TestMatrixesOfSPsRRs 
                                          : public CIDClosedFormModel::IProduct
{
public:
///////////////////////////////////////////////////////////////////////////////
//  R E F L E C T I O N   &   D L L                                          //
///////////////////////////////////////////////////////////////////////////////
    CIDClosedFormProduct2TestMatrixesOfSPsRRs(  
                          const Instrument2TestMatrixesOfSPsRRs * gInstrument);
    virtual ~CIDClosedFormProduct2TestMatrixesOfSPsRRs(){;};
private:
    /** 
      * Pricing happens inside this method ************************************
    **/
    void price( CIDClosedFormModel * gModel,
                CControl           * gControl, 
                CResults           * gResults) const;
    /** 
      * Clone constructor and operator = are not defined. Do not use them *****
    **/
    CIDClosedFormProduct2TestMatrixesOfSPsRRs(
                             const CIDClosedFormProduct2TestMatrixesOfSPsRRs&);
    void operator=(const CIDClosedFormProduct2TestMatrixesOfSPsRRs&);

    const Instrument2TestMatrixesOfSPsRRs * instrument;
};
///////////////////////////////////////////////////////////////////////////////
DRLIB_END_NAMESPACE
#endif // QLIB_CIDCLOSEDFORMPRODUCT2TESTMATRIXESOFSPSRRS_HPP
