//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CDSCreditSupport.cpp
//
//   Description : Credit support for CDS instrument
//
//   Author      : Dmytro Zhuravytsky
//
//   Date        : 20 May 2005
//
//
//----------------------------------------------------------------------------

#ifndef CDS_CREDIT_SUPPORT_HPP
#define CDS_CREDIT_SUPPORT_HPP

#include "edginc/CreditSupport.hpp"

DRLIB_BEGIN_NAMESPACE

class CredDefSwap;

class PRODUCTS_DLL ParSpreadsUnd : public CreditUnd
{
public:
    ParSpreadsUnd(
        ICDSParSpreads * spreads,
        const CVolBase * vol );

    virtual string getName() const;
    virtual string getTrueName() const;
    virtual CVolProcessed * getProcessedVol(
        const CVolRequest * volRequest ) const;
    virtual void fwdValue(
        const DateTimeArray & dateList,
        CDoubleArray        & result ) const;
    virtual double getSpot() const;

    virtual void setSpot( double spot );
protected:
    ICDSParSpreads * m_spreads;
    const CVolBase * m_vol;
    double m_spot;
};

/** CDS credit support object  */
class PRODUCTS_DLL CDSCreditSupport : virtual public CreditSupport
{
public:

    /** preprocess instrument for a given set of path dates */
    virtual void preProcess(const DateTimeArray& dates, 
                            const DoubleArray& atmFwd, 
                            const DoubleArray& atmVar);

    /** return par spreads */
    virtual CreditUndSP getUnderlier() const;

    /** calculate values for a given path */
    virtual void calcPathValues(DoubleArray& results, const DateTimeArray& dates, 
                    const double* spots, double spotRef);

    /** return model for this instrument */
    virtual IModelSP getModel();

    /** return instrument ccy ISO code */
    virtual string getInstCcyCode() const;

    /** return instrument's last exposure date */   
    virtual DateTime getInstLastExposureDate() const;

    /** interpolate values from grid */
    void InterpValue(const CDoubleArray& spots, CDoubleArray& results);    

    /** store spot-value grid to be used for interpolation
        interp coeff (spline) pre-calculated */
    void SetGrid(const CDoubleArray& spots, const CDoubleArray& values);

    CDSCreditSupport(CInstrument* inst, CMarketDataSP market);

private:

    IModelSP model;
    
    CredDefSwap* instrCDSOrig;

    smartPtr<CredDefSwap> instrCDS;

    DoubleArray        AtmFwd;

    // spline coeff
    DoubleArrayArray y2;
    // pre-processed values for spline interpolation.
    DoubleArrayArray spotGrid;
    DoubleArrayArray valueGrid;

    // par spreads vol
    CVolBaseWrapper vol;
};

DRLIB_END_NAMESPACE

#endif
