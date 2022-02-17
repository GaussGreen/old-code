//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PVDivCreditSupport.hpp
//
//   Description : Credit support object for PVDiv
//
//   Author      : Ning Shen
//
//   Date        : 10 June 2002
//
//
//----------------------------------------------------------------------------

#ifndef PVDDIV_CREDIT_HPP
#define PVDDIV_CREDIT_HPP

#include "edginc/CreditSupport.hpp"

DRLIB_BEGIN_NAMESPACE

class CPVDiv;
/** PVDiv credit support object  */
class PRODUCTS_DLL PVDivCreditSupport : virtual public CreditSupport
{
public:
 
    
    /** some comment */
    virtual void preProcess(const DateTimeArray& dates, 
							const DoubleArray& atmFwd, 
							const DoubleArray& atmVar);

    virtual	IModelSP getModel();

    /** return instrument ccy ISO code */
    virtual string getInstCcyCode() const;

    /** return instrument's last exposure date */
    virtual DateTime getInstLastExposureDate() const;

    /** return asset */
    virtual CreditUndSP getUnderlier() const;

	virtual void calcPathValues(DoubleArray& results,
                                const DateTimeArray& dates,
                                const double* spots,
                                double spotRef);

	PVDivCreditSupport(CInstrument* inst, CMarketDataSP market);
protected:
        PVDivCreditSupport();  
private:
    IModelSP model;

	CPVDiv* pvDivOriginal;
    smartPtr<CPVDiv> pvDiv; // we'll need to reset the instrument after it's been modified

    
    void getFixingDates(DateTimeArray& allFixingDates);

    bool useLinear; // when PVDiv is linear (i.e. nominal case), we use price/delta cache for computing price
    DoubleArray priceCache;
    DoubleArray deltaCache;
	
//	void zeroPaidDividends(const DateTime& creditDate);
//	void convertToNominal();
//	void convertYield2Dollar();

};

DRLIB_END_NAMESPACE

#endif
