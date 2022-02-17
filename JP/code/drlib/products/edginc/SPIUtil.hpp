//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPIUtil.hpp
//
//   Description : Synthetic Portfolio Insurance utilities
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#ifndef EDR_SPI_UTIL_HPP
#define EDR_SPI_UTIL_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/YieldCurve.hpp"

DRLIB_BEGIN_NAMESPACE

// a base class to handle the delayed validation in Pyramid
// so we trap after construction
class PRODUCTS_DLL SPIInterfaceIMS : public CObject {
public:
    static CClassConstSP const TYPE;
    SPIInterfaceIMS();// for reflection

	bool isValid();
    string errString();

protected:
	SPIInterfaceIMS(CClassConstSP clazz); // for inheritance 
	bool                 isOK;
    string               err;

private:
    SPIInterfaceIMS(const SPIInterfaceIMS& rhs); // not implemented
    SPIInterfaceIMS& operator=(const SPIInterfaceIMS& rhs); // not implemented

    static IObject* defaultSPIInterfaceIMS();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};

// for the output requests PAYMENT_DATES & KNOWN_CASHFLOWS

class PRODUCTS_DLL SPIKnownFlows {
public:
    // knowing 'today' allows any known but not yet past (i.e. unpaid) 
    // flows to be calculated and provided as additional value
    SPIKnownFlows(const DateTime& today,
                const YieldCurve* disc,
                const DateTime*   whenValued = 0);

    // MC handles its own PV so will be most convenient to get value
    // at maturity date
    // Actually we don't really want to do this inside the simulation
    // - better to add as a constant afterwards. Not sure we can do that...
    double getUnpaidValue();

	// XXX Might need to have a pair of dates - the payment("notification") date and
    // XXX the corresponding settlement date. For the moment (I think)
    // XXX this implementation assumes a flow is added only when it
    // XXX has been notified.
    void addFlow(const DateTime& when,
                 double          howMuch);

    const CashFlowArray* getFlows();

private:
    const DateTime&   today;
    CashFlowArraySP   knownFlows;
    const YieldCurve* disc; // reference
    const DateTime*   whenValued;
};
typedef refCountPtr<SPIKnownFlows> SPIKnownFlowsSP;

DRLIB_END_NAMESPACE

#endif
