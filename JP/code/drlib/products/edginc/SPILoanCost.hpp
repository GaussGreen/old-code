//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPILoanCost.hpp
//
//   Description : Loan Cost interface for SPI products
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#ifndef EDR_SPI_LOANCOST_HPP
#define EDR_SPI_LOANCOST_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Theta.hpp"
#include "edginc/YieldCurve.hpp"
#include "edginc/SPIUtil.hpp"
#include "edginc/SVGenDiscFactor.hpp"
#include "edginc/SVGenExpectedDiscFactor.hpp"

DRLIB_BEGIN_NAMESPACE

// this is the actual interface of what a loan cost object 
// needs to do internally
class PRODUCTS_DLL ILoanCostSPI {
public:

    virtual void init(const DateTime&   today,
                      const YieldCurve* ycForFixings,
                      double               Basis,
                      const DateTimeArray& rebalDates) = 0;

    virtual void refresh(const vector<SVExpectedDiscFactorSP>* dfSV) = 0;

    virtual void roll(const Theta::Util&   thetaUtil,
                      const YieldCurve*    yc,
                      double               Basis) = 0;

    virtual double getLiborRate(const YieldCurve* yc,
                                const DateTime&   startDate, 
                                const DateTime&   endDate, 
                                double            Basis) = 0;
    
    virtual double getAccrueFactor(const YieldCurve* yc,
                                   const DateTime&   startDate, 
                                   const DateTime&   endDate, 
                                   double            Basis) = 0;

    virtual double getFutureLiborRate(int iStep) = 0;
    virtual double getFutureAccrueFactor(int iStep) = 0;
    
    virtual CashFlowSP getLoanCostRateForPyramid() = 0;

    virtual ~ILoanCostSPI();
};
typedef refCountPtr<ILoanCostSPI> ILoanCostSPISP;

/*****************************************************************************/
// this is the external interface for abstraction so that users can bolt in 
// any Loan Cost type they want - note this includes the SPILoanCostWrapper
// which was necessary before we had abstraction in IMS 
// we yank out the real interface ILoanCostSPI as soon as possible
class PRODUCTS_DLL ILoanCostSPIInterface : virtual public IObject {
public:
    static CClassConstSP const TYPE;

    // get the internal interface
    virtual ILoanCostSPISP getLoanCostSPI() = 0;

private:
    static void load(CClassSP& clazz);
};
typedef smartPtr<ILoanCostSPIInterface> ILoanCostSPIInterfaceSP;

/*****************************************************************************/
class PRODUCTS_DLL SPILoanCostStd : public SPIInterfaceIMS,
                       virtual public ILoanCostSPI,
                       virtual public ILoanCostSPIInterface {
public:
    CashFlowArraySP  pastMMRates;
    double           pastSpread;      // borrowing spread added to past fixings
    double           futureSpread;    // spread added to future interpolated swap rates
    bool             spreadOnly;      // Libor is not included in the loan cost - only the spread
    
    //transient
    DateTime         today;
    const YieldCurve* ycForFixings; // $unregistered
    CashFlowSP       newFixingMMRate;  // for Pyramid processing $unregistered
    double           Basis;
    vector<int>      numDays; // between rebalDates.  $unregistered
    const vector<SVExpectedDiscFactorSP>* dfSV;

    // get the internal interface
    virtual ILoanCostSPISP getLoanCostSPI();

    void init(const DateTime&   today,
              const YieldCurve* ycForFixings,
              double            Basis,
              const DateTimeArray& rebalDates);

    void refresh(const vector<SVExpectedDiscFactorSP>* dfSV);

    // This interface says it all - this is a very poorly defined class. 
    // It should contain at least Basis and probably YieldCurve.
    void roll(const Theta::Util&   thetaUtil,
              const YieldCurve*    yc,
              double               Basis);
    
    double getLiborRate(const YieldCurve* yc,
                        const DateTime&   startDate, 
                        const DateTime&   endDate, 
                        double            Basis);

    double getAccrueFactor(const YieldCurve* yc,
                           const DateTime&   startDate, 
                           const DateTime&   endDate, 
                           double            Basis);

    // that applying from iStep-1 to iStep
    double getFutureLiborRate(int iStep);
    double getFutureAccrueFactor(int iStep);

    CashFlowSP getLoanCostRateForPyramid();

    static CClassConstSP const TYPE;
    SPILoanCostStd(); // for reflection

private:
    SPILoanCostStd(const SPILoanCostStd& rhs); // not implemented
    SPILoanCostStd& operator=(const SPILoanCostStd& rhs); // not implemented

    static IObject* defaultSPILoanCostStd();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};
typedef smartPtr<SPILoanCostStd> SPILoanCostStdSP;

/*****************************************************************************/
#define SPI_LOAN_COST_TYPE_STD   "Standard"

// the wrapper object which we needed before abstraction in IMS
class PRODUCTS_DLL SPILoanCostWrapper : public CObject,
                       virtual public ILoanCostSPIInterface {
public:
    static CClassConstSP const TYPE;

    string            SPILoanCostType;
    SPILoanCostStdSP  loanCostStd;

public:

    // get the internal interface
    virtual ILoanCostSPISP getLoanCostSPI();

    // validation
    void validatePop2Object();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
     
    // for reflection
    SPILoanCostWrapper();

    static IObject* defaultSPILoanCostWrapper();
};
typedef smartPtr<SPILoanCostWrapper> SPILoanCostWrapperSP;

DRLIB_END_NAMESPACE

#endif
