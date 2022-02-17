//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SPIAlgorithm.hpp
//
//   Description : Algorithm interface for SPI products
//
//   Date        : Dec 2005
//
//
//----------------------------------------------------------------------------
#ifndef EDR_SPI_ALGO_HPP
#define EDR_SPI_ALGO_HPP

#include "edginc/config.hpp"
#include "edginc/SPIFees.hpp"
#include "edginc/SPICutoff.hpp"
#include "edginc/SPIUtil.hpp"

DRLIB_BEGIN_NAMESPACE

// this is the actual interface of what an algorithm object 
// needs to do internally
class PRODUCTS_DLL IAlgorithmSPI {
public:

    virtual bool doRebalance(double UL,
                             double SL,
                             int    iStep) = 0;
    virtual const DoubleArray& rebalWeights(double crash) = 0;
    virtual double unbalCrash(const DoubleArray& nE,
                              const DoubleArray& E,
                              int                iStep) const = 0;
    virtual double sustCrash(double buffer, int iStep) const = 0;
    virtual double targetExp(double SE) const = 0;
    virtual double feeAtMin(IFeesSPIConstSP fees) const = 0;
    virtual double gapRiskAtStep(double B, double BF, double equityComp, int iStep) const = 0;
    virtual void crossValidate(ICutoffSPISP cutoff, IFeesSPISP fees, 
                               bool isRainbowSPI) = 0;
    virtual int getNumRiskyAssets() const = 0; // how many risky assets involved in rebalancing?
    virtual void init(const DoubleArray* bondPrices) = 0;
    virtual ~IAlgorithmSPI();
};
typedef refCountPtr<IAlgorithmSPI> IAlgorithmSPISP;

// this is the external interface for abstraction so that users can bolt in 
// any Algorithm type they want - note this includes the SPIAlgorithmWrapper
// which was necessary before we had abstraction in IMS 
// we yank out the real interface IAlgorithmSPI as soon as possible
class PRODUCTS_DLL IAlgorithmSPIInterface : virtual public IObject {
public:
    static CClassConstSP const TYPE;

    virtual IAlgorithmSPISP getAlgorithmSPI() = 0;

private:
    static void load(CClassSP& clazz);
};
typedef smartPtr<IAlgorithmSPIInterface> IAlgorithmSPIInterfaceSP;

class PRODUCTS_DLL SPIAlgorithmStd : public SPIInterfaceIMS,
                        virtual public IAlgorithmSPI,
                        virtual public IAlgorithmSPIInterface {
public:
    // bounding equity exposures. Not sure how this generalises to N
    DoubleArray          exposureUpperBound; // [numEquities], max exposure bound of equity basket to each equity
    double               overInvestBound;    // over-investment threshold
    double               underInvestBound;   // under-investment threshold
    double               equityExposureMin;
    double               equityExposureMax;
    DoubleArray          crashSizes;         // [numEquities]
    double               bondCrashSize;      // late arrival so optional default = 0.0

private:
    int                  numAlgAssets;          // taken from size of exposureUpperBound
    double               crashDenom;         // for 1 asset = 1.0, for 2 assets = c1-c0
    double               Cmax;
    double               Cmin;
    double               OIB;
    double               UIB;
    DoubleArray          pE;                 // [numEquities] - weights for actual rebalance
    DoubleArray          bondCrashAdj;       // [nbSteps] for adjustment = (Z ^ -bondCrashSize) - 1 $unregistered

public:
    static CClassConstSP const TYPE;
    SPIAlgorithmStd();// for reflection

    virtual void validatePop2Object();

    void init(const DoubleArray* bondPrices);

    virtual bool doRebalance(double UL,
                             double SL,
                             int    iStep);

    virtual const DoubleArray& rebalWeights(double crash);

    virtual double unbalCrash(const DoubleArray& nE,
                              const DoubleArray& E,
                              int                iStep) const;
    
    virtual double sustCrash(double buffer, int iStep) const;

    virtual double targetExp(double SE) const;
    
    virtual double feeAtMin(IFeesSPIConstSP fees) const;
    
    // incremental gap risk for this step
    virtual double gapRiskAtStep(double B, 
                                 double BF, 
                                 double equityComp,
                                 int    iStep) const;

    virtual void crossValidate(ICutoffSPISP cutoff, IFeesSPISP fees, 
                               bool isRainbowSPI);
    virtual int getNumRiskyAssets() const;

    // get the internal interface
    virtual IAlgorithmSPISP getAlgorithmSPI();
private:
    SPIAlgorithmStd(const SPIAlgorithmStd& rhs); // not implemented
    SPIAlgorithmStd& operator=(const SPIAlgorithmStd& rhs); // not implemented

    static IObject* defaultSPIAlgorithmStd();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
};
typedef smartPtr<SPIAlgorithmStd> SPIAlgorithmStdSP;

#define SPI_ALGORITHM_TYPE_STD   "Standard"

// the wrapper object which we needed before abstraction in IMS
class PRODUCTS_DLL SPIAlgorithmWrapper : public CObject,
                        virtual public IAlgorithmSPIInterface {
public:
    static CClassConstSP const TYPE;

    string          SPIAlgorithmType;
    SPIAlgorithmStdSP   algorithmStd;

public:

    virtual IAlgorithmSPISP getAlgorithmSPI();

    // validation
    void validatePop2Object();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
     
    // for reflection
    SPIAlgorithmWrapper();

    static IObject* defaultSPIAlgorithmWrapper();
};
typedef smartPtr<SPIAlgorithmWrapper> SPIAlgorithmWrapperSP;

DRLIB_END_NAMESPACE

#endif

