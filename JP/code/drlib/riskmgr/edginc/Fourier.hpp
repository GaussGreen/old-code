//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Fourier.hpp
//
//   Description : Interfaces for Fourier Processes & Products
//
//   Date        : 22 April 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_FOURIERPROC_HPP
#define EDR_FOURIERPROC_HPP
#include "edginc/Object.hpp"
#include "edginc/Class.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Complex.hpp"
#include "edginc/Calibrator.hpp"

DRLIB_BEGIN_NAMESPACE

class DateTime;
class TimeMetric;

/** Forward Starting Log Return Fourier Product Interface */
class RISKMGR_DLL FwdStFourierProductLogRtn { 
public:
    virtual const DateTime& getStartDate() const = 0;
    virtual ~FwdStFourierProductLogRtn(){};
};

typedef smartPtr<FwdStFourierProductLogRtn> FwdStFourierProductLogRtnSP;
typedef smartConstPtr<FwdStFourierProductLogRtn> FwdStFourierProductLogRtnConstSP;


/**ARNAUD**////////////////////////////////////////////////////////////////////
/** Forward Starting Expected Quadratic Variation Fourier Product Interface */ 
class RISKMGR_DLL FwdStFourierProductExpQuadVar { 
public:
    virtual const DateTime& getStartDate() const = 0;
    virtual ~FwdStFourierProductExpQuadVar(){};
};

typedef smartPtr<FwdStFourierProductExpQuadVar> FwdStFourierProductExpQuadVarSP;
typedef smartConstPtr<FwdStFourierProductExpQuadVar> FwdStFourierProductExpQuadVarConstSP;

/////////////////////////////////////////////////////////////////////////////////




/** Started Fourier Product Interface */
class RISKMGR_DLL StFourierProductLogRtn{
public:
    virtual ~StFourierProductLogRtn(){};
};

typedef smartPtr<StFourierProductLogRtn> StFourierProductLogRtnSP;
typedef smartConstPtr<StFourierProductLogRtn> StFourierProductLogRtnConstSP;


/** Forward Starting Fourier Product Integrated Variance Interface */
class RISKMGR_DLL FwdStFourierProductIntVar { 
public:
    virtual const DateTime& getStartDate() const = 0;
    virtual ~FwdStFourierProductIntVar(){};
};

typedef smartPtr<FwdStFourierProductIntVar> FwdStFourierProductIntVarSP;
typedef smartConstPtr<FwdStFourierProductIntVar> FwdStFourierProductIntVarConstSP;


/** Started Fourier Product Integrated Variance Interface */
class RISKMGR_DLL StFourierProductIntVar{
public:
    virtual ~StFourierProductIntVar(){};
};

typedef smartPtr<StFourierProductIntVar> StFourierProductIntVarSP;
typedef smartConstPtr<StFourierProductIntVar> StFourierProductIntVarConstSP;

/** Started Fourier Product Quadratic Variance Interface */
class RISKMGR_DLL StFourierProductQuadVar{
public:
    virtual ~StFourierProductQuadVar(){};
};

typedef smartPtr<StFourierProductQuadVar> StFourierProductQuadVarSP;
typedef smartConstPtr<StFourierProductQuadVar> StFourierProductQuadVarConstSP;

/** Forward Starting Fourier Product Quadratic Variance Interface */
class RISKMGR_DLL FwdStFourierProductQuadVar { 
public:
    virtual const DateTime& getStartDate() const = 0;
    virtual ~FwdStFourierProductQuadVar(){};
};

typedef smartPtr<FwdStFourierProductQuadVar> FwdStFourierProductQuadVarSP;
typedef smartConstPtr<FwdStFourierProductQuadVar> FwdStFourierProductQuadVarConstSP;

/** Fourier Process Object */
class FourierProduct;
class RISKMGR_DLL FourierProcess: public CObject{
public:
    static CClassConstSP const TYPE;

    FourierProcess(const CClassConstSP& clazz);

    /** Frequency range class. Contains the lower/upper bounds
        of the domain of definition. Needs to be inputted by user 
        when there is no known analytical solution to the bounds. */
    struct RISKMGR_DLL FrequencyRange: public CObject{
        static CClassConstSP const TYPE;
        friend class FourierProcessFrequencyRangeHelper;

        FrequencyRange();

        double lowerRealBound;
        double upperRealBound;
    };
    typedef smartConstPtr<FrequencyRange> FrequencyRangeConstSP;

    /** Create a MarketDataFetcher which will be used by the [Fourier] model
     *  for retrieving market data etc */
    virtual MarketDataFetcherSP marketDataFetcher() const = 0;

    /** * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this process.
     * See IModel::wantsRiskMapping(). */
    virtual IModel::WantsRiskMapping wantsRiskMapping() const = 0;

    /** Should do for now, even though it won't work if FourierProcess 
        encompasses more than 1 asset (as, then, there is more than 1 time 
        metric) */
    virtual const TimeMetric& getTimeMetric() const = 0;
    
    /** Returns the parameters of the process */
    virtual string getParameters() const;

    /** Give the process the chance to do some extra product specific 
        validation and to initialize its transient fields */
    virtual void validate(const FourierProduct* product) = 0;

protected:
    /** Extracts the parameters from an IAdjustable object */
    static string extractParameters(const Calibrator::IAdjustable* adjustable);

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<FourierProcess> FourierProcessSP;
typedef smartConstPtr<FourierProcess> FourierProcessConstSP;

/** Forward Starting Fourier Process Interface */
class RISKMGR_DLL FwdStFourierProcessLogRtn: virtual public IObject{
public:
    static CClassConstSP const TYPE;
    
    virtual Complex scalelessCumulant(const FwdStFourierProductLogRtn& product, 
                                      const Complex& z, 
                                      const DateTime& matdate) const = 0;
    virtual double lowerRealBound(const FwdStFourierProductLogRtn& product, 
                                  const DateTime& matdate) const = 0;
    virtual double upperRealBound(const FwdStFourierProductLogRtn& product, 
                                  const DateTime& matdate) const = 0;

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<FwdStFourierProcessLogRtn> FwdStFourierProcessLogRtnSP;
typedef smartConstPtr<FwdStFourierProcessLogRtn> FwdStFourierProcessLogRtnConstSP;



/**ARNAUD**/////////////////////////////////////////////////////////////
/** Forward Starting Fourier Process Interface */
class RISKMGR_DLL FwdStFourierProcessExpQuadVar: virtual public IObject{
public:
    static CClassConstSP const TYPE;
    
    virtual Complex cumulant(const FwdStFourierProductExpQuadVar& product, 
                                      const Complex& z, 
                                      const DateTime& matdate) const = 0;
    virtual double lowerRealBound(const FwdStFourierProductExpQuadVar& product, 
                                  const DateTime& matdate) const = 0;
    virtual double upperRealBound(const FwdStFourierProductExpQuadVar& product, 
                                  const DateTime& matdate) const = 0;

private:
    static void load(CClassSP& clazz);
};


typedef smartPtr<FwdStFourierProcessExpQuadVar> FwdStFourierProcessExpQuadVarSP;
typedef smartConstPtr<FwdStFourierProcessExpQuadVar> FwdStFourierProcessExpQuadVarConstSP;


////////////////////////////////////////////////////////////////////////////



/** Started Fourier Process Interface */
class RISKMGR_DLL StFourierProcessLogRtn: virtual public IObject{
public:
    static CClassConstSP const TYPE;
    
    virtual Complex scalelessCumulant(const StFourierProductLogRtn& product, 
                                      const Complex& z, 
                                      const DateTime& matdate) const = 0;
    virtual double lowerRealBound(const StFourierProductLogRtn& product, 
                                  const DateTime& matdate) const = 0;
    virtual double upperRealBound(const StFourierProductLogRtn& product, 
                                  const DateTime& matdate) const = 0;

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<StFourierProcessLogRtn> StFourierProcessLogRtnSP;
typedef smartConstPtr<StFourierProcessLogRtn> StFourierProcessLogRtnConstSP;

/** Started Fourier Process Interface */
class RISKMGR_DLL StFourierProcessIntVar: virtual public IObject{
public:
    static CClassConstSP const TYPE;
    
    virtual Complex cumulant(const StFourierProductIntVar& product, 
                             const Complex z, 
                             const DateTime& matdate) const = 0;
    virtual double lowerRealBound(const StFourierProductIntVar& product, 
                                  const DateTime& matdate) const = 0;
    virtual double upperRealBound(const StFourierProductIntVar& product, 
                                  const DateTime& matdate) const = 0;

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<StFourierProcessIntVar> StFourierProcessIntVarSP;
typedef smartConstPtr<StFourierProcessIntVar> StFourierProcessIntVarConstSP;


/** Forward Starting Fourier Process Interface */
class RISKMGR_DLL FwdStFourierProcessIntVar: virtual public IObject{
public:
    static CClassConstSP const TYPE;
    
    virtual Complex cumulant(const FwdStFourierProductIntVar& product, 
                             const Complex z, 
                             const DateTime& matdate) const = 0;
    virtual double lowerRealBound(const FwdStFourierProductIntVar& product, 
                                  const DateTime& matdate) const = 0;
    virtual double upperRealBound(const FwdStFourierProductIntVar& product, 
                                  const DateTime& matdate) const = 0;

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<FwdStFourierProcessIntVar> FwdStFourierProcessIntVarSP;
typedef smartConstPtr<FwdStFourierProcessIntVar> FwdStFourierProcessIntVarConstSP;


/** Started Fourier Process Quadratic Variance Interface */
class RISKMGR_DLL StFourierProcessQuadVar: virtual public IObject{
public:
    static CClassConstSP const TYPE;
    
    virtual Complex cumulant(const StFourierProductQuadVar& product, 
                             const Complex z, 
                             const DateTime& matdate) const = 0;
    virtual double lowerRealBound(const StFourierProductQuadVar& product, 
                                  const DateTime& matdate) const = 0;
    virtual double upperRealBound(const StFourierProductQuadVar& product, 
                                  const DateTime& matdate) const = 0;

    virtual double expectation(const StFourierProductQuadVar& product,
                               const DateTime& matdate) const = 0;

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<StFourierProcessQuadVar> StFourierProcessQuadVarSP;
typedef smartConstPtr<StFourierProcessQuadVar> StFourierProcessQuadVarConstSP;


/** Forward Starting Fourier Process Quadratic Variance Interface */
class RISKMGR_DLL FwdStFourierProcessQuadVar: virtual public IObject{
public:
    static CClassConstSP const TYPE;
    
    virtual Complex cumulant(const FwdStFourierProductQuadVar& product, 
                             const Complex z, 
                             const DateTime& matdate) const = 0;
    virtual double lowerRealBound(const FwdStFourierProductQuadVar& product, 
                                  const DateTime& matdate) const = 0;
    virtual double upperRealBound(const FwdStFourierProductQuadVar& product, 
                                  const DateTime& matdate) const = 0;

    virtual double expectation(const FwdStFourierProductQuadVar& product,
                               const DateTime& matdate) const = 0;


private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<FwdStFourierProcessQuadVar> FwdStFourierProcessQuadVarSP;
typedef smartConstPtr<FwdStFourierProcessQuadVar> FwdStFourierProcessQuadVarConstSP;



DRLIB_END_NAMESPACE
#endif
