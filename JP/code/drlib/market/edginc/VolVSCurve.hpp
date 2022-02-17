//GAD, 15/02/2006

#ifndef EDR_VOLVSCURVE_HPP
#define EDR_VOLVSCURVE_HPP

#include "edginc/VolBaseParam.hpp"
#include "edginc/Asset.hpp"
#include "edginc/VolParam.hpp"
#include "edginc/Fourier.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/IDynamicsParameter.hpp"
#include "edginc/Complex.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Expiry.hpp"


DRLIB_BEGIN_NAMESPACE

/********************/
/*** class FVSCPC ***/

//piecewise constant forward variance swap curve

class MARKET_DLL FVSCPC: public CObject{
private:
    //empty constructor
    FVSCPC();

    //wrap around empty constructor
    static IObject* defaultCtor();

    //load method
    static void load(CClassSP& clazz);
	
public:
    //constructor
    FVSCPC(const DoubleArray&  inTenors,
           const DoubleArray&  inValues);

    //destructor
    virtual ~FVSCPC();

    //find the left tenor (piecewise constant curve is assumed to be right continuous)
    int tenorLocation(double tenor)const;

    //yield using piecewise constant interpolation
    double yield(double tenor) const;
    
    //find the floor fvscpcFloor on range [0,xMax]
    smartPtr<FVSCPC> floor(double   limit,
                           double   xMax) const;

    //registered fields
    DoubleArray tenors; //tenors used
    DoubleArray values; //values between tenors

    //TYPE, used for registration
    static CClassConstSP const TYPE;
};
typedef smartPtr<FVSCPC> FVSCPCSP;
typedef smartConstPtr<FVSCPC> FVSCPCConstSP;
typedef array<FVSCPCSP, FVSCPC> FVSCPCArray;

/*** end of class FVSCPC ***/
/***************************/


/********************/
/*** class VSCAff ***/

//piecewise affine variance swap curve 

class MARKET_DLL VSCAff: public CObject{
private:
    //empty constructor
    VSCAff();

    //wrap around empty constructor
    static IObject* defaultCtor();

    //find the left tenor
    int tenorLocation(double tenor) const;

    //load method
    static void load(CClassSP& clazz);

public:
    //constructor
    VSCAff(const DoubleArray&  inTenors,
           const DoubleArray&  inValues);

    //destructor
    virtual ~VSCAff();

    //yield using linear interpolation
    double yield(double tenor) const;

    //registered fields
    DoubleArray tenors; //tenors used
    DoubleArray values; //values at tenors

    //TYPE, used for registration
    static CClassConstSP const TYPE;
};
typedef smartPtr<VSCAff> VSCAffSP;
typedef smartConstPtr<VSCAff> VSCAffConstSP;
typedef array<VSCAffSP, VSCAff> VSCAffArray;

/*** end of class VSCAff ***/
/***************************/


/************************/
/*** class VolVSCurve ***/

//variance swap curve

class MARKET_DLL VolVSCurve: public VolBaseParam,
                  public virtual Calibrator::IAdjustable{
private:
    // Throws away past varswaps
    virtual bool sensShift(Theta* shift);

    //implied vol implemented in VolVSCurve to avoid dereferencing
	void ComputeImpVol(const CLatticeDouble&	strikes,
					   const DateTimeArray&		maturities,
				       CLatticeDouble&			impV) const;
	
	//matrix vol surface implemented in VolVSCurve to avoid dereferencing
	VolSurface* spotVolSurfaceFromStrikes(const CDoubleArray& strikes) const;

    //default values for parameters
    struct MARKET_DLL DefaultVal{
		static const double correlation;
		static const double volVol;
        static const double meanReversRate;
        static const double lambda;
        static const double crashMRR;
        static const double crashRate;
        static const double crashSize;
        static const double spotcrashRate;
        static const double spotcrashSize;
        static const double spotcrashUncertainty;
    };

    //ranges for parameters
    struct MARKET_DLL RangeDef{
		static const Range correlation;
		static const Range volVol;
        static const Range meanReversRate;
        static const Range lambda;
        static const Range crashMRR;
        static const Range crashRate;
        static const Range crashSize;
        static const Range spotcrashRate;
        static const Range spotcrashSize;
        static const Range spotcrashUncertainty;
    };

    //empty constructor
    VolVSCurve();

    ///default constructor
    static IObject* defaultCtor(){
	    return new VolVSCurve();
    }

    //invoked when class is loaded
    static void load(CClassSP& clazz);

protected:
	//market data
    virtual void getMarket(const IModel* model, const MarketData* market);
    
    //build vol and cache values
    //void buildCache();

    //VolBaseParam class
    class MARKET_DLL VSCurveVolParam:public CVolParam{
    private:
		//load method
        static void load(CClassSP& clazz){
            REGISTER(VSCurveVolParam, clazz);
            SUPERCLASS(CVolParam);
        }

	public:
		//TYPE used for registration
        static CClassConstSP const TYPE;

        //constructor
        VSCurveVolParam();

        //implied vol
        virtual void ComputeImpVol(const CVolBase*          vol,
                                   const CLatticeDouble&    strikes,
                                   const DateTimeArray&     maturities,
                                   CLatticeDouble&          impV) const;

        //matrix vol surface
        virtual VolSurface* spotVolSurfaceFromStrikes(const CVolBase*       vol,
                                                      const CDoubleArray&   strikes) const;
    };

public:
    friend class VSCurveVolParam;
    friend class CalcAlphaBeta;
    friend class FVSCPC;
    friend class VSCAff;

    //registered fields
    int n; //number of factors
    bool isFwdStartingVS; //whether or not the initial curve correponds to fwd starting or fwd starting VS
    
    DoubleArray correlation; //correlation array of size n
    DoubleArray volVol; //vol of vol array of size n
    DoubleArray meanReversRate; //mean reversion speed array of size n
    
    DoubleArray lambda; //weight array of size n,must sum up to 1

    // adding positive jumps in the variance swap curve
    // in effect, we are adding a compensated Poisson process to the diffusion of the variance swap curve
    // jumps are exponentially distributed with mean crashSize and arrive at a rate crashRate
    double crashMRR;
    double crashRate;
    double crashSize;

    // adding jumps in the spot
    double spotcrashRate;
    double spotcrashSize;
    double spotcrashUncertainty;

    ExpiryArray tenorCurve; //tenor array used to build the FVSCPC and VSCAFF declared as transient IS MANDATORY
    DoubleArray valueCurve; //value array used to build the FVSCPC and VSCAFF declared as transient inputed or not, IS OPTIONAL
                            //corresponds to the variance swap IMPLIED VOL, if input in the model it requires a transformation

	//transient fields
    FVSCPCSP fvscpc; //initial fvscpc built internally using StringArray tenorCurve !!! HERE WE USE THE INTEGRATED VARIANCE
    VSCAffSP vscaff; //initial vsc built internally using StringArray tenorCurve !!! HERE WE USE THE INTEGRATED VARIANCE
    DateTime baseDate; //not used
    FVSCPCSP fvscpcFloor;

    //TYPE for registration
    static CClassConstSP const TYPE;

    void validatePop2Object();

    //method that builds a CVolParam
    virtual CVolParam* createVolParam() const{
	    return new VSCurveVolParam();
    }

    //called by FourierProcessVSCurve
    Complex scalelessCumulant(const StFourierProcessLogRtn& process,
                              const StFourierProductLogRtn& product, 
                              const Complex&       		    z,
                              const DateTime&				matDate) const;
    
	Complex scalelessCumulant(const FwdStFourierProcessLogRtn&  process,
                              const FwdStFourierProductLogRtn&	product, 
                              const Complex&       					z,
                              const DateTime&					matDate) const;
    
	Complex cumulant(const StFourierProcessIntVar&  process,
                     const StFourierProductIntVar&	product, 
                     const Complex&    				z, 
                     const DateTime&				matDate) const;
    
	Complex cumulant(const FwdStFourierProcessIntVar&   process,
                     const FwdStFourierProductIntVar&	product, 
                     const Complex&    					z,
                     const DateTime&					matDate) const;

    Complex cumulant(const StFourierProcessQuadVar& process,
                     const StFourierProductQuadVar&	product, 
                     const Complex&    				z, 
                     const DateTime&				matDate) const;
    
	Complex cumulant(const FwdStFourierProcessQuadVar&  process,
                     const FwdStFourierProductQuadVar&	product, 
                     const Complex&   					z,
                     const DateTime&					matDate) const;
    
    double expectation(const StFourierProcessQuadVar& process,
                       const StFourierProductQuadVar& product, 
                       const DateTime& matDate) const;

    double expectation(const FwdStFourierProcessQuadVar& process,
                       const FwdStFourierProductQuadVar& product, 
                       const DateTime& matDate) const;

    Complex cumulant(const FwdStFourierProcessExpQuadVar& process,              
                     const FwdStFourierProductExpQuadVar& product, 
                     const Complex& z, 
                     const DateTime& matDate) const;

    //needed for IAdjustable interface, returns market data name for vol
    virtual string getName() const;

    //called after adjustments have been made
    virtual void update();

    // Called after adjustments have been made to fields (eg calibrator) 
    virtual void fieldsUpdated(const CFieldArray& fields);

    //access to model parameters
    enum EnumParam {N,
                    ISFWDSTARTINGVS,
                    CORRELATION,
                    VOL_VOL,
			        MEAN_REVERS_RATE,
			        LAMBDA,
                    CRASH_MRR,
                    CRASH_RATE,
                    CRASH_SIZE,
                    SPOT_CRASH_RATE,
                    SPOT_CRASH_SIZE,
                    SPOT_CRASH_UNCERTAINTY};

	DoubleArray getVSCurveParam(const EnumParam param) const;
};
typedef smartPtr<VolVSCurve> VolVSCurveSP;
typedef smartConstPtr<VolVSCurve> VolVSCurveConstSP;
typedef array<VolVSCurveSP, VolVSCurve>  VolVSCurveArray;

/*** end of class VolVSCurve ***/
/*******************************/


DRLIB_END_NAMESPACE
#endif
