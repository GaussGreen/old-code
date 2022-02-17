//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolStochGarf.hpp
//
//   Description : stock diffusion, vol diffusion,  commonStockVolJump
//
//   Date        : 20 Apr 2004
//
//
//----------------------------------------------------------------------------
#ifndef EDR_VOLStochGarf_HPP
#define EDR_VOLStochGarf_HPP

#include "edginc/VolBaseParam.hpp"
#include "edginc/Asset.hpp"
#include "edginc/VolParam.hpp"
#include "edginc/Fourier.hpp"
#include "edginc/Calibrator.hpp"
#include <list>
//#include "edginc/VolBaseParamSurface.hpp"
#include "edginc/VolatilityBS.hpp"
#include "edginc/VolatilityDVF.hpp"
//#include "edginc/VolProcessedBS.hpp"

DRLIB_BEGIN_NAMESPACE

/** A parameterised CVolBase (does not support BS and DVF views of the world for now) */
class MARKET_DLL VolStochGarf: public VolBaseParam,
                    //public CVolBaseParamSurface,
                    virtual public IVolProcessed,
                    virtual public IVolatilityBS,
                    virtual public IVolatilityDVF,                    
               public virtual Calibrator::IAdjustable {
public:
    typedef list<double>::iterator LIST_ITER;
    typedef enum{START_DATE, SAMPLE_DATE, DIFF_DATE, JUMP_DATE} STEP_TYPE;
    typedef list<STEP_TYPE>::iterator STEP_TYPE_ITER;
private:
    /** Our VolParam class. This class exists primarily to keep the
        ComputeImpVol and spotVolSurfaceFromStrikes methods out of the
        CVolBaseParamSurface interface */
    class MARKET_DLL StochGarfVolParam: public CVolParam{
    public:
        static CClassConstSP const TYPE;
        // constructor
        StochGarfVolParam(): CVolParam(TYPE){}

        virtual void ComputeImpVol(const CVolBase*          vol,
                                   const CLatticeDouble&    strikes,
                                   const DateTimeArray&     maturities,
                                   CLatticeDouble&          impV) const;

        /*** Creates a "Matrix Vol Surface" on provided Strikes  ***/
        virtual VolSurface* spotVolSurfaceFromStrikes(
            const CVolBase*       vol,
            const CDoubleArray&   strikes) const;

        virtual double computeLocVol(const CVolBase* vol, 
                                     const double strike,
                                     const DateTime& maturity) const;

    private:
        static void load(CClassSP& clazz){
            REGISTER(StochGarfVolParam, clazz);
            SUPERCLASS(CVolParam);
        }
    };

    /* Methods called by StochGarfVolParam. Implemented on VolStochGarf to avoid dereferencing */ 
    void ComputeImpVol(const CLatticeDouble&      strikes,
                       const DateTimeArray&       maturities,
                       CLatticeDouble&            impV) const;

    double computeLocVol(const double spot, const DateTime date) const;

	// return ATM vol level, by taking scale.  Scale should be -1, 0, 1 for 3 status
    double   getATMVol(const double scale);

    VolSurface* spotVolSurfaceFromStrikes(const CDoubleArray&   strikes) const;

    // registered fields
    double alpha;
    double beta;
    double meanReversRate;
    double volVol;
    double correlation;
    double atmVol;
    double strikeRef;
    ExpiryArraySP   ATMVolBM; // term structure diffusion vol bench mark labels $unregistered
    StringArray     BenchMarkStrg; // String Array for BenchMarks. $unregistered
    
    // transient fields (won't appear in dd interface)
    DateTime      baseDate;

    /* Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    static IObject* defaultCtor(){
        return new VolStochGarf();
    }

protected:
    /*** Build the parameterised vol and cache any values **/
    virtual void getMarket(const IModel* model, const MarketData* market);

public:
    friend class StochGarfVolParam;
    friend class FD2DStochGarf;
    friend class CalcAlpha; // local helper functions
    friend class FD2DSolverJumps;
    friend class VolProcessedStochGarf;


    static CClassConstSP const TYPE;

    VolStochGarf();

    // constructor for each local vol state
    //VolStochGarf(const VolStochGarf& vol, double atmVols);
    VolStochGarf(const VolStochGarf& vol, double scale);

    void validatePop2Object();

    /** method that builds a CVolParam. */
    virtual CVolParam* createVolParam() const{
        return new StochGarfVolParam();
    }

    IVolProcessed*	getProcessedVol(const CVolRequest* volRequest,const CAsset* asset) const;

    // return multiple processed vol with different volatility levels.
	IVolProcessed*  getProcessedVol(double scale) const;

    /** Needed for IAdjustable interface. Returns market data name for vol */
    virtual string getName() const;

};

typedef smartPtr<VolStochGarf> VolStochGarfSP;
typedef smartConstPtr<VolStochGarf> VolStochGarfConstSP;
typedef array<VolStochGarfSP, VolStochGarf>  VolStochGarfArray;

DRLIB_END_NAMESPACE
#endif
