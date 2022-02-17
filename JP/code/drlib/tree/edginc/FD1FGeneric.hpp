//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : FD1FGeneric.hpp
//
//   Description : One factor generic finite difference base class
//
//   Author      : André Segger
//
//   Date        : 15 October 2003
//
//----------------------------------------------------------------------------

#ifndef FD1F_GENERIC_HPP
#define FD1F_GENERIC_HPP

#include "edginc/Model.hpp"
#include "edginc/Asset.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/TimeLine.hpp"
#include "edginc/FDEngineGeneric.hpp"
#include "edginc/Model1F.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL FD1FGeneric: public CModel,
                   public Model1F {
public:
    static CClassConstSP const TYPE;
    friend class FD1FGenericHelper;

    FDEngine1FGeneric fdEngine; // $unregistered

    // time line
    vector<double> stockMaxSeg; // $unregistered
    vector<double> stockMinSeg; // $unregistered
    vector<double> strikeSeg; // $unregistered

    /** fd initialisation */
    virtual void Setup(const DateTime& valDate, const DateTimeArray& segDates, 
               const vector<int>& density, const DateTimeArray* critDates,
               double minGap, bool equalTime, int numOfPrice);
    
    class TREE_DLL IProduct: virtual public Model1F::Product1F,
                    public FDPayoff1FGeneric {
    public:

        FD1FGeneric *genericFDModel;

  //      virtual double getVolFD(double stockPrice, int step){
            // proof of concept function. i.e. it's not robust. Should call fdModel->getStepVol instead.
  //          return sqrt((fdModel->Variance[step] - fdModel->Variance[step-1])/fdModel->TimePts.TradeYrFrac[step]);
 //       }

        virtual void getVolFD(int step, vector<double>& vol, const double* s, int start, int end){
	    genericFDModel->GetStepVol(step,vol,s,start,end);  
        }

        // currently returns a single coupon payment - could be made more generic to allow for spot dependent coupons
        // the coupon is assumed to be the coupon itself, ie. without the intensity (like default rate)
        virtual double getCoupon(int step, const double* s, int start, int end) = 0;

        virtual bool hasEquityLayer() = 0;

        IProduct() : genericFDModel(0){}

        virtual ~IProduct() {};
        
//        virtual double scalePremium(const double& fairValue) {return fairValue;};
        
        virtual void InitFD(Control*    control) {
            Init(control);
        };
        
        virtual void InitProdFD() {
            InitProd();
        };

        virtual void postProcessFD(const double *price, const int numLayers) {
            // do nothing by default
        };

        /** allow per instrument default grid size setting */
        virtual void InitGridSize() { 
            if (genericFDModel->stockStepsToUse == -1) { // grid hasn't been set yet
                if (genericFDModel->stockSteps == -1) {
                    genericFDModel->stockStepsToUse = 100;
                } else {
                    genericFDModel->stockStepsToUse = genericFDModel->stockSteps;
                }
                if (genericFDModel->stepsPerYear == -1 ) {
                    genericFDModel->stepsPerYearToUse = 100;
                } else {
                    genericFDModel->stepsPerYearToUse = genericFDModel->stepsPerYear;
                }
            }
        }

        virtual YieldCurveConstSP GetRiskFreeCurve() {
            return GetDiscCurveRef();
        }

    };
    
    
    /** interface that the instrument must implement */
    class TREE_DLL IIntoProduct{
    public:
        friend class FD1FGenericHelper;
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(FD1FGeneric* model) const = 0;
    };

    FD1FGeneric(const string& volType);

    /** Simple constructor - for derived classes */
    FD1FGeneric(CClassConstSP clazz);
    virtual ~FD1FGeneric();
 
    /** clean up */
    virtual void Clear();
  
    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);
    
    virtual int GetStepVol(int step, vector<double>& vol, const double* s, int start, int end) = 0;

    virtual void postProcessFD(const double* price, const int numPriceArrays) {
        product->postProcessFD(price, numPriceArrays);
    }

    virtual void preProcessGrid(const double* s, int start, int end);

	virtual FDTermStructureSP getDriftTerm(int step, const double* s, int start, int end, bool useFwdGris,
                                           const double &irPert, const double &divPert, const bool doEquityLayer) = 0;

	virtual FDTermStructureSP getDiffusionTerm(int step, const double* s, int start, int end, bool useFwdGris,
                                               const double &irPert, const double &divPert, const bool doEquityLayer) = 0;

	virtual FDTermStructureSP getCouponTerm(int step, const double* s, int start, int end, bool useFwdGris,
                                            const double &irPert, const double &divPert, const bool doEquityLayer) = 0;

	virtual FDTermStructureSP getDiscountTerm(int step, const double* s, int start, int end, bool useFwdGris,
                                              const double &irPert, const double &divPert, const bool doEquityLayer) = 0;

    virtual bool doLambdaAdjust() 
    { 
        return false; 
    }

    virtual double getLambda() 
    { 
        return 0.0; 
    }

    virtual double getLambdaAdjustedSpot(const double spot) 
    { 
        return spot; 
    }

    virtual bool isE2C() 
    { 
        return false; 
    }

    virtual bool hasEquityLayer() 
    {
        return true; 
    }

    virtual void mapToEquitySpace(double* stockArray, int numnStockSteps)
    {
        // do nothing
    }

    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingIrrelevant, because FD1FGeneric models are not
     * parametric (latent, non-observable).  See IModel::wantsRiskMapping().
     */

    virtual IModel::WantsRiskMapping wantsRiskMapping() const;

    // input variable made public for penultimate step smoothing
    CDoubleArray Variance; // $unregistered
    bool   useFwdGrid;
    int    stockSteps;
    int    stockStepsToUse; 
    int    stepsPerYear;
    int    stepsPerYearToUse;
    int    numPriceArrays; // $unregistered
    bool   DEBUG_UseCtrlVar;
    bool   DEBUG_SameGridDelta;
    bool   DEBUG_SameGridVegaRho; // true => same grid for rho, vega, etc. 

    // const refs to asset and discount curve
    YieldCurveConstSP   DiscountCurve; // $unregistered
    CAssetConstSP       Underlier; // $unregistered
    YieldCurveConstSP   riskFreeCurve; // $unregistered

    CDoubleArray inForwards; // $unregistered
    vector<double> inPVs; // $unregistered
    vector<double> inDivy; // $unregistered
    vector<double> inIr; // $unregistered
    vector<int> inSegEnd; // $unregistered
    vector<double> inSegMax; // $unregistered
    vector<double> inSegMin; // $unregistered
    vector<double> inSegStrike; // $unregistered
    vector<double> times; // $unregistered


    CDoubleArray divPert; // $unregistered
    CDoubleArray irPert; // $unregistered

protected:
    FD1FGeneric::IProduct *product; // $unregistered

    virtual void InitVol() = 0;
    /** calculate a term structure vol^2 */
    virtual void CalcV2Term(const DateTime& valDate, const DateTime& startDate,
                            const DateTime& matDate, CTermStructure& v2) = 0;
    /** set up variance array */
    virtual void PostSetup() = 0;
    
    /** calculate the stuff which specifies the grid */
    virtual void SetGridData();

    virtual bool SameGrid(const CControl* control);

	// possible output request for model itself. Used by DDE
	virtual void recordOutputRequests(CControl* control, CResults* results){}; 

    // input variables

    int    gridType;


    int    varMethod;

private:
    FD1FGeneric(const FD1FGeneric &rhs);
    FD1FGeneric& operator=(const FD1FGeneric& rhs);
};


typedef smartPtr<FD1FGeneric> FD1FGenericSP;

DRLIB_END_NAMESPACE
#endif
