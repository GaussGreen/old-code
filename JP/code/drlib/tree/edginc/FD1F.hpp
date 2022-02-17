//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : FD1F.hpp
//
//   Description : One factor finite difference base class
//
//   Author      : Tycho von Rosenvinge
//
//   Date        : November 20, 2001
//
//----------------------------------------------------------------------------

#ifndef FD1F_HPP
#define FD1F_HPP
#include "edginc/Model.hpp"
#include "edginc/Asset.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/TimeLine.hpp"
#include "edginc/FDEngine.hpp"
#include "edginc/FDEngineGeneric.hpp"
#include "edginc/Model1F.hpp"

DRLIB_BEGIN_NAMESPACE




class TREE_DLL FD1F: public CModel,
            public Model1F {
public:
    static CClassConstSP const TYPE;
    friend class FD1FHelper;

    FDEngine1F fdEngine; // $unregistered

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

        FD1F *fdModel;

  //      virtual double getVolFD(double stockPrice, int step){
            // proof of concept function. i.e. it's not robust. Should call fdModel->getStepVol instead.
  //          return sqrt((fdModel->Variance[step] - fdModel->Variance[step-1])/fdModel->TimePts.TradeYrFrac[step]);
 //       }

        virtual void getVolFD(int step, vector<double>& vol, const double* s, int start, int end){
	    fdModel->GetStepVol(step,vol,s,start,end);  
        }

        IProduct() : fdModel(0){}

        virtual ~IProduct() {};
        
//        virtual double scalePremium(const double& fairValue) {return fairValue;};
        
        virtual void InitFD(Control*    control) {
            Init(control);
        };
        
        virtual void InitProdFD() {
            InitProd();
        };

        /** allow per instrument default grid size setting */
        virtual void InitGridSize() { 
            if (fdModel->stockStepsToUse == -1) { // grid hasn't been set yet
                if (fdModel->stockSteps == -1) {
                    fdModel->stockStepsToUse = 100;
                } else {
                    fdModel->stockStepsToUse = fdModel->stockSteps;
                }
                if (fdModel->stepsPerYear == -1 ) {
                    fdModel->stepsPerYearToUse = 100;
                } else {
                    fdModel->stepsPerYearToUse = fdModel->stepsPerYear;
                }
            }
        }
    };
    
    
    /** interface that the instrument must implement */
    class TREE_DLL IIntoProduct{
    public:
        friend class FD1FHelper;
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(FD1F* model) const = 0;
    };

    FD1F(const string& volType);

    /** Simple constructor - for derived classes */
    FD1F(CClassConstSP clazz);
    virtual ~FD1F();

    /** clean up */
    virtual void Clear();
  
    /** calculate single price and store result in CResult */
    virtual void Price(CInstrument*  instrument, 
                       CControl*     control, 
                       CResults*     results);
    
    virtual int GetStepVol(int step, vector<double>& vol, const double* s, int start, int end) = 0;
    
    /**
     * Whether to enable RiskMapping when computing sensitivities for
     * instruments priced using this model
     *
     * Returns riskMappingIrrelevant, because FD1F models are not parametric
     * (latent, non-observable).  See IModel::wantsRiskMapping().
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
    FD1F::IProduct *product; // $unregistered

    virtual void InitVol() = 0;
    /** calculate a term structure vol^2 */
    virtual void CalcV2Term(const DateTime& valDate, const DateTime& startDate,
                            const DateTime& matDate, CTermStructure& v2) = 0;
    /** set up variance array */
    virtual void PostSetup() = 0;
    
    /** calculate the stuff which specifies the grid */
    virtual void SetGridData();

    virtual bool SameGrid(const CControl* control);

    // input variables

    int    gridType;


    int    varMethod;

private:
    FD1F(const FD1F &rhs);
    FD1F& operator=(const FD1F& rhs);
};


typedef smartPtr<FD1F> FD1FSP;

DRLIB_END_NAMESPACE
#endif
