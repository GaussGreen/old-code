//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : FD1FE2C.hpp
//
//   Description : One factor finite difference base class for E2C dynamic
//
//   Author      : André Segger
//
//   Date        : 16 October 2003
//
//----------------------------------------------------------------------------

#ifndef FD1FE2C_HPP
#define FD1FE2C_HPP

#include "edginc/FD1FGeneric.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/CDSHelper.hpp"
#include "edginc/VolSurface.hpp"


DRLIB_BEGIN_NAMESPACE

class TREE_DLL FD1FE2C: public FD1FGeneric {
public:
    static CClassConstSP const TYPE;
    friend class FD1FE2CHelper;

    virtual TimeMetricConstSP GetTimeMetric() const;

    FD1FE2C(const string& volType);

    /** Simple constructor */
    FD1FE2C();
    virtual ~FD1FE2C();

    /** Less simple constructor */
    FD1FE2C(int stepsPY, int stockSteps);
    FD1FE2C(CClassConstSP clazz);

    virtual int GetStepVol(int step, vector<double>& vol, const double* s, int start, int end);

    virtual int GetAssetStepVol(int step, vector<double>& vol, const double*, int, int end);

    CVolProcessedBSConstSP getProcessedVol();

    virtual void SetGridData();

    virtual void preProcessGrid(const double* s, int start, int end);

	virtual FDTermStructureSP getDriftTerm(int step, const double* s, int start, int end, bool useFwdGris,
                                           const double &irPert, const double &divPert, const bool doEquityLayer);

	virtual FDTermStructureSP getDiffusionTerm(int step, const double* s, int start, int end, bool useFwdGris,
                                               const double &irPert, const double &divPert, const bool doEquityLayer);

	virtual FDTermStructureSP getCouponTerm(int step, const double* s, int start, int end, bool useFwdGris,
                                            const double &irPert, const double &divPert, const bool doEquityLayer);

	virtual FDTermStructureSP getDiscountTerm(int step, const double* s, int start, int end, bool useFwdGris,
                                              const double &irPert, const double &divPert, const bool doEquityLayer);

    virtual bool doLambdaAdjust();

    virtual double getLambda();

    virtual double getLambdaAdjustedSpot(const double spot);

    virtual bool isE2C();

    virtual bool hasEquityLayer();

    virtual void mapToEquitySpace(double* stockArray, int numnStockSteps);

    double getFXRate();

    /** Override default createMDF in order to set the right MDF */
    virtual MarketDataFetcherSP createMDF() const;

    // additional grid information for equity layers
    vector<double> equityMinSeg; // $unregistered
    vector<double> equityMaxSeg; // $unregistered
    vector<double> equityStrikeSeg; // $unregistered

    vector<double> equityInSegMin; // $unregistered
    vector<double> equityInSegMax; // $unregistered
    vector<double> equityInSegStrike; // $unregistered

protected:

    virtual void InitVol();
    /** calculate a term structure vol^2 */
    virtual void CalcV2Term(const DateTime& valDate, const DateTime& startDate,
                            const DateTime& matDate, CTermStructure& v2);
    /** set up variance array */
    virtual void PostSetup();
    
    FD1FE2C(const FD1FE2C& rhs);
    FD1FE2C& operator=(const FD1FE2C& rhs);

    // log-normal vol info 
    CVolProcessedBSSP                   VolLN; // $unregistered
    CVolProcessedBSSP                   VolLNEquity; // $unregistered
    string                              volType;
    CDSHelper::CParSpreadDefaultRatesSP liquiditySpreadRates; // $unregistered
	CDSHelper::CFirmAssetDefaultRatesSP equitySpreadRates; // $unregistered


    // additional information for equity layer
    CDoubleArray equityVariance; // $unregistered

    // additional information for equity layer
    CDoubleArray                        inFwdEQ; // $unregistered
    vector<double>                      inDivyEQ; // $unregistered
    vector<double>                      inIrEQ; // $unregistered
    vector<double>                      inPVRiskFree; // $unregistered
    
    // additional information for asset drift
    bool                                useDivsForAssetDrift;
    double                              additionalAssetDrift;

    // grow equity at risky rate
    bool                                riskyEquityGrowth;
    CDoubleMatrix                       equitySpreads; // $unregistered

    // transient field - default barrier to speed up calculations
    double                              defaultBarrier; // $unregistered

};

typedef smartPtr<FD1FE2C> FD1FE2CSP;

DRLIB_END_NAMESPACE
#endif
