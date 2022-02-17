//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : FD1FDDE.hpp
//
//   Description : One factor finite difference base class for DDE
//
//   Author      : Qing Hou
//
//   Date        : Nov 21, 2003
//
//----------------------------------------------------------------------------

#ifndef FD1F_DDE_HPP
#define FD1F_DDE_HPP

#include "edginc/FD1FGeneric.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/SpreadEquityFunc.hpp"
#include "edginc/DDEModule.hpp"

DRLIB_BEGIN_NAMESPACE

class Control;
class Results;

class TREE_DLL FD1FDDE: public FD1FGeneric, 
				virtual public IHaveNonEqSpread,
				virtual public Asset::ICanHaveDDE
 {
public:
    static CClassConstSP const TYPE;
    friend class FD1FDDEHelper;

    virtual TimeMetricConstSP GetTimeMetric() const;

    FD1FDDE(const string& volType);

    /** Simple constructor */
    FD1FDDE();
    virtual ~FD1FDDE();

    /** Less simple constructor */
    FD1FDDE(int stepsPY, int stockSteps);
    FD1FDDE(CClassConstSP clazz);

    virtual int GetStepVol(int step, vector<double>& vol, const double* s, int start, int end);

    CVolProcessedBSConstSP getProcessedVol();

	virtual FDTermStructureSP getDriftTerm(int step, const double* s, int start, int end, bool useFwdGris,
                                           const double &irPert, const double &divPert, const bool doEquityLayer);

	virtual FDTermStructureSP getDiffusionTerm(int step, const double* s, int start, int end, bool useFwdGris,
                                               const double &irPert, const double &divPert, const bool doEquityLayer);

	virtual FDTermStructureSP getCouponTerm(int step, const double* s, int start, int end, bool useFwdGris,
                                            const double &irPert, const double &divPert, const bool doEquityLayer);

	virtual FDTermStructureSP getDiscountTerm(int step, const double* s, int start, int end, bool useFwdGris,
                                              const double &irPert, const double &divPert, const bool doEquityLayer);

	const CleanSpreadCurve *getNonEqSpreads() const;

	bool isDDE() const { return true; }

    // functions called by product to initialize default payoff to get def payoff for each step
    // product need to add a price slice "defProbs", and pass this slice into the following calls
    // at maturity and before maturity
    // the defPayoff input is default payoff at maturity and at particular step date respectively
    void setDefaultProb(int step, double *defProbs, int start, int end);
    // adjust price by the default payoff
    void adjustPriceByDefPO(int nextStep, double nextDefPayoff, double *defProbs, double *price, int start, int end);

    /** Override default createMDF in order to set the right MDF */
    virtual MarketDataFetcherSP createMDF() const;

protected:

    virtual void InitVol();
    /** calculate a term structure vol^2 */
    virtual void CalcV2Term(const DateTime& valDate, const DateTime& startDate,
                            const DateTime& matDate, CTermStructure& v2);
    /** set up variance array */
    virtual void PostSetup();
    
	void recordOutputRequests(Control* control, CResults* results);

    FD1FDDE(const FD1FDDE& rhs);
    FD1FDDE& operator=(const FD1FDDE& rhs);

    // log-normal vol info 
    CVolProcessedBSSP VolLN; // $unregistered
    string            volType;

	// DDE
	smartPtr<DDEModule>			ddeModule;
	SpreadEquityFuncConstSP		csFunc;	// credit spread function $unregistered
	string						ddeType;
    DDECalibSP                  calibParams; // calibration method and param (CloseForm, MC, tree etc)

    // internal member to make sure the defProb call is step by step in correct order
    int                         defProbStep; // $unregistered
};

typedef smartPtr<FD1FDDE> FD1FDDESP;

DRLIB_END_NAMESPACE
#endif
