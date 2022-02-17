//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : FD1FLN.hpp
//
//   Description : One factor finite difference base class for log-normal processes
//
//   Author      : André Segger
//
//   Date        : April 04 2001
//
//----------------------------------------------------------------------------

#ifndef FD1FLN_HPP
#define FD1FLN_HPP

#include "edginc/FD1F.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/VolSurface.hpp"


DRLIB_BEGIN_NAMESPACE

class TREE_DLL FD1FLN: public FD1F {
public:
    static CClassConstSP const TYPE;
    friend class FD1FLNHelper;

    virtual TimeMetricConstSP GetTimeMetric() const;

    FD1FLN(const string& volType);

    /** Simple constructor */
    FD1FLN();
    virtual ~FD1FLN();

    /** Less simple constructor */
    FD1FLN(int stepsPY, int stockSteps);
    FD1FLN(CClassConstSP clazz);

    virtual int GetStepVol(int step, vector<double>& vol, const double* s, int start, int end);

    /** Override default createMDF in order to set the right MDF */
    virtual MarketDataFetcherSP createMDF() const;

    CVolProcessedBSConstSP getProcessedVol();

protected:

    virtual void InitVol();
    /** calculate a term structure vol^2 */
    virtual void CalcV2Term(const DateTime& valDate, const DateTime& startDate,
                            const DateTime& matDate, CTermStructure& v2);
    /** set up variance array */
    virtual void PostSetup();
    
    FD1FLN(const FD1FLN &rhs);
    FD1FLN& operator=(const FD1FLN& rhs);

    // log-normal vol info 
    CVolProcessedBSSP VolLN; // $unregistered
    string            volType;
};

typedef smartPtr<FD1FLN> FD1FLNSP;

DRLIB_END_NAMESPACE
#endif
