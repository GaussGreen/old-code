//----------------------------------------------------------------------------
//
//   Group       : Convertibles DR
//
//   Filename    : FD1FLNResettable.hpp
//
//   Description : One factor finite difference base class for log-normal processes
//
//   Author      : André Segger
//
//   Date        : 13 June 2003
//
//----------------------------------------------------------------------------

#ifndef FD1FLNRESETTABLE_HPP
#define FD1FLNRESETTABLE_HPP

#include "edginc/FD1FLN.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/VolSurface.hpp"


DRLIB_BEGIN_NAMESPACE

class TREE_DLL FD1FLNResettable: public FD1FLN {
public:
    static CClassConstSP const TYPE;
    friend class FD1FLNResettableHelper;

    virtual TimeMetricConstSP GetTimeMetric() const;

    FD1FLNResettable(const string& volType);

    /** Simple constructor */
    FD1FLNResettable();
    virtual ~FD1FLNResettable();

    /** Less simple constructor */
    FD1FLNResettable(int stepsPY, int stockSteps);
    FD1FLNResettable(CClassConstSP clazz);

    virtual int GetStepVol(int step, vector<double>& vol, const double* s, int start, int end);

    CVolProcessedBSConstSP getProcessedVol();

    int getNum3DSteps() const { return num3DSteps;}

protected:

    virtual void InitVol();
    /** calculate a term structure vol^2 */
    virtual void CalcV2Term(const DateTime& valDate, const DateTime& startDate,
                            const DateTime& matDate, CTermStructure& v2);
    /** set up variance array */
    virtual void PostSetup();
    
    FD1FLNResettable(const FD1FLNResettable& rhs);
    FD1FLNResettable& operator=(const FD1FLNResettable& rhs);

    // additional fields
    int     num3DSteps;
};

typedef smartPtr<FD1FLNResettable> FD1FLNResettableSP;

DRLIB_END_NAMESPACE
#endif
