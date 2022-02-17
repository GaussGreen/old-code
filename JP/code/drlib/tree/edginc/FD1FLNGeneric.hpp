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

#ifndef FD1FLN_GENERIC_HPP
#define FD1FLN_GENERIC_HPP

#include "edginc/FD1FGeneric.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/VolSurface.hpp"


DRLIB_BEGIN_NAMESPACE

class TREE_DLL FD1FLNGeneric: public FD1FGeneric {
public:
    static CClassConstSP const TYPE;
    friend class FD1FLNGenericHelper;

    virtual TimeMetricConstSP GetTimeMetric() const;

    FD1FLNGeneric(const string& volType);

    /** Simple constructor */
    FD1FLNGeneric();
    virtual ~FD1FLNGeneric();

    /** Less simple constructor */
    FD1FLNGeneric(int stepsPY, int stockSteps);
    FD1FLNGeneric(CClassConstSP clazz);

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

    virtual bool doLambdaAdjust() 
    { 
        return false; 
    }

    /** Override default createMDF in order to set the right MDF */
    virtual MarketDataFetcherSP createMDF() const;

protected:

    virtual void InitVol();
    /** calculate a term structure vol^2 */
    virtual void CalcV2Term(const DateTime& valDate, const DateTime& startDate,
                            const DateTime& matDate, CTermStructure& v2);
    /** set up variance array */
    virtual void PostSetup();
    
    FD1FLNGeneric(const FD1FLNGeneric& rhs);
    FD1FLNGeneric& operator=(const FD1FLNGeneric& rhs);

    // log-normal vol info 
    CVolProcessedBSSP VolLN; // $unregistered
    string            volType;
};

typedef smartPtr<FD1FLNGeneric> FD1FLNGenericSP;

DRLIB_END_NAMESPACE
#endif
