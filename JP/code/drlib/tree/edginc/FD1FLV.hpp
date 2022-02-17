//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD1FLV.hpp
//
//   Description : One factor finite difference base class for local vol processes
//
//   Author      : Oleg Divinskiy
//
//   Date        : April 23 2002
//
//----------------------------------------------------------------------------

#ifndef FD1FLV_HPP
#define FD1FLV_HPP

#include "edginc/FD1F.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/VolProcessedDVF.hpp"

DRLIB_BEGIN_NAMESPACE

const double FD_MIN = 1.0e-10;

class TREE_DLL FD1FLV: public FD1F {
public:
    static CClassConstSP const TYPE;
    friend class FD1FLVHelper;

    virtual TimeMetricConstSP GetTimeMetric() const;

    FD1FLV(const string& volType);

    /** Simple constructor */
    FD1FLV();
    virtual ~FD1FLV();

    virtual int GetStepVol(int step, vector<double>& vol, const double* s_inp, int start, int end);

    /** Override default createMDF in order to set the right MDF */
    virtual MarketDataFetcherSP createMDF() const;

protected:
    virtual void InitVol();
    /** calculate a term structure vol^2 */
    virtual void CalcV2Term(const DateTime& valDate, const DateTime& startDate,
                            const DateTime& matDate, CTermStructure& v2);
    /** set up variance array */
    virtual void PostSetup();
    
private:
    FD1FLV(const FD1FLV &rhs);
    FD1FLV& operator=(const FD1FLV& rhs);

    // log-normal vol info 
    CVolProcessedDVFSP VolLV; // $unregistered
    string             volType;
};

typedef smartPtr<FD1FLV> FD1FLVSP;

DRLIB_END_NAMESPACE
#endif
