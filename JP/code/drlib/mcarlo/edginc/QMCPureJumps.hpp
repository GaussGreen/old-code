//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : QMCPureJumps.hpp
//
//   Description : Generator of jumps and retrieval of generated data
//

//
//----------------------------------------------------------------------------

#ifndef QLIB_QMCPUREJUMPS_HPP
#define QLIB_QMCPUREJUMPS_HPP

#include "edginc/IQMCPureJumps.hpp"
#include "edginc/QMCHelperBoundedDiffusion.hpp"
#include "edginc/QMCStrata.hpp"

DRLIB_BEGIN_NAMESPACE


class MCARLO_DLL QMCPureJumps : public IQMCPureJumps
{
public:

    // Common hookups for the engine

    /** finalize the timelines, allocate necessary memory */
    virtual void finalize(DateTimeArrayConstSP allDates);

    virtual double getStrataProbability(QMCStrataConstSP strata);

    /** after the path was generated, accessing the generated list of jump dates */
    virtual size_t      getTotalNumberOfJumps();
    virtual size_t      getNumberOfJumps(size_t sourceIdx);
    virtual DateTime    getJumpDateByIdx(size_t sourceIdx, size_t jumpIdx);
    virtual double      getJumpDateByIdxAsDouble(size_t sourceIdx, size_t jumpIdx);
    virtual double      getJumpIntensity(size_t sourceIdx)
    { return intensities.at(sourceIdx); }

    /** generate path across all dates. */
    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/);

    /** getting the simulation start date */
    virtual DateTime getBaseDate() { return today; }

    /** getting information about asset specific diffusion bounds */
    virtual QMCHelperBoundedDiffusion * getDiffusionBound() { return &diffusionBound;}

    QMCPureJumps() : 
        today(), intensities(), low(-999), high(-999),
        probaLow(-999.0), probaZeroLow(-999.0), 
        probaZeroHigh(-999.0), totalIntensity(-999.0), 
        lastDate(-999.0) {}
    ~QMCPureJumps() {}

    void setQMCPureJumps(
        const DateTime&       _baseDate,
        const vector<double>& _intensities);

private:
    /** given a (random) number in [0,1], give the associated Poisson (random) number
        If condionalOnAtLeastOneJump, does so using the Poisson law conditional on the Poisson rv > 0 */
    static size_t getPoissonSampleStratified(double intensity, double unif, 
                                      int _low, 
                                      double _probaLow, 
                                      double _probaInLowHigh);
    static size_t getPoissonSample(double intensity, double unif);
    static size_t getBinomialSample(double p, size_t N, double unif);

    void  modifyLowHigh(int _low, int _high);

    DateTime         today;
    vector<double>   intensities;
    double             lastDate;
    QMCHelperBoundedDiffusion diffusionBound;

    // transient fields:
    vector<size_t>   offsets;
    vector<double>   jumpDates;

    double totalIntensity;

    int low, high;
    double probaLow, probaZeroLow, probaZeroHigh;

};

/** a hookup for across-the-assets diffusion collector */
DECLARE(QMCPureJumps);


DRLIB_END_NAMESPACE
#endif //QLIB_QMCPUREJUMPS_HPP

