#ifndef	SRMRATESDETDIFF_HPP
#define	SRMRATEDETDIFF_HPP

#include "edginc/DECLARE.hpp"
#include "edginc/SRMRatesDetUtil.hpp"
#include "edginc/QMCRatesDiffuse.hpp"


DRLIB_BEGIN_NAMESPACE

class SRMRatesDetermDiffuse : public QMCRatesDiffuse
{
public:
    SRMRatesDetermDiffuse(void);
    virtual ~SRMRatesDetermDiffuse(void);

    virtual void generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/);
    virtual void finalize(DateTimeArrayConstSP allDates);

    void setSRMRatesDetermDiffuse(SRMRatesDetermUtilSP _detUtilSP, const DateTime& _today);

    size_t getDiscYCIdx(void) 
    {
        if (discYCIdx < 0)
        {
            ASSERT(detUtilSP.get()); 
            discYCIdx =  registerYCFlavor(detUtilSP->getDiscYC());
        }
        return discYCIdx;

    }
    size_t getDiffYCIdx(void) 
    {
        if (diffYCIdx < 0)
        {
            ASSERT(detUtilSP.get()); 
            diffYCIdx =  registerYCFlavor(detUtilSP->getDiffYC());
        }
        return diffYCIdx;
    }

    QMCRatesUtilSP getRatesUtil(void) { ASSERT(detUtilSP.get()); return detUtilSP;}


    // performance critical	functions
    double getLnExpectedDiscFactor(size_t idx, FwdIdx i, FwdIdx	j) 
    {
        double ycf = getYCForward(idx, i, j);
        return ycf;
        
    }

    double getExpectedDiscFactor(size_t	idx, FwdIdx	i, FwdIdx j)
    {
        return exp(this->getLnExpectedDiscFactor(idx, i, j));        
    }


    // called after origSVol externally recalibrated via ICE 
    // not defined
    virtual void recalibrate(QMCRatesUtilSP thisSRMRatesUtilSP) { return; } 


private:

    SRMRatesDetermUtilSP detUtilSP;
    vector<double>  yearFrac;     // between sim dates
    int             diffYCIdx;
    int             discYCIdx;


};

DECLARE(SRMRatesDetermDiffuse);

DRLIB_END_NAMESPACE

#endif
