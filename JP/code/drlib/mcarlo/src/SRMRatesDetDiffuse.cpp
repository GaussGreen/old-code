#include "edginc/config.hpp"
#include "edginc/QMCRatesDiffuse.hpp"
#include "edginc/SRMRatesDetDiffuse.hpp"
#include "edginc/SRMUtil.hpp"

DRLIB_BEGIN_NAMESPACE

SRMRatesDetermDiffuse::SRMRatesDetermDiffuse(void) : discYCIdx(-1),
                                                diffYCIdx(-1),
                                                detUtilSP()
{
}

SRMRatesDetermDiffuse::~SRMRatesDetermDiffuse(void)
{
}

void SRMRatesDetermDiffuse::setSRMRatesDetermDiffuse(SRMRatesDetermUtilSP _detUtilSP, const DateTime& _today)
{
    detUtilSP = _detUtilSP;
    today = _today;
}

void SRMRatesDetermDiffuse::finalize(DateTimeArrayConstSP allDatesSP)
{
    static const string method("SRMRatesDetermDiffuse::finalize");

    const size_t numSimDates = SRMUtil::getNumSimDates(today, *allDatesSP);

    DateTimeArray dfRequestedDates = getTimeLogic()->getDFDates();
    DateTimeArray edfRequestedDates = getTimeLogic()->getReqEDFDates();
    DateTimeArray  edfForwardDates = getTimeLogic()->getFwdEDFDates();

    SRM_VERIFY(DateTime::isSubset(edfForwardDates, edfRequestedDates),
        "Internal error in expected discount factor dates",
        method);

    SRM_VERIFY(detUtilSP.get() != NULL,
        "Internal error : ratesHJMUtil void",
        method);

    DateTimeArray   diffusionDates = SRMUtil::calcDiffusionDates(today, *allDatesSP);

    const size_t Nall =  SRMUtil::getNumSimDates(today, *allDatesSP); // rename later
    const size_t Nedf =  edfRequestedDates.size(); // rename later
    const size_t Ndf  =  dfRequestedDates.size();

    yearFrac.assign(Nall-1, 0.0);
    
    for (size_t i = 0; i < yearFrac.size(); i++){
        yearFrac[i] = diffusionDates[i].yearFrac(diffusionDates[i+1]);
    }


    // cached variables
    // TODO getSubIndexes should throw an error if input date is not found in allDates
    SRM_VERIFY(DateTime::isSubset((*allDatesSP), dfRequestedDates),
        "Internal error : dfRequestDates not in AllDates", 
        method); 

    dfIndexes = DateTime::getIndexes((*allDatesSP),  dfRequestedDates); // indexes for when we save discount factors [Ndf]

    SRM_VERIFY(DateTime::isSubset((*allDatesSP), edfRequestedDates),
        "Internal error : edfRequestDates not in AllDates", 
        method); 

    expDFIndexes = DateTime::getIndexes((*allDatesSP), edfRequestedDates); // indexes for when we compute expected  [Nedf]

    if (!fx && expDFIndexes.empty() && dfIndexes.empty())
    {
        throw ModelException(method, "Internal error - no "
            "results required from diffused path");
    }

    todayIdx = today.find((*allDatesSP));

    todayIndex = today.findUpper(dfRequestedDates);
    if (todayIndex >= 0 && todayIndex != dfRequestedDates.size() && dfRequestedDates[todayIndex] == today)
        ++todayIndex; // skip today as diffusion results are always in futureDates

    SRM_VERIFY( todayIndex == dfRequestedDates.size() || dfRequestedDates[todayIndex] > today,
        "Internal error : problem in dfRequestedDates", 
        method); 

    // now make life easier - add request for index off the end
    size_t errIndex = (*allDatesSP).size()+numSimDates+1;
    dfIndexes.push_back(errIndex);
    expDFIndexes.push_back(errIndex);
       

    
    SRM_VERIFY(numSimDates == diffusionDates.size(),
        "Internal error : problem in timeline",
        method);

    getDiffYCIdx();
    getDiscYCIdx();

    //Deterministic rates: discount factor gives r_t
    for(size_t i=0; i < ycForwardsDB.size(); ++i)
    {
        IYieldCurveConstSP  yc = ycForwardsDB[i].first;
        vector<double> &   fwd = ycForwardsDB[i].second;

        const DateTimeArray *dta = allDatesSP.get();
        detUtilSP->computeLogDiscFactor(edfForwardDates, fwd, yc); // initialize fwd between t0 and T
    }
    detUtilSP->computeLogDiscFactor(logDiscFactor); // [Nall]

    df.resize(Ndf, 0.0);
    int dfDatePos = todayIndex; // position in dfIndexes/df.   
    int dfDateIdx = dfIndexes[dfDatePos];
    if (dfDateIdx == 0)
    {
        // save the value
        double saveVal = 1.0; //exp(-total); 
        df[dfDatePos] = saveVal;
        dfDatePos++;
        dfDateIdx = dfIndexes[dfDatePos];
    }

    //'domLnMoney' is log of money market account, offset by 1: add an extra index so we can use the same indices for the discount factors.
    domLnMONEY = vector<double>(numSimDates-1);
    double total = 0.0;

    for (size_t i = 0;i<domLnMONEY.size();i++)
    {
        double inc = -logDiscFactor[i];
        total += inc;
        domLnMONEY[i] = total;
        if (i + 1 + todayIdx == dfDateIdx)
        {
            // save the value
            double saveVal = exp(-total); 
            df[dfDatePos] = saveVal;
            dfDatePos++;
            dfDateIdx = dfIndexes[dfDatePos];
        }
    }  
}

void SRMRatesDetermDiffuse::generatePath(IQMCRNGManagerSP rngMgr, QMCStrataConstSP /*strata*/)
{
    /*
     * Do not do much: domLnMONEY has already been set in finalize ('deterministically').
     * We need to go through the motions, though, in order to diffuse the FX.
     * Indeed, currently, the FX doesn't get diffused in generatePath. 
     * This is because the rates needs a CUPS adjust, and the FX needs the rates, 
     * so they must be diffused 'at the same time'.
     */

    int numDates = logDiscFactor.size();
    if (fx.get()){
        double sigmaFX = fx->begin(rngMgr);
        for (int i = 0; i < numDates; ++i){
            fx->moveToNextDate(domLnMONEY[i], i); // get new sigmaFX
        }
    }

    if (fx.get()){
        fx->end();
    }
}

DRLIB_END_NAMESPACE


