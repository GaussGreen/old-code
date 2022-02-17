//----------------------------------------------------------------------------
//
//   Group       : QR Credit Hybrids
//
//   Filename    : QMCDetermRates.hpp
//
//   Description : deterministic IRs asset
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/QMCDetermRates.hpp"

DRLIB_BEGIN_NAMESPACE

///////////////////////////////////////////////////////////////////////////////
QMCDetermRates::QMCDetermRates()
    : yc(0)
{;}
///////////////////////////////////////////////////////////////////////////////
QMCDetermRates::~QMCDetermRates()
{;}
///////////////////////////////////////////////////////////////////////////////
void QMCDetermRates::setQMCDetermRates(DateTime today, IYieldCurveConstSP gYC)
{
    baseDate = today;
    yc = gYC;
    QLIB_VERIFY(!(!yc), "Error: yc points to NULL");
}
///////////////////////////////////////////////////////////////////////////////
void QMCDetermRates::generatePath(IQMCRNGManagerSP, QMCStrataConstSP)  // rngMgr
{}
///////////////////////////////////////////////////////////////////////////////
void QMCDetermRates::finalize(DateTimeArrayConstSP)  // gAllDates
{
    const static string method = "QMCDetermRates::finalize";
    try
    {
        QLIB_VERIFY(!(!getTimeLogic()), 
                             "Result of getTimeLogic() points to NULL");
        finalizeForGivenDates(getTimeLogic()->getDFDates(), cachedSpotIdxDFs);
        finalizeForGivenDates(getTimeLogic()->getFwdEDFDates(), 
                                                            cachedFwdIdxDFs);
    }
    catch(exception& e)
    {
        throw ModelException(e, method);
    }
}
///////////////////////////////////////////////////////////////////////////////
void QMCDetermRates::finalizeForGivenDates(const DateTimeArray & gDates,
                                                 DoubleArray   & gResult) 
{
    const static string method = "QMCDetermRates::finalizeForGivenDates";
    try
    {
        const int _qtyOfDates = gDates.size();
        gResult.resize(_qtyOfDates);
        int _q=0; for (; _q != _qtyOfDates; _q++)
        {
            gResult[_q] = yc->pv(gDates[_q]);
        }
    }
    catch(exception& e)
    {
        throw ModelException(e, method);
    }
}
///////////////////////////////////////////////////////////////////////////////
double  QMCDetermRates::getDiscFactor(SpotIdx measurementDateIdx)
{
    const static string method = "QMCDetermRates::getDiscFactor";
    try
    {
        // check the inputs ---------------------------------------------------
        QLIB_VERIFY(measurementDateIdx >= 0,"measurementDateIdx should be>=0");
        QLIB_VERIFY(measurementDateIdx < cachedSpotIdxDFs.size(), 
                    string("measurementDateIdx should be < cachedDFs.size() ")
                    + "= " + Format::toString(cachedSpotIdxDFs.size()));
        // output the result --------------------------------------------------
        return cachedSpotIdxDFs[measurementDateIdx];
    }
    catch(exception& e)
    {
        throw ModelException(e, method);
    }
}
///////////////////////////////////////////////////////////////////////////////
double QMCDetermRates::getExpectedDiscFactor(size_t,// ycIdx
                                             FwdIdx measurementDateIdx,
                                             FwdIdx futureDateIdx)
{
    const static string method = "QMCDetermRates::getExpectedDiscFactor";
    try
    {
        // check the inputs ---------------------------------------------------
        QLIB_VERIFY(measurementDateIdx >= 0,"measurementDateIdx should be>=0");
        QLIB_VERIFY(measurementDateIdx < futureDateIdx,
                    "measurementDateIdx should be < futureDateIdx");
        QLIB_VERIFY(futureDateIdx < cachedFwdIdxDFs.size(), 
            string("futureDateIdx should be < cachedDFs.size() ")
            + "= " + Format::toString(cachedFwdIdxDFs.size()));
        // output the result --------------------------------------------------
        return (cachedFwdIdxDFs[futureDateIdx]/cachedFwdIdxDFs[measurementDateIdx]);
    }
    catch(exception& e)
    {
        throw ModelException(e, method);
    }
}
///////////////////////////////////////////////////////////////////////////////
double QMCDetermRates::getLnExpectedDiscFactor(size_t ycIdx,
                                               FwdIdx measurementDateIdx,
                                               FwdIdx futureDateIdx)
{
    return log(this->getExpectedDiscFactor(ycIdx, 
                                           measurementDateIdx, 
                                           futureDateIdx));
}
///////////////////////////////////////////////////////////////////////////////
double QMCDetermRates::getOriginalLnDiscFactor(SpotIdx futureDateIdx)
{
    return log(this->getDiscFactor(futureDateIdx));
}
///////////////////////////////////////////////////////////////////////////////
double QMCDetermRates::getOriginalLnExpectedDiscFactor(size_t ycIdx,
                                               FwdIdx measurementDateIdx,
                                               FwdIdx futureDateIdx)
{
    return this->getLnExpectedDiscFactor(ycIdx, 
                                         measurementDateIdx, 
                                         futureDateIdx);
}
///////////////////////////////////////////////////////////////////////////////
size_t QMCDetermRates::registerYCFlavor(IYieldCurveConstSP)
{
    return 0;
}
///////////////////////////////////////////////////////////////////////////////
size_t QMCDetermRates::getDiscYCIdx(void)
{
    return 0;
}
///////////////////////////////////////////////////////////////////////////////
size_t QMCDetermRates::getDiffYCIdx(void)
{
    return 0;
}
///////////////////////////////////////////////////////////////////////////////
vector<double>& QMCDetermRates::getDomLnMONEY()
{
    const static string method 
        = "QMCDetermRates::getDomLnMONEY";
    throw ModelException("Makes no sense for now YVXXX", method);
}
///////////////////////////////////////////////////////////////////////////////
const vector<double>* QMCDetermRates::getSigmaFX()
{
    return NULL;
    const static string method 
        = "QMCDetermRates::getSigmaFX";
    throw ModelException("Makes no sense for now YVXXX", method);
}
///////////////////////////////////////////////////////////////////////////////
DateTime QMCDetermRates::getBaseDate() 
{
    return baseDate;
    // return getTimeLogic()->getBaseDate();
    const static string method 
        = "QMCDetermRates::getBaseDate";
    throw ModelException("Under construction", method);
}


///////////////////////////////////////////////////////////////////////////////
const vector<double>* QMCDetermRates::getSpotFXFullPath()
{
    return NULL;
    const static string method 
        = "QMCDetermRates::getSpotFXFullPath";
    throw ModelException("Makes no sense for now JDXXX", method);
}

DRLIB_END_NAMESPACE
