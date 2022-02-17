//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : VolCalibInterp.cpp
//
//   Description : Simulates a VolCalibInterp for testing purposes!
//
//   Author      : Anwar E Sidat
//
//   Date        : 07-Oct-2006
//
//----------------------------------------------------------------------------

#include <stdio.h>
#include <math.h>
#include <vector>
#include "edginc/config.hpp"
#include "edginc/VolCalibInterp.hpp"
#include "edginc/Addin.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/B30360.hpp"
#include "edginc/ActualActual.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/Surface.hpp"
#include "edginc/IRVolBase.hpp"
#include "edginc/IRVol.hpp"
#include "edginc/VolRequestCalib.hpp"
#include "edginc/VolProcessedCalib.hpp"

DRLIB_BEGIN_NAMESPACE

// Make sure the class links
bool VolCalibInterpLoad() { return (VolCalibInterp::TYPE != 0); }

VolCalibInterp::VolCalibInterp()
    : CObject(TYPE)
{
    validatePop2Object();
}

VolCalibInterp::~VolCalibInterp() {}

VolCalibInterp::VolCalibInterp(const CClassConstSP& clazz): 
    CObject(clazz) {}

void VolCalibInterp::validatePop2Object() {}

void VolCalibInterp::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("Simulates a VolCalibInterp!");
    REGISTER(VolCalibInterp, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(vol,              "Vol!");
    FIELD(calibType,        "Calibration type eg. FIX, CMS, FMD or CAP etc.");
    FIELD(calibTenor,       "Calibration Tenor (or explicit maturity date when using FIXED_MATURITY.")
    FIELD(calibVolOverride, "Calibration spot vol override value - vol surfaces are ignored when used.")
    FIELD_MAKE_OPTIONAL(calibVolOverride);

    Addin::registerObjectMethod("VolCalibInterp",
                                Addin::MARKET,
                                "Interpolates for vol calibration.",
                                false,
                                Addin::expandMulti,
                                &VolCalibInterp::goVolCalibInterp);
}

CDoubleMatrixSP VolCalibInterp::goVolCalibInterp()
{
    static const string method = "VolCalibInterp::goVolCalibInterp";
    
    // Get Vol
    MarketObjectSP ptrMO = vol.getMO();
    MarketObject*  pMO = ptrMO.get();
    if (!pMO)
        throw ModelException(method, " - Vol object not found!");
    IRVolBase* pIRVolBase = dynamic_cast<IRVolBase*>(pMO);
    
    // Get Calib Vols
    VolRequestCalib     volReqCalib(calibType, calibTenor, calibVolOverride);
    CVolProcessedSP     volProcessed(pIRVolBase->getProcessedVol(&volReqCalib, 0));
    VolProcessedCalibSP ptrVolProcessedCalib(VolProcessedCalibSP::dynamicCast(volProcessed));
    DoubleArraySP       ptrVolArray = ptrVolProcessedCalib->getVols();
    DateTimeArraySP     ptrSwapStartDates = ptrVolProcessedCalib->getSwapStartDates();
    DateTimeArraySP     ptrSwapMaturityDates = ptrVolProcessedCalib->getSwapMaturityDates();

    CDoubleMatrixSP  ret(new DoubleMatrix(3, ptrVolArray->size()));
    DoubleMatrix&    m = *ret;
    for (int i = 0; i < m.numRows(); ++i)
    {
        m[0][i] = (*ptrVolArray)[i];
        m[1][i] = (double)(*ptrSwapStartDates)[i].toIrDate();
        m[2][i] = (double)(*ptrSwapMaturityDates)[i].toIrDate();
    }
    return (ret);
}

CClassConstSP const VolCalibInterp::TYPE = CClass::registerClassLoadMethod(
    "VolCalibInterp", typeid(VolCalibInterp), load);

/** Returns name of struct */
string VolCalibInterp::getName() const
{
    return "AnwarVolCalibInterp";
}

/** static interpVolsForCalibration */
VolProcessedCalibSP VolCalibInterp::interpVolsForCalibration(
    const ExpiryArray&    xArray,
    const ExpiryArray&    yArray,
    const CDoubleMatrix&  zMatrix,
    DateTime              baseDate,
    IRVol::CalibType      calibType,
    ExpirySP              calibTenor,
    CDoubleSP             calibVolOverride) // defaults to null smart ptr
{
    // Initialise
    static const string method = "VolCalibInterp::interpVolsForCalibration";
    B30360          dcc30360;
    DateTimeArraySP swapStartDates(new DateTimeArray());
    DateTimeArraySP swapMaturityDates(new DateTimeArray());
    DoubleArraySP   vols(new DoubleArray());
    ExpiryArraySP   selectedTenors(new ExpiryArray());
    ExpiryArraySP   selectedExpiries(new ExpiryArray());
    IntArraySP      selectedTenorIndices(new IntArray());
    IntArraySP      selectedExpiryIndices(new IntArray());
    DateTimeArraySP selectedSwapMaturityDates(new DateTimeArray());
    DateTime        swapStart;
    DateTime        swapMaturity;
    int             i, j;
    int             numExps = xArray.size();
    int             numMats = yArray.size();
    
    // Handle vol overrides (deploy if matrix size is scalar or override given explicitly)
    bool   useVolOverride = calibVolOverride.get() ? true : false;
    double volOverride = calibVolOverride.get() ? calibVolOverride->doubleValue() : 0.0;
    if (!useVolOverride && zMatrix.numRows() == 1 && zMatrix.numCols() == 1)
    {
        volOverride = zMatrix[0][0];
        useVolOverride = true;
    }

    // Select interp type
    switch (calibType)
    {
        case IRVol::CAP:
        case IRVol::CMS:
        {
            // Find smallest tenor column nearest to calibTenor
            j = 0;
            swapMaturity = calibTenor->toDate(baseDate);
            while ((j < numMats - 1) && (swapMaturity > yArray[j]->toDate(baseDate)))
                j++;

            if (swapMaturity > yArray[j]->toDate(baseDate))
                throw ModelException(method, "CMS tenor is not in volatility matrix!");

            // Populate results
            for (i = 0; i < numExps; ++i)
            {
                // Swap Start/Maturity Dates
                swapStart = xArray[i]->toDate(baseDate);     // add settleDays here if needed
                swapStartDates->push_back(swapStart);
                swapMaturity = calibTenor->toDate(swapStart);
                swapMaturityDates->push_back(swapMaturity);

                // Selected Nodes & Vols
                selectedTenors->push_back(yArray[j]);
                selectedTenorIndices->push_back(j);
                selectedExpiries->push_back(xArray[i]);
                selectedExpiryIndices->push_back(i);
                selectedSwapMaturityDates->push_back(yArray[j]->toDate(swapStart));
                vols->push_back(useVolOverride ? volOverride : zMatrix[j][i]);
            }
            break;
        }

        case IRVol::FIX:
        {
            // Convert calib tenor to year frac (not using dates)
            MaturityPeriodSP ptrFixedTenor = MaturityPeriodSP::dynamicCast(calibTenor);
            if (!ptrFixedTenor.get())
                throw ModelException(method, "FIX tenor must be expressed as a MaturityPeriod tenor (eg. '1Y'!");
            double fixedMat = ptrFixedTenor->toYears();
            
            // Setup year fraction arrays            
            DoubleArray  swapMats(numMats);
            DoubleArray  optExps(numExps);
            for (j = 0; j < numMats; ++j)
            {
                swapMaturity = yArray[j]->toDate(baseDate);
                swapMats[j] = dcc30360.years(baseDate, swapMaturity);
            }
            for (i = 0; i < numExps; ++i)
            {
                DateTime iOptionExpiry = xArray[i]->toDate(baseDate);
                optExps[i] = (iOptionExpiry.daysDiff(baseDate)) / 365.0;
            }

            // Setup output result structures
            for (i = 0; i < numExps; ++i)
            {
                double curMat = fixedMat - optExps[i];

                // This must be greater than the smallest maturity on our grid
                if (curMat > swapMats[0])
                {
                    // Find maturity tenors
                    j = 0;                                                                       
                    while ((j < numMats - 1) && (curMat >= swapMats[j]))
                        ++j;

                    // Use nearest column (no interpolation)
                    if (2 * curMat < swapMats[j-1] + swapMats[j])
                        --j;

                    // Swap Start Dates
                    swapStart = xArray[i]->toDate(baseDate);
                    swapStartDates->push_back(swapStart);

                    // Swap Maturity Dates - Use Aladdin method to back out swap maturity
                    // Previous ESL method was: swapMaturity = yArray[j]->toDate(swapStart);
                    swapMaturity = MaturityPeriod::toDate((curMat * 12.0) + 0.5, "M", swapStart);
                    swapMaturityDates->push_back(swapMaturity);

                    // Selected Nodes & Vols
                    selectedTenors->push_back(yArray[j]);
                    selectedTenorIndices->push_back(j);
                    selectedExpiries->push_back(xArray[i]);
                    selectedExpiryIndices->push_back(i);
                    selectedSwapMaturityDates->push_back(yArray[j]->toDate(swapStart));
                    vols->push_back(useVolOverride ? volOverride : zMatrix[j][i]);
                }
            }
            break;
        }

        case IRVol::FMD:
        {
            // Find exact column otherwise interpolate linearly across column.
            
            // CalibTenor must be a benchmark date
            BenchmarkDateSP calibDate = BenchmarkDateSP::dynamicCast(calibTenor);
            if (!calibDate)
                throw ModelException(method, "FMD tenor must be an explicit date!");
                
            // Setup year fraction arrays            
            DoubleArray  swapMats(numMats);
            DoubleArray  optExps(numExps);
            for (j = 0; j < numMats; ++j)
            {
                DateTime jSwapMaturity = yArray[j]->toDate(baseDate);
                swapMats[j] = dcc30360.years(baseDate, jSwapMaturity);
            }
            for (i = 0; i < numExps; ++i)
            {
                DateTime iOptionExpiry = xArray[i]->toDate(baseDate);
                optExps[i] = (iOptionExpiry.daysDiff(baseDate)) / 365.0;
            }

            ActualActual dccActAct;
            DateTime     fixedMaturityDate = calibDate->toDate();
            double       fixedMat = dccActAct.years(baseDate, fixedMaturityDate);
            double       vol;
            bool         curMatLessThanFirstTenor = false;
            const double yearFrac1D = 1.0 / 366.0;  // smallest day frac

            for (int i = 0; i < numExps; ++i)
            {
                // Determine current maturity and check if below first maturity point on grid.
                double curMat = fixedMat - optExps[i];
                if (curMat < swapMats[0])
                {
	                curMat = swapMats[0];
                    curMatLessThanFirstTenor = true;
                }

                // Find maturity tenors
                j = 0;                                                                       
                while ((j < numMats - 1) && (swapMats[j] <= curMat))
                    ++j;

                // Swap Start/Maturity Dates
                swapStart = xArray[i]->toDate(baseDate);
                swapStartDates->push_back(swapStart);
                if (curMatLessThanFirstTenor)
                    swapMaturity = xArray[i]->toDate(swapStart);
                else
                    swapMaturity = fixedMaturityDate;
                swapMaturityDates->push_back(swapMaturity);

                if (!useVolOverride)
                {
                    // Interpolate across maturity columns
                    vol = ( zMatrix[j-1][i] * (swapMats[j] - curMat)
                        +   zMatrix[j][i]   * (curMat - swapMats[j-1]) )
                        /   (swapMats[j] - swapMats[j-1]);
                }

                // Setup selected node points that were used in the interpolation.
                if (fabs(curMat - swapMats[j-1]) >= yearFrac1D) 
                {
                    selectedTenors->push_back(yArray[j]);
                    selectedTenorIndices->push_back(j);
                    selectedExpiries->push_back(xArray[i]);
                    selectedExpiryIndices->push_back(i);
                    selectedSwapMaturityDates->push_back(yArray[j]->toDate(swapStart));
                }
                if (fabs(swapMats[j] - curMat) >= yearFrac1D)
                {
                    selectedTenors->push_back(yArray[j-1]);
                    selectedTenorIndices->push_back(j-1);
                    selectedExpiries->push_back(xArray[i]);
                    selectedExpiryIndices->push_back(i);
                    selectedSwapMaturityDates->push_back(yArray[j-1]->toDate(swapStart));
                }

                vols->push_back(useVolOverride ? volOverride : vol);
            }              
            break;
        }
        default:
            throw ModelException(method, "Invalid calibration type!");
    }
    
    smartPtr<VolProcessedCalib> volProcessedCalibSP(new VolProcessedCalib(
        swapStartDates, swapMaturityDates, vols, selectedTenors, selectedTenorIndices, selectedExpiries,
        selectedExpiryIndices, selectedSwapMaturityDates));
    return volProcessedCalibSP;
}
DRLIB_END_NAMESPACE
