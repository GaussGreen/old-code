/****************************************************************
 * Module:        PenGuin
 * Submodule:        
 * File:        voldat.c
 * Function:        
 * Author:        Christian Daher - David Liu
 ****************************************************************/
#include "kvoldat.h"
#include "krate.h"


#include "kutilios.h"
extern        "C" {
#include <math.h>
#include "drlio.h"
#include "drlinter.h"
#include "drlstr.h"
#include "drltime.h"
#include "crxvnfm.h"
#include "stdio.h"

#include "dritkwrp.h"           // TDrWrapper

#include "vnfmdata.h"           // Vnfm 
};

#undef  RIDX
#define RIDX(n,i,j)     ((i) <= (j) ? \
        ((n)*((n)-1)/2 - ((n)-(i))*((n)-(i)-1)/2 + ((j)-(i)-1)) : \
        ((n)*((n)-1)/2 - ((n)-(j))*((n)-(j)-1)/2 + ((i)-(j)-1)))

#define MIN_VOL  0.00001

static
int  FreqToBasis(char freq);

//===============================================================
//

KVolCalibIndex::KVolCalibIndex()
{
        mCalibType = 0;        // nil clibration default
        mFinalAsDate = FALSE;
        mMatInt = 0e0;
        mMatDate = 0L;
}


//===============================================================
//
KVolCalibIndex::~KVolCalibIndex()
{
}


//---------------------------------------------------------------
// Scans a calibration index in a DR wrapper style.
//


KVolCalibIndex::KVolCalibIndex(const char* idxStr)
{
static        char        routine[] = "KVolCalibIndex::KVolCalibIndex(const char*)";
        char        *p;
        char        buf[256];

   try {

        strncpy(buf, idxStr, sizeof(buf));

        if ((p = strstr(buf, "nil"))) {
                // Nil clibration
                mCalibType = 0;

        } else if ((p = strstr(buf, "1m"))) {
                // Base vol
                mCalibType = 1;
                mFinalAsDate = FALSE;
                DppReadFromString("1M", mMatInt);

        } else if ((p = strstr(buf, "3m"))) {
                // Base vol
                mCalibType = 1;
                mFinalAsDate = FALSE;
                DppReadFromString("3M", mMatInt);

        } else if ((p = strstr(buf, "6m"))) {
                // Base vol
                mCalibType = 1;
                mFinalAsDate = FALSE;
                DppReadFromString("6M", mMatInt);

        } else if ((p = strstr(buf, "Cms"))) {
                // Cms calibration
                *p = '\0';
                mCalibType = 1;
                mFinalAsDate = FALSE;
                DppReadFromString(buf, mMatInt);

        } else if ((p = strstr(buf, "Fix"))) {
                // FInal calibration as interval
                *p = '\0';
                mCalibType = 2;
                mFinalAsDate = FALSE;
                DppReadFromString(buf, mMatInt);

        } else {
                // FInal calibration as date
                mCalibType = 2;
                mFinalAsDate = TRUE;
                IF_FAILED_THROW( GtoStringToDate(buf, &mMatDate));
        }


    }
    catch(KFailure) {
        throw KFailure("%s: failed reading `%s'.\n", routine, idxStr);
    }

}


//---------------------------------------------------------------
//

istream&
operator>>(istream& is, KVolCalibIndex& object)
{
static        char        routine[] = "operator>>(istream&,KVolCalibIndex&)";
        char        vs[256];

    try {
        strcpy(vs, getString(is, "calibration index"));

        object = KVolCalibIndex(vs);

        /*DppScanCalibIndex(
                vs,
                &object.mCalibType,
                &object.mFinalDate,
                &object.mMatInt,
                &object.mMatDate);*/

        return(is);
    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}

//---------------------------------------------------------------
//

ostream&
operator<<(ostream& os, const KVolCalibIndex& object)
{
        switch (object.mCalibType) {
        case 0:
                os << "nil";
                break;
        case 1:
                // CMS calibration
                os << object.mMatInt << "Cms";
                break;
        case 2:
                // Final 
                if (object.mFinalAsDate) {
                        os << DrlTDatePrint(NULL, object.mMatDate);
                } else {
                        os << object.mMatInt << "Fix";
                }
                break;
        default:
                os << format("(invalid type %d)", object.mCalibType);
        }

        return (os);
}


//---------------------------------------------------------------
//

bool
KVolCalibIndex::operator==(const KVolCalibIndex& object) const
{
        if (mCalibType != object.mCalibType)
                return (false);
        if (mCalibType == 0)
                return (true);
        if (mFinalAsDate != object.mFinalAsDate)
                return (false);
        if (mFinalAsDate) {
                if (mMatDate != object.mMatDate)
                        return (false);
        } else {
                if (mMatInt != object.mMatInt)
                        return (false);
        }
        return (true);
}








//===============================================================
// 
//
//===============================================================


//---------------------------------------------------------------
//

KVolDiag::KVolDiag()
{
        mVolType = LOGVOL;
}

//---------------------------------------------------------------
//

KVolDiag::KVolDiag(
        const KVector(TDate) &volExpDates,
        const KVector(TDate) &volMatDates,
        const KVector(int) &volFreqs,
        const KVector(double) &volRates)
{
static        char        routine[] = "KVolDiag::KVolDiag";

    try {
        mVolType = LOGVOL;

        ASSERT_OR_THROW(volExpDates.size() == volMatDates.size());
        ASSERT_OR_THROW(volExpDates.size() == volFreqs.size());
        ASSERT_OR_THROW(volExpDates.size() == volRates.size());

        mVolMats.resize(volExpDates.size());

        mVolDates = volExpDates;
        //
        // Convert vol maturity dates to 
        //
        for (int i=0; i<volMatDates.size(); i++) {
                IF_FAILED_THROW( GtoDayCountFraction(
                        volExpDates[i],
                        volMatDates[i],
                        GTO_B30_360,
                        &mVolMats[i]));
                mVolRates.push_back(volRates[i]);
        }

        mVolFreqs = volFreqs;

    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}




//---------------------------------------------------------------
//

KVolDiag::KVolDiag(
        KVolType              volType,
        const KVector(TDate)  &volExpDates,
        const KVector(TDate)  &volMatDates,
        const KVector(int)    &volFreqs,
        const KVector(double) &volRates)
{
static        char        routine[] = "KVolDiag::KVolDiag";

    try {

        ASSERT_OR_THROW(volExpDates.size() == volMatDates.size());
        ASSERT_OR_THROW(volExpDates.size() == volFreqs.size());
        ASSERT_OR_THROW(volExpDates.size() == volRates.size());

        mVolMats.resize(volExpDates.size());

        mVolType = volType;
        mVolDates = volExpDates;
        //
        // COnvert vol maturity dates to 
        //
        for (int i=0; i<volMatDates.size(); i++) {
                IF_FAILED_THROW( GtoDayCountFraction(
                        volExpDates[i],
                        volMatDates[i],
                        GTO_B30_360,
                        &mVolMats[i]));

                mVolRates.push_back(volRates[i]);
        }

        mVolFreqs = volFreqs;

    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}



//---------------------------------------------------------------
//

KVolDiag::KVolDiag(
        const KVolCalibIndex& calibIdx,                // (I) calibration index
        TDrWrapperData *drWrapData)                 // (I) DRW data 
{
static        char        routine[] = "KVolDiag::KVolDiag";
        int     numVolDates, i;
        TDate   *volDates = NULL;
        double  *volMat = NULL;
        int     *volFreq = NULL;
        double  *volRates = NULL;
        TDate        valueDate;

    try {

        // Wrapper volatility is always %
        //
        mVolType = LOGVOL;

        if (calibIdx.IsNil()) {
                // No interp with Nil calibration
                throw KFailure("%s: can't interpolate volatilities with"
                        " nil calibration index.\n", routine);
        } else {
                IF_FAILED_THROW( DriTDrWrapperDataGetInterpVol(
                        drWrapData,
                        //calibIdx.mCalibFinal,
                        calibIdx.IsFinal(),
                        calibIdx.FinalDate(),
                        calibIdx.Tenor(),
                        &numVolDates,
                        &volDates,
                        &volMat,
                        &volFreq,
                        &volRates));

 
                // Get value date 
                valueDate = drWrapData->fDiscZcCurve->fBaseDate;

                for (i=0; i<=numVolDates-1; i++)
                {
                    // Only add dates past value date (plus buffer)
                    //
                    if (volDates[i] > valueDate+2) {
                        mVolDates.push_back(volDates[i]);
                        mVolMats.push_back(volMat[i]);
                        mVolFreqs.push_back(volFreq[i]);
                        mVolRates.push_back(volRates[i]);
                    }
                }
        }


        FREE(volDates);
        FREE(volMat);
        FREE(volFreq);
        FREE(volRates);
    }
    catch (KFailure) {
        FREE(volDates);
        FREE(volMat);
        FREE(volFreq);
        FREE(volRates);
        throw KFailure("%s: failed.\n", routine);
    }
}



//---------------------------------------------------------------
//
void
KVolDiag::IRVolDiag(
        CRX_INPUT           *crxInput)                 // (I) MAW data 
{
static        char        routine[] = "KVolDiag::IRVolDiag";
        
        int      i;
        TDate    baseDate;
        double   swapMat;    // underlying swap maturity in years
        int      freq;       // underlying swap coupon freqency per year

        IR_INPUT *irInput = NULL;

    try {

        // Wrapper volatility is always %
        //
        mVolType = LOGVOL;

        // Diffuse IR curve
        //
        irInput = &(crxInput->irInput[0]);

        // Get base date 
        baseDate = irInput->AlibBaseDate;

        // Swap Freq
        if (toupper(irInput->Freq) == 'I')
        {
            IF_FAILED_THROW( GtoDayCountFraction(irInput->AlibSwapSt[0],
                                                 irInput->AlibSwapMat[0],
                                                 GTO_B30_360,
                                                 &swapMat));
            freq = (int) (1./swapMat);
        }
        else
            freq = FreqToBasis(irInput->Freq);


        for (i=0; i<=irInput->NbVol-1; i++)
        {
            // Only add dates past value date (plus buffer)
            //
            if (irInput->AlibVolDate[i] > baseDate+2) {
                IF_FAILED_THROW( GtoDayCountFraction(irInput->AlibSwapSt[i], 
                                                     irInput->AlibSwapMat[i], 
                                                     GTO_B30_360,
                                                     &swapMat));

                mVolDates.push_back(irInput->AlibVolDate[i]);
                mVolMats.push_back(swapMat);
                mVolFreqs.push_back(freq);
                mVolRates.push_back(irInput->Vol[i]);
            }
        }


    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}

//---------------------------------------------------------------
//
void
KVolDiag::IRVolDiag(
        BS_INPUT           *bsInput)                 // (I) MAW data 
{
static        char        routine[] = "KVolDiag::IRVolDiag";
        
        int      i;
        TDate    baseDate;
        double   swapMat;    // underlying swap maturity in years
        int      freq;       // underlying swap coupon freqency per year

        IR_INPUT *irInput = NULL;

    try {

        switch(bsInput->irInput->VolUnit)
        {
        case 0:
            mVolType = LOGVOL;
            break;
        case 1:
            mVolType = NORMVOL;
            break;
        default:
            throw KFailure("%s failed: Invalid vol type %d\n",
                            routine,
                            bsInput->irInput->VolUnit);
        }

        // Diffuse IR curve
        //
        irInput = &(bsInput->irInput[0]);

        // Get base date 
        baseDate = irInput->AlibBaseDate;

        // Swap Freq
        if (toupper(irInput->Freq) == 'I')
        {
            IF_FAILED_THROW( GtoDayCountFraction(irInput->AlibSwapSt[0],
                                                 irInput->AlibSwapMat[0],
                                                 GTO_B30_360,
                                                 &swapMat));
            freq = (int) (1./swapMat);
        }
        else
            freq = FreqToBasis(irInput->Freq);


        for (i=0; i<=irInput->NbVol-1; i++)
        {
            // Only add dates past value date (plus buffer)
            //
            if (irInput->AlibVolDate[i] > baseDate+2) {
                IF_FAILED_THROW( GtoDayCountFraction(irInput->AlibSwapSt[i], 
                                                     irInput->AlibSwapMat[i], 
                                                     GTO_B30_360,
                                                     &swapMat));

                mVolDates.push_back(irInput->AlibVolDate[i]);
                mVolMats.push_back(swapMat);
                mVolFreqs.push_back(freq);
                mVolRates.push_back(irInput->Vol[i]);
            }
        }


    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}

//---------------------------------------------------------------
//
void
KVolDiag::CRVolDiag(
        CRX_INPUT           *crxInput)                 // (I) MAW data 
{
static        char        routine[] = "KVolDiag::CRVolDiag";
        
        int      i;
        TDate    baseDate;
        double   swapMat;    // underlying swap maturity in years
        int      freq;       // underlying swap coupon freqency per year

        CR_INPUT *crInput = NULL;

    try {

        // Wrapper volatility is always %
        //
        mVolType = LOGVOL;

        // CR curve
        //
        crInput = &(crxInput->crInput[0]);

        // Get base date 
        baseDate = crInput->AlibBaseDate;

        // Swap Freq
        if (toupper(crInput->Freq) == 'I')
        {
            IF_FAILED_THROW( GtoDayCountFraction(crInput->AlibSwapSt[0],
                                                 crInput->AlibSwapMat[0],
                                                 GTO_B30_360,
                                                 &swapMat));
            freq = (int) (1./swapMat);
        }
        else
            freq = FreqToBasis(crInput->Freq);


        for (i=0; i<=crInput->NbVol-1; i++)
        {
            // Only add dates past value date (plus buffer)
            //
            if (crInput->AlibVolDate[i] > baseDate+2) {
                IF_FAILED_THROW( GtoDayCountFraction(crInput->AlibSwapSt[i], 
                                                     crInput->AlibSwapMat[i], 
                                                     GTO_B30_360,
                                                     &swapMat));

                mVolDates.push_back(crInput->AlibVolDate[i]);
                mVolMats.push_back(swapMat);
                mVolFreqs.push_back(freq);

                // Minimum vol set to 0.01%
                mVolRates.push_back(crInput->Vol[i]);
            }
        }


    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}

//---------------------------------------------------------------
//
void
KVolDiag::SPVolDiag(
        BS_INPUT           *bsInput)                 // (I) MAW data 
{
static        char        routine[] = "KVolDiag::SPVolDiag";
        
        int      i;
        TDate    baseDate;
        double   swapMat;    // underlying swap maturity in years
        int      freq;       // underlying swap coupon freqency per year
        double   bsVolInterp;

        SP_INPUT *spInput = NULL;
        IR_INPUT *irInput = NULL;

    try {

            switch(bsInput->spInput->VolUnit)
            {
            case 0:
                mVolType = LOGVOL;
                break;
            case 1:
                mVolType = NORMVOL;
                break;
            default:
                throw KFailure("%s failed: Invalid vol type %d\n",
                            routine,
                            bsInput->spInput->VolUnit);
            }


        // BS curve
        //
        spInput = &(bsInput->spInput[0]);
        irInput = &(bsInput->irInput[0]);

        // Get base date 
        baseDate = spInput->AlibBaseDate;

        // vol type
        if (spInput->VolUnit == 1)
            mVolType = NORMVOL;

        swapMat = 1e0 / (double)(4);
        freq    = 0;
        for (i=0; i<=irInput->NbVol-1; i++)
        {
            // Only add dates past value date (plus buffer)
            //
            if (irInput->AlibVolDate[i] > baseDate+2) {
                IF_FAILED_THROW( DrlTDateLinearInterp1d(
                    spInput->AlibVolDate,
                    spInput->Vol,
                    spInput->NbVol,
                    irInput->AlibVolDate[i],
                    &bsVolInterp))

                mVolDates.push_back(irInput->AlibVolDate[i]);
                mVolMats.push_back(swapMat);
                mVolFreqs.push_back(freq);
                mVolRates.push_back(bsVolInterp);
            }
        }


    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}

static
int  FreqToBasis(char freq)
{
    switch (freq) {
    case 'A':
        return  1;
        break;
    case 'S':
        return  2;
        break;
    case 'Q':
        return  4;
        break;
    case 'M':
        return  12;
        break;
    }

    return -1;
}



//---------------------------------------------------------------
//
void
KVolDiag::BasisVolDiag(
        const char*           bsVolType,
        const KVector(TDate)  &bsVolExpDates,
        const KVector(double) &bsVolRates,
        int                   bsVolFreq,
        KVolDiag              &irVolDiag)
{
static        char        routine[] = "KVolDiag::BasisVolDiag";
        int     i;

        double        bsVolInterp;

        KVector(TDate)  volExpDates;
        KVector(double) volRates;

 

    try {

        if (bsVolExpDates.empty())
                throw KFailure("%s: need basis vol input.\n", routine);

 
        // Check size consistency
        ASSERT_OR_THROW(bsVolExpDates.size() == bsVolRates.size());

        // Copy to non-const vector
        volExpDates = bsVolExpDates;
        volRates    = bsVolRates;
 

        for (i=0; i<irVolDiag.Size(); i++)
        {
                // interpolate basis vol on IR vol benchmark dates
                IF_FAILED_THROW(DrlTDateLinearInterp1d(
                             &volExpDates[0],
                             &volRates[0],
                             bsVolExpDates.size(),
                                 irVolDiag.mVolDates[i],
                                 &bsVolInterp));
                                
                mVolRates[i] = bsVolInterp;
                mVolFreqs[i] = bsVolFreq;
 
                // Constant maturity specified by vol curve frequency
                //
                mVolMats[i]  = 1e0/(double)(bsVolFreq == 0 ? 1 : bsVolFreq);
        }

 
        if (toupper(bsVolType[0]) == 'N')
                mVolType = NORMVOL;
        else if (toupper(bsVolType[0]) == 'L')
                mVolType = LOGVOL;
        else
                throw KFailure("%s: invalid basis vol type (%s).\n",
                                routine, bsVolType);


    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }

}




//---------------------------------------------------------------
//
KVolDiag&
KVolDiag::operator=(const KVolDiag& volDiag)
{
 
    if (this != &volDiag)
    {
        mVolType  = volDiag.mVolType;
        mVolDates = volDiag.mVolDates;
        mVolMats  = volDiag.mVolMats;
        mVolFreqs = volDiag.mVolFreqs;
        mVolRates = volDiag.mVolRates;
    }
 
    return(*this);
 
}



//---------------------------------------------------------------
//

ostream&
operator<<(ostream& os, const KVolDiag& volDiag)
{
        int        n, idx;

        os << "#KVolDiag:" << endl;

        n = volDiag.mVolDates.size();
        if ((volDiag.mVolMats.size() != n) ||
            (volDiag.mVolFreqs.size() != n) ||
            (volDiag.mVolRates.size() != n)) {
                throw KFailure("operator<<(ostream&,KVolDiag&): "
                        "inconsistent size\n");
        }
        switch (volDiag.mVolType) {
        case NORMVOL:
                os << "#Volatility type\n NORMVOL" << endl;
                break;
        case LOGVOL:
                os << "#Volatility type\n LOGVOL" << endl;
                break;
        }

        if (n == 0) {
                os << "(empty)" << endl;
        } else {
                os << format("# DATE     MAT      FREQ  VOLATILITY") << endl;
                for (idx=0; idx<n; idx++) {
                    os << format("%s %8.4f  %d  %12.8f",
                        DrlTDatePrint(NULL, volDiag.mVolDates[idx]),
                        volDiag.mVolMats[idx],
                        volDiag.mVolFreqs[idx],
                        volDiag.mVolRates[idx]) << endl;
                }
        }

        return (os);
}


//---------------------------------------------------------------
ostream&
KVolDiag::YacctionWrite(ostream& os, int indent)
{
        char buf[64];

        TDate   matDate;

        KDateInterval dtInterval = KDateInterval(mVolMats.front());

        if (IS_ALMOST_ZERO(mVolRates[0]+1e0))        // nil
        {
                os << "# FirstCalibIndex" << endl;
                os << "nil" << endl;
                os << "# SecondCalibIndex" << endl;
                os << "nil" << endl;
        }
        else if (fabs(mVolMats.front()
                -mVolMats.back()) < 0.0027)   // CMS
        {
                // Always convert to "M"onths.
                
                strcpy(buf,
                       GtoFormatDateInterval(&((TDateInterval)dtInterval)));
 
                os << "# FirstCalibIndex" << endl;
                os << buf << "Cms" << endl;
                os << "# SecondCalibIndex" << endl;
                os << buf << "Cms" << endl;
        }
        else {  // Fixed date
                // The GtoTDateAdvanceYears does more properly
                // because it adds the # years first, then the
                // fraction of year, while the TDateInterval
                // would assume 365/year.

                //matDate = mVolDates.front() + dtInterval;
                IF_FAILED_THROW(GtoTDateAdvanceYears(
                                        mVolDates.front(),
                                        mVolMats.front(),
                                        &matDate));        
 
                os << "# FirstCalibIndex" << endl;
                os << GtoFormatDate(matDate) <<  endl;
                os << "# SecondCalibIndex" << endl;
                os << GtoFormatDate(matDate) <<  endl;
        }

        return (os);

}


//---------------------------------------------------------------
ostream&
KVolDiag::BasisYacctionWrite(ostream& os, int indent)
{
        os << "# Basis vol type ('L'ognormal % vol, 'N'ormal bp vol)" << endl;
        if (mVolType == LOGVOL)
                os << "L" << endl;
        else
                os << "N" << endl;

        os << "# Number of entries" << endl;
        os << mVolDates.size() << endl;

        os << "# Basis spread volatility dates and volatilities in pct (-100%=nil)" << endl;
        for (int i=0; i < mVolDates.size(); i++)
        {
                os << GtoFormatDate(mVolDates[i]) << "\t"
                   << format("%10.5f\n", mVolRates[i]*100);
        }

        return (os);
}




//---------------------------------------------------------------

bool
KVolDiag::IsValid() const
{
    try {

        ASSERT_OR_THROW(mVolDates.size() == mVolDates.size());
        ASSERT_OR_THROW(mVolMats.size() == mVolDates.size());
        ASSERT_OR_THROW(mVolFreqs.size() == mVolDates.size());
        ASSERT_OR_THROW(mVolRates.size() == mVolDates.size());

        for (int idx=1; idx<mVolDates.size(); idx++) {
                if (mVolDates[idx-1] >= mVolDates[idx])
                        throw KFailure("date %d (%s) >= date %d (%s).\n",
                                idx-1,
                                DrlTDatePrint(NULL, mVolDates[idx-1]),
                                idx,
                                DrlTDatePrint(NULL, mVolDates[idx]));
        }


        return (true);
    }
    catch (KFailure) {
        return (false);
    }
}


//---------------------------------------------------------------
//

bool
IsConsistent(const KVolDiag& volDiag1,const KVolDiag& volDiag2)
{

        if (volDiag1.mVolDates.size() != volDiag2.mVolDates.size())
                return (false);

        for (int idx=0; idx<volDiag1.mVolDates.size(); idx++)
                if (volDiag2.mVolDates[idx] != volDiag2.mVolDates[idx])
                        return (false);

        return (true);
}

//---------------------------------------------------------------

KVolDiag&
KVolDiag::operator=(double value)
{
        for (int i=0; i<mVolRates.size(); i++) {
                mVolRates[i] = value;
        }
        return (*this);
}


//---------------------------------------------------------------
//

double
KVolDiag::VolInterp(TDate date)
{
static        char        routine[] = "KVolDiag::VolInterp";
        double        value;

        ASSERT_OR_THROW(mVolRates.size() == mVolDates.size());
        IF_FAILED_THROW( DrlTDateLinearInterp1d(
                &mVolDates[0],
                &mVolRates[0],
                mVolRates.size(),
                date,
                &value));

        return(value);
}


//===============================================================
// 
//
//===============================================================


//--------------------------------------------------------------------
//

void        
DppCalibIRVolDiag(
        const KMrParam &irMrParam,        // (I) full dimension mr info
        const KSmileParam &irSmileParam,// (I) ir smile info
        int nIRDim,                        // (I) ir dimension
        const KVolDiag &irVolDiag,        // (I) ir vol info
        const KZCurve &diffCurve,        // (I) diffuse zero curve
        KVector(TDate) &volDates,        // (O) ir spot vol dates
        KVector(KVector(double)) &factVolCurves)// (O) ir factor spot vols
{
static  char    routine[] = "DppCalibIRVolDiag";
 
        int        i, j, k, idx,
                rIdx;

        int     irDim;

        int        idxEnd;

        int        numVolDates;
        TCurve        *diffTCurve = (TCurve*)diffCurve;

        TDate        valueDate = diffTCurve->fBaseDate;
        TDate   voldateSt;

        VnfmData *vnfmData = NULL;                // Vnfm data structure

        KVector(double) factVol;

        KVolDiag        volDiag = irVolDiag;        // Make copy ir vol info

    try{

        // Check irDim valid
        if (nIRDim <= 0) {
                throw KFailure("%s: bad IR dim %d.\n", routine, nIRDim);
        }

        irDim = nIRDim;

        //-----------------------------------------------
        // Special case of the nil calibration
        // The spot curves are set to 1.0
        //-----------------------------------------------




        if (volDiag.Size() == 0) {
                throw KFailure("%s: can't handle empty vol curves.\n", routine);
        }
        ASSERT_OR_THROW(diffTCurve->fNumItems >= 1);

        //-----------------------------------------------
        // Check base date of IR diffuse zc curve and vol dates.
        // We may need to offset by 1 the input vector
        // if value date is NOT present in the vol dates
        // And insert a dummy value corresponding to valueDate
        // in mIRVolMats, mIRVolFreqs, and mIRVolRates.
        //-----------------------------------------------
 
        voldateSt = volDiag.mVolDates[0];
        if(valueDate < voldateSt){
                volDiag.mVolDates.insert(volDiag.mVolDates.begin(),
                                         valueDate);
                volDiag.mVolMats.insert(volDiag.mVolMats.begin(),
                                         volDiag.mVolMats[0]);
                volDiag.mVolFreqs.insert(volDiag.mVolFreqs.begin(),
                                         volDiag.mVolFreqs[0]);
                volDiag.mVolRates.insert(volDiag.mVolRates.begin(),
                                         volDiag.mVolRates[0]);
        }
        numVolDates = volDiag.mVolDates.size();

        if (irSmileParam.mNumIter<0) {
            /* INTERPRET VOLS AS SPOT VOLS */

            /* This is only meaningful if there is just one IR dimension */
            if (irDim!=1) {
                throw KFailure("%s: You may only calibrate "
                    "IR spot vols directly if there is only 1 IR dimension: "
                    "you should use \"nil\" vol calibration, or "
                    "reduce the IR dimensionality.\n", routine);
            }
            /*=================================================================
             * IF NEGATIVE numIter, THEN JUST TAKE GIVEN VOLS AS SPOT VOLS
             * but scale them by the weight factor first.
             *===============================================================*/
            //factVolCurves.clear();
            double alpha2 = 1e0 / irMrParam.AlphaNorm();
            alpha2 = alpha2 * alpha2;
            //volDiag.mVolDates;
            for (idx=0; idx<=volDiag.mVolDates.size()-1; idx++)
            {
                factVol.push_back(volDiag.VolInterp(volDiag.mVolDates[idx]) * alpha2);
            
            }
            //factVolCurves.insert(factVolCurves.end(), factVol);


        } else if (!IS_ALMOST_ZERO(volDiag.mVolRates[0]+1e0)) {

	        //-----------------------------------------------
	        // Regular case: set up a Vnfm, etc.
	        //-----------------------------------------------
	
	        // Create Vnfm structure and fill it with parameters
	        //
	        vnfmData = VnfmNew(numVolDates, irDim, diffTCurve->fNumItems+1);
	        ASSERT_OR_THROW(vnfmData != NULL);
	
			//
	        // Copy parameters and volatilities 
	        //
	        for(j=0; j<=vnfmData->fNf-1; j++) {
	                vnfmData->fBeta[j]  = irMrParam.mBeta[j];
	                vnfmData->fAlpha[j] = irMrParam.mAlpha[j];
	        }
	        for (i=0; i<= numVolDates-1; i++) {
	                vnfmData->fDate[i] = volDiag.mVolDates[i];
	
                // Floor vol at 0.00001
                volDiag.mVolRates[i] = MAX(volDiag.mVolRates[i], MIN_VOL);
	 
                // Compute the time line
                GtoDayCountFraction(valueDate, vnfmData->fDate[i],
                        GTO_ACT_365F, &vnfmData->fTime[i]);
 
                for(j=0; j<=vnfmData->fNf-1; j++) {
                        vnfmData->fSigma[j][i] = 0.0;
                }


 
                for(j=0;   j<=vnfmData->fNf-1; j++)
                for(k=j+1; k<=vnfmData->fNf-1; k++) {
                        rIdx = RIDX(vnfmData->fNf,j,k);
                        vnfmData->fRho[rIdx][i] = irMrParam.mRho[rIdx];
                }
        	}

	        //
	        // Setup smile
	        //
	        IF_FAILED_THROW( VnfmSetSmile(
	                vnfmData,
	                irMrParam.mBackboneQ,        // backbone q (0=normal, 1=log)
	                irSmileParam.mQ1,        // smile params 
	                irSmileParam.mQ2,
	                irSmileParam.mQF));


	        /* set the zero curve */
	        vnfmData->fZcCurve = diffTCurve;

	#ifndef VNFM_V5X
	        /* set fZTime */
	        vnfmData->fZTime[0] = 0.0;
	                        /* equals REFDATE or (vnfmData)->fDate[0] */
	        for (i=1; i<=diffTCurve->fNumItems; i++) {
	            GtoDayCountFraction(vnfmData->fDate[0], 
	                                diffTCurve->fArray[i-1].fDate,
	                                GTO_ACT_365F, 
	                                &(vnfmData->fZTime[i]));
	        }
	#endif  /*VNFM_V5X*/


        	IF_FAILED_THROW(VnfmCheckValid(vnfmData));


	        // compute coeffs
	        IF_FAILED_THROW(VnfmComputeCoeff(vnfmData));
	








        

	        // find last cali time = 1st negative maturity 
	        for (idxEnd=1; idxEnd<=vnfmData->fNDates-1; idxEnd++) {
	            if (volDiag.mVolMats[idxEnd] <= 0e0) {
	                break;
	            }
			}
        	idxEnd--;

	        // Perform the IR vol calibration
	        IF_FAILED_THROW( VnfmVolCalib1VArbitrary(
	                vnfmData,
	                0,     
	                idxEnd,    
	                &volDiag.mVolMats[0],
	                &volDiag.mVolFreqs[0],
	                &volDiag.mVolRates[0],
	                volDiag.mVolType,
	                TRUE,
	                NULL));
	


	        for (idx=0; idx<=volDiag.mVolDates.size()-1; idx++)
	            factVol.push_back(vnfmData->fSigma[0][idx]);
	
        } else {
                //-----------------------------------------------
                // NIL calibration: set spot vols to 1/Alpha,
        		// to take into account the Bbq factor (=Alpha)
                //-----------------------------------------------
                double        spotVol = 1e0 / irMrParam.AlphaNorm();

                for (idx=0; idx<=volDiag.mVolDates.size()-1; idx++)
                {
                            factVol.push_back(spotVol);
                }
        }



        // Clear factVolCurves
        //factVolCurves.clear();

        volDates = volDiag.mVolDates;
        factVolCurves.insert(factVolCurves.end(), irDim, factVol);





        // Print the vol factors
#ifdef        _SKIP
        if(debugLevel > 10) {
            dppLog << endl;
            dppLog << "============ FACTOR VOLS ============" << endl;
            dppLog << "Date        ";
            for (idx=1; idx<=irDim; ++idx)
                   dppLog << format("Factor %-3d", idx);

            dppLog << endl;

            for (idx=0; idx<=volDiag.mVolDates.size()-1; idx++)
            {
                    dppLog << format("%-12s",
                                GtoFormatDate(volDiag.mVolDates[idx]));

                    for (KVector(KVector(double))::iterator 
                        itVol=factVolCurves.begin();
                        itVol!=factVolCurves.end(); ++itVol)
                        dppLog << format("%-10.6f", (*itVol)[idx]);

                    dppLog << endl;
            }
        }
#endif


        VnfmFree(vnfmData);

    }
    catch (KFailure) {
        VnfmFree(vnfmData);
        throw KFailure("%s: failed.\n", routine);
    }
}




//--------------------------------------------------------------------
//
/*---------------------------------------------------------------
 * Calibration of spread spot volatilities given
 * the base volatilities of spread.
 * Analagous to DppCalibIRVolDiag, but with all the Q and B
 * coeff set to 1.  The zero curve in "that" is obsolete in Q and B
 * calculation, but used instead to store the forward spread curve
 * for the calculation of 2q vol correction.
 */
void        
DppCalibSpreadVolDiag(
        const KMrParam &bsMrParam,                // (I) full dimension mr info
        const KSmileParam &bsSmilePar,        // (I) basis smile info
        int nBSDim,                        // (I) basis dimension
        TDate baseDate,                        // (I) base date of ZC
        KSpd  bsType,                        // (I) SUB_SPREAD, ADD_SPREAD or PER_SPREAD
        const KZCurve        &refCurve,        // (I) reference zero curve
        const KZCurve        &bsCurve,        // (I) basis zero curve
        const KDayCc        &refDCC,        // (I) reference rate DCC
        const KDayCc        &bsDCC,                // (I) basis rate DCC
        const KVolDiag        &bsVolDiag,        // (I) basis vol info
        KVector(TDate)        &volDates,        // (O) spd spot vol dates
        KVector(KVector(double)) &factVolCurves)// (I/O) basis spot vols
{
static  char    routine[] = "DppCalibSpreadVolDiag";
 
        int        i, j, k, idx,
                rIdx;

        int     bsDim;

        int        idxEnd;

        int        numVolDates;

        TDate        voldateSt;
        KVolDiag        volDiag = bsVolDiag;        // Make copy basis vol info

        // Dummy zero curve
        TCurve        *zcCurve = NULL;
        TDate        *zcDates = NULL;
        double        *zcRates = NULL;

        VnfmData *vnfmData = NULL;                // Vnfm data structure

        KVector(double) factVol;

        double        bsFwdRate,
                refFwdRate,
                spread;

        TDate        resetDate;

        KRate        *bsRate = NULL, *refRate = NULL;

    try{

        // Check bsDim valid
        if (nBSDim <= 0) {
                throw KFailure("%s: bad basis dim %d.\n", routine, nBSDim);
        }


        bsDim = nBSDim;


        // Check base date of IR diffuse zc curve and vol dates.
        // We may need to offset by 1 the input vector
        // if value date is NOT present in the vol dates
        // And insert a dummy value corresponding to valueDate
        // in mIRVolMats, mIRVolFreqs, and mIRVolRates.
 
        voldateSt = volDiag.mVolDates[0];
        if(baseDate < voldateSt){
                volDiag.mVolDates.insert(volDiag.mVolDates.begin(),
                                         baseDate);
                volDiag.mVolRates.insert(volDiag.mVolRates.begin(),
                                         volDiag.mVolRates[0]);
                volDiag.mVolFreqs.insert(volDiag.mVolFreqs.begin(),
                                         volDiag.mVolFreqs[0]);
                volDiag.mVolMats.insert(volDiag.mVolMats.begin(),
                                        volDiag.mVolMats[0]);
        }

        numVolDates = volDiag.mVolDates.size();
        
        //-----------------------------------------------
        // Check if not "NIL" calibration (represented by -1 vols)
        //-----------------------------------------------
        if (!IS_ALMOST_ZERO(volDiag.mVolRates[0]+1e0)) {

        //-----------------------------------------------
        // Regular case: set up a Vnfm, etc.
        //-----------------------------------------------

        // Create Vnfm structure and fill it with parameters
        //
        vnfmData = VnfmNew(numVolDates, bsDim, numVolDates+1);
        ASSERT_OR_THROW(vnfmData != NULL);

        // Allocate memeory for zcCurve
        zcDates = new TDate[numVolDates];
        zcRates = new double[numVolDates];

        if (zcDates == NULL ||
            zcRates == NULL)
                throw KFailure("%s: failed to allocate memory for zcDates "
                              "or zcRates.\n",
                                routine);


        for(j=0; j<=vnfmData->fNf-1; j++) {
                vnfmData->fBeta[j]  = bsMrParam.mBeta[j];
                vnfmData->fAlpha[j] = bsMrParam.mAlpha[j];
        }

        for (i=0; i<= numVolDates-1; i++) {
                resetDate = volDiag.mVolDates[i];

                // Floor vol at 0.00001
                volDiag.mVolRates[i] = MAX(volDiag.mVolRates[i], MIN_VOL);


                vnfmData->fDate[i] = resetDate;
 
                // Compute the time line
                GtoDayCountFraction(baseDate, resetDate,
                        GTO_ACT_365F, &vnfmData->fTime[i]);
 
                for(j=0; j<=vnfmData->fNf-1; j++) {
                        vnfmData->fSigma[j][i] = 0.0;
                }

 
                for(j=0;   j<=vnfmData->fNf-1; j++)
                for(k=j+1; k<=vnfmData->fNf-1; k++) {
                        rIdx = RIDX(vnfmData->fNf,j,k);
                        vnfmData->fRho[rIdx][i] = bsMrParam.mRho[rIdx];
                }

                // Compute the forward spread
                //
                bsRate = new KRate(
                                KDateInterval(volDiag.mVolMats[i]),
                                KDateInterval(volDiag.mVolMats[i]),
                                bsDCC,                
                                KDateInterval(0e0),        // spotOffset
                                0e0,        // spread
                                1e0);        // non-zero weight for floating rate

                refRate = new KRate(
                                KDateInterval(volDiag.mVolMats[i]),
                                KDateInterval(volDiag.mVolMats[i]),
                                refDCC,                
                                KDateInterval(0e0),        // spotOffset
                                0e0,        // spread
                                1e0);        // non-zero weight for floating rate

        
                bsFwdRate  = bsRate->Forward(bsCurve, resetDate);
                refFwdRate = refRate->Forward(refCurve, resetDate);

                if (bsType == SUB_SPREAD)
                        spread = refFwdRate - bsFwdRate;    /* Basis = Libor - S */
                else if (bsType == ADD_SPREAD)          /* Basis = Libor + S */
                        spread = bsFwdRate - refFwdRate;
        else                                    /* Basis = Libor * S */
                        spread = bsFwdRate / refFwdRate;


                zcDates[i] = resetDate;
                zcRates[i] = spread;        

                delete refRate;
                refRate = NULL;

                delete bsRate;
                bsRate = NULL;
        }


        // Create the spread curve for 2q vol correction.
        //
        ASSERT_OR_THROW((zcCurve = GtoMakeTCurve(
                                        baseDate,
                                        zcDates,
                                        zcRates,
                                        numVolDates,
                                        1e0,
                                        GTO_ACT_365F)) != NULL);

        // set the zero curve
        vnfmData->fZcCurve = zcCurve;

        //
        // Setup smile
        //
        IF_FAILED_THROW( VnfmSetSmile(
                vnfmData,
                bsMrParam.mBackboneQ,        // backbone q (0=normal, 1=log)
                bsSmilePar.mQ1,        // smile params 
                bsSmilePar.mQ2,
                bsSmilePar.mQF));



#ifndef VNFM_V5X
        /* set fZTime */
        vnfmData->fZTime[0] = 0.0;
                        /* equals REFDATE or (vnfmData)->fDate[0] */
        for (i=1; i<=zcCurve->fNumItems; i++) {
            GtoDayCountFraction(vnfmData->fDate[0], 
                                zcCurve->fArray[i-1].fDate,
                                GTO_ACT_365F, 
                                &(vnfmData->fZTime[i]));
        }
#endif  /*VNFM_V5X*/


        IF_FAILED_THROW(VnfmCheckValid(vnfmData));

        // compute coeffs
        IF_FAILED_THROW(VnfmComputeCoeff(vnfmData));

        idxEnd = vnfmData->fNDates-1;

        // Perform the spread vol calibration
        // where the B or Q coefficients are 1
        //
        IF_FAILED_THROW( VnfmSpreadVolCalib1VArbitrary(
                vnfmData,
                0,     
                idxEnd,    
                &volDiag.mVolRates[0],
                volDiag.mVolType,
                TRUE,
                NULL));

        for (idx=0; idx<=volDiag.mVolDates.size()-1; idx++)
            factVol.push_back(vnfmData->fSigma[0][idx]);

        } else {
                //-----------------------------------------------
                // NIL calibration: set spot vols to 1/Alpha,
        // to take into account the Bbq factor (=Alpha)
                //-----------------------------------------------
                double        spotVol = 1e0 / bsMrParam.AlphaNorm();

                for (idx=0; idx<=volDiag.mVolDates.size()-1; idx++)
                {
                            factVol.push_back(spotVol);
                }
        }


        // Clear factVol
        //factVolCurves.clear();

        volDates = volDiag.mVolDates;
        factVolCurves.insert(factVolCurves.end(), bsDim, factVol);


        GtoFreeTCurve(zcCurve);
        VnfmFree(vnfmData);
        delete [] zcDates;
        delete [] zcRates;

        if (refRate) delete refRate;
        if (bsRate) delete bsRate;

    }
    catch (KFailure) {
        GtoFreeTCurve(zcCurve);
        VnfmFree(vnfmData);
        delete [] zcDates;
        delete [] zcRates;
        if (refRate) delete refRate;
        if (bsRate) delete bsRate;
        throw KFailure("%s: failed.\n", routine);
    }
}

//--------------------------------------------------------------------
//
/*---------------------------------------------------------------
 * Calibration of credit spread spot volatilities given
 * the CDS volatilities using 2Q VNFM.
 * Input vols are Black-Scholes implied vols unless crSmilePar.mNumIter<0
 * in which case 
 */
void        
DppCalibCRSpreadVolDiag(
    const KMrParam           &crMrParam,        // (I) CR tree info
    const KMrParam           &irMrParam,        // (I) IR tree info
    double                   ircrCorr,          // (I) IR/CR/ correlation
    const KSmileParam        &crSmilePar,       // (I) credit smile info
    TDate                    baseDate,          // (I) base date of ZC
    const KZCurve            &irCurve,          // (I) IR reference zero curve
    const KZCurve            &crCurve,          // (I) credit zero curve
    const KVolDiag           &crVolDiag,        // (I) credit vol info
    double                   recovery,          // (I) recovery rate for BM CDS
    KVector(TDate)           &volDates,         // (I/O) spd spot vol dates
    KVector(KVector(double)) &factVolCurves)    // (I/O) Comb IR/CR spot vols
{
static  char    routine[] = "DppCalibCRSpreadVolDiag";

    int             idx;
    KVector(double) factVol;
    KVector(int)    calibratePoint();
    KVector(TDate)  bmMaturityDates;
    KVolDiag        volDiag = crVolDiag;    // Make a copy of credit vol info
    TDate           volDateSt;

    if (crMrParam.mNumFact!=1)
    {
        throw KFailure("%s: credit VNFM calibration requires CR diffusion"
            "dimension 1, not %d",
            routine, crMrParam.mNumFact);
    }
    else if (irMrParam.mNumFact!=1) 
    {
        throw KFailure("%s: credit VNFM calibration requires IR diffusion"
            "dimension 1, not %d",
            routine, irMrParam.mNumFact);
    }

 try{

        /*=====================================================================
         * INITIALIZATION CHECKS
         *===================================================================*/

        // Check base date of IR diffuse zc curve and vol dates.
        // We may need to offset by 1 the input vector
        // if value date is NOT present in the vol dates
        // And insert a dummy value corresponding to valueDate
        // in mVolMats, mVolFreqs, and mVolRates.
        volDateSt = volDiag.mVolDates[0];
        if(baseDate < volDateSt){
                volDiag.mVolDates.insert(volDiag.mVolDates.begin(),
                                         baseDate);
                volDiag.mVolRates.insert(volDiag.mVolRates.begin(),
                                         volDiag.mVolRates[0]);
                volDiag.mVolFreqs.insert(volDiag.mVolFreqs.begin(),
                                         volDiag.mVolFreqs[0]);
                volDiag.mVolMats.insert(volDiag.mVolMats.begin(),
                                        volDiag.mVolMats[0]);
        }

        /* Remove any credit spot vol curves that may be there */
        while (factVolCurves.size()>1) factVolCurves.pop_back();

        if (crSmilePar.mNumIter<0) {

            /*=================================================================
             * IF NEGATIVE numIter, THEN JUST TAKE GIVEN VOLS AS SPOT VOLS
             * but scale them by the weight factor first.
             *===============================================================*/
            double alpha2 = 1e0 / crMrParam.AlphaNorm();
            alpha2 = alpha2 * alpha2;
            for (idx=0; idx<=volDates.size()-1; idx++)
            {

                factVol.push_back(MAX(volDiag.VolInterp(volDates[idx]) * alpha2, MIN_VOL));
            
            }
            factVolCurves.insert(factVolCurves.end(), factVol);

        } else {

            // insert extra dates corresponding to the credit dates,
            // if they are not there already from the rates vol dates
            // note that vols are not interpolated: they are kept flat
            // so as not to spoil the IR calibration that has already happened
            KVector(TDate)::iterator vdIdx = volDates.begin();
            KVector(double)::iterator vvIdx = factVolCurves[0].begin();
            for (idx=1; idx<volDiag.mVolDates.size(); idx++) {
                TDate cDate = volDiag.mVolDates[idx];
                while (vdIdx!=volDates.end() && *vdIdx<cDate) {
                    vdIdx++;
                    vvIdx++;
                }
                if (vdIdx==volDates.end()) {
                    // no match in whole vector - add at end
                    // and extend value
                    vdIdx = volDates.insert(vdIdx, cDate);
                    vvIdx = factVolCurves[0].insert(vvIdx,factVolCurves[0][factVolCurves[0].size()-1]);
                } else if (*vdIdx==cDate) {
                    // date there already - nothing to do
                } else {
                    // no match, but not at end: always insert
                    // next IR vol, since want to keep
                    // IR vol step function unchanged
                    vdIdx = volDates.insert(vdIdx, cDate);
                    vvIdx = factVolCurves[0].insert(vvIdx, *vvIdx);

                    /* OLD CODE THAT USED TO INTERPOLATE */
                    //TDate thisDate = *vdIdx;
                    //if (vdIdx!=volDates.begin()) vdIdx--;
                    //TDate lastDate = *vdIdx;
                    //if (thisDate==lastDate) {
                    //    vdIdx = volDates.insert(vdIdx, cDate);
                    //    vvIdx = factVolCurves[0].insert(vvIdx, *vvIdx);
                    //} else {
                        // interpolate variance
                    //    double thisVol = *vvIdx;
                     //   vvIdx--;
                    //    double lastVol = *vvIdx;
                    //    vdIdx++;
                    //    vvIdx++;
                    //    vdIdx = volDates.insert(vdIdx, cDate);
                    //    vvIdx = factVolCurves[0].insert(vvIdx, 
                    //        sqrt((thisVol*thisVol*(cDate-lastDate) + lastVol*lastVol*(thisDate-cDate))
                    //            /(thisDate-lastDate)));
                    //}
                }
            }
        
            /*=================================================================
             * DO VNFM APPROXIMATION
             *===============================================================*/
            int numVolCR = volDiag.mVolRates.size();
            KVector(int) calibratePoint(numVolCR,1); // calibrate all points
            // the first point is actually the spot date, so there should be no vol
            // calibrated, or you get errors
            calibratePoint[0] = 0;
            TDateInterval cdsFreq;
            // assume all benchmarks have same frequency
            if (FAILURE==GtoMakeDateInterval(12/volDiag.mVolFreqs[0],'M',&cdsFreq))
            {
                throw KFailure("%s: failed to create BM CDS TDateInterval.\n", routine);
            }

            // build the vector of BM maturity dates
            for (idx=0; idx<volDiag.mVolMats.size(); idx++) {
                bmMaturityDates.push_back((TDate)(volDiag.mVolDates[idx]+(TDate)(365*volDiag.mVolMats[idx])));
                
                // don't try to calibrate if no
                //if (bmMaturityDates.back()<=volDiag.mVolDates[idx]) calibratePoint[idx] = 0;
                //cout << idx << "\tvolDate=" << GtoFormatDate(volDiag.mVolDates[idx])
                //    << "\t\tBM Maturity=" << GtoFormatDate(bmMaturityDates.back()) << endl;
            }
            KVector(double) parSpreadAry(volDiag.mVolMats.size());
            KVector(double) annuityAry(volDiag.mVolMats.size());


            /*=================================================================
             * COMPUTE CREDIT VNFM APPROXIMATIONS
             *===============================================================*/
            if(FAILURE == CrxCreditVNFMBootstrapSpotVols (
                baseDate,
                numVolCR,      /**<(I) Num. of input credit vol points     */
                &volDiag.mVolDates[0],     /**<(I) Credit vol dates                    */
                &volDiag.mVolRates[0],         /**<(I/O) Input BS vols, output spot if cal */
                &calibratePoint[0],     /**<(I/O) Input/output 1 if should/did cal  */
                TRUE,    /**<(I) If TRUE, skip & go on if fail BM cal*/
                cdsFreq,       /**<(I) Underlying BM CDS fee frequency     */
                GTO_ACT_360,        /**<(I) Underlying BM CDS fee day-count - HARDCODED! */
                recovery,      /**<(I) Recovery Rate                       */
                &volDiag.mVolDates[0], /**<(I) Underlying BM CDS start dates       */
                &bmMaturityDates[0],   /**<(I) Underlying BM CDS maturity dates    */
                crMrParam.mBeta[0],        /**<(I) Credit mean-reversion               */
                crMrParam.mBackboneQ,
                crMrParam.AlphaNorm(),
                crCurve,      /**<(I) Credit clean spread curve           */
                crSmilePar.mQ1,      /**<(I) Left Credit 2Q smile (1=lognormal)  */
                crSmilePar.mQ2,     /**<(I) Right Credit 2Q smile (1=lognormal) */
                crSmilePar.mQF,
                ircrCorr,      /**<(I) IR/CR spot-vol correlation          */
                volDates.size(),      /**<(I) Num. of IR spot vol points          */
                &volDates[0],     /**<(I) IR vol dates                        */
                &factVolCurves[0][0],         /**<(I) IR spot vols                        */
                irMrParam.mBeta[0],        /**<(I) IR mean-reversion                   */
                irMrParam.mBackboneQ,
                irMrParam.AlphaNorm(),
                irCurve,       /**<(I) IR zero curve                       */
                &parSpreadAry[0],  /**<(O)Output BM par spreads, where cal     */
                &annuityAry[0]    /**<(O)Output BM annuity value, where cal   */
                ))
            {
            
                throw KFailure("%s: failed - credit VNFM bootstrap failed.\n", routine);
            }

            /*=================================================================
             * Insert into factVolCurves.
             *===============================================================*/
            for (idx=0; idx<=volDates.size()-1; idx++) {
                factVol.push_back(MAX(volDiag.VolInterp(volDates[idx]), MIN_VOL));
            }
            factVolCurves.insert(factVolCurves.end(), factVol);
        }
    }

    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
    

}

