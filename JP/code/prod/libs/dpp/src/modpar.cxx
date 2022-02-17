/****************************************************************
 * Module:        PenGuin
 * Submodule:        
 * File:        modpar.c
 * Function:        
 * test
 * Author:        Christian Daher, David Liu
 * Revision:        $Header$
 ****************************************************************/
#include "kmodpar.h"


#include "kutilios.h"
extern        "C" {
#include "drlio.h"
#include "drlstr.h"
#include "drltime.h"

#include "dritkwrp.h"           // TDrWrapper

};

#define        EQUAL_STR(s1, s2)        (!strcmp(s1, s2))

//===============================================================
//
//        KMrParam
//
//===============================================================

//---------------------------------------------------------------
//

KMrParam::KMrParam()
{

        mNumFact = 0;
        mBeta[0]  = 0.0;
        mAlpha[0] = 1.0;

        mBackboneQ = 1e0;        // Default is lognormal

        mPpy = 4;
        mSmoothFact  = 0;
        mNumStdevCut = 5;
        mCorrAdjFlag = 0;
}

//---------------------------------------------------------------
//

KMrParam::~KMrParam()
{
}


//---------------------------------------------------------------
//
bool
KMrParam::IsValid()
{
        int        k;
        double        res;

#undef                SQR
#define                SQR(x)  ((x)*(x))

  try {
        if(mNumFact == 0)
                return true;
        
        if(mNumFact < 0 ||
           mNumFact > 3) 
                throw KFailure("Number of factors (%d) is either > 3 "
                                "or < 0.\n", mNumFact);        
                
        for (k=0; k<=mNumFact-1; k++) {
            if (IS_ALMOST_ZERO(mBeta[k]))
                mBeta[k] = 1e-8;
        }
                
        if (mNumFact >= 2) {
            for (k=0; k<mNumFact*(mNumFact-1)/2; k++) {
                if (fabs(mRho[k]) > 1e0)
                        throw KFailure("Correlation #%d (%f) must satisfy "
                                       "-1 < rho < 1.\n", k, mRho[k]);

                // To avoid lattice degeneration, round to number close to 1
                //
                mRho[k] = MAX(MIN(mRho[k], 0.99), -0.99);
            }
        }

        if (mNumFact == 3) {
                res = 1e0 - SQR(mRho[0])
                          - SQR(mRho[1])
                          - SQR(mRho[2])
                          + 2e0 * mRho[0]
                                * mRho[1]
                                * mRho[2];
                if (res < 0e0) 
                        throw KFailure("Determinant of correlation matrix "
                                       "(res %f) < 0.\n", res);
        }


        if (mBackboneQ < 0e0 ||
            mBackboneQ > 1e0) 
                throw KFailure("Backbone Q (%f) is out of the range [0,1].\n",
                                mBackboneQ);

        if (mPpy < 1) 
                throw KFailure("PPY (%d) < 1.\n", mPpy);
                
        
        if (mSmoothFact < 0)
                throw KFailure("Smooth factor (%f) < 0.\n", mSmoothFact);
                
        if (mNumStdevCut < 0)
                throw KFailure("Number of stdev cut (%f) < 0.\n", 
                                  mNumStdevCut);

        return true;

    }
    catch (KFailure) {
        return false;
    }

}



//---------------------------------------------------------------
//

istream&
operator>>(istream& is, KMrParam& object)
{
static        char        routine[] = "operator>>(istream&,KMrParam&)";
        int        i;

    try {


        object.mNumFact = getInt(is,"mNumFact");
        if ((object.mNumFact < 0) || (object.mNumFact > 3)) {
                throw KFailure("%s: invalid number of fact %d.\n",
                        routine, object.mNumFact);
        }

        for (i=0; i<=object.mNumFact-1; i++) 
                object.mBeta[i] = getDouble(is, "mBeta");
        for (i=0; i<=object.mNumFact-1; i++) 
                object.mAlpha[i] = getDouble(is, "mAlpha");
        switch (object.mNumFact) {
        case 1:
                break;
        case 2:
                object.mRho[0] = getDouble(is, "mRho[0]");
                break;
        case 3:
                object.mRho[0] = getDouble(is, "mRho[0]");
                object.mRho[1] = getDouble(is, "mRho[1]");
                object.mRho[2] = getDouble(is, "mRho[2]");
                break;
        }

        object.mPpy = getInt(is, "mPpy");
        object.mSmoothFact= getInt(is, "mSmoothFact");
        object.mNumStdevCut = getInt(is, "mNumStdevCut");

        return is;
    }
    catch (KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}


//---------------------------------------------------------------
//

ostream&
operator<<(ostream& os, const KMrParam& object)
{
        int        i;
        
        os << '\n';
        os << "#mNumFact:\n\t" << object.mNumFact << '\n';

        for (i=0; i<=object.mNumFact-1; i++) 
                os << format("#mBeta[%2d]:\n\t%lf\n", i, object.mBeta[i]);
        for (i=0; i<=object.mNumFact-1; i++) 
                os << format("#mAlpha[%2d]:\n\t%lf\n", i, object.mAlpha[i]);

        switch (object.mNumFact) {
        case 1:
                break;
        case 2:
                os << "#mRho[0]:\n\t" << object.mRho[0] << '\n';
                break;
        case 3:
                os << "#mRho[0]:\n\t" << object.mRho[0] << '\n';
                os << "#mRho[1]:\n\t" << object.mRho[1] << '\n';
                os << "#mRho[2]:\n\t" << object.mRho[2] << '\n';
                break;
        }

        os << "#mPpy:\n\t" << object.mPpy << '\n';
        os << "#mBackboneQ:\n\t" << object.mBackboneQ << '\n';
        os << "#mSmoothFact:\n\t" << object.mSmoothFact << '\n';
        os << "#mNumStdevCut:\n\t" << object.mNumStdevCut << '\n';

        return os.flush();
}


//---------------------------------------------------------------
//

ostream&
KMrParam::YacctionWrite(ostream& os, int indent) 
{
        int i;


        os << "# NumFact" << endl;
        os << mNumFact;
        os << endl;

        os << "# Betas" << endl;
        for (i=0; i<mNumFact; i++)
                os << mBeta[i] << "\t";
        os << endl;
 
        os << "# Alphas" << endl;
        for (i=0; i<mNumFact; i++)
                os << mAlpha[i] << "\t";
        os << endl;
 
        os << "# Rhos" << endl;
        if (mNumFact ==1)
                os << "nil";
        else
        {
                for (i=0; i<(mNumFact-1)*mNumFact/2; i++)
                os << mRho[i] << "\t";
        }
        os << endl;

        os << "# IR Backbone Q (0=Lognormal, 1=Normal) " << endl;
        os << 1.0 - mBackboneQ << endl;
 
        os << "# Ppy" << endl;
        os << mPpy << endl;
 
        os << "# Smoothing" << endl;
        os << (mSmoothFact == 1 ? "Y" : "N") << endl;
 
        os << "# NumStdevCut" << endl;
        os << mNumStdevCut << endl;

        return (os);
}




//---------------------------------------------------------------

ostream&
KMrParam::BasisYacctionWrite(ostream& os, int indent) 
{
        os << "# Basis Backbone Q (0=Lognormal, 1=Normal) " << endl;
        os << 1.0 - mBackboneQ << endl;

        return (os);
}



//---------------------------------------------------------------

KMrParam&
KMrParam::operator=(double argument)
{
        switch (mNumFact) {
        case 3:
                mBeta[2]  = argument;
                mAlpha[2] = argument;
                mRho[1] = argument;
                mRho[2] = argument;
        case 2:
                mBeta[1]  = argument;
                mAlpha[1] = argument;
                mRho[0] = argument;
        case 1:
                mBeta[0]  = argument;
                mAlpha[0] = argument;
                break;
        }

        return(*this);
}

//---------------------------------------------------------------

KMrParam&
KMrParam::operator+=(double argument)
{
        switch (mNumFact) {
        case 3:
                mBeta[2]  += argument;
                mAlpha[2] += argument;
                mRho[1] += argument;
                mRho[2] += argument;
        case 2:
                mBeta[1]  += argument;
                mAlpha[1] += argument;
                mRho[0] += argument;
        case 1:
                mBeta[0]  += argument;
                mAlpha[0] += argument;
                break;
        }

        // WE DO NOT ADD THE PPY !

        return(*this);
}

//---------------------------------------------------------------

KMrParam&
KMrParam::operator*=(double argument)
{
        switch (mNumFact) {
        case 3:
                mBeta[2]  *= argument;
                mAlpha[2] *= argument;
                mRho[1] *= argument;
                mRho[2] *= argument;
        case 2:
                mBeta[1]  *= argument;
                mAlpha[1] *= argument;
                mRho[0] *= argument;
        case 1:
                mBeta[0]  *= argument;
                mAlpha[0] *= argument;
                break;
        }

        // WE DO NOT MULTIPLY THE PPY !

        return(*this);
}

//---------------------------------------------------------------

KMrParam&
KMrParam::operator=(const KMrParam& m)
{
        mNumFact = m.mNumFact;


        switch (mNumFact) {
        case 3:
                mBeta[2] = m.mBeta[2];
                mAlpha[2] = m.mAlpha[2];
                mRho[1] = m.mRho[1];
                mRho[2] = m.mRho[2];
        case 2:
                mBeta[1] = m.mBeta[1];
                mAlpha[1] = m.mAlpha[1];
                mRho[0] = m.mRho[0];
        case 1:
                mBeta[0] = m.mBeta[0];
                mAlpha[0] = m.mAlpha[0];
                break;
        }

        mPpy = m.mPpy;
        mSmoothFact = m.mSmoothFact;
        mNumStdevCut = m.mNumStdevCut;
        mCorrAdjFlag = m.mCorrAdjFlag;
        mBackboneQ   = m.mBackboneQ;

        return(*this);
}



//---------------------------------------------------------------

KMrParam&
KMrParam::operator+=(const KMrParam& m)
{
        if (mNumFact != m.mNumFact)
            throw KFailure(
                "KMrParam::operator=(const KMrParam&)"
                " incompatible arguments.\n");


        switch (mNumFact) {
        case 3:
                mBeta[2] += m.mBeta[2];
                mAlpha[2] += m.mAlpha[2];
                mRho[1] += m.mRho[1];
                mRho[2] += m.mRho[2];
        case 2:
                mBeta[1] += m.mBeta[1];
                mAlpha[1] += m.mAlpha[1];
                mRho[0] += m.mRho[0];
        case 1:
                mBeta[0] += m.mBeta[0];
                mAlpha[0] += m.mAlpha[0];
                break;
        }
        return(*this);
}


//--------------------------------------------------------------
//

void
KMrParam::ReadDrw(
        istream        &is,                        // (I) deal input stream
        TDrWrapperData *drWrapData)        // (I) wrapper market data 
{
static        char        routine[] = "KMrParam::ReadDrw";

        char    vs[128];
        int        nFact, idxF;

    try {

        // Number of factors
        //
        nFact = getInt(is, "number of factors");
        switch (nFact) {
        case 1:
                mNumFact   = 1;
                mBeta[0]   = drWrapData->f1Beta;
                mAlpha[0]  = drWrapData->f1Weight;
                mPpy       = drWrapData->f1Ppy;
                break;
        case 2:
                mNumFact   = 2;
                mBeta[0]   = drWrapData->f2Beta1;
                mBeta[1]   = drWrapData->f2Beta2;
                mAlpha[0]  = drWrapData->f2Weight1;
                mAlpha[1]  = drWrapData->f2Weight2;
                mRho[0]    = drWrapData->f2Corr12;
                mPpy       = drWrapData->f2Ppy;
                break;
        case 3:
                mNumFact   = 3;
                mBeta[0]   = drWrapData->f3Beta1;
                mBeta[1]   = drWrapData->f3Beta2;
                mBeta[2]   = drWrapData->f3Beta3;
                mAlpha[0]  = drWrapData->f3Weight1;
                mAlpha[1]  = drWrapData->f3Weight2;
                mAlpha[2]  = drWrapData->f3Weight3;
                mRho[0]    = drWrapData->f3Corr12;
                mRho[1]    = drWrapData->f3Corr13;
                mRho[2]    = drWrapData->f3Corr23;
                mPpy       = drWrapData->f3Ppy;
                break;
        default:
                throw KFailure("%s: bad num factors (%d).\n",
                        routine, nFact);
        }

        //
        // Other defaults
        //
        mBackboneQ  = 1e0;        // Log normal !!!
        mSmoothFact = 1;        //
        mNumStdevCut = 5;
        mCorrAdjFlag = 0;

        //
        // Now reads overwrites from data file
        //

        // betas
        strcpy(vs, getString(is, "mean reversions"));
        if (!EQUAL_STR(vs, "nil")) {
                DppReadFromString(vs, mBeta[0]);
                for (idxF=1; idxF<nFact; idxF++) {
                        mBeta[idxF] = getDouble(is);
                }
        }

        // weights
        strcpy(vs, getString(is, "volatilities"));
        if (!EQUAL_STR(vs, "nil")) {
                DppReadFromString(vs, mAlpha[0]);
                for (idxF=1; idxF<nFact; idxF++) {
                        mAlpha[idxF] = getDouble(is);
                }
        }

        // correlations
        strcpy(vs, getString(is, "correlations"));
        if (!EQUAL_STR(vs, "nil")) {
                DppReadFromString(vs, mRho[0]);
                switch (nFact) {
                case 2:
                        break;
                case 3:
                        mRho[1] = getDouble(is);
                        mRho[2] = getDouble(is);
                        break;
                }
        }

        // Backbone
        strcpy(vs, getString(is, "backboneQ"));
        if (!EQUAL_STR(vs, "nil")) {
                DppReadFromString(vs, mBackboneQ); 
                mBackboneQ = 1.0 - mBackboneQ;
        }

        
        strcpy(vs, getString(is, "ppy"));
        if (!EQUAL_STR(vs, "nil")) {
                DppReadFromString(vs, mPpy);
        }

        strcpy(vs, getString(is, "smoothing"));
        if (toupper(vs[0]) == 'Y') {
                mSmoothFact = 1;
        } else if (toupper(vs[0]) == 'N') {
                mSmoothFact = 0;
        } else {
                DppReadFromString(vs, mSmoothFact);
        }

        strcpy(vs, getString(is, "num stdev"));
        if (!EQUAL_STR(vs, "nil")) {
                DppReadFromString(vs, mNumStdevCut);
        }


    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}


//--------------------------------------------------------------
//

void
KMrParam::IRMAWMrParam(
        istream        &is,                        // (I) deal input stream
        CRX_INPUT *crxInput)        // (I) wrapper market data 
{
static        char        routine[] = "KMrParam::IRMAWMrParam";

        char    vs[128];
        int        i;

        IR_INPUT  *irInput = NULL;

    try {


        // IR diffuse curve
        //
        irInput = &(crxInput->irInput[0]);

        // Number of factors
        //
        mNumFact = irInput->NbFactor;

        // Only 1+1 is supported
        if (mNumFact > 1)
        {
                throw KFailure("%s: only 1-factor IR + 1-factor Credit "
                               "is supported.\n",
                               routine);
        } 

        // fill the MR and weight factors
        for (i=0; i<mNumFact; i++)
        {
                mBeta[i]   = irInput->Beta[i];
                mAlpha[i]  = irInput->Alpha[i];
        }

        // fill the factor correlation
        for (i=0; i<mNumFact*(mNumFact-1)/2; i++)
        {
                mRho[i]    = irInput->Rho[i];
        }


        //
        // Backbone
        //
        mBackboneQ  = 1.0 - irInput->BackboneCoeff;


        //
        // Now reads additional model params from data file
        //

        
        strcpy(vs, getString(is, "ppy"));
        if (!EQUAL_STR(vs, "nil")) {
                DppReadFromString(vs, mPpy);
        }

        strcpy(vs, getString(is, "smoothing"));
        if (toupper(vs[0]) == 'Y') {
                mSmoothFact = 1;
        } else if (toupper(vs[0]) == 'N') {
                mSmoothFact = 0;
        } else {
                DppReadFromString(vs, mSmoothFact);
        }

        strcpy(vs, getString(is, "num stdev"));
        if (!EQUAL_STR(vs, "nil")) {
                DppReadFromString(vs, mNumStdevCut);
        }

    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}


//--------------------------------------------------------------
//

void
KMrParam::IRMAWMrParam(
        istream        &is,                        // (I) deal input stream
        BS_INPUT       *bsInput)                   // (I) wrapper market data 
{
static        char        routine[] = "KMrParam::IRMAWMrParam";

        char    vs[128];
        int        i;

        IR_INPUT  *irInput = NULL;

    try {

/*
        strcpy(vs, getString(is, "mNumFact"));
        DppReadFromString(vs, mNumFact);*/

        // IR diffuse curve
        //
        irInput = &(bsInput->irInput[0]);

        // Number of factors
        //
        mNumFact = irInput->NbFactor;

        // fill the MR and weight factors
        for (i=0; i<mNumFact; i++)
        {
                mBeta[i]   = irInput->Beta[i];
                mAlpha[i]  = irInput->Alpha[i];
        }

        // fill the factor correlation
        for (i=0; i<mNumFact*(mNumFact-1)/2; i++)
        {
                mRho[i]    = irInput->Rho[i];
        }


        //
        // Backbone
        //
        mBackboneQ  = 1.0 - irInput->BackboneCoeff;


        //
        // Now set additional model params 
        //

        mPpy = bsInput->PPY;
        mNumStdevCut = bsInput->nbCutoff;

        strcpy(vs, getString(is, "smoothing"));
        if (toupper(vs[0]) == 'Y') {
                mSmoothFact = 1;
        } else if (toupper(vs[0]) == 'N') {
                mSmoothFact = 0;
        } else {
                DppReadFromString(vs, mSmoothFact);
        }


    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}




//--------------------------------------------------------------
//

void
KMrParam::CRMAWMrParam(
        istream        &is,                        // (I) deal input stream
        CRX_INPUT *crxInput)        // (I) wrapper market data 
{
static        char        routine[] = "KMrParam::CRMAWMrParam";


        CR_INPUT  *crInput = NULL;

    try {


        // CR curve
        //
        crInput = &(crxInput->crInput[0]);

        // Number of factors
        //
        mNumFact = crInput->NbFactor;

        // Only 1+1 is supported
        if (mNumFact > 1)
        {
                throw KFailure("%s: only 1-factor IR + 1-factor Credit "
                               "is supported.\n",
                               routine);
        } 

        if (mNumFact != 1)
            throw KFailure("%s: credit factor (%d) != 1.\n", routine, mNumFact);
    

        mBeta[0]   = crInput->Beta;
        mAlpha[0]  = crInput->Alpha;


        //
        // Backbone
        //
        mBackboneQ  = 1.0 - crInput->BackboneCoeff;


        //
        // Other tree params (NOT used, overwritten by IR tree params)
        //
        mPpy = 24;
        mSmoothFact = 1;
        mNumStdevCut = 5;

        mCorrAdjFlag = crInput->CorrAdjFlag;
        if (mCorrAdjFlag < 0 || mCorrAdjFlag > 2)
            throw KFailure("%s: Ir/CR CorrAdjFlag only be (0,1,2)\n.", routine);

    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}

//
void
KMrParam::SPMAWMrParam(
        istream        &is,                        // (I) deal input stream
        BS_INPUT       *bsInput)                   // (I) wrapper market data 
{
static        char        routine[] = "KMrParam::SPMAWMrParam";


        SP_INPUT  *spInput = NULL;

    try {


        // SP curve
        //
        spInput = &(bsInput->spInput[0]);

        // Number of factors
        //
        mNumFact = spInput->NbFactor;

        // Only 1+1 is supported
        if (mNumFact > 1)
        {
                throw KFailure("%s: only 1-factor IR + 1-factor spread "
                               "is supported.\n",
                               routine);
        } 

        if (mNumFact != 1)
            throw KFailure("%s: spread factor (%d) != 1.\n", routine, mNumFact);
    

        mBeta[0]   = spInput->Beta;
        mAlpha[0]  = spInput->Alpha;


        //
        // Backbone
        //
        mBackboneQ  = 1.0 - spInput->BackboneCoeff;


        //
        // Other tree params (NOT used, overwritten by IR tree params)
        //
        mPpy = 24;
        mSmoothFact = 1;
        mNumStdevCut = 5;


    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}


//---------------------------------------------------------------
//

ostream&
KMrParam::MAWYacctionWrite(ostream& os, int indent) 
{

        os << "# Ppy" << endl;
        os << mPpy << endl;
 
        os << "# Smoothing" << endl;
        os << (mSmoothFact == 1 ? "Y" : "N") << endl;
 
        os << "# NumStdevCut" << endl;
        os << mNumStdevCut << endl;

        return (os);
}


//--------------------------------------------------------------


void
KMrParam::ReadLil(const double *mrParamsL)
{
static        char        routine[] = "KMrParam::ReadLil";
        int        i, esz;

    try {
        if (ARGSIZE(mrParamsL) < 5) {
                throw KFailure("%s: expects at least array length 5 "
                        "(got %d).\n", routine, ARGSIZE(mrParamsL));
        }
        mPpy         = (int) mrParamsL[1];
        mNumStdevCut = mrParamsL[2];
        mSmoothFact  = mrParamsL[3];

        // THE CONVENTION FOR I/O IS backQ=0 FOR LOGNORMAL, backQ=1 FOR NORMAL.
        // INTERNALLY, THE "GOOD" CONVENTION HOLDS, I.E.
        // backQ=0 FOR NORMAL, backQ=1 FOR LOGNORMAL.
        // And for Vnfm, backQ=0 for Lognormal, and backQ=0.5 for normal.
        // i.e. (1-backQ)*0.5
        mBackboneQ   = 1e0 - mrParamsL[4];
        mNumFact     = (int) mrParamsL[5];

        if (mNumFact < 0) {
                throw KFailure("%s: invalid number of factors (%d).\n",
                        routine, mNumFact);
        }
        if (mNumFact > 3) {
                throw KFailure("%s: invalid number of factors (%d).\n",
                        routine, mNumFact);
        }
        if (mNumFact == 0) return;

        esz = 5 + 2*mNumFact + mNumFact*(mNumFact-1)/2;
        if (ARGSIZE(mrParamsL) != esz) {
            throw KFailure("%s: expects at array length %d+5 with %d factors "
                "(got total %d).\n", routine, esz-5, mNumFact, ARGSIZE(mrParamsL));
        }
        for (i=0; i<mNumFact; i++) {
                mBeta[i]  = mrParamsL[1+5+i];
                mAlpha[i] = mrParamsL[1+5+mNumFact+i];
        }
        for (i=0; i<mNumFact*(mNumFact-1)/2; i++) {
                mRho[i] = mrParamsL[1+5+2*mNumFact+i];
        }


    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}




//--------------------------------------------------------------
//
void
KMrParam::ReadMagnet(
        const KVector(double) &modParams,
        const KVector(int)    &treeParams)
{
static        char        routine[] = "KMrParam::ReadMagnet";
        int        idx, numFact, sizeModParams;

try {

        // Error checking
        //
        if (modParams.size() < 4)
                throw KFailure("%s: expects at least array length 4 for "
                                "modParams. (got %d).\n",
                                routine, modParams.size());
 
        numFact = (int) modParams[0];
        if (numFact > 3 || numFact < 1)
                throw KFailure("%s: invalid number of factors (%d).\n",
                                routine, numFact);
 
        mNumFact = numFact;
 
        //
        // Expected size = 1 (numFact) + numFact (betas) + numFact (alphas)
        // +  numFact*(numFact-1)/2 (rhos) + 1 (backbone)
        //
        sizeModParams = 1 + 2*numFact + numFact*(numFact-1)/2 + 1;

        bool hasIrCrCorrAdjFlag = false;
        if (sizeModParams == modParams.size())
            hasIrCrCorrAdjFlag = false;
        else if (sizeModParams + 1 == modParams.size())
            hasIrCrCorrAdjFlag = true;
        else 
            throw KFailure("%s: size of input model parameters "
                                "inconsistent with factor number.\n.", routine);
        //
        // Assign the parameters
        //
        for (idx = 0; idx < numFact; idx++)
        {
                mBeta[idx]  = modParams[1+idx];
                mAlpha[idx] = modParams[1+numFact+idx];
        }
 
        for (idx = 0; idx < numFact*(numFact-1)/2; idx++)
                mRho[idx] = modParams[1+2*numFact+idx];
 
        mBackboneQ = 1e0 - modParams[sizeModParams-1];

        if (hasIrCrCorrAdjFlag)
        {
            mCorrAdjFlag  = modParams[modParams.size()-1];
            if (mCorrAdjFlag < 0 || mCorrAdjFlag > 2)
                throw KFailure("%s: Ir/CR CorrAdjFlag only be (0,1,2)\n.", routine);
        }
        else
            mCorrAdjFlag = 0;

        //
        // Tree Parameters
        //
 
        ASSERT_OR_THROW(treeParams.size() == 3);
        mPpy          = treeParams[0];
        mSmoothFact   = treeParams[1];
        mNumStdevCut  = treeParams[2];

    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }


}



//--------------------------------------------------------------
//

double
KMrParam::AlphaNorm() const
{
static        char        routine[] = "KMrParam::AlphaNorm";
        double norm = 0e0;
        for (int idxF = 0; idxF<mNumFact; idxF++)
                norm += mAlpha[idxF] * mAlpha[idxF];
        if (norm < DBL_EPSILON)
                throw KFailure("%s: total alpha too small %g.\n",
                        routine, norm);
 
        norm = sqrt(norm);
        return (norm);
}


//--------------------------------------------------------------
//


KMrParam Correlate(
        KMrParam &mrParam1,        // (I)
        KMrParam &mrParam2,        // (I)
        double corr)                // (I)
{
static        char        routine[]= "Correlate(KMrParam&,KMrParam&,double)";
        KMrParam        newParam;

    try {
        // Trivial case: one of the dim is 0
        //
        if (mrParam1.mNumFact == 0) return (mrParam2);
        if (mrParam2.mNumFact == 0) return (mrParam1);

        // Correlation
        if(fabs(corr) > 1e0)
                throw KFailure("%s: correlation (%f) between two assets must "
                               "be within [-1, 1].\n", 
                                routine,
                                corr);

        //
        // Non trivial case
        //
        if ((mrParam1.mNumFact == 1) &&
            (mrParam2.mNumFact == 1)) {

                newParam.mNumFact = 2;
                newParam.mBeta[0]  = mrParam1.mBeta[0];
                newParam.mBeta[1]  = mrParam2.mBeta[0];
                newParam.mAlpha[0] = mrParam1.mAlpha[0];
                newParam.mAlpha[1] = mrParam2.mAlpha[0];

                newParam.mRho[0] = corr;

                newParam.mBackboneQ = mrParam1.mBackboneQ;
                newParam.mPpy = mrParam1.mPpy;
                newParam.mSmoothFact = mrParam1.mSmoothFact;
                newParam.mNumStdevCut = mrParam1.mNumStdevCut;

        } else
        if ((mrParam1.mNumFact == 1) &&
            (mrParam2.mNumFact == 2)) {
                throw KFailure("%s: sorry, not implemented dim %d x dim %d.\n",
                        routine, mrParam1.mNumFact,
                        mrParam2.mNumFact);


        } else {
                throw KFailure("%s: cannot handle dim %d x dim %d.\n",
                        routine, mrParam1.mNumFact,
                        mrParam2.mNumFact);
        }

        return (newParam);
    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}


//===============================================================
//
//        KSmileParam
//
//===============================================================



//---------------------------------------------------------------
//

KSmileParam::KSmileParam()
{
        mQ1 = 1e0;
        mQ2 = 1e0;
        mQF = 0e0;
        mNumIter = -1; // default is no VNFM or CET
}

//---------------------------------------------------------------
//

KSmileParam::~KSmileParam()
{
}



//---------------------------------------------------------------
//
bool
KSmileParam::IsValid()
{

 try{
        double        QMid, q;

        /* I have changed this so that can use mNumIter=-1 as a way of
         * signalling that should interpret vols as spot vols, and not
         * BS vols, rather like "nil" specification currently         
         * Charles Morcom 10/26/05 */
        if (mNumIter < -1)
                throw KFailure("Number of Black iteration (%d) < -1.\n",
                                 mNumIter);

        //
        // Check on mQF non-singular, i.e. mQF != -1,
        // QMid * mQF != -1, and QRight * mQF > -1, etc.
        //

        if (IS_ALMOST_ZERO(1. + mQF))
            throw KFailure("Forward shift is -1, singular in the mapping "
                           "function.\n");

        QMid = (mQ1 + mQ2) / 2.;

        if (IS_ALMOST_ZERO(1. + QMid * mQF))                
            throw KFailure("Smile parameters are singular in 2Q mapping. "
                          "1 + (qL +qR)/2 * Fsh is zero.\n");

        //
        // Ensure that calculation of Zt in CalcDiscDiff
        // would not fail.
        //
        q = (mQF > 0) ? mQ2 : mQ1;
        if ((fabs(q) > 1.e-5) &&
            (q * mQF < -1.))
            throw KFailure("Invalid forward shift: %lf. "
                           "The smile distribution at X=0 does not "        
                           "correspond to the forward.\n",
                            mQF);

        return true;        

    }
    catch(KFailure) {
        return false;
    }

}



//---------------------------------------------------------------
//

istream&
operator>>(istream& is, KSmileParam& object)
{
static        char        routine[] = "operator>>(istream&, KSmileParam&)";

    try {

        object.mQ1 = getDouble(is, "mQ1");
        object.mQ2 = getDouble(is, "mQ2");
        object.mQF = getDouble(is, "mQF");
        object.mNumIter = getInt(is, "mNumIter");
        return is;
    }
    catch (...) {
        throw KFailure("%s: failed.\n", routine);
    }
}


//---------------------------------------------------------------
//

ostream&
operator<<(ostream& os, const KSmileParam& object)
{
        os << "#mQ1:\n\t" << object.mQ1 << '\n';
        os << "#mQ2:\n\t" << object.mQ2 << '\n';
        os << "#mQF:\n\t" << object.mQF << '\n';
        os << "#mNumIter:\n\t" << object.mNumIter << '\n';
        return os.flush();
}


//---------------------------------------------------------------
//

ostream&
KSmileParam::YacctionWrite(ostream& os, int indent)
{
        os << "#DistType" << endl;
        os << 1e0 - mQ1 << "\t"
           << 1e0 - mQ2 << "\t"
           << mQF << "\t"
           << mNumIter << endl;

        return (os);
}


//---------------------------------------------------------------
//

ostream&
KSmileParam::BasisYacctionWrite(ostream& os, int indent)
{
        os << "# Basis lower smile q" << endl;
        os << 1e0 - mQ1 << endl;
        os << "# Basis higher smile q" << endl;
        os << 1e0 - mQ2 << endl;
        os << "# Basis smile forward shift" << endl;
        os << mQF << endl;

        return (os);
}





//---------------------------------------------------------------
//


void
KSmileParam::ReadDrw(
        istream        &is,                        // (I) deal input stream
        TDrWrapperData *drWrapData)        // (I) wrapper market data 
{
static        char        routine[] = "KSmileParam::ReadDrw";
        char    vs[128];

    try {

        // Distribution type:
        // L | N
        // Q-, Q+, shft, nbBlIter
        // we need it once we have the # factors.
        strcpy(vs, getString(is, "distribution type"));

        // Now overwrite distrib type
        if (EQUAL_STR(vs, "L")) {
                // lognormal 
                mQ1 = 1e0;
                mQ2 = 1e0;
                mQF = 0e0;
                mNumIter = 0;
        } else if (EQUAL_STR(vs, "N")) {
                // normal
                mQ1 = 0e0;
                mQ2 = 0e0;
                mQF = 0e0;
                mNumIter = 0; 
        } else if (EQUAL_STR(vs, "nil")) {
                // Get defaults.
                mQ1 = 1e0 - drWrapData->LQ;
                mQ2 = 1e0 - drWrapData->RQ;
                mQF = drWrapData->FwdShift;
                mNumIter = drWrapData->nbIter;
        } else {
                DppReadFromString(vs, mQ1);
                mQ1 = 1e0 - mQ1;
                mQ2 = 1e0 - getDouble(is, "mQ2");
                mQF = getDouble(is, "mQF");
                mNumIter = getInt(is, "mNumIter");
        }


    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }

}


//---------------------------------------------------------------



void
KSmileParam::IRMAWSmileParam(
        istream            &is,        // (I) deal input stream
        CRX_INPUT   *crxInput)        // (I) wrapper market data 
{
static        char        routine[] = "KSmileParam::IRMAWSmileParam";

    IR_INPUT  *irInput = NULL;

    try {


        // IR diffuse curve
        //
        irInput = &(crxInput->irInput[0]);

        mQ1 = irInput->QLeft;
        mQ2 = irInput->QRight;
        mQF = irInput->FwdShift;

        // read CET iteration from deal file
        mNumIter = getInt(is, "mNumIter");


    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }

}

void
KSmileParam::IRMAWSmileParam(
        istream            &is,        // (I) deal input stream
        BS_INPUT           *bsInput)   // (I) wrapper market data 
{
static        char        routine[] = "KSmileParam::IRMAWSmileParam";

    IR_INPUT  *irInput = NULL;

    try {


        // IR diffuse curve
        //
        irInput = &(bsInput->irInput[0]);

        mQ1 = irInput->QLeft;
        mQ2 = irInput->QRight;
        mQF = irInput->FwdShift;

        mNumIter = bsInput->nbCetIter;


    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }

}


void
KSmileParam::CRMAWSmileParam(
        istream            &is,        // (I) deal input stream
        CRX_INPUT   *crxInput)        // (I) wrapper market data 
{
static        char        routine[] = "KSmileParam::CRMAWSmileParam";

    CR_INPUT  *crInput = NULL;

    try {


    // CR curve
    //
    crInput = &(crxInput->crInput[0]);

        mQ1 = crInput->QLeft;
        mQ2 = crInput->QRight;
        mQF = crInput->FwdShift;

        mNumIter = -1; // this is read later

    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }

}

void
KSmileParam::SPMAWSmileParam(
        istream            &is,        // (I) deal input stream
        BS_INPUT           *bsInput)   // (I) wrapper market data 
{
static        char        routine[] = "KSmileParam::CRMAWSmileParam";

    SP_INPUT  *spInput = NULL;

    try 
    {
        // SP curve
        //
        spInput = &(bsInput->spInput[0]);

        mQ1 = spInput->QLeft;
        mQ2 = spInput->QRight;
        mQF = spInput->FwdShift;

        mNumIter = 0; // No CET on spread

    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }

}

//---------------------------------------------------------------
//

ostream&
KSmileParam::MAWYacctionWrite(ostream& os, int indent) 
{

        os << "# NumCETIter" << endl;
        os << mNumIter << endl;
 
        return (os);
}



//--------------------------------------------------------------


void
KSmileParam::ReadLil(const double *smileParamsL)
{
static        char        routine[] = "KSmileParam::ReadLil";

    try {
        if (ARGSIZE(smileParamsL) != 4) {
                throw KFailure("%s: expects at array length 4 "
                        "(got %d).\n", routine, ARGSIZE(smileParamsL));
        }

        // THE CONVENTION FOR I/O IS Q=0 FOR LOGNORMAL, Q=1 FOR NORMAL.
        // INTERNALLY, THE "GOOD" CONVENTION HOLDS, I.E.
        // Q=0 FOR NORMAL, Q=1 FOR LOGNORMAL.
        mQ1      = 1e0 - smileParamsL[1];
        mQ2      = 1e0 - smileParamsL[2];
        mQF      = smileParamsL[3];
        mNumIter = (int) smileParamsL[4];

    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}



//---------------------------------------------------------------

void
KSmileParam::ReadMagnet(const KVector(double) &smileParams)
{
static        char        routine[] = "KSmileParam::ReadMagnet";

    try {

        // Error checking
        //
        if (smileParams.size() != 4)
                throw KFailure("%s: need 4 smile inputs (%d).\n",
                                routine, smileParams.size());
 
        //
        // User smile conventions:  1=Normal, 0=Lognormal
        // Model smile conventions: 0=Normal, 1=Lognormal
        //
        mQ1      = 1.0 - smileParams[0];
        mQ2      = 1.0 - smileParams[1];
        mQF      = smileParams[2];
        mNumIter = (int)smileParams[3];

    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}





//---------------------------------------------------------------

KSmileParam&
KSmileParam::operator=(double argument)
{
        mQ1 = argument;
        mQ2 = argument;
        mQF = argument;
        return(*this);
}

//---------------------------------------------------------------

KSmileParam&
KSmileParam::operator+=(double argument)
{
        mQ1 += argument;
        mQ2 += argument;
        mQF += argument;
        return(*this);
}

//---------------------------------------------------------------

KSmileParam&
KSmileParam::operator*=(double argument)
{
        mQ1 *= argument;
        mQ2 *= argument;
        mQF *= argument;
        return(*this);
}

//---------------------------------------------------------------

KSmileParam&
KSmileParam::operator=(const KSmileParam& m)
{
        mQ1 = m.mQ1;
        mQ2 = m.mQ2;
        mQF = m.mQF;
        mNumIter = m.mNumIter;
        return(*this);
}



//---------------------------------------------------------------

KSmileParam&
KSmileParam::operator+=(const KSmileParam& m)
{
        mQ1 += m.mQ1;
        mQ2 += m.mQ2;
        mQF += m.mQF;
        return(*this);
}


