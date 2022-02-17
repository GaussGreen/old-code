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
#include "kvoldat.h"

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
        
            if (mBeta[k] < 0e0) 
                throw KFailure("Mean reversion #%d (%f) < 0.\n",
                                 k, mBeta[k]);
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
        if (sizeModParams != modParams.size())
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


        //
        // Tree Parameters
        //
 
        ASSERT_OR_THROW(treeParams.size() >= 3);
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

        TmxFlag  = FALSE; // default is FIX tree
        OWFlag   = FALSE; // Overwritten smiles
}

//---------------------------------------------------------------
//

KSmileParam::~KSmileParam()
{
    for (KVector(double*)::iterator iter = mMQSmile.begin();
         iter != mMQSmile.end(); iter++)
         {
             FREE((*iter));
         }
}



//---------------------------------------------------------------
//
bool
KSmileParam::IsValid()
{

 try{
        double        QMid, q;

        if (!IsTmx())
        {
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

        }
        return true;        

    }
    catch(KFailure) {
        return false;
    }

}

//---------------------------------------------------------------
//
bool
KSmileParam::IsTmx() const
{
    if (TmxFlag == TRUE)
        return true;
    else
        return false;
}

//---------------------------------------------------------------
//
bool
KSmileParam::IsOW() const
{
    if (OWFlag == TRUE)
        return true;
    else
        return false;
}

//---------------------------------------------------------------
//

istream&
operator>>(istream& is, KSmileParam& object)
{
static        char        routine[] = "operator>>(istream&, KSmileParam&)";

    try {

        if( object.IsTmx())
        {
            object.mSkew    = getDouble(is, "mSkew");
            object.mVov     = getDouble(is, "mVov");
            object.mBBV     = getDouble(is, "mBBV");
            object.mBBR     = getDouble(is, "mBBR");
            object.mDeltaL  = getDouble(is, "mDeltaL");
            object.mTauL    = getDouble(is, "mTauL");
            object.mDeltaR  = getDouble(is, "mDeltaR");
            object.mTauR    = getDouble(is, "mTauR");
            object.mCutoff  = getDouble(is, "mCutoff");
            object.mNCK     = getDouble(is, "mNCK");
            object.mNumIter = getInt(is, "mNumIter");
        }
        else
        {
            object.mQ1 = getDouble(is, "mQ1");
            object.mQ2 = getDouble(is, "mQ2");
            object.mQF = getDouble(is, "mQF");
            object.mNumIter = getInt(is, "mNumIter");
        }
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
static        char        routine[] = "operator<<(ostream&, const KSmileParam&)";
    try
    {
        if (object.IsTmx())
        {
            KVector(TDate)::const_iterator itDate;
            KVector(double*)::const_iterator itSmile = object.mMQSmile.begin();

            for ( itDate = object.mSmlDates.begin();
                  itDate != object.mSmlDates.end();
                  itDate++)
                  {
                    os << "Smile Date: " 
                         << DrlTDate2DrDate(*itDate) << endl;
                    if (itSmile == object.mMQSmile.end())
                        throw KFailure("%s failed: No MQ smile parameters "
                                       "assigned on %ld", 
                                       routine, 
                                       DrlTDate2DrDate(*itDate));

                    os << "#mSkew:\t" << (*itSmile)[0] ;
                    os << "\t#mVov:\t" << (*itSmile)[1];
                    os << "\t#mBBV:\t" << (*itSmile)[2];
                    os << "\t#mBBR:\t" << (*itSmile)[3];
                    os << "\t#mDeltaL:\t" << (*itSmile)[4];
                    os << "\t#mTauL:\t" << (*itSmile)[5];
                    os << "\t#mDeltaR:\t" << (*itSmile)[6];
                    os << "\t#mTauR:\t" << (*itSmile)[7];
                    os << "\t#mCutoff:\t" << object.mCutoff;
                    os << "\t#mNCK:\t" << object.mNCK;
                    os << "\t#mNumIter:\t" << object.mNumIter << endl;

                    itSmile++;
                  }


            for ( itDate = object.mLiqBmkDates.begin();
                  itDate != object.mLiqBmkDates.end();
                  itDate++)
              {
                    os << " Liquid Benchmark Date : " 
                       << DrlTDate2DrDate(*itDate) << endl;
              }
        }
        else
        {
            os << "#mQ1:\n\t" << object.mQ1 << '\n';
            os << "#mQ2:\n\t" << object.mQ2 << '\n';
            os << "#mQF:\n\t" << object.mQF << '\n';
            os << "#mNumIter:\n\t" << object.mNumIter << '\n';
        }
        return os.flush();
    }
    catch (...) {
        throw KFailure("%s: failed.\n", routine);
    }
}


//---------------------------------------------------------------
//

ostream&
KSmileParam::YacctionWrite(ostream& os, int indent)
{
        if (IsTmx())
        {
            os << "#TreeType\n\t" << "TMX" << endl;
            os << "#DistType\n\t" << "nil" << endl;
        }
        else
        {
            os << "#TreeType\n\t" << "FIX" << endl;
            os << "#DistType\n\t";
            os << 1e0 - mQ1 << "\t"
               << 1e0 - mQ2 << "\t"
               << mQF << "\t"
               << mNumIter << endl;
        }

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
                throw KFailure("%s: defaults (nil) supported.\n",
                        routine);
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
//

void
KSmileParam::ReadDrw_x(
        istream        &is,                        // (I) deal input stream
        TDrWrapperData *drWrapData)        // (I) wrapper market data 
{
static        char        routine[] = "KSmileParam::ReadDrw_x";
        char    vs[128];

    try {

        // Distribution type:
        // L | N
        // Fix: Q-, Q+, shft, nbBlIter
        // Tmx: skew, vov, bbv, bbr, deltaL, tauL, deltaR, tauR,
        //      cutoff, NCK, nbIter
        // we need it once we have the # factors.
        strcpy(vs, getString(is, "tree and distribution type"));

        if (EQUAL_STR(vs, "TMX"))
        {
            TmxFlag = TRUE;
            strcpy(vs, getString(is, "distribution type"));
        }
        else if (EQUAL_STR(vs, "FIX"))
        {
            TmxFlag = FALSE;
            strcpy(vs, getString(is, "distribution type"));
        }


        // Now overwrite distrib type
        if (EQUAL_STR(vs, "L")) 
        {
            // lognormal 
            if( TmxFlag == TRUE)
            {
                mSkew   = 0e0;   
                mVov    = 0e0;
                mBBV    = 0e0;
                mBBR    = 0e0;
                mDeltaL = 0.25;
                mTauL   = 0e0;
                mDeltaR = 0.75;
                mTauR   = 0e0;
                mCutoff = 8.;
                mNCK    = 1070433665.;
                mNumIter= 15.;
                OWFlag  = TRUE;

            }
            else
            {
                mQ1 = 1e0;
                mQ2 = 1e0;
                mQF = 0e0;
                mNumIter = 0;
            }
        } else if (EQUAL_STR(vs, "N")) 
        {
            if( TmxFlag == TRUE)
            {
                mSkew   = 1e0;   
                mVov    = 0e0;
                mBBV    = 0e0;
                mBBR    = 0e0;
                mDeltaL = 0.25;
                mTauL   = 1e0;
                mDeltaR = 0.75;
                mTauR   = 1e0;
                mCutoff = 8.;
                mNCK    = 1070433665.;
                mNumIter= 15.;
                OWFlag  = TRUE;
            }
            else
            {
                mQ1 = 0e0;
                mQ2 = 0e0;
                mQF = 0e0;
                mNumIter = 0;
            }
        } 
        else if (EQUAL_STR(vs, "nil")) 
        {
                if (!IsTmx())
                {
                    throw KFailure("%s: defaults (nil) not supported.\n",
                        routine);
                }
                mNumIter= 15.;
                OWFlag = FALSE;

        }
        else 
        {
            if ( IsTmx() )
            {
                DppReadFromString(vs, mSkew);
                mVov    = getDouble(is, "mVov");
                mBBV    = getDouble(is, "mBBV");
                mBBR    = getDouble(is, "mBBR");
                mDeltaL = getDouble(is, "mDeltaL");
                mTauL   = getDouble(is, "mTauL");
                mDeltaR = getDouble(is, "mDeltaR");
                mTauR   = getDouble(is, "mTauR");
                mCutoff = getDouble(is, "mCutoff");
                mNCK    = getDouble(is, "mNCK");
                mNumIter= getInt(is, "mNumIter");
                OWFlag  = TRUE;
            }
            else
            {
                DppReadFromString(vs, mQ1);
                mQ1 = 1e0 - mQ1;
                mQ2 = 1e0 - getDouble(is, "mQ2");
                mQF = getDouble(is, "mQF");
                mNumIter = getInt(is, "mNumIter");
            }
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

        TmxFlag = FALSE;

    }
    catch(KFailure) {
        throw KFailure("%s: failed.\n", routine);
    }
}


void
KSmileParam::ReadMagnet_tmx(const KVector(double) &smileParams,
                            const KVector(TDate)  &bmkLiqDates,
                            const KVector(TDate)  &smlDates)
{
static        char        routine[] = "KSmileParam::ReadMagnet";

    try {
        int i, j;
        int nbSmlDates, nbSmlPars;
        double *MQSml;

        SetSmlDates(smlDates);

        nbSmlDates = mSmlDates.size();
        nbSmlPars  = NBMQSMLPAR * nbSmlDates;

        // Error checking
        //
        if (smileParams.size() < nbSmlPars + 7)
                throw KFailure("%s: need %d  smile inputs (%d).\n",
                                routine, nbSmlPars+7, smileParams.size());
        TmxFlag = TRUE;
        OWFlag  = FALSE;

        for ( i = 0 ; i < nbSmlDates; i++)
        {
            ASSERT_OR_THROW(
                (MQSml = NEW_ARRAY(double, NBMQSMLPAR)) != NULL);
            for (j = 0; j < NBMQSMLPAR; j++)
            {
                MQSml[j] = smileParams[i+j*nbSmlDates];
            }
            mMQSmile.push_back(MQSml);
            MQSml = NULL;
        }

        MQNormT     = smileParams[nbSmlPars];
        MQNck       = smileParams[nbSmlPars+1];
        MQNbDevL    = smileParams[nbSmlPars+2];
        MQNbFwdL    = smileParams[nbSmlPars+3];
        MQBinMapStdL= smileParams[nbSmlPars+4];
        MQBinMapStdR= smileParams[nbSmlPars+5];
        mNumIter    = smileParams[nbSmlPars+6];
        

        for (KVector(TDate)::const_iterator itDate = bmkLiqDates.begin();
             itDate != bmkLiqDates.end();
             itDate++)
             {
                mLiqBmkDates.push_back(*itDate);
             }

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


void
KSmileParam::SetSmlDates(const KVector(TDate) &SmlDate)
{
    static  char routine[] = "KSmileParam::SetSmlDates";
    try
    {
        for (KVector(TDate)::const_iterator itDate = SmlDate.begin();
             itDate != SmlDate.end();
             itDate++)
             {
                mSmlDates.push_back(*itDate);
             }

    }
    catch(KFailure) 
    {
        throw KFailure("%s: failed.\n", routine);
    }

 
}


/*---------------------------------------------------------------
 * TMX MultiQ smile interpolations
 */
void
KSmileParam::SmileDiag(
             const KVolCalibIndex& calibIdx,
             TDrWrapperData *drWrapperData)
{
    static  char routine[] = "KSmileParam::SmileDiag";
    int         i,j;
    TDate      *BmkDates = NULL;
    double    **Smiles = NULL;
 
    try
    {
        if (!IsTmx())
        {
            throw KFailure("%s: FIX tree does not support MQ smile.\n",
                        routine);
        }

        TDate       baseDate = drWrapperData->fToday;
        double      bvMat,
                    calibMatYrs,
                    minMat = 0.25;

        int         nbBmks;
        double     *MQSml = NULL;

        int         calibFinal = calibIdx.IsFinal();
        TDate       calibMatDate = calibIdx.FinalDate();
        TDateInterval   calibMatInt = calibIdx.Tenor();


        if (calibIdx.IsNil())
        {
            throw KFailure("%s: cannot interpolate MQ smile with "
                " nil calibration index. \n", routine);
        }
        else
        {
            KVector(TDate)::const_iterator itDate;
            /* 
             * MultiQ Smile parameters: 
             * Smile[NBMQSMLPAR][VOLDATES]
             * 
             */

            ASSERT_OR_THROW( (Smiles = NEW_ARRAY(
                                      double*,
                                      NBMQSMLPAR)) != NULL );
            for ( i = 0; i < NBMQSMLPAR; i++)
            {
                ASSERT_OR_THROW( (Smiles[i] = NEW_ARRAY( 
                                              double,
                                              mSmlDates.size())) != NULL);
            }
            
        	/* base vol maturity */
	        bvMat = 1e0 / (double) drWrapperData->fBvCurve->fBasis;
	        minMat = bvMat;

            if ((!calibFinal) || (calibMatDate == 0L)) 
            {
                IF_FAILED_THROW( GtoDateIntervalToYears(
                &calibMatInt, &calibMatYrs));
            }
            /*
             * If final calibration or if CMS calib with
             * mat larger than base vol mat, use swaptions
             * Otherwise use base vol.
             */
            if ((calibFinal) ||
                (calibMatYrs > bvMat + 1e-4)) 
            {
                /*
                * Interp from swaption smile
                */
                TSwaptionMatrix2D *SmlIdx = NULL;

                /* If final, but not date, use interval */
                if ((calibFinal) && (calibMatDate == 0L)) 
                {
                    IF_FAILED_THROW( GtoDtFwdAny(
                                    baseDate,
                                    &calibMatInt,
                                    &calibMatDate));
                }

                /* Liquid benchmark dates */
                SmlIdx = (drWrapperData->fSwapSmile)[0];
                nbBmks = SmlIdx->table->matrix->numDim1;
                ASSERT_OR_THROW( (BmkDates = NEW_ARRAY(
                                             TDate, 
                                             nbBmks)) != NULL );

                for ( i = 0; i < nbBmks; i++)
                {
                    IF_FAILED_THROW (GtoTDateAdvanceYears(
                                        baseDate,
                                        SmlIdx->table->dim1Values[i],
                                        &(BmkDates[i])));
                }
               
                /* Interp smiles according to SmlDate*/
                //double tExp, fwdMat;
                for ( i = 0; i < NBMQSMLPAR; i++)
                {

                    SmlIdx = (drWrapperData->fSwapSmile)[i];
                    for ( itDate = mSmlDates.begin(), j = 0;
                          itDate != mSmlDates.end();
                          itDate++, j++)
                      {
/*                          IF_FAILED_THROW( GtoDayCountFraction(
                                           baseDate,
                                           (*itDate),
                                           GTO_ACT_365F,
                                           &tExp));

                          IF_FAILED_THROW( DrlTSwaptionMatrix2DInterpValue(
                                            SmlIdx,
                                            baseDate,
                                            tExp,
                                            calibFinal,
                                            calibMatDate,
                                            calibMatInt,
                                            &(Smiles[i][j]),
                                            &fwdMat,
                                            FALSE));*/

                        IF_FAILED_THROW( DrlTSwaptionMatrix2DDateInterpValue(
                                         SmlIdx,
                                         baseDate,
                                         (*itDate),
                                         calibFinal,
                                         calibMatDate,
                                         calibMatInt,
                                         &(Smiles[i][j]),
                                         FALSE));
                      }
                }
                    
            } 
            else 
            {
                /*
                 * Interp from base smile
                 */
                 TCurve  *bvCurve = NULL;
    
                 bvCurve  = drWrapperData->fBaseSmile[0];
                 nbBmks   = bvCurve->fNumItems;

                 ASSERT_OR_THROW(
                    (BmkDates = NEW_ARRAY(TDate, nbBmks)) != NULL);

                 for (j = 0; j < nbBmks; j++) 
                 {
                     BmkDates[j] = bvCurve->fArray[j].fDate;
                 }

                 for (i = 0; i < NBMQSMLPAR ; i++)
                 {
                    bvCurve  = drWrapperData->fBaseSmile[i];
                    
                    for ( itDate = mSmlDates.begin(), j = 0;
                          itDate != mSmlDates.end();
                          itDate++, j++)
                        {
                            IF_FAILED_THROW ( GtoInterpRate(
                                              (*itDate),
                                              bvCurve,
                                              GTO_LINEAR_INTERP,
                                              &(Smiles[i][j])) );
                        }
                 }
            }
            // Get value date 
            TDate valueDate = drWrapperData->fDiscZcCurve->fBaseDate;
    
            for (i = 0; i < nbBmks; i++)
            {
                if (BmkDates[i] > valueDate+2) 
                {
                    mLiqBmkDates.push_back(BmkDates[i]);
                }
            }

            for (i = 0; i < mSmlDates.size(); i++)
            {
                ASSERT_OR_THROW(
                    (MQSml = NEW_ARRAY(double, NBMQSMLPAR)) != NULL);
                
                if (IsOW())
                {
                    MQSml[0] = mSkew;
                    MQSml[1] = mVov;
                    MQSml[2] = mBBV;
                    MQSml[3] = mBBR;
                    MQSml[4] = mDeltaL;
                    MQSml[5] = mTauL;
                    MQSml[6] = mDeltaR;
                    MQSml[7] = mTauR;
                //    MQSml[8] = mCutoff;
                //    MQSml[9] = mNCK;
                    MQNormT  = mCutoff;
                    MQNck    = mNCK;

                }
                else
                {
                    for ( j = 0; j < NBMQSMLPAR; j++ )
                    {
                        MQSml[j] = Smiles[j][i];                
                    }
//                    MQSml[NBMQSMLPAR]   = drWrapperData->MQNormT;
//                    MQSml[NBMQSMLPAR+1] = drWrapperData->MQNck;
                }
                mMQSmile.push_back(MQSml);
                MQSml = NULL;
            }
        }
        FREE(BmkDates);

        for (i = 0; i < NBMQSMLPAR; i++)
        {
            FREE(Smiles[i]);
        }
        FREE(Smiles);

    }
    catch(KFailure) 
    {
        FREE(BmkDates);

        for (i = 0; i < NBMQSMLPAR; i++)
        {
            FREE(Smiles[i]);
        }
        FREE(Smiles);
        throw KFailure("%s: failed.\n", routine);
    }
}

void    
KSmileParam::CheckLiqBmkDates(
        const KVolDiag &VolDiag)
{
    static  char routine[] = "KSmileParam::SmileDiag";

    try
    {
        if (IsTmx())
        {
            if (mLiqBmkDates.size() < 2)
                throw KFailure("%s failed: Require at least 2 liquid "
                               "benchmark dates\n", routine);

            KVector(TDate)::const_iterator itDate1;
            for (KVector(TDate)::iterator itDate2 = mLiqBmkDates.begin();
                itDate2 != mLiqBmkDates.end();
                ++itDate2)
            {
                for ( itDate1  = VolDiag.mVolDates.begin();
                      itDate1 != VolDiag.mVolDates.end();
                      ++itDate1)
                      {
                          if ( (*itDate1) == (*itDate2) )
                              break;
                      }

                if (itDate1 == VolDiag.mVolDates.end() )
                    throw KFailure("%s failed: The liquid benchmark dates "
                                   " are not in ATM vol grid\n", routine);
   
            }
        }
    }
    catch(KFailure) 
    {
        throw KFailure("%s: failed.\n", routine);
    }
}
