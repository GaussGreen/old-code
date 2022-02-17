/****************************************************************
 * Module:        PenGuin
 * Submodule:        
 * File:        modpar.c
 * Function:        
 * Author:        Christian Daher - David Liu
 * Revision:        $Header$
 ****************************************************************/
#include "kmktcrv.h"
#include "kmodpar.h"
#include "kvoldat.h"


#include "kutilios.h"
extern        "C" {
#include "drlio.h"
#include "drlstr.h"
#include "drlinter.h"
#include "drltime.h"

#include "dritkwrp.h"           // TDrWrapper

};

#define        EQUAL_STR(s1, s2)        (!strcmp(s1, s2))

static  TCurve* GetTZCurve(char zcType, TDrWrapperData *drWrapData);
static  void    GetMAWTZCurveInfo(
                        string     zcName,     /* (I) curve name */
                        CRX_INPUT  *crxInput,  /* (I) wrapper data input */
                        char       *zcType,    /* (O) curve type, "ir" "cr", */
                        int        *zcId);     /* (O) curve ID   */
static  void    GetMAWTZCurveInfo(
                        string     zcName,     /* (I) curve name */
                        BS_INPUT   *bsInput,   /* (I) wrapper data input */
                        char       *zcType,    /* (O) curve type, "ir" "bs", */
                        int        *zcId);     /* (O) curve ID   */

//===============================================================
// 
// KBasisTreeModelParam methods.
//
//===============================================================



//--------------------------------------------------------------

KMarketCurves::KMarketCurves()
{
        mIsBasis = FALSE;

        mBSType = SUB_SPREAD;
        mLiborCVName = String("Default");
        mBSDiscCVName = String("Default");
        mBSDelayShift = 0e0;

        mBasisDCC = KDayCc("ACT/365");
        mLiborDCC = KDayCc("ACT/360");


        // Credit
        mIsCredit = FALSE;
        mRecovery = 0e0;
        mIRDiscCVName = String("Default");

}


//--------------------------------------------------------------

KMarketCurves::~KMarketCurves()
{
}



//--------------------------------------------------------------
//
bool
KMarketCurves::IsValid()
{
    try {

        // Check for duplicate curve types
        //
        KVector(int) copyCVType = mCVTypes;
        sort(copyCVType.begin(), copyCVType.end());
        KVector(int)::iterator iterEnd
                = unique(copyCVType.begin(), copyCVType.end());
        if (iterEnd != copyCVType.end())
                throw KFailure("Duplicate curve types in market curves.\n");
 
        // Check today <= value date
        for (KVector(int)::iterator iter=mCVTypes.begin();
                iter!=mCVTypes.end(); ++iter)
        {
                if(GetValueDate(*iter)<mToday)
                        throw KFailure("Curve %d value date (%s) "
                                     "< today (%s).\n",
                                *iter,
                                GtoFormatDate(GetValueDate(*iter)),
                                GtoFormatDate(mToday));
 
        } 

        // Check for basis
        if(IsBasis())
        {
                if (mBSType != SUB_SPREAD &&
                        mBSType != ADD_SPREAD &&
                    mBSType != PER_SPREAD ) 
                        throw KFailure("Unknown basis type (%d). "
                                       "must be either 0 for SUB_SPREAD, "
                                       "1 for PERCENT, or 2 for ADD_SPREAD.\n",
                                        mBSType);

                if(!IS_ALMOST_ZERO(mBSDelayShift))
                        throw KFailure("Only zero basis delay shift (%f) "
                                "is supported.\n",
                                 mBSDelayShift); 

        }

        return true; 

    }
    catch (KFailure)
    {
        return false;
    }

}



//--------------------------------------------------------------
//

void
KMarketCurves::InsertZc(
        const KZCurve &zcCurve,                // (I) zc curve
        int   zcType,                        // (I) zc type (KV_DIFF, etc.)
        TDate valueDate,                // (I) zc value date
        const String& zcName)                // (I) zc name
{
static        char        routine[] = "KMarketCurves::InsertZc";

    try {

        // Check type
        switch(zcType) {
        case KV_DET:
        case KV_DIFF:
        case KV_IDX1:
        case KV_IDX2:
        case KV_BASIS:
        case KV_SPREAD:
        case KV_PAR_SPREAD:
        case KV_CREDIT_RISKY:
                break;
        default:
                throw KFailure("%s: invalid curve type %d for curve `%s'.\n",
                        routine, zcType, zcName.c_str());
        }

        // Check already inserted
        if (find(mCVTypes.begin(), mCVTypes.end(), zcType) != mCVTypes.end()) {
                throw KFailure("%s: curve type %d (`%s') already inserted.\n",
                        routine, zcType, zcName.c_str());
                //continue;
        }


        // Insert zero curve data
        //
        mCVTypes.push_back(zcType);
        mCV.insert(KMap(int,KZCurve)::value_type(zcType, zcCurve));
        mCVNames.insert(KMap(int,String)::value_type(zcType, zcName));
        mValueDates.insert(KMap(int, TDate)::value_type(
                                zcType, valueDate));


    }
    catch (KFailure)
    {
        throw KFailure("%s: failed.\n", routine);
    }

}



//--------------------------------------------------------------
//

void
KMarketCurves::Insert(
        const KZCurve &zcCurve,                // (I) zc curve
        int zcType,                        // (I) zc type (KV_DIFF, etc.)
        TDate valueDate,                // (I) zc value date
        const String& zcName)                // (I) zc name
{
static        char        routine[] = "KMarketCurves::Insert";

    try {

        // Check type
        switch(zcType) {
        case KV_DET:
        case KV_DIFF:
        case KV_IDX1:
        case KV_IDX2:
                break;
        case KV_BASIS:
        case KV_SPREAD:
        case KV_PAR_SPREAD:
        case KV_CREDIT_RISKY:
                break;
                throw KFailure("%s: curve `%s' is a basis/credit curve"
                        " (use InsertBasis() or InsertCredit() instead).\n",
                        routine, zcName.c_str());
        default:
                throw KFailure("%s: invalid curve type %d for curve `%s'.\n",
                        routine, zcType, zcName.c_str());
        }


        // Insert zero curve
        InsertZc(zcCurve, zcType, valueDate, zcName);

    }
    catch (KFailure)
    {
        throw KFailure("%s: failed.\n", routine);
    }

}




//--------------------------------------------------------------
//


void
KMarketCurves::InsertBasis(
        const KZCurve &zcCurve,         // (I) zc curve
        int zcType,                     // (I) zc type
        TDate valueDate,                // (I) zc value date
        const String& zcName,           // (I) zc name
        KSpd bsType,                    // (I) SUB_SPREAD/PER_SPREAD/ADD_SPREAD
        const String& bsLiborCVName,    // (I) Libor curve name
        const String& bsDiscCVName,     // (I) Basis disc curve name
        double bsDelayShift,            // (I) Reset delay shift
        const KDayCc& bsDCC,            // (I) Basis DCC 
        const KDayCc& liborDCC)         // (I) Libor DCC 
                
{
static        char        routine[] = "KMarketCurves::InsertBasis";

    try {
        // Check type
        switch(zcType) {
        case KV_BASIS:
        case KV_SPREAD:
        case KV_PAR_SPREAD:
                break;
        case KV_DET:
        case KV_DIFF:
        case KV_IDX1:
        case KV_IDX2:
                throw KFailure("%s: curve `%s' is not a basis curve"
                        " (use Insert() instead).\n",
                        routine, zcName.c_str());
        default:
                throw KFailure("%s: invalid curve type %d for curve `%s'.\n",
                        routine, zcType, zcName.c_str());
        }

        if (!IS_ALMOST_ZERO(bsDelayShift))
                throw KFailure("%s: basis deley shift (%f) != 0!\n",
                                routine,
                                bsDelayShift);

        // Insert zero curve
        InsertZc(zcCurve, zcType, valueDate, zcName);

        // Extra basis info
        mIsBasis = true;
        mBSType = bsType;
        mLiborCVName = bsLiborCVName;
        mBSDiscCVName = bsDiscCVName;        
        
        if (bsDelayShift >= 0e0 && bsDelayShift <= 1e0)
                mBSDelayShift = bsDelayShift;
        else
                throw KFailure("%s: invalid delay shift (%f) for basis curve. "
                               "Must be within range [0, 1].\n",
                               routine, bsDelayShift);

        mBasisDCC = bsDCC;
        mLiborDCC = liborDCC;
    }
    catch (KFailure)
    {
        throw KFailure("%s: failed.\n", routine);
    }

}



//--------------------------------------------------------------
bool            
KMarketCurves::IsBasis()
{
static        char        routine[] = "KMarketCurves::IsBasis";
        //----------------------------------------------
        // Turn basis flag ON if basis curve
        // or spread is involved in the environment
        // 
        //----------------------------------------------
 
        KMap(int, KZCurve)::iterator itBS   = mCV.find(KV_BASIS);
        KMap(int, KZCurve)::iterator itSprd = mCV.find(KV_SPREAD);
        KMap(int, KZCurve)::iterator itParSprd
                                        = mCV.find(KV_PAR_SPREAD);
        if ((itBS      != mCV.end()) ||
            (itSprd    != mCV.end()) ||
            (itParSprd != mCV.end()))
        {
            // Credit can not co-exist with basis
            //
            KMap(int, KZCurve)::iterator itCR   = mCV.find(KV_CREDIT_RISKY);
            if (itCR != mCV.end()) 
                throw KFailure("%s: invalid curve dynamics (6=credit) in "
                               "basis env.\n", routine);
            else            
                mIsBasis =true;
        }
 
 
        // The basis curve name should always be present in the tree
        // in order to compute the basis index rate, even if only
        // basis (par)spreads are given in the product.
        //
        if (mIsBasis && (itBS == mCV.end()))
        {
                mCVTypes.push_back(KV_BASIS);

                // Need to find a basis curve name which
                // has not been used yet.  
                // Check if name "BASIS_FWD" has been used
                //
                if (itSprd != mCV.end() &&
                    mCVNames[KV_SPREAD] == String("BASIS_FWD"))
                {
                        throw KFailure("%s: name BASIS_FWD is reserved "
                                "for basis diffuse curve (type 3).  Please "
                                "choose another name for the spread curve.\n",
                                routine);
                }
                else if (itParSprd != mCV.end() &&
                         mCVNames[KV_PAR_SPREAD] == String("BASIS_FWD"))
                {
                        throw KFailure("%s: name BASIS_FWD is reserved "
                                "for basis diffuse curve (type 3).  Please "
                                "choose another name for the par sprd curve.\n",
                                routine);
                }
                else
                        mCVNames.insert(
                            KMap(int, String)::value_type(KV_BASIS, 
                                                          String("BASIS_FWD")));
                
 
                if (itSprd != mCV.end())
                {
                        mCV.insert(
                                KMap(int, KZCurve)::value_type(KV_BASIS,
                                (*itSprd).second));
                        InsertValueDate(
                                KV_BASIS,
                                GetValueDate(KV_SPREAD));
                }
                else    // par spread curve is present
                {
                        mCV.insert(
                        KMap(int, KZCurve)::value_type(KV_BASIS,
                                (*itParSprd).second));
                        InsertValueDate(
                                KV_BASIS,
                                GetValueDate(KV_PAR_SPREAD));
                }
 
        }

        return (mIsBasis);

}




//--------------------------------------------------------------
//


void
KMarketCurves::InsertCredit(
        const KZCurve &zcCurve,             // (I) zc curve
        int   zcType,                       // (I) zc type
        TDate valueDate,                    // (I) zc value date
        const String& zcName,               // (I) zc name
        const String& irDiscountName,       // (I) IR disc curve name
        double    recovery)                 // (I) Recovery
                
{
static        char        routine[] = "KMarketCurves::InsertCredit";

    try {
        // Check type
        if(zcType != KV_CREDIT_RISKY) {
                throw KFailure("%s: curve `%s' is not a credit curve"
                        " (use Insert() instead).\n",
                        routine, zcName.c_str());
        }

        // Insert zero curve
        InsertZc(zcCurve, zcType, valueDate, zcName);

        // Extra basis info
        mIsBasis      = true;
        mIRDiscCVName = irDiscountName;
        mRecovery     = recovery;
        
    }
    catch (KFailure)
    {
        throw KFailure("%s: failed.\n", routine);
    }

}




//--------------------------------------------------------------
bool            
KMarketCurves::IsCredit()
{
static        char        routine[] = "KMarketCurves::IsCredit";

        //----------------------------------------------
        // Turn credit flag ON if Credit curve
        // is involved in the environment
        // 
        //----------------------------------------------
 
        KMap(int, KZCurve)::iterator itCR   = mCV.find(KV_CREDIT_RISKY);
        if (itCR != mCV.end()) 
        {
            // Basis can not co-exist with credit
            //
            KMap(int, KZCurve)::iterator itBS   = mCV.find(KV_BASIS);
            if (itBS != mCV.end()) 
                throw KFailure("%s: invalid curve dynamics (3=basis) in "
                               "credit env.\n", routine);
            else            
                mIsCredit =true;
        }
 

        return (mIsCredit);

}

 

//--------------------------------------------------------------


ostream&
operator<<(ostream& os, const KMarketCurves& object)
{

        //
        // Zero curves
        //
        os << "MARKET CURVES DATA:" << endl;
        os << "TODAY:\n" << GtoFormatDate(object.mToday) << endl;

        os << "NUMBER OF CURVES:\n" <<
                object.mCVTypes.size() << endl;

        for (KVector(int)::const_iterator itCV=object.mCVTypes.begin();
                                    itCV!=object.mCVTypes.end();
                                    ++itCV) 
        {
                KMap(int, String)::const_iterator 
                                iterName = object.mCVNames.find(*itCV);
                os << "CURVE NAME: `" << iterName->second << "'\n";
                os << "CURVE TYPE: " ;
                switch (*itCV) {
                case KV_DET:
                        os << "Deterministic curve    ";
                        break;
                case KV_DIFF:
                        os << "Diffuse                ";
                        break;
                case KV_IDX1:
                        os << "Deterministic spread 1 ";
                        break;
                case KV_IDX2:
                        os << "Deterministic spread 2 ";
                        break;
                case KV_BASIS:
                        os << "Basis Rate             ";
                        break;
                case KV_SPREAD:
                        os << "Basis Spread           ";
                        break;
                case KV_PAR_SPREAD:
                        os << "Basis Par Spread       ";
                        break;
                case KV_CREDIT_RISKY:
                        os << "credit                 ";
                        break;
                }
                os << endl;
                
                KMap(int, String)::const_iterator iterFile 
                                = object.mZCFiles.find(*itCV);
                if (iterFile != object.mZCFiles.end())
                        os << "CURVE FILE: `" << iterFile->second << "'" << endl;
                
                os << "CURVE DETAILS: " << endl;

                KMap(int, KZCurve)::const_iterator iterCV 
                                                = object.mCV.find(*itCV);
                os << iterCV->second << endl;
        }

        //
        // Extra basis information
        //

        if (object.mIsBasis) {
                if (object.mBSType == SUB_SPREAD)
                        os << "#mBSType:\n\t" << "SUB SPREAD" << '\n';
                else if (object.mBSType == ADD_SPREAD)
                        os << "#mBSType:\n\t" << "ADD SPREAD" << '\n';
                else
                        os << "#mBSType:\n\t" << "PERCENT SPREAD" << '\n';
        
                os << "#mLiborCVName:\n\t" << object.mLiborCVName << '\n';
                os << "#mBSDiscCVName:\n\t" << object.mBSDiscCVName << '\n';
                os << "#mBasisDCC:\n\t" << object.mBasisDCC << '\n';
                os << "#mLiborDCC:\n\t" << object.mLiborDCC << '\n';
                os << "#mBSDelayShift:\n\t" << object.mBSDelayShift << '\n';
        }        
        else if (object.mIsCredit) {
                os << "#mIRDiscCVName:\n\t" << object.mIRDiscCVName << '\n';
                os << "#mRecovery:    \n\t" << object.mRecovery << '\n';
        } else {
                os << "#Basis factors\n\tNIL" << endl;
        }

        return(os.flush());
}


//--------------------------------------------------------------
//
ostream& 
KMarketCurves::DrWYacctionWrite( ostream& os, int indent)
{

        String zcFiles[4] = {"Zero.dat",
                             "Disczero.dat",
                             "Riskzero.dat",
                             "Riskzero.dat"}; // basis curve file

        ASSERT_OR_THROW(mCVTypes.size() <= 3);
        os << "# Number of Curves ( <= 3 ) " << endl;
        os << mCVTypes.size() << endl;
 
        os << "# ZC    Type    Name (Type 0=diff; 1=idx1; 2=idx2; 3=bs; 6=credit)"<<endl;
 
        for (KVector(int)::const_iterator itCV=mCVTypes.begin();
                        itCV!=mCVTypes.end(); ++itCV)
        {
            KMap(int, String)::const_iterator iterName
                                        = mCVNames.find(*itCV);

            if (iterName != mCVNames.end()) 
                    os  << zcFiles[*itCV] << "\t"
                    << *itCV << "\t"
                    << iterName->second << endl;
        }

        return (os);
}



//--------------------------------------------------------------
//
ostream& 
KMarketCurves::BasisYacctionWrite( ostream& os, int indent)
{
        if (mIsBasis)        
        {
                os << "# Basis Type ('S'ubtractive; 'A'dditive, 'P'ercent)" << endl;
                if (mBSType == SUB_SPREAD)
                    os << "S" << endl;
                else if (mBSType == ADD_SPREAD)
                    os << "A" << endl;
                else
                    os << "P" << endl;
 
                os << "# Reference Libor index curve name" << endl;
                os << mLiborCVName << endl;
 
                os << "# Reference discount curve name" << endl;
                os << mBSDiscCVName << endl;

                os << "# Basis rate day count convention" << endl;
                os << mBasisDCC << endl;

                os << "# Libor rate day count convention" << endl;
                os << mLiborDCC << endl;
        }

        return (os);
}


//--------------------------------------------------------------
//
ostream& 
KMarketCurves::MAWCreditYacctionWrite( ostream& os, int indent)
{
static  char    routine[] = "KMarketCurves::MAWCreditYacctionWrite";

        char    fileName[256];

        static  char irZCFile[]  = "ir_curve";
        static  char crZCFile[]  = "cr_zcurve";

        os << "# Number of Curves ( <= 4 ) " << endl;
        os << mCVTypes.size() << endl;
 
        os << "# ZC    Type    Name (Type 0=diff; 1=idx1; 2=idx2; 6=credit)"<<endl;

        for (KVector(int)::const_iterator itCV=mCVTypes.begin();
                        itCV!=mCVTypes.end(); ++itCV)
        {
            KMap(int, String)::const_iterator iterName
                                        = mCVNames.find(*itCV);

            if (*itCV <= KV_IDX2)            // IR curve
                sprintf(fileName, "%s%d_0.dat", irZCFile, *itCV); 
            else if (*itCV == KV_CREDIT_RISKY)  // Single name CR curve
                sprintf(fileName, "%s_0.dat", crZCFile); 
            else
                throw KFailure("%s: invalid curve dynamics (%d). "
                               "in credit environment.\n",
                                routine, *itCV);

            if (iterName != mCVNames.end()) 
                    os  << fileName << "\t"
                    << *itCV << "\t"
                    << iterName->second << endl;
        }

        return (os);
}


//--------------------------------------------------------------
//
 
TDate
KMarketCurves::GetValueDate(int curveIdx)
{
static  char    routine[] = "KMarketCurves::GetValueDate";
 
        KMap(int, TDate)::iterator itDate = mValueDates.find(curveIdx);
 
        if (itDate == mValueDates.end())
                throw KFailure("%s: invalid curve index (%d). "
                               "Value date for curve [%d] does not exist.\n",
                                routine, curveIdx, curveIdx);
 
        return (*itDate).second;
 
}


//---------------------------------------------------------------

KZCurve&
KMarketCurves::GetDiffuse()
{
static        char        routine[] = "KMarketCurves::GetDiffuse";
        KMap(int, KZCurve)::iterator it = mCV.find(KV_DIFF);
        if (it == mCV.end()) {
                throw KFailure("%s: diffuse zero curve [%d] "
                        "is not available.\n", routine, KV_DIFF);
        }
        return (*it).second;
}


//---------------------------------------------------------------
// GetIRDiffuse should replace the previous GetDiffuse for
// the clarity of diffusion (we have 2 diffusion process)
//

KZCurve&
KMarketCurves::GetIRDiffuse()
{
static        char        routine[] = "KMarketCurves::GetIRDiffuse";
        KMap(int, KZCurve)::iterator it = mCV.find(KV_DIFF);
        if (it == mCV.end()) {
                throw KFailure("%s: diffuse zero curve [%d] "
                        "is not available.\n", routine, KV_DIFF);
        }
        return (*it).second;
}


//---------------------------------------------------------------
// GetIRDiffuse should replace the previous GetDiffuse for
// the clarity of diffusion (we have 2 diffusion process)
//

KZCurve&
KMarketCurves::GetIRDiscount()
{
static        char        routine[] = "KMarketCurves::GetIRDiscount";
        
        int    irDiscountIdx;

        irDiscountIdx = GetCurveType(mIRDiscCVName);

        KMap(int, KZCurve)::iterator it = mCV.find(irDiscountIdx);
        if (it == mCV.end()) {
                throw KFailure("%s: diffuse zero curve [%d] "
                        "is not available.\n", routine, KV_DIFF);
        }
        return (*it).second;
}




//---------------------------------------------------------------
int
KMarketCurves::GetCurveType(String &cvName)
{
static  char    routine[] = "KMarketCurves::GetCurveType";

        for(KMap(int, String)::iterator iter = mCVNames.begin();
                iter != mCVNames.end(); ++iter)
        {
                if (cvName == (*iter).second)
                        return (*iter).first; 
        }

        throw KFailure("%s: invalid curve name %s, not found in marketCurve.\n",
                        routine, cvName.c_str());

        return -1;

}



//---------------------------------------------------------------
// Return the zero shift between today and value date
// for the discount curve.
//

double
KMarketCurves::ZeroShift(const String &discCurveName)
{
static        char        routine[] = "KMarketCurves::ZeroShift";

        int        curveIdx;
        bool       isFind = false;
        double     zeroShift;
                
    try {

        for(KVector(int)::iterator iterCV = mCVTypes.begin();
                        iterCV != mCVTypes.end(); ++iterCV)
        {
                curveIdx = *iterCV;
 
                KMap(int, String)::iterator itStr
                                = mCVNames.find(curveIdx);
                KMap(int, KZCurve)::iterator itZCV
                                = mCV.find(curveIdx);

                if(itStr  != mCVNames.end() &&
                   itZCV  != mCV.end()) {
                        // Check if it is the product discount curve
                        // --- convert all to lower case for comparison
                        //     need to be done.
                        //
                        if((*itStr).second == discCurveName)
                        {
                                isFind = true;

                                // Risky discount factor 
                                // on IR value date!!!
                                if (curveIdx == KV_CREDIT_RISKY) 
                                {
                                    int   irDiscCurve 
                                            = GetCurveType(mIRDiscCVName);
                                    TDate irValueDate 
                                            = mValueDates[irDiscCurve];

                                    zeroShift  = GetIRDiscount().ZeroShift();

                                    zeroShift *=
                                         (*itZCV).second.DiscFact(irValueDate);
                                }
                                else
                                    zeroShift = (*itZCV).second.ZeroShift();

                                break;
                        }
                }
                else
                        throw KFailure("%s: incorrect mapping of "
                                       "cvNames[%d], or cv[%d].\n",
                                       routine, curveIdx, curveIdx);
        }


        if (!isFind)
                throw KFailure("%s: product discount curve (%s) is not found "
                               "in the market.\n",
                                routine, discCurveName.c_str());

        return (zeroShift);

    }
    catch (KFailure)
    {
        throw KFailure("%s: failed.\n", routine);
    }
}


//---------------------------------------------------------------
//

bool
KMarketCurves::ReadDrw(
        istream& is,                        // (I) wrapper data stream
        TDrWrapperData *drWrapData)        // (I) wrapper market data 
{
static        char        routine[] = "KMarketCurves::ReadDrw";
                
                        // --- zero curves
        int        i;
        int     numCurves;
        int     cvtype;
	//int	zcInterp = GTO_LINEAR_INTERP;	// default for IR ZC
	int	zcInterp = GTO_FLAT_FORWARDS;	// default for IR ZC
        String  cvname, zcname;
        char    zcname0;
        TCurve  *zcurve = NULL;


    try {
        mIsBasis  = false;        // default
        mIsCredit = false;        // default

        // Today date
        //
        mToday = drWrapData->fToday;


        // Zero curves
        //

        // Number of curves
        numCurves = getInt(is, "number of curves");
 
        // Reading in curve infomation:
        //   zcname: curve name in wrapper (disc, zero, risk, etc.)
        //   cvtype: curve type (diffuse, index, basis, par, etc, )
        //   cvname: curve name in product 
        for (i=0; i<=numCurves-1; ++i)
        {
                // Read curve name in wrapper env
                zcname  = getString(is, format("curve %d", i));
                zcname0  = toupper(*(zcname.begin()));
                if (zcname0 != 'D' &&
                    zcname0 != 'Z' &&
                    zcname0 != 'R' &&
                    zcname0 != 'C' &&
                    zcname0 != 'B')
                        throw KFailure("%s: invalid zero curve name (%s) in "
                                       "row %d.\n"
                                        "Must be one of the following:\n"
                                        "(d)isczero, (z)ero, (r)iskzero, "
                                        "(c)reditzero, (b)asiszero.\n",
                                        routine, zcname.c_str(), i+1);

                // Read curve type
                cvtype = getInt(is, "curve type");
                if (cvtype != KV_DET  &&
                    cvtype != KV_DIFF &&
                    cvtype != KV_IDX1 &&
                    cvtype != KV_IDX2 && 
                    cvtype != KV_BASIS &&
                    cvtype != KV_SPREAD && 
                    cvtype != KV_PAR_SPREAD &&
                    cvtype != KV_CREDIT_RISKY)
                        throw KFailure("%s: invalid curve type %d "
                                        "in row %d. "
                                        "%d=Det, %d=Diff, %d=Idx1, %d=Idx2, "
                                        "%d=Basis, %d=Basis Spread, "
                                        "%d=Basis Par Spread, "
                                        "%d=Credit.\n",
                                        routine, cvtype, i+1,
                                        KV_DET, KV_DIFF, KV_IDX1, KV_IDX2,
                                        KV_BASIS, KV_SPREAD, KV_PAR_SPREAD,
                                        KV_CREDIT_RISKY);
                
                // Check if it is a basis curve 
                if (cvtype >= 3) 
                        mIsBasis = true;

                if (cvtype >= KV_CREDIT_RISKY) 
                        mIsCredit = true;
                        
                // Read curve name in product
                cvname = getString(is, "curve name");
                if (cvname.empty())
                        throw KFailure("%s: can not have empty curve name "
                                       "in row %d.\n",
                                        routine, i+1);        

                mCVTypes.push_back(cvtype);
                mCVNames.insert(KMap(int, String)::value_type(cvtype, cvname));


                // Assign zero curves
                //

                if ((zcurve = GetTZCurve(zcname0, drWrapData)) != NULL)
                {
                        mValueDates.insert(KMap(int, TDate)::value_type(
                                cvtype, zcurve->fBaseDate));
                        mCV.insert(KMap(int, KZCurve)::value_type(cvtype,
                                      KZCurve(zcurve, mToday, zcInterp)));
                }
                else
                        throw KFailure("%s: zero curve name (%s) in "
                                       "row %d is NOT available in "
                                       "the environment.\n",
                                        routine, zcname.c_str(), i+1);


                // Assign names of zc files
                //
                switch (zcname0) {
                case 'D': 
                        mZCFiles.insert(KMap(int, String)::value_type(
                                                cvtype,
                                                "disczero.dat"));
                        break;

                case 'Z':
                        mZCFiles.insert(KMap(int, String)::value_type(
                                                cvtype,
                                                "zero.dat"));
                        break;

                case 'R':
                        mZCFiles.insert(KMap(int, String)::value_type(
                                                cvtype,
                                                "riskzero.dat"));

                        break;
                
                case 'B':        
                        mZCFiles.insert(KMap(int, String)::value_type(
                                                cvtype,
                                                "basiszero.dat"));
                        break;
                case 'C':        
                        mZCFiles.insert(KMap(int, String)::value_type(
                                                cvtype,
                                                "creditzero.dat"));
                        break;
                }
        }

        return (mIsBasis);
    }
    catch (KFailure)
    {
        throw KFailure("%s: failed.\n", routine);
    }
}

//---------------------------------------------------------------
//

bool
KMarketCurves::ReadDrwMAW(
        istream&  is,                // (I) wrapper data stream
        TDate     today,             // (I) Today's date
        BS_INPUT *bsInput)           // (I) MAW market data 
{
static        char        routine[] = "KMarketCurves::ReadDrwMAW";
                        // --- zero curves
        int     i;
        int     numCurves;
        int     cvtype;
	    
//        int	zcInterp = GTO_LINEAR_INTERP;
    	int	zcInterp = GTO_FLAT_FORWARDS;	// default for IR ZC
        String  cvname, zcName;
        char    zcType = 'i';
        int     zcId   = 0;
        char    fileName[256];
        
        TCurve  *zcurve = NULL;

        static  char irZCFile[]  = "ir_curve";
        static  char bsZCFile[]  = "bs_zcurve";


    try {
        mIsBasis = false;        // default


        // Today's date
        mToday = today;


        // Zero curves
        //

        // Number of curves
        numCurves = getInt(is, "number of curves");
 
        // Reading in curve infomation:
        //   zcname: curve name in wrapper (ir0, ir1, ir2, bs0, etc.)
        //   cvtype: curve type (diffuse, index, credit, etc, )
        //   cvname: curve name in product 
        for (i=0; i<=numCurves-1; ++i)
        {
                // Read curve name in wrapper env
                zcName  = getString(is, format("curve %d", i));

                GetMAWTZCurveInfo(zcName,
                                  bsInput,
                                  &zcType,
                                  &zcId); 

                switch (zcType) {
                case 'i':    // IR curve family (single currency only)
                    // Check that IR curve is present
                    if (bsInput->nbIRInput == 0)
                    {
                        throw KFailure("%s: no IR curve present in the "
                                       "multi-asset wrapper env.\n",
                                       routine);
                    }

                    zcurve = bsInput->irInput[0].AlibZeroCurve[zcId];
                    break;
                case 'b':    // BS curve family
                    // Check that credit curve is present
                    if (bsInput->nbSPInput == 0)
                    {
                        throw KFailure("%s: no credit name present in the "
                                       "multi-asset wrapper env.\n",
                                       routine);
                    }

                    zcurve = bsInput->spInput[0].AlibZeroCurve;
                    break;
                }


                // Read curve type
                cvtype = getInt(is, "curve type");
                if (cvtype != KV_DET  &&
                    cvtype != KV_DIFF &&
                    cvtype != KV_IDX1 &&
                    cvtype != KV_IDX2 && 
                    cvtype != KV_BASIS)
                        throw KFailure("%s: invalid curve type %d "
                                        "in row %d. "
                                        "%d=Det, %d=Diff, %d=Idx1, %d=Idx2, "
                                        "%d=Basis.\n",
                                        routine, cvtype, i+1,
                                        KV_DET, KV_DIFF, KV_IDX1, KV_IDX2,
                                        KV_BASIS);
                
                // Check if it is a basis curve 
                if (cvtype >= KV_BASIS) 
                        mIsBasis = true;
                        
                // Read curve name in product
                cvname = getString(is, "curve name");
                if (cvname.empty())
                        throw KFailure("%s: can not have empty curve name "
                                       "in row %d.\n",
                                        routine, i+1);        

                mCVTypes.push_back(cvtype);
                mCVNames.insert(KMap(int, String)::value_type(cvtype, cvname));

                mValueDates.insert(KMap(int, TDate)::value_type(
                                   cvtype, zcurve->fBaseDate));
                mCV.insert(KMap(int, KZCurve)::value_type(cvtype,
                                      KZCurve(zcurve, mToday, zcInterp)));

                // Assign names of zc files
                //
                switch (zcType) {
                case 'i': 
                        sprintf(fileName, "%s%d_0.dat", irZCFile, zcId); 
                        mZCFiles.insert(KMap(int, String)::value_type(
                                                cvtype,
                                                String(fileName)));
                        break;

                case 'b':
                        sprintf(fileName, "%s_%d.dat", bsZCFile, zcId); 
                        mZCFiles.insert(KMap(int, String)::value_type(
                                                cvtype,
                                                String(fileName)));
                        break;
                }
        }

        return (mIsBasis);
    }
    catch (KFailure)
    {
        throw KFailure("%s: failed.\n", routine);
    }
}

//---------------------------------------------------------------
//

bool
KMarketCurves::ReadCreditMAW(
        istream& is,                // (I) wrapper data stream
        TDate     today,            // (I) Today's date
        CRX_INPUT *crxInput)        // (I) MAW market data 
{
static        char        routine[] = "KMarketCurves::ReadCreditMAW";
                
                        // --- zero curves
        int     i;
        int     numCurves;
        int     cvtype;
	int	zcInterp = GTO_LINEAR_INTERP;
        String  cvname, zcName;
        char    zcType = 'i';
        int     zcId   = 0;
        char    fileName[256];
        
        TCurve  *zcurve = NULL;

        static  char irZCFile[]  = "ir_curve";
        static  char crZCFile[]  = "cr_zcurve";


    try {
        mIsCredit = false;        // default


        // Today's date
        mToday = today;


        // Zero curves
        //

        // Number of curves
        numCurves = getInt(is, "number of curves");
 
        // Reading in curve infomation:
        //   zcname: curve name in wrapper (ir0, ir1, ir2, cr0, etc.)
        //   cvtype: curve type (diffuse, index, credit, etc, )
        //   cvname: curve name in product 
        for (i=0; i<=numCurves-1; ++i)
        {
                // Read curve name in wrapper env
                zcName  = getString(is, format("curve %d", i));

                GetMAWTZCurveInfo(zcName,
                                  crxInput,
                                  &zcType,
                                  &zcId); 

                switch (zcType) {
                case 'i':    // IR curve family (single currency only)
                    // Check that IR curve is present
                    if (crxInput->nbIRInput == 0)
                    {
                        throw KFailure("%s: no IR curve present in the "
                                       "multi-asset wrapper env.\n",
                                       routine);
                    }

                    zcurve = crxInput->irInput[0].AlibZeroCurve[zcId];
		    zcInterp = GTO_LINEAR_INTERP;
                    break;
                case 'c':    // CR curve family
                    // Check that credit curve is present
                    if (crxInput->nbCRInput == 0)
                    {
                        throw KFailure("%s: no credit name present in the "
                                       "multi-asset wrapper env.\n",
                                       routine);
                    }

                    zcurve = crxInput->crInput[0].AlibZeroCurve;
		    zcInterp = GTO_FLAT_FORWARDS;
                    break;
                }


                // Read curve type
                cvtype = getInt(is, "curve type");
                if (cvtype != KV_DET  &&
                    cvtype != KV_DIFF &&
                    cvtype != KV_IDX1 &&
                    cvtype != KV_IDX2 && 
                    cvtype != KV_CREDIT_RISKY)
                        throw KFailure("%s: invalid curve type %d "
                                        "in row %d. "
                                        "%d=Det, %d=Diff, %d=Idx1, %d=Idx2, "
                                        "%d=Credit.\n",
                                        routine, cvtype, i+1,
                                        KV_DET, KV_DIFF, KV_IDX1, KV_IDX2,
                                        KV_CREDIT_RISKY);
                
                // Check if it is a credit curve 
                if (cvtype >= KV_CREDIT_RISKY) 
                        mIsCredit = true;
                        
                // Read curve name in product
                cvname = getString(is, "curve name");
                if (cvname.empty())
                        throw KFailure("%s: can not have empty curve name "
                                       "in row %d.\n",
                                        routine, i+1);        

                mCVTypes.push_back(cvtype);
                mCVNames.insert(KMap(int, String)::value_type(cvtype, cvname));

                mValueDates.insert(KMap(int, TDate)::value_type(
                                   cvtype, zcurve->fBaseDate));
                mCV.insert(KMap(int, KZCurve)::value_type(cvtype,
                                      KZCurve(zcurve, mToday, zcInterp)));

                // Assign names of zc files
                //
                switch (zcType) {
                case 'i': 
                        sprintf(fileName, "%s%d_0.dat", irZCFile, zcId); 
                        mZCFiles.insert(KMap(int, String)::value_type(
                                                cvtype,
                                                String(fileName)));
                        break;

                case 'c':
                        sprintf(fileName, "%s_%d.dat", crZCFile, zcId); 
                        mZCFiles.insert(KMap(int, String)::value_type(
                                                cvtype,
                                                String(fileName)));
                        break;
                }
        }

        return (mIsCredit);
    }
    catch (KFailure)
    {
        throw KFailure("%s: failed.\n", routine);
    }
}





//----------------------------------------------
//
void
KMarketCurves::ReadMagnet(
        TDate                    today,           // (I) today
        const Array<const TZeroCurve*> &zcCurves, // (I) ALIB zero curve objs
        const Array<int>         &zcTypes,        // (I) zc types (KV_DIFF, ..)
        const Array<String>      &zcNames,        // (I) array of zc names
	const Array<int>	 &zcInterps)	  // (I) interp types
{
static  char routine[] = "KMarketCurves::ReadMagnet";

        int idx;

        TZeroCurve    *zcInput = NULL;
        TZeroCurve    *zcCopy  = NULL;

 try {
        // Check size consistency
        ASSERT_OR_THROW(zcCurves.size() == zcTypes.size());
        ASSERT_OR_THROW(zcCurves.size() == zcNames.size());
        ASSERT_OR_THROW(zcCurves.size() == zcInterps.size());

        //
        // Today's date
        //
        mToday = today;
 
        //
        // Add all the curves
        //
        for (idx = 0; idx < zcCurves.size(); idx++)
        {
                // Make a copy first
                zcInput = const_cast<TZeroCurve*>(zcCurves[idx]);

                ASSERT_OR_THROW (GtoCheckTCurve(
                                    zcInput->curve,
                                    -100.,
                                    const_cast<char*>(zcNames[idx].c_str()),
                                    routine) == SUCCESS);

                zcCopy = GtoZeroCurveMakeCopy(zcInput);
                ASSERT_OR_THROW(zcCopy != NULL);

                // Construct old TCurve structure with static date/rate 
                // from ALib TZeroCurve.  The ALIB TZeroCurve does NOT
                // necessarily pre-populate the static data, and in principle
                // computed on the fly.
                if (GtoEncapsulatedCurveSetRates(zcCopy->curve) 
                    == FAILURE)
                        throw KFailure("%s: failed to convert ALIB zero curve "
                                       "object #%d.\n", routine, idx+1);

                if (zcCopy->curve)
                {

                        InsertZc(KZCurve(zcCopy->curve, today, zcInterps[idx]),
                                 zcTypes[idx],
                                 zcCopy->curve->fBaseDate,
                                 zcNames[idx]);

                        mZCFiles.insert(KMap(int,String)::value_type(
                                        zcTypes[idx], 
                                        "nil"));    // not used

                }

                GtoZeroCurveDelete(zcCopy);
                zcCopy = NULL;
        }        

    }
    catch (KFailure)
    {
        if (zcCopy != NULL) 
        {
            GtoZeroCurveDelete(zcCopy);
            zcCopy = NULL;
        }
        throw KFailure("%s: failed.\n", routine);
    }
}



//----------------------------------------------
//
void
KMarketCurves::ReadMagnetBasis(
        TDate                   today,             // (I) today
        const Array<const TZeroCurve*>&zcCurves,   // (I) ALIB zero curve objs
        const Array<int>        &zcTypes,          // (I) zc types (KV_DIFF, ..)
        const Array<String>     &zcNames,          // (I) array of zc names
	const Array<int>        &zcInterps,	   // (I) interp types
        const Array<String>     &bsInfo)           // (I) basis zero curve info
{

        // Curve table
        //
        ReadMagnet(today,
                   zcCurves, 
                   zcTypes,
                   zcNames,
		   zcInterps);

        //
        // Basis curve info
        //
        mLiborCVName = bsInfo[0];
        mBSDiscCVName = bsInfo[1];
        mBasisDCC     = KDayCc(bsInfo[2].c_str());
        mLiborDCC     = KDayCc(bsInfo[3].c_str());

        switch (toupper((bsInfo[4])[0]))
        {
    case 'P':
            mBSType = PER_SPREAD;
                break;
    case 'S':
            mBSType = SUB_SPREAD;
                break;
    case 'A':
            mBSType = ADD_SPREAD;
                break;
        }

        mBSDelayShift = (bsInfo.size() > 5 ?
                        atof(bsInfo[5].c_str()) : 0e0);


}




//----------------------------------------------
//
static        TCurve*
GetTZCurve(char zcType, TDrWrapperData *drWrapData)
{
static  char    routine[] = "GetTZCurve";
 
        TCurve  *zc = NULL;
 
        zcType = toupper(zcType);
 
        switch (zcType) {
        case 'D':
                zc = drWrapData->fDiscZcCurve;
                break;
        case 'Z':
                zc = drWrapData->fZcCurve;
                break;
        case 'R':
                zc = drWrapData->fRiskZcCurve;
                break;
        case 'B':
                zc = drWrapData->fBasisZcCurve;
                break;
        default:
                throw KFailure("%s: invalid zero curve type (%c).\n", 
                                routine, zcType);
        }
 
        return(zc);
}



//----------------------------------------------
//
static  void
GetMAWTZCurveInfo(string     zcName,     /* (I) curve name */
                  CRX_INPUT  *crxInput,  /* (I) wrapper data input */
                  char       *zcType,    /* (O) curve type, "ir" "cr", .. */
                  int        *zcId)      /* (O) curve ID   */
{
static  char    routine[] = "GetMAWTZCurve";
 
        TCurve  *zc = NULL;

        String  ir_numerics("012");
        String  cr_numerics("0");        // single name credit

 
        String::size_type pos;


 try {

        *zcType = tolower(*(zcName.begin()));
 
        switch (*zcType) {
        case 'i':    // IR curve family (single currency only for time being)
                // curve id
                pos=zcName.find_first_of(ir_numerics);
                if (pos == string::npos)
                    throw KFailure("%s: invalid IR zero curve id (%s).\n", 
                                    routine, zcName.c_str());

                sscanf(&zcName[pos], "%1d", zcId);

                break;

        case 'c':    // CR curve family
                // curve id
                pos=zcName.find_first_of(ir_numerics);
                sscanf(&zcName[pos], "%1d", zcId);

                if (pos == string::npos)
                    throw KFailure("%s: invalid CR zero curve id (%s).\n", 
                                    routine, zcName.c_str());

                break;
        default:
                throw KFailure("%s: invalid zero curve name (%s).\n", 
                                routine, zcName.c_str());
        }
 
    }
    catch (KFailure)
    {
        throw KFailure("%s: failed.\n", routine);
                
    }
}

//----------------------------------------------
//
static  void
GetMAWTZCurveInfo(string     zcName,     /* (I) curve name */
                  BS_INPUT   *bsInput,   /* (I) wrapper data input */
                  char       *zcType,    /* (O) curve type, "ir" "cr", .. */
                  int        *zcId)      /* (O) curve ID   */
{
static  char    routine[] = "GetMAWTZCurve";
 
        TCurve  *zc = NULL;

        String  ir_numerics("012");
        String  bs_numerics("0");        // single basis

 
        String::size_type pos;


 try {

        *zcType = tolower(*(zcName.begin()));
 
        switch (*zcType) {
        case 'i':    // IR curve family (single currency only for time being)
                // curve id
                pos=zcName.find_first_of(ir_numerics);
                if (pos == string::npos)
                    throw KFailure("%s: invalid IR zero curve id (%s).\n", 
                                    routine, zcName.c_str());

                sscanf(&zcName[pos], "%1d", zcId);

                break;

        case 'b':    // BS curve family
                // curve id
                pos=zcName.find_first_of(bs_numerics);
                sscanf(&zcName[pos], "%1d", zcId);

                if (pos == string::npos)
                    throw KFailure("%s: invalid BS zero curve id (%s).\n", 
                                    routine, zcName.c_str());

                break;
        default:
                throw KFailure("%s: invalid zero curve name (%s).\n", 
                                routine, zcName.c_str());
        }
 
    }
    catch (KFailure)
    {
        throw KFailure("%s: failed.\n", routine);
                
    }
}



//--------------------------------------------------------------
// Reads and loads the model/market data.
//

bool
BasisTWrapperRead(
        istream& is,                        // (I) deal input stream
        char *pathDir,                        // (I) directory to read from
        char drwType,                        // (I) '2', 'B'asis

        KMarketCurves&        mktCurves,        // (O) curves and curve types
        KVolDiag&        irVolDiag,        // (O) IR volatility data.
        KMrParam&        irMrParam,        // (O) IR mr data.
        KSmileParam&        irSmileParam,        // (O) IR skew data.
        KVolDiag&        bsVolDiag,        // (O) Basis volatility data.
        KMrParam&        bsMrParam,        // (O) Basis mr data.
        KSmileParam&        bsSmileParam,        // (O) Basis skew data.
        double&                irBsCorr)        // (O) IR basis correlation.

{
static        char        routine[] = "BasisTWrapperRead";
                
        TDrWrapperData        *drWrapData = NULL;
        long                drw_type = DRI_DRW_TYPE2_2CURVES;

                        // --- zero curves

                        // --- vol indices
        char                idx1s[128], idx2s[128];
        KVolCalibIndex      volIdx1, volIdx2;


        KMrParam        treeMrParam;

        int                bsDim;
        int                bsType;
        char                bsLiborCVName[256];
        char                bsDiscCVName[256];
        double                bsDelayShift, bsBackboneQ,
                        bsQ1, bsQ2, bsQF;




                        // --- volatilities for calibration
        TDate        *volDates = NULL;
        double        *volMat = NULL;
        int        *volFreq = NULL;
        double        *volRates = NULL;

                        // --- basis volatility curve
        bool        isBasis = false;        // default
        TCurve  *bsVolCurve = NULL;
        double  bsVolInterp;
        KVector(double) bsVolTmp;

    try {

        //----------------------------------------------
        // (1) Read Dr Wrapper data
        //----------------------------------------------

        IF_FAILED_THROW( DriTDrWrapperDataGetFull(
                pathDir,
                DRI_DRW_BASIS,
                &drWrapData));



        //----------------------------------------------
        // Load zero curves
        //----------------------------------------------

        isBasis = mktCurves.ReadDrw(is, drWrapData);
        

        //----------------------------------------------
        // Calibration indices
        //----------------------------------------------
        strcpy(idx1s, getString(is, "calibration index 1"));
        strcpy(idx2s, getString(is, "calibration index 2"));

        DppReadFromString(idx1s, volIdx1);
        DppReadFromString(idx1s, volIdx2);

        if (!(volIdx1 == volIdx2)) {
                throw KFailure("%s: different volatility calibration"
                        " indices (%s and %s).\n",
                        routine,
                        DppWriteToString(volIdx1).c_str(),
                        DppWriteToString(volIdx2).c_str());
        }

        //----------------------------------------------
        // Distribution type
        //----------------------------------------------
        irSmileParam.ReadDrw(is, drWrapData);



        //----------------------------------------------
        // Read full tree mr parameters
        //----------------------------------------------
        treeMrParam.ReadDrw(is, drWrapData);

        //----------------------------------------------
        // Read extra basis parameters and fill structures.
        //----------------------------------------------

        // Read wrapper data 
        if (isBasis) {

                bsDim  = getInt(is, "basis dimension");
                bsType = getInt(is, "basis type (0=Spd, 1=Percent)");
                strcpy(bsLiborCVName, getString(is, "basis libor curve name"));
                strcpy(bsDiscCVName, getString(is, "basis disc curve name"));

                //bsDelayShift = getDouble(is, "basis delay shift");
                //
                // To allow SIMPLE/PAR stub in basis leg,
                // the delayShift can NOT be anything other than 0.
                //
                bsDelayShift = 0.0;
                if (!IS_ALMOST_ZERO(bsDelayShift))
                        throw KFailure("%s: basis deley shift (%f) != 0!\n",
                                        routine,
                                        bsDelayShift);

                // Read in the bais backbone
                bsBackboneQ = getDouble(is, "basis backbone q");

                bsQ1 = getDouble(is, "basis smile Q1");
                bsQ2 = getDouble(is, "basis smile Q2");
                bsQF = getDouble(is, "basis smile QF");

                //
                // Fill market curve basis info
                //
                switch (bsType) {
                case 0:
                        mktCurves.mBSType = SUB_SPREAD;
                        break;
                case 1:
                        mktCurves.mBSType = PER_SPREAD;
                        break;
                case 2:
                        mktCurves.mBSType = ADD_SPREAD;
                        break;
                default:
                        throw KFailure("%s: bad basis type %d (only 0, 1, or 2).\n",
                                routine, bsType);
                }
                mktCurves.mLiborCVName = String(bsLiborCVName);
                mktCurves.mBSDiscCVName = String(bsDiscCVName);

                if (bsDelayShift >= 0e0 && bsDelayShift <= 1e0)
                        mktCurves.mBSDelayShift = bsDelayShift;
                else
                        throw KFailure("%s: invalid basis delay shift (%f). "
                                       "Must be within the range [0, 1].\n",
                                        routine, bsDelayShift);

                //
                // Fill basis mr parameters
                //
                bsMrParam.mNumFact = bsDim;
                switch (bsDim) {
                case 0:
                        // No basis
                        irMrParam = treeMrParam;

                        bsMrParam    = 0;
                        bsMrParam.mNumFact = bsDim;

                        bsSmileParam = 0;

                        break;
                case 1:
                        switch (treeMrParam.mNumFact) {
                        case 2:
                                // Keep only first dim in IR
                                irMrParam = treeMrParam;
                                irMrParam.mNumFact = 1;

                                // Get second dim for basis
                                bsMrParam = treeMrParam;
                                bsMrParam.mNumFact = 1;
                                bsMrParam.mBeta[0]  = irMrParam.mBeta[1];
                                bsMrParam.mAlpha[0] = irMrParam.mAlpha[1];
                                bsMrParam.mBackboneQ = bsBackboneQ;
        
                                // IR basis correlation.
                                irBsCorr = treeMrParam.mRho[0];

                                break;
                        default:
                                throw KFailure("%s: tree dim %d and bs dim %d "
                                        "(only support 1+1 mode).\n", 
                                        routine, treeMrParam.mNumFact,
                                        bsDim);
                        }
                        break;
                default:
                        throw KFailure("%s: only supports basis dimension of 0 or 1.\n",
                                routine);
                }

                //
                // Basis smile parameters
                //
                bsSmileParam.mQ1 = 1e0 - bsQ1;
                bsSmileParam.mQ2 = 1e0 - bsQ2;
                bsSmileParam.mQF = bsQF;
                bsSmileParam.mNumIter = 0;

        } else {
                //
                // No basis present (dummy params)
                //
                irMrParam = treeMrParam;

                bsMrParam    = 0;
                bsMrParam.mNumFact = 0;
                bsSmileParam = 0;
                irBsCorr = 0;
        }



        if (!volIdx1.IsNil()) {

            //----------------------------------------------
            // IR Volatilities: interpolate at calib index
            //----------------------------------------------

            irVolDiag = KVolDiag(volIdx1, drWrapData);

        } else {
            //----------------------------------------------
            // "NIL" calibration: we set up the dates, but basence
            // of vol is represented by -1 (UGLY !!! MUST BE CHANGED)
            // $$$
            //----------------------------------------------
            KVolCalibIndex dummyVolIdx("3m");
            irVolDiag = KVolDiag(dummyVolIdx, drWrapData);
            irVolDiag = -1e0;
        }



        //----------------------------------------------
        // Basis volatilities: use same timeline as IR
        //----------------------------------------------
        bsVolDiag = irVolDiag;
        bsVolDiag = 0e0;

        // Basis volatility curve
        //
        if (isBasis) {
                bsVolCurve = drWrapData->fBSVolCurve;
                if (bsVolCurve == NULL) {
                    // If no basis curve (must be a type II)
                    // set vols to 1 to use the weight.
                    if (drwType == 'B') {
                        throw KFailure("%s: basis volatility curve "
                               "(basisvol.dat) is NOT available "
                               "from the environment.\n",
                                routine);
                    }
                    bsVolDiag = 1e0;

                } else {
                    // Basis vol curve available, interpolate
                    for (int i=0; i<bsVolDiag.mVolRates.size(); i++)
                    {
                        // interpolate basis vol on these benchmark dates
                        IF_FAILED_THROW(GtoInterpRate(
                                        irVolDiag.mVolDates[i],
                                        bsVolCurve,
                                        GTO_LINEAR_INTERP,
                                        &bsVolInterp));

                        bsVolDiag.mVolRates[i] = bsVolInterp;
                        bsVolDiag.mVolFreqs[i] = 0;
                        
                        // Constant maturity specified by vol curve basis
                        //
                        bsVolDiag.mVolMats[i]  = 1e0/(double)(bsVolCurve->fBasis);
                    }
                }
        } else {
                // No basis
                bsVolDiag = 0e0;
        }



        //----------------------------------------------
        // Logging
        //----------------------------------------------

        dppLog << "===========================================================================" << endl;
        dppLog << routine << ": INPUT MARKET AND MODEL DATA " << endl;
        dppLog << "===========================================================================" << endl;

        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: mktCurves:\n" << mktCurves << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: irVolDiag:\n" << irVolDiag << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: irMrParam:\n" << irMrParam << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: irSmileParam:\n" << irSmileParam << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: bsVolDiag:\n" << bsVolDiag << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: bsMrParam:\n" << bsMrParam << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: bsSmileParam:\n" << bsSmileParam << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: irBsCorr:\n" << irBsCorr << endl;
        dppLog << "===========================================================================" << endl;

        // Free memory
        DriTDrWrapperDataFree(drWrapData);

        return(isBasis);

    }
    catch (KFailure)
    {
        // Free memory
        DriTDrWrapperDataFree(drWrapData);

        throw KFailure("%s: failed (dir=`%s').\n",
                routine, (pathDir ? pathDir : "."));
                
    }
}


bool
BasisTWrapperRead_Mod(
        istream& is,                       // (I) deal input stream
        char *pathDir,                     // (I) directory to read from
        char drwType,                      // (I) '2', 'B'asis

        KMarketCurves&   mktCurves,        // (O) curves and curve types
        KVolDiag&        irVolDiag,        // (O) IR volatility data.
        KMrParam&        irMrParam,        // (O) IR mr data.
        KSmileParam&     irSmileParam,     // (O) IR skew data.
        KVolDiag&        bsVolDiag,        // (O) Basis volatility data.
        KMrParam&        bsMrParam,        // (O) Basis mr data.
        KSmileParam&     bsSmileParam,     // (O) Basis skew data.
        double&          irBsCorr)         // (O) IR basis correlation.
{
static        char        routine[] = "BasisTWrapperRead_Mod";

    bool    isMAWrapper = false;

    try
    {
        ifstream    MASummary ("summary.dat", ios::in);

        if ( MASummary)
            isMAWrapper = true;

        if (isMAWrapper)
        {
            return BasisMAWrapperRead_Mod(
                    is,
                    pathDir,
                    drwType,
                    mktCurves,
                    irVolDiag,
                    irMrParam,
                    irSmileParam,
                    bsVolDiag,
                    bsMrParam,
                    bsSmileParam,
                    irBsCorr);
        }
        else
        {
            return BasisDRW2WrapperRead_Mod(
                    is,
                    pathDir,
                    drwType,
                    mktCurves,
                    irVolDiag,
                    irMrParam,
                    irSmileParam,
                    bsVolDiag,
                    bsMrParam,
                    bsSmileParam,
                    irBsCorr);
        }

    }
    catch (KFailure)
    {
    	throw KFailure ("%s failed\n", routine);
    }
}


//--------------------------------------------------------------
// Reads and loads the model/market data.
// The basis spread vol is given inside the wrapper file
//

bool
BasisDRW2WrapperRead_Mod(
        istream& is,                       // (I) deal input stream
        char *pathDir,                     // (I) directory to read from
        char drwType,                      // (I) '2', 'B'asis

        KMarketCurves&   mktCurves,        // (O) curves and curve types
        KVolDiag&        irVolDiag,        // (O) IR volatility data.
        KMrParam&        irMrParam,        // (O) IR mr data.
        KSmileParam&     irSmileParam,     // (O) IR skew data.
        KVolDiag&        bsVolDiag,        // (O) Basis volatility data.
        KMrParam&        bsMrParam,        // (O) Basis mr data.
        KSmileParam&     bsSmileParam,     // (O) Basis skew data.
        double&          irBsCorr)         // (O) IR basis correlation.

{
static        char        routine[] = "BasisDRW2rapperRead_Mod";
                
        TDrWrapperData        *drWrapData = NULL;
        long                drw_type = DRI_DRW_TYPE2_2CURVES;

                        // --- zero curves

                        // --- vol indices
        char                idx1s[128], idx2s[128];
        KVolCalibIndex        volIdx1, volIdx2;


        KMrParam        treeMrParam;

        int                bsDim;
        String                bsType, bsVolType;
        char                bsLiborCVName[256];
        char                bsLiborDCC[256];
        char                bsBasisDCC[256];
        char                bsDiscCVName[256];
        double                bsDelayShift,
                        bsQ1, bsQ2, bsQF, bsBackboneQ;




                        // --- volatilities for calibration
        TDate        *volDates = NULL;
        double        *volMat = NULL;
        int        *volFreq = NULL;
        double        *volRates = NULL;

                        // --- basis volatility curve
        bool        isBasis = false;        // default
        double  bsVolInterp;
        KVector(double) bsVolTmp;

        

    try {

        //----------------------------------------------
        // (1) Read Dr Wrapper data
        //----------------------------------------------

        IF_FAILED_THROW( DriTDrWrapperDataGetFull(
                    pathDir,
                    DRI_DRW_BASIS,
                    &drWrapData));



        //----------------------------------------------
        // Load zero curves
        //----------------------------------------------

        isBasis = mktCurves.ReadDrw(is, drWrapData);
        

        //----------------------------------------------
        // Calibration indices
        //----------------------------------------------
        strcpy(idx1s, getString(is, "calibration index 1"));
        strcpy(idx2s, getString(is, "calibration index 2"));

        DppReadFromString(idx1s, volIdx1);
        DppReadFromString(idx1s, volIdx2);

        if (!(volIdx1 == volIdx2)) {
                throw KFailure("%s: different volatility calibration"
                        " indices (%s and %s).\n",
                        routine,
                        DppWriteToString(volIdx1).c_str(),
                        DppWriteToString(volIdx2).c_str());
        }


        if (!volIdx1.IsNil()) {

            //----------------------------------------------
            // IR Volatilities: interpolate at calib index
            //----------------------------------------------

            irVolDiag = KVolDiag(volIdx1, drWrapData);

        } else {
            //----------------------------------------------
            // "NIL" calibration: we set up the dates, but basence
            // of vol is represented by -1 (UGLY !!! MUST BE CHANGED)
            // $$$
            //----------------------------------------------
            KVolCalibIndex dummyVolIdx("3m");
            irVolDiag = KVolDiag(dummyVolIdx, drWrapData);
            irVolDiag = -1e0;
        }


        //----------------------------------------------
        // Distribution type
        //----------------------------------------------
        irSmileParam.ReadDrw(is, drWrapData);



        //----------------------------------------------
        // Read full tree mr parameters
        //----------------------------------------------
        treeMrParam.ReadDrw(is, drWrapData);

        // Check if have NORM/LOGN flag.  If not, volType defaults to lognormal.

        String irVolType;
        int c = is.peek();
        if (!is.good()) {
            irVolDiag.mVolType = LOGVOL;
        } else {
            irVolType = getString(is, "ir vol type ('L'ognormal, 'N'ormal)");
            if (!::isdigit(irVolType[0])) {
                if (toupper(irVolType[0]) == 'N')
                    irVolDiag.mVolType = NORMVOL;
                else
                    irVolDiag.mVolType = LOGVOL;
            }
       }

        //----------------------------------------------
        // Read extra basis parameters and fill structures.
        //----------------------------------------------

        // Read wrapper data 
        if (isBasis) {

                if (::isdigit(irVolType[0])) // looks like we read basis dim instead of vol type
                    bsDim  = ::atoi(&irVolType[0]);
                else
                    bsDim  = getInt(is, "basis dimension");

                bsType = getString(is, "basis type ('S'ubtractive spd, 'A'dditive sprd, or 'P'ercent spd)");
                strcpy(bsLiborCVName, getString(is, "basis libor curve name"));
                strcpy(bsDiscCVName, getString(is, "basis disc curve name"));

                strcpy(bsBasisDCC,    getString(is, "basis dcc"));
                strcpy(bsLiborDCC,    getString(is, "libor dcc"));

                //bsDelayShift = getDouble(is, "basis delay shift");
                bsDelayShift = 0e0;

                //
                // To allow SIMPLE/PAR stub in basis leg,
                // the delayShift can NOT be anything other than 0.
                //
                if (!IS_ALMOST_ZERO(bsDelayShift))
                        throw KFailure("%s: basis deley shift (%f) != 0!\n",
                                        routine,
                                        bsDelayShift);

                // Read in the bais backbone
                bsBackboneQ = getDouble(is, "basis backbone q");
                bsBackboneQ = 1.0 - bsBackboneQ;

                bsQ1 = getDouble(is, "basis smile Q1");
                bsQ2 = getDouble(is, "basis smile Q2");
                bsQF = getDouble(is, "basis smile QF");


                //
                // Fill market curve basis info
                //
                switch (toupper(bsType[0])) {
                case 'S':
                        mktCurves.mBSType = SUB_SPREAD;
                        break;
                case 'A':
                        mktCurves.mBSType = ADD_SPREAD;
                        break;
                case 'P':
                        mktCurves.mBSType = PER_SPREAD;
                        break;
                default:
                        throw KFailure("%s: bad basis type %s (only S or P).\n",
                                routine, bsType.c_str());
                }

                mktCurves.mLiborCVName = String(bsLiborCVName);
                mktCurves.mBSDiscCVName = String(bsDiscCVName);

                mktCurves.mBasisDCC = KDayCc(bsBasisDCC);
                mktCurves.mLiborDCC = KDayCc(bsLiborDCC);

                mktCurves.mBasisDCC = KDayCc(bsBasisDCC);
                mktCurves.mLiborDCC = KDayCc(bsLiborDCC);

                if (bsDelayShift >= 0e0 && bsDelayShift <= 1e0)
                        mktCurves.mBSDelayShift = bsDelayShift;
                else
                        throw KFailure("%s: invalid basis delay shift (%f). "
                                       "Must be within the range [0, 1].\n",
                                        routine, bsDelayShift);

                //
                // Fill basis mr parameters
                //
                bsMrParam.mNumFact = bsDim;
                switch (bsDim) {
                case 0:
                        // No basis
                        irMrParam = treeMrParam;

                        bsMrParam    = 0;
                        bsMrParam.mNumFact = 0;

                        bsSmileParam = 0;

                        break;
                case 1:
                        switch (treeMrParam.mNumFact) {
                        case 2:
                                // Keep only first dim in IR
                                irMrParam = treeMrParam;
                                irMrParam.mNumFact = 1;

                                // Get second dim for basis
                                bsMrParam = treeMrParam;
                                bsMrParam.mNumFact = 1;
                                bsMrParam.mBeta[0]  = irMrParam.mBeta[1];
                                bsMrParam.mAlpha[0] = irMrParam.mAlpha[1];
                                bsMrParam.mBackboneQ = bsBackboneQ;        

                                // IR basis correlation.
                                irBsCorr = treeMrParam.mRho[0];

                                break;
                        default:
                                throw KFailure("%s: tree dim %d and bs dim %d "
                                        "(only support 1+1 mode).\n", 
                                        routine, treeMrParam.mNumFact,
                                        bsDim);
                        }
                        break;
                default:
                        throw KFailure("%s: only supports basis dimension of 0 or 1.\n",
                                routine);
                }

                //
                // Basis smile parameters
                //
                bsSmileParam.mQ1 = 1e0 - bsQ1;
                bsSmileParam.mQ2 = 1e0 - bsQ2;
                bsSmileParam.mQF = bsQF;
                bsSmileParam.mNumIter = 0;

        } else {
                //
                // No basis present (dummy params)
                //
                irMrParam = treeMrParam;

                bsMrParam    = 0;
                bsMrParam.mNumFact = 0;
                bsSmileParam = 0;
                irBsCorr = 0;
        }




        //----------------------------------------------
        // Basis volatilities: use same timeline as IR
        //----------------------------------------------
        bsVolDiag = irVolDiag;
        bsVolDiag = 0e0;

        // Basis volatility curve
        //
        if (isBasis) {
            int    i;
            bsVolType = getString(is, "basis vol type ('L'ognormal, 'N'ormal)");
            int    numBSVols   = getInt(is, "number of spread vols");
            TDate  *bsVolDates = new TDate[numBSVols];
            double *bsVolRates = new double[numBSVols];

            if (toupper(bsVolType[0]) == 'L')
                    bsVolDiag.mVolType = LOGVOL;
            else if (toupper(bsVolType[0]) == 'N')
                bsVolDiag.mVolType = NORMVOL;
            else
                throw KFailure("%s: invalid basis vol type (%s).\n",
                                routine, bsVolType.c_str());

            for (i=0; i<=numBSVols-1; i++)
            {
                bsVolDates[i] = getTDate(is, "Spread Vol Date");        
                bsVolRates[i] = getDouble(is, "Spread Vol Rate")/100e0;
            }

            // Basis vol curve available, interpolate
            for (i=0; i<bsVolDiag.mVolRates.size(); i++)
            {
                // interpolate basis vol on these benchmark dates
                // interpolate basis vol on IR vol benchmark dates
                IF_FAILED_THROW(DrlTDateLinearInterp1d(
                             bsVolDates,
                             bsVolRates,
                             numBSVols,
                             irVolDiag.mVolDates[i],
                             &bsVolInterp));

                bsVolDiag.mVolRates[i] = bsVolInterp;
                bsVolDiag.mVolMats[i]  = 1e0/(double)(4);  /* 3M */
                bsVolDiag.mVolFreqs[i] = 0;
            }

            delete [] bsVolDates;
            delete [] bsVolRates;

        } else {
                // No basis
                bsVolDiag = 0e0;
        }



        //----------------------------------------------
        // Logging
        //----------------------------------------------

        dppLog << "===========================================================================" << endl;
        dppLog << routine << ": INPUT MARKET AND MODEL DATA " << endl;
        dppLog << "===========================================================================" << endl;

        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: mktCurves:\n" << mktCurves << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: irVolDiag:\n" << irVolDiag << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: irMrParam:\n" << irMrParam << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: irSmileParam:\n" << irSmileParam << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: bsVolDiag:\n" << bsVolDiag << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: bsMrParam:\n" << bsMrParam << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: bsSmileParam:\n" << bsSmileParam << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: irBsCorr:\n" << irBsCorr << endl;
        dppLog << "===========================================================================" << endl;

        // Free memory
        DriTDrWrapperDataFree(drWrapData);

        return(isBasis);

    }
    catch (KFailure)
    {
        // Free memory
        DriTDrWrapperDataFree(drWrapData);

        throw KFailure("%s: failed (dir=`%s').\n",
                routine, (pathDir ? pathDir : "."));
                
    }
}


//--------------------------------------------------------------
// Reads and loads the model/market data.
// The basis spread vol is given inside the wrapper file
//

bool
BasisMAWrapperRead_Mod(
        istream& is,                       // (I) deal input stream
        char *pathDir,                     // (I) directory to read from
        char drwType,                      // (I) '2', 'B'asis

        KMarketCurves&   mktCurves,        // (O) curves and curve types
        KVolDiag&        irVolDiag,        // (O) IR volatility data.
        KMrParam&        irMrParam,        // (O) IR mr data.
        KSmileParam&     irSmileParam,     // (O) IR skew data.
        KVolDiag&        bsVolDiag,        // (O) Basis volatility data.
        KMrParam&        bsMrParam,        // (O) Basis mr data.
        KSmileParam&     bsSmileParam,     // (O) Basis skew data.
        double&          irBsCorr)         // (O) IR basis correlation.

{
static        char        routine[] = "BasisMArapperRead_Mod";
                
        BS_INPUT           *bsInput = NULL;
        long                DrToday;
        TDate               Today;

                        // --- volatilities for calibration
        TDate               *volDates = NULL;
        double              *volMat = NULL;
        int                 *volFreq = NULL;
        double              *volRates = NULL;

                        // --- credit volatility curve
        bool                isBasis = false;        // default
    try {

        //----------------------------------------------
        // (1) Read Dr Wrapper data
        //----------------------------------------------

        bsInput = (BS_INPUT*)calloc(1, sizeof(BS_INPUT));
        ASSERT_OR_THROW (bsInput != NULL);

        IF_FAILED_THROW( BSMAWReadInput(
                            &DrToday,
                            bsInput,
                            pathDir));

        // Convert to TDate from London long date
        Today = CrxDrDate2TDate(DrToday);


        //----------------------------------------------
        // Load zero curves
        //----------------------------------------------

        isBasis = mktCurves.ReadDrwMAW(is, Today, bsInput);
        
        //----------------------------------------------
        // IR calibration vol curve
        //----------------------------------------------
        irVolDiag.IRVolDiag(bsInput); 

         //----------------------------------------------
        // IR model parameters + tree params
        //----------------------------------------------
        irMrParam.IRMAWMrParam(is, bsInput);

        //----------------------------------------------
        // IR Smile param + CET iterations
        //----------------------------------------------
        irSmileParam.IRMAWSmileParam(is, bsInput);

        //----------------------------------------------
        // Read extra basis parameters and fill structures.
        //----------------------------------------------
        // Read wrapper data 
        if (isBasis) 
        {
            // Check the MAW is configured to contain 1-IR and 1-SP
            if (bsInput->nbSPInput == 0)
            {
                throw KFailure("%s: no spread exists in the "
                               "multi-asset wrapper env!\n",
                               routine);
            }

            // Check IR is 1 factor
            if (bsInput->irInput->NbFactor > 1)
            {
                throw KFailure("%s: IR should be 1 factor only\n",
                               routine);
            }

            //----------------------------------------------
            // SP calibration vol curve
            //----------------------------------------------
            bsVolDiag.SPVolDiag(bsInput); 
            
            //----------------------------------------------
            // SP model parameters 
            //----------------------------------------------
            bsMrParam.SPMAWMrParam(is, bsInput);


            //----------------------------------------------
            // SP Smile param + CET iterations
            //----------------------------------------------
            bsSmileParam.SPMAWSmileParam(is, bsInput);


            //----------------------------------------------
            // IR & CR correlation
            //----------------------------------------------
            irBsCorr = bsInput->corrInput[0][1];

            // Fill basis info
            switch (toupper(bsInput->spInput->BasisType)){
            case 'S':
                mktCurves.mBSType = SUB_SPREAD;
                break;
            case 'A':
                mktCurves.mBSType = ADD_SPREAD;
                break;
            case 'P':
                mktCurves.mBSType = PER_SPREAD;
                break;
            default:
                throw KFailure("%s: bad basis type %c (only S or P).\n",
                                routine, 
                                bsInput->spInput->BasisType);
            }

            /* reference LIBOR curve*/
            char CurveName[MAXBUFF];
            strcpy (CurveName, "ir_curve");
            sprintf(CurveName, "%s%d_%d.dat",
                    CurveName,
                    bsInput->spInput->LiborCrv,
                    bsInput->spInput->BaseIR_id);

            KMap(int, String)::iterator pZCFile;
            pZCFile = mktCurves.mZCFiles.begin();
            while( pZCFile != mktCurves.mZCFiles.end())
            {
                if (strcmp(CurveName, pZCFile->second.c_str()) == 0)
                    break;
                pZCFile++;
            }
            if (pZCFile == mktCurves.mZCFiles.end())
                throw KFailure("%s failed: %s is an invalid curve file name.\n",
                                routine,
                                CurveName);

            KMap(int, String)::iterator pZCLabel = mktCurves.mCVNames.find(
                                    pZCFile->first);
            mktCurves.mLiborCVName = pZCLabel->second;

            /* discounting curve */
            strcpy (CurveName, "ir_curve");
            sprintf(CurveName, "%s%d_%d.dat",
                    CurveName,
                    bsInput->spInput->DiscCrv,
                    bsInput->spInput->BaseIR_id);

            pZCFile = mktCurves.mZCFiles.begin();
            while( pZCFile != mktCurves.mZCFiles.end())
            {
                if (strcmp(CurveName, pZCFile->second.c_str()) == 0)
                    break;
                pZCFile++;
            }
            if (pZCFile == mktCurves.mZCFiles.end())
                throw KFailure("%s failed: %s is an invalid curve file name.\n",
                                routine,
                                CurveName);

            pZCLabel = mktCurves.mCVNames.find( pZCFile->first);
            mktCurves.mBSDiscCVName = pZCLabel->second;           

            /* Basis DCC */
            char CurveDCC[MAXBUFF];
            switch(bsInput->spInput->BasisDCC)
            {
            case '3':
                strcpy(CurveDCC, "30/360");
                break;
            case '5':
                strcpy(CurveDCC, "ACT/365");
                break;
            case '0':
                strcpy(CurveDCC, "ACT/360");
                break;
            default:
                throw KFailure("%s failed: Invalid DCC\n", routine);
            }

            mktCurves.mBasisDCC = KDayCc(CurveDCC);

            /* Libor DCC */
            switch(bsInput->spInput->LiborDCC)
            {
            case '3':
                strcpy(CurveDCC, "30/360");
                break;
            case '5':
                strcpy(CurveDCC, "ACT/365");
                break;
            case '0':
                strcpy(CurveDCC, "ACT/360");
                break;
            default:
                throw KFailure("%s failed: Invalid DCC\n", routine);
            }

            mktCurves.mLiborDCC = KDayCc(CurveDCC);

            mktCurves.mBSDelayShift = 0e0;

 
        } 
        else
        {
            //
            // No basis present (dummy params)
            //
            bsVolDiag           = irVolDiag;
            bsVolDiag           = 0e0;
            bsMrParam           = 0;
            bsMrParam.mNumFact  = 0;
            bsSmileParam        = 0;
            irBsCorr            = 0;
        }

        //----------------------------------------------
        // Logging
        //----------------------------------------------

        dppLog << "===========================================================================" << endl;
        dppLog << routine << ": INPUT MARKET AND MODEL DATA " << endl;
        dppLog << "===========================================================================" << endl;

        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: mktCurves:\n" << mktCurves << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: irVolDiag:\n" << irVolDiag << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: irMrParam:\n" << irMrParam << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: irSmileParam:\n" << irSmileParam << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: bsVolDiag:\n" << bsVolDiag << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: bsMrParam:\n" << bsMrParam << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: bsSmileParam:\n" << bsSmileParam << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: irBsCorr:\n" << irBsCorr << endl;
        dppLog << "===========================================================================" << endl;

        BSMAWFreeInput(bsInput);
        free(bsInput);
        return(isBasis);

    }
    catch (KFailure)
    {

        BSMAWFreeInput(bsInput);
        free(bsInput);
        throw KFailure("%s: failed (dir=`%s').\n",
                routine, (pathDir ? pathDir : "."));
                
    }
}


//--------------------------------------------------------------
// Reads and loads the model/market data.
// The credit spread vol is given inside the wrapper file
//

bool
CreditMAWrapperRead(
        istream&         is,               // (I) deal input stream
        char             *pathDir,         // (I) directory to read from
        char             drwType,          // (I) '2', 'B'asis

        KMarketCurves&   mktCurves,        // (O) curves and curve types
        KVolDiag&        irVolDiag,        // (O) IR volatility data.
        KMrParam&        irMrParam,        // (O) IR mr data.
        KSmileParam&     irSmileParam,     // (O) IR skew data.
        KVolDiag&        crVolDiag,        // (O) Credit volatility data.
        KMrParam&        crMrParam,        // (O) Credit mr data.
        KSmileParam&     crSmileParam,     // (O) Credit skew data.
        double&          irCrCorr)         // (O) IR credit correlation.

{
static        char        routine[] = "CreditMAWrapperRead";
                
        CRX_INPUT           *crxInput = NULL;

        long                DrToday;
        TDate               Today;;


                        // --- vol indices
        KVolCalibIndex      volIdx1, volIdx2;


        KMrParam            treeMrParam;





                        // --- volatilities for calibration
        TDate               *volDates = NULL;
        double              *volMat = NULL;
        int                 *volFreq = NULL;
        double              *volRates = NULL;

                        // --- credit volatility curve
        bool                isCredit = false;        // default
 //       double              crVolInterp;
        KVector(double)     crVolTmp;

    try {

        //----------------------------------------------
        // (1) Read Dr Wrapper data
        //----------------------------------------------

        crxInput = (CRX_INPUT*)calloc(1, sizeof(CRX_INPUT));
        ASSERT_OR_THROW (crxInput != NULL);

        IF_FAILED_THROW( CrxMAWReadInput(
                            &DrToday,
                            crxInput,
                            pathDir));

        // Convert to TDate from London long date
        Today = CrxDrDate2TDate(DrToday);


        //----------------------------------------------
        // Load zero curves
        //----------------------------------------------

        isCredit = mktCurves.ReadCreditMAW(is, Today, crxInput);
        

        //----------------------------------------------
        // IR calibration vol curve
        //----------------------------------------------
        irVolDiag.IRVolDiag(crxInput); 

        //----------------------------------------------
        // IR model parameters + tree params
        //----------------------------------------------
        irMrParam.IRMAWMrParam(is, crxInput);

        //----------------------------------------------
        // IR Smile param + CET iterations
        //----------------------------------------------
        irSmileParam.IRMAWSmileParam(is, crxInput);


        //----------------------------------------------
        // Read extra credit parameters and fill structures.
        //----------------------------------------------

        // Read wrapper data 
        if (isCredit) {

                // Check the MAW is configured to contain 1-IR and 1-CR
                if (crxInput->nbCRInput == 0)
                {
                        throw KFailure("%s: no credit name exists in the "
                                       "multi-asset wrapper env!\n",
                                        routine);
                }


                //----------------------------------------------
                // CR calibration vol curve
                //----------------------------------------------
                crVolDiag.CRVolDiag(crxInput); 

                //----------------------------------------------
                // CR model parameters 
                //----------------------------------------------
                crMrParam.CRMAWMrParam(is, crxInput);

                //----------------------------------------------
                // CR Smile param + CET iterations
                //----------------------------------------------
                crSmileParam.CRMAWSmileParam(is, crxInput);

                //----------------------------------------------
                // Recovery
                //----------------------------------------------
                mktCurves.mRecovery = crxInput->crInput[0].RecoveryRate;

                if (mktCurves.mRecovery > 1e0 ||
                    mktCurves.mRecovery < 0e0)
                        throw KFailure("%s: recovery rate (%f) outside of "
                                       "[0, 1]!\n",
                                        routine,
                                        mktCurves.mRecovery);

                //----------------------------------------------
                // IR & CR correlation
                //----------------------------------------------
                irCrCorr = crxInput->corrInput[0][1];

                
                mktCurves.mIRDiscCVName = 
                       getString(is, "IR discount curve name");

                //----------------------------------------------
                // Credit CET iterations - not part of CR smile
                // to keep backwards compatibility with deal
                // file
                //----------------------------------------------
                if (checkSilent(is)) {
                    crSmileParam.mNumIter = getInt(is,"CR CET Iterations");
                } else {
                    // if there is nothing there, set to no CET, no VNFM
                    crSmileParam.mNumIter = -1;
                }

        } else {
                //
                // No credit present (dummy params)
                //
                irMrParam = treeMrParam;

                crMrParam    = 0;
                crMrParam.mNumFact = 0;
                crSmileParam = 0;
                crVolDiag = 0e0;
                irCrCorr = 0;
        }




        //----------------------------------------------
        // Logging
        //----------------------------------------------

        dppLog << "===========================================================================" << endl;
        dppLog << routine << ": INPUT MARKET AND MODEL DATA " << endl;
        dppLog << "===========================================================================" << endl;

        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: mktCurves:\n" << mktCurves << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: irVolDiag:\n" << irVolDiag << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: irMrParam:\n" << irMrParam << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: irSmileParam:\n" << irSmileParam << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: crVolDiag:\n" << crVolDiag << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: crMrParam:\n" << crMrParam << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: crSmileParam:\n" << crSmileParam << endl;
        dppLog << "---------------------------------------------------------------------------" << endl;
        dppLog << "DATA: irCrCorr:\n" << irCrCorr << endl;
        dppLog << "===========================================================================" << endl;

        // Free memory
        CrxMAWFreeInput(crxInput);
        free(crxInput);

        return(isCredit);

    }
    catch (KFailure)
    {
        // Free memory
        CrxMAWFreeInput(crxInput);
        free(crxInput);

        throw KFailure("%s: failed (dir=`%s').\n",
                routine, (pathDir ? pathDir : "."));
                
    }
}


