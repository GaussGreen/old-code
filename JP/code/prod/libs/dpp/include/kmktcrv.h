/***************************************************************
 * Module:       PenGuin
 * Submodule:       
 * File:       kmktcrv.h
 * Function:       
 * Author:       Christian Daher, David Liu
 ***************************************************************/
#ifndef       _kmktcrv_H
#define       _kmktcrv_H
#include "kstdinc.h"
#include "ktypes.h"
#include "kzcurve.h"
#include "kmodpar.h"
#include "kvoldat.h"
#include "crxwrpio.h"


extern       "C" {
#include "dritkwrp.h"           // TDrWrapper
};


#define KV_DET                  -1
#define KV_DIFF                  0
#define KV_IDX1                  1
#define KV_IDX2                  2

#define KV_BASIS                 3
#define KV_SPREAD                4
#define KV_PAR_SPREAD            5

#define KV_CREDIT_RISKY          6
#define KV_CREDIT_DEFPROB        7
#define KV_CREDIT_PROT_DEFRECOV  8
#define KV_CREDIT_PROT_BINRECOV  9
 
#define MAX_CURVES               15


//--------------------------------------------------------------
/**
 * Class for market curve environment with possibly one basis curve.
 * It contains * the zero curves with name (tring by which it is indentified
 * and curve type in the model (i.e. diffuse, index, basis, etc.)
 */

class       KMarketCurves {
public:
       /** Default constructor. */
              KMarketCurves();

       /** Destructor. */
              ~KMarketCurves();

       /** Write to a stream. */
friend  ostream& operator<<(ostream& os, const KMarketCurves& object);

       /**
        * Writes in yacction format.
        */
virtual ostream& DrWYacctionWrite( ostream& os, int indent=FALSE);

       /**
        * Writes the basis info in yacction format.
        */
ostream& BasisYacctionWrite( ostream& os, int indent=FALSE);

       /**
        * Writes the credit info in yacction format.
        */
ostream& MAWCreditYacctionWrite( ostream& os, int indent=FALSE);

       /**
        * Adds a nonbasis zero curve to the environment.
        */
       void       Insert(
              const KZCurve &zcCurve,              // (I) zc curve
              int zcType,                     // (I) zc type (KV_DIFF, etc.)
              TDate valueDate,              // (I) zc value date
              const String& zcName);              // (I) zc name

       /** 
        * Insert value date for curve index idx
        */
       void       InsertValueDate(
              int       zcType,                     // (I) zc type (KV_DIFF, etc.)
              TDate       valueDate)              // (I) value date for zcurve
       {mValueDates.insert(KMap(int, TDate)::value_type(
                            zcType, valueDate));}

       /**
        * Adds a basis zero curve to the environment.
        */
       void       InsertBasis(
              const KZCurve &zcCurve,              // (I) zc curve
              int zcType,                     // (I) zc type (KV_DIFF, etc.)
              TDate valueDate,              // (I) zc value date
              const String& zcName,              // (I) zc name
              KSpd bsType,                     // (I) SUB_SPREAD/ADD_SPREAD/PER_SPREAD
              const String& bsLiborCVName,       // (I) Libor curve name
              const String& bsDiscCVName,       // (I) Basis disc curve name
              double bsDelayShift,              // (I) Reset delay shift
              const KDayCc& bsDCC,              // (I) Basis DCC
              const KDayCc& liborDCC);       // (I) Libor DCC


       /**
        * Adds a credit zero curve to the environment.
        */
        void     InsertCredit(
        const KZCurve &zcCurve,             // (I) zc curve
        int   zcType,                       // (I) zc type
        TDate valueDate,                    // (I) zc value date
        const String& zcName,               // (I) zc name
        const String& irDiscountName,       // (I) IR disc curve name
        double    recovery);                // (I) Recovery




       /**
        * Add a zero curve regardless of its type.
        */
       void       InsertZc(
              const KZCurve &zcCurve,              // (I) zc curve
              int zcType,                     // (I) zc type (KV_DIFF, etc.)
              TDate valueDate,              // (I) zc value date
              const String& zcName);              // (I) zc name



       /**
        * Reads from the input deal data file (input stream is)
        * the specifications for the zero curves 
        * loads them from the DR wrapper files
        * (located in the pathDir directory or current if NULL).
        * The recognized format for the input file is of the form:
        * <PRE>
        * #Number of curves ( <= 4 )
        * 1
        * # ZC    Type    Name (Type 0=diffuse; 1=index1; 2=index2; 3=basis)
        * Disc    0       Curve2
        * </PRE>
        * For each zero curve, the following are specified:<BR>
        * 1. Zero curve designation in the DR Wrapper (D=discount,
        *    Z=index, R=risky, B=basis) <BR>
        * 2. The type assigned to the curve in the model
        *    (0=diffuse, 1=index1, 2=index2, 3=basis).<BR>
        * 3. The name assigned to the curve in the model (arbitrary).
        * <BR>
        * Return true if basis curves are required,
        * false if pure IR products.
        */
       bool       ReadDrw(
              istream& is,                     // (I) wrapper data stream
              TDrWrapperData *drWrapData);     // (I) wrapper market data 

       /**
        * Reads from the multi-asset env
        * credit env.
        */
       bool       ReadDrwMAW(
              istream&  is,               // (I) wrapper data stream
              TDate     today,            // (I) Today's date
              BS_INPUT *bsInput);         // (I) MAW market data

       /**
        * Reads from the input deal data file (input stream is)
        * credit env.
        */
       bool       ReadCreditMAW(
              istream& is,                // (I) wrapper data stream
              TDate     today,            // (I) Today's date
              CRX_INPUT *crxInput);       // (I) MAW market data


       /**
        * Reads the curves from Magnet interface.
        */
virtual void       ReadMagnet(
       TDate                    today,           // (I) today
       const Array<const TZeroCurve*> &zcCurves, // (I) ALIB zero curve objs
       const Array<int>         &zcTypes,        // (I) zc types (KV_DIFF, ..)
       const Array<String>      &zcNames,        // (I) array of zc names
       const Array<int>         &zcInterps);     // (I) array of zc interps


       /**
        * Reads the basis curves from Magnet interface.
        */
virtual void       ReadMagnetBasis(
       TDate              today,                  // (I) today
       const Array<const TZeroCurve*> &zcCurves,  // (I) ALIB zero curve objs
       const Array<int>        &zcTypes,          // (I) zc types (KV_DIFF, ..)
       const Array<String>     &zcNames,          // (I) array of zc names
       const Array<int>        &zcInterps,        // (I) array of zc interps
       const Array<String>     &bsInfo);          // (I) basis zero curve info

       /**
        * Returns a reference to the diffused curve,
        * throws an exception if not found.
        */
       KZCurve&       GetDiffuse();

       /**
        * Return curve type given the curve name
        */
       int               GetCurveType(String &cvName);
        
       /**
        * Check validity of class members.
        */
       bool              IsValid();

       /**
        * Returns true if contains a basis curve.
        */
       bool              IsBasis();

       /**
        * Returns true if contains a credit curve.
        */
       bool              IsCredit();


       /**
        * Return today's date
        */
       TDate              Today() 
                     { return mToday;}

       /** Return the value date of given curve
        */
       TDate              GetValueDate(int curveIdx);

       /** Return the zero shift between value date and today
        *  for the discount curve.
        */
       double              ZeroShift(const String &discCurveName);


                /** Credit market */
       /**
        * Returns a reference to the IR diffused curve,
        * throws an exception if not found.
        */
       KZCurve&       GetIRDiffuse();

       /**
        * Returns a reference to the IR discout curve,
        * throws an exception if not found.
        */
       KZCurve&       GetIRDiscount();

public:
                            /** Today's date */
       TDate                  mToday;
                            /** curve types */
       KVector(int)           mCVTypes;
                            /** zero curves */
       KMap(int, KZCurve)  mCV;
                            /** curve names */
       KMap(int, String)   mCVNames;
                            /** store the zero curve file name */
       KMap(int, String)   mZCFiles;
                            /** store the zero curve value date */
       KMap(int, TDate)   mValueDates;

                            /** Set to true if basis parameters are needed*/
       bool              mIsBasis;       
                            /** Basis type. SUB_SPREAD/ADD_SPREAD/PER_SPREAD */
       KSpd              mBSType;
                            /** Libor curve name */
       String              mLiborCVName;
                            /** Basis discount curve name*/
       String              mBSDiscCVName;       
                            /** Reset delay shift */
       KDayCc              mBasisDCC;
                            /** Libor day count convention */
       KDayCc              mLiborDCC;
                            /** Reset delay shift */
       double              mBSDelayShift;       


                /** Credit market */

                            /** Set to true if pricing credit */
       bool                mIsCredit;       
                            /** Recovery rate */
       double              mRecovery;
                            /** IR discount curve name. */
       String              mIRDiscCVName;
};





//--------------------------------------------------------------
// Reads and loads the model/market data.
//
bool
BasisTWrapperRead(
       istream& is,                     // (I) deal input stream
       char *pathDir,                     // (I) directory to read from
       char drwType,                     // (I) '2', 'B'asis

       KMarketCurves&       mktCurves,       // (O) today, curves and curve types
       KVolDiag&       irVolDiag,       // (O) IR volatility data.
       KMrParam&       irMrParam,       // (O) IR mr data.
       KSmileParam&       irSmileParam,       // (O) IR skew data.
       KVolDiag&       bsVolDiag,       // (O) Basis volatility data.
       KMrParam&       bsMrParam,       // (O) Basis mr data.
       KSmileParam&       bsSmileParam,       // (O) Basis skew data.
       double&              irBsCorr);       // (O) IR basis correlation.

//--------------------------------------------------------------
// Reads and loads the model/market data.
// Libor and basis rate DCCs and the basis spread vol are 
// provided inside the wrapper file
//
bool
BasisTWrapperRead_Mod(
       istream& is,                     // (I) deal input stream
       char *pathDir,                     // (I) directory to read from
       char drwType,                     // (I) '2', 'B'asis

       KMarketCurves&       mktCurves,       // (O) today, curves and curve types
       KVolDiag&       irVolDiag,       // (O) IR volatility data.
       KMrParam&       irMrParam,       // (O) IR mr data.
       KSmileParam&       irSmileParam,       // (O) IR skew data.
       KVolDiag&       bsVolDiag,       // (O) Basis volatility data.
       KMrParam&       bsMrParam,       // (O) Basis mr data.
       KSmileParam&       bsSmileParam,       // (O) Basis skew data.
       double&              irBsCorr);       // (O) IR basis correlation.

bool
BasisDRW2WrapperRead_Mod(
       istream& is,                    // (I) deal input stream
       char *pathDir,                  // (I) directory to read from
       char drwType,                   // (I) '2', 'B'asis

       KMarketCurves&  mktCurves,       // (O) today, curves and curve types
       KVolDiag&       irVolDiag,       // (O) IR volatility data.
       KMrParam&       irMrParam,       // (O) IR mr data.
       KSmileParam&    irSmileParam,    // (O) IR skew data.
       KVolDiag&       bsVolDiag,       // (O) Basis volatility data.
       KMrParam&       bsMrParam,       // (O) Basis mr data.
       KSmileParam&    bsSmileParam,    // (O) Basis skew data.
       double&         irBsCorr);       // (O) IR basis correlation.

bool
BasisMAWrapperRead_Mod(
       istream& is,                     // (I) deal input stream
       char *pathDir,                   // (I) directory to read from
       char drwType,                    // (I) '2', 'B'asis

       KMarketCurves&  mktCurves,       // (O) today, curves and curve types
       KVolDiag&       irVolDiag,       // (O) IR volatility data.
       KMrParam&       irMrParam,       // (O) IR mr data.
       KSmileParam&    irSmileParam,    // (O) IR skew data.
       KVolDiag&       bsVolDiag,       // (O) Basis volatility data.
       KMrParam&       bsMrParam,       // (O) Basis mr data.
       KSmileParam&    bsSmileParam,    // (O) Basis skew data.
       double&         irBsCorr);       // (O) IR basis correlation.



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
        double&          irCrCorr);        // (O) IR credit correlation.


#endif


