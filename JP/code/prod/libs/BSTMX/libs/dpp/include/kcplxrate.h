/****************************************************************
 * Module:  Complex rate
 * Submodule:	
 * File:    kcplxrate.h
 * Function:
 * Author:  Changhong He
 * Revision:$Header$
 *****************************************************************/
#ifndef	_KCPLXRATE_H
#define	_KCPLXRATE_H

#include "krate.h"
#include "kvpatom.h"

class KCplxRate : public KVPAtom{
public:
    
    KCplxRate();
    
    KCplxRate ( const String&          name,         // (I) index name
                const KVector(KRate*)& rateInstr,    // (I) instruments
                const String&          formula);     // (I) formula

    //******************** Constructs from a single rate ***************
    // Compatible to KRate constructors

    KCplxRate ( const KRate& rate);     // (I) CplxRate is a single rate

    KCplxRate& operator = (const KRate& cplxRate);

    KCplxRate(
        const KDateInterval &mat,       // (I) rate maturity
        const KDateInterval &freq,      // (I) rate payment freq
        KDayCc dayCc,                   // (I) rate day count conv
        const KDateInterval &spotOffset,// (I) spot offset
        double spread,                  // (I) spread
        double weight);                 // (I) weight

    /**
     * General constructor with curve name.
     */
    KCplxRate(
        const String  &curveName,       // (I) curve name
        const KDateInterval &mat,       // (I) rate maturity
        const KDateInterval &freq,      // (I) rate payment freq
        KDayCc dayCc,                   // (I) rate day count conv
        const KDateInterval &spotOffset,// (I) spot offset
        double spread,                  // (I) spread
        double weight);                 // (I) weight

    /**
     * General constructor for a simple rate between two
     * specified dates.
     */
    KCplxRate(
        const String  &curveName,       // (I) curve name
        const TDate   startDate,        // (I) start date
        const TDate   endDate,          // (I) end date
        KDayCc dayCc,                   // (I) rate day count conv
        const KDateInterval &spotOffset,// (I) spot offset
        double spread,                  // (I) spread
        double weight);                 // (I) weight


    /**
     * General constructor with debugging and curve name.
     */
    KCplxRate(
        const String  &name,            // (I) debug name
        const String  &curveName,       // (I) curve name
        const KDateInterval &mat,       // (I) rate maturity
        const KDateInterval &freq,      // (I) rate payment freq
        KDayCc dayCc,                   // (I) rate day count conv
        const KDateInterval &spotOffset,// (I) spot offset
        double spread,                  // (I) spread
        double weight);                 // (I) weight

    //**************************************************************

    // Copy constructor
    KCplxRate ( const KCplxRate& cplxRate); 

    // Copy operator    
    KCplxRate& operator = (const KCplxRate& cplxRate);

    // Destructor
    ~KCplxRate();
    
    // Type Name
    virtual const char* TypeName() const{return ("KCplxRate");}
    
    /**
     * Reads the object from a stream.
     */
    virtual istream& Get(istream& is, int drwFmt = FALSE);

    /**
     * Writes the object to a stream.
     */
    virtual ostream& Put(ostream& os, int indent=FALSE) const;

    /**
     * Writes in yacction format with the name.
     */
    virtual ostream& YacctionWrite( ostream& os, int indent=FALSE);

    /**
     * Writes in yacction format without the name.
     */
    virtual ostream& WriteSimple( ostream& os);

    /**
     * Reads from a stream.
     */
    friend  istream& operator>>(istream& is, KCplxRate& cplxRate)
        { cplxRate.Get(is); return (is);}

    // 
    int	IsAllFloating() const;

    /**
     * Get the curve name of the rate
     */
    const String&   CurveName() const
        { return mCurveName;}

    /**
     * Set the curve name.
     */
    void    SetCurveName(const String curveName)
        { mCurveName = curveName;}

    void    SetCurveName(const char *curveName)
        { mCurveName = curveName;}

    // Sets the name of KCplxRate object
    void SetName(const char* name)
        {   mName = name;    }

    void SetName(const String name)
        {   mName = name;    }

    // Gets the name of KCplxRate object
    virtual const char* GetName() const
        { return mName.c_str(); }

    int nbInstr() const
        { return mRates.size(); }

    //KRate* instrIdx (int) const;

    KRate& instrIdx (int) const;

    const String& getFormula() const;

protected:
    KVector(KRate*) mRates;
   
    String mformula;
                /** name for geometry of the tree slice */           
    String mCurveName;


};



class KCplxRateReset : public KVPAtom
{
public:

    KCplxRateReset( 
        const KVector(KRateReset*)& resetIdx,    // (I) reset instruments
        const KCplxRate&            cplxRate);   // (I) complex rate

    KCplxRateReset (
        TDate               resetDate,    //  (I) observation dates
        const KCplxRate&    cplxRate);     // (I) rate description

    KCplxRateReset (
        TDate   resetDate,    // (I) observation dates
        TDate   effDate,      // (I) effective dates
        const KCplxRate&        cplxRate);     // (I) rate description

    KCplxRateReset (
        TDate   resetDate,     // (I) observation dates
        TDate   effDate,       // (I) effective dates
        TDate   endDate,       // (I) end dates
        const KCplxRate &cplxRate        // (I) rate description
        );

    //******************** Constructs object from KRateReset ***************
    // Compatible to KRateReset constructors
    KCplxRateReset( KRateReset& );

    // Constructs object from KRate
    KCplxRateReset(
        TDate resetDate,        // (I) observation date
        const KRate& rate);     // (I) rate description

    KCplxRateReset(
        TDate resetDate,        // (I) observation date
        TDate effDate,          // (I) effective date
        const KRate& rate);     // (I) rate description

    KCplxRateReset(
        const String& curveName,    // (I) curve name
        TDate resetDate,            // (I) observation date
        TDate effDate,              // (I) effective date
        const KRate& rate);         // (I) rate description

    /**
     * Creates a simple stub rate reset between effDate and endDate,
     * using the same day count convention, spread, etc. as "rate".
     */
    KCplxRateReset(
        TDate resetDate,        // (I) observation date
        TDate effDate,          // (I) effective date
        TDate endDate,          // (I) End date
        const KRate& rate);     // (I) rate description

    //**************************************************************

    // Copy constructor
    KCplxRateReset( const KCplxRateReset& );

    // Destructor
    ~KCplxRateReset();
    
    //
    bool operator < (const KCplxRateReset& cplxRateReset) const;

    virtual const char* TypeName() const { return ("KCplxRateReset");}
    
    // returns the instrument
    const KRate& RateIndex(int idx) const
        { return mcplxRate.instrIdx(idx); }

    // return the 1st instrument
    const KRate& Rate() const
        { return RateIndex(0);    }

    // returns RateReset
    const KRateReset& RateResetIndex(int idx) const
        { return *(mRateResets[idx]);    }

    // returns the underlying
    const KCplxRate&  CplxRate() const
        { return mcplxRate;    }

    KCplxRate& CplxRate() 
        { return mcplxRate;    }

    // returns the reset date
    TDate ResetDate() const
        { return mResetDate;   }

    // returns the effective date
    TDate EffDate() const
        { return mEffDate;     }

    // generates a series of zero reset dates from rate reset

//    void GetZeroResetList(KVector(KZeroReset)  &) const;

    /**
     * Writes to a stream.
     */
    virtual ostream& Put(ostream& os, int indent=FALSE) const;

    /**
     * Writes in yacction format.
     */
    virtual ostream& YacctionWrite( ostream& os, int indent=FALSE)
        {Put(os, indent); return os;}

private:

    KCplxRate            mcplxRate;

    KVector(KRateReset*) mRateResets;

    TDate       mResetDate;

    TDate       mEffDate;
};


#endif
