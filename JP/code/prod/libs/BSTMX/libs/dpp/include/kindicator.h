/****************************************************************
 * Module:  Indicator
 * Submodule:	
 * File:    kindicator.h
 * Function:
 * Author:  Changhong He
 * Revision:$Header$
 *****************************************************************/
#ifndef	_KINDICATOR_H
#define	_KINDICATOR_H

#include "kcplxrate.h"
#include "kvpatom.h"
#include "ktypes.h"


class KIndicator : public KVPAtom{
public:
    KIndicator();

    KIndicator( const String&          name,         // (I) indicator name
                const KVector(KRate*)& rateInstr,    // (I) instruments
                const String&          formula,      // (I) formula
                double                 lbarrier,     // (I) lower barrier
                double                 hbarrier,     // (I) higher barrier
                const KKnockIO        &ioWindow,     // (I) Inside or Outside
                const KSmooth         &smooth);      // (I) Node smoothing flag

    KIndicator( 
        const String&          name,         // (I) indicator name
        const KCplxRate       *rateIndex,// (I) complex rate index
        double                 lbarrier,     // (I) lower barrier
        double                 hbarrier,     // (I) higher barrier
        const KKnockIO        &ioWindow,     // (I) Inside or Outside
        const KSmooth         &smooth);      // (I) Node smoothing flag

    KIndicator( 
        const String&          name,         // (I) indicator name
        const KRate           *rateIndex,// (I) complex rate index
        double                 lbarrier,     // (I) lower barrier
        double                 hbarrier,     // (I) higher barrier
        const KKnockIO        &ioWindow,     // (I) Inside or Outside
        const KSmooth         &smooth);      // (I) Node smoothing flag

    // copy constructor

    KIndicator (const KIndicator& Indicator);
   
    // Copy operator    
    KIndicator& operator = (const KIndicator& Indicator);

    ~KIndicator();

     // Type Name
    virtual const char* TypeName() const{return ("KIndicator");}
    
    /**
     * Get the curve name of the rate
     */
    const String&   CurveName() const
        { return mCplxRate->CurveName();}

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
    friend  istream& operator>>(istream& is, KIndicator& Indicator)
        { Indicator.Get(is); return (is);}
   
    // Sets the name of KIndicator object
    void SetName(const char* name)
        {   mName = name;    }

    void SetName(const String name)
        {   mName = name;    }

    // Gets the name of KIndicator object
    virtual const char* GetName() const
        { return mName.c_str(); }

    // returns the underlying complex rate
    const KCplxRate&  CplxRate() const
        { return *mCplxRate;    }

    void setLBarrier(double lb)
        { mLbarrier = lb; }

    void setHBarrier(double hb)
        { mHbarrier = hb; }

    void setIoWindow(KKnockIO ioFlag)
        { mIoWindow = ioFlag; }

    void setSmooth(KSmooth smoothFlag)
        { mSmooth = smoothFlag; }

    double getLBarrier() const
        { return mLbarrier ;}

    double getHBarrier() const
        { return mHbarrier ;}

    KKnockIO getIoWindow() const
        { return mIoWindow ;}
    
    KSmooth getSmooth() const
        { return mSmooth ;}

protected:
    
    KCplxRate   *mCplxRate;
    
    double      mLbarrier;

    double      mHbarrier;

    KKnockIO    mIoWindow;

    KSmooth     mSmooth;

};



#endif