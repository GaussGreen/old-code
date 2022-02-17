/****************************************************************
 * Module:  VirtualProduct
 * Submodule:   
 * File:    vpkio.h
 * Function:    
 * Author:  Christian Daher
 * Revision:    $Header$
 *****************************************************************/
#ifndef _vpkio2Idx_H
#define _vpkio2Idx_H

#include "vpbase.h"
#include "krate.h"
#include "vpkio.h"


//--------------------------------------------------------------
/**
 * Class for knock in/out definition. 
 */


class KVPKnockIO2Idx : public KVPKnockIO {
public:

    /**
     * Convenience constructor given observation freqency.
     * Time-independent dual ko rate index.
     */
    KVPKnockIO2Idx(
        const char  *name,                    // (I) name
        const KKnockIO  &ioType,              // (I) Knock in/out type
        const KKnockIO  &ioWindow1,           // (I) Pay in/out 1
        const SharedPointer<KRate> rateIndex1,// (I) Knock in/out rate index 1
        const KKnockIO  &ioWindow2,           // (I) Pay in/out 2
        const SharedPointer<KRate> rateIndex2,// (I) Knock in/out rate index 2
        const KSmooth   &smooth,              // (I) Node smoothing flag
        TDate           startDate,            // (I) Not included
        TDate           matDate,              // (I) Maturity date
        const TDateInterval &freq,            // (I) Obs frequency as interval
        TBoolean        stubAtEnd,            // (I) stub
        const KDateInterval   &notDays,       // (I) # of notif days
        const KVector(TDate)  &barrierDates,  // (I) barrier dates
        const KVector(double) &barrierLo1,    // (I) Low barrier  1
        const KVector(double) &barrierHi1,    // (I) High barrier 1
        const KVector(double) &barrierLo2,    // (I) Low barrier  2
        const KVector(double) &barrierHi2,    // (I) High barrier 2
        const KVector(double) &rebates,       // (I) Rebates
        const char* discZcName);              // (I) discount curve


    /**
     * Convenience constructor for arbitrary observation and
     * settlement dates.  Information on effective observation dates 
     * implied from the rateIndex spot offset. 
     * Time-independent dual ko rate index.
     */
    KVPKnockIO2Idx(
        const char *name,           // (I) name
        const KKnockIO &ioType,     // (I) Knock in/out type
        const KKnockIO &ioWindow1,  // (I) Pay in/out
        const SharedPointer<KRate> &rateIndex1,// (I) Ko rate index
        const KKnockIO &ioWindow2,  // (I) Pay in/out
        const SharedPointer<KRate> &rateIndex2,// (I) Ko rate index
        const KSmooth &smooth,      // (I) Node smoothing flag

        const KVector(TDate)  &obsDates,      // (I) Observ dates
        const KVector(TDate)  &settleDates,   // (I) Settlement dates
        const KVector(double) &barrierLo1,    // (I) Low barrier  1
        const KVector(double) &barrierHi1,    // (I) High barrier 1
        const KVector(double) &barrierLo2,    // (I) Low barrier  2
        const KVector(double) &barrierHi2,    // (I) High barrier 2
        const KVector(double) &rebates,       // (I) Rebates
        const char*     discZcName);          // (I) discount curve


    /**
     * Convenience constructor for arbitrary observation and
     * settlement dates.  Information on effective observation dates 
     * is given explicitly. Time-independent dual ko rate index.
     */
    KVPKnockIO2Idx(
        const char *name,           // (I) name
        const KKnockIO &ioType,     // (I) Knock in/out type
        const KKnockIO &ioWindow1,  // (I) Pay in/out
        const SharedPointer<KRate> &rateIndex1,// (I) Ko rate index
        const KKnockIO &ioWindow2,  // (I) Pay in/out
        const SharedPointer<KRate> &rateIndex2,// (I) Ko rate index
        const KSmooth &smooth,      // (I) Node smoothing flag

        const KVector(TDate)  &obsDates,      // (I) Observ dates
        const KVector(TDate)  &obsEffDates,   // (I) Observ eff dates
        const KVector(TDate)  &settleDates,   // (I) Settlement dates
        const KVector(double) &barrierLo1,    // (I) Low barrier  1
        const KVector(double) &barrierHi1,    // (I) High barrier 1
        const KVector(double) &barrierLo2,    // (I) Low barrier  2
        const KVector(double) &barrierHi2,    // (I) High barrier 2
        const KVector(double) &rebates,       // (I) Rebates
        const char*     discZcName);          // (I) discount curve 


    /**
     * Convenience constructor for arbitrary observation dates
     * settlement dates.  Effective observation dates are given
     * explicitly. One pair of knock rate indices per observation
     * 
     */
    KVPKnockIO2Idx(
        const char *name,           // (I) name
        const KKnockIO &ioType,     // (I) Knock in/out type
        const KKnockIO &ioWindow1,  // (I) Pay in/out
        const KVector(SharedPointer<KRate>) &rateIndex1,// (I) Ko rate index
        const KKnockIO &ioWindow2,  // (I) Pay in/out
        const KVector(SharedPointer<KRate>) &rateIndex2,// (I) Ko rate index
        const KSmooth &smooth,      // (I) Node smoothing flag

        const KVector(TDate)  &obsDates,      // (I) Observ dates
        const KVector(TDate)  &obsEffDates,   // (I) Observ eff dates
        const KVector(TDate)  &settleDates,   // (I) Settlement dates
        const KVector(double) &barrierLo1,    // (I) Low barrier  1
        const KVector(double) &barrierHi1,    // (I) High barrier 1
        const KVector(double) &barrierLo2,    // (I) Low barrier  2
        const KVector(double) &barrierHi2,    // (I) High barrier 2
        const KVector(double) &rebates,       // (I) Rebates
        const char*     discZcName);          // (I) discount curve


    /**
     * Destructor
     */
    ~KVPKnockIO2Idx();


    /**
     * Include only the events with observation date >= today.
     */
virtual void    ValidEvents(TDate today);   // today's date

    /**
     * Type name.
     */
virtual const char* TypeName() const {return("KVPKnockIO2Idx");}

    /**
     * Reads the object from a stream.
     */
virtual istream& Get(istream& is, int drw=FALSE);

    /**
     * Writes the object to a stream.
     */
virtual ostream& Put(ostream& os, int indent = FALSE) const;

    /**
     * Writes in yacction format.
     */
virtual ostream& YacctionWrite( ostream& os, int indent=FALSE);

    /**
     * Reads from a stream.
     */
friend  istream& operator>>(istream& is, KVPKnockIO2Idx& asset)
        { asset.Get(is); return (is);}

    /**
     * Writes to a stream.
     */
friend  ostream& operator<<(ostream& os, const KVPKnockIO2Idx& asset)
        { asset.Put(os); return (os);}


    /*
     * Data
     */
private:
                    
    KKnockIO             mIOWindow2;  /** 2nd ko index in / out */

    KVector(KRateReset*) mRateIndex2; /** 2nd ko index */

    KVector(TDate)       mObsDates2;  /** duplicate list -- needed in
                                        * ValidEvents function */

    KVector(double)      mBarrierLo2; /** low barriers 2 */
                    
    KVector(double)      mBarrierHi2; /** high barriers 2 */

                    
friend  class   KVPToolKnockIO2Idx;   /* So can access */

};


#endif  /* _vpkio2Idx_H */

