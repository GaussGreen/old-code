/****************************************************************
 * Module:	VirtualProduct
 * Submodule:	
 * File:	vpdefprotect.h
 * Function:	
 * Author:	David Liu, 6/7/1/2005
 *****************************************************************/
#ifndef	_vpdefprotect_H
#define	_vpdefprotect_H

#include "vpbase.h"

/**     
 * Default protection type     
 */     
enum DefType {
    DEF_EXPOSURE,    // default exposure         
    DEF_KNOCKIN      // default knock-in
};


//--------------------------------------------------------------
/**
 * Class for default defprotect of general underlying definition. 
 */


class KVPDefProtect : public KVPInstr {
public:

	/**
	 * Convenience constructor given observation freqency.
	 */
	KVPDefProtect(
		const char 	*name,		// (I) name
		TDate           startDate,      // (I) start date of period
		TDate           endDate,        // (I) end date of period
		TDateInterval   &freq,	        // (I) Obs frequency as interval
		DefType         type,	        // (I) default protection type
		double          recovery,	// (I) default recovery rate
		const char* discZcName);	// (I) discount curve

	/**
	 * Convenience constructor for arbitrary start and end schedules
	 */
	KVPDefProtect(
		const char            *name,         // (I) name
		const KVector(TDate)  &startDates,   // (I) Start dates
		const KVector(TDate)  &endDates,     // (I) End dates
		const KVector(TDate)  &settleDates,  // (I) Settlement dates
		DefType               type,	     // (I) default protect type
		double                recovery,	     // (I) default recovery 
		const char*           discZcName);   // (I) discount curve

	/**
	 * Convenience constructor with rebates for arbitrary start 
         * and end schedules
	 */
	KVPDefProtect(
		const char            *name,         // (I) name
		const KVector(TDate)  &startDates,   // (I) Start dates
		const KVector(TDate)  &endDates,     // (I) End dates
		const KVector(TDate)  &settleDates,  // (I) Settlement dates
		const KVector(double) &rebates,      // (I) Rebates
		DefType               type,	     // (I) default protect type
		double                recovery,	     // (I) default recovery 
		const char*           discZcName);   // (I) discount curve


	/**
	 * Destructor
	 */
	~KVPDefProtect() {};


	/**
	 * Check validity of structure
	 */
void	CheckIsValid() const;

	/**
	 * Type name.
	 */
virtual	const char*	TypeName() const {return("KVPDefProtect");}

	/**
	 * Default type.
	 */
virtual	const DefType	GetDefType() const {return mDefType;}

	/**
	 * Reads the object from a stream.
	 */
virtual	istream& Get(istream& is, int drw=FALSE);

	/**
	 * Writes the object to a stream.
	 */
virtual	ostream& Put(ostream& os, int indent = FALSE) const;

	/**
	 * Writes in yacction format.
	 */
virtual ostream& YacctionWrite( ostream& os, int indent=FALSE);

	/**
	 * Reads from a stream.
	 */
friend	istream& operator>>(istream& is, KVPDefProtect& asset)
		{ asset.Get(is); return (is);}

	/**
	 * Writes to a stream.
	 */
friend	ostream& operator<<(ostream& os, const KVPDefProtect& asset)
		{ asset.Put(os); return (os);}


	/*
	 * Data
	 */
private:
				/** Default protection start dates */
	KVector(TDate)	mStartDates;
				/** Default protection end dates */
	KVector(TDate)	mEndDates;
				/** Default protection settelement dates */
	KVector(TDate)	mSettleDates;
				/** Default protection rebates */
	KVector(double)	mRebates;
				/** Default protection type */
	DefType	        mDefType;
				/** Default recovery rate.  Only for exposure
                                 *  calculation. 
                                 */
	double	        mRecovery;
				

				/* So can access */
friend	class 	KVPToolDefProtect;

};



#endif	/* _vpdefprotect_H */

