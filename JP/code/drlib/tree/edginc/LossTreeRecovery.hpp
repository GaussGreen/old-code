//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : LossTreeRecovery.hpp
//
//   Description : Recovery model for spread loss tree
//					recovery can be afunction of number of defaults that have incurred so far
//					Note that this has to be independent of the time of default
//
//   Author      : Matthias Arnsdorf
//
//   Date        : June 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_LOSS_TREE_RECOVERY_HPP
#define QLIB_LOSS_TREE_RECOVERY_HPP

#include "edginc/Object.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE


class LossTreeRecovery: 
	public CObject
{

public:
    static CClassConstSP const TYPE;
    
	virtual ~LossTreeRecovery();

	/** Called immediately after object constructed */
    virtual void validatePop2Object();    

	/** set the number of names in the portfolio */
	void setNumNames(int numberOfNames);
	


	/** recovery rate for next default given the number of defaults that have happend so far */
	double recovery(
		int numberOfDefaults // number of defaults incurred up to now
		) const;
	/** returns the total loss incurred up to date given the number of defaults  
	as proportion of unit portfolio notional*/
	double totalLoss(
		int numberOfDefaults // number of defaults incurred up to now
		) const;

	/** returns the total recovery  up to date given the number of defaults 
	 as proportion of unit portfolio notional*/
	double totalRecovery(
		int numberOfDefaults // number of defaults incurred up to now
		) const;
	
	/** array of total losses for each default  as proportion of unit portfolio notional*/
	DoubleArrayConstSP getTotalLosses() const {return totalLosses;} 

	/** array of total recovories for each default  as proportion of unit portfolio notional*/
	DoubleArrayConstSP getTotalRecoveries() const {return totalRecoveries;} 

	/** array of recovery for each default */
	DoubleArrayConstSP getRecoveries() const {return recoveries;} 

private:
    // For reflection
    static void load (CClassSP& clazz);

	/** private constructor */
	LossTreeRecovery();

	static IObject* defaultLossTreeRecovery();

	/** set up recovery and loss arrays */
	void initialise();

	//INPUTS
	/** array of recovories per default. recoveries[i] gives the recovery for the next default given we have
	i defaults now*/
	DoubleArraySP recoveries;


	// TRANSIENT FIELDS

	/** number of names in portfolio */
	int numNames;
	
	
	/** array of losses per default as proportion of unit portfolio notional*/
	DoubleArraySP totalLosses;

	
	/** array of total recoveries per default as proportion of unit portfolio notional */
	DoubleArraySP totalRecoveries;

	// INPUT FIELDS
};

DECLARE(LossTreeRecovery);

DRLIB_END_NAMESPACE

#endif
