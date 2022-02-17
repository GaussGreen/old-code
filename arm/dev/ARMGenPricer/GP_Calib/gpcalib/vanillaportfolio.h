/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file vanillaportfolio.h
 *
 *  \brief vanilla portfolio
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */


#ifndef _INGPCALIB_VANILLAPORTFOLIO_H
#define _INGPCALIB_VANILLAPORTFOLIO_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/rootobject.h"

/// gpcalib
#include "vanillaarg.h"
#include "typedef.h"


/// gpinfra
#include "gpbase/port.h"
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/gpvector.h"

/// ARM Kernel
#include <glob/expt.h>

/// STL
#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )

#ifdef __GP_STRICT_VALIDATION
	#define CHECKRANGE(i)  CheckRange(i);
#else
	#define CHECKRANGE(i)  
#endif



struct ARM_VanillaPortfolio: public ARM_VanillaArg
{
public:
	ARM_VanillaPortfolio(
		const vector<ARM_VanillaArgPtr>& assets,
        const ARM_GP_Vector& weights,
        const ARM_GP_Vector& mktprices,
        const ARM_GP_Vector& vegas );
	
	ARM_VanillaPortfolio( const ARM_VanillaPortfolio& rhs );
	ARM_VanillaPortfolio& operator=( const ARM_VanillaPortfolio& rhs );
	virtual ~ARM_VanillaPortfolio();
	
	/// accessors
	inline ARM_VanillaArgPtr GetAsset( size_t i ) const { CHECKRANGE(i); return itsElems[i]->GetAsset(); }
	inline double GetMktPrice( size_t i ) const { CHECKRANGE(i); return itsElems[i]->GetMktPrice(); }
	inline double GetWeight( size_t i ) const { CHECKRANGE(i); return itsElems[i]->GetWeight(); }
	inline double GetPrecisiont( size_t i ) const { CHECKRANGE(i); return itsElems[i]->GetPrecisiont(); }
	inline size_t size() const { return itsElems.size(); }
	
	/// for the calibration validation
    bool HasSameAssets() const;
    bool AreAssetsGrowingByExpiry() const;
    ARM_VanillaPortfolio* GetSubPortfolio(size_t begin, size_t end) const;
    vector<ARM_VanillaPortfolio*> SplitPortfolioByMaturity() const;
	
	/// vanilla arg support
	virtual double Price(ARM_PricingModel* model ) const;
    virtual double ImpliedVol(ARM_PricingModel* model) const;
    virtual double GetMktPrice( ) const;
	virtual VanillaArgType GetType() const { return ARM_VanillaArg::VANILLA_PORTFOLIO; }
	
	/// standard ARM support
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual ARM_Object* Clone() const { return new ARM_VanillaPortfolio(*this);}
	
private:
	/// forward declaration
	struct ARM_PortfolioElem
	{
	private:
		ARM_VanillaArgPtr itsAsset;
		double itsWeight;
		double itsPrecision;
	public:
		ARM_PortfolioElem( ARM_VanillaArgPtr asset, double weight, double precision )
		:	itsAsset( asset), itsWeight(weight), itsPrecision(precision) {}
		double GetMktPrice() const;
		
		inline ARM_VanillaArgPtr GetAsset() const { return itsAsset; }
		inline double GetWeight() const { return itsWeight; }
		inline double GetPrecisiont() const { return itsPrecision; }
		ARM_PortfolioElem* Clone() const; 
	};
	
	vector<ARM_PortfolioElem*> itsElems;
	
	void CheckRange(size_t i) const
	{	if( i>= itsElems.size() )
	ARM_THROW( ERR_INVALID_ARGUMENT, "Out of Range in the portfolio!" );
	}
	
	/// for GetSubPortfolio and SplitPortfolioByMaturity easy construction
	ARM_VanillaPortfolio( const vector< ARM_PortfolioElem* >& elems )
	:	ARM_VanillaArg(), itsElems( elems ) {}
	
public:
	/// iterator support
	typedef CC_NS(std,vector)< ARM_PortfolioElem* >::const_iterator const_iterator;
	typedef CC_NS(std,vector)< ARM_PortfolioElem* >::iterator iterator;
	
	/// iterator support
	/// const version
	inline const_iterator begin() const { return itsElems.begin(); }
	inline const_iterator end() const { return itsElems.end(); }
	
	/// non const version
	inline iterator begin() { return itsElems.begin(); }
	inline iterator end() { return itsElems.end(); }
};

#undef CHECKRANGE

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
