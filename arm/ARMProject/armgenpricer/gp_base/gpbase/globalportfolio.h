/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file globalportfolio.h
 *
 *  \brief global portfolio will alow us to build 
 /// a global portfolio with all type of security ( ARM_Security, ARM_GC_Calculator, ARM_GenSecurity)
 *	\author  E. Ezzine
 *	\version 1.0
 *	\date May 2007
 */


#ifndef _INGPCALIB_GlOBALPORTFOLIO_H
#define _INGPCALIB_GlOBALPORTFOLIO_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "env.h"
#include "rootobject.h"

/// gpcalib

/// gpinfra
#include "port.h"
#include "gplinalgtypedef.h"
#include "assignop.h"

/// ARM Kernel
#include "expt.h"
#include "armglob.h"


/// STL
#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )

#ifdef __GP_STRICT_VALIDATION
	#define CHECKRANGE(i)  CheckRange(i);
#else
	#define CHECKRANGE(i)  
#endif


struct ARM_GlobalPortfolio: public ARM_RootObject
{
public:
	ARM_GlobalPortfolio(
		const vector<ARM_Object*>& assets,
        const std::vector<double>& weights,
        const std::vector<double>& mktprices);
	
	ARM_GlobalPortfolio( const ARM_GlobalPortfolio& rhs );
	ASSIGN_OPERATOR(ARM_GlobalPortfolio)
	virtual ~ARM_GlobalPortfolio();
	
	/// accessors
	inline ARM_Object* GetAsset( size_t i ) const { /*CHECKRANGE(i);*/ return itsElems[i]->GetAsset(); }
	inline double GetMktPrice( size_t i ) const { /*CHECKRANGE(i);*/ return itsElems[i]->GetMktPrice(); }
	inline double GetWeight( size_t i ) const { /*CHECKRANGE(i);*/ return itsElems[i]->GetWeight(); }
	inline size_t size() const { return itsElems.size(); }
	
	/// for the calibration validation
    bool HasSameAssets() const;
	
	/// vanilla arg support
	double Price(ARM_Object* model = NULL) const;
    double GetMktPrice( ) const;
	
	/// standard ARM support
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual ARM_Object* Clone() const { return new ARM_GlobalPortfolio(*this);}
	
private:
	/// forward declaration
	struct ARM_PortfolioElem
	{
	private:
		ARM_Object* itsAsset;
		double itsWeight;
		double itsPrice;
	public:
		ARM_PortfolioElem(const ARM_PortfolioElem& rhs);
		ASSIGN_OPERATOR(ARM_PortfolioElem)
		 ~ARM_PortfolioElem();

		ARM_PortfolioElem( const ARM_Object& asset, double weight, double price );
		double GetMktPrice() const;
		
		inline ARM_Object* GetAsset() const { return itsAsset; }
		inline double GetWeight() const { return itsWeight; }
		inline double GetPrice() const { return itsPrice; }
		ARM_PortfolioElem* Clone() const {	return new ARM_PortfolioElem(*this);} 
	};
	
	vector<ARM_PortfolioElem*> itsElems;
	
	void CheckRange(size_t i) const
	{	if( i>= itsElems.size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "Out of Range in the portfolio!" );
	}
	

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
