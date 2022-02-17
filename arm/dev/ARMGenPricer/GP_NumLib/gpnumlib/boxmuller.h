/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file boxmuller.h
 *
 *  \brief General file for the composite generator
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _INGPNUMLIB_BOXMULLER_H
#define _INGPNUMLIB_BOXMULLER_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "manipulator.h"
#include "typedef.h"


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////
/// Box Muller method
////////////////////////////////////////////////
struct ARM_BoxMuller : public ARM_ManipulatorMethod 
{
	ARM_BoxMuller(const ARM_RandomGeneratorPtr& randomGen );
	virtual double operator()() const;
	virtual string toString(const string& indent="", const string& nextIndent="") const;
//	virtual void reset( size_t dim, const ARM_GP_T_Vector<size_t>& nbOfPointsList, size_t factorsNb );

	/// copy constructor, assignment operator, destructor
	ARM_BoxMuller( const ARM_BoxMuller& rhs );
	ARM_BoxMuller& operator =( const ARM_BoxMuller& rhs );
	virtual ~ARM_BoxMuller();

	/// clone support
	virtual ARM_ManipulatorMethod* Clone() const;
//	ARM_RandomGeneratorPtr GetBaseGen() const { return itsRandomGen;}

protected:
	virtual void resetwork(const ARM_GP_T_Vector<size_t>& nbOfPointsList, const ARM_GP_T_Vector<size_t>& factorNb);

private:
//	ARM_RandomGeneratorPtr itsRandomGen;
	CC_IS_MUTABLE bool itsIset;
	CC_IS_MUTABLE double itsGset;
};



CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
