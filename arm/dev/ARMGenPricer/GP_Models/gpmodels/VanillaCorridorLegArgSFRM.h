/*
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 * $Log: VanillaCorridorLegArgSFRM.h,v $
 * Revision 1.1  2004/05/06 10:06:19  emezzine
 * Initial revision
 *
 */

/*! \file VanillaCorridorLegArgSFRM.h
 *
 *  \brief this files enables to store various data to fast calibrate
 *      caps for the SFRM model
 *	\author  E.Ezzine
 *	\version 1.0
 *	\date May 2004
 */


#ifndef _INGPMODEL_VANILLACORRIDORLEGARGSFRM_H
#define _INGPMODEL_VANILLACORRIDORLEGARGSFRM_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix

#include "gpinfra/typedef.h"
#include "gpcalib/vanillacorridor.h"

#include "gpbase/port.h"
#include <vector>
CC_USING_NS(std,vector)


CC_BEGIN_NAMESPACE( ARM )

class ARM_SFRM;

///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaCorridorLegArgSFRM
/// \brief swaption argument simple struct
///////////////////////////////////////////////////////////////
struct ARM_VanillaCorridorLegArgSFRM: public ARM_VanillaCorridorLegArg
{
	ARM_VanillaCorridorLegArgSFRM();

	ARM_VanillaCorridorLegArgSFRM( const ARM_VanillaCorridorLegArg& rhs );

	/// copy constructor
	ARM_VanillaCorridorLegArgSFRM( const ARM_VanillaCorridorLegArgSFRM& rhs );
    ARM_VanillaCorridorLegArgSFRM& operator=(const ARM_VanillaCorridorLegArgSFRM& rhs);
	virtual ~ARM_VanillaCorridorLegArgSFRM();
	virtual ARM_Object* Clone();

	/// accessors
	inline void SetLibors( const vector<ARM_VectorPtr>& libors )    { itsLibors = libors; }
	inline void SetZCPays( const vector<ARM_VectorPtr>& ZCPays )    { itsZCPays = ZCPays; }

    inline void SetRef_Libors( const vector<ARM_VectorPtrVector>& ref_libors )    { itsRef_Libors = ref_libors;  }
	inline void SetRef_ZCEnds( const vector<ARM_VectorPtrVector>& ref_ZCPays )    { itsRef_ZCEnds = ref_ZCPays;  }
	inline void SetRef_DfFwds( const vector<ARM_VectorPtrVector>& ref_ZCEnds)     { itsRef_DfFwds = ref_ZCEnds;  }

private :
	friend class ARM_SFRM;
	ARM_VectorPtrVector   itsLibors;
	ARM_VectorPtrVector   itsZCPays;

	vector< ARM_VectorPtrVector >  itsRef_Libors;
    vector< ARM_VectorPtrVector >  itsRef_DfFwds;
	vector< ARM_VectorPtrVector >  itsRef_ZCEnds;

    void CopyNoCleanUp(const ARM_VanillaCorridorLegArgSFRM& rhs);
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
