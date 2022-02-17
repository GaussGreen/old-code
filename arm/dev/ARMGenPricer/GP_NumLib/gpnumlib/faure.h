/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file faure.h
 *
 *  \brief Faure low discrepancy sequence
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#ifndef _INGPNUMLIB_FAURE_H
#define _INGPNUMLIB_FAURE_H

#include "gpbase/port.h"
#include "gpbase/env.h"
#include "quasirandom.h"

CC_BEGIN_NAMESPACE( ARM )


//////////////////////////////////////////////////////////////////////////////////
///  Faure multidimentional low-discrepancy sequence generator
//////////////////////////////////////////////////////////////////////////////////
class ARM_FaureSeq : public ARM_QuasiRandom
{
public:
	ARM_FaureSeq(int firstSimulations);
	ARM_FaureSeq(const ARM_FaureSeq& rhs);
	ARM_FaureSeq& operator=(const ARM_FaureSeq& rhs );

	virtual ~ARM_FaureSeq();
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual void SetDim( size_t dim );

	/// standard ARM support
	virtual ARM_Object* Clone() const;

	/// validation function
	virtual void ValidateResult( double result ) {};


private:
	virtual void DrawAll();
	void Init();

	/// member variables for number generations!
	int itsQs;				/// prime number greater than the dimension
	int itsSeed;			/// seed
	int itsHisum_save;
	int itsHisum;
	vector<int> itsCoeff;
	vector<int> itsYtemp;

};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
