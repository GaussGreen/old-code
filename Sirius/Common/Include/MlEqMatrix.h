//	MlEqMatrix.h :			 Matrix of doubles class
//
//	Author :				 David Cuin
/////////////////////////////////////////////////////////////////////////////

#ifndef _MLEQMATRIX_H_
#define _MLEQMATRIX_H_

#include "smart.h"
#include "cmatrix.h"

class MlEqMatrix : public CMatrix, public RCObject
{
public:
	/*operator GDA::HDElement(void) const
	{
		GDA::HDElement hdeArray(getsize(), 1);
		for (long n = 0; n < getsize(); n++){
			hdeArray(n, 0) = (*this)[(int)n];
		}			
		return hdeArray;
	}
			
	// ToDo - replace all this once the compiler supports templatised operators.
	// (the bool cast may be an exception)
	operator GVector<int>(void) const
	{
		GVector<int> gv(getsize());
		for (long n = 0; n < getsize(); n++){
			gv[n] = (int)(*this)[n];
		}		
		return gv;
	}
	operator GVector<long>(void) const
	{
		GVector<long> gv(getsize());
		for (long n = 0; n < getsize(); n++){
			gv[n] = (long)(*this)[n];
		}		
		return gv;
	}
	operator GVector<bool>(void) const
	{
		GVector<bool> gv(getsize());
		for (long n = 0; n < getsize(); n++){
			gv[n] = (*this)[n] ? true : false;
		}		
		return gv;
	}*/
	
	MlEqMatrix& operator=(const CMatrix& m)
	{		
		*(dynamic_cast<CMatrix*>(this)) = m;
		return *this;
	}
};

#endif