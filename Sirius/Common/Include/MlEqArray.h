//	MlEqArray.h :			 Array of doubles class
//
//	Author :				 David Cuin
/////////////////////////////////////////////////////////////////////////////

#ifndef _MLEQARRAY_H_
#define _MLEQARRAY_H_

#include "smart.h"
#include "cmatrix.h"

#undef min

class MlEqArray : public CVector, public RCObject
{
public:
	MlEqArray(void) {}
	MlEqArray(const CVector& a) : CVector(a) {}
	
	operator GDA::HDElement(void) const
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
	}
	MlEqArray& operator=(const CVector& v)
	{		
		*(dynamic_cast<CVector*>(this)) = v;						
		return *this;
	}

	// this is a preserving resize
	void PutSize(long nSize)
	{
		if (nSize == getsize()) return;
		CVector v(*this);
		resize(nSize);
		for (long n = 0; n < std::min(v.getsize(), getsize()); n++){
			(*this)[(int)n] = v[(int)n];
		}
	}
};

#endif