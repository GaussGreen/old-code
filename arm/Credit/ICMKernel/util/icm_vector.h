
#ifndef _ICM_VECTOR_H_
#define _ICM_VECTOR_H_

/**
	This class acts as a std::vector<> with controlled access
**/ 

#include "ICMKernel/util/icm_macro.h"
#include <vector>

template <class T> 
class ICM_Vector : public ARM_Object
{
	std::vector<T> itsValues; 
public:
	ICM_Vector()
	{}
	ICM_Vector(unsigned int size) : itsValues(size) 
	{}
	virtual ~ICM_Vector() 
	{}
	ICM_Vector(const ICM_Vector&ref) : itsValues(ref.itsValues) 
	{}
	ICM_Vector& operator=(const ICM_Vector&ref)
	{
		if (this!=&ref) 
		{
			this->~ICM_Vector(); 
			new(this)ICM_Vector(ref); 
		}
		return *this; 
	}
	virtual ARM_Object* Clone() 
	{ 
		return new ICM_Vector(*this); 
	}
	//
	const T& operator[](unsigned int i) const
	{
		if (i>=itsValues.size()) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Vector("<<i<<"): out of bounds"); 
		return itsValues[i]; 
	}
	T& operator[](unsigned int i)
	{
		if (i>=itsValues.size()) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Vector("<<i<<"): out of bounds"); 
		return itsValues[i]; 
	}
	unsigned int size() const
	{
		return itsValues.size(); 
	}
	void resize(unsigned int size)
	{
		itsValues.resize(size); 
	}
	void erase(unsigned int i) 
	{
		if (i>=itsValues.size()) 
			ICMTHROW(ERR_INVALID_ARGUMENT,"ICM_Vector::remove("<<i<<"): out of bounds"); 
		itsValues.erase(itsValues.begin()+i); 
	}
	void push_back(const T&item)
	{
		itsValues.push_back(item) ;
	}
} ;

#endif //_ICM_VECTOR_H_
