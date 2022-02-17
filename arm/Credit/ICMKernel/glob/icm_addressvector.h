/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_ADDRESSVECTOR.H
	PROJECT:	GLOB
	
	DESCRIPTION:	Vector of Pointor

  -----------------------------------------------------------------

 	CREATION:	October 13, 2004

	LAST MODIF:	October 13, 2004
  -----------------------------------------------------------------
   
	ICM Library

		version 1.0
		developped by Laurent Jacquel

  -----------------------------------------------------------------
  ----------------------------------------------------------------- */

# ifndef _ICM_ADDRESSVECTOR_H
# define _ICM_ADDRESSVECTOR_H

#include "ICMKernel\glob\icm_types.h"
#include "ICMKernel\glob\icm_enums.h"

class AddressVector
{

public:

	// Constructors and Destructors
	AddressVector(IndexId Data = 0);				// Constructor with Size
	AddressVector(const AddressVector& other);	// Copy Constructor
	~AddressVector();

private:

	void	Destroy();

public:

	void	Reset();

	void	GetCount(IndexId & Elts) const {Elts = nbel;}
	void**	GetVector() {return vector;}
	
	void*& operator [] (IndexId Id) 
	{
		if (Id>= nbel)
			ICMTHROW(ERR_INVALID_ARGUMENT,"AddressVector :: Get: out of size"); 
		return vector[Id];
	}
	void*& operator [] (IndexId Id) const 
	{
		if (Id>= nbel)
			ICMTHROW(ERR_INVALID_ARGUMENT,"AddressVector :: Get: out of size"); 
		return vector[Id];
	}

	ReturnCode	Resize(IndexId Nbel);

	// -----------------------------------------------------------
	// -----------------------------------------------------------
	// OPERATORS

	AddressVector& operator =	(const AddressVector& entree);	// assign
	
	// -----------------------------------------------------------
	// -----------------------------------------------------------
	
	void Set(IndexId Num, void* Data) { (*this)[Num]=Data; } 
	/*ReturnCode	*/ void Get(IndexId Num, void*& Data) { Data=(*this)[Num]; }
	/* ReturnCode	*/ void GetLast(void*& Data) const { Data= GetLast() ; }
	/* ReturnCode	*/ void GetLastMinusOne(void*& Data) const { Data= GetLastMinusOne();  }

	void*		Get(IndexId Num) const { return (*this)[Num]; } 
	void*		GetLast() const;
	void*		GetLastMinusOne() const;

	ReturnCode Append(IndexId Num,void* Data);
	ReturnCode Append(void* Data);
	
	// -----------------------------------------------------------
	// -----------------------------------------------------------
	
	bool Find(IndexId& Num,void* Data);

// 	ReturnCode RAZ(void* Data);

	
	// ***********************************************************
	// THE DATA
	// ***********************************************************

protected:

	IndexId	nbel;
	void*	*vector;

};

# endif
