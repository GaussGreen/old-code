/* -----------------------------------------------------------------
   -----------------------------------------------------------------

	FILE:		ICM_ADDRESSVECTOR.CPP
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


# include "icm_addressvector.h"
# include "icm_constants.h"



// -----------------------------------------------------------
// CONSTRUCTORS and DESTRUCTORS
// -----------------------------------------------------------

AddressVector :: AddressVector(IndexId NbEl)
{
	if (NbEl > 0)
	{
		nbel = NbEl;
		vector	=	(void* *)malloc(nbel*sizeof(void*));
		if (vector != NULL)
			memset(vector, 0, nbel*sizeof(void*));

	}
	else
	{
		nbel = 0;
		vector = NULL;
	}
}


AddressVector :: AddressVector(const AddressVector& Data)
{
	nbel = Data.nbel;

	vector = (void**)malloc(nbel * sizeof(void*));
	if (vector != NULL)
		memcpy(vector, Data.vector, nbel * sizeof(void*));
}


AddressVector :: ~AddressVector()
{
	Destroy();
}


void AddressVector :: Destroy()
{
	if (vector != NULL)
		free(vector);
}

void AddressVector :: Reset()
{
	Resize(0);
	vector = NULL;
}

ReturnCode AddressVector :: Resize(IndexId Nbel)
{
	IndexId nbelOld = nbel;
	if ((Nbel == 0) && (vector == NULL)) return RetOk;
	nbel = Nbel;
	vector = (void**)realloc(vector, nbel * sizeof(void*));
	if (vector == NULL) return RetNotOk;
	if (nbelOld < nbel)
		memset(vector+nbelOld, 0, (nbel-nbelOld)*sizeof(void*));
	return RetOk;
}


// -----------------------------------------------------------
// -----------------------------------------------------------

AddressVector& AddressVector :: operator  = (const AddressVector& other)
{
	if (this != &other)
	{
		if (other.nbel != nbel)
			Resize(other.nbel);
		memcpy(vector, other.vector, nbel * sizeof(void*));
	}

	return *this;
}



/*void AddressVector :: Set(IndexId Num, void* Data)
{
	if (Num >= nbel)
	{
		ICMTHROW(ERR_INVALID_ARGUMENT,"AddressVector :: Set: out of size"); 
		// return RetNotOk;
	}
	vector[Num] = Data;
	// return RetOk;
}
*/
//ReturnCode 
/*
void AddressVector :: Get(IndexId Num,void*& Data)
{
	if (Num >= nbel)
	{
		Data = 0;
		ICMTHROW(ERR_INVALID_ARGUMENT,"AddressVector :: Get: out of size"); 
		// return RetNotOk;
	}
	Data = vector[Num];
	// return RetOk;
}

*/
/*void* AddressVector :: Get(IndexId Num) const
{
	void* data;
	if (Num >= nbel)
		data = 0;
	data = vector[Num];
	return data;
}
*/
/** ReturnCode **/ 
/*
void AddressVector :: GetLast(void*& Data) const
{
	if (!nbel)
	{
		Data = 0;
		ICMTHROW(ERR_INVALID_ARGUMENT,"AddressVector :: GetLast: empty"); 
		// return RetNotOk;
	}
	Data = vector[nbel-1];
	//return RetOk;
}
*/

/* ReturnCode */ 
/*void * AddressVector :: GetLastMinusOne(void*& Data) const
{
	if (nbel < 2)
	{
		Data = 0;
		ICMTHROW(ERR_INVALID_ARGUMENT,"AddressVector :: GetLastMinusOne: 1 elt"); 
		// return RetNotOk;
	}
	Data = vector[nbel-2];
	// return RetOk;
}
*/
void* AddressVector :: GetLast() const
{
	void* data;
	if (!nbel)
		data = 0;
	data = vector[nbel-1];
	return data;
}

void* AddressVector :: GetLastMinusOne() const
{
	void* data;
	if (nbel < 2)
		data = 0;
	data = vector[nbel-2];
	return data;
}

ReturnCode AddressVector :: Append(IndexId Num, void* Data)
{
	ReturnCode Err;
	if (Num+1 > nbel)
	{
		Err = Resize(Num+1);
		if (Err != RetOk) return Err;
	}
	vector[Num] = Data;
	return RetOk;
}

ReturnCode AddressVector :: Append(void* Data)
{
	ReturnCode Err;
	Err = Resize(nbel+1);
	if (Err != RetOk) return Err;
	vector[nbel-1] = Data;
	return RetOk;
}

bool AddressVector :: Find(IndexId& Num, void* Data)
{
	void* * pointeur = vector;
	for (Num = 0; Num < nbel; Num++)
		if (Data == *(pointeur++)) return true;
	return false;
}


/** ReturnCode AddressVector :: RAZ(void* Data)
{
	IndexId i;
	void* * pointeur = vector;
	for (i = 0; i < nbel; i++)	*(pointeur++) = Data;
	return RetOk;
}**/ 
