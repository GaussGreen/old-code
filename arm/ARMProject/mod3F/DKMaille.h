#ifndef __DG_DKMAILLE_H__
#define __DG_DKMAILLE_H__

#include <stdlib.h>
#include <stdio.h>

#pragma warning(disable : 4251)
#pragma warning(default : 4251)

unsigned long Index3dToIndex1d(unsigned long uli,unsigned long ulj,unsigned long ulk,
                               unsigned long ulIMAX,unsigned long ulJMAX,unsigned long ulKMAX);
                               


template <class T>
class DKMaille 
{
  public:
  DKMaille();
  DKMaille(unsigned N);
  DKMaille(unsigned noElements, const T *array);
  DKMaille(unsigned noElements, T init);
  DKMaille(const DKMaille<T>& c);
  DKMaille& operator=(const DKMaille<T>& c);

  void print(FILE* fic);

  ~DKMaille() { dealloc(); }

#ifdef DG_FULL_CHECKING
  const T& at(unsigned i, int dummy=0) const 
  { 
    if (  (dummy == 0) && ((i>=entries()) || (i<0) ))
    {
      throw("ERROR: DKMaille write error in const at(...)");
    }
    else
      return myData[i]; 
  }

  T& at(unsigned i, int dummy=0) 
  { 
    if ((dummy == 0) && ((i >= entries() ) || ( i < 0 ) ))
    {
      throw("ERROR: DKMaille write error in at(...)");
    }
    else
      return myData[i]; 
  }
  const T& operator()(unsigned i) const 
  { 
    if ((i >= entries() ) || ( i < 0 ) )
    {
      throw("ERROR: DKMaille write error in (...)");
    }
    else
      return myData[i]; 
  }

  T& operator()(unsigned i)
  { 
    if ((i >= entries() ) || ( i < 0 ) )
    {
      throw("ERROR: DKMaille write error in (...)");
    }
    else
      return myData[i]; 
  }

  const T& operator[](unsigned i) const 
  { 
    if ((i >= entries() ) || ( i < 0 ) )
    {
      throw("ERROR: DKMaille write error in []");
    }
    else
      return myData[i]; 
  }

  T& operator[](unsigned i)
  { 
    if ((i >= entries() ) || ( i < 0 ) )
    {
      throw("ERROR: DKMaille write error in []");
    }
    else
      return myData[i]; 
  }

#else
  const T& at(unsigned i, int dummy=0) const { return myData[i]; }
  T& at(unsigned i, int dummy=0) { return myData[i]; }
  const T& operator()(unsigned i) const { return myData[i]; }
  T& operator()(unsigned i) { return myData[i]; }
  const T& operator[](unsigned i) const { return myData[i]; }
  T& operator[](unsigned i) { return myData[i]; }
#endif 

  void removeAt(unsigned p)
  {
    --myEntries;
    for(unsigned i = p; i < entries(); ++i)
    {
      myData[i] = myData[i+1];
    }
  }

  void removeAll(const T& t)
  {
    while(contains(t))
    {
      remove(t);
    }
  }

  void remove(const T& t)
  {
    int i = index(t);
    if(i != -1)
    {
      removeAt((unsigned)i);
    }
  }
  
  
  void removeBefore(const T& t)
  {
    sort();
    unsigned i = 0;

    while( i < entries() -1)
    {
      if (at(i,999) < t)
      {
        removeAt(i);
      }
      else
      {
        ++i;
      }
    }
  }

 
  bool operator==(const DKMaille<T>& rhs) const { return this == &rhs; }
  const T* asCArray() const;
  void clearAndDestroy();
  // We need to add the qsort to this
  void sort();
  void SortWithCompare( int CompareFunc(const void *pvElem1, const void *pvElem2 ) );
  void combine(const DKMaille<T>& in);
  void removeDuplicates();
  void fill(const T& in);

  void clear()
  {
    dealloc();
    myEntries = 0;
    mySize = 0;
  }
  void dealloc();
  T *alloc(unsigned N);

#ifdef DG_ARRAY_GLOBAL
  T *allocG(unsigned N);
  void deallocG();
#elif defined(DG_ARRAY_BASIC_GLOBAL)
  T *allocB(unsigned N);
  void deallocB();
#else
  T *allocL(unsigned N);
  void deallocL();
#endif


  const T* data() const { return myData; }

  void grow(size_t N)
  {
    if(N <= mySize)
    {
      return;
    }
    else
    {
      T *replacement = alloc(N);
      for(unsigned i = 0; i < entries(); ++i)
      {
        replacement[i] = myData[i];
      }
      dealloc();
      mySize = N;
      myData = replacement;
    }
  }


  // Take a particular datum for initialising the array
  void grow(size_t N, T t)
  {
    if(N <= mySize)
    {
      return;
    }
    else
    {
      T *replacement = alloc(N);
      unsigned iLimit = entries();
      for(unsigned ii = 0; ii < N; ++ii)
      {
        if (ii < iLimit)
          replacement[ii] = myData[ii];
        else
          replacement[ii] = t;
      }
      dealloc();
      mySize = N;
      myData = replacement;
    }
  } // grow 2-ary



  void  resize(size_t N)
  {
    if(N > mySize)
    {
      grow(N);
    }
    myEntries = N;
  }


  
  // Take a particular datum for initialising the array
  void  resize(size_t N, T t)
  {
    if(N > mySize)
    {
      grow(N, t);
    }
    myEntries = N;
  } // resize 2-ary



  void append(const T& t)
  {
    if(myEntries == mySize)
    {
      grow(mySize * 2 + 1);
    }
    at(myEntries,999) = t;
    ++myEntries;
  }

  
  void insert(const T& t) { append(t); }

  unsigned length() const { return myEntries; }
  unsigned entries() const { return myEntries; }
  int contains(const T& t) const
  {
    for(unsigned i = 0; i < entries(); ++i)
    {
      if(t == at(i,999)) return 1;
    }
    return 0;
  }

  int index(const T& t) const
  {
    for(unsigned i = 0; i < entries(); ++i)
    {
      if(t == at(i,999)) return i;
    }
    return -1;
  }

  int bsearch(const T& value, bool exact = false) const;

  T getStddev() const;
  T getAverage() const;
  const T& getMax() const;
  const T& getMin() const;
  T getSum() const;

  static int compare(const void* elem1, const void* elem2);

  DKMaille& operator-=(const DKMaille& rhs);
  DKMaille& operator+=(const DKMaille& rhs);
  DKMaille& operator*=(const T& rhs);

  T dot(const DKMaille& rhs);

protected:
  T *myData;
  unsigned mySize;
  unsigned myEntries;
};


// RWDECLARE_PERSIST_TEMPLATE_IO(DKMaille,RWvistream,RWvostream) // we don't want the RWFile stuff


//----------------------------------------------------------------------------------
template <class T>
inline DKMaille<T>::DKMaille() : myData(NULL), mySize(0), myEntries(0)
{
}
//----------------------------------------------------------------------------------

//----------------------------------------------------------------------------------
template <class T>
inline DKMaille<T>::DKMaille(unsigned N) 
 : myData(alloc(N)), mySize(N), myEntries(N)
{
}
//----------------------------------------------------------------------------------

//----------------------------------------------------------------------------------
template <class T>
DKMaille<T>::DKMaille(const DKMaille<T>& c) 
  : myData(alloc(c.mySize)), mySize(c.mySize), myEntries(c.myEntries)
{
  for(unsigned i = 0; i < entries(); ++i)
  {
    myData[i] = c.myData[i];
  }
}
//----------------------------------------------------------------------------------

//----------------------------------------------------------------------------------
template <class T>
DKMaille<T>& DKMaille<T>::operator=(const DKMaille<T>& c)
{
  dealloc();
  myData    = alloc(c.mySize);
  mySize    = c.mySize;
  myEntries = c.myEntries;
  for(unsigned i = 0; i < entries(); ++i)
  {
    myData[i] = c.myData[i];
  }
  return *this;
}
//----------------------------------------------------------------------------------

//----------------------------------------------------------------------------------
template <class T>
inline DKMaille<T>::DKMaille(unsigned noElements, const T *array) 
  : myData(alloc(noElements)), mySize(noElements), myEntries(noElements)
{
  for(unsigned i = 0; i < entries(); ++i)
  {
    myData[i] = array[i];
  }
}
//----------------------------------------------------------------------------------


//----------------------------------------------------------------------------------
template <class T>
inline DKMaille<T>::DKMaille(unsigned N, T init) 
 : myData(alloc(N)), mySize(N), myEntries(N)  
{
  for(unsigned i = 0; i < N; ++i)
  {
    myData[i] = init;
  }
}
//----------------------------------------------------------------------------------


//----------------------------------------------------------------------------------
template <class T>
void DKMaille<T>::dealloc()
{ 
#ifdef DG_ARRAY_GLOBAL
  deallocG();
#elif  defined(DG_ARRAY_BASIC_GLOBAL)
  deallocB();
#else
  deallocL();
#endif
}
//----------------------------------------------------------------------------------

//----------------------------------------------------------------------------------
template <class T>
T *DKMaille<T>::alloc(unsigned N)
{
#ifdef DG_ARRAY_GLOBAL
  return allocG(N);
#elif  defined(DG_ARRAY_BASIC_GLOBAL)
  return allocB(N);
#else
  return allocL(N);
#endif

}
//----------------------------------------------------------------------------------

#ifdef DG_ARRAY_GLOBAL
//----------------------------------------------------------------------------------
template <class T>
T *DKMaille<T>::allocG(unsigned N)
{
  if(N == 0) return NULL;

  T *ret = (T*)mt_malloc(N * sizeof(T));
  for(unsigned i = 0; i < N; ++i)
  {
#pragma warning(disable : 4291)
    T *tmp = new ((void*)(ret + i)) T;   
#pragma warning(default : 4291)
  }
  return ret;
}
//----------------------------------------------------------------------------------
#elif  defined(DG_ARRAY_BASIC_GLOBAL)
//----------------------------------------------------------------------------------
template <class T>
T *DKMaille<T>::allocB(unsigned N)
{
  return CAST_MALLOC(T, N);
}
//----------------------------------------------------------------------------------
#else
//----------------------------------------------------------------------------------
#ifdef DG_FULL_CHECKING
template <class T>
T *DKMaille<T>::allocL(unsigned N)
{
  T *pT = NULL;
  try
  {
  pT = new T[N];
  }
  catch(...)
  {
    throw("ERROR: memory allocation failed");
  }

  return pT;
}
#else // it's release
template <class T>
T *DKMaille<T>::allocL(unsigned N)
{
  return new T[N];
}
#endif // DG_FULL_CHECKING
//----------------------------------------------------------------------------------
#endif


#ifdef DG_ARRAY_GLOBAL
//----------------------------------------------------------------------------------
template <class T>
void DKMaille<T>::deallocG()
{ 
  if(myData == NULL) return;

  for(unsigned i = 0; i < mySize; ++i)
  {
    (myData+i)->~T();     
  }
  mt_free((void*)myData);
}
//----------------------------------------------------------------------------------
#elif  defined(DG_ARRAY_BASIC_GLOBAL)
//----------------------------------------------------------------------------------
template <class T>
void DKMaille<T>::deallocB()
{
  CAST_FREE(myData);
}
//----------------------------------------------------------------------------------
#else
//----------------------------------------------------------------------------------
template <class T>
void DKMaille<T>::deallocL()
{ 
/*
  if (myData)
    delete[] myData; 
    */
  if (mySize) // RDJ 04/02/02 - just deleting doesn't set to NULL - check size instead!
    delete[] myData; 

}
//----------------------------------------------------------------------------------
#endif


//----------------------------------------------------------------------------------
template <class T>
void inline DKMaille<T>::clearAndDestroy()
{
  for(unsigned i = 0; i < entries(); i++)
  {
    delete at(i,999);
  }
  dealloc();
  mySize = 0;
  myEntries = 0;
}

//----------------------------------------------------------------------------------
template <class T>
inline const T *DKMaille<T>::asCArray() const
{
	return data();
}

//----------------------------------------------------------------------------------

template <class T>
void DKMaille<T>::sort()
{
  qsort((void*)data(), myEntries, sizeof(T), DKMaille::compare);
}



//----------------------------------------------------------------------------------

template <class T>
void DKMaille<T>::SortWithCompare( int CompareFunc(const void *pvElem1, const void *pvElem2 ) )
{
  qsort((void*)data(), myEntries, sizeof(T), CompareFunc);
}


//----------------------------------------------------------------------------------
template <class T>
void DKMaille<T>::fill(const T& in)
{
  for (unsigned i = 0; i < entries(); ++i)
  {
    at(i,999) = in;
  }
}
//----------------------------------------------------------------------------------
// append the array to the end of this one
template <class T>
void DKMaille<T>::combine(const DKMaille<T>& in)
{
  unsigned oldN = entries();
  unsigned inN  = in.entries();

  resize(oldN + inN);

  unsigned i = oldN;
  unsigned j = 0;

  for (; j < inN; ++i, ++j)
  {
    at(i,999) = in(j);
  }
}


//----------------------------------------------------------------------------------
// Remove Duplicates - sorts the array in the process
template <class T>
void DKMaille<T>::removeDuplicates()
{
  sort();

  unsigned i = 0;

  while( i < entries() -1)
  {
    if (at(i,999) == at(i+1,999))
    {
      removeAt(i+1);
    }
    else
    {
      ++i;
    }
  }
}

//----------------------------------------------------------------------------------
template <class T>
int DKMaille<T>::compare(const void* elem1, const void* elem2)
{
  T *d1 = (T*)elem1;
  T *d2 = (T*)elem2;
  if(*d1 < *d2)
  {
    return -1;
  }
  else
  if(*d2 < *d1)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}


//----------------------------------------------------------------------------------
template <class T>
int DKMaille<T>::bsearch(const T& value, bool exact) const
{
  // binary search to find the value in the y grid
  int top = entries()-1, bot = 0, mid;

  if ((at(top,999) < value || value < at(bot,999)) == false)       // check we've been given a reasonable value
  {
	  while ( (top - bot) > 1)   
	  {
	    mid = (top + bot)/2;
      if ( at(mid,999) < value)
      {
		    bot = mid;
      }
	    else 
      if ( value < at(mid,999))
      {
		    top = mid;
      }
	    else           // we've got an exact match 
      {
		    top = bot = mid;
        break;
      }  
    }

	// check for the special case where value is the last value in array
	if (exact == true && value == at(entries()-1,999))
	{
		bot = entries()-1;
	}

    return bot;
	}
  else                   // error checking
	{
    return RW_NPOS;
  }
}

/*
template <class T>
void DKMaille<T>::print(const RWCString &title) const
{
    std::cout << std::endl;
    std::cout << title << std::endl;
    for(unsigned i = 0; i < entries(); ++i)
    {
	    std::cout << "\t" << i << "\t" << at(i,999) << std::endl;
	}
}
*/
template <class T>
T DKMaille<T>::getStddev() const
{
	T	ave  = 0.;
	T	vari = 0.;

	for (unsigned i=0; i < entries(); ++i)	{
		ave  += at(i,999);
		vari += at(i,999)*at(i,999);
	}

	ave  = ave /(double)entries();
	vari = vari/(double)entries();

	return 	sqrt(vari - ave*ave); // What's going on if "T" != double ??????
}

// --------------------------------------------------------------------------------
template <class T>
T DKMaille<T>::getAverage() const
{
	T	ave  = 0.;
	for (unsigned i=0; i < entries(); ++i)	{
		ave  += at(i,999);
	}
	return (ave /(double)entries());
}

// --------------------------------------------------------------------------------
template <class T>
const T& DKMaille<T>::getMax() const
{
	const T *high = &at(0,999);
	if (entries() > 0)
	{
		for (unsigned i = 1; i < entries(); ++i)
		{
			if (*high < at(i,999))
			{
				high = &at(i,999);
			}
		}
	}
	return *high;
}

// --------------------------------------------------------------------------------
template <class T>
const T& DKMaille<T>::getMin() const
{
	const T *low = &at(0,999);
	if (entries() > 0)
	{
		for (unsigned i = 1; i <entries(); ++i)
		{
			if (at(i,999) < *low)
			{
				low = &at(i,999);
			}
		}
	}
	return *low;
}

// --------------------------------------------------------------------------------
template <class T>
T DKMaille<T>::getSum() const
{
	T sum = at(0,999);
	for (unsigned i = 1; i <entries(); ++i)
	{
		sum += at(i,999);
	}
	return sum;
}

#ifdef WIN32
// --------------------------------------------------------------------------------
// --------------------------------------------------------------------------------
// Global Helper Functions

template <class Ret, class In>
DKMaille<Ret> convertArray(const DKMaille<In>& inArray)
{
	DKMaille<Ret> retArray(inArray.entries());
	for (int i=0; i < inArray.entries(); ++i) 
	{
		retArray(i) = inArray(i);
	}
	return retArray;
}
#endif

template <class T>
DKMaille<T> elementMultiply(const DKMaille<T>& lhs, const DKMaille<T>& rhs)
{
  if(lhs.entries() != rhs.entries())
  {
    MESSAGE("Cannot multiply vectors by element ");
    msg.stream() << "lhs.entries() = " << lhs.entries() << "\t";
    msg.stream() << "rhs.entries() = " << rhs.entries();
    throw msg;
  }
  DKMaille<T> retArray = lhs;
  for(unsigned i = 0; i < retArray.entries(); ++i)
  {
    retArray(i) *= rhs(i);
  }
  return retArray;
}

//================================================================================
template <class T>
DKMaille<T>& DKMaille<T>::operator-=(const DKMaille& rhs)
{
  if(entries() != rhs.entries())
  {
    MESSAGE("Cannot subtract vectors by element ");
    msg.stream() << "entries() = " << entries() << "\t";
    msg.stream() << "rhs.entries() = " << rhs.entries();
    throw msg;
  }
  for(unsigned i = 0; i < entries(); ++i)
  {
     operator()(i) -= rhs(i);
  }
  return *this;
}

//================================================================================
template <class T>
DKMaille<T>& DKMaille<T>::operator+=(const DKMaille& rhs)
{
  if(entries() != rhs.entries())
  {
    MESSAGE("Cannot add vectors by element ");
    msg.stream() << "entries() = " << entries() << "\t";
    msg.stream() << "rhs.entries() = " << rhs.entries();
    throw msg;
  }
  for(unsigned i = 0; i < entries(); ++i)
  {
     operator()(i) += rhs(i);
  }
  return *this;
}

//================================================================================
template <class T>
DKMaille<T>& DKMaille<T>::operator*=(const T& rhs)
{
	for (unsigned i = 0; i < entries(); ++i)
	{
		at(i,999) *= rhs;
	}

	return *this;
}

//================================================================================
template <class T>
T DKMaille<T>::dot(const DKMaille& rhs)
{
  if(entries() != rhs.entries())
  {
    MESSAGE("Cannot multiply vectors by element ");
    msg.stream() << "entries() = " << entries() << "\t";
    msg.stream() << "rhs.entries() = " << rhs.entries();
    throw msg;
  }

  T ret(0.0);
  for(unsigned i = 0; i < entries(); ++i)
  {
     ret += operator()(i) * rhs(i);
  }
  return ret;
}

//================================================================================
template <class T>
DKMaille<T> operator+(const DKMaille<T>& lhs, const DKMaille<T>& rhs)
{
  if(lhs.entries() != rhs.entries())
  {
    MESSAGE("Cannot Add vectors by element ");
    msg.stream() << "lhs.entries() = " << lhs.entries() << "\t";
    msg.stream() << "rhs.entries() = " << rhs.entries();
    throw msg;
  }
  DKMaille<T> retArray(lhs);
  retArray += rhs;
  return retArray;
}

//================================================================================
template <class T>
DKMaille<T> operator-(const DKMaille<T>& lhs, const DKMaille<T>& rhs)
{
  if(lhs.entries() != rhs.entries())
  {
    MESSAGE("Cannot multiply vectors by element ");
    msg.stream() << "lhs.entries() = " << lhs.entries() << "\t";
    msg.stream() << "rhs.entries() = " << rhs.entries();
    throw msg;
  }
  DKMaille<T> retArray(lhs);
  retArray -= rhs;
  return retArray;
}

//================================================================================
template <class T>
DKMaille<T> operator* (const T& lhs, const DKMaille<T>& rhs)
{
	DKMaille<T> tmp(rhs.entries());
	for (unsigned i = 0; i < tmp.entries(); ++i)
	{
		tmp(i) = lhs * rhs(i);
	}

	return tmp;
}

//================================================================================
template <class T>
DKMaille<T> operator* (const DKMaille<T>& lhs, const T& rhs)
{
	DKMaille<T> tmp = lhs;
	return tmp *= rhs;
}


//================================================================================
template <class Type>
void DKMaille<Type>::print(FILE* fic)
{
	
	unsigned int i;

	fprintf(fic, "\n");

	for(i=0;i<entries();i++)
	{
  		fprintf(fic, "%lf ", at(i));
		
		fprintf(fic, "\n");
	}
}


#endif

