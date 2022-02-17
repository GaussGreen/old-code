#ifndef __drmatrix__h
#define __drmatrix__h

#include "drvalarray.h"

//k Simple matrix class.  Lacks a lot of functionality such as matrix addition or
//k multiplication.  Essentially, just represents data in 2D form with accessors.
template <class T> 
class DRMatrix {
public:
	DRMatrix (int length1=0, int length2=0); //
	virtual ~DRMatrix ();	//
	
	DRMatrix (const DRMatrix&); //
	DRMatrix& operator=(const DRMatrix&); //
	DRMatrix& operator=(const T&); //
	
	void resize (int, int); // Resizes the matrix.  Values in the matrix are undefined.
	int size (int dim) const; // Returns the size of the corresponding dimension
	
	operator T**();						// Cast into a T**
	operator T**() const;				// Cast into a T**
	
	DRValarray<T> get_array (int loc, int dim) const; // Returns a 1D array
	// when dim = 1 or 2, and loc tells you which row/col
	// eg. get_array(0, 1) means set the dim=1 (row) to be 0
	// In other words, get the first column

	DRMatrix& store_array (const DRValarray<T>&, int loc, int dim); // Stores 1D Array 
	// same notation as get_array
	
	friend ostream &operator<<(ostream& s, const DRMatrix& a) //
	{
		s << "Matrix(" << a.m_length1 << "," << a.m_length2 << ")" << endl;
		
		for (int i = 0; i< a.m_length1; i++) {
			for (int j = 0; j< a.m_length2; j++) s << a.m_ptr[i][j] << " ";
			s << endl;
		}
		s << endl;	
		return s;
	}
	
protected:
	void Alloc (int length1, int length2);
	void Free ();
	int m_length1;
	int m_length2;
	T** m_ptr;
};

// Common typedefs
typedef DRMatrix<double>	DMatrix;


// Member function implementation
template <class T> inline 
DRMatrix<T>::DRMatrix (int length1, int length2)
{
	Alloc(length1, length2);
}

template <class T> inline 
DRMatrix<T>::~DRMatrix ()
{
	if (m_ptr) Free();
}

template <class T> inline 
DRMatrix<T>::DRMatrix (const DRMatrix& a)
{
	Alloc (a.m_length1, a.m_length2);
	
	for (int i = 0 ; i < m_length1; i++)
		for (int j = 0; j < m_length2; j++) m_ptr[i][j] = a.m_ptr[i][j];
}

template <class T> inline 
DRMatrix<T>& DRMatrix<T>::operator= (const DRMatrix<T>& a)
{
	if (this == &a) return *this;
	
	resize(a.m_length1, a.m_length2);
	
	for (int i = 0 ; i < m_length1; i++) {
		for (int j = 0; j < m_length2; j++) 
			m_ptr[i][j] = a.m_ptr[i][j];
	}
	
	return *this;
}

template <class T> inline 
DRMatrix<T>& DRMatrix<T>::operator= (const T& a)
{
	for (int i = 0 ; i < m_length1; i++) {
		for (int j = 0; j < m_length2; j++) 
			m_ptr[i][j] = a;
	}
	
	return *this;
}

template <class T> inline 
void DRMatrix<T>::resize (int length1, int length2)
{
	if (length1 != m_length1 || length2 != m_length2) {
		Free();		
		Alloc(length1, length2);
	}
}

template <class T> inline 
int DRMatrix<T>::size (int dim) const
{
	if (dim == 1) return m_length1;
	if (dim == 2) return m_length2;
	else throw DRException ("Illegal matrix dimension");
	return 0;
}

template <class T> inline 
DRMatrix<T>::operator T**() const
{
	return m_ptr;
}

template <class T> inline 
DRMatrix<T>::operator T**() 
{
	return m_ptr;
}

template <class T> inline 
DRValarray<T> DRMatrix<T>::get_array (int loc, int dim) const
{
	if (dim != 1 && dim != 2) 
		throw DRException ("Illegal dimension ") << dim;

	int len;
	if (dim == 1) len = m_length2;
	else len = m_length1;

	DRValarray<T> ans(len);
	for (int i=0; i< len; i++) {
		if (dim == 1) {
			ans[i] = m_ptr[loc][i];
		}
		else {
			ans[i] = m_ptr[i][loc];
		}
	}

	return ans;
}

template <class T> inline 
DRMatrix<T>& DRMatrix<T>::store_array (const DRValarray<T>& a, int loc, int dim) 
{
	if (dim != 1 && dim != 2) 
		throw DRException ("Illegal dimension ") << dim;

	int len;
	if (dim == 1) len = m_length2;
	else len = m_length1;

	if (len != a.size())
		throw DRException ("Wrong-size array");

	for (int i = 0; i< len; i++) {
		if (dim == 1) {
			m_ptr[loc][i] = a[i];
		}
		else {
			m_ptr[i][loc] = a[i];
		}
	}

	return *this;
}


template <class T> inline 
void DRMatrix<T>::Alloc (int length1, int length2)
{
	m_length1 = length1;
	m_length2 = length2;
	
	if (length1 < 0 || length2 < 0)
		throw DRException ("DRMatrix can't make an array with a negative dimension");
	
	if (m_length1 == 0) {
		m_ptr = NULL;
	}
	else {
		m_ptr = new T* [m_length1];
		
		if (!m_ptr) 
			throw DRException ("Memory error in DRMatrix");
		
		for (int i=0; i< m_length1; i++) {
			if (m_length2 == 0) {
				m_ptr[i] = NULL;
			}
			else {
				m_ptr[i]=new T [m_length2]; 
				if (!m_ptr[i]) 
					throw DRException ("Memory error in DRMatrix");
			}
		}
	}
}

template <class T> inline 
void DRMatrix<T>::Free ()
{
	if (m_ptr) 
	{
		for (int i=0; i<m_length1; i++) {
			if (m_ptr[i]) delete [] m_ptr[i];
		}

		delete [] m_ptr;
	}
}

#endif

