// DRValarray.h: interface for the drvalarray class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_DRVALARRAY_H__14A1C2E1_226B_11D2_97B2_00C04FD8EB9A__INCLUDED_)
#define AFX_DRVALARRAY_H__14A1C2E1_226B_11D2_97B2_00C04FD8EB9A__INCLUDED_

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000

#include "drexception.h"
#include IOSTREAM_H

//k  Emulates STL valarray without the slices.

template <class T>
class DRValarray  
{
public:
	DRValarray (); //
	explicit DRValarray (size_t n); //
	DRValarray (const T&, size_t); //
	DRValarray (const T*, size_t); //
	virtual ~DRValarray(); //
	
	DRValarray (const DRValarray&);		// copy constructor
	
	operator T*();	//
	operator const T*() const; //
	
	DRValarray& operator=(const DRValarray&);	// = operator 
	DRValarray& operator=(const T&); // same as fill
	
	int size() const; //
	void resize(int size); //
	
	DRValarray apply (T fn(T)) const; //
	DRValarray apply (T fn(const T&)) const; //
	
	T min() const;		// Gets smallest element
	T max() const;		// Gets biggest element
	T sum() const;		// Computes the sum of the array
	T mean() const;		// Computes the mean
	double variance () const; // Computes the variance and returns a double
	
	friend T dot_product (const DRValarray<T>& x, const DRValarray<T>& y) // Computes dot product
	{
		x.CheckSameSize(y); 
		T sum = 0;
		
		for (int i = 0; i < x.m_length; i++) 
			sum += x.m_ptr[i] * y.m_ptr[i];
		
		return sum;
	}
	
	friend double covariance (const DRValarray<T>&, const DRValarray<T>&) //
	{throw DRException ("Never been used"); return 0;}
	
	friend double correlation (const DRValarray<T>&, const DRValarray<T>&) //
	{throw DRException ("Never been used"); return 0;}
	
	friend ostream &operator<<(ostream& s, const DRValarray &a) //
	{
		s << "Length(" << a.m_length << "): " <<endl;
		
		for (int i = 0; i< a.m_length; i++) 
			s << a.m_ptr[i] << " ";
		
		s << endl; 
		return s;
	}
	
	// the math operations +,-,/,*,+=,-=,/=,*= act only on arrays of the same size
	DRValarray operator/(const DRValarray&) const; //
	DRValarray operator/(const T&) const; //
	DRValarray& operator/=(const DRValarray&); //
	DRValarray& operator/=(const T&); //
	
	DRValarray operator*(const DRValarray&) const; //
	DRValarray operator*(const T&) const; //
	DRValarray& operator*=(const DRValarray&); //
	DRValarray& operator*=(const T&); //
	
	DRValarray operator+(const DRValarray&) const; //
	DRValarray operator+(const T&) const; //
	DRValarray& operator+=(const DRValarray&); //
	DRValarray& operator+=(const T&); //
	
	DRValarray operator-(const DRValarray&) const; //
	DRValarray operator-(const T&) const; //
	DRValarray& operator-=(const DRValarray&); //
	DRValarray& operator-=(const T&); //

	bool operator==(const DRValarray&) const;

	void print_short (ostream&) const; //
	
protected:
	int m_length, m_capacity;
	T* m_ptr;
	void Alloc(int size);
	void CheckSameSize(const DRValarray&) const;
};

template <class T> inline
void DRValarray<T>::print_short (ostream& s) const
{ 
	if (size() < 3) {
		for (int i = 0; i < size(); i++) s << m_ptr[i] << "\t";
		s << endl;
	}
	else {
		for (int i = 0; i < 2; i++) s << m_ptr[i] << "\t";
		s << "...\t" << m_ptr[size() -1] << endl;
	}
}

template <class T> inline
DRValarray<T> MakeCountedArray (const DRValarray<T>& a) //f Creates a counted array
{ 
	DRValarray<T> ans(a.size() + 1);
	ans[0] = a.size();
	for (int i = 0; i < a.size(); i++)  ans[i+1] = a[i];
	return ans;
}


// Common typedefs

typedef DRValarray<int>		IArray;
typedef DRValarray<double>	DArray;
typedef DRValarray<float>	FArray;
typedef DRValarray<short>	SArray;
typedef DRValarray<long>	LArray;
typedef DRValarray<long>	MbsTDateArray;


// Member function implementation
template <class T> inline
DRValarray<T>::DRValarray () : m_ptr (NULL), m_length(0), m_capacity(0) {}

template <class T> inline
DRValarray<T>::DRValarray (size_t n)
{
	Alloc(n);
}

template <class T> inline
DRValarray<T>::DRValarray (const T* p, size_t n) 
{
	Alloc(n);
	for (int i = 0 ; i < m_length; i++) m_ptr[i] = p[i];
}

template <class T> inline
DRValarray<T>::DRValarray (const T& t, size_t n)
{
	Alloc (n); 
	for (int i = 0 ; i < m_length; i++) m_ptr[i] = t;
}

template <class T> inline
DRValarray<T>::~DRValarray() 
{
	if (m_ptr) delete [] m_ptr;
}


template <class T> inline
DRValarray<T>::DRValarray (const DRValarray& a)
{
	Alloc (a.m_length);
	for (int i = 0 ; i < m_length; i++) m_ptr[i] = a.m_ptr[i];
}


template <class T> inline
DRValarray<T>::operator T*()
{
	return m_ptr;
}

template <class T> inline
DRValarray<T>::operator const T*() const
{
	return m_ptr;
}

template <class T> inline
DRValarray<T>& DRValarray<T>::operator=(const DRValarray& a)
{
	if (m_length != 0) CheckSameSize(a);
	else Alloc(a.m_length);
	
	for (int i = 0; i < m_length; i++) m_ptr[i] = a.m_ptr[i];
	return *this;
}

template <class T> inline
DRValarray<T>& DRValarray<T>::operator=(const T& t)
{
	for (int i = 0; i < m_length; i++) m_ptr[i] = t;
	return *this;
}

template <class T> inline
int DRValarray<T>::size() const
{
	return m_length;
}

template <class T> inline
bool DRValarray<T>::operator==(const DRValarray& a) const
{
	if (m_length != a.m_length) return false;

	for (int i = 0; i < m_length; i++) if (m_ptr[i] != a.m_ptr[i]) return false;

	return true;
}

template <class T> inline
void DRValarray<T>::resize(int size)
{
	if (m_capacity < size) {
		if (m_ptr) delete [] m_ptr;
		Alloc(size);
	}
	else {
		m_length = size;
	}
}

template <class T> inline
DRValarray<T> DRValarray<T>::apply (T fn(T)) const
{
	DRValarray ans (m_length);
	for (int i = 0; i < m_length; i++)
		ans.m_ptr[i] = fn (m_ptr[i]);
	return ans;
}

template <class T> inline
DRValarray<T> DRValarray<T>::apply (T fn(const T&)) const
{
	DRValarray ans (m_length);
	for (int i = 0; i < m_length; i++)
		ans.m_ptr[i] = fn (m_ptr[i]);
	return ans;
}

template <class T> inline
T DRValarray<T>::min() const
{
	if (m_length < 1) 
		throw DRException ("Min must have at least one element");
	
	T ans = m_ptr[0];
	for (int i = 1; i < m_length; i++) 
		if (m_ptr[i] < ans) ans = m_ptr[i];
		
		return ans;
}

template <class T> inline
T DRValarray<T>::max() const
{
	if (m_length < 1) 
		throw DRException ("Max must have at least one element");
	
	T ans = m_ptr[0];
	for (int i = 1; i < m_length; i++) 
		if (m_ptr[i] > ans) ans = m_ptr[i];
		
		return ans;
}

template <class T> inline
T DRValarray<T>::sum() const
{
	T ans = 0;
	
	for (int i = 0; i < m_length; i++) 
		ans += m_ptr[i];
	
	return ans;
}

template <class T> inline
T DRValarray<T>::mean() const
{ 
	return sum() / m_length;
}

template <class T> inline
double DRValarray<T>::variance () const
{
	double x2 = 0;
	for (int i = 0; i < m_length; i++)
		x2 += m_ptr[i] * m_ptr[i];
	
	double m = mean();
	
	return x2 / m_length - m*m;
}

template <class T>
inline void DRValarray<T>::CheckSameSize (const DRValarray& a) const 
{
	if (m_length != a.m_length)
		throw DRException("Arrays must be of the same size!");
}

template <class T> inline
void DRValarray<T>::Alloc(int size)
{
	if (size < 0)
		throw DRException ("Cannot Allocate an Array of negative size");
	
	if (size != 0) {
		m_ptr = new T [size];
		if (!m_ptr) 
			throw DRException ("Out of Memory in TSimpleArray");
	}
	else {
		m_ptr = NULL;
	}
	
	m_length = size;
	m_capacity = size;
}

template <class T> 
inline DRValarray<T>& DRValarray<T>::operator/=(const T& a)
{
	for (int i = 0; i < m_length; i++) m_ptr[i] /= a;
	return *this;
}

template <class T> 
inline DRValarray<T>& DRValarray<T>::operator/=(const DRValarray& a)
{
	CheckSameSize(a);
	for (int i = 0; i<m_length; i++) m_ptr[i] /= a.m_ptr[i];
	return *this;
}

template <class T> 
inline DRValarray<T>& DRValarray<T>::operator*=(const T& a)
{
	for (int i = 0; i < m_length; i++) m_ptr[i] *= a;
	return *this;
}

template <class T> 
inline DRValarray<T>& DRValarray<T>::operator*=(const DRValarray& a)
{
	CheckSameSize(a);
	for (int i = 0; i<m_length; i++) m_ptr[i] *= a.m_ptr[i];
	return *this;
}

template <class T> 
inline DRValarray<T>& DRValarray<T>::operator+=(const DRValarray& a)
{
	CheckSameSize(a);
	for (int i = 0; i<m_length; i++) m_ptr[i] += a.m_ptr[i];
	return *this;
}

template <class T> 
inline DRValarray<T>& DRValarray<T>::operator+=(const T& a)
{
	for (int i = 0; i<m_length; i++) m_ptr[i] += a;
	return *this;
}

template <class T> 
inline DRValarray<T>& DRValarray<T>::operator-=(const DRValarray& a)
{
	CheckSameSize(a);
	for (int i = 0; i<m_length; i++) m_ptr[i] -= a.m_ptr[i];
	return *this;
}

template <class T> 
inline DRValarray<T>& DRValarray<T>::operator-=(const T& a)
{
	for (int i = 0; i<m_length; i++) m_ptr[i] -= a;
	return *this;
}

template <class T> 
inline DRValarray<T> DRValarray<T>::operator*(const DRValarray& a) const
{
	CheckSameSize(a);
	DRValarray ans(a.m_length);
	for (int i=0; i<a.m_length; i++) ans.m_ptr[i] = m_ptr[i] * a.m_ptr[i];
	return ans;
}

template <class T> 
inline DRValarray<T> DRValarray<T>::operator*(const T& a) const
{
	DRValarray ans(m_length);
	for (int i=0; i<m_length; i++) ans.m_ptr[i] = m_ptr[i] * a;
	return ans;
}

template <class T> 
inline DRValarray<T> DRValarray<T>::operator-(const DRValarray& a) const
{
	CheckSameSize(a);
	DRValarray ans(a.m_length);
	for (int i=0; i<a.m_length; i++) ans.m_ptr[i] = m_ptr[i] - a.m_ptr[i];
	return ans;
}

template <class T> 
inline DRValarray<T> DRValarray<T>::operator-(const T& a) const
{
	DRValarray ans(m_length);
	for (int i=0; i<m_length; i++) ans.m_ptr[i] = m_ptr[i] - a;
	return ans;
}

template <class T> 
inline DRValarray<T> DRValarray<T>::operator/(const DRValarray& a) const
{
	CheckSameSize(a);
	DRValarray ans(a.m_length);
	for (int i=0; i<a.m_length; i++) ans.m_ptr[i] = m_ptr[i] / a.m_ptr[i];
	return ans;
}

template <class T> 
inline DRValarray<T> DRValarray<T>::operator/(const T& a) const
{
	DRValarray ans(m_length);
	for (int i=0; i<m_length; i++) ans.m_ptr[i] = m_ptr[i] / a;
	return ans;
}

template <class T> 
inline DRValarray<T> DRValarray<T>::operator+(const DRValarray& a) const
{
	CheckSameSize(a);
	DRValarray ans(a.m_length);
	for (int i=0; i<a.m_length; i++) ans.m_ptr[i] = m_ptr[i] + a.m_ptr[i];
	return ans;
}

template <class T> 
inline DRValarray<T> DRValarray<T>::operator+(const T& a) const
{
	DRValarray ans(m_length);
	for (int i=0; i<m_length; i++) ans.m_ptr[i] = m_ptr[i] + a;
	return ans;
}


#endif // !defined(AFX_DRVALARRAY_H__14A1C2E1_226B_11D2_97B2_00C04FD8EB9A__INCLUDED_)
