// KValarray.h: interface for the KValarray class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_KVALARRAY_H__14A1C2E1_226B_11D2_97B2_00C04FD8EB9A__INCLUDED_)
#define AFX_KVALARRAY_H__14A1C2E1_226B_11D2_97B2_00C04FD8EB9A__INCLUDED_

/**@#-*/
#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000
/**@#+*/

#include "kexception.h"
#include IOSTREAM_H

//k  Emulates STL valarray without the slices.

template <class T>
class KValarray  
{
public:
	KValarray (); //
	explicit KValarray (size_t n); //
	KValarray (const T&, size_t); //
	KValarray (const T*, size_t); //
	virtual ~KValarray(); //
	
	KValarray (const KValarray&);		// copy constructor
	
	operator T*();	//
	operator const T*() const; //
	
	KValarray& operator=(const KValarray&);	// = operator 
	KValarray& operator=(const T&); // same as fill
	
	int size() const; //
	void resize(int size); //
	
	KValarray apply (T fn(T)) const; //
	KValarray apply (T fn(const T&)) const; //
	
	T min() const;		// Gets smallest element
	T max() const;		// Gets biggest element
	T sum() const;		// Computes the sum of the array
	T mean() const;		// Computes the mean
	double variance () const; // Computes the variance and returns a double
	
	friend T dot_product (const KValarray<T>& x, const KValarray<T>& y) // Computes dot product
	{
		x.CheckSameSize(y); 
		T sum = 0;
		
		for (int i = 0; i < x.m_length; i++) 
			sum += x.m_ptr[i] * y.m_ptr[i];
		
		return sum;
	}
	
	friend double covariance (const KValarray<T>&, const KValarray<T>&) //
	{throw KException ("Never been used"); return 0;}
	
	friend double correlation (const KValarray<T>&, const KValarray<T>&) //
	{throw KException ("Never been used"); return 0;}
	
	friend std::ostream &operator<<(std::ostream& s, const KValarray &a) //
	{
		s << "Length(" << a.m_length << "): " <<std::endl;
		
		for (int i = 0; i< a.m_length; i++) 
			s << a.m_ptr[i] << " ";
		
		s << std::endl; 
		return s;
	}
	
	// the math operations +,-,/,*,+=,-=,/=,*= act only on arrays of the same size
	KValarray operator/(const KValarray&) const; //
	KValarray operator/(const T&) const; //
	KValarray& operator/=(const KValarray&); //
	KValarray& operator/=(const T&); //
	
	KValarray operator*(const KValarray&) const; //
	KValarray operator*(const T&) const; //
	KValarray& operator*=(const KValarray&); //
	KValarray& operator*=(const T&); //
	
	KValarray operator+(const KValarray&) const; //
	KValarray operator+(const T&) const; //
	KValarray& operator+=(const KValarray&); //
	KValarray& operator+=(const T&); //
	
	KValarray operator-(const KValarray&) const; //
	KValarray operator-(const T&) const; //
	KValarray& operator-=(const KValarray&); //
	KValarray& operator-=(const T&); //


	bool operator==(const KValarray&) const;

	void print_short (std::ostream&) const; //
	friend	KValarray	operator *(const T &t, const KValarray &v){return v*t;}
protected:
	int m_length, m_capacity;
	T* m_ptr;
	void Alloc(int size);
	void CheckSameSize(const KValarray&) const;
};

template <class T> inline
void KValarray<T>::print_short (std::ostream& s) const
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
KValarray<T> MakeCountedArray (const KValarray<T>& a) //f Creates a counted array
{ 
	KValarray<T> ans(a.size() + 1);
	ans[0] = a.size();
	for (int i = 0; i < a.size(); i++)  ans[i+1] = a[i];
	return ans;
}


// Common typedefs

typedef KValarray<int>		IArray;
typedef KValarray<double>	DArray;

// Member function implementation
template <class T> inline
KValarray<T>::KValarray () : m_ptr (NULL), m_length(0), m_capacity(0) {}

template <class T> inline
KValarray<T>::KValarray (size_t n)
{
	Alloc(n);
}

template <class T> inline
KValarray<T>::KValarray (const T* p, size_t n) 
{
	Alloc(n);
	for (int i = 0 ; i < m_length; i++) m_ptr[i] = p[i];
}

template <class T> inline
KValarray<T>::KValarray (const T& t, size_t n)
{
	Alloc (n); 
	for (int i = 0 ; i < m_length; i++) m_ptr[i] = t;
}

template <class T> inline
KValarray<T>::~KValarray() 
{
	if (m_ptr) delete [] m_ptr;
}


template <class T> inline
KValarray<T>::KValarray (const KValarray& a)
{
	Alloc (a.m_length);
	for (int i = 0 ; i < m_length; i++) m_ptr[i] = a.m_ptr[i];
}


template <class T> inline
KValarray<T>::operator T*()
{
	return m_ptr;
}

template <class T> inline
KValarray<T>::operator const T*() const
{
	return m_ptr;
}

template <class T> inline
KValarray<T>& KValarray<T>::operator=(const KValarray& a)
{
	if (m_length != 0) CheckSameSize(a);
	else Alloc(a.m_length);
	
	for (int i = 0; i < m_length; i++) m_ptr[i] = a.m_ptr[i];
	return *this;
}

template <class T> inline
KValarray<T>& KValarray<T>::operator=(const T& t)
{
	for (int i = 0; i < m_length; i++) m_ptr[i] = t;
	return *this;
}

template <class T> inline
int KValarray<T>::size() const
{
	return m_length;
}

template <class T> inline
bool KValarray<T>::operator==(const KValarray& a) const
{
	if (m_length != a.m_length) return false;

	for (int i = 0; i < m_length; i++) if (m_ptr[i] != a.m_ptr[i]) return false;

	return true;
}

template <class T> inline
void KValarray<T>::resize(int size)
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
KValarray<T> KValarray<T>::apply (T fn(T)) const
{
	KValarray ans (m_length);
	for (int i = 0; i < m_length; i++)
		ans.m_ptr[i] = fn (m_ptr[i]);
	return ans;
}

template <class T> inline
KValarray<T> KValarray<T>::apply (T fn(const T&)) const
{
	KValarray ans (m_length);
	for (int i = 0; i < m_length; i++)
		ans.m_ptr[i] = fn (m_ptr[i]);
	return ans;
}

template <class T> inline
T KValarray<T>::min() const
{
	if (m_length < 1) 
		throw KException ("Min must have at least one element");
	
	T ans = m_ptr[0];
	for (int i = 1; i < m_length; i++) 
		if (m_ptr[i] < ans) ans = m_ptr[i];
		
		return ans;
}

template <class T> inline
T KValarray<T>::max() const
{
	if (m_length < 1) 
		throw KException ("Max must have at least one element");
	
	T ans = m_ptr[0];
	for (int i = 1; i < m_length; i++) 
		if (m_ptr[i] > ans) ans = m_ptr[i];
		
		return ans;
}

template <class T> inline
T KValarray<T>::sum() const
{
	T ans = 0;
	
	for (int i = 0; i < m_length; i++) 
		ans += m_ptr[i];
	
	return ans;
}

template <class T> inline
T KValarray<T>::mean() const
{ 
	return sum() / m_length;
}

template <class T> inline
double KValarray<T>::variance () const
{
	double x2 = 0;
	for (int i = 0; i < m_length; i++)
		x2 += m_ptr[i] * m_ptr[i];
	
	double m = mean();
	
	return x2 / m_length - m*m;
}

template <class T>
inline void KValarray<T>::CheckSameSize (const KValarray& a) const 
{
	if (m_length != a.m_length)
		throw KException("Arrays must be of the same size!");
}

template <class T> inline
void KValarray<T>::Alloc(int size)
{
	if (size < 0)
		throw KException ("Cannot Allocate an Array of negative size");
	
	if (size != 0) {
		m_ptr = new T [size];
		if (!m_ptr) 
			throw KException ("Out of Memory in TSimpleArray");
	}
	else {
		m_ptr = NULL;
	}
	
	m_length = size;
	m_capacity = size;
}

template <class T> 
inline KValarray<T>& KValarray<T>::operator/=(const T& a)
{
	for (int i = 0; i < m_length; i++) m_ptr[i] /= a;
	return *this;
}

template <class T> 
inline KValarray<T>& KValarray<T>::operator/=(const KValarray& a)
{
	CheckSameSize(a);
	for (int i = 0; i<m_length; i++) m_ptr[i] /= a.m_ptr[i];
	return *this;
}

template <class T> 
inline KValarray<T>& KValarray<T>::operator*=(const T& a)
{
	for (int i = 0; i < m_length; i++) m_ptr[i] *= a;
	return *this;
}

template <class T> 
inline KValarray<T>& KValarray<T>::operator*=(const KValarray& a)
{
	CheckSameSize(a);
	for (int i = 0; i<m_length; i++) m_ptr[i] *= a.m_ptr[i];
	return *this;
}

template <class T> 
inline KValarray<T>& KValarray<T>::operator+=(const KValarray& a)
{
	CheckSameSize(a);
	for (int i = 0; i<m_length; i++) m_ptr[i] += a.m_ptr[i];
	return *this;
}

template <class T> 
inline KValarray<T>& KValarray<T>::operator+=(const T& a)
{
	for (int i = 0; i<m_length; i++) m_ptr[i] += a;
	return *this;
}

template <class T> 
inline KValarray<T>& KValarray<T>::operator-=(const KValarray& a)
{
	CheckSameSize(a);
	for (int i = 0; i<m_length; i++) m_ptr[i] -= a.m_ptr[i];
	return *this;
}

template <class T> 
inline KValarray<T>& KValarray<T>::operator-=(const T& a)
{
	for (int i = 0; i<m_length; i++) m_ptr[i] -= a;
	return *this;
}

template <class T> 
inline KValarray<T> KValarray<T>::operator*(const KValarray& a) const
{
	CheckSameSize(a);
	KValarray ans(a.m_length);
	for (int i=0; i<a.m_length; i++) ans.m_ptr[i] = m_ptr[i] * a.m_ptr[i];
	return ans;
}

template <class T> 
inline KValarray<T> KValarray<T>::operator*(const T& a) const
{
	KValarray ans(m_length);
	for (int i=0; i<m_length; i++) ans.m_ptr[i] = m_ptr[i] * a;
	return ans;
}

template <class T> 
inline KValarray<T> KValarray<T>::operator-(const KValarray& a) const
{
	CheckSameSize(a);
	KValarray ans(a.m_length);
	for (int i=0; i<a.m_length; i++) ans.m_ptr[i] = m_ptr[i] - a.m_ptr[i];
	return ans;
}

template <class T> 
inline KValarray<T> KValarray<T>::operator-(const T& a) const
{
	KValarray ans(m_length);
	for (int i=0; i<m_length; i++) ans.m_ptr[i] = m_ptr[i] - a;
	return ans;
}

template <class T> 
inline KValarray<T> KValarray<T>::operator/(const KValarray& a) const
{
	CheckSameSize(a);
	KValarray ans(a.m_length);
	for (int i=0; i<a.m_length; i++) ans.m_ptr[i] = m_ptr[i] / a.m_ptr[i];
	return ans;
}

template <class T> 
inline KValarray<T> KValarray<T>::operator/(const T& a) const
{
	KValarray ans(m_length);
	for (int i=0; i<m_length; i++) ans.m_ptr[i] = m_ptr[i] / a;
	return ans;
}

template <class T> 
inline KValarray<T> KValarray<T>::operator+(const KValarray& a) const
{
	CheckSameSize(a);
	KValarray ans(a.m_length);
	for (int i=0; i<a.m_length; i++) ans.m_ptr[i] = m_ptr[i] + a.m_ptr[i];
	return ans;
}

template <class T> 
inline KValarray<T> KValarray<T>::operator+(const T& a) const
{
	KValarray ans(m_length);
	for (int i=0; i<m_length; i++) ans.m_ptr[i] = m_ptr[i] + a;
	return ans;
}


#endif // !defined(AFX_KValarray_H__14A1C2E1_226B_11D2_97B2_00C04FD8EB9A__INCLUDED_)
