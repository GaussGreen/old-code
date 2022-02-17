/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file specialsort.h
 *  \brief special sort can sort two vector in once
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date May 2004
 */

#ifndef _INGPBASE_SPECIALSORT_H
#define _INGPBASE_SPECIALSORT_H

#include "port.h"
#include "gpvector.h"

CC_BEGIN_NAMESPACE( ARM )

template <typename T = double, typename U = T>
	struct ARM_T_Sort
{
private:
// FIXMEFRED: mig.vc8 (21/05/2007 10:54:04): typedef typename
#if _MSC_VER >= 1310	// Visual C++ 2005 & .net 2003
	typedef typename CC_NS(ARM,ARM_GP_T_Vector)<T>::iterator t_iterator;
	typedef typename CC_NS(ARM,ARM_GP_T_Vector)<T>::const_iterator t_const_iterator;
	typedef typename CC_NS(ARM,ARM_GP_T_Vector)<U>::iterator u_iterator;
	typedef typename CC_NS(ARM,ARM_GP_T_Vector)<U>::const_iterator u_const_iterator;
#else					// Visual C++ 6
	typedef CC_NS(ARM,ARM_GP_T_Vector)<T>::iterator t_iterator;
	typedef CC_NS(ARM,ARM_GP_T_Vector)<T>::const_iterator t_const_iterator;
	typedef CC_NS(ARM,ARM_GP_T_Vector)<U>::iterator u_iterator;
	typedef CC_NS(ARM,ARM_GP_T_Vector)<U>::const_iterator u_const_iterator;
#endif 

public:
	/// sort global constant
	static const size_t _INSERTION_SORT_SIZE;

	/// static
	static void sortTwoVectorsWithSameSize( ARM_GP_T_Vector<T>& vec1, ARM_GP_T_Vector<U>& vec2 );

	/// part on sort inspired by the STL
	static void _Sort(t_iterator _F, t_iterator _L, ARM_GP_T_Vector<T>& vec1, ARM_GP_T_Vector<U>& vec2);
	static t_iterator _Unguarded_partition(t_iterator _F, t_iterator _L, T _Piv, ARM_GP_T_Vector<T>& vec1, ARM_GP_T_Vector<U>& vec2);
	static void _Insertion_sort(t_iterator _F, t_iterator _L, ARM_GP_T_Vector<T>& vec1, ARM_GP_T_Vector<U>& vec2);
	static void _Unguarded_insert(t_iterator _L, ARM_GP_T_Vector<T>& vec1, ARM_GP_T_Vector<U>& vec2);
	static inline T _Median(T _X, T _Y, T _Z)
	{	if (_X < _Y) 
			return (_Y<_Z ? _Y: _X<_Z ? _Z:_X);
		else  
			return (_X<_Z ? _X: _Y<_Z ? _Z:_Y);
	}
};

template <typename T, typename U> 
	const size_t ARM_T_Sort<T,U>::_INSERTION_SORT_SIZE  = 16;

/// sort

template <typename T, typename U> 
	inline void ARM_T_Sort<T,U>::sortTwoVectorsWithSameSize( ARM_GP_T_Vector<T>& vec1, ARM_GP_T_Vector<U>& vec2 )
{	
	t_iterator begin = vec1.begin(), end = vec1.end();

	if (end - begin <= _INSERTION_SORT_SIZE)
		ARM_T_Sort<T,U>::_Insertion_sort(begin, end, vec1, vec2);
	else
	{
		ARM_T_Sort<T,U>::_Sort(begin, end, vec1, vec2);
		ARM_T_Sort<T,U>::_Insertion_sort(begin, begin + _INSERTION_SORT_SIZE, vec1, vec2);
		for( begin += _INSERTION_SORT_SIZE; begin != end; ++begin )
			ARM_T_Sort<T,U>::_Unguarded_insert(begin, vec1, vec2); 
	}
}

template <typename T, typename U> 
	inline void ARM_T_Sort<T,U>::_Sort(t_iterator _F, t_iterator _L, ARM_GP_T_Vector<T>& vec1, ARM_GP_T_Vector<U>& vec2 )
{
	for (; _INSERTION_SORT_SIZE < _L - _F; )
	{
		t_iterator _M = _Unguarded_partition(_F, _L, _Median(*_F,*(_F+(_L-_F)/2),*(_L-1)), vec1, vec2 );
		if (_L - _M <= _M - _F)
			ARM_T_Sort<T,U>::_Sort(_M, _L, vec1, vec2 ), _L = _M;
		else
			ARM_T_Sort<T,U>::_Sort(_F, _M, vec1, vec2 ), _F = _M; 
	}
}

template <typename T, typename U> 
// FIXMEFRED: mig.vc8 (21/05/2007 11:10:35): t_iterator instead of the described type
inline typename ARM_T_Sort<T, U>::t_iterator ARM_T_Sort<T,U>::_Unguarded_partition(t_iterator _F, t_iterator _L, T _Piv, ARM_GP_T_Vector<T>& vec1, ARM_GP_T_Vector<U>& vec2)
{
	for (; ; ++_F)
	{
		for (; *_F < _Piv; ++_F) ;
		for (; _Piv < *--_L; ) ;
		if (_L <= _F)
			return (_F);
		CC_NS(std,iter_swap)(vec2.begin()+(_F-vec1.begin()),vec2.begin()+(_L-vec1.begin()) );
		CC_NS(std,iter_swap)(_F, _L);
	}
}

template <typename T, typename U> 
	inline void ARM_T_Sort<T,U>::_Insertion_sort(t_iterator _F, t_iterator _L, ARM_GP_T_Vector<T>& vec1, ARM_GP_T_Vector<U>& vec2)
{
	if (_F != _L)
		for (t_iterator _M = _F; ++_M != _L; )
		{
			T _V = *_M;
			if (!(_V < *_F))
				ARM_T_Sort<T,U>::_Unguarded_insert(_M,vec1,vec2);
			else
			{
				u_iterator _FII = vec2.begin()+(_F-vec1.begin()),
					_MII = vec2.begin()+(_M-vec1.begin());
				U _VII = *_MII;
				CC_NS(std,copy_backward)(_F, _M, _M+1);
				*_F = _V; 
				CC_NS(std,copy_backward)(_FII, _MII, _MII+1);
				*_FII = _VII; 
			}
		}
}

template <typename T, typename U>
	inline void ARM_T_Sort<T,U>::_Unguarded_insert(t_iterator _L, ARM_GP_T_Vector<T>& vec1, ARM_GP_T_Vector<U>& vec2)
{
	u_iterator _LII = vec2.begin()+(_L-vec1.begin());
	T _V = *_L;
	U _VII = *_LII;
	t_iterator _M = _L;
	u_iterator _MII = _LII;
	for(; --_MII, _V<*--_M; _L=_M, _LII=_MII)
	{
		*_L = *_M;
		*_LII = *_MII; 
	}
	*_L	= _V;
	*_LII = _VII;
}

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
