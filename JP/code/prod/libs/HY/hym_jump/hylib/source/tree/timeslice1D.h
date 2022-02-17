/*********************************************************
*
*	timeslice1D.h
*
**********************************************************/

#ifndef _TIME_SLICE_1D_H_
#define	_TIME_SLICE_1D_H_

#include "kplatdep.h"

template <class T>
class TimeSlice1D	
{
	int			_startIdx;
	KVector(T)	_data;
public:
	TimeSlice1D(){_startIdx = 0;}
	TimeSlice1D(int lim, const T & t = T()):_startIdx(-lim), _data(2*lim+1, t){}
	~TimeSlice1D(void){}
	void	resize(int lim, const T & t = T())
	{_data.clear(); _startIdx = -lim; _data.resize(2*lim+1, t);}
	T	&operator[](int i){return _data[i-_startIdx];}
	T	operator[](int i)const{return _data[i-_startIdx];}
	//equal to a constant, assuming that TimeSlice has been initialized
	void	equal(const T&t, int lim);
	void	apply (const TimeSlice1D& m, KBinaryFunction<T, T, T> &f, int lim);//*this = f(*this, m)
	void	apply (const T& t, KBinaryFunction<T, T, T> &f, int lim); //*this = f(*this, t)
	void	apply (KBinaryFunction<T, T, T> &f, const T& t, int lim); //*this = f(t, *this)

	int size(){return _data.size();}
	void	print(std::ostream &, int lim)const;
};

template <class T> inline
void	TimeSlice1D<T>::equal(const T&t, int lim)
{
	for (int i = -lim; i <= lim; i++) 
		(*this)[i] = t;
}

template <class T> inline
void TimeSlice1D<T>::apply (const TimeSlice1D<T>& m, KBinaryFunction<T, T, T> &f, int lim)
{
	for (int i = -lim; i <= lim; i++)
	{
		T	&tmp = (*this)[i];
		tmp = f(tmp, m[i]);
	}
}

template <class T> inline
void TimeSlice1D<T>::apply (const T& t, KBinaryFunction<T, T, T> &f, int lim)
{
	for (int i = -lim; i <= lim; i++)
	{
		T	&tmp = (*this)[i];
		tmp = f(tmp, t);
	}
}

template <class T> inline
void TimeSlice1D<T>::apply (KBinaryFunction<T, T, T> &f, const T& t, int lim)
{
	for (int i = -lim; i <= lim; i++)
	{
		T	&tmp = (*this)[i];
		tmp = f(t, tmp);
	}
}

template <class T> inline
void	TimeSlice1D<T>::print(std::ostream &out, int lim)const
{
	for (int i = -lim; i <= lim; i++)
		out <<"["<<i<<", "<<(*this)[i]<<"]\t";
}

#endif

