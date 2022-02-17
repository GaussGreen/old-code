/*********************************************************
*
*	timeslice2D.h
*
**********************************************************/

#ifndef _TIME_SLICE_2D_H_
#define	_TIME_SLICE_2D_H_

#include "kplatdep.h"

template <class T>
class TimeSlice2D	
{
	int			_startIdx0;
	int			_startIdx1;
	KVector(T)	_data;
public:
	TimeSlice2D(){_startIdx0 = _startIdx1 = 0;}
	TimeSlice2D(int lim1, int lim2, const T&t=T());
	~TimeSlice2D(void){}
	void	resize(int lim1, int lim2, const T&t=T());
	T	&operator()(int i, int j );
	T	operator()(int i, int j )const;
	//equal to a constant, assuming that TimeSlice has been initialized
	void	equal(const T&t, int lim1, int lim2);
	void	apply (const TimeSlice2D& m, KBinaryFunction<T, T, T> &f, int lim1, int lim2);	//*this = f(*this, m)
	void	apply (const T& t, KBinaryFunction<T, T, T> &f, int lim1, int lim2); //*this = f(*this, t)
	void	apply (KBinaryFunction<T, T, T> &f, const T& t, int lim1, int lim2); //*this = f(t, *this)
	void	print(std::ostream &, int lim1, int lim2)const;
};


template <class T> inline
TimeSlice2D<T>::TimeSlice2D<T>(int lim1, int lim2, const T&t)
{
	_startIdx0 = -lim1;
	_startIdx1 = -lim2;
	_data.clear();
  	_data.resize((2 * lim1 + 1)*(2 * lim2 + 1), t);
}
template <class T> inline
void	TimeSlice2D<T>::resize(int lim1, int lim2, const T&t)
{
	_startIdx0 = -lim1;
	_startIdx1 = -lim2;
	_data.clear();
  	_data.resize((2 * lim1 + 1)*(2 * lim2 + 1), t);
}

template <class T>	inline
T	&TimeSlice2D<T>::operator()(int i, int j )
{
	return _data[(i - _startIdx0) * ( -2 * _startIdx1 + 1)  +  j - _startIdx1];
}
template <class T>	inline
T	TimeSlice2D<T>::operator()(int i, int j )const
{
	return _data[(i - _startIdx0) * ( -2 * _startIdx1 + 1)  +  j - _startIdx1];
}

template <class T> inline
void	TimeSlice2D<T>::equal(const T&t, int lim1, int lim2)
{
	int	i, j;
	for (i = -lim1; i <= lim1; i++) 
		for (j = -lim1; j <= lim1; j++) 
			(*this)(i, j) = t;
}

template <class T> inline
void TimeSlice2D<T>::apply (const TimeSlice2D<T>& m, KBinaryFunction<T, T, T> &f, int lim1, int lim2)
{
	int	i, j;
	for (i = -lim1; i <= lim1; i++)
		for (j = -lim2; j <= lim2; j++)
		{
			T	&tmp = (*this)(i, j);
			tmp = f(tmp, m(i, j));
		}
}
template <class T> inline
void TimeSlice2D<T>::apply (const T& t, KBinaryFunction<T, T, T> &f, int lim1, int lim2)
{
	int	i, j;
	for (i = -lim1; i <= lim1; i++)
		for (j = -lim2; j <= lim2; j++)
		{
			T	&tmp = (*this)(i, j);
			tmp = f(tmp, t);
		}
}
template <class T> inline
void TimeSlice2D<T>::apply (KBinaryFunction<T, T, T> &f, const T& t, int lim1, int lim2)
{
	int	i, j;
	for (i = -lim1; i <= lim1; i++)
		for (j = -lim2; j <= lim2; j++)
		{
			T	&tmp = (*this)(i, j);
			tmp = f(t, tmp);
		}
}

template <class T> inline
void	TimeSlice2D<T>::print(std::ostream &out, int lim1, int lim2)const
{
	int	i, j;
	for (i = -lim1; i <= lim1; i++)
	{
		for (j = -lim2; j <= lim2; j++)
			out <<(*this)(i, j)<<"\t";
		out<<std::endl;
	}
}
#endif

