/*********************************************************
*
*	timeslice.h
*
**********************************************************/

#ifndef _TIME_SLICE_H
#define	_TIME_SLICE_H
#include "kfunction.h"
#include "timeslice1D.h"
#include "timeslice2D.h"

class	BaseTree;
template<class T>
class TimeSlice  
{
public:
	TimeSlice(void){_ptree = NULL; t_type = EMPTY;}
	TimeSlice(const BaseTree *, const T& t=T());
	TimeSlice(const TimeSlice&);
	void	resize(const BaseTree *, const T& t=T());
	int	 dim()const{return t_type;} // 

	//for efficiency, no type checking
	T	&operator[](int i){return t_1d[i];}
	T	operator[](int i)const{return t_1d[i];}
	//for efficiency, no type checking
	T	&operator()(int i, int j){return t_2d(i, j);}
	T	operator()(int i, int j)const{return t_2d(i, j);}

	int size(){
		if(t_type == ONE_DIM)
		{
			return t_1d.size();
		}
		else
		{
			return 0;
		}
	}


	TimeSlice	&operator= (const TimeSlice&);
	TimeSlice	&operator = (const T&t);
	void	operator += (const TimeSlice &);
	void	operator -= (const TimeSlice &);
	void	operator *= (const TimeSlice &);
   	void	operator /= (const TimeSlice &);

	void	operator += (const T &t);
	void	operator -= (const T &t);
	void	operator *= (const T &t);
   	void	operator /= (const T &t);

	TimeSlice	operator + (const TimeSlice &m)const
	{	TimeSlice ans(*this); ans += m; return ans;}
	TimeSlice	operator - (const TimeSlice &m)const
	{	TimeSlice ans(*this); ans -= m; return ans;}
	TimeSlice	operator * (const TimeSlice &m)const
	{	TimeSlice ans(*this); ans *= m; return ans;}
	TimeSlice	operator / (const TimeSlice &m)const
	{	TimeSlice ans(*this); ans /= m; return ans;}

	TimeSlice	operator + (const T &t)const
	{	TimeSlice ans(*this); ans += t; return ans;}
	TimeSlice	operator - (const T &t)const
	{	TimeSlice ans(*this); ans -= t; return ans;}
	TimeSlice	operator * (const T &t)const
	{	TimeSlice ans(*this); ans *= t; return ans;}
   	TimeSlice	operator / (const T &t)const
	{	TimeSlice ans(*this); ans /= t; return ans;}

	friend	TimeSlice	operator + (const T &t, const TimeSlice &m){return m + t;}
	friend	TimeSlice	operator - (const T &t, const TimeSlice &m);
	friend	TimeSlice	operator * (const T &t, const TimeSlice &m){return m * t;}
   	friend	TimeSlice	operator / (const T &t, const TimeSlice &m);
	
	friend	std::ostream	& operator<<(std::ostream &, const TimeSlice&);
protected:
	const BaseTree	*_ptree; 
	enum TimeSliceType {EMPTY = 0, ONE_DIM, TWO_DIM};

	TimeSliceType	t_type;
	TimeSlice1D<T>	t_1d;
	TimeSlice2D<T>	t_2d;
};


typedef	TimeSlice<double>	DTimeSlice;
typedef	TimeSlice<KDoubleFunction>	FTimeSlice;

template	<class T> inline
TimeSlice<T>& TimeSlice<T>::operator= (const TimeSlice<T>& m)
{
	if(this != &m)
	{
		if (_ptree == NULL) 
			resize(m._ptree);
		if(m.t_type == ONE_DIM)
		{
			t_type = m.t_type; 
			KValarray<int>	&lim = _ptree->get_curr_limit();
			for(int i=-lim[0]; i<= lim[0]; i++)
				t_1d[i] = m.t_1d[i];
		}
		else if(m.t_type == TWO_DIM)
		{
			t_type = m.t_type; 
			KValarray<int>	&lim = _ptree->get_curr_limit();
			for(int i=-lim[0]; i<= lim[0]; i++)
				for(int j= -lim[1]; j<= lim[1]; j++)
					t_2d(i, j) = m.t_2d(i, j);
		}
		else
			throw KException("TimeSlice is not defined!");
	}
	return *this;
}
template	<class T> inline
TimeSlice<T>::TimeSlice<T>(const TimeSlice<T>& m)
{
	if (m._ptree != NULL)
	{
		resize(m._ptree);
		if(t_type == ONE_DIM)
		{
			KValarray<int>	&lim = _ptree->get_curr_limit();
			for(int i=-lim[0]; i<= lim[0]; i++)
				t_1d[i] = m.t_1d[i];
		}
		else if(t_type == TWO_DIM)
		{
			KValarray<int>	&lim = _ptree->get_curr_limit();
			for(int i=-lim[0]; i<= lim[0]; i++)
				for(int j= -lim[1]; j<= lim[1]; j++)
					t_2d(i, j) = m.t_2d(i, j);
		}
		else
			throw KException("TimeSlice is not defined!");
	}
}

template <class T> inline
TimeSlice<T>::TimeSlice<T>( const BaseTree *tree, const T&t)
{
	_ptree = tree;
	KValarray<int>	&lim = _ptree->get_max_limit();
	if(lim.size() == 1)
	{
		t_type = ONE_DIM;
		t_1d.resize(lim[0], t);
	}
	else if(lim.size() == 2)
	{
		t_type = TWO_DIM;
		t_2d.resize(lim[0], lim[1], t);
	}
	else
		throw KException("Higher dimensional TimeSlice is not yet implemented!"); 
}

	

template <class T> inline
void	TimeSlice<T>::resize(const BaseTree *tree, const T&t)
{
	if(tree == NULL)
		throw KException("Null pointer!"); 
	_ptree = tree;
	KValarray<int>	&lim = _ptree->get_max_limit();
	if(lim.size() == 1)
	{
		t_type = ONE_DIM;
		t_1d.resize(lim[0], t);
	}
	else if(lim.size() == 2)
	{
		t_type = TWO_DIM;
		t_2d.resize(lim[0], lim[1], t);
	}
	else
		throw KException("Higher dimensional TimeSlice is not yet implemented!"); 
}

template <class T>
std::ostream	& operator<<(std::ostream &out, const TimeSlice<T>&ts)
{
	if(ts.dim() == 1)
	{
		KValarray<int>	&lim = ts._ptree->get_curr_limit();
		ts.t_1d.print(out, lim[0]);
	}
	else if(ts.dim() == 2)
	{
		KValarray<int>	&lim = ts._ptree->get_curr_limit();
		ts.t_2d.print(out, lim[0], lim[1]);
	}
	else
		throw KException("TimeSlice is not defined!");
	return out;
}
template <class T> inline
TimeSlice<T>	&TimeSlice<T>::operator = (const T &t)
{
	if(t_type == ONE_DIM)
	{
		KValarray<int>	&lim = _ptree->get_curr_limit();
		t_1d.equal(t, lim[0]);
	}
	else if(t_type == TWO_DIM)
	{
		KValarray<int>	&lim = _ptree->get_curr_limit();
		t_2d.equal(t, lim[0], lim[1]);
	}
	else
		throw KException("TimeSlice is not defined!");
	return *this;
}

template <class T> inline
void	TimeSlice<T>::operator += (const TimeSlice<T> &m)
{ 
	if(t_type == ONE_DIM)
	{
		KValarray<int>	&lim = _ptree->get_curr_limit();
		t_1d.apply(m.t_1d, KPlus<T>(), lim[0]);
	}
	else if(t_type == TWO_DIM)
	{
		KValarray<int>	&lim = _ptree->get_curr_limit();
		t_2d.apply(m.t_2d, KPlus<T>(), lim[0], lim[1]);
	}
	else
		throw KException("TimeSlice is not defined!");
}

template <class T> inline
void	TimeSlice<T>::operator -= (const TimeSlice<T> &m)
{ 
	if(t_type == ONE_DIM)
	{
		KValarray<int>	&lim = _ptree->get_curr_limit();
		t_1d.apply(m.t_1d, KMinus<T>(), lim[0]);
	}
	else if(t_type == TWO_DIM)
	{
		KValarray<int>	&lim = _ptree->get_curr_limit();
		t_2d.apply(m.t_2d, KMinus<T>(), lim[0], lim[1]);
	}
	else
		throw KException("TimeSlice is not defined!");
}

template <class T>
void	TimeSlice<T>::operator *= (const TimeSlice<T> &m)
{ 
	if(t_type == ONE_DIM)
	{
		KValarray<int>	&lim = _ptree->get_curr_limit();
		t_1d.apply(m.t_1d,  KMultiplies<T>(), lim[0]);
	}
	else if(t_type == TWO_DIM)
	{
		KValarray<int>	&lim = _ptree->get_curr_limit();
		t_2d.apply(m.t_2d,  KMultiplies<T>(), lim[0], lim[1]);
	}
	else
		throw KException("TimeSlice is not defined!");
}

template <class T>
void	TimeSlice<T>::operator /= (const TimeSlice<T> &m)
{ 
	if(t_type == ONE_DIM)
	{
		KValarray<int>	&lim = _ptree->get_curr_limit();
		t_1d.apply(m.t_1d,  KDivides<T>(), lim[0]);
	}
	else if(t_type == TWO_DIM)
	{
		KValarray<int>	&lim = _ptree->get_curr_limit();
		t_2d.apply(m.t_2d,  KDivides<T>(), lim[0], lim[1]);
	}
	else
		throw KException("TimeSlice is not defined!");
}

template <class T> inline
void	TimeSlice<T>::operator += (const T &t)
{ 
	if(t_type == ONE_DIM)
	{
		KValarray<int>	&lim = _ptree->get_curr_limit();
		t_1d.apply(t, KPlus<T>(), lim[0]);
	}
	else if(t_type == TWO_DIM)
	{
		KValarray<int>	&lim = _ptree->get_curr_limit();
		t_2d.apply(t, KPlus<T>(), lim[0], lim[1]);
	}
	else
		throw KException("TimeSlice is not defined!");
}

template <class T> inline
void	TimeSlice<T>::operator -= (const T &t)
{ 
	if(t_type == ONE_DIM)
	{
		KValarray<int>	&lim = _ptree->get_curr_limit();
		t_1d.apply(t, KMinus<T>(), lim[0]);
	}
	else if(t_type == TWO_DIM)
	{
		KValarray<int>	&lim = _ptree->get_curr_limit();
		t_2d.apply(t, KMinus<T>(), lim[0], lim[1]);
	}
	else
		throw KException("TimeSlice is not defined!");
}

template <class T> inline
void	TimeSlice<T>::operator *= (const T &t)
{ 
	if(t_type == ONE_DIM)
	{
		KValarray<int>	&lim = _ptree->get_curr_limit();
		t_1d.apply(t,  KMultiplies<T>(), lim[0]);
	}
	else if(t_type == TWO_DIM)
	{
		KValarray<int>	&lim = _ptree->get_curr_limit();
		t_2d.apply(t,  KMultiplies<T>(), lim[0], lim[1]);
	}
	else
		throw KException("TimeSlice is not defined!");
}

template <class T> inline
void	TimeSlice<T>::operator /= (const T &t)
{ 
	if(t_type == ONE_DIM)
	{
		KValarray<int>	&lim = _ptree->get_curr_limit();
		t_1d.apply(t,  KDivides<T>(), lim[0]);
	}
	else if(t_type == TWO_DIM)
	{
		KValarray<int>	&lim = _ptree->get_curr_limit();
		t_2d.apply(t,  KDivides<T>(), lim[0], lim[1]);
	}
	else
		throw KException("TimeSlice is not defined!");
}

template <class T> inline
TimeSlice<T>	operator - (const T &t, const TimeSlice<T> &m)
{
	TimeSlice<T> ans(m);
	if(ans.dim() == 1)
	{
		KValarray<int>	&lim = ans._ptree->get_curr_limit();
		ans.t_1d.apply(KMinus<T>(), t, lim[0]);
	}
	else if(ans.dim() == 2)
	{
		KValarray<int>	&lim = ans._ptree->get_curr_limit();
		ans.t_2d.apply(KMinus<T>(), t, lim[0], lim[1]);
	}
	else
		throw KException("TimeSlice is not defined!");
	return ans;
}
template <class T>	inline
TimeSlice<T>	operator / (const T &t, const TimeSlice<T> &m)
{
	TimeSlice<T> ans(m);
	if(ans.dim() == 1)
	{
		KValarray<int>	&lim = ans._ptree->get_curr_limit();
		ans.t_1d.apply(KDivides<T>(), t, lim[0]);
	}
	else if(ans.dim() == 2)
	{
		KValarray<int>	&lim = ans._ptree->get_curr_limit();
		ans.t_2d.apply(KDivides<T>(), t, lim[0], lim[1]);
	}
	else
		throw KException("TimeSlice is not defined!");
	return ans;
}



#endif



