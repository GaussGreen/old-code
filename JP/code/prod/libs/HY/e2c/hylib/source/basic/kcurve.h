#ifndef __K_DATE_CURVE_H__
#define __K_DATE_CURVE_H__
#include "kplatdep.h"
#include <map>
#include <vector>

template<class Key, class T>
class KCurve
{
public:
	typedef	KMap(Key, T)::iterator iterator;
	typedef	KMap(Key, T)::const_iterator const_iterator;
	iterator begin(){ return _m.begin(); }
	const_iterator begin() const { return _m.begin(); }
	iterator end(){ return _m.end(); }
	const_iterator end() const { return _m.end(); }
	size_t	size()const{return _m.size();}
public:
	KCurve(){}

	KCurve(const KVector(Key) &k, const KVector(T) &t)
	  {
	    if(k.size() != t.size())
	      throw KException("Constructor failed!");
	    KVector(Key)::const_iterator iter1;
	    KVector(T)::const_iterator iter2;
	    for(iter1 = k.begin(), iter2 = t.begin(); iter1 != k.end(); iter1++, iter2++)
	      insert(*iter1, *iter2);
	  }

	~KCurve(){}
	//curve operations
	//assume that T type has operations of +, +=
	KCurve	operator + ( const KCurve &c)const
				{KCurve ans(*this); ans += c; return ans;}
	//merge two curves
	void	operator +=( const KCurve &c);
	
	std::pair<T, bool>	get_value(const Key&d)const;
	//if d is in the curve return &f(d),else insert d and return &f(d)
	T	&operator[](const Key &d){return _m[d];}
	//f(d) += t;
	void	insert(const Key&d, const T&t = T());
	void	erase(const	Key &d){_m.erase(d);}
	void	erase(iterator first, iterator last){_m.erase(first, last);}
	//searching handles
	iterator find(const Key& d){return _m.find(d);}
	const_iterator find(const Key& d) const{return _m.find(d);}
    iterator lower_bound(const Key& d){return _m.lower_bound(d);}
    const_iterator lower_bound(const Key& d) const{return _m.lower_bound(d);}
    iterator upper_bound(const Key& d){return _m.upper_bound(d);}
    const_iterator upper_bound(const Key& d) const{return _m.upper_bound(d);}
	friend	std::ostream & operator<<(std::ostream &, const KCurve&);
private:
	KMap(Key, T)	_m;
} ;
	
template<class Key, class T>	
void	KCurve<Key, T>::operator +=( const KCurve &c)
{
	KMap(Key, T)::const_iterator	iter;
	KMap(Key, T)::iterator iter1;
	for(iter = c._m.begin(); iter != c._m.end(); iter++)
	{
		iter1 = _m.find(iter->first);
		if(iter1 == _m.end())
			_m.insert(KMap(Key, T)::value_type(iter->first, iter->second));
		else
			iter1->second += iter->second; 
	}
}

template <class Key, class T>	
std::pair<T, bool>	KCurve<Key, T>::get_value(const Key&d)const
{
	std::pair<T, bool> ans;
	ans.second = false;
	KMap(Key, T)::iterator	iter;
	iter  = _m.find(d);
	if(iter != _m.end())
	{
		ans.first = iter->second;
		ans.second = true;
	}
	return ans;
}
template <class Key, class T>	
void	KCurve<Key, T>::insert(const Key&d, const T&t)
{
	KMap(Key, T)::iterator	iter;
	iter  = _m.find(d);
	if(iter == _m.end())
		_m.insert(KMap(Key, T)::value_type(d, t));
	else
		iter->second += t; 
}
template<class Key, class T>
std::ostream & operator<<(std::ostream &out, const KCurve<Key, T>&c)
{
	KCurve<Key, T>::const_iterator	iter;
	for(iter = c.begin(); iter != c.end(); iter++)
		out<<iter->first<<"\t"<<iter->second<<std::endl;
	return out;
}
#endif


