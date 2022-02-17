//----------------------------------------------------------------------------
//
//   Group       : Cross Asset QR
//
//   Filename    : CondensedVector.hpp
//
//   Description : Utility tool for SRM and other piecewise-constant functions
//                 accessed consecuetively
//
//   Date        : 
//
//----------------------------------------------------------------------------

#ifndef QLIB_CONDENSED_VECTOR_HPP
#define QLIB_CONDENSED_VECTOR_HPP

#include "edginc/DECLARE.hpp"
#include <vector>
#include <algorithm>
//#include <cassert>

DRLIB_BEGIN_NAMESPACE

template <typename T, typename IDX = unsigned short>
class CondensedVector {
public:
    CondensedVector(const std::vector<T>& v)
    {
        if (v.empty())
            return;
        data.push_back(v[0]);
        ends.push_back(1);
        for(size_t i=1; i<v.size(); ++i)
            if(
               max(v[i], data.back())-
               min(v[i], data.back()) < 1e-6
              )
                ++ends.back();
            else
            {
                data.push_back(v[i]);
                ends.push_back(ends.back()+1);
            }
    }

    // purposefully no operator[] provided, as this is a possible 
    // yet inefficient way of accessing stored data
    const T& at(size_t pos) const
    {
        return data[findpos(pos)];
    }


    operator std::vector<T> () const 
    {
        std::vector<T> v;
        if (data.empty())
            return v;
        v.reserve(ends.back());
        for(const_iterator i=begin(); i!=end(); ++i)
            v.push_back(*i);
        return v;

    }

	class const_iterator;
	friend class const_iterator;

	class const_iterator
		{	// iterator for nonmutable vector
	public:

        const_iterator() : myvec(0), condensed_pos(0), real_pos(0) {}	// construct with null pointer

        const_iterator(const CondensedVector<T,IDX> *anothervec, IDX cpos = 0, IDX rpos = 0) : 
        myvec(anothervec), condensed_pos(cpos), real_pos(rpos) {} // construct with pointer _Ptr

		const T& operator*() const
        {	// return designated object
		    return (myvec->data[condensed_pos]);
		}

		const T* operator->() const
		{	// return pointer to class object
			return (&**this);
		}

		const_iterator& operator++()
		{	// preincrement
            ++real_pos;
            if (real_pos >= myvec->ends[condensed_pos])
                ++condensed_pos;
			return (*this);
		}

		const_iterator operator++(int)
		{	// postincrement
			const_iterator _Tmp = *this;
			++*this;
			return (_Tmp);
		}

		const_iterator& operator--()
		{	// predecrement
            --real_pos;
            if (condensed_pos && real_pos < myvec->ends[condensed_pos-1])
                --condensed_pos;
			return (*this);
		}

		const_iterator operator--(int)
		{	// postdecrement
			const_iterator _Tmp = *this;
			--*this;
			return (_Tmp);
		}

		int operator-(const const_iterator& _Right) const
		{	// return difference of iterators
			return (real_pos - _Right.real_pos);
		}

		bool operator==(const const_iterator& _Right) const
		{	// test for iterator equality
			return (myvec == _Right.myvec && real_pos == _Right.real_pos);
		}

		bool operator!=(const const_iterator& _Right) const
		{	// test for iterator inequality
			return (!(*this == _Right));
		}

		bool operator<(const const_iterator& _Right) const
		{	// test if this < _Right
			return (real_pos < _Right.real_pos);
		}

		bool operator>(const const_iterator& _Right) const
		{	// test if this > _Right
			return (_Right < *this);
		}

		bool operator<=(const const_iterator& _Right) const
		{	// test if this <= _Right
			return (!(_Right < *this));
		}

		bool operator>=(const const_iterator& _Right) const
		{	// test if this >= _Right
			return (!(*this < _Right));
		}

		const CondensedVector<T,IDX>* myvec;	// condensed vector
        IDX condensed_pos;
        IDX real_pos; 

	};

    const_iterator begin() const
    {
        return const_iterator(this);
    }
    const_iterator end() const
    {
        if (empty())
            return begin();

        return const_iterator(this, (IDX) ends.size(), (IDX) ends.back());
    }
    const_iterator position(size_t pos) const
    {
        if (pos<size())
            return begin();

        return const_iterator(this, findpos(pos), (IDX) pos);
    }

    size_t  size() const
    {
        if (empty())
            return 0;
        return ends.back();
    }
    bool empty() const 
    {
        return data.empty();
    }

    size_t memsize() const 
        // approximate measurement as structure 
        // fields can be aligned
    {
        return data.size()*(sizeof(T)+sizeof(IDX))+sizeof(CondensedVector);
    }
    size_t memsizeorig() const 
        // approximate measurement as structure 
        // fields can be aligned -- for comparison to see if memsize < memsizeorig
    {
        return size()*sizeof(T)+sizeof(vector<T>);
    }

    size_t condsize() const
    {
        return data.size();
    }

    int memsavings() const { return memsizeorig() - memsize();}
private:
    std::vector<T>   data;
    std::vector<IDX> ends;

    IDX findpos(size_t pos) const
    {
        return (IDX) (std::upper_bound(ends.begin(), ends.end(), pos)-ends.begin());
    }
};

template <typename T>
        CondensedVector<typename T::value_type> makeCondensedVector(const T& t)
{
    return CondensedVector<typename T::value_type>(t);
}

#ifndef STUDY
#define STUDY(x) { cerr << "Condensing "<<#x<<" saves " << makeCondensedVector(x).memsavings() << endl;}
#endif

DRLIB_END_NAMESPACE
#endif


