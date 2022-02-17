#include <algorithm>
#include <iterator>
#include "edginc/ModelException.hpp"

#ifndef SMALL_LINEAR_INTERPOLATOR_HPP
#define SMALL_LINEAR_INTERPOLATOR_HPP

DRLIB_BEGIN_NAMESPACE


template <typename T, typename F, typename D = const double> 
class func_iterator : public iterator<random_access_iterator_tag, D, long,
			D*, const D>
            //public _Ranit<D, size_t, D*, const D>
{
public:
    typedef random_access_iterator_tag iterator_category;
    typedef D value_type;
    typedef long difference_type;
    typedef size_t size_type;
    typedef D* pointer;
    typedef const D reference;
    typedef const D const_reference;

    func_iterator(T* _p, F _f, size_t _pos = 0) : p(_p), fun(_f), pos(_pos)
    {}    // construct with pointer _Ptr
    
    reference operator*() const
    {    // return designated object
        return (p->*fun)(pos);//calling p->function(pos);
    }

    func_iterator& operator++()
    {    // preincrement
        ++pos;
        return (*this);
    }

    func_iterator operator++(int)
    {    // postincrement
        func_iterator _Tmp = *this;
        ++pos;
        return (_Tmp);
    }

    func_iterator& operator--()
    {    // predecrement
        --pos;
        return (*this);
    }

    func_iterator operator--(int)
    {    // postdecrement
        func_iterator _Tmp = *this;
        --pos;
        return (_Tmp);
    }

    func_iterator& operator+=(difference_type _Off)
    {    // increment by integer
        pos += _Off;
        return (*this);
    }

    func_iterator operator+(difference_type _Off) const
    {    // return this + integer
        func_iterator _Tmp = *this;
        return (_Tmp += _Off);
    }

    func_iterator& operator-=(difference_type _Off)
    {    // decrement by integer
        pos -= _Off;
        return (*this);
    }

    func_iterator operator-(difference_type _Off) const
    {    // return this - integer
        func_iterator _Tmp = *this;
        return (_Tmp -= _Off);
    }

    difference_type operator-(const func_iterator& _Right) const
    {    // return difference of iterators
        return (pos - _Right.pos);
    }

    const_reference operator[](difference_type _Off) const
    {    // subscript
        return (*(*this + _Off));
    }

    bool operator==(const func_iterator& _Right) const
    {    // test for iterator equality


        return (p == _Right.p && pos == _Right.pos);
    }

    bool operator!=(const func_iterator& _Right) const
    {    // test for iterator inequality
        return (!(*this == _Right));
    }

    bool operator<(const func_iterator& _Right) const
    {    // test if this < _Right
        return (pos < _Right.pos);
    }

    bool operator>(const func_iterator& _Right) const
    {    // test if this > _Right
        return (_Right < *this);
    }

    bool operator<=(const func_iterator& _Right) const
    {    // test if this <= _Right
        return (!(_Right < *this));
    }

    bool operator>=(const func_iterator& _Right) const
    {    // test if this >= _Right
        return (!(*this < _Right));
    }

    friend func_iterator operator+(difference_type _Off,
        const func_iterator& _Right)
    {    // return iterator + integer
        return (_Right + _Off);
    }
    func_iterator() : p(NULL), fun(NULL), pos(0) {}

private:

    T*        p;
    F         fun;
    size_type pos;
};

template <typename T, typename F> 
func_iterator<T, F> make_iterator(T* ptr, F fun, size_t offset)
{
    return func_iterator<T, F>(ptr, fun, offset);
}
template <typename T, typename F> 
func_iterator<T, F> make_begin_iterator(T* ptr, F fun)
{
    return func_iterator<T, F>(ptr, fun, 0);
}

template <typename T, typename F> 
func_iterator<T, F> make_end_iterator(T* ptr, F fun, size_t nElems)
{
    return func_iterator<T, F>(ptr, fun, nElems);
}

template <typename RndIterX, typename RndIterY>
double smallLinearInterpolation(
                                RndIterX xBegin, 
                                RndIterX xEnd, 
                                RndIterY yBegin, 
                                RndIterY yEnd, 
                                double xx)
{
    QLIB_VERIFY(xBegin < xEnd, "Interpolation called on empty array.");
    QLIB_VERIFY (xEnd-xBegin <= yEnd-yBegin, "The ranges of x and y have different sizes. ");
    if (xx <= *xBegin)
        return *yBegin;  // no extrapolation
    if (*(xEnd-1) < xx)
        return *(yEnd-1); // no extrapolation, but !!!NOTE!!! - y is allowed to be longer than x
                          // in any case, the last value in y array will be returned.
    RndIterX xptr = lower_bound(xBegin, xEnd, xx);
    double xnext = *xptr;
    double ynext = yBegin[xptr-xBegin];
    double xprev = *--xptr;
    double yprev = yBegin[xptr-xBegin];
    double dxx = xx-xprev;
    double dx = xnext-xprev;
    double dy = ynext-yprev;

    return yprev + (dy/dx)*dxx;

}


/*******

    A tutorial on using this interpolator class

    if you have a class that has two methods inside,

    class MyData {
    public:
        double funcX(int i);
        double funcY(int i);
        size_t nMax();
        // ...
    private:
        // ...
    };
    with allowed range for i to be between [0, nMax-1] inclusive, 
    then to find the interpolated value of y that corresponds to some x

    and an object of this class

    MyData *myObject;

    then call 
    double y = smallLinearInterpolation(
                make_begin_iterator(myObject, &MyData::funcX),
                make_end_iterator(myObject, &MyData::funcX, myObject->nMax()),
                make_begin_iterator(myObject, &MyData::funcY),
                make_end_iterator(myObject, &MyData::funcY, myObject->nMax()),
                x);

    and if any of the arrays (x or y or both) is a regular vector<double> - use
    the usual semantics 

    vector<double> vx, vy;
    // vx, vy populated
    double y = smallLinearInterpolation(
                vx.begin(),
                vx.end(),
                vy.begin(),
                vy.end(),
                x);

    this interpolation method is not accessing more funcX(i) than necessary -- 
    it is checking the funcX() on the boundaries first and if the value is inside 
    uses bisection to converge to the interval where x lives. And then it calls 
    funcY() exactly two times.

    **********/




DRLIB_END_NAMESPACE

#endif //SMALL_LINEAR_INTERPOLATOR_HPP
