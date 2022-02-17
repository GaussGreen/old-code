//----------------------------------------------------------------------------
//
//   Group       : QR&D Core Analytics Team
//
//   Filename    : Equals.hpp
//
//   Description : Templates for equality tests
//
//   Author      : Richard Appleton
//
//   Date        : 1tth November 2005
//
//----------------------------------------------------------------------------

#ifndef EQUALS_HPP
#define EQUALS_HPP

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"


DRLIB_BEGIN_NAMESPACE


template <class T>
bool equalsValue(const smartConstPtr<T>& lhs, const smartConstPtr<T>& rhs)
{
    if (lhs.get() == NULL && rhs.get() == NULL)
        return true;
    else if (lhs.get() == NULL || rhs.get() == NULL)
        return false;
    else
        return lhs.get()->equals(rhs.get());
}


template<class T>
bool equalsClass(const smartConstPtr<T>& lhs, const smartConstPtr<T>& rhs)
{
    if (lhs.get() == NULL && rhs.get() == NULL)
        return true;
    else if (lhs.get() == NULL || rhs.get() == NULL)
        return false;
    else
        return lhs.get()->getClass() == rhs.get()->getClass();
}

template<class T>
bool NumericArrayEquals( 
    const smartConstPtr<T>& lhs, 
    const smartConstPtr<T>& rhs )
{
    if (lhs.get() == NULL && rhs.get() == NULL)
        return true;
    else if (lhs.get() == NULL || rhs.get() == NULL)
        return false;
    else {
        if ( lhs->size() != rhs->size() ) { 
            return false;
        }

        int len = lhs->size();
        for (int i = 0; i < len; i++) {
            if ( ! Maths::equals( (*lhs)[i], (*rhs)[i] ) ) {
                return false;
            }
        }
    }
    return true;
}

template<class T>
bool StringArrayEquals( 
    const smartConstPtr<T>& lhs, 
    const smartConstPtr<T>& rhs )
{
    if (lhs.get() == NULL && rhs.get() == NULL)
        return true;
    else if (lhs.get() == NULL || rhs.get() == NULL)
        return false;
    else {
        if ( lhs->size() != rhs->size() ) { 
            return false;
        }

        int len = lhs->size();
        for (int i = 0; i < len; i++) {
            if ( (*lhs)[i] != (*rhs)[i] ) {
                return false;
            }
        }
    }
    return true;
}

template<class T>
bool GenericArrayEquals( 
    const smartConstPtr<T>& lhs, 
    const smartConstPtr<T>& rhs )
{
    // This method will give a compilation error if it is used
    // with types that do not support equals(T*) method.
    // For instance, the DoubleArrayEquals above should be used for
    // array of doubles.
    if (lhs.get() == NULL && rhs.get() == NULL)
        return true;
    else if (lhs.get() == NULL || rhs.get() == NULL)
        return false;
    else {
        if ( lhs->size() != rhs->size() ) { 
            return false;
        }

        int len = lhs->size();
        for (int i = 0; i < len; i++) {
            if ( !(*lhs)[i]->equals((*rhs)[i].get() )) {
                return false;
            }
        }
    }
    return true;
}

DRLIB_END_NAMESPACE

#endif // EQUALS_HPP
