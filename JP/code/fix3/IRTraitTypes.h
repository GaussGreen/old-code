// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2000 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// 1.0 08/27/2003 Victor Paskhaver
//
// $Header$
//

#pragma once

#if ! defined(_LEGO_IR_TYPES_)
#define _LEGO_IR_TYPES_

#include <vector>
#include "assert.h"

using namespace std;

namespace IR {

template<class T>
class Matrix {
    int m_rows;
    int m_cols;
    std::vector<T> m_data;
public:
    typedef T value_type;
    Matrix() {
        m_rows = 0; 
        m_cols = 0;
    }
    int rows() const { return m_rows; }
    int cols () const  { return m_cols; }
    void resize(int rows, int cols ) {
        m_data.resize(rows*cols);
        m_rows = rows;
        m_cols = cols;
    }
    T &at(int row, int col)
    {
        return m_data.at(col + row*m_cols);
    }
    
};


struct Date {

    // following spills over month end, modified following
    // stays in current month if a bad day
    enum EMonthEndRule { FOLLOWING, MODIFIEDFOLLOWING };
    static char noLeapYear[13]; // = {0,31,28,31,30,31,30,31,31,30,31,30,31};
    static char leapYear[13]; // =  {0,31,29,31,30,31,30,31,31,30,31,30,31};

    short year;
    char  month;
    char  day;

    Date() {}
    
    Date( long date ) {
        day = date % 100;
        month = date % 10000 / 100;
        year = date / 10000;
    }
    Date (const Date& date)
    {
        year = date.year;
        month = date.month;
        day = date.day;
    }

    void AddMonths(int numMonths, EMonthEndRule rule)
    {


    }

    void AddYears(int numYears, EMonthEndRule rule)
    {
        year+=numYears;
        char daysInMonth;

        // check leap year hasn't caused day to become invalid
        if (IsLeapYear())
            daysInMonth = leapYear[month];
        else
            daysInMonth = noLeapYear[month];

        if (day > daysInMonth)
        {
            if (rule == MODIFIEDFOLLOWING)
            {
                // move backwards to last good day of month
                day = daysInMonth;    
            }
            else  // move to start of next month if bad day
            {
                day = 1;
                char nextMonth = month + 1;
                if (nextMonth > 12)
                {
                    month = 1;
                    year++;
                }
            }
        }
    }

    bool IsLeapYear() const
    {
        if (year%4 != 0)
            return false;  // not divisible by 4
        if (year%100 != 0)
            return true;  // divisible by 4 but not 100
        if (year%400 != 0)
            return false;  // divisible by 100 but not 400

        return true;  // divisible by 400, so is a leap year
    }

    bool IsValid() const
    {
        if (month < 1 || month > 12)
            return false;
        if (year < 0 || year > 3000)
            return false;

        char numDays = IsLeapYear() ? leapYear[year] : noLeapYear[year];
        if (day < 1 || day > numDays)
            return false;

        return true;
    }
    
    operator long() const {
         return ((unsigned long) year ) * 10000 + 
                ((unsigned long) month) * 100   + 
                  (unsigned long) day;
    }

    bool operator < (const Date& date)
    {
        if (year < date.year)
            return true;
        if (year == date.year && month < date.month)
            return true;
        if (year == date.year && month == date.month && day < date.day)
            return true;

        return false;
    }
};


// Execute returns the result of the calculation.  Therefore, if an exception occurs,
// the API layer needs to throw an exception.  

// With finalization, this can be overcome by modifying the Finalize(...) function
// on the IRTraits.h to return an int and then modifying
//      static void ClassTraits::finalize( T &that ) {
//                Finalize( that );
//      }
// to check the result of Finalize and throw an exception on error.

class Error {
protected:
public:
    Error(const char * str) : m_what(str) {}
    Error(const std::string &str) : m_what(str) {}

    std::string what() const { return m_what; }
private:
    std::string m_what;
    Error(){}
};

// This version of smart pointer stores the reference count outside of the object.  
// That way, non-reference counted objects can be reference counted without modifying 
// the IR / Hybrids library.  However, when defining new objects, one needs to be
// careful.  For example, suppose you have an object A which includes an object B.  
// There are three ways to handle this without any problems.  These are:
// 1. Include one object into another (ie. copy it).
// 2. Use reference counted pointers, ie. A stores a shared_ptr to B.
// 3. Carefully set up the C style pointers.  The current implementation doesn't
//      support this yet but with additional methods, C users can manually
//      increase and decrease the reference count of the object.

class shared_count {
    int m_count;
protected:
    virtual void destroy() = 0;
public:
    shared_count() {
        m_count = 0;
    }
    void addRef() {
        ++m_count;
    }
    void release() {
        if( --m_count == 0 )
            destroy();
    }
};

template <class T>
class shared_count_typed : public shared_count
{
    T *m_ptr;
public:
    shared_count_typed( T* ptr) : m_ptr(ptr) {}

    void destroy() {
        delete m_ptr;
        delete this;
    }
};

struct shared_ptr_base {
    void *m_ptr;
    shared_count *m_count;
};
template<class T>
class shared_ptr : public shared_ptr_base
{
public:
    typedef T value_type;
    shared_ptr() {
        m_count = 0;
        m_ptr = 0;
    }
    template <class U>
    shared_ptr(U* ptr) {
        m_ptr = static_cast<T*>(ptr);
        m_count = new shared_count_typed<T>(ptr);
        m_count->addRef();
    }

    shared_ptr(int value)  // only for the meantime
    {
        m_count = 0;
        m_ptr = 0;
    }

    void operator = (const shared_ptr &src)
    {
        if(m_count)
            m_count->release();
        m_ptr = src.m_ptr;
        m_count = src.m_count;
        if(m_count)
            m_count->addRef();

    }
    

    template <class U>
    void operator = (U *ptr) 
    {
        if(m_count)
            m_count->release();
        m_ptr = static_cast<T*>(ptr);
        m_count = new shared_count_typed<T>(ptr);
        if(m_count)
            m_count->addRef();
    }

    ~shared_ptr() {
        if(m_count)
            m_count->release();
    }

    T* get() const {
        return (T*) m_ptr;
    }
    T* operator->() const {
        return (T*) m_ptr;
    }
    bool IsNull() const
    {
        if (m_count == 0 && m_ptr == 0)
            return true;
        else
            return false;
    }
    
};

// The following are utilities to convert between C arrays to std::vectors

template< int N>
struct CArrayCheck {
	enum { arraySize=N ? N : 1};
	template <class T>
		static bool providedTypeIsNotAnArray(T (&array)[arraySize]);
};
template <typename T>
inline int CArraySize(const T&arr)
{
	enum {arraySize = sizeof(arr) / sizeof(arr[0])};
	sizeof(CArrayCheck<arraySize>::providedTypeIsNotAnArray(arr));
    return arraySize;
}

template<typename S, typename T> 
void DynamicCArrayToVector( const S* src, int N, std::vector<T>& target ) {
    target.clear();
    target.resize( N );
    for ( int i = 0; i < N; ++i ) {
        target[i] = src[i];
    }
}

template<typename S, typename T> 
void CArrayToVector( const S &src, int N, std::vector<T>& target ) {
	enum {arraySize = sizeof(src) / sizeof(src[0])};
	sizeof(CArrayCheck<arraySize>::providedTypeIsNotAnArray(src));
	assert ( arraySize >= N );
    DynamicCArrayToVector( src, N, target );
}

template<typename S, typename T> 
void VectorToCArray( const std::vector<S>& src, T& target ) {
	enum {arraySize = sizeof(target) / sizeof(target[0])};
	sizeof(CArrayCheck<arraySize>::providedTypeIsNotAnArray(target));
	if ( arraySize < src.size() ) 
		throw Error("Static array size is less than source size");
	
    for ( int i = 0; i < src.size(); ++i ) {
        target[i] = src[i];
    }
}

template<typename S, typename T> 
void VectorToDynamicCArray( const std::vector<S>& src, T* &target ) {
    target = (T*)malloc(src.size()*sizeof(T));
    for ( int i = 0; i < src.size(); ++i ) {
        target[i] = src[i];
    }
}

inline std::string operator <<(std::string dst, int src)
{
    char buffer[128];
    sprintf(buffer, "%i", src);
    return buffer;
}

inline std::string operator <<(std::string dst, const char *src)
{
    dst += src;
    return dst;
}

inline std::string operator <<(std::string dst, std::string src)
{
    dst += src;
    return dst;
}


} // IR


#endif
