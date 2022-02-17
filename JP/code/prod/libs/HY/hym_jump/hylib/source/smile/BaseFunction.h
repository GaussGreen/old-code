// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 1999 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// 1.0 10/7/99 Afshin Bayrooti
// $Header$
//

#if ! defined(BASEFUNCTION)
#define BASEFUNCTION




// Helper class for the root solver
class AbstractFunctorWrapper {
public:
    virtual double operator()(double) const = 0;
};


template<class Functor>
class ConcreteFunctorWrapper : public AbstractFunctorWrapper
{
public:
    ConcreteFunctorWrapper(Functor& function) : m_function( function ) {};

    double operator()(double x) const {
        return m_function(x);
    }
        
private:
    Functor& m_function;
};



#endif
