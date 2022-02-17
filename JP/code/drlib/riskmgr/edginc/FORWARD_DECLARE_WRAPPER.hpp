/**
 * @file FORWARD_DECLARE_WRAPPER.hpp
 */

#ifndef QLIB_FORWARD_DECLARE_WRAPPER_H
#define QLIB_FORWARD_DECLARE_WRAPPER_H

#include "edginc/MarketObject.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

#define FORWARD_DECLARE_WRAPPER(T)        \
    FORWARD_DECLARE(T);                   \
    typedef MarketWrapper<T> T##Wrapper;

DRLIB_END_NAMESPACE

#endif
