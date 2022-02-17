/**
 * @file FORWARD_DECLARE.hpp
 */

#ifndef QLIB_FORWARD_DECLARE_H
#define QLIB_FORWARD_DECLARE_H

#include "edginc/DECLARE.hpp"

// NB You must use these macros inside a DRLIB_BEGIN_NAMESPACE etc

#define FORWARD_DECLARE(T)                              \
    class T;                                            \
    DECLARE(T)

#define FORWARD_DECLARE_REF_COUNT(T)                    \
    class T;                                            \
    DECLARE_REF_COUNT(T)

#endif
