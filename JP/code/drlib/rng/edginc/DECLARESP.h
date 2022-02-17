//
// C++ Interface: DECLARESP
//
// Description:
//
//
// Author: Vladimir A Grebinskiy <vladimir.a.grebinskiy@jpmorgan.com>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef DECLARESP_H
#define DECLARESP_H

#include "edginc/coreConfig.hpp"
#include "edginc/refCountPtr.hpp"

CORE_BEGIN_NAMESPACE

#define DECLARESP(T) \
    typedef refCountPtr<T> T##SP;

#define FORWARD_DECLARESP(T) \
    class T; \
    DECLARESP(T);
    
CORE_END_NAMESPACE

#endif

