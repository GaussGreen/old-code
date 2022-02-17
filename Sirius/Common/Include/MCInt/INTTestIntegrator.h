
#ifndef __INTTESTINTEGRATOR_H__
#define __INTTESTINTEGRATOR_H__

//#include <eqspmcint/StdAfx.h>

#include "INTSTL.h"

// Perform single integration
 void testintegrator(void);

// Parameters k: dimension, class of option in {1 = Asian, 2 = Basket}, Q: number of instances
 void testsystematic(int k, int optionclass, int Q, double corr = 0.0){};



#endif
