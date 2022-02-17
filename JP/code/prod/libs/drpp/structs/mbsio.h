#if !defined(mbsio_h)
#define mbsio_h

#if _MSC_VER >= 1000
#pragma once
#endif // _MSC_VER >= 1000
#pragma warning(disable:4786)

#include "drinput.h"
#include "cdate.h"

typedef KInput DRInput;
typedef KInputMap DRInputMap;
typedef KInputVector DRInputVector;

typedef DRInput MbsIO;
typedef DRInputMap MbsIOMap;
typedef DRInputVector MbsIOVector;

//#include "drdate.h"

void MakeCounted (DRInputVector& vec, TDate* &dates);
void MakeCounted (DRInputVector& vec, double* &values);
void MakeCounted (DRInputVector& vec, int* &values);
void MakeCounted (DRInputVector& vec, char* &values);

double toNumber (DRInput& input);

#endif 
