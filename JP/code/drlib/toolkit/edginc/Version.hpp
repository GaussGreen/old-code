#ifndef EDG_VERSION_H
#define EDG_VERSION_H
#include <string>
using namespace std;
DRLIB_BEGIN_NAMESPACE

class TOOLKIT_DLL CVersion
{
public:
    static const int    DRLIB_BUILD_DAYCOUNT;
    static const int    DRLIB_VER_MAJOR;
    static const int    DRLIB_VER_MINOR;
    static const string DRLIB_PLATFORM;

    static string DRLibVersion(void);

};


DRLIB_END_NAMESPACE
#endif
