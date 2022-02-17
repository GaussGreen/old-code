#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/Version.hpp"
#include "edginc/buildnum.hpp"

DRLIB_BEGIN_NAMESPACE

// buildnum.hpp contains a running total of all the build's we've ever made
// it needs to be in a separate file so that we get a different number even
// when we rebuild the same version

const int CVersion::DRLIB_VER_MAJOR       = 6;
const int CVersion::DRLIB_VER_MINOR       = 6;
/* to reset DRLIB_BUILD_DAYCOUNT leave it blank between the '=' and the ';' 
   However, if you are about to run the overnight script using 'attachlabel'
   then you must put in zero (ie a value) here and then relabel the file
 */
const int CVersion::DRLIB_BUILD_DAYCOUNT  = 33;

#ifdef WIN32
#if defined (_MSC_VER) && (_MSC_VER == 1310)
const string CVersion::DRLIB_PLATFORM = "NT vc 7.1";
#else
const string CVersion::DRLIB_PLATFORM = "NT";
#endif
#elif defined(LINUX)
const string CVersion::DRLIB_PLATFORM = "Linux";
#else
const string CVersion::DRLIB_PLATFORM = "Solaris";
#endif

string CVersion::DRLibVersion(void)
{
    return Format::toString("%d.%d.%d.0 devel (Build %d) [%s]",
                            CVersion::DRLIB_VER_MAJOR,
                            CVersion::DRLIB_VER_MINOR,
                            CVersion::DRLIB_BUILD_DAYCOUNT,
                            EDRLIB_BUILD_NUMBER,
                            CVersion::DRLIB_PLATFORM.c_str());
}

DRLIB_END_NAMESPACE
