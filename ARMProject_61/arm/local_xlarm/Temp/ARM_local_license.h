#ifndef _ARM_LOCAL_LICENSE_H
#define _ARM_LOCAL_LICENSE_H


#include <libCCTools++\CCString.h>


#ifdef __cplusplus
extern "C"
{
#endif	/* __ cplusplus */
extern int isAuthorized ();
#ifdef __cplusplus
}
#endif	/* __cplusplus */


VECTOR<CCString> GetIpAddress(CCString hostname);
CCString PcName();
VECTOR<CCString> GetPcIp();
int GetLicense();


#endif /* _ARM_LOCAL_LICENSE_H */
