/* this file contains the actual definitions of */
/* the IIDs and CLSIDs */

/* link this file in with the server and any clients */


/* File created by MIDL compiler version 5.01.0164 */
/* at Fri May 18 15:56:40 2007
 */
/* Compiler settings for .\local_DLLARM.idl:
    Oicf (OptLev=i2), W1, Zp8, env=Win32, ms_ext, c_ext
    error checks: allocation ref bounds_check enum stub_data 
*/
//@@MIDL_FILE_HEADING(  )
#ifdef __cplusplus
extern "C"{
#endif 


#ifndef __IID_DEFINED__
#define __IID_DEFINED__

typedef struct _IID
{
    unsigned long x;
    unsigned short s1;
    unsigned short s2;
    unsigned char  c[8];
} IID;

#endif // __IID_DEFINED__

#ifndef CLSID_DEFINED
#define CLSID_DEFINED
typedef IID CLSID;
#endif // CLSID_DEFINED

const IID IID_IARMModule = {0xCA2D4B2C,0x934D,0x4643,{0xBA,0x08,0xCA,0xED,0xAB,0x80,0xB9,0xB8}};


const IID LIBID_LOCAL_DLLARMLib = {0x09044F72,0x94F4,0x460E,{0xAA,0x8F,0xDE,0x46,0x2B,0xB5,0x2D,0x7C}};


const CLSID CLSID_ARMModule = {0x47CC0BC7,0x4EDF,0x4853,{0x83,0x8C,0x89,0xB4,0x7E,0x57,0x9A,0xD4}};


#ifdef __cplusplus
}
#endif

