#include "imsl_inc.h"

#if defined( _MSC_VER )
#pragma warning( once : 4101 4102 4244 4305 )
#endif
/*
#define COMPUTER_Mac2ms
*/
    /* MACHINE DEPENDENT CODE */

    /* OPERATING SYSTEM VERSION */
/* MAR: 2002/5/3 - this is just a pain in the neck so maintain so admit 
   defeat now */
static Mchar *lv_os_version = 
#ifdef UNIX
    "UNIX";
#else
    "Windows";
#endif
#if 0
#if defined(COMPUTER_Mac2ms)
    "Apple System 7.5";
#endif
#if defined(COMPUTER_SUN)
    "SunOS 4.1.1";
#endif
#if defined(COMPUTER_SUN4)
#if defined(COMPUTER_SLRSXS)
    "Solaris 2.3";
#else
    "SunOS 4.1.3_U1";
#endif
#endif
#if defined(COMPUTER_GSUN4) 
     "SunOS 4.1.3_U1";
#endif
#if defined(COMPUTER_SUN5)
    "SunOS Version 5.5.1 (Solaris)";
#endif
#if defined(COMPUTER_PMXUX)
    "Ultrix V4.1";
#endif
#if defined(COMPUTER_DECOSF)
    "DEC OSF/1 V1.2";
#endif
#if defined(COMPUTER_MIPS)
    "RISC/os (UMIPS) 4.52";
#endif
#if defined(COMPUTER_RTXLXS)
    "AIX Version 3.2";
#endif
#if defined(COMPUTER_ALFAC_IEEE)
    "DEC System ALFAC AXP";
#endif
#if defined(COMPUTER_VAX)
    "VAX/VMS V5.3-1";
#endif
#if defined(COMPUTER_SGRUXS)
    "IRIX Release 5.2";
#endif
#if defined(COMPUTER_HP93C)
    "HP-UX C Release 7.0.3";
#endif
#if defined(COMPUTER_HP97C)
    "HP-UX Release 9.01";
#endif
#if defined(COMPUTER_HP98C)
    "HP-UX Release 7.0";
#endif
#if defined(COMPUTER_APLC)
    "Apollo Domain/IX Release 10.2";
#endif
#if defined(IMSL_MACHINE_NT)
    "Microsoft Windows NT Version 4.0";
#endif
#if defined(COMPUTER_IRIX_N32)
    "SGI IRIX 6.5";
#endif
#if defined(COMPUTER_LINUX)
    "Linux 2.2.5-22";
#endif
#endif

    /* C COMPILER VERSION */
/* MAR: 2002/5/3 - this is just a pain in the neck so maintain so admit 
   defeat now */
static Mchar *lv_compiler_version = "EDR Default Compiler";
#if 0
#if defined(COMPUTER_Mac2ms)
    "Metrowerks CodeWarrior Gold CW4 C/C++ 68K Release 1.1";
#endif
#if defined(COMPUTER_SUN)
    "SunOS 4.1.1 cc";
#endif
#if defined(COMPUTER_SUN4)
#  if defined(COMPUTER_SLRSXS)
      "SPARCompiler C 3.0";
#  elif defined(COMPUTER_SUN4ACC)
      "Sun ANSI C 2.01";
#  elif defined(COMPUTER_SUN4SC)
      "SPARCompiler C 2.0.1";
#  else
#    if !defined(COMPUTER_GSUN4)
      "SunOS 4.1.3 cc";
#    endif
#  endif
#endif
#if defined(COMPUTER_SUN5)
#  if defined(COMPUTER_SLRSACC)
       "Sun ANSI C Compiler for Solaris";
#  elif defined(COMPUTER_SLRSGCC)
       "GNU C 2.7.2.0";
#  endif
#endif
#if defined(COMPUTER_GSUN) || defined(COMPUTER_GSUN4)
   "GNU C 2.7.2.0";
#endif
#if defined(COMPUTER_PMXUX)
    "Ultrix V2.0 cc";
#endif
#if defined(COMPUTER_DECOSF)
    "DEC OSF/1 C compiler V1.2";
#endif
#if defined(COMPUTER_MIPS)
    "MIPS C compiler 2.20";
#endif
#if defined(COMPUTER_RTXLXS)
    "XL C 1.3.0.41";
#endif
#if defined(COMPUTER_ALFAC_IEEE)
    "DEC C Compiler Version 1.2 with/IEEE";
#endif
#if defined(COMPUTER_VAX)
			  /* it's a VAX using the /NOG_FLOAT option */
#if defined(COMPUTER_VAX) && !defined(COMPUTER_VAXG)
    "VAX C V3.1-051";
#else                     /* it's a VAX using the /G_FLOAT option */
    "VAX C V3.1-051 with /G_FLOAT";
#endif
#endif
#if defined(COMPUTER_SGRUXS)
    "IRIX C 3.18";
#endif
#if defined(COMPUTER_HP93C)
    "HP-UX Release 68.8.1.2 cc";
#endif
#if defined(COMPUTER_HP97C)
    "HP A.09.01, c89 Compiler";
#endif
#if defined(COMPUTER_HP98C)
    "HP-UX C Release 7.10";
#endif
#if defined(COMPUTER_APLC)
    "C compiler 68K Rev 6.7(316)";
#endif
#if defined(COMPUTER_NTLNT)
    "cl386 Version 8.00.3190a";
#endif
#if defined(COMPUTER_MIPNT)
    "mcl Version 8.00.081";
#endif
#if defined(COMPUTER_ALIBNT)
    "Microsoft CL, Version 10.00.5270";
#endif
#if defined(COMPUTER_ALFNT)
    "claxp Version 8.00.9B";
#endif
#if defined(COMPUTER_IRIX_N32)
    "MIPSpro Compiler Version 7.3.1.1m";
#endif
#if defined(COMPUTER_LINUX)
    "gcc version egcs-2.91.66 19990314/Linux (egcs-1.1.2 release)";
#endif
#endif

#define VERSION_BUF_LEN   255
static Mchar lv_library_version[VERSION_BUF_LEN]; 


static Mchar   *PROTO(l_new_string,(Mchar *string));
static Mchar   *PROTO(l_format_library_version,(void));

#ifdef ANSI
Mchar *imsl_version(Mint code)
#else
Mchar *imsl_version(code)
   Mint          code;
#endif
{
    switch (code) {
	case IMSL_OS_VERSION:
	    return l_new_string(lv_os_version);
	case IMSL_COMPILER_VERSION:
	    return l_new_string(lv_compiler_version);
	case IMSL_LIBRARY_VERSION:
	    return l_new_string(l_format_library_version());
	case IMSL_LICENSE_NUMBER:
	    return l_new_string(imsl_find_message(IMSL_LICENSE_NUMBER_ENTRY));
	default:
	    imsl_e1psh ("imsl_version");
	    imsl_e1sti (1, code);
	    imsl_e1sti (2, 1);
	    imsl_ermes (IMSL_TERMINAL, IMSL_UNKNOWN_OPTION);
	    imsl_e1pop ("imsl_version");
	    return NULL;
    }
}


#ifdef ANSI
static Mchar *l_new_string(Mchar *string)
#else
static Mchar *l_new_string(string)
    Mchar   *string;
#endif
{
    return strcpy((Mchar*)imsl_malloc(strlen(string)+1), string);
}

#ifdef ANSI
static Mchar * l_format_library_version (void) 
#else
static Mchar * l_format_library_version () 
#endif
{ 
#if !defined(IMSL_MAJOR_REV) || !defined(IMSL_MINOR_REV)
/* MAR: 2002/5/3 - I believe that IMSL_MAJOR_REV and IMSL_MINOR_REV would come
   from the version control label. Don't particular want a separate label for
   this code so manually define version here */
#define IMSL_MAJOR_REV 1
#define IMSL_MINOR_REV 0
#endif
    sprintf (lv_library_version, "EDR IMSL C/Math/Library Version %d.%d", 
             IMSL_MAJOR_REV, IMSL_MINOR_REV); 
    return (lv_library_version); 
}

