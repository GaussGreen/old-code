/***************************************************************************
 *	SCCS Keyword Information
 *	------------------------
 *	Module name	:  mbsbase.c
 *	Company Name	:  JP Morgan Securities Inc.
 *	Author		:  Davis (Shuenn-Tyan) Lee
 *			   Derivatives Research
 *	Version		:  1.55
 *	Extracted	:  6/6/97 at 17:55:44
 *	Last Updated	:  6/6/97 at 17:55:42
 ***************************************************************************
 *      Base library for DR/MBS analytics (err handling, etc.)
 *      NB: error/mesg handling now uses GTO functions exclusively
 *
 *      Copyright 1996 J.P. Morgan & Co. Incorporated. All rights reserved.
 ***************************************************************************/
#include "mbsbase.h"
#include "mbsutils.h"

/****************************************************************************
 *  STATIC vars/funcs
 ***************************************************************************/
static TBoolean mbsInitDone = FALSE;         /* set when lib init done */
static TBoolean mbsEnableStdOutMsgs = TRUE;  /* enables stdout msgs */
static TBoolean mbsEnableLogMsgs = FALSE;    /* enables logfile msgs */
static char mbsAppLabel[LINESTRLEN] = DEF_MBSBASE_LIBNAME;
                                             /* label for msgs */
static char mbsLastErrMesg[ERRSTRLEN] = "";  /* most recent error msg */
static char mbsRunIdStr[L_tmpnam] = "";      /* unique string for run */
static TBoolean mbsEnableGtoTrapping = TRUE;     /* en/disables error trapping */
static TBoolean mbsNowInGtoTrap  = FALSE;    /* signals that MbsMesg() has
                                              * been called from trap routine*/

PRIVATE void EnsureMbsInit(void);
PRIVATE
TBoolean MbsTrapGtoErr
    (char    *errMesg,
     void    *callBackData);
PRIVATE
int GetUniqueTimeString
   (long    maxStrLen,          /* (I) max length of uniqueStr */
    char   *uniqueStr);         /* (O) unique string will be copied to 
                                 * this string (must already exist;
                                 * should have length >= L_tmpnam) */

/****************************************************************************
 * Vars shared w/mbsutils module
 ****************************************************************************/
char mbsPrimaryDataDir[MAXPATHLEN] = "";  /* first dir to search
						      for data	*/
char mbsFimsDataDir[MAXPATHLEN] = "";  /* dir to search for IMS
						   data		*/
char mbsDrDataDir[MAXPATHLEN]   = "";  /* dir to search for DR
						   data		*/

/****************************************************************************
 *  PUBLIC FUNCTIONS
 ***************************************************************************/

/****************************************************************************
 *  MbsBaseInit
 *  Newer (8/96) version of init function for mbsbase library
 *  Calling code should call this lib when starting up
 *  Notes:
 *   - can turn on/off msgs to stdout & log file independently
 *   - uses our code to send msgs to stdout
 *   - uses GTO funcs to send msgs to log file
 *   (and therefore this func will turn on/off GTO error logging)
 ***************************************************************************/
void MbsBaseInit
   (char     *appLabel,	        /* (I) brief label (for output msgs) */
    TBoolean  enableStdOutMsgs, /* (I) TRUE=enable output msgs to stdout */
    TBoolean  enableLogMsgs,    /* (I) TRUE=enable output msgs to logfile */
    char     *logFileName)      /* (I) Log file name (if logging enabled);
                                 * if NULL or empty, current/default
                                 * GTO log file name will be used;
                                 * SPECIAL CASE: if name begins with '@',
                                 * will use subsequent label together with
                                 * unique run id to make log file name */
{
    F_INIT("MbsBaseInit");
    char    fullLogName[MAXPATHLEN];
    char    tmpDirName[MAXPATHLEN];
    char    defOutputPath[MAXPATHLEN];
    long    iPos;


    /* Can only be called once */
    if( mbsInitDone )
    {
        return;
    }

    /* set flag to show init done 
     * (blocks infinite loop if any other base library func 
     * used, since that func might also try to init library) */
    mbsInitDone = TRUE;

    /* Get id string for run */
#ifdef UNIX_SKIP
    XONFAIL(GetUniqueTimeString(L_tmpnam,
                                mbsRunIdStr));
#endif
#ifdef DOS
    /* @@ don't yet know how this works on PC */
    strcpy(mbsRunIdStr,"");
#endif

    /* enable GTO error trapping */
    mbsEnableGtoTrapping = TRUE;
    mbsNowInGtoTrap = FALSE;
//    GtoErrMsgAddCallback(MbsTrapGtoErr,
//                         TRUE,
//                         NULL);

    /* enable msgs to stdout? */
    mbsEnableStdOutMsgs = enableStdOutMsgs;
    
    /* Enable msgs to log file? */
    mbsEnableLogMsgs = enableLogMsgs;
    /* (either way, alter GTO logging flag accordingly) */
    if( mbsEnableLogMsgs )
    {
        GtoErrMsgOn();
    }
    else
    {
        GtoErrMsgOff();
    }

    /* (if logging) Set/change log file name? */
    if(mbsEnableLogMsgs &&
       logFileName ISNT NULL &&
       logFileName[0] ISNT '\0')
    {
        /* if user wants us to create log file name using
         * (unique) run id string: */
        if( logFileName[0] IS '@' )
        {
            XONFAIL(MbsGetDefOutputPath(MAXPATHLEN-1,defOutputPath));
#ifdef UNIX
            sprintf(fullLogName,
                    "%s%c%s.%s",
                    defOutputPath,
                    MBS_PATH_DELIM,
                    mbsRunIdStr,
                    logFileName+1);
#endif
#ifdef DOS
            sprintf(fullLogName,
                    "%s%c%s%s",
                    defOutputPath,
                    MBS_PATH_DELIM,
                    mbsRunIdStr,
                    logFileName+1);
#endif
        }
        else 
        {
            
            /* if no explicit path in name, append home dir */
            str_npos(logFileName,MBS_PATH_DELIM,1,&iPos);
            if(iPos < 0)
            {
                XONFAIL(MbsGetDefOutputPath(MAXPATHLEN-1,defOutputPath));
#ifdef UNIX
                sprintf(fullLogName,"%s/%s",
                        defOutputPath,
                        logFileName);
#endif
#ifdef DOS
                sprintf(fullLogName,"%s%c%s",
                        defOutputPath,
                        MBS_PATH_DELIM,
                        logFileName);
#endif
            }
            else
            {
                strncpy(fullLogName,logFileName,MAXPATHLEN-2);
            }
        }
        GtoErrMsgFileName(fullLogName,TRUE); 
    }


    /* save global run id */
    if ( appLabel )
    {
	strncpy(mbsAppLabel,appLabel,LINESTRLEN-1);
    }

    /* note if data directories defined 
     * (FIMS data dir only avail on UNIX) */
#ifdef UNIX
    if( getenv("DATA") )
    {
	strncpy(tmpDirName,getenv("DATA"),MAXPATHLEN-1);
        XONFAIL(MbsSetFimsDataDir(tmpDirName));
    }
#endif

    /* For unix or PC: */
    if ( getenv("DRDATA") ) 
    {
	strncpy(tmpDirName,getenv("DRDATA"),MAXPATHLEN-1);
        XONFAIL(MbsSetDrDataDir(tmpDirName));
    }

    /* @@ later, I'd prefer that all functions refuse to run
     * if we haven't found a valid DRDATA directory */

    /* reset behavior linked-list ptr */
    p_behavior_start = NULL;  

    status = SUCCESS;
  done:
    return;
}


/*****************************************************************************
 *  MbsMesg
 *  Newer version (2/97) of mesg function
 *  Now uses GTO msg functions exclusively
 *****************************************************************************/
void MbsMesg
   (long      msgStatus,         /* (I) MINFO=info mesg(log), MWARN=warning
				   (stderr), MERROR=error(stderr) */
    TBoolean timeStamp,         /* (I) false=no timestamp, true=also print
				   timestamp */
    char    *mesg)              /* (I) msg to print */
{
    char   currTime[64],
	   locMesg[ERRSTRLEN] = "",
	   tmpMesg[ERRSTRLEN] = "",
	   fullMesg[ERRSTRLEN] = "";

    /* !! Note that we deliberately don't call EnsureMbsInit() here,
       since in some cases we may want to use this function before init
       was done */

    /*
     * Formatting of message 
     */

    /* make local copy of msg */
    if( mesg )
    {
	strncpy(locMesg,mesg,ERRSTRLEN);
    }
    /* must always have something... */
    if ( locMesg[0] IS '\0' ) 
    {
	strcpy(locMesg,"unknown error --- there is no error message");
    }
    /* trim any trailing CR */
    if ( locMesg[strlen(locMesg)-1] IS '\n' )
    {
	locMesg[strlen(locMesg)-1] = '\0';
    }
    /* format fuller mesg */
    switch (msgStatus) {
    case MINFO:
	sprintf(fullMesg, "%s", locMesg);
	break;
    case MWARN:
	sprintf(fullMesg, "WARNING from %s:\n%s", mbsAppLabel, locMesg);
	break;
    case MERROR:
	sprintf(fullMesg, "ERROR in %s:\n%s", mbsAppLabel, locMesg);
	break;
    case MERRTRACE:
	sprintf(fullMesg, "Error traceback: [%s]", locMesg);
	break;
    default:
	sprintf(fullMesg, "Message from %s:\n%s", mbsAppLabel, locMesg);
	break;
    }
    /* add timestamp */
    if (timeStamp)
    {
        GetTime(currTime);
	strncpy(tmpMesg,fullMesg,ERRSTRLEN-1);
	sprintf(fullMesg,"%s\n%s",currTime,tmpMesg);
    }

    /* Save most recent error msg 
     */
    if ( msgStatus IS MERROR ) 
    {
	strncpy(mbsLastErrMesg,fullMesg,ERRSTRLEN-1);
    }
    /* for error trace messages, append to end of message */
    else if ( msgStatus IS MERRTRACE )
    {
        strncat(mbsLastErrMesg,
                "\n",
                MAX(0,ERRSTRLEN-3-strlen(mbsLastErrMesg)));
        strncat(mbsLastErrMesg,
                fullMesg,
                MAX(0,ERRSTRLEN-3-strlen(mbsLastErrMesg)));
    }

    /* Output msg to stdout? 
     * (disabled for PC use)
     */
#ifndef DOS
    if ( mbsEnableStdOutMsgs ) 
    {
        fprintf(stdout,"%s\n",fullMesg);
    }
#endif

    /* Output msg to log file? 
     * (but don't do so if MbsMesg() called from trap, to prevent
     * infinite loop) */
    if(mbsEnableLogMsgs &&
       (!mbsNowInGtoTrap) )
    {
        /* to prevent infinite loop, disable GTO error trapping
         * (since here, our error handler is calling GTO's,
         * so we don't want GTO to then call ours...) */
        mbsEnableGtoTrapping = FALSE;
        GtoErrMsg("%s\n",fullMesg);
        /* when done, turn on trapping again */
        mbsEnableGtoTrapping = TRUE;
    }
}


/****************************************************************************
 * MbsTrapGtoErr
 * Traps calls to GTO error handler, to allow us to 
 * save the (last) message, and to output to stdout,
 * if needed
 ****************************************************************************/
PRIVATE
TBoolean MbsTrapGtoErr
    (char    *errMesg,
     void    *callBackData)
{
    /* only do this if trapping enabled 
     * (i.e., only if code really called GtoErrMesg()) */
    if( mbsEnableGtoTrapping )
    {
        /* signal that we're calling MbsMesg() from trap */
        mbsNowInGtoTrap = TRUE;
        /* call our err mesg handler (saves mesg, does stdout if needed) */
        MbsMesg(MINFO,FALSE,errMesg);
        /* done with trap */
        mbsNowInGtoTrap = FALSE;
    }
    /* always return TRUE, to let GTO error handler do its normal thing */
    return TRUE;
}



/****************************************************************************
 *  mbs_init
 *  Older version of mbsbase init function, for backward compatibility
 ***************************************************************************/
void mbs_init(char *app_label,	/* (I) brief label for app/run (for output
				   msgs) */
	      long enable_mesgs,	/* (I) true enables mbs_mesg() output to
				   stdout/stderr, false prevents any
				   stdout/stderr output			*/
	      char *logfile)	/* (I) name for log file (NULL for none,
				   "" for default) */
{
    /* Just call the newer function */
    MbsBaseInit(app_label, 
                enable_mesgs, 
                (logfile ISNT NULL) && (logfile[0] ISNT '\0'),
                logfile);
}


/*****************************************************************************
 *  mbs_mesg
 *  Older function to output message (either normal or error msg)
 *  Included for backward compatibility
 *****************************************************************************/
void mbs_mesg(long msgStatus,	/* (I) MINFO=info mesg(log), MWARN=warning
				   (stderr), MERROR=error(stderr) */
	      long timestamp,	/* (I) false=no timestamp, true=also print
				   timestamp */
	      char *mesg)	/* (I) msg to print */
{
    /* Just call new function */
    MbsMesg(msgStatus, (TBoolean) timestamp, mesg);
}


/*****************************************************************************
 *  MbsDbg
 *  Outputs string only in debug mode
 *  (Simplified, special-purpose version of mbs_mesg)
 *****************************************************************************/
/*ARGSUSED*/
void MbsDbg(char *msg)
{
#ifdef DEBUG_MBS
    mbs_mesg(MINFO,0,msg);
#endif
}

/*****************************************************************************
 *  MbsGetLastErrMesg
 *  Copies most recent err mesg to output string
 *****************************************************************************/
int MbsGetLastErrMesg(long maxlen,            /* (I) max # chars to copy */
			 char *last_errmesg)	/* (O) mesg copied to this
						   string */
{
    F_INIT("MbsGetLastErrMesg");

    /* make sure init has been done */
    EnsureMbsInit();

    XONTRUE( last_errmesg IS NULL, "Failed to copy err msg--null ptr" );
    strncpy(last_errmesg,mbsLastErrMesg,maxlen);

    F_END;
}

/*****************************************************************************
 *  MbsSetPrimaryDataDir
 *  Sets/overrides path to Primary data directory 
 *****************************************************************************/
int MbsSetPrimaryDataDir
   (char *dataDir)          /* (I) name of directory */
{
    F_INIT("MbsSetPrimaryDataDir");

    /* make sure init has been done */
    EnsureMbsInit();

    XONTRUE(dataDir IS NULL ||
            dataDir[0] IS '\0',
            "Failed to supply path for Primary data directory");
    strncpy(mbsPrimaryDataDir,dataDir,MAXPATHLEN-1);
    /* ensure trailing (back)slash */
    XONFAIL( MbsAddTrailSlash(mbsPrimaryDataDir) );

    F_END;
}


/*****************************************************************************
 *  MbsSetDrDataDir
 *  Sets/overrides path to DRDATA data directory 
 *****************************************************************************/
int MbsSetDrDataDir
   (char *dataDir)          /* (I) name of directory */
{
    F_INIT("MbsSetDrDataDir");
    char    tmpFileName[MAXPATHLEN];

    /* make sure init has been done */
    EnsureMbsInit();

    XONTRUE(dataDir IS NULL ||
            dataDir[0] IS '\0',
            "Failed to supply path for DR DATA directory");
    strncpy(mbsDrDataDir,dataDir,MAXPATHLEN-1);
    /* ensure trailing (back)slash */
    XONFAIL( MbsAddTrailSlash(mbsDrDataDir) );

    /* Now, we check for signature file in this specific DRDATA dir,
     * and declare an error if not found */
    XONFAIL(MbsBuildFileName
            (mbsDrDataDir,
             DRDATA_SIG_FILE,
             tmpFileName));
    XONTRUE( ! mbs_filefound(tmpFileName),
            "DRDATA data directory is not valid--signature file not found");

    F_END;
}


/*****************************************************************************
 *  MbsSetFimsDataDir
 *  Sets/overrides path to FIMS $DATA data directory 
 *****************************************************************************/
int MbsSetFimsDataDir
   (char *dataDir)          /* (I) name of FIMS $DATA directory */
{
/*
    F_INIT("MbsSetFimsDataDir");
*/
/*    char   tmpDirName[MAXPATHLEN]; */
	
	printf ("Procedure not implemented\n");
	return FAILURE;

    /* make sure init has been done */
/*    EnsureMbsInit();

    XONTRUE(dataDir IS NULL ||
            dataDir[0] IS '\0',
            "Failed to supply path for FIMS DATA directory");

    strncpy(tmpDirName,dataDir,MAXPATHLEN-1);
    XONFAIL(MbsRemoveTrailSlash(tmpDirName,FALSE));
    if(strlen(tmpDirName) > MAXPATHLEN-1-strlen(FIMS_DATADIR_SUFFIX))
    {
        tmpDirName[MAX(0,MAXPATHLEN-1-strlen(FIMS_DATADIR_SUFFIX))] = '\0';
    }
    strcat(tmpDirName,FIMS_DATADIR_SUFFIX);
    strncpy(mbsFimsDataDir,tmpDirName,MAXPATHLEN-1);
*/
    /* ensure trailing (back)slash */
/*    XONFAIL(MbsAddTrailSlash(mbsFimsDataDir));
*/
/*
	F_END;
*/
}


/***************************************************************
 * MbsGetRunIdString
 * Returns pointer to run id string for this job
 ***************************************************************/
int EXPORT MbsGetRunIdString
   (char  **runIdStr)             /* (O) returns ptr to run id string */
{
    /* make sure init has been done */
    EnsureMbsInit();

    *runIdStr = (char*) &mbsRunIdStr;
    return(SUCCESS);
}



/*****************************************************************************
 *  LOCAL FUNCTIONS
 ****************************************************************************/

/****************************************************************************
 *  EnsureMbsInit
 *  Local function called by every function in module to
 *  guarantee that initialization has been done
 ****************************************************************************/
PRIVATE
void EnsureMbsInit(void)
{
    /* if this function called, it means application 
       forgot to explicitly init base lib;
       given this, be conservative and init the
       lib in "quiet" mode (no output of any kind) */
    MbsBaseInit("mbsbase",FALSE,FALSE,NULL);
}


/**********************************************************************
 * GetUniqueTimeString
 * Returns unique string based on current time
 * Useful for unique run labels
 **********************************************************************/
PRIVATE
int GetUniqueTimeString
   (long    maxStrLen,          /* (I) max length of uniqueStr */
    char   *uniqueStr)          /* (O) unique string will be copied to 
                                 * this string (must already exist;
                                 * should have length >= L_tmpnam) */
{
    F_INIT("GetUniqueTimeString");
    static char dumbTempNameDos[] = "\temp.fil";
    static char dumbTempNameUnix[] = "/temp.fil";
    char   *ptempnam = NULL;

    XONTRUE(maxStrLen < L_tmpnam,
            "Output string too short--must have length = L_tmpnam");

    /* tempnam uses malloc() to store created name */
#ifdef UNIX
#ifdef SUN
    ptempnam = tempnam("/", "A");	    
#endif
#ifdef HP
    /* For now, have problems w/temp name on HP */
    ptempnam = dumbTempNameUnix;
#endif
#endif

#ifdef DOS
    /* For now, doubt this works on DOS */
    ptempnam = dumbTempNameDos;
#endif
    XONTRUE(ptempnam IS NULL,
            "Bizarre: Call to UNIX tempnam() function failed");
    /* create unique identifier for filenames */
    strncpy(uniqueStr, ptempnam + 2, maxStrLen);   

    status = SUCCESS;
  done:
    /* DON'T USE GTO FREE() here: de-alloc temp space created 
     * by tempnam() using malloc() */
#ifdef UNIX
#ifdef SUN
    if(ptempnam ISNT NULL)
    {
        free(ptempnam);
    }
#endif
#endif
    return status;
}

