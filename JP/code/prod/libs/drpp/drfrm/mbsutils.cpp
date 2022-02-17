/***************************************************************************
 *	SCCS Keyword Information
 *	------------------------
 *	Module name	:  mbsutils.c
 *	Company Name	:  JP Morgan Securities Inc.
 *	Author		:  Davis (Shuenn-Tyan) Lee
 *			   Derivatives Research
 *	Version		:  1.73
 *	Extracted	:  6/18/97 at 13:06:13
 *	Last Updated	:  6/18/97 at 13:06:11
 ***************************************************************************
 *	Misc utilities for DR/MBS analytics base lib 
 *      (file i/o, memory, etc)
 *
 *      Copyright 1996 J.P. Morgan & Co. Incorporated. All rights reserved.
 ***************************************************************************/

#include "mbsutils.h"
#include "mbsbase.h"
#define DOS


/****************************************************************************
 *  EXTERN vars (shared from mbsbase)
 ****************************************************************************/
extern char mbsPrimaryDataDir[MAXPATHLEN]; /* first dir to search
						for data	*/
extern char mbsFimsDataDir[MAXPATHLEN]; /* dir to search for IMS
					     data		*/
extern char mbsDrDataDir[MAXPATHLEN];	/* dir to search for DR
					   data		*/

/****************************************************************************
 *  STATIC vars/funcs
 ***************************************************************************/
static
void verify_mbs_init(void);
static
int find_behavior_flag(char *behavior_label,     /* (I) label to search fo */
		       BEHAVIOR_STRUC **pfound); /* (O) ptr to matching item */
static
int clean_up_behavior_label(char *behavior_label, /* (I) initial label */
			    char *clean_label);   /* (O) cleaned-up label */
static
int find_perl_line(FILE *f,            /* (I) ptr to (open) file */
		   char *perl_lbl,     /* (I) identifier label */
		   long *ndata,         /* (O) # data fragments in line */
		   char **p_data);     /* (O) pointer to start of 1st data frag */


/***************************************************************************
 *  "Behavior control" funcs
 **************************************************************************/

/*****************************************************************
 *  set_behavior
 *  Turns specified behavior "on"
 *  (part of run-time "behavior control" system)
 *****************************************************************/
int set_behavior(char *behavior_label)	    /* (I) some label, e.g. "smooth oas" */
{
    F_INIT("set_behavior");
    BEHAVIOR_STRUC *pfound;

    /* make sure init has been done */
    verify_mbs_init();

    /* search for this label in linked list */
    XONFAIL( find_behavior_flag(behavior_label,&pfound) );

    /* if not already in list, add it to start of list */
    if( pfound IS NULL ) {
	/* allocate new item */
	XONFAIL( mbs_malloc(sizeof(BEHAVIOR_STRUC),(char**) &pfound) );
	/* clean up (copy of) label, and store in item */
	XONFAIL( clean_up_behavior_label(behavior_label,pfound->behavior_label) );

	/* make this item the new start of the list */
	pfound->prev = NULL;
	pfound->next = (void*) p_behavior_start;
	p_behavior_start = pfound;
    }

    F_END;
}

/*****************************************************************
 *  unset_behavior
 *  Turns specified behavior "off"
 *  (part of run-time "behavior control" system)
 *****************************************************************/
int unset_behavior(char *behavior_label)    /* (I) some label, e.g. "smooth oas" */
{
    F_INIT("unset_behavior");
    BEHAVIOR_STRUC *pfound;

    /* make sure init has been done */
    verify_mbs_init();

    /* search for this label in linked list */
    XONFAIL( find_behavior_flag(behavior_label,&pfound) );

    /* if found, remove it */
    if( pfound ) {
	/* connect the two adjacent elements of linked list */
	if( pfound->prev )
	    ((BEHAVIOR_STRUC*) (pfound->prev))->next = pfound->next;
	if( pfound->next )
	    ((BEHAVIOR_STRUC*) (pfound->next))->prev = pfound->prev;
	/* if this item also happens to be first in list, alter start ptr */
	if( pfound IS p_behavior_start )
	    p_behavior_start = ((BEHAVIOR_STRUC*) (pfound->next));
	/* now de-alloc item */
        mbs_free((char*) pfound);
    }

    F_END;
}

/*****************************************************************
 *  is_behavior
 *  BOOLEAN function: returns TRUE if behavior is "on"
 *  (part of run-time "behavior control" system)
 *****************************************************************/
int is_behavior(char *behavior_label)	    /* (I) some label, e.g. "smooth oas" */
{
    F_INIT("is_behavior");
    BEHAVIOR_STRUC *pfound;

    /* make sure init has been done */
    verify_mbs_init();

    /* search for this label in linked list */
    XONFAIL( find_behavior_flag(behavior_label,&pfound) );

    status = SUCCESS;
  done:
    if( status IS SUCCESS )
	return (pfound ISNT NULL);
    else
	return FALSE;
}



/***************************************************************************
 *  File I/O funcs
 **************************************************************************/

/****************************************************************************
 * MbsBuildFileName
 * File name (string) utility: prepends directory name to (raw) file
 * name, if that file name currently lacks a path
 * NOTE: calling code must have allocated fullFileName as string
 * of length MAXPATHLEN
 * NOTE: OK to use same file name string for input & output (raw/fullFileName)
 ****************************************************************************/
int MbsBuildFileName
   (char    *pathName,          /* (I) name of path (dir) to prepend */
    char    *rawFileName,       /* (I) raw name of file */
    char    *fullFileName)      /* (O) full name for file */
{
    F_INIT("MbsBuildFileName");
    char     locPathName[MAXPATHLEN] = "";
    char     locFileName[MAXPATHLEN] = "";
    long     iChar;
    TBoolean foundSlash;

    /* make sure init has been done */
    verify_mbs_init();

    XONTRUE(fullFileName IS NULL, "No ptr to output string" );
    XONTRUE(rawFileName IS NULL ||
            rawFileName[0] IS '\0', "No input file supplied" );

    /* search raw file name for slash char */
    foundSlash = FALSE;
    for(iChar=0; iChar< (long) strlen(rawFileName); iChar++)
    {
        if(rawFileName[iChar] IS MBS_PATH_DELIM)
        {
            foundSlash = TRUE;
        }
    }

    /* In some cases, no prepending needed: no path given, 
     * or file name already has path */
    if(pathName IS NULL ||
       pathName[0] IS '\0' ||
       foundSlash)
    {
        strncpy(locFileName,rawFileName,MAXPATHLEN-1);
    }
    else
    {
        strncpy(locPathName,pathName,MAXPATHLEN-1);
        XONFAIL(MbsAddTrailSlash(locPathName) );
        strncpy(locFileName,locPathName,MAXPATHLEN-1);
        strncat(locFileName,rawFileName,MAXPATHLEN-1-strlen(locFileName));
    }
    /* either way, copy to output string */
    strncpy(fullFileName,locFileName,MAXPATHLEN-1);

    F_END;
}


/****************************************************************************
 * MbsGetDefOutputPath
 * Determines the default output directory 
 * (at present, this is $HOME on UNIX, and c:\ on PC)
 * Calling code must have already allocated output string
 ****************************************************************************/
int MbsGetDefOutputPath
   (long    maxNumChars,        /* (I) max # chars to copy to outputPath */
    char   *outputPath)         /* (O) default output path */
{
    F_INIT("MbsGetDefOutputPath");

    XONTRUE(outputPath IS NULL, "No pointer to output string" );
#ifdef UNIX
    strncpy(outputPath,getenv("HOME"),MAX(0,maxNumChars-1));
#endif
#ifdef DOS
    strncpy(outputPath,DOS_HOME_DIR,MAX(0,maxNumChars-1));    
#endif

    F_END;
}




/****************************************************************************
 *  MbsAddTrailSlash
 *  (DOS or UNIX) Adds trailing fwd (UNIX) or bkwd (DOS) slash to
 *  file/path spec string
 ****************************************************************************/
int MbsAddTrailSlash
   (char *pathSpec)             /* (I/O) path spec to alter */
{
    F_INIT("MbsAddTrailSlash");
    long    sLen;

    /* make sure init has been done */
    verify_mbs_init();

    XONTRUE(pathSpec IS NULL,
            "Failed to supply path name string" );

    /* remove leading/trailing spaces */
    XONFAIL(remove_end_blanks(pathSpec));

    sLen = strlen(pathSpec);
    /* if string is empty or lacks trailing slash, add it */
    if(sLen IS 0 ||
       pathSpec[MAX(0,sLen-1)] ISNT MBS_PATH_DELIM)
    {
        pathSpec[sLen] = MBS_PATH_DELIM;
        pathSpec[sLen+1] = '\0';
    }

    F_END;
}



/****************************************************************************
 *  MbsRemoveTrailSlash
 *  (DOS or UNIX) Removes trailing fwd (UNIX) or bkwd (DOS) slash 
 *  from file/path spec string; 
 *  Can handle case of root path ("\" or "/") in two ways
 ****************************************************************************/
int MbsRemoveTrailSlash
   (char     *pathSpec,             /* (I/O) path spec to alter */
    TBoolean  keepRootPath)         /* (I) if input path is root path
                                     * (i.e., just single slash), then if:
                                     *   TRUE: keep the trailing (only) slash
                                     *   FALSE: remove trailing slash (will
                                     *   then return empty string for path) */
{
    F_INIT("MbsRemoveTrailSlash");
    long    sLen;

    /* make sure init has been done */
    verify_mbs_init();

    XONTRUE(pathSpec IS NULL,
            "Failed to supply path name string" );
    /* remove leading/trailing spaces */
    XONFAIL(remove_end_blanks(pathSpec));
    sLen = strlen(pathSpec);
    /* check if we have just root dir name, and user wants
     * to keep trailing slash in this case */
    if(sLen IS 1 &&
       pathSpec[0] IS MBS_PATH_DELIM &&
       keepRootPath)
    {
        /* in this case, do nothing */
    }
    /* else if string isn't empty and has trailing slash, remove it */
    else if(sLen > 0 &&
            pathSpec[MAX(0,sLen-1)] IS MBS_PATH_DELIM)
    {
        pathSpec[sLen-1] = '\0';
    }

    F_END;
}

/****************************************************************************
 *  mbs_filefound
 *  BOOLEAN func to check if file exists
 *    If rootname contains a path delim, is ordinary "file-find" function
 *    Otherwise, looks in all the "usual" data locations ($DATA/LATEST, etc);
 *  Returns boolean true if found
 ****************************************************************************/
int mbs_filefound(char *rootname)    /* (I) root name (w/o path) of file */
{
    F_INIT("mbs_filefound");
    char fullname[MAXPATHLEN];

    /* make sure init has been done */
    verify_mbs_init();

    /* try to find file (if not found, sets fullname to empty) */
    XONFAIL( mbs_findfile(rootname,MAXPATHLEN-1,fullname) );

    status = SUCCESS;
  done:
    return (fullname[0] ISNT '\0');
}


/****************************************************************************
 *  mbs_fopen
 *  Finds and opens data file for input (w/"smart" file-find funcs)
 *  Returns pointer to file (or NULL, if failed)
 *  NB: if rootname contains slash char, it's assumed to be full pathspec
 ****************************************************************************/
FILE *mbs_fopen(char *rootname,      /* (I) root name of file */
		char *mode)          /* (I) file mode ("r", "a", etc.) */
{
    FILE *file = NULL;
    char fullname[MAXPATHLEN];
    F_INIT("mbs_fopen");

    /* make sure init has been done */
    verify_mbs_init();

    /* try to find file (if not found, sets fullname to empty) */
    XONFAIL( mbs_findfile(rootname,MAXPATHLEN-1,fullname) );
    sprintf(outmesg,"Failed to find file: %s",rootname);
    XONTRUE( fullname[0] IS '\0', outmesg );

    /* open it */
    file = fopen(fullname, mode);

    status = SUCCESS;
  done:
    if( status IS SUCCESS )
	return file;
    else
	return NULL;
}


/****************************************************************************
 *  mbs_findfile
 *  Fundamental function to find data file
 *    If rootname contains a path delim, is ordinary "file-find" function
 *    Otherwise, looks in all the "usual" data locations ($DRDATA, etc);
 *  Sets fullname to full name of file (if found), or empty string (if not found)
 *  "Usual" dirs: user-defined [if set], DATA/LATEST, DR_DATA, 
 *                /opt/fims/1.0.0.0/data/LATEST
 *  NB: function itself does not fail (return FAILURE) if file not found;
 *  only fails if some fatal error occurs
 ****************************************************************************/
int mbs_findfile(char *rootname,       /* (I) root name of file */
		 long maxnamelen,       /* (I) max # chars to copy to fullname */
		 char *fullname)       /* (O) (if found) full name of file */
{
    F_INIT("mbs_findfile");
    char tmpname[MAXPATHLEN];
    FILE *file = NULL;
    long ipos;

    /* make sure init has been done */
    verify_mbs_init();

    /* reset return arg */
    fullname[0] = '\0';

    /* special case: if rootname contains path delim, 
     * we assume it's a full pathspec */
    XONFAIL( str_pos(rootname,MBS_PATH_DELIM,&ipos) );
    if ( ipos >= 0 ) 
    {
        file = fopen(rootname, "r");
        if ( file ) 
        {
            fclose(file);
            strncpy(fullname,rootname,maxnamelen-1);
        }
        else 
        {
            fullname[0] = '\0';
        }
        status = SUCCESS;
	goto done;
    }

    /* First: the user-defined primary (i.e., first-searched) data dir */
    if ( (!file) && mbsPrimaryDataDir[0] ISNT '\0' ) 
    {
        XONFAIL(MbsBuildFileName(mbsPrimaryDataDir, rootname, tmpname) );
        file = fopen(tmpname, "r");
    }

    /* DR "DRDATA" dir */
    if ( (!file) && mbsDrDataDir[0] ISNT '\0' ) 
    {
        XONFAIL(MbsBuildFileName(mbsDrDataDir, rootname, tmpname) );
        file = fopen(tmpname, "r");
    }

#ifdef DOS
    /* Hard-coded PC DRDATA location (when DRDATA not defined) */
    if (!file ) 
    {
        XONFAIL(MbsBuildFileName(PC_DRDATADIR, rootname, tmpname) );
        file = fopen(tmpname, "r");
    }
#endif

#ifdef UNIX_DONT_DO_IT
    /* FIMS "DATA/LATEST" directory */
    if ( (!file) && mbsFimsDataDir[0] ISNT '\0' ) 
    {
        XONFAIL(MbsBuildFileName(mbsFimsDataDir, rootname, tmpname) );
        file = fopen(tmpname, "r");
    }
    /* hard-coded trading FIMS data dir */
    if (!file ) 
    {
        XONFAIL(MbsBuildFileName(TRADING_FIMSDATADIR, rootname, tmpname) );
        file = fopen(tmpname, "r");
    }
    /* alternate hard-coded trading FIMS data dir */
    if (!file ) 
    {
        XONFAIL(MbsBuildFileName(TRADING_FIMSDATADIR2, rootname, tmpname) );
        file = fopen(tmpname, "r");
    }
#endif

    /* as final possibility, look in current directory */
    if (!file ) 
    {
	strcpy(tmpname,rootname);
	file = fopen(tmpname,"r");
    }

    /* Found? */
    if ( file )
        strncpy(fullname,tmpname,maxnamelen-1);
    else
        fullname[0] = '\0';

    status = SUCCESS;
  done:
    if ( file )
        fclose(file);
    return status;
}




/***************************************************************
 * mbs_delfile
 * Returns SUCCESS/FAILURE
 ***************************************************************/
int mbs_delfile(char* fname)      /* (I) name of file to delete	*/
{
    F_INIT("mbs_delfile");

    /* make sure init has been done */
    verify_mbs_init();

    /* returns error code 0=success, non-zero=failure */
    XONFAIL( remove(fname) );	

    F_END;
}


/***************************************************************************
 *  Arith. funcs
 ***************************************************************************/

/**************************************************
 *  mbs_pow
 *  Faster version of pow() (uses exp/ln)
 **************************************************/
double mbs_pow(double x, double y)
{
    return exp(y*log(x));
}

/**************************************************
 *  mbs_sgn
 *  Returns -1 for x<0, 0 for x=0, 1 for x>0
 **************************************************/
double mbs_sgn(double x)
{
    if ( x < 0. )
	return -1.;
    else if ( x > 0. )
	return 1.;
    else
	return 0.;
}

/**************************************************
 *  mbs_round
 *  Returns rounded value 
 **************************************************/
long mbs_round(double x)
{
  long i;

  i = (long) floor(x);
  if ( x-i >= 0.5 )
      i++;
  return i;
}

/**************************************************
 *  mbs_lininterp1
 *  1-d linear interpolation (assumes ascending xarr)
 *  error if not in ascending order (3/28/96 - DL) 
 *  Returns boolean true if successful
 *  NB: does NOT extrapolate
 **************************************************/
int mbs_lininterp1(long n,	      /* (I) # pts */
		   double xarr[],     /* (I) array of indep var */
		   double yarr[],     /* (I) array of dep var */
		   double x,	      /* (I) value at which to interpolate */
		   double *y)	      /* (O) interpolated value */
{
    F_INIT("mbs_lininterp1");
    long i,found;

    /* make sure init has been done */
    verify_mbs_init();

    *y = 0.;

    XONTRUE( n < 2, "Not enough pts to interp" );

    for ( i = 0; i < n-1; i++ ) {
	XONTRUE( xarr[i] > xarr[i+1], "xarr not in ascending order" );
    }

    for ( i = 0; i < n-1; i++ ) {
	XONTRUE( xarr[i] IS xarr[i+1], "adjacent x numbers should not be same");
    }

    if ( x <= xarr[0] ) {
	*y = yarr[0];
    }
    else if ( x >= xarr[n-1] ) {
	*y = yarr[n-1];
    }
    else {
	found = FALSE;
	for ( i = 0; i < n-1; i++ ) 
        {
	    if ( xarr[i] <= x && x <= xarr[i+1] ) 
            {
		*y = yarr[i] + (yarr[i+1] - yarr[i]) /
			(xarr[i+1] - xarr[i]) * (x - xarr[i]);
		found = TRUE;
	    }
	}
        XONTRUE( !found, "Failed to interpolate" );
    }

    F_END;
}

/**************************************************
 *  mbs_lininterp2
 *  1-d linear interpolation (assumes ascending xarr)
 *  error if not in ascending order (3/28/96 - DL) 
 *  Returns boolean true if successful
 *  NB: does extrapolate
 **************************************************/
int mbs_lininterp2(long n,	      /* (I) # pts */
		   double xarr[],     /* (I) array of indep var */
		   double yarr[],     /* (I) array of dep var */
		   double x,	      /* (I) value at which to interpolate */
		   double *y)	      /* (O) interpolated value */
{
    F_INIT("mbs_lininterp2");
    long i,found;

    /* make sure init has been done */
    verify_mbs_init();

    *y = 0.;

    XONTRUE( n < 2, "Not enough pts to interp" );

    for ( i = 0; i < n-1; i++ ) {
	XONTRUE( xarr[i] > xarr[i+1], "xarr not in ascending order" );
    }

    for ( i = 0; i < n-1; i++ ) {
	XONTRUE( xarr[i] IS xarr[i+1], "adjacent x numbers should not be same");
    }

    if ( x < xarr[0] ) {
	*y = yarr[0] + (yarr[i] - yarr[0]) /
		(xarr[1] - xarr[0]) * (x - xarr[0]);
    }
    else if ( x IS xarr[0] ) {
	*y = yarr[0];
    }
    else if ( x IS xarr[n-1] ) {
	*y = yarr[n-1];
    }
    else if ( x > xarr[n-1] ) {
	*y = yarr[n-2] + (yarr[n-1] - yarr[n-2]) /
		(xarr[n-1] - xarr[n-2]) * (x - xarr[n-2]);
    }
    else {
	found = FALSE;
	for ( i = 0; i < n-1; i++ ) {
	    if ( xarr[i] <= x && x <= xarr[i+1] ) {
		*y = yarr[i] + (yarr[i+1] - yarr[i]) /
			(xarr[i+1] - xarr[i]) * (x - xarr[i]);
		found = TRUE;
	    }
	    XONTRUE( !found, "Failed to interpolate" );
	}
    }

    F_END;
}


/***************************************************************************
 *  Date funcs (supplementing GTO funcs)
 ***************************************************************************/

/************************************************************
 *  GetDayCountStr
 *  Converts GTO day count integer into string
 *  @@ someday, find GTO macro/func that does this
 *  and use it inside this func (couldn't find it today...)
 ************************************************************/
int EXPORT GetDayCountStr
   (long      dayCount,         /* (I) day count int */
    char     *dayCountStr)      /* (O) string rep */
{
    switch (dayCount)
    {
      case  GTO_ACT_365:
        strcpy(dayCountStr,LOC_GTO_ACT_365_STR);
        break;
      case  GTO_ACT_365F:
        strcpy(dayCountStr,LOC_GTO_ACT_365F_STR);
        break;
      case  GTO_ACT_360:
        strcpy(dayCountStr,LOC_GTO_ACT_360_STR);
        break;
      case  GTO_B30_360:
        strcpy(dayCountStr,LOC_GTO_B30_360_STR);
        break;
      case  GTO_B30E_360:
        strcpy(dayCountStr,LOC_GTO_B30E_360_STR);
        break;
      case  GTO_ACT_365FJ:
        strcpy(dayCountStr,LOC_GTO_ACT_365FJ_STR);
        break;
      case  GTO_B30E_360I:
        strcpy(dayCountStr,LOC_GTO_B30E_360I_STR);
        break;
      case  GTO_B30_360_FIXED:
        strcpy(dayCountStr,LOC_GTO_B30_360F_STR);
        break;
      case  GTO_B30EP_360:
        strcpy(dayCountStr,LOC_GTO_B30EP_360_STR);
        break;
      default:
        strcpy(dayCountStr,"Unknown day count type");
        break;
    }
    return(SUCCESS);
}


/************************************************************
 * SimplifyDL()
 * Modifies date list (may actually shorten list):
 * sorts in order of increasing date,
 * and removes any repeated or zero dates 
 ************************************************************/
int EXPORT SimplifyDL
   (TDateList **DL)   /* (I/O) ptr to ptr for date list */
{
    F_INIT("SimplifyDL");
    TBoolean stillSorting;
    long dtIdx;

    TDateList *newDL = NULL;
    long idx;
    long idxNew;
    long nValid;
    
    /* check inputs */
    XONTRUE( DL IS NULL || *DL IS NULL,
            "Failed to supply ptr to date list");

    /* Sort the list */
    do {
        stillSorting = FALSE;
        for(dtIdx=0; dtIdx<(*DL)->fNumItems-1; dtIdx++)
            if( (*DL)->fArray[dtIdx] > (*DL)->fArray[dtIdx+1] )
            {
                SwitchDates(&((*DL)->fArray[dtIdx]),
                          &((*DL)->fArray[dtIdx+1]));
                stillSorting = TRUE;
            }
    } while(stillSorting);

    /* Zero out any repeated dates */
    for(dtIdx=0; dtIdx<(*DL)->fNumItems-1; dtIdx++)
        if( (*DL)->fArray[dtIdx] IS (*DL)->fArray[dtIdx+1] )
            (*DL)->fArray[dtIdx] = 0;

    /* Count # valid dates (non-zero) */
    nValid = 0;
    for(dtIdx=0; dtIdx<(*DL)->fNumItems; dtIdx++)
        if( (*DL)->fArray[dtIdx] ISNT 0 )
            nValid++;
    XONTRUE( nValid IS 0,
            "Date list contained no valid dates" );

    /* If # valid dates fewer than total # dates,
     * create new, shorter date list
     */
    if( nValid < (*DL)->fNumItems )
    {
        newDL = GtoNewEmptyDateList(nValid);
        idxNew = 0;
        for(idx=0; idx < (*DL)->fNumItems; idx++)
        {
            if( (*DL)->fArray[idx] ISNT 0 )
            {
                newDL->fArray[idxNew] = (*DL)->fArray[idx];
                idxNew++;
            }
        }
        /* de-alloc old list */
        GtoFreeDateList(*DL);
        /* set output pointer to new list */
        *DL = newDL;
    }

    status = SUCCESS;
  done:
    if( status IS FAILURE )
        GtoFreeDateList(newDL);
    return (status);
}


/************************************************************
 * Simply exchanges the values of the two date vars
 ************************************************************/
int EXPORT SwitchDates
   (TDate *date1,               /* (I/O) */
    TDate *date2)               /* (I/O) */
{
    TDate tmpDate;
    
    tmpDate = *date1;
    *date1 = *date2;
    *date2 = tmpDate;
    return(SUCCESS);
}


/***************************************************************************
 * MbsWala
 * Utility to compute WALA of MBS, using the 
 * (time-independent) eff orig term (WARM+WALA) and maturity date
 * Note: this function can and should return a WALA < 0
 * in cases where the MBS is to be issued in future
 ***************************************************************************/
int EXPORT MbsWala
   (TDate   currDate,             /* (I) date for which WALA needed */
    TDate   mbsMatDate,           /* (I) maturity of MBS */
    long    mbsEffOrigTerm,       /* (I) WARM+WALA of mbs, in months */
    long   *wala)                 /* (O) WALA, in months */
{
    F_INIT("MbsWala");
    long   warm;

    /* compute WARM */
    XONFAIL( MbsWarm(currDate,mbsMatDate,&warm) );
    *wala = mbsEffOrigTerm - warm;

    F_END;
}


/***************************************************************************
 * MbsWarm
 * Utility to compute WARM of MBS, using the 
 * (time-independent) maturity date
 * Note: this function can and should return a WARM > MBS term
 * in cases where the MBS is to be issued in future
 ***************************************************************************/
int EXPORT MbsWarm
   (TDate   currDate,             /* (I) date for which WARM needed */
    TDate   mbsMatDate,           /* (I) maturity of MBS */
    long   *warm)                 /* (O) WARM, in months */
{
    F_INIT("MbsWarm");

    /* # months to maturity */
    XONFAIL( MonsDiff(currDate, mbsMatDate, warm) );

    F_END;
}


/***************************************************************************
 *  GetTime
 *  Copies string containing current time
 ***************************************************************************/
int EXPORT GetTime
   (char   *currTime)           /* (O) string w/current time */
{
    time_t  tNow;

    /* I think this works on both UNIX and DOS */
    time(&tNow);
    strcpy(currTime, ctime(&tNow));
    /* trim trailing CR-LF at end of ctime() string */
    currTime[strlen(currTime) - 1] = '\0';  
    return(SUCCESS);
}


/***********************************************************
 *  GetHolidayFileName
 *  Determines location/name of holiday file
 ***********************************************************/
int GetHolidayFileName
   (char   *holidayFile)        /* (O) full name (w/path) of holiday file */
{
    F_INIT("GetHolidayFileName");

    /* Init outputs */
    holidayFile[0] = '\0';

    /* Find our std holiday file */
    XONFAIL( mbs_findfile(HOLIDAY_FILE,MAXPATHLEN-1,holidayFile) );
    /* If not found, use GTO "NONE" (i.e., weekends only) */
    if( holidayFile[0] IS '\0' )
    {
        strcpy(holidayFile,"NONE");
    }

    F_END;
}



/***********************************************************
 *  make_TDate
 *  Utility to set TDate to explicit month/day/year
 ***********************************************************/
int make_TDate(long month,	/* (I) 1-12 */
	       long day,		/* (I) 1-31 */
	       long year,	/* (I) should be 4-digit yyyy form */
	       TDate *date)	/* (O) set to corresponding TDate */
{
    F_INIT("make_TDate");
    TMonthDayYear mdy;

    XONTRUE( year < 100, "Must supply 4-digit year number" );
    mdy.month = month;
    mdy.day = day;
    mdy.year = year;
    GtoMDYToDate(&mdy,date);

    F_END;
}


/***********************************************************
 *  TDateOf
 *  Converts YYYYMMDD date to GTO TDate
 ***********************************************************/
int TDateOf(long idate,         /* (I) YYYYMMDD date */
	    TDate *Date)        /* (O) TDate */
{
    int status = FAILURE;
    TMonthDayYear tmdy;

    tmdy.day = idate % 100;
    tmdy.month = (idate/100) % 100;
    tmdy.year = idate/10000;
    if( GtoMDYToDate(&tmdy, Date) ISNT SUCCESS )
        goto done;

    status = SUCCESS;
  done:
    return status;
}


/***********************************************************
 *  IDateOf
 *  Converts GTO TDate to YYYYMMDD date
 ***********************************************************/
int IDateOf(TDate date,
	    long *iDate)
{
    int status = FAILURE;
    TMonthDayYear mdy;

    if( GtoDateToMDY(date,&mdy) ISNT SUCCESS)
        goto done;
    *iDate = mdy.day+100*(mdy.month+100*mdy.year);
    
    status = SUCCESS;
  done:
    return status;
}


/***************************************************************************
 *  FindDateInDatelist()
 *  Determines if specified date is in date list, and if so,
 *  returns index in date list
 ***************************************************************************/
int FindDateInDatelist
   (TDate       date,           /* (I) date to find in date list */
    TDateList  *dateList,       /* (I) datelist to search */
    TBoolean   *dateFound,      /* (O) TRUE if date found */
    long       *dlIdx           /* (O) (if found) index of this date
                                 * in datelist; -1 if not found */
    )
{
    int status = FAILURE;
    static char routine[] = "FindDateInDatelist";
    long idx;

    /* reset outputs */
    *dateFound = FALSE;
    *dlIdx = -1;

    /* check inputs */
    if( dateList IS NULL )
    {
        GtoErrMsg("%s: Did not supply datelist\n",routine);
        goto done;
    }
    
    for(idx=0; idx<dateList->fNumItems; idx++)
        if( dateList->fArray[idx] IS date )
        {
            *dateFound = TRUE;
            *dlIdx = idx;
        }

    status = SUCCESS;
  done:
    return status;
}


/************************************************************************
 * AdvDate()
 * Advances given date by nIntervals
 ************************************************************************/
int AdvDate(TDate *date,                  /* (I/O) date to alter */
            long nIntervals,               /* (I) # intervals */
            TDateInterval *dateInterval)  /* (I) (single) interval */
{
    int status = FAILURE;
    TDate tmpDate;		       /* for date arith */
    TDateInterval totInterval;

    /* copy interval */
    totInterval = *dateInterval;
    totInterval.prd *= nIntervals;

    /* advance the date */
    tmpDate = *date;
    if (GtoDtFwdAny(tmpDate, &totInterval, date) ISNT SUCCESS)
	goto done;

    status = SUCCESS;
  done:
    return (status);
}


/************************************************************************
 * NxtMth()
 * Advances date by nMonths months;
 * Note: OK to use on same date (e.g., NxtMth(date,&date))
 * Note: if day-of-month beyond last valid day in new month,
 * function sets day-of-month to last valid day;
 * Note also that function tries to preserve the day-of-month
 * of startDate--if we start at Jan 31, and advance 2 months
 * (i.e., past Feb), the result would be Mar 31
 ************************************************************************/
int NxtMth(TDate startDate,
           long nMonths,
           TDate * endDate
           )
{
    int status = FAILURE;
    TDateInterval monthInterval;	       /* 1M */

    /* create a month interval */
    if (GtoMakeDateInterval(1, 'M', &monthInterval) ISNT SUCCESS)
	goto done;
    /* advance the date */
    *endDate = startDate;
    if( AdvDate(endDate,nMonths,&monthInterval) ISNT SUCCESS)
	goto done;

    status = SUCCESS;
  done:
    return (status);
}

/****************************************************************************
 * IDateNxtMth
 * Computes new YYYYMMDD date nMonths fwd/back from start date
 ****************************************************************************/
int EXPORT IDateNxtMth
   (long     startIDate,        /* (I) start date (YYYYMMDD) */
    long     nMonths,           /* (I) # months to advance (or bkwd) */
    TBoolean preserveEndOfMon,  /* (I) TRUE: if start date is last day of mon
                                 * (e.g., 2/28), then ensure that newIDate
                                 * is also last day of new month (e.g., 4/31);
                                 * else, if FALSE, always preserve day-of-month*/
    long    *newIDate)          /* (O) new date (YYYYMMDD)  */
{
    F_INIT("IDateNxtMth");
    long     startMon;
    long     nxtMon;
    TDate    startDate;
    TDate    newDate;
    TBoolean forceToLastDayOfMonth;


    /* convert to TDate */
    XONFAIL(TDateOf(startIDate,&startDate));
    /* determine month of start date, and day after */
    XONFAIL(GetMOY(startDate,&startMon));
    XONFAIL(GetMOY(startDate+1,&nxtMon));
    /* will we need to force resulting day to last day of new month? */
    forceToLastDayOfMonth = (preserveEndOfMon && (startMon ISNT nxtMon));
    
    /* advance to new date */
    XONFAIL(NxtMth(startDate,nMonths,&newDate));
    /* if needed, force day-of-month to last day of new month */
    if(forceToLastDayOfMonth)
    {
        /* go to 1st of next month, then back 1 day */
        XONFAIL(NxtMth(newDate,1,&newDate));
        XONFAIL(SetDOM(1,FALSE,&newDate));
        newDate--;
    }

    /* convert back to YYYYMMDD */
    XONFAIL(IDateOf(newDate,newIDate));

    F_END;
}

/************************************************************************
 * MonthNumOfTDate
 * Computes "monthnum" from TDate: 12*year+month
 * (useful in some monthly date arithmetic)
 ************************************************************************/
int MonthNumOfTDate(TDate    date,     /* (I) input date */
                    long    *monthNum) /* (O) monthnum: 12*yr+month */
{
    F_INIT("MonthNumOfTDate");
    TMonthDayYear mdy;

    /* reset output */
    *monthNum = 0;

    XONFAIL( GtoDateToMDY(date, &mdy) );
    *monthNum = 12*mdy.year + (mdy.month-1);

    F_END;
}


/************************************************************************
 * TDateOfMonthNum
 * Computes TDate from "monthnum" (12*year+month),
 * where lost day-of-month information is supplied
 * as additional input
 ************************************************************************/
int TDateOfMonthNum(long     monthNum, /* (I) monthnum (12*yr+mon) */
                    long     dom,      /* (I) day-of-mon (1-31) to
                                        * which to set output date */
                    TDate   *date)     /* (O) output date */
{
    F_INIT("TDateOfMonthNum");
    TMonthDayYear mdy;

    /* reset output */
    *date = 0;

    mdy.day = dom;
    mdy.month = 1 + (monthNum % 12);
    mdy.year = monthNum/12;
    XONFAIL( GtoMDYToDate(&mdy, date) );

    F_END;
}


/************************************************************************
 * MonsDiff()
 * Computes the number of months between two dates (date2 - date1)
 * (day-of-month for dates ignored)
 ************************************************************************/
int MonsDiff(TDate date1,       /* (I) date */
             TDate date2,       /* (I) date */
             long *nmons)       /* (O) difference, in months */
{
    F_INIT("MonsDiff");
    long iMonNum1;
    long iMonNum2;

    /* reset output */
    *nmons = 0;

    /* compute 12*year + mon for each date */
    XONFAIL( MonthNumOfTDate(date1,&iMonNum1) );
    XONFAIL( MonthNumOfTDate(date2,&iMonNum2) );
    /* compute difference in months */
    *nmons = iMonNum2 - iMonNum1;

    F_END;
}


/************************************************************************
 * GetDOM()
 * Returns the day-of-month of a TDate
 ************************************************************************/
int GetDOM(TDate  date,         /* (I) date */
           long  *dom)          /* (O) day-of-month (1-31)*/
{
    F_INIT("GetDOM");
    TMonthDayYear tmpMdy;

    /* decompose date */
    XONFAIL(GtoDateToMDY(date, &tmpMdy) );
    *dom = tmpMdy.day;

    F_END;
}


/************************************************************************
 * SetDOM()
 * Sets the day-of-month of a TDate to specified value
 ************************************************************************/
int SetDOM(long dom,            /* (I) day-of-month (1-31)*/
           TBoolean endOfMonth, /* (I) determines handling of a day-of-month
                                 * later than last day in mon (e.g., 
                                 * dom=31 for February):
                                 * TRUE=set dom to last valid day,
                                 * FALSE=declare error */
           TDate * date)        /* (I/O) date to alter */
{
    F_INIT("SetDOM");
    TMonthDayYear tmpMdy;
    TDate saveDate;

    /* remember date */
    saveDate = *date;
    /* check dom */
    XONTRUE(dom < 1 || dom > 31, "Invalid day-of-month" );
    /* decompose date */
    XONTRUE( GtoDateToMDY(*date, &tmpMdy) ISNT SUCCESS,
            "Failed to decompose TDate into m,d,y--bad date" );
    /* alter day-of-month */
    tmpMdy.day = dom;
    /* try to convert m/d/y back to TDate */
    if( GtoMDYToDate(&tmpMdy, date) ISNT SUCCESS )
    {
        /* if failure is not invalid-dom, then quit now */
        XONTRUE( dom <= 28, "Failed to set day-of-mon" );
        /* otherwise, then clearly this day-of-mon
         * is beyond end of this month */
        if( endOfMonth )
        {
            /* if user requests it, go to last valid day:
             * go to 1st of next month, then back 1 day */
            *date = saveDate;
            XONFAIL( NxtMth(*date,1,date) );
            XONFAIL( SetDOM(1,TRUE,date) );
            (*date)--;
        }
        else 
        {
            /* otherwise, declare an error */
            sprintf(outmesg,
                    "Cannot set day-of-mon of %s to %ld",
                    GtoFormatDate(saveDate),
                    dom);
            XONTRUE( TRUE, outmesg );
        }

    }

    F_END;
}


/************************************************************************
 * ForceToClosestDOM()
 * Forces the date to the nearest date having the desired day-of-month
 ************************************************************************/
int ForceToClosestDOM
   (long desiredDom,            /* (I) desired day-of-month (1-31)*/
    TBoolean endOfMonth,        /* (I) determines handling of a day-of-month
                                 * later than last day in mon (e.g., 
                                 * dom=31 for February):
                                 * TRUE=set dom to last valid day,
                                 * FALSE=declare error */
    TDate * date)               /* (I/O) date to alter */
{
    F_INIT("ForceToClosestDOM");
    long   dateDom;
    TDate  currMonDate;
    TDate  nextMonDate;
    TDate  prevMonDate;

    /* date in same month with this d-o-m */
    currMonDate = *date;
    XONFAIL( SetDOM(desiredDom,TRUE,&currMonDate) );

    XONFAIL( GetDOM(*date,&dateDom) );
    if( dateDom IS desiredDom )
    {
        /* done--no change */
    }
    else if( dateDom > desiredDom )
    {
        /* compare w/date in next month */
        XONFAIL( NxtMth(currMonDate,1,&nextMonDate) );
        XONFAIL( SetDOM(desiredDom,TRUE,&nextMonDate) );
        if( nextMonDate-(*date) > (*date)-currMonDate )
        {
            *date = currMonDate;
        }
        else
        {
            *date = nextMonDate;
        }
    }
    else
    {
        /* compare w/date in prev month */
        XONFAIL( NxtMth(currMonDate,-1,&prevMonDate) );
        XONFAIL( SetDOM(desiredDom,TRUE,&prevMonDate) );
        if( currMonDate-(*date) > (*date)-prevMonDate )
        {
            *date = prevMonDate;
        }
        else
        {
            *date = currMonDate;
        }
    }


    F_END;
}



/************************************************************************
 * GetMOY()
 * Returns the month-of-year of a TDate
 ************************************************************************/
int GetMOY(TDate  date,         /* (I) date */
           long  *moy)          /* (O) month-of-year (1-12)*/
{
    F_INIT("GetMOY");
    TMonthDayYear tmpMdy;

    /* decompose date */
    XONFAIL(GtoDateToMDY(date, &tmpMdy) );
    *moy = tmpMdy.month;

    F_END;
}




/***************************************************************************
 *  Misc utility functions
 ***************************************************************************/


/***************************************************************************
 *  CopyDouble
 *  Copies array of doubles from one array to another
 ***************************************************************************/
int EXPORT CopyDouble
   (long      nValues,          /* (I) # values to copy */
    double   *inpArray,         /* (I) input array */
    double   *outArray)         /* (O) output array */
{
    F_INIT("CopyDouble");
    long   idx;

    XONTRUE( inpArray IS NULL, "No input array specified" );
    XONTRUE( outArray IS NULL, "No output array specified" );
    for(idx=0; idx<nValues; idx++)
    {
        outArray[idx] = inpArray[idx];
    }

    F_END;
}



/***************************************************************************
 *  Misc GTO-interface utility functions
 ***************************************************************************/


/***************************************************************************
 *  SetSingleFloatArray
 *  Constructor for TFloatRateArray to define a single
 *  floating rate (i.e., weight=1.0)
 ***************************************************************************/
int EXPORT SetSingleFloatArray
   (long     cpnsPerYear,           /* (I) cpns per year (0=zro cpn rate) */
    long     matInMonths,           /* (I) maturity in months */
    long     dayCountConv,          /* (I) day count convention */
    long     curveIndex,            /* (I) which curve (e.g.,
                                       GTO_CURVE_DISCOUNT) */
    long     numSettleDays,         /* (I) (aka spotOffsetDays) # bus. days
                                     * between trade & settle for this rate;
                                     * equiv. to # bus. days from observation of 
                                     * rate to day on which it is true spot */
    TFloatRateArray **floatIndex)   /* (O) ptr to new TFloatRateArray */
{
    F_INIT("SetSingleFloatArray");
    TFloatRateArray *fltIdx = NULL;
    double     defaultWgt = 1.0;
    double     defaultSprd = 0.;
    long       cpnIntInMonths;
	TDateAdjIntvl adjInterval;

    /* check inputs */
    XONRANGE(cpnsPerYear,0,12,"num cpns/year");
    XONTRUE( matInMonths <= 0, "Maturity must be > 0 months" )

    /* call simple GTO constructor */
    XONTRUE( (fltIdx = GtoFloatRateArrayMakeSimple
              (&cpnsPerYear,    /* payments/yr */
               &matInMonths,    /* mat in months */
               &dayCountConv,
               &defaultWgt,     /* weight */
               1)) IS NULL,     /* # index rates */
            "Failed to alloc TFloatRateArray");

    /* for simple rates, we also force the cpn interval 
     * to be no longer than the maturity interval (even if 
     * GTO allows this), since we need this to be true
     * for some tree pricing code 
     */
    /* for zero coupon, always force pay pd to mat pd */
    if( cpnsPerYear IS 0 )
    {
        fltIdx->defs[0].payInterval =
            fltIdx->defs[0].matInterval;
    }
    /* for cpn-bearing, check if pay pd > mat pd */
    else
    {
        cpnIntInMonths = 12 / cpnsPerYear;
        if( cpnIntInMonths > matInMonths )
        {
            fltIdx->defs[0].payInterval =
                fltIdx->defs[0].matInterval;
        }
    }
    /* also set hidden member data */
	GTO_SET_ADJ_INTERVAL_DAYS (adjInterval, numSettleDays)
    fltIdx->defs[0].spotOffset = adjInterval;
//    fltIdx->defs[0].spotOffset.numDays = numSettleDays;
    fltIdx->defs[0].spread = defaultSprd;
    fltIdx->curveIndices[0] = curveIndex;
#ifndef CLIBVER_7_4
    fltIdx->rateInfo = NULL;
#endif

    *floatIndex = fltIdx;
    status = SUCCESS;
  done:
    if( status ISNT SUCCESS ) 
    {
        GtoFloatRateArrayFree(fltIdx);
        *floatIndex = NULL;
    }
    return (status);
}



/***************************************************************************
 *  Rate funcs
 ***************************************************************************/

/********************************************************************************
 *  ResetArmCpn
 *  Function to reset ARM cpn (one reset)
 ********************************************************************************/
int EXPORT ResetArmCpn
   (double    uncappedNewCpn,     /* (I) rate+sprd, before caps */
    double    prevCpn,            /* (I) previous cpn */
    double    pdCapSprd,          /* (I) sticky cap spread */
    double    pdFlrSprd,          /* (I) sticky floor spread */
    double    lifeCap,            /* (I) life cap on cpn */
    double    lifeFlr,            /* (I) life floor on cpn */
    double    roundingMultiple,   /* (I) rounding multiple (e.g., 0.00125) */
    TBoolean  useSmoothing,       /* (I) TRUE to use smoothing of MAX/MIN,
                                   * FALSE to use ordinary (discontin) MAX/MIN */
    double   *maxDiffs,           /* (I) Max diff from other nodes */
    double   *newCpn)             /* (O) new (capped/rounded) cpn */
{
    double   tmpCpn;
    double   tmpCpn2;
    double   capStrike;
    double   flrStrike;

    tmpCpn = uncappedNewCpn;
    if(!IS_ALMOST_ZERO(roundingMultiple))
    {
        tmpCpn = roundingMultiple*ROUND(tmpCpn/roundingMultiple);
    }
    capStrike = MIN(lifeCap, prevCpn+pdCapSprd);
    flrStrike = MAX(lifeFlr, MAX(0.,prevCpn-fabs(pdFlrSprd)));

    if(useSmoothing)
    {
        /* Apply cap/floor to new ARM cpn:
         *   cap: smoothed form of MIN(newNetCpn,capstrike),
         *   flr: smoothed form of MAX(newNetCpn,flrStrike) */
        GtoMinSmooth(tmpCpn, *maxDiffs, capStrike, &tmpCpn2);
        GtoMaxSmooth(tmpCpn2,*maxDiffs, flrStrike, newCpn);
    }
    else
    {
        tmpCpn2 = MIN(tmpCpn,capStrike);
        *newCpn = MAX(tmpCpn2,flrStrike);
    }
    return(SUCCESS);
}


/***********************************************************************
 * GetIndexRate
 * Given rate definition, reset date and zero curves, 
 * computes index rate (incl wgt/spread, if any)
 ***********************************************************************/
int EXPORT GetIndexRate
   (TDate             effResetDate,     /* (I) eff. reset date */
    TFloatRateArray  *indexRateDefn,    /* (I) defn of index rate */
    TCurve           *discZeroCurve,    /* (I) disc zero curve */
    TCurve          **indexZeroCurves,  /* (I) index zero curves */
    char             *holidayFile,      /* (I) holiday file */
    long              interpType,       /* (I) interp for zero curve;
                                         * e.g., GTO_LINEAR_INTERP */
    long              stubType,         /* (I) for coupon-bearing;
                                         * e.g., GTO_STUB_BOND */
    TBoolean          stubAtEnd,        /* (I) for coupon-bearing;
                                         * eg., FALSE */
    double           *indexRate)        /* (O) index rate (w/wgt,sprd) */
{
    F_INIT("GetIndexRate");
    TDate     lkupDate;
    TDate     firstCpnDate;
    TDate     matDate;
    long      iRate;
    double    fundRate;

    /* for each component of index rate defn */
    *indexRate = 0.;
    for(iRate=0; iRate<indexRateDefn->numRates; iRate++)
    {
        /* get lookup date */
        XONFAIL(GtoDateFromBusDaysOffset
                (effResetDate,
                 -indexRateDefn->defs[iRate].spotOffset.interval.prd,
                 holidayFile,
                 &lkupDate) );
        /* first cpn date */
        XONFAIL( GtoDtFwdAny(lkupDate,
                             &(indexRateDefn->defs[iRate].payInterval),
                             &firstCpnDate) );
        /* maturity date */
        XONFAIL( GtoDtFwdAny(lkupDate,
                             &(indexRateDefn->defs[iRate].matInterval),
                             &matDate) );
        /* get fundmental rate */
        if(firstCpnDate >= matDate )    /* simple rate */
        {
            XONFAIL( GtoZerosToSimplePoint
                    (((indexRateDefn->curveIndices[iRate]
                       IS GTO_CURVE_DISCOUNT) ?
                      discZeroCurve : 
                      indexZeroCurves[indexRateDefn->curveIndices[iRate]-1]),
                     interpType,
                     lkupDate,
                     matDate,
                     indexRateDefn->defs[iRate].dayCountConv,
                     &fundRate) );
        }
        else                            /* cpn-bearing rate */
        {
            XONFAIL( GtoZerosToCouponsPoint
                    (((indexRateDefn->curveIndices[iRate] 
                       IS GTO_CURVE_DISCOUNT) ?
                      discZeroCurve : 
                      indexZeroCurves[indexRateDefn->curveIndices[iRate]-1]),
                     interpType,
                     lkupDate,
                     &(indexRateDefn->defs[iRate].payInterval),
                     matDate,
                     indexRateDefn->defs[iRate].dayCountConv,
                     stubType,
                     stubAtEnd,
                     &fundRate) );
        }
        /* apply wgt/spread, & add to total */
        *indexRate += indexRateDefn->defs[iRate].weight * fundRate
            + indexRateDefn->defs[iRate].spread;
    }

    F_END;
}


/***************************************************************************
 *  Misc funcs
 ***************************************************************************/


/***************************************************************************
 *  MbsCopyDArray
 *  Creates copy of double array
 *  Sets copyArray to NULL if no values or orig is NULL
 *  Special behavior: if origArray is NULL, but nValues > 0
 *  (user wants array copy, but no input array exists),
 *  create output array of size nValues and fill w/def value
 ***************************************************************************/
int EXPORT MbsCopyDArray
   (double  *origArray,         /* (I) original array to copy from */
    long     nValues,           /* (I) # elements to fill in copyArray */
    double   defValue,          /* (I) default value (used if origArray=NULL)*/
    double **copyArray)         /* (O) copy */
{
    F_INIT("MbsCopyDArray");
    long   idx;

    if(nValues <= 0)
    {
        *copyArray = NULL;
    }
    else if(nValues > 0 &&
           origArray IS NULL)
    {
        XONTRUE( (*copyArray = NEW_ARRAY(double,nValues)) IS NULL,
                "Failed to alloc copy array" );
        for(idx=0; idx<nValues; idx++)
        {
            (*copyArray)[idx] = defValue;
        }
    }
    else
    {
        XONTRUE( (*copyArray = NEW_ARRAY(double,nValues)) IS NULL,
                "Failed to alloc copy array" );
        for(idx=0; idx<nValues; idx++)
        {
            (*copyArray)[idx] = origArray[idx];
        }
    }

    F_END;
}




/***************************************************************************
 *  MbsCopyLArray
 *  Creates copy of long int array
 *  Sets copyArray to NULL if no values or orig is NULL
 ***************************************************************************/
int EXPORT MbsCopyLArray
   (long     nValues,           /* (I) # elements in origArray */
    long    *origArray,         /* (I) original */
    long   **copyArray)         /* (O) copy */
{
    F_INIT("MbsCopyDArray");
    long   idx;

    if(nValues <= 0 ||
       (origArray IS NULL))
    {
        *copyArray = NULL;
    }
    else
    {
        XONTRUE( (*copyArray = NEW_ARRAY(long,nValues)) IS NULL,
                "Failed to alloc copy array" );
        for(idx=0; idx<nValues; idx++)
        {
            (*copyArray)[idx] = origArray[idx];
        }
    }

    F_END;
}



/***************************************************************************
 *  Fundamental memory alloc funcs
 ***************************************************************************/

/***********************************************************
 *  mbs_calloc
 *  Central "calloc" function
 *  Examples:
 *       XONFAIL( mbs_calloc(2,1000,(char**) &somevec) );
 *    XONFAILMSG( mbs_calloc(2,1000,(char**) &somevec), "failed to alloc somevec" );
 ***********************************************************/
int mbs_calloc(long nvars,	     /* (I) # vars/objs to alloc */
	       long varsize, 	     /* (I) size of each var/obj */
	       char** ptr)	     /* (O) ptr to allocated var/obj */
{
    F_INIT("mbs_calloc");

    /* make sure init has been done */
    verify_mbs_init();

    XONTRUE( (*ptr = (char *) calloc(nvars,varsize)) IS NULL, 
	     "Failed to mbs_calloc() var/struc" );

    F_END;
}

/**************************************************
 *  mbs_malloc
 *  Central "malloc" function
 *  Examples:
 *       XONFAIL( mbs_malloc(1000,(char**) &somevec) );
 *    XONFAILMSG( mbs_malloc(1000,(char**) &somevec), "failed to alloc somevec" );
 **************************************************/
int mbs_malloc(long varsize,	     /* (I) size of var/obj */
	       char** ptr)	     /* (O) ptr to allocated var/obj */
{
    F_INIT("mbs_malloc");

    /* make sure init has been done */
    verify_mbs_init();

    XONTRUE( (*ptr = (char *) malloc(varsize)) IS NULL, 
	     "Failed to mbs_malloc() var/struc" );

    F_END;
}

/**************************************************
 *  mbs_free
 *  Central memory de-allocation function
 **************************************************/
void mbs_free(char* ptr)
{
    /* make sure init has been done */
    verify_mbs_init();

    /* for now, just call plain C func */
    if( ptr )
	free(ptr);

    return;
}


/***************************************************************************
 *  Port of Numerical Recipes memory functions
 *  Altered functions to use our central malloc func,
 *  and to use our return code scheme
 ***************************************************************************/

int mbs_vector(long nl, long nh, float** p)
{
    F_INIT("vector");
    float *pp;

    XONFAILMSG( mbs_malloc((unsigned) (nh-nl+1)*sizeof(float), (char**) &pp), 
		"Failure in mbs_vector()" );
    *p = pp-nl;

    F_END;
}

int mbs_ivector(long nl, long nh, long **p)
{
    F_INIT("ivector");
    long *pp;

    XONFAILMSG( mbs_malloc((unsigned) (nh-nl+1)*sizeof(long), (char**) &pp), 
		"Failed in mbs_ivector()" );
    *p = pp-nl;

    F_END;
}

int mbs_dvector(long nl, long nh, double **p)
{
    F_INIT("dvector");
    double *pp;

    XONFAILMSG( mbs_malloc((unsigned) (nh-nl+1)*sizeof(double), (char**) &pp), 
		"Failed in mbs_dvector()" );
    *p = pp-nl;

    F_END;
}

int mbs_matrix(long nrl, long nrh, long ncl, long nch, float ***p)
{
    F_INIT("mbs_matrix");
    long i;
    float **m;

    XONFAILMSG( mbs_malloc((unsigned) (nrh-nrl+1)*sizeof(float*), (char**) &m), 
		"Failed in mbs_matrix()/1" );
    m -= nrl;
    for( i = nrl; i <= nrh; i++ ) {
	XONFAILMSG( mbs_malloc((unsigned) (nch-ncl+1)*sizeof(float), (char**) &(m[i])), 
		    "Failed in mbs_matrix()/2" );
        m[i] -= ncl;
    }
    *p = m;

    F_END;
}

int mbs_dmatrix(long nrl, long nrh, long ncl, long nch, double ***p)
{
    F_INIT("mbs_dmatrix");
    long i;
    double **m;

    XONFAILMSG( mbs_malloc((unsigned) (nrh-nrl+1)*sizeof(double*), 
                           (char**) &m), 
               "Failed in mbs_dmatrix()/1" );
    m -= nrl;
    for( i = nrl; i <= nrh; i++ ) {
	XONFAILMSG( mbs_malloc((unsigned) (nch-ncl+1)*sizeof(double), 
                               (char**) &(m[i])), 
                   "Failed in mbs_dmatrix()/2" );
        m[i] -= ncl;
    }
    *p = m;

    F_END;
}

int mbs_imatrix(long nrl, long nrh, long ncl, long nch, long ***p)
{
    F_INIT("mbs_imatrix");
    long i;
    long **m;

    XONFAILMSG( mbs_malloc((unsigned) (nrh-nrl+1)*sizeof(long*), (char**) &m), 
		"Failed in mbs_imatrix()/1" );
    m -= nrl;
    for( i = nrl; i <= nrh; i++ ) {
	XONFAILMSG( mbs_malloc((unsigned) (nch-ncl+1)*sizeof(long), (char**) &(m[i])), 
		    "Failed in mbs_imatrix()/2" );
        m[i] -= ncl;
    }
    *p = m;

    F_END;
}

void mbs_free_vector(float *v, long nl, long nh)
{
    mbs_free((char*) (v+nl));

    return;
}

void mbs_free_ivector(long *v, long nl, long nh)
{
    mbs_free((char*) (v+nl));

    return;
}

void mbs_free_dvector(double *v, long nl, long nh)
{
    mbs_free((char*) (v+nl));

    return;
}

void mbs_free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
{
    long i;

    for ( i = nrh; i >= nrl; i-- ) {
        mbs_free((char*) (m[i]+ncl));
    }
    mbs_free((char*) (m+nrl));

    return;
}

void mbs_free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
    long i;

    for ( i = nrh; i >= nrl; i-- ) {
        mbs_free((char*) (m[i]+ncl));
    }
    mbs_free((char*) (m+nrl));

    return;
}

void mbs_free_imatrix(long **m, long nrl, long nrh, long ncl, long nch)
{
    long i;

    for ( i = nrh; i >= nrl; i-- ) {
        mbs_free((char*) (m[i]+ncl));
    }
    mbs_free((char*) (m+nrl));

    return;
}

/***************************************************************************
 *  Misc. string functions
 ***************************************************************************/

/******************************************************************************
 * remove_blanks
 * Removes all blanks (or tabs) from string
 ******************************************************************************/
int remove_blanks(char *s)
{
    F_INIT("remove_blanks");
    long i, j, slen;
    char *stmp = NULL;

    /* make sure init has been done */
    verify_mbs_init();

    /* only do something if non-null & not empty */
    if( (s ISNT NULL) && (s[0] ISNT '\0') ) {
	slen = strlen(s);
	/* create temp. string of sufficient size */
	XONFAIL( mbs_calloc(slen + 1, sizeof(char), (char**) &stmp) );
	j = 0;
	for (i = 0; i < (long) strlen(s); i++)
	    /* only copy char if not: blank, tab */
	    if (s[i] ISNT ' ' && s[i] ISNT '\t')
		stmp[j++] = s[i];
	stmp[j] = '\0';
	strcpy(s, stmp);
    }

    status = SUCCESS;
  done: 
    if( stmp ) {
	mbs_free(stmp);
    }
    return status;
}

/******************************************************************************
 * remove_lead_blanks
 * Removes leading blanks (or tabs) from string
 ******************************************************************************/
int remove_lead_blanks(char *s)
{
    char *p;

    /* make sure init has been done */
    verify_mbs_init();

    /* if ptr non-null */
    if ((p = s)) {			    
	/* skip past any leading blanks or tabs */
	while ((p[0] ISNT '\0') && (p[0] IS ' ' || p[0] IS '\t'))	
	    p++;
	/* shift s backward to eliminate leading blanks */
	strcpy(s, p);			    
    }
    return(SUCCESS);
}

/******************************************************************************
 * remove_trail_blanks
 * Removes trailing blanks (or tabs) from string
 ******************************************************************************/
int remove_trail_blanks(char *s)
{
    long slen;

    /* make sure init has been done */
    verify_mbs_init();

    slen = strlen(s);
    while (slen > 0 && (s[slen - 1] IS ' ' || s[slen - 1] IS '\t'))
	s[--slen] = '\0';
    return(SUCCESS);
}

/******************************************************************************
 * remove_end_blanks
 * Removes leading and trailing blanks (or tabs) from string
 ******************************************************************************/
int remove_end_blanks(char *s)
{
    F_INIT("remove_end_blanks");

    /* make sure init has been done */
    verify_mbs_init();

    XONFAIL( remove_lead_blanks(s) );
    XONFAIL( remove_trail_blanks(s) );
    F_END;
}

/************************************************************************************
 * remove_trail_cr
 * Removes trailing carriage return from string;
 ************************************************************************************/
int remove_trail_cr(char *s)
{
    long slen;

    /* make sure init has been done */
    verify_mbs_init();

    /* if ptr <> NULL */
    if (s) {
	slen = strlen(s);
	if (slen > 1 && s[slen - 1] IS '\n')
	    s[slen - 1] = '\0';
    }
    return(SUCCESS);
}

/******************************************************************************
 * dncase
 * Simple function to convert string to lowercase
 * (created this local routine in case we need special handling)
 ******************************************************************************/
int dncase(char *s)
{
    long i, slen;

    /* make sure init has been done */
    verify_mbs_init();

    slen = strlen(s);
    for (i = 0; i < slen; i++)
	s[i] = tolower(s[i]);
    return(SUCCESS);
}

/******************************************************************************
 * upcase
 * Simple function to convert string to uppercase
 * (created this local routine in case we need special handling)
 ******************************************************************************/
int upcase(char *s)
{
    long i, slen;

    /* make sure init has been done */
    verify_mbs_init();

    slen = strlen(s);
    for (i = 0; i < slen; i++)
	s[i] = toupper(s[i]);
    return(SUCCESS);
}


/******************************************************************************
 * str_npos
 * Finds position of nth occurence of char c in string
 * (or -1 if not found)
******************************************************************************/
int str_npos(char *s,		/* (I) string to search */
	     char c,		/* (I) char to find */
             long n,             /* (I) occurence to search for */
	     long *ipos)		/* (O) position of char */
{
    F_INIT("str_npos");
    long i, nocc;

    /* make sure init has been done */
    verify_mbs_init();

    /* reset return arg */
    *ipos = -1;    

    /* check inputs */
    XONTRUE( n < 1, "Cannot find 0th occurrence" );
    XONTRUE( s IS NULL, "No ptr to string supplied" );

    /* only bother checking if s not empty 
       (but empty string should not cause failure of func) */
    if( s[0] ISNT '\0' ) {
        i = 0;
        nocc = 0;
        while( s[i] ISNT '\0' && *ipos < 0 ) {
            if( s[i] IS c ) {
                nocc++;
                if( nocc IS n )
                    *ipos = i;
            }                
            i++;
        }
    }

    F_END;
}


/******************************************************************************
 * str_pos
 * Finds position of 1st occurence of char c in string
 * (or -1 if not found)
******************************************************************************/
int str_pos(char *s,		/* (I) string to search */
	    char c,		/* (I) char to find */
	    long *ipos)		/* (O) position of char */
{
    F_INIT("str_pos");

    /* just call more general func */
    XONFAIL( str_npos(s,c,1,ipos) );

    F_END;
}


/******************************************************************************
 * str_match
 * Determines if strings match (case-sensitive)
 * NB: differs from strcmp() in that it can handle strings of
 * different length--compares shorter string to the same-length
 * starting piece of longer string
******************************************************************************/
int str_match(char *s1,		/* (I) string 1 */
	      char *s2,		/* (I) string 2 */
	      long *match)	/* (O) boolean true if matching */
{
    long ilen, i, OK;

    /* make sure init has been done */
    verify_mbs_init();

    ilen = MIN(strlen(s1),strlen(s2));
    /* return false (no match) if smaller string is empty */
    if (ilen IS 0)
	OK = FALSE;
    else {
	i = 0;
	do {
	    OK = (s1[i] IS s2[i]);
	    i++;
	} while (OK && i < ilen);
    }
    *match = OK;
    return(SUCCESS);
}


/******************************************************************************
 * parse_bond_price
 * Converts string representation of bond price to floating point #;
 * works with either "nn-mm+" (tick format) or "nn.mm" (decimal) forms
 * Returns boolean (true if successful)
 @@@# leftoff on cleaning up this function, ...
******************************************************************************/
int parse_bond_price(char *sz,	/* (I) string to parse */
		     double *x)	/* (O) converted decimal value */
{
    F_INIT("parse_bond_price");
    long len, ipos, ifrac;
    char *stmp = NULL;
    char *ticks;

    /* make sure init has been done */
    verify_mbs_init();

    /* reset output arg */
    *x = 0.;

    /* check input string */
    XONFAIL( remove_blanks(stmp) );
    XONTRUE( sz IS NULL || sz[0] IS '\0', "Cannot parse null or empty string" );
    len = strlen(sz);
    XONTRUE( len < 3, "price string not long enough to parse" );

    /* create local copy of string */
    XONFAIL( mbs_malloc(len+1,(char**) &stmp) );
    strcpy(stmp,sz);

    /* find hyphen */
    XONFAIL( str_pos(stmp,'-',&ipos) );
    /* if no hyphen found, assume decimal form */
    if( ipos < 0 ) 
	*x = atof(stmp);
    else {
	XONTRUE( ipos IS 0 || ipos IS len - 1,
		 "Hypen in price string cannot be at ends" );
	/* cheat: create two strings, for a moment */
	stmp[ipos] = '\0';			    
	*x = (double) atoi(stmp);
	stmp[ipos] = '-';	    /* restore original string */
	ticks = &(stmp[ipos + 1]);  /* ptr to 2nd part of string */

	/* special case: trailing "+" */
	if(ticks[strlen(ticks)-1] IS '+') {
	    *x += 0.5/32.;
	    ticks[strlen(ticks)-1] = '\0';
	}

	/* must have fractional part in proper format */
	XONTRUE( ticks[0] IS '\0' || strlen(ticks) > 2,
		 "Bad format for fractional part of price string" );
	/* step past leading zero, if any */
	if (strlen(ticks) > 1 && ticks[0] IS '0')	
	    ticks++;
	ifrac = atoi(ticks);
	XONTRUE( ifrac >= 32, "Cannot have >= 32 ticks in price" );
	*x += (double) ifrac / 32.;
    }

    status = SUCCESS;
  done:
    mbs_free(stmp);
    return status;
}


/***************************************************************************
 *  "Perl" data-file parsing functions
 ***************************************************************************/

/***************************************************************************
  Rules for lines in "perl-formatted" ascii flat files:
     -Lines must be: empty, a comment (starts with '#'), or data line
     -Data lines have the format:
         label|N|data1|data2|...|dataN
      as, e.g.,
         Run label | 1 | Fixed Rate Bonds
         Libor|6|6.025|6.025|6.025|6.025|6.025|6.025
     -Fields are delimited by '|' char, though this delimiter
      should NOT be used at start or end of line
     -The first field is a unique string label (case insensitive)
     -Second "count" field N MUST match the number
      of data elements to follow (note, therefore, that
      the total number of fields is N+2)
     -For now, each data line should contain only ONE type
      of data: long integer, double float, or string;
      do not mix these types in a single data line
     -Limit on character length of data lines: 16k bytes
 ****************************************************************************/

/*********************************************************************************
 *  read_perl_long
 *  Gets long int data from matching line in "perl" data file
 *********************************************************************************/
int read_perl_long(FILE *f,           /* (I) file (already open) to read from */
		   char *perl_lbl,    /* (I) label of data line to get */
		   long nvalues,       /* (I) # of values to get */
		   long *outarr)      /* (O) values copied to this array */
{
    F_INIT("read_perl_long");
    long
        i,
        ndata;
    char
        *pdata;

    /* check inputs */
    XONTRUE( f IS NULL, "No file ptr supplied" );    
    XONTRUE( perl_lbl IS NULL || perl_lbl[0] IS '\0',
        "No label supplied" );
    XONTRUE( nvalues < 1, "Bad nvalues" );
    XONTRUE( outarr IS NULL, "No output ptr supplied" );
    
    /* find line in data file w/matching label */
    XONFAIL( find_perl_line(f,perl_lbl,&ndata,&pdata) );
    XONTRUE( ndata < nvalues, "Not enough data in perl line" );

    /* convert each data fragment to long int */
    for(i=0; i<MIN(ndata,nvalues); i++) 
    {
        outarr[i] = atol(pdata);
        pdata += 1+strlen(pdata);
    }

    F_END;
}

/*********************************************************************************
 *  read_perl_dbl
 *  Gets double data from matching line in "perl" data file
 *********************************************************************************/
int read_perl_dbl(FILE *f,           /* (I) file (already open) to read from */
	          char *perl_lbl,    /* (I) label of data line to get */
		  long nvalues,       /* (I) # of values to get */
		  double *outarr)    /* (O) values copied to this array */
{
    F_INIT("read_perl_dbl");
    long
        i,
        ndata;
    char
        *pdata;

    /* check inputs */
    XONTRUE( f IS NULL, "No file ptr supplied" );    
    XONTRUE( perl_lbl IS NULL || perl_lbl[0] IS '\0',
        "No label supplied" );
    XONTRUE( nvalues < 1, "Bad nvalues" );
    XONTRUE( outarr IS NULL, "No output ptr supplied" );
    
    /* find line in data file w/matching label */
    XONFAIL( find_perl_line(f,perl_lbl,&ndata,&pdata) );
    XONTRUE( ndata < nvalues, "Not enough data in perl line" );

    /* convert each data fragment to double */
    for(i=0; i<MIN(ndata,nvalues); i++) 
    {
        outarr[i] = atof(pdata);
        pdata += 1+strlen(pdata);
    }

    F_END;
}

/*********************************************************************************
 *  read_perl_str
 *  Gets strings from matching line in "perl" data file
 *  NB: strings are copied consecutively to outstr, 
 *  separated by null chars
 *********************************************************************************/
int read_perl_str(FILE *f,           /* (I) file (already open) to read from */
		  char *perl_lbl,    /* (I) label of data to get */
		  long nvalues,       /* (I) # of strings to get */
                  long maxlen,        /* (I) max # chars to copy to outstr */
                  long trim_ends,     /* (I) boolean: if true, removes
                                        lead/trail blanks from strings */
		  char *outstr)      /* (O) values copied to this string */
{
    F_INIT("read_perl_str");
    long
        i,
        space,
        ndata;
    char
        *pdata,
        *pnext,
        *pout;

    /* check inputs */
    XONTRUE( f IS NULL, "No file ptr supplied" );    
    XONTRUE( perl_lbl IS NULL || perl_lbl[0] IS '\0',
        "No label supplied" );
    XONTRUE( nvalues < 1, "Bad nvalues" );
    XONTRUE( maxlen < 2, "Bad maxlen" );
    XONTRUE( outstr IS NULL, "No output string ptr supplied" );
    
    /* find line in data file w/matching label */
    XONFAIL( find_perl_line(f,perl_lbl,&ndata,&pdata) );
    XONTRUE( ndata < nvalues, "Not enough data in perl line" );

    /* avail space in outstr */
    space = maxlen-1;
    /* ptr to current location in outstr */
    pout = outstr;
    for(i=0; i<MIN(ndata,nvalues); i++) {
        /* before altering fragment, point to start
           of next frag */
        pnext = pdata+strlen(pdata)+1;
        /* shorten string, if not enough space */
        if( (long) strlen(pdata)+1 > space )
            pdata[space-1] = '\0';
        /* now we can safely copy the output */
        strcpy(pout,pdata);
        /* trim blanks, if requested */
        if( trim_ends ) {
            XONFAIL( remove_end_blanks(pout) );
        }
        /* update space */
        space -= strlen(pout)+1;
        /* advance ptrs */
	pdata = pnext;
        pout += strlen(pout)+1;
    }

    F_END;
}

/* @# later, add perl "write" functions */


/*****************************************************************************
 *  LOCAL FUNCTIONS
 ****************************************************************************/

/****************************************************************************
 *  verify_mbs_init
 *  Local function called by every function in module to
 *  guarantee that initialization has been done
 ****************************************************************************/
static
void verify_mbs_init(void)
{
    /* if this function called, it means application 
       forgot to explicitly init base lib;
       given this, be conservative and init the
       lib in "quiet" mode (no output of any kind) */
    mbs_init("mbsbase",FALSE,NULL);
}


/*****************************************************************
 *  find_behavior_flag
 *  internal func to find behavior label in linked list
 *****************************************************************/
static
int find_behavior_flag(char *behavior_label,     /* (I) label to search fo */
		       BEHAVIOR_STRUC **pfound)	 /* (O) ptr to matching item */
{
    F_INIT("find_behavior_flag");
    char lbl[LINESTRLEN];
    BEHAVIOR_STRUC *p;

    *pfound = NULL;

    /* only bother if there's something in linked list */
    if( p_behavior_start ) {
	/* clean up (copy of) label */
	XONFAIL( clean_up_behavior_label(behavior_label,lbl) );

	strncpy(lbl,behavior_label,LINESTRLEN-1);
	XONFAIL( remove_end_blanks(lbl) );
	XONFAIL( upcase(lbl) );

	p = p_behavior_start;
	while (p && (*pfound IS NULL)) {
	    if( !strcmp(lbl, p->behavior_label) )
		*pfound = p;
	    p = (BEHAVIOR_STRUC*) p->next;
	}
    }

    F_END;
}

/*****************************************************************
 *  clean_up_behavior_label
 *  internal func to put label string in std form
 *****************************************************************/
static
int clean_up_behavior_label(char *behavior_label, /* (I) initial label */
			    char *clean_label)    /* (O) cleaned-up label */
{
    F_INIT("clean_up_behavior_label");

    strncpy(clean_label,behavior_label,LINESTRLEN-1);
    XONFAIL( remove_end_blanks(clean_label) );
    XONFAIL( upcase(clean_label) );

    F_END;
}


/******************************************************************
 *  find_perl_line
 *  Internal func to locate line in perl-formatted flat file
 *  If successful, returns # data frags, ptr to 1st data frag
 ******************************************************************/
static
int find_perl_line(FILE *f,            /* (I) ptr to (open) file */
		   char *perl_lbl,     /* (I) identifier label */
		   long *ndata,         /* (O) # data fragments in line */
		   char **p_data)      /* (O) pointer to start of 1st data frag */
{
  F_INIT("find_perl_line");
  long
    ipass,
    ichar,
    slen,
    delim_pos,
    match,
    keep_looking,
    nvalues_found,
    nvalues_claimed;
  long
    orig_filepos;
  char
    sought_lbl[LINESTRLEN],
    *curr_lbl,
    *pfrag;
  /* static long string for storing current line */
  static char perl_line[MAX_PERL_LEN];

  /* reset return params */
  *ndata = 0;
  *p_data = NULL;

    /* Don't bother to check inputs,
       since they were already checked in public funcs) */

    /* create cleaned-up copy of label */
    strncpy(sought_lbl,perl_lbl,LINESTRLEN-1);
    XONFAIL( remove_end_blanks(sought_lbl) );
    XONFAIL( upcase(sought_lbl) );

  /* we do two passes through file: first starts from current point,
     then we rewind and start from the top (but only goes up to original
     starting point) */
  ipass = 1;
  match = FALSE;
  keep_looking = TRUE;
  /* remember point where we started */
  orig_filepos = ftell(f);    

  /* For each pass through file: */
  while( (!match) && ipass <= 2 ) {

    /* look for matching line */
    while( keep_looking && (!match) 
    && (fgets(perl_line,MAX_PERL_LEN-1,f) ISNT NULL) ) {
        /* strip leading blanks & trailing CR */
        XONFAIL( remove_lead_blanks(perl_line) );
        XONFAIL( remove_trail_cr(perl_line) );

        /* only consider lines that aren't comments or empty lines */
        if( !( perl_line[0] IS '\0' || perl_line[0] IS '\n' 
        || perl_line[0] IS '#') ) {
            /* ASAP: determine if label matches */
            /* isolate label by replacing first delim w/null */
            XONFAIL( str_npos(perl_line,'|',1,&delim_pos) );
            XONTRUE( delim_pos < 0, "No delimiters in perl line" );
            perl_line[delim_pos] = '\0';
            /* convenience ptr to label & start of 2nd field */
            curr_lbl = perl_line;  
            pfrag = perl_line+delim_pos+1;
            /* clean up label */
            XONFAIL( remove_end_blanks(curr_lbl) );
            XONFAIL( upcase(curr_lbl) );
            /* match? */
            if( !strcmp(curr_lbl,sought_lbl) ) {
                match = TRUE;
                /* convert all subseq. delims to nulls */
                slen = strlen(pfrag);
                nvalues_found = 0;
                for(ichar=0; ichar<slen; ichar++) {
                    if( pfrag[ichar] IS '|' ) {
                        nvalues_found++;
                        pfrag[ichar] = '\0';
                    }
                }
  	        /* Must always have > 1 frags */
                XONTRUE( nvalues_found < 1, "Too few fields in perl line" );
                /* determine # values claimed to be in this line */
                nvalues_claimed = atol(pfrag); 
                XONTRUE( nvalues_claimed < 1, "Perl line must have some data" );
                XONTRUE( nvalues_claimed ISNT nvalues_found,
                    "Mismatch in # data values in perl line" );
                /* set outputs */
                *ndata = nvalues_found;
                /* then point to start of 3rd field (1st data field) */
                *p_data = pfrag+strlen(pfrag)+1;
            }  /* if match */

        }  /* if-not-comment */

        /* On first pass, keep reading until EOF;
           on second pass, only go up to point at which 
           we started on first pass */
        if(ipass IS 1)
	    keep_looking = TRUE;
        else
   	    keep_looking = (ftell(f) <= orig_filepos);

        }  /* while-keep-looking */

    /* if we didn't get it on first pass, start again at the top */
    if(ipass IS 1 && !match)
      rewind(f);
    ipass++;
    }  /* end of curr pass thru file */            
  
  status = SUCCESS;
 done:
  return status;
}

