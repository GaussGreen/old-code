#include "drlptable.h"
#include "drlgetop.h"		/* getopt() */
#include "drlstr.h"
#include "drlio.h"
#include "drltime.h"

#include "drical.h"

#include <errno.h>
#include <math.h>
#include <string.h>

#define	__PROGRAM__	"caldrw_t"
#define	__DEBUG__
#undef	__DEBUG__

static	void	usage(char *s, ...);

#define READ_DATA(fp, type, ptr, str) \
        { if (DrlFScanVType(fp, type, (void*) ptr) != SUCCESS) \
        { GtoErrMsg("%s: can't read %s.\n", routine, str); goto done;}}

/*----------------------------------------------------------------------
 *
 */

int
main(int argc, char **argv)
{
static	char	routine[256];
	int	status = FAILURE;
	int	c;

	FILE		*fpIn =NULL;
	FILE		*fpCur = NULL;

	TDrWrapperData	*drwData = NULL;
	char	dataFnam[256];
	char	pathdir[256];
	char	tagName[256];

	int		doLiverate = FALSE;
	int		doDate = FALSE;
	int		doPrintEnv = FALSE;
	char		curFnam[256];		/* Currency config file */
	char		curCode[] = "xxx";
	char		logFnam[256];
	char		logLink[256];
        int		mmDen = 360;
        int		bvFreq = 4;
        int		swapFreq = 2;
        TDayCount	swapDcc = GTO_B30_360;

	char		yldcrvBnam[256];	/* yield curve file basename */
	char		baseAtmBnam[256];	/* ATM base vol file basename */
	char		swoAtmBnam[256];	/* ATM swaption vol file */

	strncpy(routine, argv[0], sizeof(routine));
	GtoErrMsgFilePointer(stdout);
	GtoErrMsgOn();


	/*
	 * Parse options
	 */
	sprintf(dataFnam, "%s.dat", __PROGRAM__);
	sprintf(logFnam, "");
	sprintf(logLink, "");
	strcpy(pathdir, "");

#ifdef	__DEBUG__
		/*!!! DEBUGGING */
		strcpy(pathdir, "tmp.liverate");
		strcpy(curFnam, "config/config.usd");
		doLiverate = TRUE;
#endif

	while ((c = getopt(argc, argv, "c:d:ehi:l:tvL:R:")) != EOF)
	switch (c) {
	case 'c':
		ASSERT_OR_DONE(sscanf(optarg, "%s", curCode) == 1);
		IF_FAILED_DONE( DrlStrSubsEnv(curCode));
		break;
	case 'd':
		ASSERT_OR_DONE(sscanf(optarg, "%s", pathdir) == 1);
		IF_FAILED_DONE( DrlStrSubsEnv(pathdir));
		break;
	case 'i':
		ASSERT_OR_DONE(sscanf(optarg, "%s", dataFnam) == 1);
		IF_FAILED_DONE( DrlStrSubsEnv(dataFnam));
		break;
	case 'l':
		ASSERT_OR_DONE(sscanf(optarg, "%s", curFnam) == 1);
		IF_FAILED_DONE( DrlStrSubsEnv(curFnam));
		doLiverate = TRUE;
		break;
	case 'e':
		doPrintEnv = TRUE;
		break;
	case 't':
		doDate = TRUE;
		break;
	case 'v':
		GtoLoggingSet(TRUE);
		break;
	case 'L':
		ASSERT_OR_DONE(sscanf(optarg, "%s", logLink) == 1);
		IF_FAILED_DONE( DrlStrSubsEnv(logLink));
		break;
	case 'R':
		ASSERT_OR_DONE(sscanf(optarg, "%s", logFnam) == 1);
		IF_FAILED_DONE( DrlStrSubsEnv(logFnam));
		break;
	default:
	case 'h':
		usage(NULL);
		status = SUCCESS;
		goto done;
	}
	argv += optind - 1 ;
	/*if (argv[1] == NULL) usage(NULL);*/


	/*
	 * Read market data
	 */
	if (doLiverate) {
		if ((fpCur = fopen(curFnam, "r")) == NULL) {
			GtoErrMsg("%s: can't open `%s' (%s).\n",
				routine, curFnam, strerror(errno));
			goto done;
		}
		READ_DATA(fpCur, DRL_CHAR_ARRAY_T, curCode, "currency code")
		READ_DATA(fpCur, DRL_INT_T, &mmDen, "mm denominator")
		READ_DATA(fpCur, DRL_INT_T, &bvFreq, "bv Freq")
		READ_DATA(fpCur, DRL_INT_T, &swapFreq, "swap Freq")
		READ_DATA(fpCur, DRL_TDAYCOUNT_T, &swapDcc, "swap dcc")
 
		READ_DATA(fpCur, DRL_CHAR_ARRAY_T, 
				yldcrvBnam,  "yield curve file basename")
		READ_DATA(fpCur, DRL_CHAR_ARRAY_T, 
				baseAtmBnam, "ATM base vol file basename")
		READ_DATA(fpCur, DRL_CHAR_ARRAY_T, 
				swoAtmBnam,  "ATM swaption vol file basename")

		IF_FAILED_DONE(DriTDrWrapperDataGetLiverate(
        		pathdir,
        		curCode,
			mmDen,
			bvFreq,
			swapFreq,
			swapDcc,
			yldcrvBnam,
			baseAtmBnam,
			swoAtmBnam,
        		&drwData));
	
	} else {
		/* Wrapper style */

		IF_FAILED_DONE( DriTDrWrapperDataGet(
			(pathdir && pathdir[0] ? pathdir : NULL),
			&drwData));
	}



	/* If just need dates */
	if (doDate) {
		fprintf(stdout, "%04d-%02d-%02d",
			DrlYEAR(drwData->fToday),
			DrlMONTH(drwData->fToday),
			DrlDAY(drwData->fToday));
		status = SUCCESS;
		goto done;
	}

	if (doPrintEnv) {
		DriTDrWrapperDataFpWrite(
			drwData,
			stdout);
		status = SUCCESS;
		goto done;
	}



	/* Set up tag name */
	sprintf(tagName, "%s",
		(strcmp(dataFnam, __PROGRAM__".dat") ?
			DrlGetRootName(NULL, dataFnam) :
			"MRCALIBRATION" ));


	if ((fpIn = fopen(dataFnam, "r")) == NULL) {
		GtoErrMsg("%s: can't open `%s' (%s).\n",
			routine, dataFnam, strerror(errno));
		goto done;
	}


	fprintf(stdout, "logFnam = `%s'\n", logFnam);
	fprintf(stdout, "curCode = `%s'\n", curCode);
	fprintf(stdout, "tagName = `%s'\n", tagName);


	while (DrlFAdvanceToNextChar(fpIn) == SUCCESS) {
		IF_FAILED_DONE( DriCalibrateVnfmParams(
			logFnam,
			logLink,
			curCode,
			tagName,
			fpIn,
			drwData));
	}




	/* OK */
	status = SUCCESS;
done:
	if (fpIn) fclose(fpIn);
	if (fpCur) fclose(fpCur);
	DriTDrWrapperDataFree(drwData);

	if (status != SUCCESS)
		GtoErrMsg("%s: failed (-h for help).\n", routine);
	return(status);
}


/*------------------------------------------------------
 *
 */

static	void
usage(char *s, ...)
{
	va_list	ap;
	char	buf[255];

	if (s != NULL) {
		va_start(ap, s);
		vsprintf(buf, s, ap);
		va_end(ap);

		fprintf(stderr, "ERROR: %s (-h for help)\n", buf);
		exit(1);
	} else {
		fprintf(stdout, "\
USAGE\n\
    %s [options] [<input dat file>]\n\
DESCRIPTION\n\
    Performs a batch MR calibration and generates a detailed report.\n\
    Supports DRWrapper or liverate format market data input.\n\
OPTIONS\n\
    -c <code>	Use <code> for currency code (in DRW mode)\n\
    -t          Just scans and returns the the environment read\n\
                    (do not perform calibration, for debugging)\n\
    -d <dir>	Use <dir> to get DRW data or liverate\n\
    -i <file>	Use <file> as input file (default is "__PROGRAM__".dat)\n\
    -l <file>	Use liverate format (default is DRW) and <file> to get currency information\n\
    -R <file>	Use <file> for detailed output report (default is stdout)\n\
    -L <link>	Add html <link> to the 1-line report (default is none)\n\
    -t          Just scans and returns the base date of the environment read\n\
                    (do not perform calibration)\n\
    -v          Verbose\n\
INPUT FILE FORMAT\n\
\n", __PROGRAM__);
		fprintf(stdout, "    Created %s  %s\n",
			__DATE__,__TIME__);
		fflush(stdout);
		exit(0);
	}
}




