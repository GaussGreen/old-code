#include "cerror.h"
#include "dritridx.h"



/*----------------------------------------------------------------------
 *
 */

static const char* versionInfo = "VERSION 4.5 COMPILED " __DATE__ " " __TIME__;

int
main(int argc, char **argv)
{
static	char	routine[] = "main";
	int	status = FAILURE;

	GtoErrMsgFilePointer(stdout);
	GtoErrMsgOn();


	if (DriIdxForwardSwapCapFloorW("dritridx_w.dat") != SUCCESS)
		goto done;


	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}


