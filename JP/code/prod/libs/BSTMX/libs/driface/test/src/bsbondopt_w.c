#include "cerror.h"
#include "bsbondopt.h"



/*----------------------------------------------------------------------
 *
 */

static const char* versionInfo = "VERSION 1.0 COMPILED" __DATE__ " " __TIME__;

int
main(int argc, char **argv)
{
static	char	routine[] = "main";
	int	status = FAILURE;

	GtoErrMsgOn();
	GtoErrMsgFilePointer(stdout);

	if (DriOutBSBondOptionW("bsbondopt_w.dat") != SUCCESS)
		goto done;


	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}


