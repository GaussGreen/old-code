#include "cerror.h"
#include "pfmopt.h"



/*----------------------------------------------------------------------
 *
 */

int
main(int argc, char **argv)
{
static	char	routine[] = "main";
	int	status = FAILURE;

	GtoErrMsgOn();
	GtoErrMsgFilePointer(stdout);

	if (DriOutPerformanceIdxOptionW("pfmopt_w.dat") != SUCCESS)
		goto done;


	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}


