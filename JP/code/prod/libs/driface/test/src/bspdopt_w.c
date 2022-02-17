#include "cerror.h"
#include "bspdopt.h"

/*----------------------------------------------------------------------
 *
 */

int
main(int argc, char **argv)
{
static	char	routine[] = "main";
	int	status = FAILURE;

	GtoErrMsgFilePointer(stdout);
	GtoErrMsgOn();


	if (DriBasisSpreadOptW("bspdopt_w.dat") != SUCCESS)
		goto done;


	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}


