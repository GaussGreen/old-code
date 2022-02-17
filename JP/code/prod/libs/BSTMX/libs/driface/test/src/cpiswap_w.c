#include "cerror.h"
#include "cpiswap.h"


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


	if (DriCpiSwapW("cpiswap_w.dat") != SUCCESS)
		goto done;


	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}


