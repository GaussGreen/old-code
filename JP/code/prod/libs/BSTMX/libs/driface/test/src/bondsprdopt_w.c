#include "cerror.h"
#include "bondsprdopt.h"



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

	if (DriOutSpreadStrikeptionW("bondsprdopt_w.dat") != SUCCESS)
		goto done;


	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}


