#include "driwbsp.h"

#include <math.h>
#include <string.h>
        extern int      optind;
        extern char     *optarg;



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


	if (DriSpreadBarrierBinaryW("wbinswst_t.dat") != SUCCESS)
		goto done;


	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}


