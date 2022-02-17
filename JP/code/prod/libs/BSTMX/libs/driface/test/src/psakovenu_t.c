#include "drlgetop.h"
#include "dripsavenu.h"

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
	int	c;


	/*GtoErrMsgFileName("psakovenu_t.log", TRUE);*/
	GtoErrMsgFilePointer(stdout);
	GtoErrMsgOn();


  	while ((c = getopt(argc, argv, "hv")) != EOF)
  	switch (c) {
	case 'v':
		GtoLoggingSet(1);
		break;
	case 'h':
		fprintf(stderr, "usage: %s [-hv]\n", argv[0]);
		exit(0);
		break;
	default:
		break;
	}
	argv += optind - 1 ;




	if (DriPsaSwapKoVenuW("psakovenu_t.dat") != SUCCESS)
		goto done;


	/* OK */
	status = SUCCESS;
done:
	if (status != SUCCESS)
		GtoErrMsg("%s: failed.\n", routine);
	return(status);
}


