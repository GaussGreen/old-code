#include "grf_h_all.h"

Err grf_f_princomll(COMLL_PTR com, FILE *out)
{
	int i;

	if(!com)return NULL;
	com=comll_gototop(com);	
	while(com){
		fprintf(out,"COMLL\n");
		fprintf(out,"%.4lf\n",com->dval);
		for(i=0;i<com->dfindlen;i++)
		{
			fprintf(out,"%d %.4lf\n",com->dfind[i],com->cvg[i]);
		}
		com=com->next;
	}
	return NULL;
}



Err grf_f_prinevent(GrfnEvent *ev, FILE *out)
{
	int ix=0;
	Err err;
	if(!ev)
	{
	  fprintf(out,"NULLEVENT\n");
	  return NULL;
	}
	fprintf(out,"EVENT time %.4lf\n",ev->t);	
	err=grf_f_princomll(ev->icom, out);
	return err;
}


