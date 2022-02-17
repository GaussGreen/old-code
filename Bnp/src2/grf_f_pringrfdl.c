#include "grf_h_all.h"



Err grf_f_pringrfrng(
		FILE *out,
		long numgrng,
		GrfnRng *grng 
			)
{
   fprintf(out,"%s\n",GRF2020LINE);
   fprintf(out,"%s\n",GRFBEGGRN);
   fprintf(out,"%s\n",GRF2020LINE);


   fprintf(out,"%s\n",GRF2020LINE);
   fprintf(out,"%s\n",GRFENDGRN);
   fprintf(out,"%s\n",GRF2020LINE);
   return NULL;
}


Err grf_f_prinauxrng(
		FILE *out,
		long auxwidth,
		long *auxlen,
		double **aux
			)
{
   long i,j,maxlen;

   fprintf(out,"%s\n",GRF2020LINE);
   fprintf(out,"%s\n",GRFBEGAUX);
   fprintf(out,"%s\n",GRF2020LINE);

   if(auxwidth>0){
	fprintf(out,"%d",auxwidth);
	fprintf(out,"%c\n",GRFSEPCHAR);
	for(i=0;i<auxwidth;i++){
	  fprintf(out,"%d",auxlen[i]);
	  fprintf(out,"%c",GRFSEPCHAR);
	}
	fprintf(out,"\n");
	fprintf(out,"%s\n",GRF2020LINE);
	maxlen=0;
	for(i=0;i<auxwidth;i++)maxlen=IMAX(maxlen,auxlen[i]);
	for(j=0;j<maxlen;j++){
	  for(i=0;i<auxwidth;i++){
		if(j<auxlen[i]){
		  fprintf(out,"%.10lf",aux[i][j]);
		}else{
		  fprintf(out,"%c",GRFBLNCHAR);
		}
		fprintf(out,"%c",GRFSEPCHAR);
	  }
	  fprintf(out,"\n");
	}
   }
   fprintf(out,"%s\n",GRF2020LINE);
   fprintf(out,"%s\n",GRFENDAUX);
   fprintf(out,"%s\n",GRF2020LINE);
   return NULL;
}

Err grf_f_princells(
		FILE *out, 
		Date *eventdates,
		GrfnCell **gcells,
		long nr, 
		long nc
		)
{
   long i,j;

   fprintf(out,"%s\n",GRF2020LINE);
   fprintf(out,"%s\n",GRFBEGEV);
   fprintf(out,"%s\n",GRF2020LINE);

   if(nr>0 && gcells != NULL){
	fprintf(out,"%d%c%d%c\n",
		nr,GRFSEPCHAR,nc,GRFSEPCHAR);
	fprintf(out,"%s\n",GRF2020LINE);
	for(i=0;i<nr;i++){
	  fprintf(out,"%d",eventdates[i]);
	  fprintf(out,"%c",GRFSEPCHAR);
	  for(j=0;j<nc;j++){
		switch(gcells[i][j].type){
		  case GRFNSCELL:
			fprintf(out,"%s",gcells[i][j].sval);
			break;
		  case GRFNDCELL:
			fprintf(out,"%.10lf",gcells[i][j].dval);
			break;
		  default:
			fprintf(out,"%c",GRFBLNCHAR);
			break;
		}
		fprintf(out,"%c",GRFSEPCHAR);
	  }
	  fprintf(out,"\n");
	}
   }

   fprintf(out,"%s\n",GRF2020LINE);
   fprintf(out,"%s\n",GRFENDEV);
   fprintf(out,"%s\n",GRF2020LINE);
   return NULL;
}



/** PRINT A GRFN DEAL TO A FLAT FILE FROM WHICH IT CAN BE READ BY MAD **/

Err grf_f_pringrfdl(FILE *out, GrfnDeal *gd)
{
   	Err 	err;


	err = grf_f_pringrfrng(
		out,
		gd->num_grng,
		gd->grng
		);
	if(err)return err;

	err = grf_f_prinauxrng(
		out,
		gd->auxwidth,
		gd->auxlen,
		gd->aux
		);
	if(err)return err;

	err = grf_f_princells(
		out, 
		gd->event_dates,
	 	gd->gcells,
		gd->sslength, 
		gd->sswidth
		);
	if(err)return err;

   	return NULL;
}

Err grf_f_prgrdl(
		String		filename,
		int			numeventdates,
		Date		*eventdates,
		long		nrows,
		long		ncols,
    	GrfnCell	**sprdsht,
		long		numgrng,
		GrfnRng		*grng,
		long   		auxwidth,
		long		*auxlen,
		double		**aux)
{
   	GrfnDeal 	gd;
   	Err 		err=NULL;
   	FILE 		*f;
   	int 		i;	

/*** ERROR CHECK ***/
   	if(numeventdates < 1 || numeventdates != nrows || ncols < 1)
		return GRERR_BAD_EVDIM;
   	for(i=0;i<numeventdates-1;i++)
   	{
     		if(eventdates[i]>=eventdates[i+1])
     		{
       			return serror("GRF:need incr. event dates.");
     		}
   	}	

/*** PUT DATA IN GRFN DEAL STRUCTURE **********************************/
   	grfn_copy_to_GrfnDeal(&gd,nrows,ncols,sprdsht,eventdates,
		numgrng,grng,auxwidth,auxlen,aux);

	f = fopen(filename,"w");

	grf_f_pringrfdl(f,&gd);

	fclose(f);

	return NULL;
}




