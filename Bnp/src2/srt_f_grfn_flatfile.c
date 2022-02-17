/*****************************************************************************
   FUNCTION	: srt_f_grfn_flatfile.c
    	
   AUTHOR	: E.AULD	 20jul94

   DESCRIPTION	:reads a srt model and inputs to grfn from 
	a flat file, calls grfn.
 	
  		
    AMENDMENTS:
 	Reference	:
 	Author          :J.Malhi
 	Date            :27 Oct 1994
 	Description     :amend to read ts from ssheet

******************************************************************************/

#include "srt_h_all.h"
#include "grf_h_all.h"

Err srt_f_grfn_from_flat_file(
		String 		filename,
		SrtGrfnParam *grfnparam,
		SrtUndPtr 	und,
		double 		*answer)
{
FILE 		*in;

Date 		*eventdates;
long 		nrow;
long 		ncol;
GrfnCell 	**sprdsht;
long 		numgrng		= 0;
GrfnRng 	*grng;
long 		auxwidth	= 0;
long 		*auxlen;
double 		**aux;
Err err;
	
	
    grng	= NULL;
    aux		= NULL;
	auxlen	= NULL;

	in  = fopen(filename,"r");
    if(!in)
		return serror("couldn't open %s\n",filename);


	err=grf_f_readgrfrng(	in,
							&numgrng,
							&grng);
	if(err)
		return err;
	
	err = grf_f_readauxrng( in,
							&auxwidth,
							&auxlen,
							&aux);
	if(err)
		return err;
	
	err = grf_f_readcells(	in,
							&eventdates,
							&sprdsht,
							&nrow, 
							&ncol);
	if(err)
		return err;


	fclose(in);

   	err = srt_f_grfn(	und,
				grfnparam,
				nrow,&eventdates,
				&nrow,&ncol,&sprdsht,
        		numgrng,grng,
				auxwidth,auxlen,aux,
				answer,
				0, 0);  
   
   	if(aux)
	{
   		if(aux[0])
			srt_free(aux[0]);
		srt_free(aux);
   	}
   	
	if(auxlen)
		srt_free(auxlen);
   	
	if(sprdsht)
	{
		grfn_free_GrfnCellmatrix(sprdsht,nrow,ncol);
	}

	if(eventdates)
		srt_free(eventdates);
   	
	if(grng)
		srt_free(grng);/*FIX*/

   	return err;
}
