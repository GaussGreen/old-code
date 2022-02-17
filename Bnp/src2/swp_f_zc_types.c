/* ===========================================================================
   
	 FILENAME:     swp_f_zc_types.c

     PURPOSE:      Some function for the structures used in stripping

   =========================================================================== */

#include "swp_h_all.h"

ZC_Obj 	*new_ZC_Obj()
{
	ZC_Obj *new_zc_obj;
	new_zc_obj = (ZC_Obj *) srt_calloc(1,sizeof(ZC_Obj));
        if(!new_zc_obj) return NULL;	

	new_zc_obj->data = new_Field_List();
	if(!new_zc_obj->data) return NULL;

	Field_List_set_Dwidth(new_zc_obj->data, (int)ZC_DWIDTH);
	Field_List_set_Iwidth(new_zc_obj->data, (int)ZC_IWIDTH);
	return new_zc_obj;
}

YC_Obj 	*new_YC_Obj()
{
	YC_Obj *new_yc_obj;
	new_yc_obj = (YC_Obj *) srt_calloc(1,sizeof(YC_Obj));
        if(!new_yc_obj) return NULL;	

	new_yc_obj->data = new_Field_List();
	if(!new_yc_obj->data) return NULL;

	Field_List_set_Dwidth(new_yc_obj->data, (int)YC_DWIDTH);
	Field_List_set_Iwidth(new_yc_obj->data, (int)YC_IWIDTH);
	return new_yc_obj;
}

Leg_Obj	*new_Leg()
{
	Leg_Obj *new_leg;
	new_leg = (Leg_Obj *) srt_calloc(1,sizeof(Leg_Obj));
        if(!new_leg) return NULL;	

	new_leg->data = new_Field_List();
	if(!new_leg->data) return NULL;

	Field_List_set_Dwidth(new_leg->data, (int)LEG_DWIDTH);
	Field_List_set_Iwidth(new_leg->data, (int)LEG_IWIDTH);
	return new_leg ;
}

Swap_Obj*new_Swap()
{
	Swap_Obj *new_swap;
	new_swap = (Swap_Obj *) srt_calloc(1,sizeof(Swap_Obj));
	return new_swap;
}

Arg_Obj *new_Arg()
{
	Arg_Obj *new_arg;
	new_arg = (Arg_Obj *) srt_calloc(1,sizeof(Arg_Obj));
	return new_arg;	
}

/**************** free functions *********************************/

void free_Swap(Swap_Obj *s)
{
	int i,nl = Swap_Iget(s,SWAP_NUM_LEG);
	for(i=0;i<nl;i++)
	{
		free_Leg(Swap_field_get(s,SWAP_LEG,i));
	}
	srt_free (s);
}

void free_Arg(Arg_Obj *ga)
{
	srt_free(ga);
}

void free_YC_Obj(YC_Obj *yco)
{
	free_Field_List(yco->data);
	srt_free(yco);
}

void free_ZC_Obj(ZC_Obj *zco)
{
	free_Field_List(zco->data);
	srt_free(zco);
}

void free_Leg(Leg_Obj *leg)
{
	free_Field_List(leg->data);
	srt_free(leg);
}

/**
	initialize a Arg according to the 
	sort of information it needs to carry, and
	where it needs to carry it.
**/

/**
	create an Arg and
	populate it with defaults appropriate
	to the purpose for which it will be used
	as shown by m
**/


Err 	init_Arg(Arg_Obj **ga, Message m)
{
	*ga = new_Arg();
	if(!*ga) return "new_Arg() failed!";

	switch(m){
	case GENERATE_SWAP:
	case COMPUTE_FRA:
	case COMPUTE_DISC_FACTOR:
	case COMPUTE_ZERO_RATE:
	default:
		break;
	}
	return 0;
}
	
/****** allocate a leg AND give it size pay_size ***********/
	
void leg_allocate(Leg_Obj **leg, int pay_size)
{       
	Leg_Obj *new_leg;
 	new_leg = new_Leg();
	
	Leg_field_set_Ilength(new_leg, pay_size);
	Leg_field_set_Dlength(new_leg, pay_size);

	*leg = new_leg;
}

/********************* field list ***************************/

Err Field_List_set_Iwidth(Field_List *fl, int w)
{
        /* set number of fields */
	
	fl->Iwidth = w;
	if(!(fl->Ilengths = (int *) srt_calloc(w, sizeof(int))))
		return "srt calloc failed.!";
	if(!(fl->Ifields = (long **) srt_calloc(w, sizeof(long *))))
		return "srt calloc failed.!";
	return 0;
}
Err Field_List_set_Dwidth(Field_List *fl, int w)
{
        /* set number of fields */
	
	fl->Dwidth = w;
	if(!(fl->Dlengths = (int *) srt_calloc(w, sizeof(int))))
		return "srt calloc failed.!";
	if(!(fl->Dfields = (double **) srt_calloc(w, sizeof(double *))))
		return "srt calloc failed.!";
	return 0;
}

Err Field_List_set_Dlength(Field_List *fl, int l)
{
        /* set length of fields */

	int i;	
	for(i=0;i<fl->Dwidth;i++){
		fl->Dlengths[i] = l;
		if(!(fl->Dfields[i] = (double *)srt_calloc(l, sizeof(double))))
			return "srt calloc failed.!";
	}
	return 0;
}
Err Field_List_set_Ilength(Field_List *fl, int l)
{
        /* set length of fields */
	int i;	
	for(i=0;i<fl->Iwidth;i++){
		fl->Ilengths[i] = l;
		if(!(fl->Ifields[i] = (long *)srt_calloc(l, sizeof(long))))
			return "srt calloc failed.!";
	}	
	return 0;
}

Field_List *new_Field_List(void){
	Field_List *new_field_list;
	new_field_list = (Field_List *) srt_calloc(1, sizeof(Field_List));
	return new_field_list;
}

void free_Field_List(Field_List *fl)
{
	int i;

	for(i=0;i<fl->Dwidth;i++)
		srt_free(fl->Dfields[i]);
	srt_free(fl->Dfields);
	srt_free(fl->Dlengths);

	for(i=0;i<fl->Iwidth;i++)
		srt_free(fl->Ifields[i]);
	srt_free(fl->Ifields);
	srt_free(fl->Ilengths);

	srt_free(fl);
}
