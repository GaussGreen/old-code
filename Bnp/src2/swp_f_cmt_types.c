	
/*
	swp_f_cmt_types.c
*/
#include "swp_h_all.h"
#include "swp_h_cmt_types.h"

/*******************  allocation functions ****/

FwdT_Obj * new_FwdT_Obj()
{
	FwdT_Obj * obj;
	obj = (FwdT_Obj *) srt_calloc(1,sizeof(FwdT_Obj));
    if(!obj) 
		return NULL;	

	obj->data = new_Field_List();
	if(!obj->data) 
		return NULL;

	Field_List_set_Dwidth(obj->data, (int)FwdT_DWIDTH);
	Field_List_set_Iwidth(obj->data, (int)FwdT_IWIDTH);
	
	return obj;
}

CMT_Obj * new_CMT_Obj()
{
	CMT_Obj * obj;
	obj = (CMT_Obj *) srt_calloc(1,sizeof(CMT_Obj));
    if(!obj) 
		return NULL;	

	obj->data = new_Field_List();
	if(!obj->data) 
		return NULL;

	Field_List_set_Dwidth(obj->data, (int)CMT_DWIDTH);
	Field_List_set_Iwidth(obj->data, (int)CMT_IWIDTH);
	
	return obj;
}

void free_CMT_Obj(CMT_Obj * cmto)
{
	free_Field_List(cmto->data);
	srt_free(cmto);
}

void free_FwdT_Obj(FwdT_Obj * fwdTo)
{
	free_Field_List(fwdTo->data);
	srt_free(fwdTo);
}

 
