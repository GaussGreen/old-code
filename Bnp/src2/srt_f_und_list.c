/* =============================================================================
   FILENAME:   srt_f_und_list.c

   PURPOSE:    functions to store and handle underlyings in a list...
   ============================================================================= */

/* ------------------------------------------------------------------------------- 
                              Include Statements  
   ------------------------------------------------------------------------------- */

#include "srt_h_all.h"
#include "srt_h_und_list.h"

/* ------------------------------------------------------------------------------- 
                            Static Declarations  
   ------------------------------------------------------------------------------- */

/* ------ The Static Pointer used to refer to the Underlying List ---------------- */ 
static SrtUndListPtr _srt_underlyings_list = NULL;


/* ------------------------------------------------------------------------------- 
                   Functions to operate on the underlying list  
   ------------------------------------------------------------------------------- */

/* Allocate memory for the list where all the underlyings will be stored */
Err create_underlying_list(String und_list_name)
{
	Err err ;

	err=srt_f_lstcreate(&_srt_underlyings_list,und_list_name);
	if (err)
		return serror("Could not create underlying list %s",und_list_name) ;

	return NULL ;
}

/* -------------------------------------------------------------------------- */

/* Destroy the list where all the underlyings are stored, as well as the underlyings */
Err destroy_all_underlyings()
{
Err err;

	err = srt_f_lstfree(_srt_underlyings_list,SRT_YES);
	if (err)
		return serror("%s in destroy_all_underlyings",err) ;
	
	srt_free(_srt_underlyings_list);

	return NULL;	
}

/* -------------------------------------------------------------------------- */

/* An easy way to access the static pointer to the underlying list */
SrtUndListPtr get_underlying_list(void)
{
 	return _srt_underlyings_list;
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- 
                    Functions to operate on a single underlying 
   -------------------------------------------------------------------------- */

/* Check if a list has go an underlying */
SRT_Boolean srt_f_isundinlist(SrtUndListPtr und_list, String und_name)
{
  SRT_Boolean isthere = SRT_NO ;

  isthere = srt_f_lsthas(*und_list,und_name,0) ;
  
  return (isthere) ;
}/* END SRT_Boolean srt_f_isundinlist(...) */

/* -------------------------------------------------------------------------- */

/* Get the Underlying (structure) corresponding to the name in a list */
Err srt_f_getundinlist(SrtUndListPtr und_list,String und_name, SrtUndPtr *und ) 
{
int               found = 0 ;
Err               err = NULL;
SrtListObject    *obj;

/* Get the underlying == SrtListObject in the underlying_list */
	err  = srt_f_lstgetobj(*und_list, und_name, 0.0, &obj);
	if (err)
		return serror ("Could not find %s underlying in list",und_name);

/* The underlying is the pval of the SrtObject */
	(*und) = (SrtUndPtr) (obj->val.pval);

/* Return a success message */
	return NULL;

} /* END Err srt_f_getundinlist(...) */
				
/* --------------------------------------------------------------------------------- */
/* Function to free an underlying object in a linked list  */
Err srt_f_unddescfree(SrtUndDesc *spec_desc)
{
Err	  err  = NULL;

	if (!spec_desc)
		return NULL;

/* Free the TermStructure attached first (SrtUndDesc = SrtUndPtr) */
	err = free_underlying_ts(spec_desc); 
	if (err)
		return err;

/* Free the remainder of the spec_desc */
	srt_free(spec_desc->spec_desc);
	srt_free(spec_desc);

	return NULL;

}/* END Err srt_f_unddescfree(...) */

/* --------------------------------------------------------------------------------- */
/* Function to pass to srt_f_lstins to free an underlying object in a linked list  */
Err srt_f_undvalfree(void *unddesc)
{
SrtUndDesc *spec_desc = (SrtUndDesc*)(unddesc);
Err   err = NULL;

	err = srt_f_unddescfree(unddesc);
	
	return err;
}
/* --------------------------------------------------------------------------------- */

/* Add an underlying to a list (overwrite the previous und if exist ) */

Err srt_f_addundtolist (
			SrtUndListPtr  und_list,
			String		   und_name,       /* == und_name in SrtUndDesc */
 			String		   und_lbl,        /* FX_UND, IR_UND, EQ_UND */
			String         und_ccy,
			String		   mdl_lbl,
			String		   crv_name1,      /* == discount curve */
			String		   crv_name2,      /* == dividend curve */
			String		   crv_name3,      /* == repo curve */
			TermStruct	   *ts,
			double		   spot)
{
SrtUndDesc	    *und ;
SrtIrDesc	    *undint ;
SrtEqDesc	    *undeq ;
SrtFxDesc	    *undfx ;
SrtErr	        err  ;

/* Create space for SrtUndDesc */

  und = (SrtUndDesc *) srt_calloc(1,sizeof(SrtUndDesc)) ;

/* Checks the underlying type */ 
	err = srt_f_interp_under(und_lbl, &(und->underl_type));
	if (err) 	
		return err;
	strcpy(und->underl_lbl,und_lbl) ;

/* Sets underlying name in the SrtUndDesc: uppercase it  */
	strcpy(und->underl_name,und_name) ;
	strupper(und->underl_name);
	strip_white_space(und->underl_name);

/* Store the currency */
	und->underl_ccy = und_ccy;

/* Set up the elements in the underlying object val */
	switch(und->underl_type)
	{
		case INTEREST_RATE_UND :
			undint = (SrtIrDesc *) srt_calloc(1,sizeof(SrtIrDesc)) ;

			if (mdl_lbl!=0)
			{
				if (err = srt_f_interp_model(mdl_lbl, 
						&(undint->mdl_type), &(undint->mdl_dim) ) )
					return err;
				strcpy(undint->mdl_lbl,mdl_lbl) ;
			}
			else
			{
				undint->mdl_type = NONE ;
				strcpy(undint->mdl_lbl,"NONE") ;
			}
			undint->ts		= ts		;
			strcpy(undint->yc_name,crv_name1) ;
			und->spec_desc		= undint ;
		break ;            

		case EQUITY_UND :
			undeq = (SrtEqDesc *) srt_calloc(1,sizeof(SrtEqDesc)) ;
			if (mdl_lbl!=0)
			{
				if (err = srt_f_interp_model(mdl_lbl, 
							&(undeq->mdl_type), &(undeq->mdl_dim)))
					return err;
				strcpy(undeq->mdl_lbl,mdl_lbl) ;
			}
			else
			{
				undeq->mdl_type = NONE ;
				strcpy(undeq->mdl_lbl,"NONE") ;
			}
			undeq->ts		= ts		;
			undeq->spot		= spot		;
			strcpy(undeq->disc_name,crv_name1) ;
			strcpy(undeq->dvd_name,crv_name2) ;
			strcpy(undeq->repo_name,crv_name3) ;
			und->spec_desc		= undeq ;
		break ;              

		case FOREX_UND :	
			undfx = (SrtFxDesc *) srt_calloc(1,sizeof(SrtFxDesc)) ;
			if (mdl_lbl!=0)
			{
				if (err = srt_f_interp_model(mdl_lbl, 
						&(undfx->mdl_type),&(undfx->mdl_dim)))
					return err;
				strcpy(undfx->mdl_lbl,mdl_lbl) ;
			}
			else
			{
				undfx->mdl_type = NONE ;
				strcpy(undfx->mdl_lbl,"NONE") ;
			}
			undfx->ts		= ts		;
			undfx->spot		= spot		;
			strcpy(undfx->dom_name,crv_name1) ;
			strcpy(undfx->for_name,crv_name2) ;
		
			und->spec_desc		= undfx ;
		break ;

		case BOND_UND :	
		default :
			return serror("Unknown underlying label %s",und_lbl);
		break ;
	}

/* Insert (overwrite & update ticker) the object in the underlying list (with the same name) */
	err = srt_f_lstins(und_list, und_name, 0.0, OBJ_PTR_UND, 
		und, &srt_f_undvalfree, &(und->underl_ticker));
	if (err)
		return (serror("Error in initialising %s object",und->underl_lbl)) ;

/* Return a success message */
	return NULL ;

} /* END Err srt_f_addundtolist(...) */

/* --------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------
                      FUNCTION TO LOOKUP AN UNDERLYING OR DESTROY IT                
   --------------------------------------------------------------------------------- */

/* Get the Underlying  corresponding to the name in the underlying list used */
SrtUndPtr srt_f_lookup_und(String und_name)
{
Err             err ;
String          copy_name;
int             len;
SrtUndListPtr   und_list;
SrtUndPtr       undptr = NULL;

/* Get the underlying list */
	und_list = get_underlying_list();

/* Make a copy of the string not to modify it when UpperCasing it */	
	len = strlen(und_name);
	copy_name = (char *) srt_malloc(sizeof(char) * ( len + 1) );
	strncpy(copy_name, und_name, len);
	copy_name[len] = '\0';

/* Remove any extra information from the name and make it uppercase */
	strupper(copy_name);
	strip_white_space(copy_name);
	rem_tick_string(copy_name,copy_name);

/* Get the underlying object in the list (if it is there) */
	err = srt_f_getundinlist(und_list, copy_name, &undptr);
		
/* Free the copy_name memory */
	srt_free(copy_name);

/* Return whatever has to be returned */
	if (err)
	{
		return NULL ;
	}
	else
	{
		return undptr ;
	}

} /* END SrtUndPtr srt_f_lookup_und(...) */

/* Get the FX Underlying  corresponding to the domname and forname in the underlying list used */
SrtUndPtr srt_f_lookup_fxund(String dom_name, String for_name)
{
Err             err ;
String          copy_dom_name;
String          copy_for_name;
String          temp_dom;
String          temp_for;
int             len;
SrtUndListPtr   und_list;
SrtUndPtr       tempund = NULL;
SrtListObject   *obj;
SrtLst			*ls;

/* Get the underlying list */
	und_list = get_underlying_list();

/* Make a copy of the string not to modify it when UpperCasing it */	
	/*dom*/
	len = strlen(dom_name);
	copy_dom_name = (char *) srt_malloc(sizeof(char) * ( len + 1) );
	strncpy(copy_dom_name, dom_name, len);
	copy_dom_name[len] = '\0';
	/*for*/
	len = strlen(for_name);
	copy_for_name = (char *) srt_malloc(sizeof(char) * ( len + 1) );
	strncpy(copy_for_name, for_name, len);
	copy_for_name[len] = '\0';

/* Remove any extra information from the name and make it uppercase */
	/*dom*/
	strupper(copy_dom_name);
	strip_white_space(copy_dom_name);
	rem_tick_string(copy_dom_name,copy_dom_name);
	/*for*/
	strupper(copy_for_name);
	strip_white_space(copy_for_name);
	rem_tick_string(copy_for_name,copy_for_name);

/* Get the underlying object in the list (if it is there) */

/* Check that the list is not empty */
	if ( (und_list->head == und_list->tail) && (und_list->head->element== NULL) ) 
	{
		srt_free(copy_dom_name);
		srt_free(copy_for_name);
		return NULL;
	}

/* Go though the list to search for the object name */
	ls = und_list->head;
	while (ls!=NULL)
	{
		if (ls->element->type == OBJ_PTR_UND)
 		{
			obj = ls->element;
			tempund =  (SrtUndPtr) (obj->val.pval);

			if (tempund->underl_type == FOREX_UND)
			{
				err = get_fx_underlying_currencies(tempund, &temp_dom, &temp_for);
				if (err)
				{
					srt_free(copy_dom_name);
					srt_free(copy_for_name);
					return NULL;
				}
				
				if (strcmp(temp_dom, copy_dom_name) == 0 && strcmp(temp_for, copy_for_name) == 0)
				{
					srt_free(copy_dom_name);
					srt_free(copy_for_name);
					return tempund;
				}
			}
		}
		ls = ls->next;
	}
		
/* Free the copy_name memory */
	srt_free(copy_dom_name);
	srt_free(copy_for_name);
/* Return whatever has to be returned */
	return NULL;

} /* END SrtUndPtr srt_f_lookup_fxund(...) */

/* ------------------------------------------------------------------------- */
/* Destroy the underlying corresponding to the UndName in the underlying list */
Err srt_f_destroy_und (String undName)
{
SrtUndPtr       undptr;
Err             err   = NULL;

SrtUndListPtr   und_list;

/* Get THE underlying list */
	und_list = get_underlying_list();

/* Get the underlying in the list */
	undptr = srt_f_lookup_und (undName);
	
/* If the underlying does not exist: nothing to do */
	if (!undptr)
	{
		return NULL;
	}

/* Removes and frees the object from the list */
	err = srt_f_lstremobj( und_list, undName, 0.0);
	if (err)
		return err;

/* Return a success message */
	return NULL;

} /* END Err srt_f_destroy_und (...) */

/* ------------------------------------------------------------------------------- */


