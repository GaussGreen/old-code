/* ==============================================================
   FILENAME:      utlist.c

   PURPOSE:       All the functions to build a double linked
                  list that allows the storage of any
                                  information (int, double, string, pointers...)
   ============================================================== */
#include "utlist.h"

#include "math.h"
#include "utallhdr.h"
#include "utliststruct.h"

/* ------------------------------------------------------------
            TOP LEVEL FUNCTIONS TO WORK ON THE LIST
   ------------------------------------------------------------ */

/* ------------------------------------------------------------------------- */

Err srt_f_lstcreate(SrtList** Ptr, char* name)
{
    SrtList* li;
    SrtLst*  ls;

    li = (SrtList*)srt_calloc(1, sizeof(SrtList));
    if (li == NULL)
    {
        *Ptr = NULL;
        return serror("Memory allocation error (1) in srt_f_lst_create");
    }

    strncpy(li->name, name, 32);
    ls = (SrtLst*)srt_calloc(1, sizeof(SrtLst));
    if (!ls)
        return serror("Memory allocation error (2) in srt_f_lst_create");
    ls->next     = NULL;
    ls->previous = NULL;
    ls->element  = NULL;

    li->tail = ls;
    li->head = ls;

    *Ptr = li;
    return NULL;
} /* END Err srt_f_lstcreate(...) */

/* ------------------------------------------------------------------------- */
/* Insert an object in a list: update ticker and delete previous if necessary */
Err srt_f_lstsetobj(SrtList* l, SrtObject* obj, long* ticker)
{
    SrtLst *ls, *lt;
    Err     err = 0;

    if (!l)
        return serror("List does not exist in srt_f_lstsetobj");

    /* By default, the version number is 1 */
    *ticker = 1;

    /* The list is empty */
    if ((l->tail == l->head) && (l->tail->element == NULL))
    {
        l->tail->element = obj;
        return err;
    }

    /* Attaches the Object as the element of the SrtLst */
    ls = (SrtLst*)srt_calloc(1, sizeof(SrtLst));
    if (!ls)
        return serror("Memory allocation (1) in srt_f_lstsetobj");
    ls->element = obj;

    /* Move along the SrtLst until a higher element is found */
    lt = l->head;
    while ((lt != NULL) && (srt_f_objcmp(obj, (lt->element)) > 0))
        lt = lt->next;

    /* Element must be inserted at the end of the list */
    if (lt == NULL)
    {
        l->tail->next = ls;
        ls->next      = NULL;
        ls->previous  = l->tail;
        l->tail       = ls;
        return err;
    }

    /* Check if the object already exists in the list */
    if (srt_f_objequal(lt->element, obj) == SRT_YES)
    {
        /* Increment by 1 the SrtObject ticker to keep track of the Object version */
        obj->ticker = lt->element->ticker + 1;
        *ticker     = obj->ticker;

        if (lt->element != obj)
            err = srt_f_objfree(lt->element);

        srt_free(ls);
        lt->element = obj;
        return err;
    }

    /* Element must be inserted at the beginning of the list */
    if (lt == l->head)
    {
        ls->previous = NULL;
        ls->next     = lt;
        lt->previous = ls;
        l->head      = ls;
        return err;
    }

    ls->previous       = lt->previous;
    ls->next           = lt;
    lt->previous       = ls;
    ls->previous->next = ls;
    return err;

} /* END Err srt_f_lstsetobj( ... ) */

/* ------------------------------------------------------------------------- */
/* Remove an element in the list: free it and re-establish connectivity  */
Err srt_f_lstremobj(SrtList* l, char* objname, double objkey)
{
    SrtLst *   ls, *lt;
    Err        err = NULL;
    SrtObject* obj;

    /* Check that the list exists */
    if (!l)
        return serror("List does not exist in srt_f_lstset");

    /* If the list is empty, nothing to do */
    if ((l->tail == l->head) && (l->tail->element == NULL))
    {
        return serror("Cannot delete %s object: the list is empty", objname);
    }

    /* Go though the list to search for the object name */
    ls         = l->head;
    l->current = NULL;
    while (ls != NULL)
    {
        if ((strcmp(ls->element->name, objname) == 0) &&
            (fabs(ls->element->key - objkey) <= SRT_EPSILON))
        {
            obj        = ls->element;
            l->current = ls;
        }
        ls = ls->next;
    }

    /* Check that l->current points to something */
    if (l->current == NULL)
        return serror("Cannot delete %s object: not in the list", objname);

    /* Free the Object itself (and the pointer to the SrtUndDesc ) */
    err = srt_f_objfree(obj);
    if (err)
    {
        return err;
    }
    obj = NULL;

    /* There is only one element that has just been deleted */
    if (l->head == l->tail)
    {
        l->current = NULL;
        return err;
    }

    /* Restablish connectivity with the following node in the SrtList */
    if (l->current->next)
    {
        l->current->next->previous = l->current->previous;
    }
    else
    {
        l->tail = l->current->previous;
    }

    /* Restablish connectivity with the previous node in the SrtList */
    if (l->current->previous)
    {
        l->current->previous->next = l->current->next;
        lt                         = l->current->previous;
    }
    else
    {
        l->head = l->current->next;
        lt      = l->current->next;
    }

    /* Free the current SrtLstAtom */
    srt_free(l->current);

    /* Make current point to the closest SrtListAtom (lt) */
    l->current = lt;

    /* Return a success message */
    return err;

} /* END Err srt_f_lstremobj( ... ) */

/* ------------------------------------------------------------------------- */
/* Create an object and insert it in the list from its components (ticker update) */
Err srt_f_lstins(
    SrtList* l,
    char*    name,
    double   key,
    int      type,
    void*    Ptr,
    Err (*ObjValPvalFreeFct)(void*),
    long* ticker)
{
    SrtObject* objPtr;
    SrtObjVal* objval;
    Err        err = 0;

    /* If list does not exist, return an error message */
    if (!l)
        return serror("List does not exist in srt_f_lstins");

    /* Allocate memory for a SrtObjVal and set its content */
    switch (type)
    {
    case OBJ_INT:
        objval = i2obj(*(int*)Ptr);
        break;
    case OBJ_DOUBLE:
        objval = d2obj(*(double*)Ptr);
        break;
    case OBJ_STRING:
        objval = s2obj((char*)Ptr);
        break;
    default:
        objval = val2obj(Ptr);
        break;
    }

    /* Allocate memory and fills in the SrtObjectPtr  */
    srt_f_objset(name, key, type, objval, &objPtr, SRT_YES, ObjValPvalFreeFct);

    /* Attaches the SrtObject to a (new) element in the list */
    srt_f_lstsetobj(l, objPtr, ticker);

    return err;

} /* END Err srt_f_lstins(...) */

/* ------------------------------------------------------------------------- */

/* Get an object in a list through its name AND its key */
Err srt_f_lstgetobj(SrtList l, char* objname, double key, SrtObject** ret)
{
    SrtLst* ls;
    Err     err = 0;

    /* Check that the list is not empty */
    if ((l.head == l.tail) && (l.head->element == NULL))
    {
        *ret = NULL;
        return serror("Object %s not found: list is empty", objname);
    }

    /* Go though the list to search for the object name */
    ls = l.head;
    while (ls != NULL)
    {
        if ((strcmp(ls->element->name, objname) == 0) &&
            (fabs(ls->element->key - key) <= SRT_EPSILON))
        {
            *ret = ls->element;
            return NULL;
        }
        ls = ls->next;
    }

    /* The object is not in the list...*/
    *ret = NULL;
    return serror("Object %s not found", objname);

} /* END Err srt_f_lstgetobj(...) */

/* ------------------------------------------------------------------------- */

SRT_Boolean srt_f_lsthas(SrtList l, char name[], double key)
{
    SrtObject* ret;

    if (srt_f_lstgetobj(l, name, key, &ret) == NULL)
    {
        return SRT_YES;
    }

    return SRT_NO;

} /* END SRT_Boolean srt_f_lsthas(...) */

/* ------------------------------------------------------------------------- */

Err srt_f_lstfree(SrtList* l, int free_obj)
{
    SrtLst *ls, *lt;
    Err     err = 0;

    /* If list does not exist, no need to free it */
    if (!l)
        return NULL;

    ls = l->head;
    while (ls != NULL)
    {
        lt = ls;
        ls = lt->next;
        if (free_obj == SRT_YES && lt->element != NULL)
        {
            err         = srt_f_objfree(lt->element);
            lt->element = NULL;
        }
        free(lt);
        lt = NULL;
    } /* END of loop on element of l list */

    /* OVE + DEN: DO not free here: free after lstfree
            srt_free(l);
    */

    l->head    = NULL;
    l->tail    = NULL;
    l->current = NULL;

    return err;

} /* END Err srt_f_lstfree(...) */

/* ------------------------------------------------------------------------- */

/* ------------------------------------------------------------
          FUNCTIONS TO WORK WITH AN OBJECT OF THE LIST
   ------------------------------------------------------------ */

int srt_f_objcmp(SrtObject* obj1, SrtObject* obj2)
{
    if (fabs(obj1->key - obj2->key) <= SRT_EPSILON)
        return (strcmp(obj1->name, obj2->name));

    return (obj1->key > obj2->key) ? 1 : -1;
}

/* ------------------------------------------------------------------------- */

SRT_Boolean srt_f_objequal(SrtObject* obj1, SrtObject* obj2)
{
    int res;
    res = (((strcmp(obj1->name, obj2->name) == 0) ? 1 : 0) && (obj1->key == obj2->key));
    if (res == 1)
        return SRT_YES;
    else
        return SRT_NO;
}

/* 2 objects are equal if name1 = name2 and key1 = key2 */

/* ------------------------------------------------------------------------------ */

/* Create a new SrtObject and fill it in with the input values  */
Err srt_f_objset(
    char*       name,
    double      key,
    int         type,
    SrtObjVal*  val,
    SrtObject** Ptr,
    SRT_Boolean ObjValPtr_allocated,
    Err (*ObjValPvalFreeFct)(void*))
{
    SrtObject* obj;
    Err        err = 0;

    /* Allocate memory for a SrtObject */
    obj = (SrtObject*)srt_malloc(sizeof(SrtObject));
    if (!obj)
    {
        *Ptr = NULL;
        return serror("Allocation error in srt_f_objset");
    }

    /* Set the name of the object */
    strncpy(obj->name, name, 32);

    /* Set the SrtObject version ticker to 1 (first one) */
    obj->ticker = 1;

    /* Copy the inputs as members of the object */
    obj->key                 = key;
    obj->type                = type;
    obj->val                 = *val;
    obj->ObjValPtr_allocated = ObjValPtr_allocated;
    obj->ObjValPvalFreeFct   = ObjValPvalFreeFct;
    srt_free(val);
    *Ptr = obj;

    return err;

} /* END Err srt_f_objset(...) */

/* ------------------------------------------------------------------------- */

/* Frees the Object and everything it contains */
Err srt_f_objfree(SrtObject* obj)
{
    Err err = 0;

    /* If object does not exist, no need to free it */
    if (!obj)
        return NULL;

    /* Check on the type of the SrtObjVal to see if a val.pval exists */
    switch (obj->type)
    {
        /* Here, the use of a free function is required */
    case (OBJ_PTR_IRM_TermStruct):
    case (OBJ_PTR_IRM2f_TermStruct):
    case (OBJ_PTR_FX_TermStruct):
    case (OBJ_PTR_EQ_TermStruct):
    case (OBJ_PTR_IO):
    case (OBJ_PTR_UND):
    case (OBJ_PTR_CorrelStruct):
    case (OBJ_PTR_HIST):
    case (OBJ_PTR_CURVE):
        /* Call the FreePVal function associated with the object */
        err = obj->ObjValPvalFreeFct(obj->val.pval);
        if (err)
            return err;
        break;
        /* In the other case, nothing to do */
    default:
        break;
    }

    /* Free the SrtObject itself */
    srt_free(obj);

    /* Return a success message */
    return err;

} /* END Err srt_f_objfree(...) */

/* ------------------------------------------------------------------------- */

/* ------------------------------------------------------------
          FUNCTIONS TO CREATE AN OBJECTVAL FOR AN OBJECT
   ------------------------------------------------------------ */

SrtObjVal* d2obj(double d)
{
    SrtObjVal* ov;
    ov       = (SrtObjVal*)srt_malloc(sizeof(SrtObjVal));
    ov->dval = d;
    return ov;
}

/* ------------------------------------------------------------------------- */

SrtObjVal* i2obj(int i)
{
    SrtObjVal* ov;
    ov       = (SrtObjVal*)srt_malloc(sizeof(SrtObjVal));
    ov->ival = i;
    return ov;
}

/* ------------------------------------------------------------------------- */

SrtObjVal* s2obj(char s[])
{
    SrtObjVal* ov;
    ov = (SrtObjVal*)srt_malloc(sizeof(SrtObjVal));
    strncpy(ov->sval, s, 32);
    return ov;
}

/* ------------------------------------------------------------------------- */

SrtObjVal* val2obj(void* ptr)
{
    SrtObjVal* ov;
    ov       = (SrtObjVal*)srt_malloc(sizeof(SrtObjVal));
    ov->pval = ptr;
    return (ov);
}

/*---------------------------------

SrtObjVal* histdata2obj(SrtHistData* shd)
{
        SrtObjVal* ov;
        ov = (SrtObjVal*) srt_malloc (sizeof(SrtObjVal));
        ov->pval = shd;
        return (ov);
}
                        -------------------------------------- */

/* ------------------------------------------------------------------------- */