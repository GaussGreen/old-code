/**********************************************************************
 *      Name: srt_f_iostruct.c                                        *
 *  Function:                                                         *
 * Copyright: (C) Paribas Capital Markets Ltd.                        *
 *--------------------------------------------------------------------*
 *    Author: Arnaud SAHAGUET                                         *
 *      Date: 19/12/94                                                *
 *--------------------------------------------------------------------*
 *    Inputs:                                                         *
 *   Returns:                                                         *
 *   Globals:                                                         *
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 **********************************************************************/

#include "srt_h_all.h"

/* -------------------------------------------------------------------------- */
/* ------------------- Type number to Type name function  ------------------- */
/* -------------------------------------------------------------------------- */

Err srt_f_IOtype2string(int type, char name[32])
{
    Err err = NULL;

    strncpy(name, "", 32);
    switch (type)
    {
    case IO_PREMIUM:
        strncpy(&(name[0]), "PREM", strlen("PREM"));
        break;
    case IO_RATE_SHIFT:
        strncpy(&(name[0]), "R_SH", strlen("R_SH"));
        break;
    case IO_SIGMA_SHIFT:
        strncpy(&(name[0]), "S_SH", strlen("S_SH"));
        break;
    case IO_TAU_SHIFT:
        strncpy(&(name[0]), "T_SH", strlen("T_SH"));
        break;
    case IO_STDEV:
        strncpy(&(name[0]), "STDEV", strlen("STDEV"));
        break;
    case IO_COLPVS:
        strncpy(&(name[0]), "COLPVS", strlen("COLPVS"));
        break;
    default:
        strncpy(&(name[0]), "", 0);
        break;
    }
    return err;
}

/* -------------------------------------------------------------------------- */
/* ----------------------- Name encoding function --------------------------- */
/* -------------------------------------------------------------------------- */

Err srt_f_IOencodename(
    char id_name[32], int type, int bucket_no, double shift_value, int shift_type, char name[32])
{
    char name_type[32];
    Err  err = 0;

    srt_f_IOtype2string(type, name_type);
    sprintf(
        name,
        "%s-%.1d-%.5lf-%d-%s",
        name_type, /* PREMIUM, TAU, ... */
        bucket_no,
        shift_value,
        shift_type,
        id_name);

    return err;
}

/* -------------------------------------------------------------------------- */
/* ----------------------------- show function ------------------------------ */
/* -------------------------------------------------------------------------- */

Err srt_f_IOvalshow(SrtIOVal* ioval, FILE* file_out)
{
    Err  err = 0;
    int  type;
    char name_type[32];

    type = ioval->type;
    srt_f_IOtype2string(type, name_type);
    fprintf(file_out, "Result type = %s\n", name_type);
    if ((type == IO_SIGMA_SHIFT) || (type == IO_TAU_SHIFT))
        fprintf(file_out, "Bucket start =%ld \n", ioval->bucket_start);
    fprintf(file_out, "Bucket end =%ld \n", ioval->bucket_end);

    if ((type == IO_SIGMA_SHIFT) || (type == IO_TAU_SHIFT) || (type == IO_RATE_SHIFT))
    {
        if (ioval->shift_type == SH_ABSOLUTE)
        {
            fprintf(file_out, "Absolute ");
        }
        else
        {
            fprintf(file_out, "Relative ");
        }
        fprintf(file_out, "shift = %f \n", ioval->shift_value);
    }
    if (ioval->done == SRT_YES)
    {
        fprintf(file_out, "Job done\n");
    }
    else
    {
        fprintf(file_out, "Job not done...\n");
    }
    fprintf(file_out, "Computed Value = %.10lf\n", ioval->dval);
    fprintf(
        file_out,
        "Function use for computation : %s\n",
        ((strcmp(ioval->computation_origin, "") == 0) ? "None" : ioval->computation_origin));
    fprintf(file_out, "ID_NAME = %s\n", ioval->id_name);
    if (ioval->pval == NULL)
    {
        fprintf(file_out, "PTR NULL\n");
    }
    else
    {
        fprintf(file_out, "PTR not NULL\n");
    }
    return err;
}

/* -------------------------------------------------------------------------- */
/* ----------------------- free function for IOVal -------------------------- */
/* -------------------------------------------------------------------------- */

Err srt_f_IOvalfree(void* iovalptr)
{
    Err       err   = NULL;
    SrtIOVal* ioval = (SrtIOVal*)iovalptr;

    if (ioval != NULL)
    {
        if (ioval->type == IO_COLPVS)
        {
            if (ioval->pval != NULL)
                free_dvector((double*)ioval->pval, 0, ioval->lval - 1);
            ioval->pval = NULL;
        }
        srt_free(ioval);
    }
    else
    {
        return "Srt_f_iovalfree() failed";
    }
    return err;
}

/* -------------------------------------------------------------------------- */
/* ----------------------- free function for IOStruct ----------------------- */
/* -------------------------------------------------------------------------- */

Err srt_f_IOstructfree(SrtIOStruct** l)
{
    Err err = NULL;

    if (err = srt_f_lstfree(*l, SRT_YES))
        return err;

    srt_free(*l);

    return err;
}

/* -------------------------------------------------------------------------- */
/* ------------------------ conversion function ----------------------------- */
/* -------------------------------------------------------------------------- */

SrtObjVal* srt_f_IOval2obj(SrtIOVal* ioval)
{
    SrtObjVal* ov;

    ov       = (SrtObjVal*)srt_malloc(sizeof(SrtObjVal));
    ov->pval = ioval;
    return (ov);
}

/* -------------------------------------------------------------------------- */
/* --------------------------- Generic set function ------------------------- */
/* -------------------------------------------------------------------------- */

Err srt_f_IOstructset(
    SrtIOStruct* l,
    char         id_name[32],
    int          type,
    int          bucket_no,
    Date         bucket_start,
    Date         bucket_end,
    double       shift_value,
    int          shift_type,
    SRT_Boolean  done,
    char         comp_origin[32],
    long         lval,
    double       dval,
    void*        pval)
{
    char      name[32];
    SrtIOVal* ioval;
    Err       err = 0;
    long      ticker;

    /* Allocate memory for the IOVAl */
    ioval = (SrtIOVal*)srt_malloc(sizeof(SrtIOVal));
    memset(ioval, 0, sizeof(SrtIOVal));

    /* Transfers all the input information into the SrtIOVal */
    ioval->type         = type;
    ioval->bucket_start = bucket_start;
    ioval->bucket_end   = bucket_end;
    ioval->shift_value  = shift_value;
    ioval->shift_type   = shift_type;
    ioval->done         = done;
    ioval->lval         = lval;
    ioval->dval         = dval;
    ioval->pval         = pval;

    /* Puts the computation origin string */
    strncpy(ioval->computation_origin, "", 32);
    strncpy(ioval->computation_origin, comp_origin, strlen(comp_origin));

    /* Puts the id name of the job in the IOVal */
    strncpy(ioval->id_name, "", 32);
    strncpy(ioval->id_name, id_name, strlen(id_name));

    /* Builds the name of the object to put into the IO list (using the ioval->type) */
    srt_f_IOencodename(
        id_name, ioval->type, bucket_no, ioval->shift_value, ioval->shift_type, name);

    /* Inserts the object (with its name) in the IO list */
    srt_f_lstins(l, name, ioval->shift_value, OBJ_PTR_IO, (void*)ioval, &srt_f_IOvalfree, &ticker);

    return err;
}

/* -------------------------------------------------------------------------- */
/* --------------------------- Generic get function ------------------------- */
/* -------------------------------------------------------------------------- */

Err srt_f_IOstructget(
    SrtIOStruct l,
    char        id_name[32],
    int         type,
    int         bucket_no,
    double      shift_value,
    int         shift_type,
    SrtIOVal**  ioval_pp)
{
    char       name[32];
    Err        err = 0;
    double     key;
    SrtObject* obj;

    srt_f_IOencodename(id_name, type, bucket_no, shift_value, shift_type, name);
    key = shift_value;

    err = srt_f_lstgetobj(l, name, key, &obj);
    if (err)
    {
        *ioval_pp = NULL;
        err       = "Error in srt_f_IOstructget(): object not found";
        return err;
    }
    else
        *ioval_pp = (SrtIOVal*)(obj->val.pval);
    return err;
}

/* -------------------------------------------------------------------------- */
/* --------------------------- Premium functions ---------------------------- */
/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgetpremiumval(SrtIOStruct l, double* premium)
{
    Err       err = 0;
    SrtIOVal* ioval_p;

    err = srt_f_IOstructgetpremium(l, &ioval_p);
    if (err)
    {
        *premium = 0.0;
    }
    else
    {
        *premium = ioval_p->dval;
    }
    return err;
}

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgetpremium(SrtIOStruct l, SrtIOVal** ioval_pp)
{
    Err err = 0;

    err = srt_f_IOstructget(l, "", IO_PREMIUM, 0, 0, 0, ioval_pp);
    return err;
}

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructsetpremium(SrtIOStruct* l, SRT_Boolean done, double result, char comp_origin[32])
{
    Err err = 0;

    err = srt_f_IOstructset(l, "", IO_PREMIUM, 0, 0, 0, 0, 0, done, comp_origin, 0, result, NULL);
    return err;
}

/* -------------------------------------------------------------------------- */
/* ------------------------ Standard Deviation functions -------------------- */
/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgetstdevval(SrtIOStruct l, double* stdev)
{
    Err       err = 0;
    SrtIOVal* ioval_s;

    err = srt_f_IOstructgetstdev(l, &ioval_s);
    if (err)
    {
        *stdev = 0.0;
    }
    else
    {
        *stdev = ioval_s->dval;
    }
    return err;
}

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgetstdev(SrtIOStruct l, SrtIOVal** ioval_ss)
{
    Err err = 0;

    err = srt_f_IOstructget(l, "", IO_STDEV, 0, 0, 0, ioval_ss);
    return err;
}

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructsetstdev(SrtIOStruct* l, SRT_Boolean done, double result, char comp_origin[32])
{
    Err err = 0;

    err = srt_f_IOstructset(l, "", IO_STDEV, 0, 0, 0, 0, 0, done, comp_origin, 0, result, NULL);
    return err;
}

/* -------------------------------------------------------------------------- */
/* Function to put the PV's of all the columns in the tableau */
/* -------------------------------------------------------------------------- */

Err srt_f_IOstructsetcolpvs(
    SrtIOStruct* l, SRT_Boolean done, double* pv_cols, long num_pv, char comp_origin[32])
{
    Err     err = 0;
    double* copy;

    /* Allocate space for the copy of the pvs */
    copy = dvector(0, num_pv - 1);
    memcpy(&copy[0], &pv_cols[0], num_pv * sizeof(double));

    /* Here, uses the lval field to store the dimension of the tableau to store */
    err = srt_f_IOstructset(
        l, "", IO_COLPVS, 0, 0, 0, 0, 0, done, comp_origin, num_pv, 0, (void*)copy);

    return err;
}

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgetcolpvs(SrtIOStruct l, double** pvs, long* num_pvs)
{
    Err       err = 0;
    SrtIOVal* ioval_s;

    err = srt_f_IOstructget(l, "", IO_COLPVS, 0, 0, 0, &ioval_s);
    if (err)
    {
        *pvs     = NULL;
        *num_pvs = 0;
    }
    else
    {
        *pvs     = (double*)ioval_s->pval;
        *num_pvs = ioval_s->lval;
    }
    return err;
}

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* FOR ALL SORTS OF GREEKS COMPUTATIONS BY A PARAMTER SHIFT */
/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* ----------------------------- rate functions ------------------------------ */
/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgetrshift(
    SrtIOStruct l, char id_name[32], double shift_value, int shift_type, SrtIOVal** ioval_pp)
{
    Err err = NULL;

    if (err = srt_f_IOstructget(l, "", IO_RATE_SHIFT, 0, shift_value, shift_type, ioval_pp))
    {
        return err;
    }

    return NULL;
}

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgetrshiftval(
    SrtIOStruct l, char id_name[32], double shift_value, int shift_type, double* value)
{
    Err       err = NULL;
    SrtIOVal* ioval_p;

    if (err = srt_f_IOstructgetrshift(l, id_name, shift_value, shift_type, &ioval_p))
    {
        return err;
    }

    *value = ioval_p->dval;

    return NULL;
}

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructsetrshift(
    SrtIOStruct* l,
    char         id_name[32],
    double       shift_value,
    int          shift_type,
    SRT_Boolean  done,
    double       value,
    char         comp_origin[32],
    void*        ptr)
{
    Err err = NULL;

    if (err = srt_f_IOstructset(
            l,
            id_name,
            IO_RATE_SHIFT,
            0,
            0,
            0,
            shift_value,
            shift_type,
            done,
            comp_origin,
            0,
            value,
            ptr))
    {
        return err;
    }

    return NULL;
}

/* -------------------------------------------------------------------------- */
/* ------------------------------- sigma functions -------------------------- */
/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgetsigmashift(
    SrtIOStruct l, int bucket_no, double shift_value, int shift_type, SrtIOVal** ioval_pp)
{
    Err err = NULL;

    if (err =
            srt_f_IOstructget(l, "", IO_SIGMA_SHIFT, bucket_no, shift_value, shift_type, ioval_pp))
    {
        return err;
    }

    return NULL;
}

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgetsigmashiftval(
    SrtIOStruct l, int bucket_no, double shift_value, int shift_type, double* value)
{
    Err       err = NULL;
    SrtIOVal* ioval_p;

    if (err = srt_f_IOstructgetsigmashift(l, bucket_no, shift_value, shift_type, &ioval_p))
    {
        return err;
    }

    *value = ioval_p->dval;

    return NULL;
}

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructsetsigmashift(
    SrtIOStruct* l,
    int          bucket_no,
    Date         bucket_start,
    Date         bucket_end,
    double       shift_value,
    int          shift_type,
    SRT_Boolean  done,
    double       value,
    char         comp_origin[32])
{
    Err err = NULL;

    if (err = srt_f_IOstructset(
            l,
            "",
            IO_SIGMA_SHIFT,
            bucket_no,
            bucket_start,
            bucket_end,
            shift_value,
            shift_type,
            done,
            comp_origin,
            0,
            value,
            NULL))
    {
        return err;
    }

    return NULL;
}

/* -------------------------------------------------------------------------- */
/* ---------------------------- tau functions ------------------------------- */
/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgettaushift(
    SrtIOStruct l, int bucket_no, double shift_value, int shift_type, SrtIOVal** ioval_pp)
{
    Err err = NULL;

    if (err = srt_f_IOstructget(l, "", IO_TAU_SHIFT, bucket_no, shift_value, shift_type, ioval_pp))
    {
        return err;
    }

    return NULL;
}

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructgettaushiftval(
    SrtIOStruct l, int bucket_no, double shift_value, int shift_type, double* value)
{
    Err       err = 0;
    SrtIOVal* ioval_p;

    err = srt_f_IOstructgettaushift(l, bucket_no, shift_value, shift_type, &ioval_p);
    if (err)
        *value = 0.0;
    else
        *value = ioval_p->dval;

    return err;
}

/* -------------------------------------------------------------------------- */

Err srt_f_IOstructsettaushift(
    SrtIOStruct* l,
    int          bucket_no,
    Date         bucket_start,
    Date         bucket_end,
    double       shift_value,
    int          shift_type,
    SRT_Boolean  done,
    double       value,
    char         comp_origin[32])
{
    Err err = 0;
    err     = srt_f_IOstructset(
        l,
        "",
        IO_TAU_SHIFT,
        bucket_no,
        bucket_start,
        bucket_end,
        shift_value,
        shift_type,
        done,
        comp_origin,
        0,
        value,
        NULL);
    return err;
}

/* -------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------- */
/* ------------------------- Create I/O list   ------------------------------ */
/* -------------------------------------------------------------------------- */

Err srt_f_IOstructcreate(SrtIOStruct** l, char* name)
{
    Err err = NULL;

    if (err = srt_f_lstcreate(l, name))
    {
        return err;
    }

    return NULL;
}
