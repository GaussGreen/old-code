/* =================================================================================
   FILENAME :     SrtInit.c

   PURPOSE:       The Function to call to set up the SORT library environment
   ================================================================================= */

#include "SrtAccess.h"

/* --------------------------------------------------------------------- */

char* SrtVersion()
{
    return "10.5.0.0";
}

/* --------------------------------------------------------------------- */

int SrtInit()
{
    char* err = NULL;

    /* Create an empty list to store underlyings */
    err = create_underlying_list("The Underlying List");
    if (err)
    {
        return 1;
    }

    /* Create an empty list to store products */
    err = create_product_list("The Product List");
    if (err)
    {
        return 1;
    }

    /* Create an empty list to store curves */
    err = create_curve_list("The Curve List");
    if (err)
    {
        return 1;
    }

    /* Create an empty list to store the Grfn requests */
    err = create_request_list("The Request List");
    if (err)
    {
        return 1;
    }

    /* Create an empty correlation list for Grfn */
    err = create_correlation_list("The Correlation List");
    if (err)
    {
        return 1;
    }

    /* Create an empty list to store Grfn histograms */
    err = create_histogram_list("The Histogram List");
    if (err)
    {
        return 1;
    }

    return 0;

} /* END int SrtInit(...) */

/* --------------------------------------------------------------------- */

int SrtClose()
{
    Err err = NULL;

    /* Destroy all the underlyings */
    err = destroy_all_underlyings();
    if (err)
    {
        return 1;
    }

    /* Destroy all the products */
    err = destroy_all_products();
    if (err)
    {
        return 1;
    }

    /* Destroy all the curves */
    err = destroy_all_curves();
    if (err)
    {
        return 1;
    }

    /* Destroy all the Grfn requests */
    err = destroy_request_list();
    if (err)
    {
        return 1;
    }

    /* Destroy the correlation list for Grfn multi-underlyings */
    err = destroy_correlation_list();
    if (err)
    {
        return 1;
    }

    /* Destroy all the histograms for Grfn */
    err = destroy_histogram_list();
    if (err)
    {
        return 1;
    }

    return 0;

} /* END int SrtClose(...) */

/* --------------------------------------------------------------------------------- */
