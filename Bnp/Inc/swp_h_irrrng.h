#ifndef SWP_H_IRRRNG_H
#define SWP_H_IRRRNG_H

Err irr_range(
    Date*          d,
    double*        cf,
    int            len,
    DateListType   list_type,
    Date           prev_date,
    double         pv,
    SrtBasisCode   basis_code,
    SrtCompounding compounding,
    double*        irr);

Err pv_range_with_irr(
    Date*          d,
    double*        cf,
    int            len,
    DateListType   list_type,
    Date           prev_date,
    SrtBasisCode   basis_code,
    SrtCompounding compounding,
    double         irr,
    double*        pv);

Err duration_range(
    Date*          d,
    double*        cf,
    int            len,
    DateListType   list_type,
    Date           prev_date,
    SrtBasisCode   basis_code,
    SrtCompounding compounding,
    double         irr,
    double*        duration);

Err convexity_range(
    Date*          d,
    double*        cf,
    int            len,
    DateListType   list_type,
    Date           prev_date,
    SrtBasisCode   basis_code,
    SrtCompounding compounding,
    double         irr,
    double*        convexity);

Err thirdmoment_range(
    Date*          d,
    double*        cf,
    int            len,
    DateListType   list_type,
    Date           prev_date,
    SrtBasisCode   basis_code,
    SrtCompounding compounding,
    double         irr,
    double*        thirdmoment);

Err modified_duration_range(
    Date*          d,
    double*        cf,
    int            len,
    DateListType   list_type,
    Date           prev_date,
    SrtBasisCode   basis_code,
    SrtCompounding compounding,
    double         irr,
    double*        mod_duration);

Err modified_convexity_range(
    Date*          d,
    double*        cf,
    int            len,
    DateListType   list_type,
    Date           prev_date,
    SrtBasisCode   basis_code,
    SrtCompounding compounding,
    double         irr,
    double*        mod_convexity);

#endif
