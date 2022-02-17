/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/version.h
// Purpose:     helpers for version tracking, used by freeway/version.h
// Author:      Vadim Zeitlin
// Created:     25.04.03
// RCS-ID:      $Id: version.h,v 1.7 2004/10/05 09:13:35 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/version.h
    @brief Helper macros for defining version strings.

    This tiny helper header defines common macros used in all
    ito33/product/version.h files (for all values of product). It is not
    supposed to be included directly by the user code, only by the other
    version.h files.
 */

#ifndef _ITO33_VERSION_H_
#define _ITO33_VERSION_H_

#ifndef ITO33_REALLY_NEED_VERSION_H
  #error "Please don't include this file directly."
#else
  #undef ITO33_REALLY_NEED_VERSION_H
#endif

/// helper for ITO33_MAKE_VERSION_STRING
#define ITO33_STR(x) #x

/**
    Produce a collated version string from the version components.

    For example, for the version 2.3.5 the returned value is "235".

    @param x major version component
    @param y minor version component
    @param z release number
    @return string containing the concatenation of the parameters
 */
#define ITO33_MAKE_VERSION_STRING(x, y, z) \
  ITO33_STR(x) ITO33_STR(y) ITO33_STR(z)

/**
    Produce a series string with dot from major the version components
    and minor version component.

    For example, for the version 2.3 the returned value is "2.3".
    This is mainly useful for the RC files which must embed the version info in
    this format.

    @param x major version component
    @param y minor version component
    @return string containing the concatenation of the parameters
 */
#define ITO33_MAKE_SERIES_DOT_STRING(x, y) \
  ITO33_STR(x) "." ITO33_STR(y)

/**
    Produce a version string with dots from the version components.

    For example, for the version 2.3.5 the returned value is "2.3.5".
    This is mainly useful for the RC files which must embed the version info in
    this format.

    @param x major version component
    @param y minor version component
    @param z release number
    @return string containing the concatenation of the parameters
 */
#define ITO33_MAKE_VERSION_DOT_STRING(x, y, z) \
  ITO33_STR(x) "." ITO33_STR(y) "." ITO33_STR(z)

/**
    Produce full version string with dots from the version components.

    For example, for the version 2.3.5.10 the returned value is "2.3.5.10".
    This is mainly useful for the RC files which must embed the version info in
    this format.

    @param x major version component
    @param y minor version component
    @param z release number
    @param b build number
    @return string containing the concatenation of the parameters
 */
#define ITO33_MAKE_FULLVERSION_DOT_STRING(x, y, z, b) \
  ITO33_STR(x) "." ITO33_STR(y) "." ITO33_STR(z) "." ITO33_STR(b)

/**
    Checks if the current product version is equal or greater than the given
    one.

    This macro allows to test for the current version in #if preprocessor
    statements and so write code working with different, possibly incompatible
    versions of the same product.

    Note that it is never used directly from the user code, please use
    ITO33_PRODUCT_CHECK_VERSION instead.

    @param product product name (the macros ITO33_PRODUCT_VERSION_MAJOR, MINOR
                   and RELEASE are supposed to be defined)
    @param x major component of the version to compare with
    @param y minor version component
    @param z release number
    @return true if x.y.z >= product_major.product_minor.product_release
 */
#define ITO33_COMPARE_VERSIONS(product, x, y, z) \
  (ITO33_ ## product ## _VERSION_MAJOR > (x) || \
  (ITO33_ ## product ## _VERSION_MAJOR == (x) && \
    ITO33_ ## product ## _VERSION_MINOR > (y)) || \
  (ITO33_ ## product ## _VERSION_MAJOR == (x) && \
    ITO33_ ## product ## _VERSION_MINOR == (y) && \
      ITO33_ ## product ## _VERSION_RELEASE >= (z)))

#endif // _ITO33_VERSION_H_

