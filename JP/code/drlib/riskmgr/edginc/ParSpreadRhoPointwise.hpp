/**
 * @file ParSpreadRhoPointwise.hpp
 */

#ifndef PS_RHOPOINTWISE_HPP
#define PS_RHOPOINTWISE_HPP

#include "edginc/PerNameRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Opaque constructor for ParSpreadRhoPointwise
 *
 * See ParSpreadRhoPointwise.cpp.
 *
 * @deprecated This is provided purely for the benefit of FirmAsset, which for
 * some reason wants to make its own ParSpreadRhoParallel and friends.
 */

RISKMGR_DLL PerNameRiskPropertySensitivity<ExpiryWindow>* newParSpreadRhoPointwise(double shiftSize);

DRLIB_END_NAMESPACE

#endif
