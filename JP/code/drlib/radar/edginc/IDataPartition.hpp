#ifndef IDATA_PARTITION_H
#define IDATA_PARTITION_H

#if defined(_MSC_VER)
#pragma warning(disable:4786)
#endif

#include "edginc/DECLARE.hpp"
#include "edginc/Format.hpp"
#include "edginc/RadarRepUtil.hpp"

#include <vector>
#include <map>
#include <numeric>

DRLIB_BEGIN_NAMESPACE

using namespace std;

// a class which partitions the fitting variables (before transformations have been applied?)
class RADAR_DLL IDataPartition {
public:
	// input:	a vector of fitting variables
	// output:	the corresponding element of the partition
	virtual size_t classify( const FittingArray & x ) const = 0;

	// get the total number of elements in the partition
	virtual size_t getSize() const = 0;
	virtual ~IDataPartition() {}
};

DECLARE_REF_COUNT(IDataPartition);

DRLIB_END_NAMESPACE
#endif
