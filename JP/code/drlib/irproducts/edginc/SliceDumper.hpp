#ifndef SLICE_DUMPER_HPP
#define SLICE_DUMPER_HPP

#include "edginc/TreeSliceOper.hpp"
#include <fstream>
#include <string>

DRLIB_BEGIN_NAMESPACE

class SliceDumper : public SliceMarker<SliceDumper> {
public:
    SliceDumper(TreeSliceSP slice, string fileName) : slice(slice), fileName(fileName) {}
    
    //const int sliceCount() const { return 1; }
	static const int sliceCount = 1;
    
    template< typename S >
        const S** listSlices(const S** list) const {
            return slice->listSlices(list);
    }

    inline double calc() const {
        ofstream out;
        out.open(fileName.c_str(), ios_base::app);
        out << slice->calc() << endl;
        out.close();

        return 0;
    }

    void printDebug(char *s) const {
        strcat(s, "(SliceDumper)");
    }
    
private:
    TreeSliceSP slice;
    string fileName;
};

DRLIB_END_NAMESPACE
#endif
