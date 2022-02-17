#include "Util.h"

#include <fstream>

static std::ofstream logger("c:\\temp\\log.txt", std::ios::out);

void dump(const Array<double>& data, char* prefix)
{
#ifdef _DEBUG
    logger << prefix << "\n";
    for (int i = 0; i < data.size(); ++i)
        logger << i << " " << data[i] << "\n";
    logger << std::flush;
#endif
}

void dump(const Array<int>& data, char* prefix)
{
#ifdef _DEBUG
    logger << prefix << "\n";
    for (int i = 0; i < data.size(); ++i)
        logger << i << " " << data[i] << "\n";
    logger << std::flush;
#endif
}

void dump(double num, char* prefix)
{
#ifdef _DEBUG
    logger << prefix << " " << num << "\n";
#endif
}


void dump(const Matrix<double>& data, int rowId, int columnId, char* prefix)
{
#ifdef _DEBUG
    logger << prefix << "\n";
    if (rowId < 0 && columnId < 0) {
        int numRows = data.rowSize();
        int numCols = data.colSize();
        for (int i = 0; i < numRows; ++i) {
            for (int j = 0; j < numCols; ++j)
                logger << data[i][j] << " ";
            logger << "\n";
        }
        logger << std::flush;
    } else if (rowId < 0) {  
        int numRows = data.rowSize();
        for (int i = 0; i < numRows; ++i)
            logger << i << " " << data[i][columnId] << "\n";
        logger << std::flush;
    } else if (columnId < 0) {
        int numCols = data.colSize();
        for (int i = 0; i < numCols; ++i)
            logger << i << " " << data[rowId][i] << "\n";
        logger << std::flush;
    }
#endif
}

void err(char* prefix)
{
    logger << prefix << std::flush;
}

void err(double num)
{
    logger << num << std::flush;
}
