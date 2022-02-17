#include "Magnet/Magnet.h"
#include "General/General.h"

using namespace CM;

void dump(const Array<double>& data, char* prefix);

void dump(const Array<int>& data, char* prefix);

void dump(const Matrix<double>& data, int rowId, int columnId, char* prefix);

void dump(double num, char* prefix);

void err(char* prefix);

void err(double num);
