
#include <windows.h>
#include "xlcall.h"
#include "framewrk.h"

void cwCenter(HWND, int);
INT_PTR CALLBACK DIALOGMsgProc(HWND hWndDlg, UINT message, WPARAM wParam, LPARAM lParam);
BOOL GetHwnd(HWND * pHwnd);
int lpwstricmp(LPWSTR s, LPWSTR t);

//
// identifier for controls
//
#define FREE_SPACE                  104
#define EDIT                        101
#define TEST_EDIT                   106
