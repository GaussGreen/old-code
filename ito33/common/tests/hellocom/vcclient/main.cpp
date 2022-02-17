#include <stdio.h>

#import "../hello.tlb"

int main()
{
  CoInitialize(NULL);

  try
  {
    Ito33HelloLib::IHelloPtr ptr("Ito33HelloLib.Ito33Hello");

    ptr->SayHello("C++ client");
  }
  catch ( const _com_error& ce )
  {
    printf("COM exception: %#08x", ce.Error());

    IErrorInfo *pErrorInfo;
    if ( SUCCEEDED(::GetErrorInfo(0, &pErrorInfo)) && pErrorInfo )
    {
      BSTR bstr;
      pErrorInfo->GetDescription(&bstr);
      if ( *bstr )
        printf(":\n\n%s", (char *)_bstr_t(bstr));

      pErrorInfo->Release();
    }

    putchar('\n');
  }

  CoUninitialize();

  return 0;
}
