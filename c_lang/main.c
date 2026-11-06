#include <locale.h>

int main()
{
    #ifdef _WIN32
        setlocale(LC_ALL, ".UTF8");
    #else
        setlocale(LC_ALL, "");
    #endif

    return 0;
}
