#include "detec.h"

#ifdef _WIN32
# define SLEEP _sleep
#else
# include <unistd.h>
# define SLEEP sleep
#endif

int main(int argc, char *argv[])
{
    if(argc>1)
    {
        Detec *detec = new Detec(argc, argv);
        detec->start();
        while(detec->isRunning()) SLEEP(100);
        delete detec;
    }
    return(0);
}
