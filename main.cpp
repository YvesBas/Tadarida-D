
#include "deteclaunch.h"


class DetecLaunch;

int main(int argc, char *argv[])
{
    if(argc<1) return(-1);
    DetecLaunch *dl= new DetecLaunch();
    dl->treat(argc,argv);
    delete dl;
    exit(0);
}

