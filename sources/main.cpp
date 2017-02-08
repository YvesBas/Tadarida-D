
#include "deteclaunch.h"
#include <QCoreApplication>


class DetecLaunch;

int main(int argc, char *argv[])
{
    if(argc<1) return(-1);
	// Deteclaunch : main class of tadaridaD (non graphic class)
    DetecLaunch *dl= new DetecLaunch();
    // Treat : this method manages the whole treatment
	dl->Treat(argc,argv);
    delete dl;
    SLEEP(3000);
    exit(0);
}

