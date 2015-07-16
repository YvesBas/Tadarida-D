#ifndef DETECLAUNCH_H
#define DETECLAUNCH_H

#ifdef _WIN32
# define SLEEP(time) _sleep(time)
#else
# include <unistd.h>
# define SLEEP(time) usleep(time * 1000)
#endif

#define PARAMNO 0
#define PARAMEXPANSION 1
#define PARAMCOMPRESS 2
#define PARAMHELP 3
#define PARAMNTHREADS 4
#define PARAMNPROCESS 5
#define PARAMNCALLED 6
#define PARAMNVERSION 7


#include <QObject>
#include <QString>
#include <QStringList>
#include <QProcess>
#include <QFile>
#include <QTextStream>
#include "fftw3.h"

#define	FFT_HEIGHT_MAX 4096


class DetecLaunch : public QObject
{
Q_OBJECT
public:
    explicit DetecLaunch(QObject *parent=0);
    ~DetecLaunch();
    bool treat(int,char **);
    bool _withTimeCsv;
    bool IDebug;
    fftwf_complex*	             _complexInput[10];
    fftwf_complex*	             _fftRes[10];
    fftwf_plan		             Plan[10][6];


private slots:
    void processFinished(int,QProcess::ExitStatus);
    void processStarted();
    void processError(QProcess::ProcessError);


private :
    bool createTxtFile(QString);
    int _modeDirFile;
    QString _wavPath;
    QStringList   _wavFileList;
    QStringList   _wavFileListProcess;
    QStringList   _wavRepList;
    int _timeExpansion;
    int _nwf;
    int _nwfp;
    int _nbThreads;
    int _nbProcess;
    int _nbPec;
    int _nCalled;
    QTextStream  _logText;
    QFile _logFile;
    int _paramVersion;
};

#endif // DETECLAUNCH_H
