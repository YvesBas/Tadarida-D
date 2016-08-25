#ifndef DETEC_H
#define DETEC_H

#include <float.h>
#include <stdint.h>
#if QT_VERSION >= 0x050000
#define VERQT 4
#else
#define VERQT 5
#endif
#ifdef __linux
 #define LINWIN 0
#else
#define LINWIN 1
#endif
typedef int64_t  __int64;
#include <iostream>
#include <QtCore/qmath.h>
#include <QDate>
#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QSettings>
#include <QTextStream>
#include <QThread>
#include "detectreatment.h"
#define NODETERMINED  0
#define DIRECTORYMODE 1
#define FILESMODE     2

#include "deteclaunch.h"


class Detec : public QThread
{
    Q_OBJECT
public:
    Detec(DetecLaunch *pdl,QString,int,QString,int,QString,QStringList,QStringList,int,bool,int,bool,bool,bool,int);
    Detec(int);
    ~Detec();
    void                                  run();
    //
    bool                                  ErrorFileOpen;
    QTextStream                   ErrorStream;
    QTextStream                   LogStream;
    bool                                  MustCompress;
    float                                  NumtE;
    int                                     NumVer;
    DetecLaunch                  *PMainWindow;
    bool                                  ReprocessingMode;
    int                                     TE;
    bool                                  TimeFileOpen;
    QTextStream                   TimeStream;
    bool                                  IDebug;
    int                                      IThread;

private:
    bool                                  createTxtFile(QString);
    void                                  endDetec();
    bool                                  initializeDetec();
    void                                  treatOneFile(QString,QString);
    void                                  wavCut();

    int                                     _detectionThreshold;
    DetecTreatment             *_detecTreatment;
    QFile                               _errorFile;
    bool                                 _fileProblem;
    QStringList                      _firstList;
    float                                  _freqCallMin;
    int                                     _freqMin;
    int                                     _highThreshold;
    int                                     _highThreshold2;
    int                                     _jumpThreshold;
    QFile                                _logFile;
    int                                     _lowThreshold;
    int                                     _lowThreshold2;
    int                                     _modeDirFile;
    int                                     _modeFreq;
    int	                                 _nbo;
    int                                     _numThread;
    int                                     _nProcess;
    int                                     _paramVersion;
    QString                            _processSuffixe;
    QString                            _threadSuffixe;
    QString                            _resultSuffix;
    int                                     _stopThreshold;
    int	                                 _timeExpansion;
    QFile                               _timeFile;
    QString                            _txtPath;
    bool                                  _useValflag;
    QStringList                      _wavFileList;
    QString                            _wavPath;
    QStringList                     _wavRepList;
    int                                    _qN;
    int                                    _qR;
    int                                    *_tabr1;
    int                                   _widthBigControl;
    int                                   _widthLittleControl;
    bool                               _withSox;
    bool                               _withTimeCsv;
};

#endif
