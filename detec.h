#ifndef DETEC_H
#define DETEC_H
#include <float.h>
#include <iostream>
#include <QDate>
#include <QDateTime>
#include <QDir>
#include <QFile>
#include <QSettings>
#include <QtCore/qmath.h>
// #include <QtGui/QMainWindow>
#include <QTextStream>
#include <QThread>
#include "detectreatment.h"

#define DIRECTORYMODE 1
#define FILESMODE            2

#include "deteclaunch.h"

class Detec : public QThread
{
    Q_OBJECT
public:
    Detec(DetecLaunch *pdl,QString,int,QString,int,QString,QStringList,QStringList,int,bool,int,bool);
    Detec(int);
    ~Detec();
    void                 run();
    bool                 InitializeDetec();
    bool                createTxtFile(QString);
    //
    DetecLaunch *PDL;
    DetecTreatment  *_detecTreatment;
    int                     _logVersion;
    bool                 MustCompress;
    int                     _userVersion;
    QString           ResultSuffix;
    QTextStream  _timeStream;
    bool                 _timeFileOpen;
    bool                 _imageData;

    bool                  _errorFileOpen;
    QTextStream   _errorStream;
    QTextStream  _logText;
    bool                  _withTimeCsv;
    bool                  ReprocessingMode;
    float _numtE;
    int _tE;
    int _numVer;
    bool _xmoitie;
    bool IDebug;
    int IThread;

private:
    int                   _numThread;
    int                   _nProcess;
    QString          _threadSuffixe;
    QString          _processSuffixe;
    void                endDetec();
    void                treatOneFile(QString,QString);

    int                   _detectionThreshold;
    QFile             _errorFile;
    bool               _fileProblem;
    QStringList   _firstList;
    float               _freqCallMin;
    int                  _freqMin;
    QFile             _logFile;
    int                   _modeDirFile;
    int	              _nbo;
    int                  _stopThreshold;
    int	              _timeExpansion;
    QFile             _timeFile;
    QString         _txtPath;
    bool               _useValflag;
    QStringList   _wavFileList;
    QString         _wavPath;
    QStringList   _wavRepList;
    // for new parameters
    int                  _highThreshold;
    int                  _jumpThreshold;
    int                 _lowThreshold;
    int                  _highThreshold2;
    int                 _lowThreshold2;
    int                 _qN;
    int                 _qR;
    int                 *_tabr1;
    int                 _widthBigControl;
    int                _widthLittleControl;
    //bool            _withNewParams;
    int                _paramVersion;
};

#endif
