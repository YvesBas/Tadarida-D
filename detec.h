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
#include <QtGui/QMainWindow>
#include <QTextStream>
#include <QThread>
#include "detectreatment.h"

#define DIRECTORYMODE 1
#define FILESMODE            2



class Detec : public QThread
{
    Q_OBJECT
public:
    Detec(QString);
    Detec(int,char **);
    ~Detec();
    void            run();
    DetecTreatment  *_detecTreatment;
    bool                  _errorFileOpen;
    QTextStream   _errorStream;
    bool                  InitializeDetec();
    int                     _logVersion;
    bool                 MustCompress;
    int                     _userVersion;
    QTextStream  _logText;
    QString           ResultSuffix;

private:
    bool                createTxtFile(QString);
    void                endDetec();
    void                treatOneFile(QString,QString);

    int                   _detectionThreshold;
    QFile             _errorFile;
    bool               _fileProblem;
    QStringList   _firstList;
    float               _freqCallMin;
    int                  _freqMin;
    bool               _initPassed;
    QFile             _logFile;
    int                   _modeDirFile;
    int	              _nbo;
    int                  _stopThreshold;
    int	              _timeExpansion;
    QString         _txtPath;
    bool               _useValflag;
    QStringList   _wavFileList;
    QString         _wavPath;
    QStringList   _wavRepList;
    // for new parameters
    int                  _highThreshold;
    int                  _jumpThreshold;
    int                 _lowThreshold;
    int                 _qN;
    int                 _qR;
    int                 *_tabr1;
    int                 _widthBigControl;
    int                _widthLittleControl;
};

#endif
