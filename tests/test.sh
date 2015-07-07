#! /bin/sh

WORK_DIR=`mktemp -d`
BASE_DIR=`readlink -f \`dirname $0\``
TADARIDAD=$BASE_DIR/../tadaridaD
OPTS=$@
WAVES_DIR=$BASE_DIR/waves

echo "Bootstrap env in $WORK_DIR..."
cp -R $WAVES_DIR $WORK_DIR
cd $WORK_DIR
echo "Running tadaridaD with '$OPTS'"
$TADARIDAD $OPTS waves
exit $?
