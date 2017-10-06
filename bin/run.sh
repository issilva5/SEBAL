#!/bin/bash
# ${SAPS_MOUNT_POINT}/${TASK_ID}/$INPUT_DIR
INPUT_DIR_PATH=$1
# ${SAPS_MOUNT_POINT}/$OUTPUT_DIR
OUTPUT_DIR_PATH=$2
# ${SAPS_MOUNT_POINT}/$PREPROCESS_DIR
PREPROCESS_DIR_PATH=$3

# Global variables
SANDBOX=$(pwd)
SAPS_DIR_PATH=$SANDBOX/SEBAL
CONF_FILE=sebal.conf
LIBRARY_PATH=/usr/local/lib
BOUNDING_BOX_PATH=example/boundingbox_vertices
TMP_DIR_PATH=/tmp

R_EXEC_DIR=$SAPS_DIR_PATH/workspace/R
R_ALGORITHM_VERSION=Algoritmo-completo-v26062017.R
R_RASTER_TMP_DIR=/mnt/rasterTmp
MAX_TRIES=2

SCRIPTS_DIR=scripts
LOG4J_PATH=$SAPS_DIR_PATH/log4j.properties

# Required info
IMAGE_NAME=
IMAGE_MTL_PATH=
IMAGE_STATION_FILE_PATH=

# This function get required file names from input path
function getFileNames {
  files=( $INPUT_DIR_PATH/*_MTL.txt )
  for file in "${files[@]}"
  do
    filename="${file##*/}"
    filenameWithoutExtension="${filename%_MTL*}"
    IMAGE_NAME="$filenameWithoutExtension"
  done

  IMAGE_MTL_PATH=$INPUT_DIR_PATH/$IMAGE_NAME"_MTL.txt"
  IMAGE_STATION_FILE_PATH=$INPUT_DIR_PATH/$IMAGE_NAME"_station.csv"
  
  echo "ImageName: $IMAGE_NAME;  ImageMTLPath: $IMAGE_MTL_PATH; ImageStationPath: $IMAGE_STATION_FILE_PATH"
}

# This function verifies if all R dependencies are installed
function verifyRScript {
  echo "Verifying dependencies for R script"

  bash -x ${SAPS_DIR_PATH}/$SCRIPTS_DIR/verify-dependencies.sh ${SAPS_DIR_PATH}/$SCRIPTS_DIR
}


# This function clean environment by deleting raster temp dir if exists
function cleanRasterEnv {
  if [ -d $R_RASTER_TMP_DIR ]
  then
    sudo rm -r $R_RASTER_TMP_DIR
  fi

}

# This function prepare a dados.csv file
function creatingDadosCSV {
  echo "Creating dados.csv for image $TASK_ID"

  cd $R_EXEC_DIR

  echo "File images;MTL;File Station Weather;Path Output" > dados.csv
  echo "$INPUT_DIR_PATH;$IMAGE_MTL_PATH;$IMAGE_STATION_FILE_PATH;$OUTPUT_DIR_PATH" >> dados.csv
}

# This function creates a raster tmp dir if not exists and start scripts to collect CPU and memory usage
function prepareEnvAndCollectUsage {
  # check if raster temporary dir exists
  if [ ! -d $R_RASTER_TMP_DIR ]
  then
    sudo mkdir $R_RASTER_TMP_DIR
  else
    count=`ls -1 $R_RASTER_TMP_DIR/r_tmp* 2>/dev/null | wc -l`
    if [ $count != 0 ]
    then
      sudo rm -r $R_RASTER_TMP_DIR/r_tmp*
    fi
  fi

  echo "Starting CPU and Memory collect..."
  sudo bash $SAPS_DIR_PATH/$SCRIPTS_DIR/collect-cpu-usage.sh | sudo tee $OUTPUT_DIR_PATH/$IMAGE_NAME"_cpu_usage.txt" > /dev/null &
  sudo bash $SAPS_DIR_PATH/$SCRIPTS_DIR/collect-memory-usage.sh | sudo tee $OUTPUT_DIR_PATH/$IMAGE_NAME"_mem_usage.txt" > /dev/null &
}

# This function executes R script
function executeRScript {
  for i in `seq $MAX_TRIES`
  do
    cleanRasterEnv
    sudo bash $SAPS_DIR_PATH/$SCRIPTS_DIR/executeRScript.sh $R_EXEC_DIR/$R_ALGORITHM_VERSION $R_EXEC_DIR $TMP_DIR_PATH
    PROCESS_OUTPUT=$?

    echo "executeRScript_process_output=$PROCESS_OUTPUT"
    if [ $PROCESS_OUTPUT -eq 0 ]
    then
      echo "NUMBER OF TRIES $i"
      break
    elif [ $PROCESS_OUTPUT -eq 124 ] && [ $i -ge $MAX_TRIES ]
    then
      exit 124
    else
      if [ $i -ge $MAX_TRIES ]
      then
	echo "NUMBER OF TRIES $i"
        exit 1
      fi
    fi
  done
}

# This function moves dados.csv to image results dir
function mvDadosCSV {
  sudo mv dados.csv $OUTPUT_DIR_PATH
  cd ../..
}

function killCollectScripts {
  echo "Killing collect CPU and Memory scripts"
  ps -ef | grep collect-cpu-usage.sh | grep -v grep | awk '{print $2}' | xargs sudo kill
  ps -ef | grep collect-memory-usage.sh | grep -v grep | awk '{print $2}' | xargs sudo kill
}

function checkProcessOutput {
  PROCESS_OUTPUT=$?

  if [ $PROCESS_OUTPUT -ne 0 ]
  then
    finally
  fi
}

# This function ends the script
function finally {
  exit $PROCESS_OUTPUT
}

getFileNames
checkProcessOutput
verifyRScript
checkProcessOutput
creatingDadosCSV
checkProcessOutput
prepareEnvAndCollectUsage
checkProcessOutput
executeRScript
checkProcessOutput
mvDadosCSV
killCollectScripts
cleanRasterEnv
checkProcessOutput
finally
