#!/bin/bash
#------------------------------------------------------------------------
# Copyright 2009-2013 by The Regents of the University of California
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# you may obtain a copy of the License from
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ------------------------------------------------------------------------
#set -x

GENOMIX_HOME="$( dirname "$( cd "$(dirname "$0")" ; pwd -P )" )"  # script's parent dir's parent
cd "$GENOMIX_HOME"

if [[ $# != 1 || ("$1" != "HYRACKS" && "$1" != "PREGELIX") ]]; then
    echo "please provide a cluster type as HYRACKS or PREGELIX! (saw $1)" 1>&2
    exit 1
fi
if [ $1 == "HYRACKS" ]; then
    MY_HOME="$GENOMIX_HOME/hyracks"
else 
    MY_HOME="$GENOMIX_HOME/pregelix"
fi

PID_FILE="$MY_HOME/cc.pid"

if [ ! -e $PID_FILE ] ; then
  echo "No pid file for CC! (expected $PID_FILE)" 
  exit 1
else
  PID=`cat $PID_FILE`
  kill -9 $PID
  rm $PID_FILE
  echo "Stopped CC on master: "`hostname`$'\t'$PID
fi
