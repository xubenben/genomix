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
set -e
set -o pipefail
#set -x

# args are cluster type and number of partitions per machine
if [ $# != 2 ]; then
    echo "Please specify exactly 2 arguments, (HYRACKS | PREGELIX) THREADS_PER_MACHINE " 1>&2
    exit 1;
fi
if [ "$1" == "HYRACKS" ]; then
    CLUSTER_TYPE="hyracks"
elif [ "$1" == "PREGELIX" ]; then
    CLUSTER_TYPE="pregelix"
else
    echo "unknown cluster type $1" 1>&2
    exit 1
fi
THREADS_PER_MACHINE="$2"

GENOMIX_HOME="$( dirname "$( cd "$(dirname "$0")" ; pwd -P )" )"  # script's parent dir's parent
cd "$GENOMIX_HOME"


# get generic cluster properties to templatize
. conf/cluster.properties

# make lengths of $IO_DIRS and $store equal to THREADS_PER_MACHINE but
# distribute to all the available $DISKS by round-robinning them
IO_DIRS=""
store=""
delim=""
IFS=',' read -ra DEV_ARRAY <<< "$DISKS"  # separate $DISKS into an array
for i in `seq 1 $THREADS_PER_MACHINE`; do
    dev_index=$(( ($i - 1) % ${#DEV_ARRAY[@]} ))
    device=${DEV_ARRAY[$dev_index]}
    IO_DIRS+="$delim""$device/$CLUSTER_TYPE/io_dir-"$i
    store+="$delim""$device/$CLUSTER_TYPE/store-"$i
    delim=","
done

# write the generated conf file to (hyracks|pregelix)/conf/cluster.properties
CONF_HOME="$GENOMIX_HOME/$CLUSTER_TYPE/conf"
mkdir -p $CONF_HOME

cat > "$CONF_HOME/cluster.properties" <<EOF
CCTMP_DIR="$WORKPATH/$CLUSTER_TYPE/cc"
CCLOGS_DIR="$WORKPATH/$CLUSTER_TYPE/cc/logs"
NCTMP_DIR="$WORKPATH/$CLUSTER_TYPE/nc"
NCLOGS_DIR="$WORKPATH/$CLUSTER_TYPE/nc/logs"
IO_DIRS="$IO_DIRS"
EOF

if [ "$CLUSTER_TYPE" == "hyracks" ]; then
    cat >> "$CONF_HOME/cluster.properties" <<EOF
CC_CLIENTPORT=$HYRACKS_CC_CLIENTPORT
CC_CLUSTERPORT=$HYRACKS_CC_CLUSTERPORT
CC_HTTPPORT=$HYRACKS_CC_HTTPPORT
CC_DEBUG_PORT=$HYRACKS_CC_DEBUG_PORT
NC_DEBUG_PORT=$HYRACKS_NC_DEBUG_PORT
EOF
elif [ "$CLUSTER_TYPE" == "pregelix" ]; then
    cat >> "$CONF_HOME/cluster.properties" <<EOF
CC_CLIENTPORT=$PREGELIX_CC_CLIENTPORT
CC_CLUSTERPORT=$PREGELIX_CC_CLUSTERPORT
CC_HTTPPORT=$PREGELIX_CC_HTTPPORT
CC_DEBUG_PORT=$PREGELIX_CC_DEBUG_PORT
NC_DEBUG_PORT=$PREGELIX_NC_DEBUG_PORT
EOF
    cat > "$CONF_HOME/stores.properties" <<EOF
store="$store"
EOF
else
    echo "unrecognized cluster-type $CLUSTER_TYPE" 1>&2 
    exit 1
fi
