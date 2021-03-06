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

#get the OS
os_name=`uname -a|awk '{print $1}'`
linux_os='Linux'

if [ $os_name = $linux_os ];
then
    #Get IP Address
    #Prefer Infiniband connection
#    IPADDR=`/sbin/ifconfig ib0 2> /dev/null | grep "inet " | awk '{print $2}' | cut -f 2 -d ':'`
#    if [ "$IPADDR" = "" ]
#    then
        ipaddr=`/sbin/ifconfig eth0 | grep "inet " | awk '{print $2}' | cut -f 2 -d ':'`
        if [ "$ipaddr" = "" ]
        then
            ipaddr=`/sbin/ifconfig lo | grep "inet " | awk '{print $2}' | cut -f 2 -d ':'`
        fi 
#    fi 
else
        ipaddr=`/sbin/ifconfig en1 | grep "inet " | awk '{print $2}' | cut -f 2 -d ':'`
	if [ "$ipaddr" = "" ]
        then
                ipaddr=`/sbin/ifconfig lo0 | grep "inet " | awk '{print $2}' | cut -f 2 -d ':'`
        fi

fi
echo $ipaddr
