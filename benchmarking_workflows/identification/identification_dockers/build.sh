#!/bin/bash

type -a docker > /dev/null

if [ $? -ne 0 ] ; then
	echo "UNCONFIGURED: No docker executable" 1>&2
	exit 1
fi

if [ $# -gt 0 ]; then
	tag_id="$1"

	for docker_name in i_validation i_metrics i_consolidation ; do
		docker build -t cjh4zavolab/"$docker_name":"$tag_id" "$docker_name"
	done
else
	echo "Usage: $0 tag_id" 1>&2
	exit 1
fi
