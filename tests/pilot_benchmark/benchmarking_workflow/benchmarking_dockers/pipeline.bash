#!/bin/bash

#Run the tcga visualizer pipeline

REALPATH="$(realpath "$0")"
BASEDIR="$(dirname "$REALPATH")"
case "$BASEDIR" in
	/*)
		true
		;;
	*)
		BASEDIR="${PWD}/$BASEDIR"
		;;
esac

TCGA_DIR="${BASEDIR}"/TCGA_full_data
TAG=0.3.1

if [ $# -gt 1 ] ; then 
	input="$1"
	RESDIR="$2"
	shift 2
	
	if [ $# -gt 0 ] ; then
		PARTICIPANT="$1"
		shift
	else
		PARTICIPANT=NEW_PARTICIPANT
	fi
	
	if [ $# -gt 0 ] ; then
		CANCER_TYPES="$@"
	else
		CANCER_TYPES="ACC BRCA"
	fi
	
	cat <<EOF
* Using version $TAG
* Running parameters

  Input: $input
  Results: $RESDIR
  Participant: $PARTICIPANT
  Cancer types: $CANCER_TYPES
EOF
	
	if [ ! -f "$input" ] ; then
		echo "ERROR: input file does not exist" 1>&2
		exit 1
	fi
	echo "* Deriving input directory"
	inputRealPath="$(realpath "$input")"
	inputBasename="$(basename "$input")"
	INPUTDIR="$(dirname "$inputRealPath")"
	case "$INPUTDIR" in
		/*)
			true
			;;
		*)
			INPUTDIR="${PWD}/$INPUTDIR"
			;;
	esac
	
	echo "* Creating $RESDIR (if it does not exist)"
	mkdir -p "$RESDIR"
	
	# REMEMBER: We need absolute paths for docker
	RESDIRreal="$(realpath "$RESDIR")"
	case "$RESDIRreal" in
		/*)
			true
			;;
		*)
			RESDIRreal="${PWD}/$RESDIRreal"
			;;
	esac
	
	ASSESSDIR="${TCGA_DIR}"/data
	METRICS_DIR="${TCGA_DIR}"/metrics_ref_datasets
	PUBLIC_REF_DIR="${TCGA_DIR}"/public_ref

	echo "=> Validating input" && \
	docker run --rm -u $UID -v "${INPUTDIR}":/app/input:ro -v "${PUBLIC_REF_DIR}":/app/ref:ro tcga_validation:"$TAG" \
		-i /app/input/"${inputBasename}" -r /app/ref/ && \
	echo "=> Computing metrics" && \
	docker run --rm -u $UID -v "${INPUTDIR}":/app/input:ro -v "${METRICS_DIR}":/app/metrics:ro -v "${RESDIRreal}":/app/results:rw tcga_metrics:"$TAG" \
		-i /app/input/"${inputBasename}" -c $CANCER_TYPES -m /app/metrics/ -p "${PARTICIPANT}" -o /app/results/ && \
	echo "=> Assessing metrics" && \
	docker run --rm -u $UID -v "${ASSESSDIR}":/app/assess:ro -v "${RESDIRreal}":/app/results:rw tcga_assessment:"$TAG" \
		-b /app/assess/ -p /app/results/ -o /app/results/ && \
	echo "* Pipeline has finished properly"


	#Build de imagenes:

	#docker build -t tcga_validation .
	#docker build -t tcga_metrics .
	#docker build -t manage_assessment_data .
else
	echo "Usage: $0 input_file results_dir [participant_id [cancer_type]*]"
fi
