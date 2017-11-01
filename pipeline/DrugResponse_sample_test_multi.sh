#!/bin/bash
# DrugResponse_sample_test_multi.sh: Example of piping DrugResponse package to test cancer samples
# Input: folder with multiple samples
# Usage:
# DrugResponse_sample_test_URL.sh path.to.folder

E_FILE_NOT_EXIST=70
E_DOWNLOAD_FAIL=71
E_MV_FAIL=72
E_UNZIP_FAIL=73
E_TEST_FAIL=74

INPUT=$1
PROJECTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/.." && pwd )"
WORKDIR="$( cd "$INPUT" && pwd )"
DATAPATH="$PROJECTDIR/KEADrugResponse/data"
TDATAPATH=$(cd $(dirname "$DATAPATH") && pwd -P)/$(basename "$DATAPATH")

if [ ! -e "$WORKDIR" ]; then
	echo "File \""$WORKDIR"\" does not exist."
	exit $E_FILE_NOT_EXIST
fi

echo "Data folder $INPUT checked"

cd $WORKDIR

for i in *.[Cc][Ee][Ll]
do
	(
		cd "$PROJECTDIR"
		if R -q -e "library(KEADrugResponse);library(affy);DrugResponse.predict('$WORKDIR/$i','cel','$i','$TDATAPATH','$WORKDIR')"; then
			echo "test for $i successful"
		else
			echo "Fail to process $i"
		fi
	)
done

echo "done"
