#!/bin/bash

#PDFs=($(ls ./rawPWHG/))
PDFs=("wm_pwgevents_14000.lhe")

for PDF in "${PDFs[@]}"; do

	# renaming the files in the lhe file list
	filename="${PDF%.lhe}"

	echo "showering $filename ..."
	# changing the input files to look at the correct PDFs
    sed -i "/TFile\* outFile = new TFile(\".*\", \"RECREATE\");/s/\"[^\"]*\"/\"${filename}.root\"/" assignment7_parallel.cc
	sed -i "/^Beams:LHEF .*/s/=.*/=\.\/rawPWHG\/${filename}.lhe/" assignment7.cmnd

	#actual showering
	make assignment7_parallel
	./assignment7_parallel > assignment7${filename}.log

	#moving files
	if [[ -f ${filename}.root && -f assignment7${filename}.log ]]; then
        mv assignment7${filename}.log ${filename}.root ./output
    else
        echo "output files missing for ${filename}. \*fucking dies\*"
		break
    fi
done

