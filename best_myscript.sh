#!/bin/sh
#SBATCH --account=physics
#SBATCH --partition=ada
#SBATCH --time=7-01:00:00
#SBATCH --nodes=1 --ntasks=40
#SBATCH --job-name="P_owh_EGGING"
#SBATCH --mail-user=BRWJOS004@myuct.ac.za
#SBATCH --mail-type=ALL

set -e

#PDF_IDs=({331100..331199} {14000..14058} {27400..27464})

PDF_IDs=({14000..14058})
max_jobs=40
mainID=14000

PATH=$PATH:/home/brwjos004/physics/lhapdf/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/brwjos004/physics/pythia8312/lib


#for boson in wp wm z; do
for boson in wm z; do

	#enter file of corresponding boson
	if [[ "$boson" == "wp" || "$boson" == "wm" ]];then
	cd /scratch/brwjos004/POWHEG-BOX-V2/W || { echo "failed to enter W"; exit 1; }
	elif [[ "$boson" == "z" ]] ;then 
	cd /scratch/brwjos004/POWHEG-BOX-V2/Z || { echo "failed to enter Z"; exit 1; }
	fi


	#copy files to generate events in seperate folders, preventing interference between different pdf events
	for PDF_ID in "${PDF_IDs[@]}"; do
		cp -r ${boson}_cooking ${boson}_cooking_${PDF_ID}
	done


	#prevent stupid pc from stopping event generation in one file to go work on another before completion
	for PDF_ID in "${PDF_IDs[@]}"; do

		#wait until the number of running background jobs is less than max_jobs.
		echo $(pgrep -P $$ | wc -l)
		while (( $(pgrep -P $$ | wc -l) >= max_jobs )); do
		   sleep 1  # Pause execution until a job finishes
		done
		echo "number of jobs is $(pgrep -P $$ | wc -l)"
		#parallemisation
		(

		#enter pdf file and edit shit in file
		cd ${boson}_cooking_${PDF_ID}
		sed -i "s/^lhans1 .*/lhans1 ${PDF_ID}/" powheg.input
		sed -i "s/^lhans2 .*/lhans2 ${PDF_ID}/" powheg.input

		#execute thingy in outer file and actually generate events
		echo "Running POWHEG with PDF ID: ${PDF_ID}"
		../pwhg_main

		#moving and renaming file to appropriate pdf/boson folder
		mv pwgevents.lhe "/scratch/brwjos004/POWHEG_output/${PDF_IDs[0]}_output/${boson}_output/${boson}_pwgevents_${PDF_ID}.lhe"

		#cleaning up garbage
		cd ..
		rm -r ${boson}_cooking_${PDF_ID}


		echo "SUCKCESSfully generated powheg files ${boson}_pwgevents_${PDF_ID} now showering"

#------------------------------------------------------------------------------------------------#

	#moving to pythia showering location
	#	cd /home/brwjos004/physics/pythia8312/assignment7
		cd /scratch/brwjos004/assignment7

	#making copies of .cc and .cmnd file for each pdf and boson for easier paramellisation
		cp /home/brwjos004/physics/pythia8312/assignment7/assignment7.cc /scratch/brwjos004/assignment7/assignment7_${boson}_pwgevents_${PDF_ID}.cc
		cp /home/brwjos004/physics/pythia8312/assignment7/assignment7.cmnd /scratch/brwjos004/assignment7/assignment7_${boson}_pwgevents_${PDF_ID}.cmnd

		#renaming the files in the lhe file list
		filename="${boson}_pwgevents_${PDF_ID}"
		echo "showering $filename ..."

		#changing the input files to look at the correct PDFs
		sed -i "/TFile\* outFile = new TFile(\".*\", \"RECREATE\");/s/\"[^\"]*\"/\"${filename}.root\"/" assignment7_${filename}.cc
		sed -i "s#pythia.readFile(\"assignment7.cmnd\")#pythia.readFile(\"assignment7_${filename}.cmnd\")#" assignment7_${filename}.cc
		sed -i "/^Beams:LHEF .*/ s#=.*#=\/scratch\/brwjos004\/POWHEG_output\/${mainID}_output\/${boson}_output\/${filename}.lhe#" assignment7_${filename}.cmnd

		if ! grep -q "assignment7_${filename}.cmnd" assignment7_${filename}.cc ||
		   ! grep -q "${filename}.root" assignment7_${filename}.cc ||
		   ! grep -q "${filename}.lhe" assignment7_${filename}.cmnd; then
			echo "The file name wasn't changed bro... \*dies\*." >&2
			exit 1;
		fi


		#actual showering
		make assignment7_${filename}
		./assignment7_${filename} > assignment7_${filename}.log

		#moving files
		if [[ -f ${filename}.root && -f assignment7_${filename}.log ]]; then
			mkdir -p /scratch/brwjos004/assignment7/${mainID}_pythia_output
			mkdir -p /scratch/brwjos004/assignment7/${mainID}_pythia_output_logs

			mv ${filename}.root /scratch/brwjos004/assignment7/${mainID}_pythia_output
			mv assignment7_${filename}.log /scratch/brwjos004/assignment7/${mainID}_pythia_output_logs

			#removing now redundant pythia files
			rm assignment7_${filename}.cc
			rm assignment7_${filename}.cmnd
			rm assignment7_${filename}

			#removing redundant powheg files
			if [[ -f /scratch/brwjos004/POWHEG_output/${mainID}_output/${boson}_output/${filename}.lhe ]]; then
				rm /scratch/brwjos004/POWHEG_output/${mainID}_output/${boson}_output/${filename}.lhe
			else 
				echo "can't find corresponding ${filename}.lhe file to brutally murder. EPIC FAIL"
				echo "here's the path tho"
				echo "/scratch/brwjos004/POWHEG_output/${mainID}_output/${boson}_output/${filename}.lhe"
				exit 1;
			fi


		else
			echo "output files missing for ${filename}. \*fucking dies\*"
			exit 1;
		fi
		) &

	done

done

wait
