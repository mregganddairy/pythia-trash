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

PDF_IDs=({27400..27464})
max_jobs=40
mainID=27400

PATH=$PATH:/home/brwjos004/physics/lhapdf/bin
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/home/brwjos004/physics/pythia8312/lib



#for boson in wp wm z; do
for boson in wp wm z; do

	#enter file of corresponding boson
	if [[ "$boson" == "wp" || "$boson" == "wm" ]];then
	cd /scratch/brwjos004/POWHEG-BOX-V2/W || { echo "failed to enter W"; exit 1; }
	elif [[ "$boson" == "z" ]] ;then 
	cd /scratch/brwjos004/POWHEG-BOX-V2/Z  || { echo "failed to enter Z"; exit 1; }
	fi

#########################################################################################################
####creating main folder which will be used for the reweighting after main LHE file has been generated###
#########################################################################################################
	cp -r ${boson}_cooking ${boson}_cooking_${mainID}
	cd ${boson}_cooking_${mainID}

	#enter pdf file and edit shit in file
	sed -i "s/^lhans1 .*/lhans1 ${mainID}/" powheg.input
	sed -i "s/^lhans2 .*/lhans2 ${mainID}/" powheg.input

	#edditing associated with reweighting
	sed -i "s/^lhrwgt_id .*/lhrwgt_id '0'/" powheg.input
	sed -i "s/^lhrwgt_descr .*/lhrwgt_descr '${mainID}'/" powheg.input

	#MAKE SURE THE REWEIGHTING LINES ARE AT LINES 82, 83 AND 85 OR IT WON'T WORK
	sed -i "82s/.*/# lhrwgt_group_name '${mainID}' /" powheg.input
	sed -i "83s/.*/# lhrwgt_group_combine 'hessian'/" powheg.input
	sed -i "85s/.*/# compute_rwgt 1 ! calculates the new lhe file using weights instead of raw dogging it./" powheg.input

	#execute thingy in outer file and actually generate events
	echo "Running MAIN POWHEG with PDF ID: ${mainID}"
	../pwhg_main
	cd ..

echo " MAIN POWHEG file with PDF ID: ${mainID} probably SUCKCESSfully generated."
#------------------------------------------------------------------------------------------------#



###################################################################################################
###################################################################################################
###################################################################################################

	

###################################################################################################
#####performing reweighting procedure in parallel##################################################
###################################################################################################

	#copy files to generate events in seperate folders, preventing interference between different pdf events
	for PDF_ID in "${PDF_IDs[@]}"; do
		if [[ "${boson}_cooking_${PDF_ID}" != "${boson}_cooking_${mainID}" ]]; then
			cp -r ${boson}_cooking_${mainID} ${boson}_cooking_${PDF_ID}
		fi
	done

	#performing reweighting procedure in parallel
	for PDF_ID in "${PDF_IDs[@]}"; do

		#prevent stupid pc from stopping event generation in one file to go work on another before completion
		#wait until the number of running background jobs is less than max_jobs.
		echo $(pgrep -P $$ | wc -l)
		while (( $(pgrep -P $$ | wc -l) >= max_jobs )); do
		   sleep 1  # Pause execution until a job finishes
		done
		echo "number of jobs is $(pgrep -P $$ | wc -l)"
		#parallemisation
		(

		#calculating new weight for each PDF_ID folder storing it seperately in each file
		#<weight id ='1'> should correspond to the weight of pdf $PDF_ID
		cd ${boson}_cooking_${PDF_ID} 

		filename="${boson}_pwgevents_${PDF_ID}"
		if [[ "${boson}_cooking_${PDF_ID}" != "${boson}_cooking_${mainID}" ]]; then
			sed -i "s/^lhans1 .*/lhans1 ${PDF_ID}/" powheg.input
			sed -i "s/^lhans2 .*/lhans2 ${PDF_ID}/" powheg.input

			sed -i "s/^lhrwgt_id .*/lhrwgt_id '1'/" powheg.input
			sed -i "s/^lhrwgt_descr .*/lhrwgt_descr '${PDF_ID}'/" powheg.input

			sed -i "82s/.*/lhrwgt_group_name '${mainID}' /" powheg.input
			sed -i "83s/.*/lhrwgt_group_combine 'gaussian'/" powheg.input
			sed -i "85s/.*/compute_rwgt 1 ! calculates the new lhe file using weights instead of raw dogging it./" powheg.input

			echo "Reweighting POWHEG with PDF ID: ${PDF_ID}"
			../pwhg_main
			mv pwgevents-rwgt.lhe pwgevents.lhe

###################################################################################################
###################################################################################################
###################################################################################################

####################################################################################
######looping over events to replace original weight with new weight################
####################################################################################
		
			echo "rearranging ${boson}_pwgevents_${PDF_ID}"
			#finding number of weights and line numbers and storing them in seperate text files 
			awk '/<wgt id='\''1'\''>/ { print $3 }' pwgevents.lhe > weights.txt 
			awk '/<event>/ {print FNR}' pwgevents.lhe > event_lines.txt

			#replacing input file name in C++ file event reweighter with boson and pdf name
			echo "attempting to rename input files to ${boson}_pwgevents_${PDF_ID}.lhe in Reorganiser.cpp" 
			sed -i "s/const std::string input_file = .*/const std::string input_file = \"pwgevents.lhe\";/" Reorganiser.cpp
			sed -i "s/const std::string output_file = .*/const std::string output_file = \"${boson}_pwgevents_${PDF_ID}.lhe\";/" Reorganiser.cpp

			if ! grep -q "const std::string input_file = \"pwgevents.lhe\";" Reorganiser.cpp ||
			   ! grep -q "const std::string output_file = \"${boson}_pwgevents_${PDF_ID}.lhe\";" Reorganiser.cpp; then
				echo "Failed to change the name of of input and/or output file in the Reorganiser.cpp... \*dies\*." >&2
				exit 1;
			fi

			#making executable for the Reorganiser and executing

			filename="${boson}_pwgevents_${PDF_ID}"
			g++ -O2 -march=native -std=c++17 Reorganiser.cpp -o Reorganiser
			./Reorganiser > ${filename}.lhe


			echo "${filename}.lhe successfully rearranged"
			echo "File size is $(ls -l --block-size=M)"
		fi

		if [[ "${boson}_cooking_${PDF_ID}" == "${boson}_cooking_${mainID}" ]]; then
			mv pwgevents.lhe "/scratch/brwjos004/POWHEG_output/${PDF_IDs[0]}_output/${boson}_output/${filename}.lhe"
		else
			mv ${filename}.lhe "/scratch/brwjos004/POWHEG_output/${PDF_IDs[0]}_output/${boson}_output/${filename}.lhe"
		fi
		#moving and renaming file to appropriate pdf/boson folder

#######################################################################################
#######################################################################################
#######################################################################################

		#cleaning up garbage
		cd ..
		rm -r ${boson}_cooking_${PDF_ID}
		echo "REMOVAL OF FILES COMMENTED OUT FOR DEBUGGING PURPOSES."
		
		echo "SUCKCESSfully generated powheg file ${boson}_pwgevents_${PDF_ID}. Proceed to showering"

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------# 
		#moving to pythia showering location
	#	cd /home/brwjos004/physics/pythia8312/assignment7
		cd /scratch/brwjos004/assignment7

	#making copies of .cc and .cmnd file for each pdf and boson for easier paramellisation
		cp /home/brwjos004/physics/pythia8312/assignment7/assignment7.cc /scratch/brwjos004/assignment7/assignment7_${filename}.cc
		cp /home/brwjos004/physics/pythia8312/assignment7/assignment7.cmnd /scratch/brwjos004/assignment7/assignment7_${filename}.cmnd

		#renaming the files in the lhe file list
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
			#	echo "file removal commented out. change this"
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
		echo "successfully showered ${filename}"

		) &

	done

done

wait
