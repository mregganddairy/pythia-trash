#!/bin/sh
#SBATCH --account=physics
#SBATCH --partition=ada
#SBATCH --time=7-01:00:00
#SBATCH --nodes=1 --ntasks=3
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


#prevent stupid pc from stopping event generation in one file to go work on another before completion
#wait until the number of running background jobs is less than max_jobs.
#echo $(pgrep -P $$ | wc -l)
#while (( $(pgrep -P $$ | wc -l) >= max_jobs )); do
#   sleep 1  # Pause execution until a job finishes
#done
#echo "number of jobs is $(pgrep -P $$ | wc -l)"
##parallemisation
#(

for boson in wp wm z; do

	#enter file of corresponding boson
	if [[ "$boson" == "wp" || "$boson" == "wm" ]];then
	cd /scratch/brwjos004/POWHEG-BOX-V2/W || { echo "failed to enter W"; exit 1; }
	elif [[ "$boson" == "z" ]] ;then 
	cd /scratch/brwjos004/POWHEG-BOX-V2/Z || { echo "failed to enter Z"; exit 1; }
	fi


	#copy files to generate events in seperate folders, preventing interference between different pdf events
	cp -r ${boson}_cooking ${boson}_cooking_${mainID}
	cd ${boson}_cooking_${PDF_ID}


	ID=0
	for PDF_ID in "${PDF_IDs[@]}"; do


		#enter pdf file and edit shit in file
		sed -i "s/^lhans1 .*/lhans1 ${PDF_ID}/" powheg.input
		sed -i "s/^lhans2 .*/lhans2 ${PDF_ID}/" powheg.input

		#edditing associated with reweighting
		sed -i "s/^lhrwgt_id .*/lhrwgt_id '${ID}'/" powheg.input
		sed -i "s/^lhrwgt_descr .*/lhrwgt_descr '${PDF_ID}'/" powheg.input

		if (( ID == 0 )); then
			sed -i "82s/.*/# lhrwgt_group_name '${mainID}' /" powheg.input
			sed -i "83s/.*/# lhrwgt_group_combine 'hessian'/" powheg.input
			sed -i "85s/.*/# compute_rwgt 1 ! calculates the new lhe file using weights instead of raw dogging it./" powheg.input
		elif (( ID > 0 )); then 
			sed -i "82s/.*/lhrwgt_group_name '${mainID}' /" powheg.input
			sed -i "83s/.*/lhrwgt_group_combine 'hessian'/" powheg.input
			sed -i "85s/.*/compute_rwgt 1 ! calculates the new lhe file using weights instead of raw dogging it./" powheg.input
		fi

		#execute thingy in outer file and actually generate events
		if (( ID == 0 )); then
		echo "Running MAIN POWHEG with PDF ID: ${PDF_ID}"
		../pwhg_main
		elif (( ID > 0 )); then 
		echo "Reweighting POWHEG with PDF ID: ${PDF_ID}"
		../pwhg_main 
		mv -f pwgevents-rwgt.lhe pwgevents.lhe
		fi
		ID=$((ID + 1))

	done

		#moving and renaming file to appropriate pdf/boson folder
		mv pwgevents.lhe "/scratch/brwjos004/POWHEG_output/${PDF_IDs[0]}_output/${boson}_output/${boson}_pwgevents_${mainID}.lhe"

		#cleaning up garbage
		cd ..
		rm -r ${boson}_cooking_${mainID}



done

echo "SUCKCESSfully generated powheg files ${boson}_pwgevents_${PDF_ID}. Proceed to the showering phase."

#) &

wait
