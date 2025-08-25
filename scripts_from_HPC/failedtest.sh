#!/bin/bash


#determine number of events
NoofEventsThingy=$(grep -m 1 "numevts" wp_pwgevents_14058_tests.lhe)
NoOfEventsThingyArray=(${NoofEventsThingy})
NoOfEvents=${NoOfEventsThingyArray[1]}

weightString=$(awk  '/<wgt id='\''1'\''>/ { print $3 }' wp_pwgevents_14058_tests.lhe)	#storing all weights in a single variable
weightArray=(${weightString}) #converting that variable into an array of weights

eventLineNoList=$(awk '/<event>/ {print FNR}' wp_pwgevents_14058_tests.lhe) #stores the relevant line numbers in a single variable
eventLineNoArray=(${eventLineNoList}) #converting that variable into an array of weight

echo $NoOfEvents

#looping over events to replace original weight with new weight
awk -v eventLines="$eventLineNoList" -v Weights="$weightString" '
	BEGIN {
		split(eventLines,eventLineArray,\n)
		split(Weights,WeightArray,\n)
	}

	{
		for (i = 1; i <= (length(eventLineArray)+1); i++) 
		{
			#check if the line number (NR) is the line with the XWGTUP and replace it with the new weight.
			if (NR == (eventLineArray[i]+1)) { $3 = WeightArray[i]}
		}
			print $0
	}
' wp_pwgevents_14058_tests.lhe > testOutput.out

	#obtaining relevant weight and line number
	#weight=${weightArray[${eventNo} - 1]}
	#eventLineNo=$((${eventLineNoArray[${eventNo} -1]} + 1))

#	sed -Ei "${eventLineNo}s/[^ ]+/${weight}/3" wp_pwgevents_14058_tests.lhe #replace the number at line "eventLineNo", position 3 with "weight"
#	awk -v line=$eventLineNo -v word=3 -v value="$weight" 'NR==line{$word=value}1' wp_pwgevents_14058_tests.lhe


	#echo $weight


