#!/bin/sh
folder="D:/generatedMolecules/xxxxx/deep/ledock/Zinc6000"
output_file="D:/generatedMolecules/xxxxxx/deep/ledock/ligands1"
softfiles=$(ls $folder)
for sfile in ${softfiles}
do 
	echo "$folder/${sfile}" >> $output_file
done
