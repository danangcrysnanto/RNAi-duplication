
#I use this short script to automate BEAST analysis for each xml files 

for file in *xml
do
beast -beagle_GPU -beagle_double -overwrite $file
done
