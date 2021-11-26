#! /bin/bash

baseDir=$(dirname $BASH_SOURCE)

# assuming the configuration file is in E1.yaml 
# in the same directory as run.E1.sh

currentPath=$(grep 'working_dir' $baseDir/E1.yaml | sed 's/working_dir://g')

# check to see if the path specified as working directory
# exists. If it does, proceed with the analysis. Otherwise,
# create a new configuration file where the path has been 
# substituted with the path to this shell script

echo "The specified working directory:  $currentPath"

if [ -d $currentPath ]; then
    echo "The working directory specified in E1.yaml exists."
    echo "proceeding to analysis"
    # run schism in sequential mode
    runSchism analyze --config $baseDir/E1.yaml
else
    echo "The working directory specified in E1.yaml does not exists."
    echo "Using the path to this shell script as the working directory."

    echo "# adjusted path to working directory" >> $baseDir/E1.path.yaml
    echo "working_dir: $baseDir" >> $baseDir/E1.path.yaml
    echo "" >> $baseDir/E1.path.yaml
    grep -v "working_dir" $baseDir/E1.yaml >> $baseDir/E1.path.yaml
    # run schism in sequential mode
    runSchism analyze --config $baseDir/E1.path.yaml
fi
