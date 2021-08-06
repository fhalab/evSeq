#!/bin/bash

# If there are no arguments given
if [ "$#" -eq "0" ]
then
    # Iterate through all notebooks
    notebooks="*.ipynb"
else
    # Iterate only through arguments
    notebooks="$@"
fi

path=../docs/
for nb in $notebooks
do
    # Pretty printing:
    # Get length of upper and lower ### set
    N=$((${#nb}+25))
    head -c $N < /dev/zero | tr '\0' '#'
    echo
    echo "###### Processing $nb ######"
    head -c $N < /dev/zero | tr '\0' '#'
    echo

    # Get the basename
    base="${nb%.*}"
    
    # Specify output
    new_path="$path$base.html"

    # Rerun notebook
    jupyter nbconvert --to notebook --execute --clear-output --inplace $nb

    # Touch new html doc with jekyll front matter
    echo -e "---\nlayout: default\n---\n" > $new_path

    # Convert to html and add to html file
    jupyter nbconvert $nb --to html --stdout >> $new_path

    # Update the relative path links to assets directory
    sed -i '' -e 's&../docs/assets/&assets/&g' $new_path
done
