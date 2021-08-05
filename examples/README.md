# `evSeq` Examples

This directory holds two things:

1. Jupyter Notebooks that are used for the docs, but can also be run by themselves.
2. Examples files for testing `evSeq` and running the notebooks.

### For docs:
When you make a change to an example notebook that is used for the docs and want to publish that change, you first need to convert the file to replace the old .html file in the docs folder. This is most easily done with the script `convert_to_docs.sh`.

First, you will need to make the script executable with
```
chmod +x convert_to_docs.sh
```
Then execute the script with
```
./convert_to_docs.sh [FILE1] [FILE2] ...
```
If no arguments are given, the script will convert _all_ notebooks in this directory. Otherwise it will only convert those passed as arguments.