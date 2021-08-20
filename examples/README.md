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
If no arguments are given, the script will convert _all_ notebooks in this directory (except the demo notebook, since that will take a long time; you must explicitly pass the demo notebook as an argument to run it). Otherwise it will only convert those passed as arguments.

If you want a cell's input to be hidden (e.g., like the `env_check()` function definition), add a tag (in the gears tab) called `'hide_input'` and the input will be hidden in the final HTML.