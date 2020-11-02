# Conformal Prediction evalulation of MS study

## Requirements
* [Nonconformist](https://github.com/donlnz/nonconformist) using commit 91fca869b7421c0658fd93590a6d84d23a96072d which implemented mondrian TCP and supports later versions of Scikit learn. The easiest way to use this is by cloning this github repo and adding the path to your PYTHONPATH.
* Packages specified in the `env.yaml` file in this repo (can be installed by running `conda env create --file env.yaml` if you have anaconda installed)
* **Optional** [plot_utils](https://github.com/pharmbio/plot_utils) plotting package used for generating plots (not identical with the ones in the manuscript, they are generated separatetly from R-code). 

## Notes
The code is organized with a `common_utils` file which loads data, performes scaling and other utility methods. The scripts that should be runned for doing the analysis is files called `stepX`. Note that the optimal SVM parameters found in Step2 has been hardcoded in the `common_utils`-file. 
