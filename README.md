

# Connectivity Analysis of Julich and Lyon Regions

## Description
This repository contains MATLAB scripts for analyzing connectivity between the Julich and Lyon brain regions in different primates. The scripts load in connectivity matrices, perform PCA and Fisher r-to-z transformations, and finally visualize the data using plots.

## Table of Contents
1. [Installation](#Installation)
2. [Usage](#Usage)
3. [Contributing](#Contributing)
4. [License](#License)

## Installation
To install, clone the repository to your local machine:
```
git clone https://github.com/yourusername/yourrepository.git
```
You will also need MATLAB to run these scripts. The scripts were developed and tested on MATLAB R2019a. 

## Usage
1. Open MATLAB and navigate to the directory containing the scripts.
2. The main script to run is `run_functional_connectivity_frontal.m`.
3. The accompanying data is available on EBRAINS and BALSA, and neuroimaging data (from Mars et al at Oxford, on PRIME-DE)
4. The constants used in the scripts can be modified in the `connectivity_analysis.m` file.
5. After running the script, you should see several `.xlsx` files generated in the same directory, along with a series of plots displayed in the MATLAB figure window.

## Further details  & references
See the paper for further details and references:
Rapan, Lucija, Sean Froudist-Walsh, Meiqi Niu, Ting Xu, Ling Zhao, Thomas Funck, Xiao Jing Wang, Karin Amunts, and Nicola Palomero-Gallagher. "Cytoarchitectonic, receptor distribution and functional connectivity analyses of the macaque frontal lobe." eLife (2023, in press)



## License
[MIT](https://choosealicense.com/licenses/mit/)

