# plasma-source

## Repository for plasma source development.

### Setup Anaconda and the virtual enviornment

Instructions updated for conda 4.6.14. To check the version use the command
```
conda --version
```

First step is to download the repository. Depending on what operating system you are you using you can do this by
either downloading GitHub desktop (Mac, Windows) or by installing git and then downloading the repository. 

Once you have the repository you should install the Anaconda Distribution which will install python and python tools.

Anaconda will handle setting up your virtual enviornment and managing packages. Open a terminal (on Windows open Anaconda Propmt)
go into therepository directory, it should contain a requirements.txt file. Inside of the directory use the following command
```
conda create --name CU-PWFA
```
This will create a virtual enviornment called "CU-PWFA". To activate the enviornment:
```
conda activate CU-PWFA
```
To install all required packages use the command
```
```


### Update the virtual enviornment

The point of the enviornment is that it contains all the packages and libraries needed to run the code in the repository. When
installing the plasma source code on a new computer the required packages need to be downloaded by updating the enviornment.
Whenever you update your local repository using the pull command, check to see if the requirement.txt file changed. If it changed
you need to update your enviornment using
```
conda env update -f requirements.txt
```

When you install a new package into your enviornment you need to add it to the master package list. Then everyone
else can use the master list to update their enviornment by installing the package.

To save your current enviornment to file, first activate the enviornment:
```
conda activate CU-PWFA
```
Next, use the command:
```
conda env export > requirements.txt
```


### Setup Spyder and IPython

To use spyder inside of the virtual enviornment, activate the virtual enviornment and launch spyder from the terminal.

To use IPython in a Jupyter-Notebook you need to install the IPython kernel in the enviornment. The package should have been downloaded 
when the virtual enviornment was created. To install use the command:
```
python -m ipykernel install --user --name=CU-PWFA
```
Now when you create a Jupyter-Notebook it should give you the option to use the CU-PWFA kernel.


### Compile the beam package

The beam package is a collection of simulation modules for lasers and particle beams. It is written in cython and needs to be compiled.
On linux this is straightforward, use the terminal to navigate to the plasma-source/python folder and run the command
```
python beam/calc/setup.py build_ext --inplace
```
If everything works each file should be successfully compiled. The package needs to be recompiled anytime the cython files change.

On windows it is a bit more complicated, you need a compiler so download the Visual Studio build tools and then use the same command
as above. In windows the code is compiled to c++ because the visual studio compiler doesn't follwo the complex specification in C.


### VSim enviornment

Currently no one uses this. We generally analyze VSim results by directly loading data from the hdf5 files.
It is no longer maintained or updated.

VSim's analysis scripts are written in python 2.7. Unfortunatly, everything we do is written in python 3.4. This means 
VSim scripts won't run in our virtual enviornment so we have a second virtual enviornment for VSim where you can write analyzer 
scripts or scripts to extract data from VSim results files. Unless you are writing an analyzer, I would suggest extracting
data using VSim's internal scripts and saving it to a numpy file in Python 2.7 and doing all analysis in Python 3.4.

There is a seperate VSim file in the repository that contains the requirements file for the VSim virtual enviornment. Place all
scripts meant to run in this enviornment in this file.

To set up the VSim virtual enviornment, first cd into the VSim folder. Next, create the enviornment using 
```
conda create --name VSim --file requirements.txt
```
To add the enviornment to Ipython and Jupyter-Notebook, follow the same steps as above. To update the enviornment, follow the
steps above.

