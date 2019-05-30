# plasma-source

## Repository for plasma source development.

### Setup Anaconda and the virtual environment

Instructions updated for conda 4.6.14. To check the version use the command
```
conda --version
```

First step is to download the repository. Depending on what operating system you are you using you can do this by
either downloading GitHub desktop (Mac, Windows) or by installing git and then downloading the repository. 

Once you have the repository you should install the Anaconda Distribution which will install python and python tools.

Anaconda will handle setting up your virtual environment and managing packages. Open a terminal (on Windows open Anaconda Prompt)
go into the repository directory, it should contain a requirements.txt file. Inside of the directory use the following command
```
conda create --name CU-PWFA
```
This will create a virtual environment called "CU-PWFA". To activate the environment:
```
conda activate CU-PWFA
```
To install all required packages, navigate into the top directory of the repository and use the command for your operating system
```
conda env update -f requirements_linux.yml
```
```
conda env update -f requirements_windows.yml
```

If the environment won't build, install each package manually following the commands in the requirements.txt file.


### Update the virtual environment

The point of the environment is that it contains all the packages and libraries needed to run the code in the repository. When
installing the plasma source code on a new computer the required packages need to be downloaded by updating the environment.
Whenever you update your local repository using the pull command, check to see if the requirement.yml file changed. If it changed
you need to update your environment using the command for your operating system
```
conda env update -f requirements_linux.yml
```
```
conda env update -f requirements_windows.yml
```

When you install a new package into your environment you need to add it to the master package list. Then everyone
else can use the master list to update their environment by installing the package. Add the package name and version
to the requirements.txt file, this is a master list of all top level packages needed to rebuild the environment.

To save your current environment to file, first activate the environment and navigate to the top directory of the repository:
```
conda activate CU-PWFA
```
Next, use the command:
```
conda env export --no-builds > requirements_windows.yml
```
```
conda env export --no-builds > requirements_linuxs.yml
```


### Setup Spyder and IPython

To use spyder inside of the virtual environment, activate the virtual environment and launch spyder from the terminal.

To use IPython in a Jupyter-Notebook you need to install the IPython kernel in the environment. The package should have been downloaded 
when the virtual environment was created. To install use the command:
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


### VSim environment

Currently no one uses this. We generally analyze VSim results by directly loading data from the hdf5 files.
It is no longer maintained or updated.

VSim's analysis scripts are written in python 2.7. Unfortunately, everything we do is written in python 3.4. This means 
VSim scripts won't run in our virtual environment so we have a second virtual environment for VSim where you can write analyzer 
scripts or scripts to extract data from VSim results files. Unless you are writing an analyzer, I would suggest extracting
data using VSim's internal scripts and saving it to a numpy file in Python 2.7 and doing all analysis in Python 3.4.

There is a separate VSim file in the repository that contains the requirements file for the VSim virtual environment. Place all
scripts meant to run in this environment in this file.

To set up the VSim virtual environment, first cd into the VSim folder. Next, create the environment using 
```
conda create --name VSim --file requirements.txt
```
To add the environment to Ipython and Jupyter-Notebook, follow the same steps as above. To update the environment, follow the
steps above.
