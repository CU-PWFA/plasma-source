# plasma-source

## Repository for plasma source development.

### Setup Anaconda and the virtual enviornment

First step is to download the repository. Depending on what operating system you are you using you can do this by
either downloading GitHub desktop (Mac, Windows) or by setting up a git repository (Linux). 

Once you have the repository you should install Anaconda3 which will install python and all of the scientific packages.

Anaconda will handle setting up your virtual enviornment. cd into the repository directory, it should contain a
requirements.txt file. Inside of the directory use the following command
```
conda create --name CU-PWFA
```
This will create a virtual enviornment called "CU-PWFA" with all the modules specified in the requirements.txt file.
To activate the enviornment:

(Mac, Linux):
```
source activate CU-PWFA
```
(Windows)
```
activate CU-PWFA
```

### Setup Spyder and IPython

To use spyder inside of the virtual enviornment, enter the virtual enviornment and then launch spyder from the command
line. 

To use IPython in a Jupyter-Notebook you need to install the IPython kernel in the enviornment. The package should have been downloaded when
the virtual enviornment was created. To install use the command:
```
python -m ipykernel install --user --name=CU-PWFA
```
Now when you create an Jupyter-Notebook it should give you the option to use the CU-PWFA kernel.

### Update the virtual enviornment

When you install a new package into your enviornment you need to add it to the master package list. Then everyone
else can use the master list to update their enviornment by installing the package.

To save your current enviornment to file, first activate the enviornment:
```
source activate CU-PWFA
```
Next, use the command:
```
conda env export > requirements.txt
```

Whenever you update your local repository using the pull command, check to see if the requirement.txt file changed. If it changed
you need to update your enviornment using
```
conda env update -f requirements.txt
```

### VSim enviornment

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

