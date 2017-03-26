# plasma-source

## Repository for plasma source development.

### Setup Anaconda and the virtual enviornment

First step is to download the repository. Depending on what operating system you are you using you can do this by
either downloading GitHub desktop (Mac, Windows) or by setting up a git repository (Linux). 

Once you have the repository you should install Anaconda3 which will install python and all of the scientific packages.

Anaconda will handle setting up your virtual enviornment. cd into the repository directory, it should contain a
requirements.txt file. Inside of the directory use the following command
```
conda create --name CU-PWFA --file requirements.txt
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

To use IPython you need to install the IPython kernel in the enviornment. The package should have been downloaded when
the virtual enviornment was created. To install use the command:
```
python -m ipykernel install --user --name=CU-PWFA
```
Now when you create an Jupyter-Notebook it should give you the option to use the CU-PWFA kernel.