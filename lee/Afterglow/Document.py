#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 16:09:52 2021

@author: valentinalee
"""

print('''
Document for running afterglow scan

Check USim simulation by USim software. Once it runs fine in USim GUI and is realy to run scan, save all the USim files in a dir, called, let's say "FinalCheck".
Go to USim4.0 folder to open "RunFromCommendLine.py" and change the dir to where "FinalCheck" is.
USim4.0 folder is located at ~/Softwares/USim-4.0 if running from my laptop. 
Once you run "RunFromCommendLine.py", newly generated USim Scan (dir) will locate in "FinalCheck" dir, having names like n50_GP5_PW70.
You can cp the newly generated dirs to another dir, let's say "AfterglowScanFinal"
Now, you need to build up the folders that is required for scan. 
AfterglowFinal-|ParameterScanPlotting.py
               |ExpansionPlots 
               |n50_GP5_PW70
               |...
               |ScanResults------|AfterglowPlots
                                 |AfterglowScan.py
                                 |AfterglowScan
             
ParameterScanPlotting.py and AfterglowScan.py can be found at /media/valentinalee/Transcend1/ValentinaLee/Scripts/Afterglow

Once the env is built, run ParameterScanPlotting.py
This will generate the python results from USim and plot it. The plot will be saved in "ExpansionPlots".
The py results will be saved in "ScanResults". 
Note that if "ScanResults" contains a given result, it will be skipped when running "ParameterScanPlotting.py". 
It is, when there are many "n50_GP5_PW70" type of USim scan in the folder, the program will only run through the one that has not been run before.
Thus, if you would like to run some of them for the second time because you updated its "n50_GP5_PW70" file, remember to delete the generated scan results. 
(Plots does not need to be deleted. It over write every time)
Check the scan range in "ParameterScanPlotting.py" before running it. 

Once there are some files in "ScanResults", you can go ahead and run  "AfterglowScan.py", which will calculate the afterglow and save the results in "AfterglowScan" and the plot in "AfterglowPlots".
"AfterglowScan.py" works same as "ParameterScanPlotting.py" that it only goes through the one that has not be process beofore. 
Check the scan range in "AfterglowScan.py" before running it. 

If you want to change the recomination coeff, it is in "AfterglowScan.py". (Obviously)

***My code goes back to look at USim's .pre file to see what is the time and space bound and grid. Do not change the space and character number when change the grid size. 
***"RunFromCommendLine.py" makes the script being run in parallel. It is currently seup as 6 core. If you would like to change the core number, go to the last line of the code, change the number 6 to whatever you would like.

Plotting final reuslts:
The results are those "n5_GP3_PW30_map.npy"s and "n5_GP3_PW30_x.npy" in "AfterglowScan" dir. You can plot them by "AfterglowCalculation.py" at /media/valentinalee/Transcend1/ValentinaLee/Scripts/Afterglow



''')
