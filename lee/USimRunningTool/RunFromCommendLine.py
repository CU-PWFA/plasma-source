#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 14:05:37 2020

@author: valentinalee
"""


#This script needs to run in USim-4.0 dir
import os   
import shutil

for density in range (1, 3, 2):
    
    for power in range (5, 6):
        
        for pw in range (30, 50, 10):

            BasicPrePath= "../../Littlesun_LT/USim/EPExchange/"
            FileName= "cap_2D"
            NewDir= 'n'+str(density)+'_GP'+str(power)+'_PW'+str(pw)+'/'
            
            path= BasicPrePath+NewDir
            try:
                os.mkdir(path)
            except OSError:
                print ("Creation of the directory %s failed" % path)
            else:
                print ("Successfully created the directory %s " % path)
            
            BasicFile= BasicPrePath+ FileName+ '.pre'
            shutil.copy(BasicFile, path)
            
            NewPreFile= path+FileName+'.pre'
            reading_file = open(NewPreFile, "r")
            
            new_file_content = ""
            for line in reading_file:
                stripped_line = line.strip()
                new_line1 = stripped_line.replace("exp(-(r^2/(30e-6)^2)^2.6)", "exp(-r^"+str(power)+"/(2*"+str(pw)+"^"+str(power)+"))")
                new_file_content += new_line1+"\n"

            reading_file.close()
            writing_file = open(NewPreFile, "w")
            writing_file.write(new_file_content)
            writing_file.close()

            reading_file = open(NewPreFile, "r")
            
            new_file_content = ""
            for line in reading_file:
                stripped_line = line.strip()
#                new_line1 = stripped_line.replace("$ N_FRAME = 30", "$ N_FRAME = 2")
                new_line1 = stripped_line.replace("$ N0 = 1e23", "$ N0 = "+str(density)+"e22")
                new_file_content += new_line1+"\n"

            reading_file.close()
            writing_file = open(NewPreFile, "w")
            writing_file.write(new_file_content)
            writing_file.close()
            
            
            os.system('./Contents/engine/bin/txpp.py '+NewPreFile)
            NewInFile= path+FileName+'.in'
            os.system('./Contents/engine/bin/mpiexec -np 6 ./Contents/engine/bin/ulixes -i '+NewInFile)