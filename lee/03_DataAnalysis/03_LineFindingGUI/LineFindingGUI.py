#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 21:06:22 2021

@author: valentinalee
"""

#%%
import tkinter as tk
import numpy as np
from PIL import Image, ImageTk
from matplotlib import cm
import os
import glob
from os import path

os.chdir(os.path.expanduser('~/Dropbox/01_Research/01_Project/01_PlasmaDianostics/01_Shadowgraphy/03_DataAnalysis/04_Functions/'))
from PeakFinding import FindPeaksFromSeed

os.chdir(os.path.expanduser('~/Dropbox/01_Research/01_Project/01_PlasmaDianostics/01_Shadowgraphy/03_DataAnalysis/03_LineFindingGUI/'))


root= tk.Tk()
root.title('Peaks Finding Clicky Tool')
root.maxsize(1100, 600)

left_frame= tk.Frame(root, width= 500, height= 500, bg='white')
left_frame.grid(row=0, column=0, padx= 1, pady= 1)
right_frame= tk.Frame(root, width= 500, height= 500, bg='white')
right_frame.grid(row=0, column=1, padx= 1, pady= 1)
bottom_frame= tk.Frame(root, width= 1100, height= 100, bg='white')
bottom_frame.grid(row=1, column=0, padx= 1, pady= 1, columnspan=2)
left_canvas= tk.Canvas(left_frame, width= 500, height= 500)
right_canvas= tk.Canvas(right_frame, width= 500, height= 500)
bottom_canvas= tk.Canvas(bottom_frame, width= 1100, height= 600)

#TODO Make names for directory from the top 
#TODO Make a new function so that when finish everything in the folder, tell user you are done
#canvas= tk.Canvas(root, width= 1100, height= 600)
left_canvas.pack()
right_canvas.pack()
bottom_canvas.pack()


coords= {'x0':0,'y0':0,'x1':0,'y1':0}
# keep a reference to all lines by keeping them in a list 
lines= []

def click(e):
    # define start point for line
    coords['x0'] = e.x
    coords['y0'] = e.y

    # create a line on this point and store it in the list
    lines.append(left_canvas.create_line(coords['x0'], coords['y0'], coords['x0'], coords['y0'], tags= 'createdLine'))

def drag(e):
    # update the coordinates from the event
    coords['x1'] = e.x
    coords['y1'] = e.y

    # Change the coordinates of the last created line to the new coordinates
    left_canvas.coords(lines[-1], coords['x0'],coords['y0'],coords['x1'],coords['y1'])

left_canvas.bind("<ButtonPress-1>", click)
left_canvas.bind("<B1-Motion>", drag) 

FirstPeaks1= np.array([])
FirstPeaks2= np.array([])
SkipList= []

def printCoord(event):
    bottom_canvas.delete('coordTxt')
    bottom_canvas.create_text(820, 50, fill="darkblue", font= "Times 16 italic bold", text='Coord:'+str(coords), tags= 'coordTxt')
    
def clearAll(event):
    left_canvas.delete('createdLine')
    try:
        right_canvas.delete('Peaks1')
    except:
        pass
    try:
        right_canvas.delete('Peaks2')
        global FirstPeaks1
        global FirstPeaks2
        FirstPeaks1= np.array([])
        FirstPeaks2= np.array([])
    except:
        pass
    try:
        bottom_canvas.delete('coordTxt')
    except:
        pass
    try:
        bottom_canvas.delete('PFStxt')
    except:
        pass
    try:
        bottom_canvas.delete('PFSerror')
    except:
        pass

def FindPeaksAndPlot1(event):
    printCoord(event)
    global FirstPeaks1
    FirstPeaks1= FindPeaksFromSeed(coords, SignalArray)
    for xy in FirstPeaks1:
        right_canvas.create_oval(xy[0], xy[1], xy[0], xy[1], fill= 'Black', tags= 'Peaks1')
    
def FindPeaksAndPlot2(event):
    printCoord(event)
    global FirstPeaks2
    FirstPeaks2= FindPeaksFromSeed(coords, SignalArray)
    for xy in FirstPeaks2:
        right_canvas.create_oval(xy[0], xy[1], xy[0], xy[1], fill= 'Black', tags= 'Peaks2')

def SavePeak(event):
    bottom_canvas.delete('coordTxt')
    if FirstPeaks1.size!= 0 and FirstPeaks2.size!= 0:
        
        LongerLen= max([FirstPeaks1.shape[0], FirstPeaks2.shape[0]])
        CombinedResult= np.zeros((LongerLen, 4))
        CombinedResult[0:FirstPeaks1.shape[0], 0:2]= FirstPeaks1
        CombinedResult[0:FirstPeaks2.shape[0], 2:4]= FirstPeaks2
        np.save(os.path.join(os.path.expanduser('~/Dropbox/01_Research/01_Project/01_PlasmaDianostics/01_Shadowgraphy/03_DataAnalysis/05_ParameterScan/01_Ar/01_PressureScan/03_PeaksFoundByClickyTool/'), DataNumber+ '_PeaksInfo.npy'), \
                CombinedResult)
        bottom_canvas.create_text(820, 50, fill="darkblue", font= "Times 16 italic bold", text='Peaks Found Saved', tags= 'PFStxt')
    else:
        print('Error: FirstPeaks1 and/or FirstPeaks2 not defined.')
        bottom_canvas.create_text(820, 50, fill="darkblue", font= "Times 16 italic bold", text='Error: FirstPeaks1 and/or FirstPeaks2 not defined.', tags= 'PFSerror')
        return

def LoadSignal(event):
    global SignalArray
    global img
    global DataNumber
    
    os.chdir(os.path.expanduser('~/Dropbox/01_Research/01_Project/01_PlasmaDianostics/01_Shadowgraphy/03_DataAnalysis/05_ParameterScan/01_Ar/01_PressureScan/02_MultishotsAvgResults/01_npy'))
    for name in glob.glob('??????????_NorSignal.npy'):
        DataNumber= name[0:10]
        CheckRun= os.path.expanduser('~/Dropbox/01_Research/01_Project/01_PlasmaDianostics/01_Shadowgraphy/03_DataAnalysis/05_ParameterScan/01_Ar/01_PressureScan/03_PeaksFoundByClickyTool/')+ \
                    DataNumber+ '_PeaksInfo.npy'
        if path.exists(CheckRun)== False:
            if DataNumber not in SkipList:
                print(name)
                SignalArray= np.load(name)
                img= ImageTk.PhotoImage(image= Image.fromarray(np.uint8(cm.bwr(SignalArray)*255)))
                left_canvas.create_image(0, 0, anchor='nw', image= img)
                right_canvas.create_image(0, 0, anchor='nw', image= img)
                break
    try:
        bottom_canvas.delete('PFStxt')
    except:
        pass

def Skip(event):
    SkipList.append(DataNumber)
    LoadSignal(event)
    
#button0
buttonBG= bottom_canvas.create_rectangle(10, 10, 110, 40, fill= 'white', outline= 'grey60')
buttonTXT= bottom_canvas.create_text(60, 25, text= 'LoadSignal')
bottom_canvas.tag_bind(buttonBG, "<Button-1>", LoadSignal) ## when the square is clicked runs function "clicked".
bottom_canvas.tag_bind(buttonTXT, "<Button-1>", LoadSignal) ## same, but for the text.

#button0.1
buttonBG= bottom_canvas.create_rectangle(120, 10, 220, 40, fill= 'white', outline= 'grey60')
buttonTXT= bottom_canvas.create_text(170, 25, text= 'Skip')
bottom_canvas.tag_bind(buttonBG, "<Button-1>", Skip) ## when the square is clicked runs function "clicked".
bottom_canvas.tag_bind(buttonTXT, "<Button-1>", Skip) ## same, but for the text.

#button1
buttonBG= bottom_canvas.create_rectangle(10, 50, 110, 80, fill= 'white', outline= 'grey60')
buttonTXT= bottom_canvas.create_text(60, 65, text= 'ConfirmLine1')
bottom_canvas.tag_bind(buttonBG, "<Button-1>", FindPeaksAndPlot1) ## when the square is clicked runs function "clicked".
bottom_canvas.tag_bind(buttonTXT, "<Button-1>", FindPeaksAndPlot1) ## same, but for the text.

#button2
buttonBG= bottom_canvas.create_rectangle(120, 50, 220, 80, fill= 'white', outline= 'grey60')
buttonTXT= bottom_canvas.create_text(170, 65, text= 'ConfirmLine2')
bottom_canvas.tag_bind(buttonBG, "<Button-1>", FindPeaksAndPlot2) ## when the square is clicked runs function "clicked".
bottom_canvas.tag_bind(buttonTXT, "<Button-1>", FindPeaksAndPlot2) ## same, but for the text.

#button3
buttonBG= bottom_canvas.create_rectangle(230, 50, 330, 80, fill= 'white', outline= 'grey60')
buttonTXT= bottom_canvas.create_text(280, 65, text= 'ClearAll')
bottom_canvas.tag_bind(buttonBG, "<Button-1>", clearAll) ## when the square is clicked runs function "clicked".
bottom_canvas.tag_bind(buttonTXT, "<Button-1>", clearAll) ## same, but for the text.

#button3
buttonBG= bottom_canvas.create_rectangle(340, 50, 440, 80, fill= 'white', outline= 'grey60')
buttonTXT= bottom_canvas.create_text(390, 65, text= 'SavePeaks')
bottom_canvas.tag_bind(buttonBG, "<Button-1>", SavePeak) ## when the square is clicked runs function "clicked".
bottom_canvas.tag_bind(buttonTXT, "<Button-1>", SavePeak) ## same, but for the text.


root.mainloop()


