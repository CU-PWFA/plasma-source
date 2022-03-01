#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 14:43:14 2022

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
import threading 
#load a signal 
#decide if that is ok
#if not 
 #delete the "signal" file
 #load one signal plot and sum with one axis plot
 #click the peak
 #let it run 
os.chdir(os.path.expanduser('~/Dropbox/01_Research/01_Project/01_PlasmaDianostics/01_Shadowgraphy/03_DataAnalysis/04_Functions/'))
from DataToSignal import DataToSignal
from MultishotAvgSeed import MultishotAvgSeed

os.chdir(os.path.expanduser('~/Dropbox/01_Research/01_Project/01_PlasmaDianostics/01_Shadowgraphy/03_DataAnalysis/03_LineFindingGUI/'))


root= tk.Tk()
root.title('Single Peak Finding Clicky Tool')
root.maxsize(1400, 600)

left_frame= tk.Frame(root, width= 500, height= 500, bg='white')
left_frame.grid(row=0, column=0, padx= 1, pady= 1)
right_frame= tk.Frame(root, width= 800, height= 500, bg='white')
right_frame.grid(row=0, column=1, padx= 1, pady= 1)
bottom_frame= tk.Frame(root, width= 1100, height= 100, bg='white')
bottom_frame.grid(row=1, column=0, padx= 1, pady= 1, columnspan=2)
left_canvas= tk.Canvas(left_frame, width= 500, height= 500)
right_canvas= tk.Canvas(right_frame, width= 800, height= 500)
bottom_canvas= tk.Canvas(bottom_frame, width= 1400, height= 600)

left_canvas.pack()
right_canvas.pack()
bottom_canvas.pack()


def draw_dot(event):
    global x1
    global y1
    x1=event.x
    y1=event.y
    x2=event.x
    y2=event.y
    # Draw an oval in the given co-ordinates
    left_canvas.create_oval(x1,y1,x2,y2,fill="black", width=5)

left_canvas.bind('<Button-1>', draw_dot)


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
    global LO
    global DataNumber
    
#    os.chdir(os.path.expanduser('~/Dropbox/01_Research/01_Project/01_PlasmaDianostics/01_Shadowgraphy/03_DataAnalysis/05_ParameterScan/01_Ar/01_PressureScan/02_MultishotsAvgResults/01_npy'))
    for name in glob.glob(os.path.expanduser('~/Dropbox/01_Research/01_Project/01_PlasmaDianostics/01_Shadowgraphy/03_DataAnalysis/05_ParameterScan/01_Ar/01_PressureScan/02_MultishotsAvgResults/01_npy/')+'??????????_NorSignal.npy'):
        DataNumber= name[172:182]
        if (DataNumber not in SkipList):
            SignalArray= np.load(name)
            img= ImageTk.PhotoImage(image= Image.fromarray(np.uint8(cm.bwr(SignalArray)*255)))
            LO= ImageTk.PhotoImage(Image.open(os.path.expanduser('~/Dropbox/01_Research/01_Project/01_PlasmaDianostics/01_Shadowgraphy/03_DataAnalysis/05_ParameterScan/01_Ar/01_PressureScan/02_MultishotsAvgResults/02_Plots/'+ \
                                                              DataNumber+'_ManualCheck.png')))
            left_canvas.create_image(0, 0, anchor='nw', image= img)
            right_canvas.create_image(0, 0, anchor='nw', image= LO)
            break
    try:
        bottom_canvas.delete('Confirm')
    except:
        pass

def Skip(event):
    SkipList.append(DataNumber)
    LoadSignal(event)
    
def LoadRawSignal(event):
    global SF
    global SignalFiltered
    global MediaPath
    global avgBG
    avgBG= np.load(os.path.expanduser('~/Dropbox/01_Research/01_Project/01_PlasmaDianostics/01_Shadowgraphy/03_DataAnalysis/05_ParameterScan/01_Ar/01_PressureScan/02_MultishotsAvgResults/01_npy/')+ DataNumber+ '_BG.npy')
    MediaPath= '/media/valentinalee/Elements/PlasmaDiagnosticsData/day_29/'
    OneShotData= MediaPath+ DataNumber+ '/17529184_'+DataNumber+'_0010.tiff'
#load one shot and do the finding signal thing
    SignalFiltered= DataToSignal(avgBG, OneShotData)
#show the one shot and the sum lineout    
    SF= ImageTk.PhotoImage(image= Image.fromarray(np.uint8(cm.bwr((SignalFiltered-np.amin(SignalFiltered))/np.amax(SignalFiltered-np.amin(SignalFiltered)))*255)))
    left_canvas.create_image(0, 0, anchor='nw', image= SF)
    bottom_canvas.delete('coordTxt')
    bottom_canvas.create_text(820, 50, fill="darkblue", font= "Times 16 italic bold", \
                              text='SingleShotSignal', tags= 'SSS')
    
    
def ReMatch(event):
    global ResultDir
    DataDir= MediaPath+ DataNumber +'/'
    ResultDir= '~/Dropbox/01_Research/01_Project/01_PlasmaDianostics/01_Shadowgraphy/03_DataAnalysis/05_ParameterScan/01_Ar/01_PressureScan/02_MultishotsAvgResults/'
    seed= y1
    threading.Thread(target= MultishotAvgSeed, args=(DataDir, avgBG, ResultDir, seed, )).start()
    try:
        bottom_canvas.delete('SSS')
    except:
        pass
    bottom_canvas.create_text(820, 50, fill="darkblue", font= "Times 16 italic bold", \
                              text='Dot Confirmed. Do not pressed Comfirmed Dot again.', tags= 'Confirm')
    SkipList.append(DataNumber)


#button0
buttonBG= bottom_canvas.create_rectangle(10, 10, 110, 40, fill= 'white', outline= 'grey60')
buttonTXT= bottom_canvas.create_text(60, 25, text= 'LoadSignal')
bottom_canvas.tag_bind(buttonBG, "<Button-1>", LoadSignal) ## when the square is clicked runs function "clicked".
bottom_canvas.tag_bind(buttonTXT, "<Button-1>", LoadSignal) ## same, but for the text.

#button0.1
buttonBG= bottom_canvas.create_rectangle(120, 10, 220, 40, fill= 'white', outline= 'grey60')
buttonTXT= bottom_canvas.create_text(170, 25, text= 'Good')
bottom_canvas.tag_bind(buttonBG, "<Button-1>", Skip) ## when the square is clicked runs function "clicked".
bottom_canvas.tag_bind(buttonTXT, "<Button-1>", Skip) ## same, but for the text.

#button0.2
buttonBG= bottom_canvas.create_rectangle(230, 10, 330, 40, fill= 'white', outline= 'grey60')
buttonTXT= bottom_canvas.create_text(280, 25, text= 'Bad')
bottom_canvas.tag_bind(buttonBG, "<Button-1>", LoadRawSignal) ## when the square is clicked runs function "clicked".
bottom_canvas.tag_bind(buttonTXT, "<Button-1>", LoadRawSignal) ## same, but for the text.

#button0.3
buttonBG= bottom_canvas.create_rectangle(340, 10, 440, 40, fill= 'white', outline= 'grey60')
buttonTXT= bottom_canvas.create_text(390, 25, text= 'Hopeless')
bottom_canvas.tag_bind(buttonBG, "<Button-1>", Skip) ## when the square is clicked runs function "clicked".
bottom_canvas.tag_bind(buttonTXT, "<Button-1>", Skip) ## same, but for the text.

#button1
buttonBG= bottom_canvas.create_rectangle(10, 50, 110, 80, fill= 'white', outline= 'grey60')
buttonTXT= bottom_canvas.create_text(60, 65, text= 'ConfirmDot')
bottom_canvas.tag_bind(buttonBG, "<Button-1>", ReMatch) ## when the square is clicked runs function "clicked".
bottom_canvas.tag_bind(buttonTXT, "<Button-1>", ReMatch) ## same, but for the text.

#button2
#buttonBG= bottom_canvas.create_rectangle(120, 50, 220, 80, fill= 'white', outline= 'grey60')
#buttonTXT= bottom_canvas.create_text(170, 65, text= 'ConfirmLine2')
#bottom_canvas.tag_bind(buttonBG, "<Button-1>", FindPeaksAndPlot2) ## when the square is clicked runs function "clicked".
#bottom_canvas.tag_bind(buttonTXT, "<Button-1>", FindPeaksAndPlot2) ## same, but for the text.

#button3
buttonBG= bottom_canvas.create_rectangle(230, 50, 330, 80, fill= 'white', outline= 'grey60')
buttonTXT= bottom_canvas.create_text(280, 65, text= 'ClearAll')
bottom_canvas.tag_bind(buttonBG, "<Button-1>", clearAll) ## when the square is clicked runs function "clicked".
bottom_canvas.tag_bind(buttonTXT, "<Button-1>", clearAll) ## same, but for the text.

#button3
#buttonBG= bottom_canvas.create_rectangle(340, 50, 440, 80, fill= 'white', outline= 'grey60')
#buttonTXT= bottom_canvas.create_text(390, 65, text= 'SavePeaks')
#bottom_canvas.tag_bind(buttonBG, "<Button-1>", SavePeak) ## when the square is clicked runs function "clicked".
#bottom_canvas.tag_bind(buttonTXT, "<Button-1>", SavePeak) ## same, but for the text.

#if it's a good signal, then skip
root.mainloop()


