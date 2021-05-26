'''
Module for creating figure and axis handles
'''

from cycler import cycler;
import matplotlib.pyplot as plt;
import numpy as np

def makefig(x = 6, y = 4, xlab = '', ylab = '', title = '', fs = 12, ts = 12):
	'''
	Function to create custom figure and axis handles. 
	Parameters:
	-----------
	x : int, optional
	    Horizontal size of the plot, in. 
	y : int, optional 
	    Vertical size of the plot, in. 
	xlab : str, optional
	       xlabel of the plot
	ylab : str, optional
	       ylabel pf the plot
	title : str, optional
	        Title of the plot
	fs    : int, optional
	        Fontsize for axes labels
	ts    : Fontsize for title

	Returns:
	--------
	fig : object
	      pyplot figure handle
	ax  : object
	      pyplot axes handle
	'''

	# Custom colors for plotting
	plot_colors = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933',\
	               '#DDCC77', '#CC6677', '#882255', '#AA4499'];
	cy  = cycler('color', plot_colors);
	fig = plt.figure(figsize = (x, y), dpi = 600);
	ax  = fig.gca();
	ax.set_prop_cycle(cy);
	ax.set_xlabel(xlab, fontsize = fs);
	ax.set_ylabel(ylab, fontsize = fs);
	ax.set_title(title, fontsize = ts, fontweight = 'bold');
	ax.tick_params(labelsize = 'large');

	return fig, ax
	

def plot_signal(E, te, sig, tsig):
    fig = plt.figure(figsize = (4,4), dpi = 200)
    ax1 = fig.gca()
    ax1.spines['left'].set_color("red")
    ax1.set_ylabel("E [V/m]", color = "red")
    ax1.yaxis.label.set_color("red")
    ax1.tick_params(axis="y", colors="red")
    ax1.set_xlabel('t [ps]')
    ax1.set_xlim([min(tsig*1e12), max(tsig*1e12)])
    ax1.plot(te*1e12, E, '-r')
    
    ax2 = ax1.twinx()
    ax2.spines['left'].set_color("blue")
    ax2.set_ylabel("Signal [AU]", color = "blue")
    ax2.yaxis.label.set_color("blue")
    ax2.tick_params(axis="y", colors="blue")
    ax2.set_xlabel('t [ps]')
    ax2.plot(tsig*1e12, sig/max(sig), '--b')
    
    # creat appropriate limits 
    x1 = te[np.argmax(E)]*1e12 - 0.15
    x2 = 0.7
    #ax1.set_xlim([x1, x2])
    
    plt.show()