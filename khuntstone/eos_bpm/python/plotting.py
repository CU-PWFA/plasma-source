'''
Module for creating figure and axis handles
'''

from cycler import cycler;
import matplotlib.pyplot as plt;

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
	
