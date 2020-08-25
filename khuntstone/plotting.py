'''
Module containing often used plotting functions
'''
import matplotlib.pyplot as plt

def makefig(x = 8, y = 6, dpi = 200, xlab = '', ylab ='', titl =''):
    """
    Function to make a figure and axis handle
    
    Parameters:
    -----------
    x    : int
           Width of the plot (in)
    y    : int 
           height of the plot (in)
    dpi  : int
           Plot resolution 
    xlab : str, optional
           x axis label
    ylab : str, optional
           y axis label
    titl : str, optional
           plot title
    
    Returns:
    --------
    fig, ax : the figure and axis handles
    """
    
    fig = plt.figure(figsize = (x,y), dpi = dpi)
    ax  = fig.gca()
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(titl)
    
    return fig, ax
