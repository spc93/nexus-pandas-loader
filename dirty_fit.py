#import __main__
from __main__ import gca, plot, axis, xlim

def fit(func, aXis=None):
    '''
    fit(func): fit first line in current axis using quickfit function func
    fit(func, axis): do the same for the line in a specified axis
    fit is designed to be simple with limited functionality
    Use fit function directly for a more sophistcated fit
    Need to import pyplot first: from matplotlib.pyplot import *
    '''
    if aXis == None:
        aXis = gca()

    [xData, yData] = aXis.get_lines()[0].get_xydata().transpose()
    xL, xU = aXis.get_xlim()
    iI = (xData >= xL) & (xData <= xU)
    xData, yData = xData[iI], yData[iI]
    print func.fit(xData, yData)
    plot(xData, func(xData),'r.',hold=True); axis('tight'); xlim(xL, xU);

