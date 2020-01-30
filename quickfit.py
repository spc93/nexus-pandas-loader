from __main__ import gca, plot, axis, xlim
from lmfit import Model
import numpy as np

def peak(xdat, ydat, nbgpts=1):
    '''
    [centre, fwhm_sd, fwhm_area, sum, height, area, m, c]=peak(x,y,nbgpts=1)
    Returns centre, fwhm_sd, fwhm_area, sum, height, area after subtracting a sloping
    linear background from the fist and nbgpts data points and the background parameters (m, c)
    fwhm_sd is from standard deviation (correct for gaussian)
    fwhm_area is from area/height (correct for gaussian) - more robust than sd with noisy data
    xdat, ydat are 1d arrays or lists
    '''
    x = np.array(xdat); 
    y = np.array(ydat); 
    npts = len(x); 
    xspan= x[-1] - x[0]			#copy data into arrays
    m = (np.mean(y[npts-nbgpts:npts]) - np.mean(y[0:nbgpts]))/(np.mean(x[npts-nbgpts:npts]) - np.mean(x[0:nbgpts]))    #slope for linear b/g
    c = np.mean(y[0:nbgpts]) - np.mean(x[0:nbgpts])*m	#intercept
    y = y - m * x - c;			#subtract background
    sumdat = sum(y);			#peak sum
    area = sumdat * xspan/npts;		#peak integral
    centre = sum(x*y)/sumdat;		#centroid calc.
    height = max(y);			#max y value after linear b/g
    fwhm_sd = np.sqrt(sum((x - centre)**2*y)/sumdat) * np.sqrt(8*np.log(2))
    fwhm_area = area/height * 0.3989 * np.sqrt(8*np.log(2))
    return [centre, fwhm_sd, fwhm_area, sumdat, height, area, m, c]    
 
### some pre-defined peak and background functions

def gauss(x, area, cen, fwhm):
    return area/fwhm*np.sqrt(4*np.log(2)/np.pi)*np.exp(-4*np.log(2)*((x-cen)/fwhm)**2)

def lorentz(x, area, cen, fwhm):
    return area/fwhm/(np.pi/2)/(1+4*((x-cen)/fwhm)**2)

def pvoigt(x, area, cen, fwhm, lfrac):
    return area/fwhm/(lfrac*np.pi/2+(1-lfrac)*np.sqrt(np.pi/4/np.log(2)))*(lfrac/(1+4*((x-cen)/fwhm)**2)+(1-lfrac)*np.exp(-4*np.log(2)*((x-cen)/fwhm)**2))

def poly2(x, m, c):
   return m * x + c

def const(x, c):
   return c

g_c = Model(gauss) + Model(const)
g_lin = Model(gauss) + Model(poly2)
lor_c = Model(lorentz) + Model(const)
lor_lin = Model(lorentz) + Model(poly2)
pv_c = Model(pvoigt) + Model(const)
pv_lin = Model(pvoigt) + Model(poly2)



class fit():
    '''
    create fit instance (i.e. do a fit) from peak-like plot data
    variables cen, fwhm, sum, amp, area, m, c are given initial values from peak function; others are zero by default
    optionally supply axis and/or parameters
    .result and .params contain results and parameters
    see lmfit documentation for information about constraining parameters etc
    for anything else use lmfit directly
    example:
      ff = fit(pv_c)                      # do a fit of current axis data using pv_c lmfit model
      pin = pv_c.make_params()            # create fresh parameters from model
      pin = ff.params                     # create parameters from fit
      pin['c'].set(50000, vary = False)   # set parameter to value and keep it fixed
      ff = fit(pv_c, params = pin)        # use these parameters (with one fixed) in fit_old
      ff.result                           # show full result (rich display in ipython)
      ff.params                           # return lmfit parameters (rich display in ipython
      pv_c.fit(y, pin, x=x)               # fit directly using lmfit model
    '''
    
    def __init__(self, func, aXis = None, params = None):

        if aXis == None:
            aXis = gca()

        [xData, yData] = aXis.get_lines()[0].get_xydata().transpose()
        xL, xU = aXis.get_xlim()
        iI = (xData >= xL) & (xData <= xU)
        xData, yData = xData[iI], yData[iI]
        
        pk_prms = {}
        [pk_prms['cen'], fwhm_sd, pk_prms['fwhm'], pk_prms['sum'], pk_prms['amp'], pk_prms['area'], pk_prms['m'], pk_prms['c']] = peak(xData, yData)
        
        if params == None:
          self.params = func.make_params()
          # assign fit parameter value to parameters that match the parameters from peak with others defaulting to zero
          for key in self.params.keys():
              try:
                  self.params[key].value = pk_prms[key]
              except:
                  self.params[key].value = 0
        else:
          # use parameters supplied (params) if given
          self.params = params

        self.result =  func.fit(yData, x=xData, params = self.params)   #do the fit 

        outstr = func.name+'\n\n'
        for pname in self.result.params:
            prm = self.result.params[pname]
            try:
                stderr = prm.stderr
                _n_dec_places = max(0, int(-np.round(np.log10(stderr))) + 1)
            except:
                stderr = 0
                _n_dec_places = 4 # in case parameter not varied then stderr is zero

            _fmt='%.' + str(_n_dec_places) + 'f'
            outstr+='%10s:  ' % prm.name + '%15s' % (_fmt % prm.value)+' +/- '+ '%-10s' % str(_fmt % stderr) +'\n'

        print(outstr+'\n')
    
        plot(xData, func.eval(x=xData, params=self.result.params),'r.'); axis('tight'); xlim(xL, xU);

