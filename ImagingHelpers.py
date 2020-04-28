from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import matplotlib.pyplot as plt

class Image:
        
    def __init__(self, filename):
        self.filename = filename
        fts = fits.open(filename)
        self.hdu = fts[0]
        
        # If the header has the deprecated keyword, replace it
        # This will change it in both self.header and self.hdu
        try:
            self.hdu.header.rename_keyword('RADECSYS','RADESYSa')
        except:
            pass
        self.data = np.array(self.hdu.data.copy(),dtype='float64')
        self.header = self.hdu.header.copy()
        self.wcs = WCS(self.header)
        fts.close()

def rgb_scale(a, log=False, vmin=None, vmax=None):
    
    scaled = a.copy()
      
    if log:
        scaled = np.log10(scaled)

    if vmin is None:
        vmin = np.nanmin(scaled)
    else:
        scaled[scaled <= vmin] = vmin
    if vmax is None:
        vmax = np.nanmax(scaled)
    else:
        scaled[scaled >= vmax] = vmax  
    
    # Scale into [0,1]
    scaled = (scaled - vmin)/(vmax - vmin)
    
    # Deal with nan's
    scaled[~np.isfinite(scaled)] = 0
    # This should not be required given the "Scale into [0,1]" line, but sometimes there are roundoff errors
    scaled[scaled < 0] = 0.
    scaled[scaled > 1] = 1.
        
    return scaled

def rgb_stack(r, g, b):
    return np.stack((r, g, b), axis=2)

def minmax(array):
    print(array.min(), array.max())
    print(np.nanmin(array), np.nanmax(array))
    
def ColorSwatch(rgb):
    from ImagingHelpers import rgb_stack
    block = np.ones([25,25])
    image = rgb_stack(block*rgb[0], block*rgb[1], block*rgb[2])
    plt.imshow(image)
    ax = plt.gca()
    ax.axes.xaxis.set_visible(False)
    ax.axes.yaxis.set_visible(False)
    
    return
