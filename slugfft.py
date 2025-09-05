import sys
import math
import numpy

class SLUG:
  """Model a single pulse as a smoothed trapezoid
Arguments;
       t0:  Absolute time of start of pulse, start of low to high ramp
       lo:  low value of pulse value, default 0.0
       hi:  high value of pulse value, default 100.0
      dur:  duriation of pulse including ramps, default 20s
     ramp:  duration of each ramp (lo to hi at start; hi to lo at end), defult 3s
  nsmooth:  size of moving average, default 6samples
"""
  def __init__(self,t0, lo=0.0,hi=100.0,dur=20.0, ramp=3.0,nsmooth=6, **kwargs):
    """Initialize values"""
    (self.t0,self.lo,self.hi,self.dur,self.rampup
    ,) = map(float,(t0, lo, hi, dur, ramp,))
    self.nsmooth = int(nsmooth)
    self.range = self.hi - self.lo 
    ### Last smoothed value calculated; relative time of last update
    self.last,self.tlast = self.lo,0.0
    ### Initialize moving average array
    ### - .isample is Next sample number
    self.samples = [self.lo / self.nsmooth]*self.nsmooth
    self.isample = 0
    ### Time when ramp down (hi to lo) starts
    self.rampdown = self.dur - self.rampup
    ### Time when pulse no longer exists
    self.end = self.dur + self.nsmooth + 1

  def val(self,t=None):
    """Update to time t; return the latest value"""
    if not (t is None): self.update(t)
    return self.last

  def active(self):
    """Return True when active; false when inactive"""
    return self.tlast <= self.end

  def update(self,t):
    """Update to new absolute time, t"""
    ### Convert to relative time; do nothing if no change
    trel = t - self.t0
    if trel == self.tlast: return
    ### Save relative time in .tlast
    ### Place new raw value into moving average array
    ### Increment the index into the moving average array, with wraparound
    self.tlast = trel
    self.samples[self.isample] = self.nextvalraw(trel) / self.nsmooth
    self.isample = (self.isample + 1) % self.nsmooth
    ### Calculate smoothed value (moving average)
    self.last = sum(self.samples)

  def nextvalraw(self,trel):
    """Get the unsmoothed value at a relative time"""
    ### Use .lo before the pulse starts
    if trel <= 0.0           : return self.lo
    ### Calculate ramp-up value from 0 to .rampup
    if trel <= self.rampup   : return self.lo + (self.range * trel / self.rampup)
    ### Use self.hi from .rampup to .rampdown
    if trel <= self.rampdown : return self.hi
    ### Calculate ramp-down value from .rampdown to .dur
    if trel <= self.dur      : return self.hi - (self.range * (trel - self.rampdown)  / self.rampup)
    ### Use .lo after the pulse ends
    return self.lo

  def __gt__(self,other): 
    """Supply greater-than mehtod for For max(SLUGS)"""
    return self.last > other.last


class SLUGS:
  """Maintain an array of periodic pulse models;
  - return highest value of overlapping pulses
Arguments
   period:  time betwween start of successive pulses
  seconds:  number of 1Hz samples that will be calculated
"""
  def __init__(self,period=25,seconds=512,**kwargs):
    """Create list to hold pulse models, save period and seconds
  - calculate number of pulse models
  - initialzze .last value for when no slugs are left
"""
    self.slugs = list()
    self.period,self.seconds = int(period),int(seconds)
    self.pulses = self.seconds // self.period
    self.last = 0.0

  def addslugs(self,**kwargs):
    """Populate list of pulse models based on .period and .seconds"""
    offset = (self.seconds - (self.pulses * self.period)) / 2.0
    for dt in range(self.pulses): self.addslug((dt*self.period)+offset,**kwargs)
    return self

  def addslug(self,t0,**kwargs):
    """Append one pulse model to list"""
    self.slugs.append(SLUG(float(t0),**kwargs))

  def update(self,t):
    """Update all pulse models to absolute time t"""
    for slug in self.slugs: slug.update(t)
    ### Keep only active models
    self.slugs = [slug for slug in self.slugs if slug.active()]

  def val(self,t=None):
    """Return the highest value of all active models"""
    if not (t is None): self.update(t)
    ### If pulse model list is empty, return the last calculated value
    if not self.slugs: return self.last
    ### Save and return the highest value from all active pulse models
    self.last = max(self.slugs).last
    return self.last


try   : import matplotlib.pyplot as plt
except: pass
def plotdata(ax,ts,vals):
  """Plot the time domain data"""
  ax.plot(ts,vals)
  ax.set_title('Data')
  ax.set_ylabel('time')
  ax.set_xlabel('values')

def plotfft(ax,freqs,fft):
  """Plot the FFT data"""
  magn = numpy.sqrt((fft.real*fft.real) + (fft.imag*fft.imag))
  ax.plot(freqs,fft.real,ls='dotted',label='real')
  ax.plot(freqs,fft.imag,ls='dotted',label='imag')
  ax.plot(freqs,magn,label='magn')
  ax.set_xlabel('frequency')
  ax.set_ylabel('FFT value')
  ax.legend(loc='best')

def plot(ts,vals,title):
  """Calculate FFT; plot all data"""

  freqs = numpy.fft.fftfreq(ts.shape[-1])
  fft = numpy.fft.fft(vals)

  fig,(valax,fftax,) = plt.subplots(2,1,height_ratios=[1,6])

  plotdata(valax,ts,vals)
  plotfft(fftax,freqs,fft)

  valax.set_title(title)

  plt.show()

if "__main__" == __name__:
  ### Convert --name=value command-line arguments to keyword dict
  kwargs = dict()
  for arg in sys.argv:
    if not arg.startswith('--'): continue
    toks = arg[2:].split('=')
    kwargs[toks[0]] = toks[1:] and '='.join(toks[1:]) or True

  ### Create SLUGS object and list of pulst models
  slugs = SLUGS(**kwargs).addslugs(**kwargs)

  ### Build the plot title
  title = f"Samples={slugs.seconds}s; PulseDuration={slugs.slugs[0].dur}s; Period={slugs.period}s\nRamp={slugs.slugs[0].rampup}s; NSmooth={slugs.slugs[0].nsmooth}s"

  ### build array of times, use SLUGs object to calculate values from those times
  ts = numpy.arange(slugs.seconds,dtype=numpy.float64)
  vals = numpy.array([slugs.val(t) for t in ts])

  ### Plot the time- and FFT-data
  plot(ts,vals,title)
