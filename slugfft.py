import sys
import math
import numpy

class SLUG:
  def __init__(self,t0, lo=0.0,hi=100.0,dur=20.0, ramp=3.0,nsmooth=6, **kwargs):
    (self.t0,self.lo,self.hi,self.dur,self.rampup
    ,) = map(float,(t0, lo, hi, dur, ramp,))
    self.nsmooth = int(nsmooth)
    self.range = self.hi - self.lo 
    self.last,self.tlast = self.lo,0.0
    self.samples = [self.lo / self.nsmooth]*self.nsmooth
    self.isample = 0
    self.rampdown = self.dur - self.rampup
    self.end = self.dur + (self.rampup * 3.0)

  def val(self,t=None):
    if not (t is None): self.update(t)
    return self.last

  def notexpired(self):
    return self.tlast <= self.end

  def update(self,t):
    t1 = t - self.t0
    if t1 == self.tlast: return
    self.tlast = t1
    self.samples[self.isample] = self.nextvalraw(t1) / self.nsmooth
    self.isample = (self.isample + 1) % self.nsmooth
    self.last = sum(self.samples)

  def nextvalraw(self,t1):
    if t1 <= 0.0           : return self.lo
    if t1 <= self.rampup   : return self.lo + (self.range * t1 / self.rampup)
    if t1 <= self.rampdown : return self.hi
    if t1 <= self.dur      : return self.hi - (self.range * (t1 - self.rampdown)  / self.rampup)
    return self.lo

  def __gt__(self,other): 
    return self.last > other.last


class SLUGS:
  def __init__(self,period=25,seconds=512,**kwargs):
    self.slugs = list()
    self.period,self.seconds = int(period),int(seconds)
    self.pulses = self.seconds // self.period
    self.last = 0.0

  def addslugs(self,**kwargs):
    offset = (self.seconds - (self.pulses * self.period)) / 2.0
    for dt in range(self.pulses): self.addslug((dt*self.period)+offset,**kwargs)
    return self

  def addslug(self,t0,**kwargs):
    self.slugs.append(SLUG(float(t0),**kwargs))

  def update(self,t):
    for slug in self.slugs: slug.update(t)
    self.slugs = [slug for slug in self.slugs if slug.notexpired()]

  def val(self,t=None):
    if not (t is None): self.update(t)
    if not self.slugs: return self.last
    self.last = max(self.slugs).last
    return self.last


try   : import matplotlib.pyplot as plt
except: pass
def plotdata(ax,ts,vals):
  ax.plot(ts,vals)
  ax.set_title('Data')
  ax.set_ylabel('time')
  ax.set_xlabel('values')

def plotfft(ax,freqs,fft):
  magn = numpy.sqrt((fft.real*fft.real) + (fft.imag*fft.imag))
  ax.plot(freqs,fft.real,ls='dotted',label='real')
  ax.plot(freqs,fft.imag,ls='dotted',label='imag')
  ax.plot(freqs,magn,label='magn')
  ax.set_xlabel('frequency')
  ax.set_ylabel('FFT value')
  ax.legend(loc='best')

def plot(ts,vals,title):

  fig,(valax,fftax,) = plt.subplots(2,1,height_ratios=[1,6])

  plotdata(valax,ts,vals)

  freqs = numpy.fft.fftfreq(ts.shape[-1])
  fft = numpy.fft.fft(vals)
  plotfft(fftax,freqs,fft)

  valax.set_title(title)

  plt.show()

if "__main__" == __name__:
  kwargs = dict()
  for arg in sys.argv:
    if not arg.startswith('--'): continue
    toks = arg[2:].split('=')
    kwargs[toks[0]] = toks[1:] and '='.join(toks[1:]) or True

  slugs = SLUGS(**kwargs).addslugs(**kwargs)
  title = f"Samples={slugs.seconds}s; PulseDuration={slugs.slugs[0].dur}s; Period={slugs.period}s\nRamp={slugs.slugs[0].rampup}s; NSmooth={slugs.slugs[0].nsmooth}s"
  ts = numpy.arange(slugs.seconds,dtype=numpy.float64)
  vals = numpy.array([slugs.val(t) for t in ts])
  plot(ts,vals,title)
