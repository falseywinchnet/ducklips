import sys, wave, numpy
name,t0,t1=sys.argv[1],float(sys.argv[2]),float(sys.argv[3])
h=wave.open(name); d=numpy.frombuffer(h.readframes(h.getnframes()),dtype='<i2').astype(float)
seg=d[int(t0*48000):int(t1*48000)]
S=numpy.abs(numpy.fft.rfft(seg*numpy.hanning(len(seg))))**2; f=numpy.fft.rfftfreq(len(seg),1/48000)
def band(a,b): return 10*numpy.log10(S[(f>=a)&(f<b)].sum()+1e-9)
top=band(4000,12000)
print(name.split('/')[-1], " ".join("%d-%d:%.0f"%(a,b,band(a,b)-top) for a,b in [(0,1000),(1000,2500),(2500,4000),(4000,6000),(6000,9000),(9000,14000)]))
