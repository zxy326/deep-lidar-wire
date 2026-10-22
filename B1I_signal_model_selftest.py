"""Fractional-delay and symbol-boundary regression, no MATLAB required."""
from pathlib import Path
import numpy as np
import B1I_simulate as sim

def run():
    code=sim.primary(19,Path(__file__).resolve().parent)
    nh=1-2*np.array(list('00000100110101001110'),int)
    tab=sim.filtered_code_tables(code,nh,4092000)
    bits=np.ones(20,dtype=int)
    t=.04+np.arange(4092)/4092000
    a=sim.filtered_spread(t+.10/sim.RATE,bits,tab)
    b=sim.filtered_spread(t+.40/sim.RATE,bits,tab)
    rms=np.sqrt(np.mean((a-b)**2))
    assert rms>.05, 'Sub-chip delays are still indistinguishable'
    worst=0.
    for p in (-1,1):
        for c in (-1,1):
            for n in (-1,1):
                for after in (-1,1):
                    left=2*int(p*c<0)+int(n*c<0)
                    right=2*int(c*n<0)+int(after*n<0)
                    worst=max(worst,abs(c*tab[left,-1]-n*tab[right,0]))
    assert worst<1e-5, f'Symbol boundary discontinuity {worst}'
    y=tab[0,:-1]; spectrum=np.fft.rfft(y)
    freq=np.fft.rfftfreq(len(y),1/(8*sim.RATE))
    leak=np.sum(abs(spectrum[freq>0.90*4092000/2+1])**2)/np.sum(abs(spectrum)**2)
    assert leak<1e-10, f'Unexpected out-of-band energy {leak}'
    print(f'PASS: delay waveform RMS difference={rms:.6f}, boundary error={worst:.3g}, out-of-band fraction={leak:.3g}')

if __name__=='__main__': run()
