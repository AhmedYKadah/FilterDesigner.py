import numpy as np
from IIR import _pz2iir
"""
Implement the Direct-Form Notch Filter
"""

def _NotchDF(fc,fs,B):

    z=np.exp(fc/fs*2*np.pi*1j)

    b,a=_pz2iir([z*(1-B/fs),np.conj(z)*(1-B/fs)],
            [z,np.conj(z)])
    return b,a

    pass

