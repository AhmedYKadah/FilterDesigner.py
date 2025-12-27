import numpy as np
from IIR import _pz2iir
"""
Implement the Direct-Form Notch Filter
"""

def _NotchDF(fc,fs,B):

    z=np.exp(fc/fs*2*np.pi*1j)

    b,a=_pz2iir([zz*(1-B/fs),np.conj(zz)*(1-B/fs)],
            [zz,np.conj(zz)])
    return b,a

    pass

