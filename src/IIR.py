import numpy as np

"""
1- Implement a pole-zero placement function 
"""
def _pz2iir(p=[],z=[],g=1):
    b=np.polynomial.polynomial.polyfromroots(z)
    a=np.polynomial.polynomial.polyfromroots(p)
    order=max(len(a),len(b))-1
    if order!=0:
        b = g*np.pad(b, (order-len(b), 0), mode='constant', constant_values=0)[::-1]
        a = np.pad(a, (order-len(a), 0), mode='constant', constant_values=0)[::-1]

        b=b*a[0]
        a=a[1:]/a[0]
    return b.real ,a.real
    pass

"""
2- Implement a stability checker
"""

def _checkStability(a):
    poles=np.roots(([1, *a]))
    return np.linalg.norm(poles,np.inf)>1.0
    pass

