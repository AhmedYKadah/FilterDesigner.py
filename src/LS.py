import numpy as np
"""
1- Implementation of Least Squares
"""
def _LS(order, Hd, weights=0):
    if weights!=0:
        return _WLS(order, Hd,n,weights)
    if order%2==0:
        order+=1
    L=(order+1)/2 # Number of unique coefficients on one side of a symmetric filter
    samples=len(Hd) # Number of samples

    n=np.arange(0,L,1) # Iterator over h[n] values in the DTFT summation
    ws=np.linspace(0,np.pi,samples+1)[:samples] # Evenly spaced sample points in the frequency domain

    F=2*np.cos(ws.reshape(-1,1) @ n.reshape(1,-1)) # Generating the Fourier matrix
    F[:,0]=F[:,0]/2 # Setting first column to all ones (h[0] is not repeated)
    FT=np.transpose(F) 

    h= (np.linalg.inv(FT @F) @ FT )@ Hd.reshape(-1,1) # Least-Squares exact solution

    return np.concatenate((h[::-1][:-1],h))

    pass

"""
2- Implementation of Weighted Least Squares ## UNTESTED
"""
def _WLS(order, Hd, weights):
    if order%2==0:
        order+=1
    L=(order+1)/2 # Number of unique coefficients on one side of a symmetric filter
    samples=len(Hd) # Number of samples

    n=np.arange(0,L,1) # Iterator over h[n] values in the DTFT summation
    ws=np.linspace(0,np.pi,samples+1)[:samples] # Evenly spaced sample points in the frequency domain

    F=2*np.cos(ws.reshape(-1,1) @ n.reshape(1,-1)) # Generating the Fourier matrix
    F[:,0]=F[:,0]/2 # Setting first column to all ones (h[0] is not repeated)
    FT=np.transpose(F) 

    h= (np.linalg.inv(FT @ np.diag(weights) @ F) @ FT ) @ np.diag(weights) @ Hd.reshape(-1,1) # Least-Squares exact solution

    return np.concatenate((h[::-1][:-1],h))

    pass
    pass
