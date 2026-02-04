from SanathananKoerner import _SanathananKoerner
from LS import _LS, _WLS
from ParksMcClellan import _parksMcClellan
import numpy as np

"""
Filter Designer 
"""
def FilterDesigner(N, Fs, Fpass, Fstop, dp, ds, fType="lowpass", algorithm="LS"):
    """
    N is the order of the filter 
    Fs is the sampling frequency
    Fpass is the passband frequency
    Fstop is the stopband frequency
    dp is the pass band ripple
    ds is the stop band attenuation
    ftype is the filter type
    algorithm is the filter design method
    """
    pass

def _DesiredResponseGenerator(N, Fs, ds, Fpass, Fstop, BW=0,Fpass2=0,Fstop2=0, fType="lowpass"):
    """
    N is the number of samples 
    """

    Fs=Fs/2 # Gives the maximum frequency measurable given the sampling frequency
    attenuation=np.power(10.0,-ds)

    if (Fpass>Fs)|( Fstop>Fs )|(Fpass2>Fs)|(Fstop2>Fs):
        raise ValueError('Passband frequency must be lower than half the sampling frequency.')
        pass
    if (Fpass<0)|( Fstop<0 )|(Fpass2<0)|(Fstop2<0)|(BW<0):
        raise ValueError('Frequencies must be positive.')
        pass
    if fType.lower()=='lowpass':
        if Fpass>Fstop:
            raise ValueError('Passband frequency must be lower than stop band frequency for lowpass filters.')
        passBand=int(N*Fpass/Fs)
        transBand=int(N*(Fstop-Fpass)/Fs)
        stopBand=N-passBand-transBand
        return np.concatenate((np.ones(passBand),np.linspace(1,attenuation,transBand),attenuation*np.ones(stopBand)))
        pass
    elif fType.lower()=='highpass':
        if Fpass<Fstop:
            raise ValueError('Passband frequency must be higher than stop band frequency for highpass filters.')
        stopBand=int(N*Fstop/Fs)
        transBand=int(N*(Fpass-Fstop)/Fs)
        passBand=N-stopBand-transBand
        return np.concatenate((attenuation*np.ones(stopBand),np.linspace(attenuation,1,transBand),np.ones(passBand)))
        pass
    elif fType.lower()=='bandpass':
        if (BW==0) & (Fpass2==0 or Fstop2==0):
            raise ValueError('No parameters given for second cut off frequencies. Give a bandwidth value or add second cut offs')
        if BW==0:
            if (Fpass<Fstop)|(Fpass2>Fstop2):
                raise ValueError('Passband frequency must be higher than stop band frequency highpass filters.')
            stopBand1=int(N*Fstop/Fs)
            transBand1=int(N*(Fpass-Fstop)/Fs)
            passBand=int(N*(Fpass2-Fpass)/Fs)
            transBand2=int(N*(Fstop2-Fpass2)/Fs)
            stopBand2=min(N-transBand2-stopBand1-transBand1-passBand,0)
            return np.concatenate((attenuation*np.ones(stopBand1),np.linspace(attenuation,1,transBand1),np.ones(passBand),
                                   np.linspace(1,attenuation,transBand2),attenuation*np.ones(stopBand2)))

        else:
            if (Fstop<(Fpass+BW/2)) & (Fstop>(Fpass-BW/2)):
                raise ValueError('Fstop must be outside the passband')
            if BW>=2*Fs:
                return np.ones(N)

            trans=max(Fstop-(Fpass+BW/2),(Fpass-BW/2)-Fstop)
            slope=(1-attenuation)/trans

            stopBand1=max(int(N*(Fpass-BW/2-trans)/Fs) ,0)

            transBand1=int(N*(trans)/Fs) if (Fpass-BW/2-trans)>=0 else max(int(N*(Fpass-BW/2)/Fs),0)

            passBand=min(int(N*(BW-max(BW/2-Fpass,0))/Fs), N - transBand1 - stopBand1)

            transBand2=min(int(N*(trans)/Fs), N-stopBand1-transBand1-passBand)
            stopBand2=min(N-transBand2-stopBand1-transBand1-passBand,0)

            sp1= max(1-(Fpass-BW/2)*slope,attenuation)
            sp2= max(1-(Fs-(Fpass+BW/2))*slope,attenuation)

            return np.concatenate((attenuation*np.ones(stopBand1),np.linspace(sp1,1,transBand1),np.ones(passBand),
                                   np.linspace(1,sp2,transBand2),attenuation*np.ones(stopBand2)))
        pass
    elif fType.lower()=='bandstop':
        return np.ones(N)-_DesiredResponseGenerator(N, Fs, ds, Fstop, Fpass, BW ,Fstop2,Fpass2, "bandpass")
        pass
    else:
        raise ValueError('fType: Invalid filter type.')
        pass
    pass
