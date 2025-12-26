# Filter Designer
## Desired Response Generator

# Finite Impulse Response (FIR) Filters 
## Least Squares FIR filter design
We are trying to obtain the filter coefficients of a symmetric filter $h[n]$ for the problem 

$\mathbf{\hat{H}_d}(\omega)=\mathbf{F}(\omega)\mathbf{\hat{h}}[n]$

which is the DTFT of the impulse response, where we are trying to minimize the square error 

$\min (\mathbf{H_d}-\mathbf{\hat{H}_d})^2$


Defining vectors for frequency samples and creating the Fourier matrix corresponding to 

$\mathbf{\Omega}= 
\begin{bmatrix} 
\omega_0 \\ 
\omega_1 \\ 
\vdots \\
\omega_n
\end{bmatrix}
\begin{bmatrix} 
0 &  
1 & 
\dots & 
n
\end{bmatrix}$
$
\mathbf{E_y}=
\begin{bmatrix} 
1 \\
1 \\
\dots \\ 
1
\end{bmatrix}$

to give

$\mathbf{F}=
\begin{bmatrix} 
\mathbf{E_y}|2\cos(\mathbf{\Omega})
\end{bmatrix} 
$

and we solve the least-squares problem exactly using 

$
\mathbf{\hat{h}}=
(\mathbf{F}^T \mathbf{F})^{-1}\mathbf{F}^T \mathbf{H_d}
$

## Weighted Least Squares FIR filter design
We do the same as the ordinary least squares problem but with the addition of a weight function for the error $W$

$\min (W^{1/2}(\mathbf{H_d}-\mathbf{\hat{H}_d}))^2$

$
\mathbf{\hat{h}}=
(\mathbf{F}^T \mathbf{W}\mathbf{F})^{-1}\mathbf{F}^T \mathbf{W}\mathbf{H_d}
$

## Parks-Mcclellan Exchange FIR filter design
### Minimax Polynomials and Approximations 
### The Remez Algorithm
### Chebyshev Polynomials

# Infinite Impulse Response (IIR) Filters 
## Manual Pole-Zero Placement

### Stability check function

## Direct-form Notch filter

## Sanathanan-Koerner Least-Squares IIR filter design

### 
