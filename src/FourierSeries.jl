__precompile__(true)
module FourierSeries

    export
        fourierSeriesStepReal,
        fourierSeriesSampledReal,
        fourierSeriesSynthesisReal,
        repeatPeriodically

    """
    # Function call

    `fourierSeriesStepReal(t,u,T,hMax)`

    # Description

    Calculates the real Fourier coefficients of a step function based on the
    data pairs `t` and `u` according to the following variables

    ```math
    \mathtt{a[k]} =
    \frac{2}{\mathtt{T}} \int_{\mathtt{t[1]}}^{\mathtt{t[1]+T}}
    \mathtt{u(t)} \cos(k\omega \mathtt{t}) \operatorname{d}\mathtt{t} \\
    \mathtt{b[k]} =
    \frac{2}{\mathtt{T}} \int_{\mathtt{t[1]}}^{\mathtt{t[1]+T}}
    \mathtt{u(t)} \sin(k\omega \mathtt{t}) \operatorname{d}\mathtt{t}
    ```

    # Function arguments

    `t` Time vector instances where a step occurs

    `u` Data vector of which Fourier coefficients shall be determined of, i.e.,
    `u[k]` is the data value from the rigth at `t[k]`

    `T` Time period of `u`

    `hMax` Maximum harmonic number to be determined

    # Return values

    `h` Vector of harmonic numbers start at, i.e., `h[1] = 0` (dc value of `u`),
    `h[hMax+1] = hMax`

    `f` Frequency vector associated with harmonic numbers, i.e. `f = h / T`

    `a` Cosine coefficients of Fourier series, where `a[1]` = dc value of `u`

    `b` Sine coefficients of Fourier series, where `b[1] = 0`

    # Examples

    Copy and paste code:

    ```julia
    # Fourier approximation of square wave form
    ts = [0;1]                  # Sampled time data
    us = [-1,1]                 # Step function of square wave
    T = 2                       # Period
    hMax = 7                    # Maximum harmonic number
    (h,f,a,b) = fourierSeriesStepReal(ts,us,2,hMax)
    (t,u)=fourierSeriesSynthesisReal(f,a,b)
    using PyPlot
    # Extend (t,u) by one element to right periodically
    (tx,ux) = repeatPeriodically(ts,us,right=1)
    step(tx,ux,where="post")    # Square wave function
    plot(t,u)                   # Fourier approximation up to hMax = 7
    ```
    """
    function fourierSeriesStepReal(t,u,T,hMax)
        # Check if t and have equal lengths
        if length(t)!=length(u)
            error("module FourierSeries: function fourierSeriesStepReal:\n
    Vectors t and u have different lengths")
        end
        # Check if vector u is NOT of type Complex
        if u[1] isa Complex
            error("module FourierSeries: function fourierSeriesStepReal:\n
    The analyzed functions `u` must not be of Type ::Complex")
        end
        # Determine length of function u
        N = length(u)
        # Initialization of real result vectors
        a=zeros(hMax+1)
        b=zeros(hMax+1)
        # Cycle through loop to determine coefficients
        i=collect(1:N)
        # Time vector, extended by period T
        τ=cat(1,t,[T+t[1]])
        # DC value
        a[1]=sum(u.*(τ[i+1]-τ[i]))/T
        b[1]=0
        for k in collect(1:hMax)
            a[k+1]=sum(u.*(+sin.(k*τ[i+1]*2*pi/T)
                           -sin.(k*τ[i]*2*pi/T)))/(k*pi)
            b[k+1]=sum(u.*(-cos.(k*τ[i+1]*2*pi/T)
                           +cos.(k*τ[i]*2*pi/T)))/(k*pi)
        end
        # Number of harmonics
        h = collect(0:hMax)
        # Frequencies
        f = h/T
        return (h,f,a,b)
    end

    function fourierSeriesSampledReal(t,u,hMax::Int64=typemax(Int64))
        # Check if vector u is NOT of type Complex
        if u[1] isa Complex
            error("module FourierSeries: function fourierSeriesSampledReal:\n
    The analyzed functions `u` must not be of Type ::Complex")
        end
        # Determine length of function u
        N = length(u)
        # Determine period of function u
        T = t[end]-t[1]+t[2]-t[1]
        # Apply fast Fourier transform
        c = 2*fft(u)/N
        # Correct DC value by factor 1/2
        c[1]=c[1]/2;
        if mod(N,2)==1
            # Number of sampled data, N, is odd
            hMaxMax=div(N-1,2)
        else
            # Number of sampled data, N, is even
            hMaxMax=div(N,2)
        end
        # Number of harmonics
        h = collect(0:hMaxMax)
        # Frequencies
        f = h/T
        # Real and imaginary coefficients
        a = +real(c[1:min(hMax,hMaxMax)+1])
        b = -imag(c[1:min(hMax,hMaxMax)+1])
        # Return vectors
        return (h,f,a,b)
    end

    function fourierSeriesSynthesisReal(f,a,b;hMax=length(a)-1,N=1000,t0=0)
        # Check if vectors f, a and b have equal lengths
        if length(f)!=length(a)
            error("module Fourier: function fourierSeriesSynthesisReal:\n
    Vectors f and a have different lengths")
        end
        if length(a)!=length(b)
            error("module Fourier: function fourierSeriesSynthesisReal:\n
    Vectors a and b have different lengths")
        end
        # Determine Period of synthesized function
        T = 1/f[2]
        # Initialization of synthesis function f
        u=fill(a[1],N)
        # time samples of synthesis function, considering start time t0
        t=collect(0:N-1)/N*T+fill(t0,N)
        # hMax may either be a scalar of vector
        if hMax isa Array
            # If hMax is an array, then synthesize only harmonic numbers
            # indicated by hMax
            kRange = hMax
        else
            # If hMax is a scalar, then treat hMax as the maximum harmonic
            # number, which may not exceed length(a)-1, as a[1] equals the dc
            # component (harmonic 0) and a[hMax-1] represents harmonic number
            # hMax
            kRange = collect(1:min(hMax,length(a)-1))
        end
        # Calculate superosition
        for k in kRange
            u = u + a[k+1]*cos.(k*t*2*pi/T) + b[k+1]*sin.(k*t*2*pi/T)
        end
        return (t,u)
    end

    function repeatPeriodically(t,u;left=0,right=1)
        # Check if vectors t and u have equal lengths
        if length(t)!=length(u)
            error("module Fourier: function repeatPeriodically:\n
    Vectors t and u have different lengths")
        end
        # Determine length of arrays t and u
        N = length(t)
        # Determine period of time axis
        T = t[end]-t[1]+t[2]-t[1]
        # Determine how many times the vectors t and u shall be repeated
        outerRight = Int(floor((right-1)/N)+1)
        outerLeft  = Int(floor((left -1)/N)+1)
        # Repeat vector t by full periods and consider period T
        tx = repeat(t,outer=outerRight+outerLeft+1) +
             repeat(collect(-outerLeft:outerRight),inner=N)*T
        # Repeat vector u by full periods
        ux = repeat(u,outer=outerLeft+outerRight+1)
        # Limit output vectors elements
        tx = tx[outerLeft*N-left+1:(outerLeft*N)+N+right]
        ux = ux[outerLeft*N-left+1:(outerLeft*N)+N+right]
        return (tx,ux)
    end

end
