module FourierSeries

    export fourierSeriesStep,
        fourierSeriesSampled,
        fourierSeriesSynthesis,
        fourierComplexToReal,
        fourierRealToComplex

    function fourierSeriesStep(t,u,hMax)
        # Check if t and have equal lengths
        if length(t)!=length(u)
            error("module FourierSeries: function fourierSeriesStep:\n
    Vectors t and u have different lengths")
        end
        # Check if vector u is NOT of type Complex
        if u[1] isa Complex
            error("module FourierSeries: function fourierSeriesStep:\n
    The analyzed functions `u` must not be of Type ::Complex")
        end
        # Determine length of function u
        N=length(u)
        # Determine period of function u
        T = t[end]-t[1]+t[2]-t[1]
        # Initialization of complex result vectors
        c=1im*zeros(hMax+1)
        # Cycle through loop to determine coefficients
        i=collect(1:N)
        c[1]=sum(u)/(2*pi)
        for k in collect(1:hMax)
            c[k+1]=sum(u.*(-exp.(-1im*k*i*2*pi/N)+exp.(-1im*k*(i-1)*2*pi/N)))/(1im*k*pi)
        end
        # Number of harmonics
        h = collect(0:hMax)
        # Frequencies
        f = h/T
        return (h,f,c)
    end

    function fourierSeriesSampled(t,u,hMax::Int64=typemax(Int64))
        # Check if vector u is NOT of type Complex
        if u[1] isa Complex
            error("module FourierSeries: function fourierSeriesSampled:\n
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
        # Limit output
        return (h,f,c[1:min(hMax,hMaxMax)+1]);
    end

    function fourierSeriesSynthesis(f,c;hMax=length(c)-1,N=1000)
        # Check if vectors f and c have equal lengths
        if length(f)!=length(c)
            error("module Fourier: function fourierSeriesSynthesis:\n
    Vectors f and c have different lengths")
        end
        # Determine Period of synthesized function
        T = 1/f[2]
        # Initialization of synthesis function f
        u=fill(real(c[1]),N)
        # Indices of synthesis function
        t=collect(0:N-1)/N*T
        # hMax may either be a scalar of vector
        if hMax isa Array
            # If hMax is an array, then synthesize only harmonic numbers indicated by
            # hMax
            kRange = hMax
        else
            # If hMax is a scalar, then treat hMax as the maximum harmonic number,
            # which may not exceed length(c)-1, as c[1] equals the dc component
            # (harmonic 0) and c[hMax-1] represents h number hMax
            kRange = collect(1:min(hMax,length(c)-1))
        end
        # Calculate superposition
        for k in kRange
            u = u + real(c[k+1]*exp.(1im*k*t*2*pi/T))
        end
        return (t,u)
    end

    function fourierComplexToReal(c)
        return (real(c),-imag(c))
    end

    function fourierRealToComplex(a,b)
        return a-1im*b
    end
end
