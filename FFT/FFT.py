import math, random, cmath, timeit
pi = math.pi
def e(n):
    return cmath.exp(complex(0, 2*pi*n))


def fast_pow(a, n):
    n_bin = bin(n)
    result = 1
    k = 2  # top
    while True:
        result *= result
        if n_bin[k] == '1':
            result*= a
        k = k + 1
        if k >= len(n_bin):
            break
    return result

def DFT(x, N, T=1):
    Xs = []
    for k in range(N):
        terms = [x(n*T) * cmath.exp(complex(0, -2*pi*n*k/N)) * T
                 for n in range(N)]
        Xs.append(sum(terms))
    return Xs

def DFT2(samples):
    Xs = []
    N = len(samples)
    for k in range(N):
        terms = [samples[n] * e(-n*k/N) for n in range(N)]
        Xs.append(sum(terms))
    return Xs

def IDFT(X, N, T=1):
    xs = []
    for n in range(N):
        terms = [X(2*pi*k/(N*T)) * cmath.exp(complex(0, 2*pi*n*k/N)) / T
                 for k in range(N)]
        xs.append(sum(terms) / N)
    return xs

def IDFT2(samples):
    xs = []
    N = len(samples)
    for n in range(N):
        terms = [samples[k] * e(n*k/N) for k in range(N)]
        xs.append(sum(terms) / N)
    return xs

def merge(s1, s2):
    merged = s2.copy()
    for i in range(len(s2)):
        merged.insert(2*i, s1[i])
    return merged

def step_merge(S, s, step):
    loc = step
    for e in s:
        S.insert(loc, e)
        loc += step+1
    return S

def FFT(samples):
    N = len(samples)
    if N <= 2:
        return DFT2(samples)
    else:
        evens = FFT([samples[n] + samples[n+N//2]
                      for n in range(N//2)])
        odds  = FFT([(samples[n] - samples[n+N//2]) * e(-n/N)
                      for n in range(N//2)])
        return merge(evens, odds)

def FFT2(samples, Q=2):
    N = len(samples)
    if N <= Q:
        return DFT2(samples)
    else:
        P = N // Q
        Fs = []
        for r in range(Q):
            f1_r_samples = [e(-p*r/(P*Q)) * sum(samples[q*P+p]*e(-r*q/Q)
                                             for q in range(Q))
                          for p in range(P)]
            F_r = FFT2(f1_r_samples)
            Fs = step_merge(Fs, F_r, r)
        return Fs

def IFFT(samples):
    N = len(samples)
    samples_conj = [sample.conjugate() for sample in samples]
    f_conjs = FFT(samples_conj)
    return [f_conj.conjugate() / N for f_conj in f_conjs]
    


if __name__ == "__main__":
    from matplotlib import pyplot as plt
    def half_exp(t):
        return math.exp(-0.1*t) if 0<=t<=100  else 0
    xs1 = [x for x in range(-10,110)]
    ys1 = [half_exp(t) for t in xs1]
    plt.plot(xs1, ys1)
    plt.title("Exponential")
    plt.show()

    def transformed(w):
        return 1/complex(0.1, w)
    xs2 = [x/10 for x in range(100)]
    ys2 = [transformed(w) for w in xs2]
    ys_real2 = [c.real for c in ys2]
    ys_imag2 = [c.imag for c in ys2]
    plt.plot(xs2, ys_real2, label="real")
    plt.plot(xs2, ys_imag2, label="imag")
    plt.legend()
    plt.title("Actual Fourie transformation of Exponential")
    plt.show()

    N, T = 1024, 1
    X = DFT(half_exp, N, T)
    xs3 = [2*pi*k/(N*T) for k in range(N)]
    ys_real3 = [c.real for c in X]
    ys_imag3 = [c.imag for c in X]
    plt.plot(xs3, ys_real3, label="real")
    plt.plot(xs3, ys_imag3, label="imag")
    plt.legend()
    plt.title("DFT of Exponential(N = 1024)")
    plt.show()

    user_input = input("More example? Input 'y' or 'yes' if you want\n")
    if user_input == 'y' or user_input == 'yes':
        def gauss(t):
            return math.exp(-0.1*((t-50)**2)) if 0<=t<=100  else 0
        ys4 = [gauss(t) for t in xs1]
        plt.plot(xs1, ys4)
        plt.title("Gaussian")
        plt.show()

        print("calculation time of DFT:")
        print(timeit.timeit("DFT(gauss, N, T)",number=1,globals=globals()))
        X2 = DFT(gauss, N, T)
        ys_real5 = [c.real for c in X2]
        ys_imag5 = [c.imag for c in X2]
        plt.plot(xs3, ys_real5, label="real")
        plt.plot(xs3, ys_imag5, label="imag")
        plt.legend()
        plt.title("DFT of Gaussian(N = 1024)")
        plt.show()

        sin_samples = [math.sin(n) for n in range(N)]
        print("calculation time of FFT:")
        print(timeit.timeit("FFT(sin_samples)",number=1,globals=globals()))
        X3 = FFT(sin_samples)
        ys_real6 = [c.real for c in X3]
        ys_imag6 = [c.imag for c in X3]
        plt.plot(xs3, ys_real6, label="real")
        plt.plot(xs3, ys_imag6, label="imag")
        plt.legend()
        plt.title("FFT of sin(x) (N = 1024)")
        plt.show()
    
