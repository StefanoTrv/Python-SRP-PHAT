import math
import numpy as np
import soundfile as sf
from scipy.fft import fft, ifft, fftshift

def main():
    # FFT data size
    L = 2048
    # overlap
    O = L
    # speed of sound
    V = 340
    # number of microphones
    M = 8
    # inter-microphone distance ULA
    D = 0.2
    # room setup
    room = [4, 3]
    # source position
    src_pos = [1, 1]
    
    hanning_window = np.hanning(L)

    wh = np.tile(hanning_window.reshape(L, 1), (1, M))
    
    P = int(M*(M-1)/2)

    mic_pos = np.zeros((M, 2))
    mic_pos[0, :] = [1.1, 0.2]
    for i in range(1,M):
        mic_pos[i, :] = [mic_pos[i-1,0]+D, 0.2];

    xsig, fs = sf.read('SRP_PHAT_8ch_rn.wav')
    W = math.floor((len(xsig) - L) / O) if (len(xsig) - L) / O >= 0 else math.ceil((len(xsig) - L) / O)

    RES = 0.01
    DX = len(np.arange(0.5, 3.5 + RES, RES))
    DY = len(np.arange(0.5, 2.5 + RES, RES))
    NG = DX * DY * P
    POS_TDOA = np.zeros((int(NG), 4))
    DSM = np.zeros((M, 1))
    idx = 0
    ix = 1

    for i in np.arange(0.5, 3.5 + RES, RES):
        iy = 1
        for j in np.arange(0.5, 2.5 + RES, RES):
            for z in range(0,M):
                DSM[z] = math.sqrt((mic_pos[z, 0] - i) ** 2 + (mic_pos[z, 1] - j) ** 2)
            g = 1
            for ii in range(0,M-1):
                for jj in range(ii+1,M):
                    TDOA = round((DSM[ii,0] - DSM[jj,0]) * fs / V)
                    POS_TDOA[idx, 0] = ix
                    POS_TDOA[idx, 1] = iy
                    POS_TDOA[idx, 2] = g
                    POS_TDOA[idx, 3] = TDOA
                    idx += 1
                    g += 1
            iy += 1
        ix += 1
    
    ZERO = int(L/2+1)

    GCC = np.zeros((P, L))
    k = 0
    for i in range(0,W):
        X = fft(wh * xsig[k:k+L, :], axis=0)
        g = 0
        for ii in range(0,M-1):
            for jj in range(ii+1,M):
                CS = X[:,ii] * np.conj(X[:, jj])
                RP = CS / np.abs(CS)
                GCC[g,:] = np.real(fftshift(ifft(RP)))
                g += 1
        map = np.zeros((DY,DX))
        for j in range(0,NG):
            map[int(POS_TDOA[j, 1])-1, int(POS_TDOA[j, 0])-1] = map[int(POS_TDOA[j, 1])-1, int(POS_TDOA[j, 0])-1] + GCC[int(POS_TDOA[j,2])-1,ZERO+int(POS_TDOA[j,3])-1]

        mv = np.max(map)
        ys, xs = np.where(map == mv)
        xs = xs*RES+0.5-RES
        ys = ys*RES+0.5-RES
        print((xs, ys))

        k = k+O






if __name__ == "__main__":
    main()