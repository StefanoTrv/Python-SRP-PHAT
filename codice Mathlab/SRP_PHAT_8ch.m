clear variables
close all
clc

% FFT data size
L = 2048;
% overlap
O = L; 
% speed of sound
V = 340;
% number of microphones
M = 8;
% inter-microhpne distance ULA
D = 0.2;
% room setup
SetupStruc.room = [4 3];  
% source position
SetupStruc.src_pos = [1 1];

wh = repmat(hann(L),1,M);

P = M*(M-1)/2;

SetupStruc.mic_pos = zeros(M,2);
SetupStruc.mic_pos(1,:) = [1.1 0.2];
for i = 2:M
    SetupStruc.mic_pos(i,:) = [SetupStruc.mic_pos(i-1,1)+D 0.2];
end

figure(1);
hold on
for i = 1:M
    plot(SetupStruc.mic_pos(i,1),SetupStruc.mic_pos(i,2),'o','Color',[0 0 0],'MarkerSize',7,'LineWidth',2);
end
plot(SetupStruc.src_pos(1),SetupStruc.src_pos(2),'+','Color',[0 0 0],'MarkerSize',8,'LineWidth',2);
hold off;
axis equal;
xlim([0 SetupStruc.room(1)]);
ylim([0 SetupStruc.room(2)]);
box on;

[xsig, fs] = audioread('SRP_PHAT_8ch_rn.wav');
W = fix((length(xsig)-L)/O);

RES = 0.01;
DX = length(0.5:RES:3.5);
DY = length(0.5:RES:2.5);
NG = DX*DY*P;
POS_TDOA = zeros(NG,4);
DSM = zeros(M,1);
idx = 1;
ix = 1;
for i = 0.5:RES:3.5
    iy = 1;
    for j = 0.5:RES:2.5
        for z = 1:M
            DSM(z) = sqrt((SetupStruc.mic_pos(z,1)-i)^2+(SetupStruc.mic_pos(z,2)-j)^2);
        end
        g = 1;
        for ii = 1:M-1
            for jj = ii+1:M
                TDOA = round((DSM(ii)-DSM(jj))*fs/V);
                POS_TDOA(idx,1) = ix;
                POS_TDOA(idx,2) = iy;
                POS_TDOA(idx,3) = g;
                POS_TDOA(idx,4) = TDOA;
                idx = idx+1;
                g = g+1;
            end
        end
        iy = iy+1;
    end
    ix = ix+1;
end

ZERO = L/2+1;

GCC = zeros(P,L);
k = 1;
for i = 1:W
    X = fft(wh.*xsig(k:k+L-1,:)); 
    g = 1;
    for ii = 1:M-1
        for jj = ii+1:M
            CS = X(:,ii).*conj(X(:,jj));
            RP = CS./abs(CS);
            GCC(g,:) = real(fftshift(ifft(RP)));
            g = g+1;
         end
    end
    map = zeros(DY,DX);
    for j = 1:NG
        map(POS_TDOA(j,2),POS_TDOA(j,1)) = map(POS_TDOA(j,2),POS_TDOA(j,1))+GCC(POS_TDOA(j,3),ZERO+POS_TDOA(j,4));
    end
    figure(3);
    imagesc(map)
    set(gca,'YDir','normal') 
   
    mv = max(map(:));
    [ys,xs]= find(map==mv);
    xs = xs*RES+0.5-RES;
    ys = ys*RES+0.5-RES;
    disp([xs ys]);

    pause;

    k = k+O;
end

