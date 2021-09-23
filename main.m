clear all; close all; clc;
%[xClean, Fs] = audioread('melodiak5_part.wav');
[xNoisy, Fs] = audioread('noisyMelody.wav');

[x1, Fs] = audioread('la4k5.wav');
[x2, Fs] = audioread('fa4k5.wav');
[x3, Fs] = audioread('do5k5.wav');
[x4, Fs] = audioread('sol4k5.wav');
[x5, Fs] = audioread('re#5k5.wav');
[x6, Fs] = audioread('si4k5.wav');
[x7, Fs] = audioread('la#4k5.wav');
[x8, Fs] = audioread('re5k5.wav');

[xNoise, Fs2] = audioread('mainsBrum50Hz.wav');
%xNoise = resample(xNoise,2,1);


segmentTime = 0.035;
segmentLength = round(segmentTime*Fs);
specWindow = hann(segmentLength, 'periodic');
% specWindow = ones(segmentLength,1);
nDft = 4096;
nOverlap = 0.5*segmentLength;
maxDynRange = 60;
lambda = 0.01;

%%
%[sparseS, sparseF, sparseT] = sparseStft(xClean, specWindow, nOverlap, nDft, Fs, lambda);
[sparseSNoisy, sparseFNoisy, sparseTNoisy] = sparseStft(xNoisy, specWindow, nOverlap, nDft, Fs, lambda);

[sparseS1, sparseF1, sparseT1] = sparseStft(x1, specWindow, nOverlap, nDft, Fs, lambda);
[sparseS2, sparseF2, sparseT2] = sparseStft(x2, specWindow, nOverlap, nDft, Fs, lambda);
[sparseS3, sparseF3, sparseT3] = sparseStft(x3, specWindow, nOverlap, nDft, Fs, lambda);
[sparseS4, sparseF4, sparseT4] = sparseStft(x4, specWindow, nOverlap, nDft, Fs, lambda);
[sparseS5, sparseF5, sparseT5] = sparseStft(x5, specWindow, nOverlap, nDft, Fs, lambda);
[sparseS6, sparseF6, sparseT6] = sparseStft(x6, specWindow, nOverlap, nDft, Fs, lambda);
[sparseS7, sparseF7, sparseT7] = sparseStft(x7, specWindow, nOverlap, nDft, Fs, lambda);
[sparseS8, sparseF8, sparseT8] = sparseStft(x8, specWindow, nOverlap, nDft, Fs, lambda);

[sparseSNoise, sparseFNoise, sparseTNoise] = sparseStft(xNoise, specWindow, nOverlap, nDft, Fs2, lambda);



%%
[D1,A1, cost1] = nnmf(abs(sparseS1),1);
[D2,A2, cost2] = nnmf(abs(sparseS2),1);
[D3,A3, cost3] = nnmf(abs(sparseS3),1);
[D4,A4, cost4] = nnmf(abs(sparseS4),1);
[D5,A5, cost5] = nnmf(abs(sparseS5),1);
[D6,A6, cost6] = nnmf(abs(sparseS6),1);
[D7,A7, cost7] = nnmf(abs(sparseS7),1);
[D8,A8, cost8] = nnmf(abs(sparseS8),1);

[DNoise,ANoise, costNoise] = nnmf(abs(sparseSNoise),1);

D = [D1 D2 D3 D4 D5 D6 D7 D8 DNoise];
N = length(sparseTNoisy);
K = 9;
A = rand(K,N);

%%


n_iter = 15;
beta = 0;

[D, A, cost] = beta_nmf_mu(abs(sparseSNoisy).^2, n_iter, D, A, beta);

figure;
subplot(3,3,1);
plot(abs(A(1,:))); title('A1 = Note C5');
subplot(3,3,2);
plot(abs(A(2,:))); title('A1 = Note F4');
subplot(3,3,3);
plot(abs(A(3,:))); title('A1 = Note A#4');
subplot(3,3,4);
plot(abs(A(4,:))); title('A1 = Note B4');
subplot(3,3,5);
plot(abs(A(5,:))); title('A1 = Note D#5');
subplot(3,3,6);
plot(abs(A(6,:))); title('A1 = Note D5');
subplot(3,3,7);
plot(abs(A(7,:))); title('A1 = Note B4');
subplot(3,3,8);
plot(abs(A(8,:))); title('A1 = Note G4');
subplot(3,3,9);
plot(abs(A(9,:))); title('A9 = Noise');

K1 = ((D1*A(1,:)).*(D*A).^(-1)).*sparseSNoisy;
K2 = ((D2*A(2,:)).*(D*A).^(-1)).*sparseSNoisy;
K3 = ((D3*A(3,:)).*(D*A).^(-1)).*sparseSNoisy;
K4 = ((D4*A(4,:)).*(D*A).^(-1)).*sparseSNoisy;
K5 = ((D5*A(5,:)).*(D*A).^(-1)).*sparseSNoisy;
K6 = ((D6*A(6,:)).*(D*A).^(-1)).*sparseSNoisy;
K7 = ((D7*A(7,:)).*(D*A).^(-1)).*sparseSNoisy;
K8 = ((D8*A(8,:)).*(D*A).^(-1)).*sparseSNoisy;

KCombined = K1 + K2 + K3 + K4 + K5 + K6 + K7 + K8;

iy = istft(KCombined, Fs, 'Window', specWindow, 'OverlapLength', nOverlap, 'FFTLength', nDft);
iy = iy(10:1000000);
%plot(abs(iy));
out = abs(iy);

soundsc(out,Fs);

% -20dB for audiowrite
R_dB = -20;

R = 10^(R_dB/20);

a = sqrt((length(out)*R^2)/(sum(out.^2)));

out = out*a;
audiowrite('newFile.wav',out, Fs);

%%
figure(1)
imagesc(sparseTNoisy,sparseFNoisy,...
    10*log10(dynamicRangeLimiting(abs(sparseSNoisy).^2, maxDynRange)))
set(gca,'YDir','normal')
xlabel('time [s]')
ylabel('frequency [Hz]')
ylim([0 1.2*10^4]);
colorbar

figure(2)
imagesc(sparseTNoisy,sparseFNoisy,...
    10*log10(dynamicRangeLimiting(abs(KCombined).^2, maxDynRange)))
set(gca,'YDir','normal')
xlabel('time [s]')
ylabel('frequency [Hz]')
ylim([0 1.2*10^4]);
colorbar