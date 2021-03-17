%% 
%figure(1)
[y, Fs] = audioread('Floyd.m4a');
y = y(1:end-1,1);
p = 0;
len = length(y);
y1 = y((len/8)*p+1:(len/8)*(p+1));
tr_gnr = length(y)/Fs; % record time in seconds
tr_gnr1 = tr_gnr/8;
%plot((1:length(y))/Fs,y);
%xlabel('Time [sec]'); ylabel('Amplitude');
%title('Comfortably Numb');
%p8 = audioplayer(y,Fs); playblocking(p8);
% t2 = linspace((tr_gnr1)*p,(tr_gnr1)*(p+1), (len+1)/8);
% t = t2(1:(len+1)/8);
t2 = linspace(0,tr_gnr, (len+1));
t = t2(1:len);
%tau = 3; % centre of window
%a = 3; % window size
%g = exp(-a*(t-tau).^2);
%Sf = g.*y';
%Sft = fft(Sf);
% k = (1/(tr_gnr1))*[0:(len/8)/2-1 -(len/8)/2:-1];
k = (1/tr_gnr)*[0:len/2-1 -len/2:-1];
ks = fftshift(k);
%figure(2)
%subplot(2,1,1) % Time domain
%plot(t,Sf,'k','Linewidth',2); axis([0 15 -.05 .05])
%set(gca,'Fontsize',16), xlabel('Time (t)'), ylabel('S(t)*g(t-\tau)')
%subplot(2,1,2) % Fourier domain 
%plot(ks,abs(fftshift(Sft))/max(abs(Sft)),'r','Linewidth',2); axis([80 1200 0 1])
%set(gca,'Fontsize',16)
%xlabel('frequency, Hz (k)'), ylabel('FFT(S(t)*g(t-\tau))')
%% 

% figure(3)
% a = 0.3;
% % tau = (tr_gnr1)*p:.5:(tr_gnr1)*(p+1);
% tau = 0:2:tr_gnr;
% Sgt_spec = [];
% for j = 1:length(tau)
%     g = exp(-a*(t - tau(j)).^2);
%     Sg = g.*y';
%     Sgt = fft(Sg);
%     Sgt_spec(:,j) = fftshift(abs(Sgt));
% end
% pcolor(tau,ks,Sgt_spec)
% shading interp
% set(gca,'ylim',[60 150],'Fontsize',16)
% colormap(hot)
% %colorbar
% xlabel('time (t)'), ylabel('frequency (k)')
% title(['a = ',num2str(a)],'Fontsize',16)

%% 
[y, Fs] = audioread('Floyd.m4a');
y = y(1:end-1,1);
p = 0;
len = length(y);
y1 = y((len/8)*p+1:(len/8)*(p+1));
tr_gnr = length(y)/Fs; % record time in seconds
tr_gnr1 = tr_gnr/8;
t2 = linspace((tr_gnr1)*p,(tr_gnr1)*(p+1), (len+1)/8);
t = t2(1:(len+1)/8);
% t2 = linspace(0,tr_gnr, (len+1));
% t = t2(1:len);
k = (1/(tr_gnr1))*[0:(len/8)/2-1 -(len/8)/2:-1];
% k = (1/tr_gnr)*[0:len/2-1 -len/2:-1];
ks = fftshift(k);

figure(5)
a = 0.3;
tau = (tr_gnr1)*p:.5:(tr_gnr1)*(p+1);
% tau = 0:1:tr_gnr;
Sgt_spec = [];
for j = 1:length(tau)
    g = exp(-a*(t - tau(j)).^2);
    Sg = g.*y1';
    Sgt = fft(Sg);
    filtered = bandpass(Sgt, [250 800], Fs);
    Sgt_spec(:,j) = fftshift(abs(filtered));
end
pcolor(tau,ks,Sgt_spec)
shading interp
set(gca,'ylim',[250 800],'Fontsize',16)
colormap(hot)
%colorbar
xlabel('time (t)'), ylabel('frequency (k)')
title(['a = ',num2str(a)],'Fontsize',16)
