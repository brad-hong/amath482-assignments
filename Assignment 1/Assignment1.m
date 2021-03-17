% Clean workspace
clear all; close all; clc

load subdata.mat % Imports the data as the 262144x49 (space by time) matrix called subdata

L = 10; % spatial domain
n = 64; % Fourier modes
x2 = linspace(-L,L,n+1); x = x2(1:n); y =x; z = x;
k = (2*pi/(2*L))*[0:(n/2 - 1) -n/2:-1]; ks = fftshift(k);

[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

%plot(x, X, 'Linewidth', 2)

%ut = fft(Kx);
%subplot(2,1,1)
%plot(fftshift(X), fftshift(abs(ut)), 'r', 'Linewidth', 2)
avg = zeros(n,n,n);
for j=1:49
Un(:,:,:)=reshape(subdata(:,j),n,n,n);
%M = max(abs(Un),[],'all');
ut = fft(Un);
%Mt = max(abs(ut),[],'all');
avg = avg + ut;
%isosurface(Kx,Ky,Kz, fftshift(abs(ut)/Mt,1), 0.7)
%close all, isosurface(X,Y,Z,abs(Un)/M,0.7)
%axis([-20 20 -20 20 -20 20]), grid on, drawnow
%pause(1)
end
avg = fftshift(abs(avg), 1)/49;
[val,idx] = max(avg(:));
dims = size(avg);
smax = cell(size(dims));
[smax{:}] = ind2sub(dims,idx);
smax
isosurface(Kx,Ky,Kz, avg, 0.7)
xlabel('frequency (kx)') % -8 3
ylabel('frequency (ky)') % -7
zlabel('frequency (kz)') % -5 5
%axis([-12 12 -12 12 -12 12]), grid on, drawnow