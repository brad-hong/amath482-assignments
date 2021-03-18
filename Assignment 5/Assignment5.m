vid = VideoReader('ski_drop_low.mp4');
dt = 1/vid.FrameRate;
t = 0:dt:vid.Duration; 
nFrames = vid.NumFrames;
height = vid.Height;
width = vid.Width;
X = zeros(height*width, nFrames);

i = 0;
while hasFrame(vid)
   i = i + 1;
   f = readFrame(vid);
   fg = rgb2gray(f);
   sz = size(fg,1) * size(fg,2);
   image = reshape(fg(:,:),sz,1);
   X(:,i) = image;
end
%%
X1 = X(:, 1:end-1);
X2 = X(:, 2:end);
%%
[U,S,V] = svd(X1, 'econ');
numModes = 3;
U_modes = U(:, 1:numModes);
S_modes = S(1:numModes, 1:numModes);
V_modes = V(:, 1:numModes);
%%
S_tilda = U_modes' * X2 * V_modes * diag(1./diag(S_modes));
[eV, D] = eig(S_tilda);
mu = diag(D);
omega = log(mu)/dt;
Phi = U_modes*eV;
%%
y0 = Phi\X1(:,1);
u_modes = zeros(length(y0), length(t));
for iter = 1:length(t)
    u_modes(:, iter) = y0.*exp(omega*t(iter));
end
u_dmd = Phi*u_modes;
%%
for i = 1:nFrames
    image = reshape(u_dmd(:,i),height,width);
    r_image = uint8(real(image));
    imshow(r_image); 
    drawnow
end

%%
u_dmd_fg = X - abs(u_dmd);
for i = 1:nFrames
    image = reshape(u_dmd_fg(:,i),height,width);
    r_image = uint8(real(image));
    imshow(r_image); 
    drawnow
end