[height1, width1, rgb1, num_frames1] = size(vidFrames1_1);
[height2, width2, rgb2, num_frames2] = size(vidFrames2_1);
[height3, width3, rgb3, num_frames3] = size(vidFrames3_1);
x1 = zeros(1, num_frames1-4);
y1 = zeros(1, num_frames1-4);
x2 = zeros(1, num_frames1-4);
y2 = zeros(1, num_frames1-4);
x3 = zeros(1, num_frames1-4);
y3 = zeros(1, num_frames1-4);
for j=1:num_frames1-4
    pixels1 = zeros(height1, width1);
    pixels2 = zeros(height1, width1);
    for i=1:rgb1    
        pixels1 = pixels1 + im2double(vidFrames1_1(:,:,i,j));
        pixels2 = pixels2 + im2double(vidFrames1_1(:,:,i,j+4));
    end
    [pixdiffmax, indices] = max(abs(pixels1 - pixels2), [], 'linear');
    x1(j) = indices(1);
    y1(j) = indices(2);
end
for j=1:num_frames1-4
    pixels1 = zeros(height2, width2);
    pixels2 = zeros(height2, width2);
    for i=1:rgb2    
        pixels1 = pixels1 + im2double(vidFrames2_1(:,:,i,j));
        pixels2 = pixels2 + im2double(vidFrames2_1(:,:,i,j+4));
    end
    [pixdiffmax, indices] = max(abs(pixels1 - pixels2), [], 'linear');
    x2(j) = indices(1);
    y2(j) = indices(2);
end
for j=1:num_frames1-4
    pixels1 = zeros(height3, width3);
    pixels2 = zeros(height3, width3);
    for i=1:rgb3    
        pixels1 = pixels1 + im2double(vidFrames3_1(:,:,i,j));
        pixels2 = pixels2 + im2double(vidFrames3_1(:,:,i,j+4));
    end
    [pixdiffmax, indices] = max(abs(pixels1 - pixels2), [], 'linear');
    y3(j) = indices(1);
    x3(j) = indices(2);
end
x1 = x1 - mean(x1);
y1 = y1 - mean(y1);
x2 = x2 - mean(x2);
y2 = y2 - mean(y2);
x3 = x3 - mean(x3);
y3 = y3 - mean(y3);
%% 
[height1, width1, rgb1, num_frames1] = size(vidFrames1_2);
[height2, width2, rgb2, num_frames2] = size(vidFrames2_2);
[height3, width3, rgb3, num_frames3] = size(vidFrames3_2);
x1 = zeros(1, num_frames1);
y1 = zeros(1, num_frames1);
x2 = zeros(1, num_frames1);
y2 = zeros(1, num_frames1);
x3 = zeros(1, num_frames1);
y3 = zeros(1, num_frames1);
%%
for j=1:num_frames1
    pixels1 = zeros(height1, width1);
    for i=1:rgb1
        pixels1 = pixels1 + im2double(vidFrames1_2(:,:,i,j));
    end
    [pixmax, indices] = max(pixels1, [], 'linear');
    x1(j) = indices(1);
    y1(j) = indices(2);
end
for j=1:num_frames1
    pixels1 = zeros(height2, width2);
    for i=1:rgb2    
        pixels1 = pixels1 + im2double(vidFrames2_2(:,:,i,j));
    end
    [pixdiffmax, indices] = max(pixels1, [], 'linear');
    x2(j) = indices(1);
    y2(j) = indices(2);
end
for j=1:num_frames1
    pixels1 = zeros(height3, width3);
    for i=1:rgb3    
        pixels1 = pixels1 + im2double(vidFrames3_2(:,:,i,j));
    end
    [pixdiffmax, indices] = max(pixels1, [], 'linear');
    y3(j) = indices(1);
    x3(j) = indices(2);
end
% x1 = x1(1:length(x1)-23);
% y1 = y1(1:length(y1)-23);
x1 = x1 - mean(x1);
y1 = y1 - mean(y1);
% x2 = x2(21:length(x2)-3);
% y2 = y2(21:length(y2)-3);
x2 = x2 - mean(x2);
y2 = y2 - mean(y2);
% x3 = x3(23:length(x3)-1);
% y3 = y3(23:length(y3)-1);
x3 = x3 - mean(x3);
y3 = y3 - mean(y3);
%% 
M = [x1;y1;x2;y2;x3;y3];
x = linspace(0,width1,length(x1));
t = linspace(0,height1,length(x1));
plot(t, y3, 'r.','Markersize',10);
xlabel('time')
ylabel('position')
%%
[U,S,V] = svd(M,'econ');
X_rank1 = S(1,1)*U(:,1)*V(:,1)';
plot(X_rank1(1,:), X_rank1(2,:),'r.','MarkerSize',10)

%%
subplot(1,1,1)
% plot(x,U(:,1),'b',x,U(:,2),'--r',x,U(:,3),':k','Linewidth',2)
% xlabel('x')
% set(gca,'Fontsize',16)
% subplot(2,1,2)
plot(t,V(:,1),'b',t,V(:,2),'--r',t,V(:,3),':k','Linewidth',2)
legend('mode 1','mode 2','mode 3')
xlabel('t')
set(gca,'Fontsize',16)
% X_rank2 = S(1:2,1:2)*U(:,1:2)*V(:,1:2)';
% plot(X_rank2(1,:),X_rank2(2,:),'r.','MarkerSize',10)

%% 

% for j=1:num_frames
% X=vidFrames1_1(:,:,:,j);
% imshow(X); drawnow
% end