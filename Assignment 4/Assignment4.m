[images, testlabels] = mnist_parse('train-images.idx3-ubyte', 'train-labels.idx1-ubyte');
sz = size(images);
images = reshape(images, [sz(1)*sz(2), sz(3)]);
%%
im_wave = a4_wavelet(images);
[U,S,V] = svd(im_wave, 'econ');

%%
figure(2)
subplot(2,1,1)
plot(diag(S),'ko','Linewidth',2)
set(gca,'Fontsize',16,'Xlim',[0 80])
subplot(2,1,2)
semilogy(diag(S),'ko','Linewidth',2)
set(gca,'Fontsize',16,'Xlim',[0 80]) % 10 modes will have a lot of info

%%

figure(3)
scatter3(V(:,2),V(:,3),V(:,5),1,testlabels);

%%

figure(4)
val1 = 0;
val2 = 7;
X = V(testlabels == val1 | testlabels == val2, :);
lim_labels = testlabels(testlabels == val1 | testlabels == val2);
scatter3(X(:,2),X(:,3),X(:,5),1,lim_labels);

%%
im1 = images(:, testlabels == val1);
im2 = images(:, testlabels == val2);
im1_wave = a4_wavelet(im1);
im2_wave = a4_wavelet(im2);
[U1,S1,V1] = svd([im1_wave im2_wave], 'econ');
%%
feature = 10;
n1 = size(im1_wave,2);
n2 = size(im2_wave,2);
digits = S1*V1'; % projection onto principal components: X = USV' --> U'X = SV'
dig1 = digits(1:feature,1:n1);
dig2 = digits(1:feature,n1+1:n1+n2);
m1 = mean(dig1,2);
m2 = mean(dig2,2);
Sw = 0;
for k = 1:n1
    Sw = Sw + (dig1(:,k) - m1)*(dig1(:,k) - m1)';
end
for k = 1:n2
    Sw =  Sw + (dig2(:,k) - m2)*(dig2(:,k) - m2)';
end
Sb = (m1-m2)*(m1-m2)';

[V2, D] = eig(Sb,Sw);
[lambda, ind] = max(abs(diag(D)));
w = V2(:,ind);
w = w/norm(w,2);
v1 = w'*dig1;
v2 = w'*dig2;

if mean(v1) > mean(v2)
    w = -w;
    v1 = -v1;
    v2 = -v2;
end
%%

figure(5)
plot(v1,zeros(n1),'ob','Linewidth',2)
hold on
plot(v2,ones(n2),'dr','Linewidth',2)
ylim([0 1.2])

%%

sort1 = sort(v1);
sort2 = sort(v2);

t1 = length(sort1);
t2 = 1;
while sort1(t1) > sort2(t2)
    t1 = t1-1;
    t2 = t2+1;
end
threshold = (sort1(t1) + sort2(t2))/2;

%%

[testimages, testlabels] = mnist_parse('t10k-images.idx3-ubyte', 't10k-labels.idx1-ubyte');
sz = size(testimages);
testimages = reshape(testimages, [sz(1)*sz(2), sz(3)]);

%%

TestNum = size(testimages,2);
Test_wave = a4_wavelet(testimages); % wavelet transform
TestMat = U1'*Test_wave; % PCA projection
pval = w'*TestMat; 

%%
Mdl = fitcsvm(images,labels);
testlabels = predict(Mdl,test);