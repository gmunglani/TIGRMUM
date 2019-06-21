clear all
close all

fname = '/home/gm/Documents/Work/Images/DIC_tubes/WT1_1.tif';

A = imread(fname,10);
histogram(A)

%v = VideoWriter('/home/gm/Documents/pollen_channel.avi','Uncompressed AVI');
%v.FrameRate=10;
%open(v);
 %for k=1:200      % assumes 10 images to write to file
 %   A = uint8(double(imread(fname,k))*255/4096);
 %   writeVideo(v,A);
 %end
 %close(v);


thresh = multithresh(A,2);
thresh(2) = thresh(2) + 300;

G = imquantize(A,thresh);
H = zeros(size(G,1),size(G,2));
F = zeros(size(G,1),size(G,2));
Fo = bwareaopen(F,50);

F(G == 1) = 1;
H(G == 3) = 1;


SE1 = strel('line',20,90);
SE2 = strel('rectangle',[200 3]);

Hn = imerode(H,SE1);
Hn = bwareaopen(Hn,100);

Fn = imerode(F,SE1);
Fn = bwareaopen(Fn,100);

for s = 1:size(H,2)
    Hnz(s) = nnz(Hn(:,s));
    Fnz(s) = nnz(Fn(:,s));
end

[tmp,Mhl] =findpeaks(Hnz,'MinPeakProminence',100);
[tmp,Mfl] =findpeaks(Fnz,'MinPeakProminence',100);
hl = Mhl(1);
fl = Mfl(1);

off = [zeros(80,1) [1:80]'];
K = graycomatrix(Fn,'offset',off);
Kg = graycoprops(K);

[Psl,Lsl] = findpeaks(Kg.Correlation(10:end));
[tmp,Mlsl] = max(Psl);
sl = Lsl(Mlsl) + 9;

count = 0;
hh = fspecial('sobel');
hv = hh';

A1 = A(1:40,hl:fl); % Top black
A2 = A(921:960,hl+sl*17:fl+sl*17); % Pollen grains
A3 = A(781:820,hl:fl); % Round noise
A4 = A(381:420,hl+sl*14:fl+sl*14); % Smear noise 
A5 = A(141:180,hl+sl*3:fl+sl*3); % Top exit
A6 = A(601:640,hl+sl*3:fl+sl*3); % Pollen tube
A7 = A(740:780,hl+sl*2:fl+sl*2); % Pollen tip
A8 = A(201:240,hl:fl); % Empty

W = zeros(size(A,1),size(A,2));
W(1:40,hl:fl) = 1;
W(921:960,hl+sl*17:fl+sl*17) = 1;
W(781:820,hl:fl) = 1;
W(381:420,hl+sl*14:fl+sl*14) = 1;
W(141:180,hl+sl*3:fl+sl*3) = 1;
W(601:640,hl+sl*3:fl+sl*3) = 1;
W(740:780,hl+sl*2:fl+sl*2) = 1;
W(201:240,hl:fl) = 1;

figure;
AS = uint8(double(A).*255/4096);
imshowpair(AS,W,'method','blend')


%Abox = A(201:240,hl:fl);
%Abox2 = A(740:780,hl+sl*2:fl+sl*2);
%imshow(Abox2,[]) 

figure
subplot(1,4,1)
imshow(A1,[])
A1h = imfilter(A1,hh);
A1v = imfilter(A1,hv);

A1st = [sum(sum(A1v.*A1v)) sum(sum(A1v.*A1h)); sum(sum(A1h.*A1v)) sum(sum(A1h.*A1h))];
A1p = eig(A1st);
F1 = Fo(1:40,hl:fl);
F1r = regionprops(F1,'Eccentricity','Orientation','Solidity');

if exist('F1r.Eccentricity','var')
    F1ecc = F1r.Eccentricity; 
else
    F1ecc = 0; 
end
if exist('F1r.Orientation','var')
    F1or = F1r.Orientation; 
else
    F1or = 0; 
end
if exist('F1r.Solidity','var')
    F1so = F1r.Solidity; 
else
    F1so = 0; 
end

A1t = imquantize(A1,thresh);
A1med = median(A1p);
A1stdev = std(A1p);
A1p = [A1p; A1med; A1stdev; length(find(A1t == 1)); length(find(A1t == 2)); length(find(A1t == 3)); F1ecc; F1or; F1so];

subplot(1,4,2)
imshow(A2,[])
A2h = imfilter(A2,hh);
A2v = imfilter(A2,hv);
A2st = [sum(sum(A2v.*A2v)) sum(sum(A2v.*A2h)); sum(sum(A2h.*A2v)) sum(sum(A2h.*A2h))];
A2p = eig(A2st);
F2 = Fo(921:960,hl+sl*17:fl+sl*17);
F2r = regionprops(F1,'Eccentricity','Orientation','Solidity');

if exist('F2r.Eccentricity','var')
    F2ecc = F2r.Eccentricity; 
else
    F2ecc = 0; 
end
if exist('F2r.Orientation','var')
    F2or = F2r.Orientation; 
else
    F2or = 0; 
end
if exist('F2r.Solidity','var')
    F2so = F2r.Solidity; 
else
    F2so = 0; 
end

size(A2p)
A2t = imquantize(A2,thresh);
A2med = median(A2p);
A2stdev = std(A2p);
A2p = [A2p; A2med; A2stdev; length(find(A2t == 1)); length(find(A2t == 2)); length(find(A2t == 3)); F2ecc; F2or; F2so];


Ap = [A1p A2p];

Hn = imdilate(Hn,SE2);
Fn = imdilate(Fn,SE2);

%figure
%X = fft2(Hn);
%imshow(fftshift(log(abs(X) + 1)), [])

%imshow(A, []);
h = fspecial('log',[51 51],2);
B = imfilter(A,h);


% C = zeros(size(B,1),size(B,2));
% C(B > 30) = 1;
% Cn = bwareaopen(C,100);
% 
% 
% K = Cn - Hn;
% K(K < 0) = 0;
% figure
% imshow(K,[]);
% 
% level = 1;
% [c,s] = wavedec2(A,level,'coif2');
% [chd1,cvd1,cdd1] = detcoef2('all',c,s,level);
% 
% 
% E = edge(A,'sobel','vertical');
% 
% 
% SE1 = strel('line',40,90);
% SE2 = strel('line',5,90);
% 
% J = imopen(C,SE1);
% K = imopen(E,SE2);
% 
% 
% figure
% subplot(1,2,1)
% imshow(J)
% subplot(1,2,2)
% imshow(K)
%figure
%subplot(1,2,1)
%imshow(B,[])
%subplot(1,2,2)
%imshow(C)



%figure
%imshowpair(E,C);


% FFT
%Bfreq = fftshift(fft2(A));
%Bamp = log(abs(Bfreq));
%C = Bamp > 12;
%C(:,x:1280-x) = 0;

%Bfreq(C) = 0;
%Dfreq = log(abs(Bfreq));

% h = mat2gray(fspecial('gaussian',1280,60));
% h(1121:1280,:) = [];
% h(1:160,:) = [];
% Dfreq = Bfreq.*h;
% E = uint16(abs(ifft2(fftshift(Dfreq))));
% 
% figure
% subplot(1,2,1)
% imshow(E,[]);
% subplot(1,2,2)
% imshow(A-E,[])
% 
% figure
% subplot(1,2,1)
% imshow(log(abs(Dfreq)),[])
% subplot(1,2,2)
% imshow(log(abs(Bfreq)),[])

% HOUGH
% [D, theta, rho] = hough(C);
% peaks = houghpeaks(D,5,'threshold',ceil(0.01*max(D(:))));
% figure
% imshow(D,[],'XData',theta,'YData',rho);
% axis on, axis normal, hold on;
% 
% lines = houghlines(D,theta,rho,peaks,'FillGap',5);
% x = theta(peaks(:,2)); y = rho(peaks(:,1));
% plot(x,y,'s','color','white');

% figure
% imshow(C,[]);
% hold on
% max_len = 0;
% for k = 1:length(lines)
%   xy = [lines(k).point1; lines(k).point2];
%   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
% 
%    % Plot beginnings and ends of lines
%    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
% 
%    % Determine the endpoints of the longest line segment
%    len = norm(lines(k).point1 - lines(k).point2);
%    if ( len > max_len)
%       max_len = len;
%       xy_long = xy;
%    end
% end
% 

