clear all
close all

fname = '~/Documents/MATLAB/Pollen_tube_diameter/Hannes_images/Lily_19.jpg';

info = imfinfo(fname);
num_images = numel(info);
A = imread(fname,'jpg');
[B,pos] = imcrop(A);
C = imadjust(mat2gray(rgb2gray(B)));  
    
%temp = fspecial('gaussian', [200 200]);
%D = imfilter(C,temp,'replicate');
D = imsharpen(C);
    
%e1 = strel('rectangle',[6 6]);
%e2 = strel('disk',3);

e1 = strel('rectangle',[8 8]);
e3 = strel('rectangle',[6 6]);
e2 = strel('disk',7);
E = imdilate(D,e1);
E = imerode(E,e2);
F = imcomplement(imbinarize(D,0.72));
G = imbinarize(imcomplement(D),0.5);

J = edge(E,'Canny',0.1);
H = imcomplement(F) + G;
M = imdilate(H,e3);
M = imerode(M,e3);
%Q = M + P;
stats = regionprops(M,'Orientation','MajorAxisLength', 'BoundingBox', ...
    'MinorAxisLength', 'Eccentricity', 'Centroid','Area','FilledImage');
L = stats(1).FilledImage;
I = bwboundaries(L,'holes');
P = bwmorph(L,'skel',inf);

%[yco,xco] = ind2sub(size(P),find(P));

% xco(xco == max(xco)) = [];yco(xco == max(xco)) = [];
% xco(xco == min(xco)) = [];yco(xco == min(xco)) = [];
% yco(yco == max(yco)) = [];xco(yco == max(yco)) = [];
% yco(yco == min(yco)) = [];xco(yco == min(yco)) = [];

%[f,fitco] = polyfit(xco,yco,3);

%temp = fspecial('gaussian', [20 20]);
%T = imfilter(imcomplement(Q),temp,'replicate');

%O = zeros(size(
%blot = 0.7;
%while(isequal(N,L) != 0)
%    N = nlfilter(O,[3 3],'mean2');
%    N(N>blot) = 1.00;
%    N(N<blot) = 0.00;
%    O = N;
%end

for x = 1:numel(I)
    temp1(x) = size(I{x},1);
end
[pos posI] = max(temp1);

bound = I{posI};
maxx = size(L,1); maxy = size(L,2); 

hullo = convhull(bound(:,1),bound(:,2));
ver = [bound(hullo,1) bound(hullo,2)];

xo1 = find(ver(:,1) == maxx); if(size(xo1,1) > 1) [diam vero start] = edge_quant(ver,xo1); end
yo1 = find(ver(:,2) == maxy); if(size(yo1,1) > 1) [diam vero start] = edge_quant(ver,yo1); end 
xo2 = find(ver(:,1) == 0); if(size(xo2,1) > 1) [diam vero start] = edge_quant(ver,xo2); end
yo2 = find(ver(:,2) == 0); if(size(yo2,1) > 1) [diam vero start] = edge_quant(ver,yo2); end

veron = circshift(vero,-start);
%xo1 = find(bound(:,1) == maxx); if(size(xo1,1) > 1) edgex = [bound(xo1,1) bound(xo1,2)]; end
%yo1 = find(bound(:,2) == maxy); if(size(yo1,1) > 1) edgey = [bound(yo1,1) bound(yo1,2)]; end 
%xo2 = find(bound(:,1) == 0); if(size(xo2,1) > 1) edgex = [bound(xo2,1) bound(xo2,2)]; end
%yo2 = find(bound(:,2) == 0); if(size(yo2,1) > 1) edgey = [bound(yo2,1) bound(yo2,2)]; end

% boundc = bound;
% if (isempty(edgex)) diam = norm(max(edgey) - min(edgey)); bound(min(edgey(:,1))+1:max(edgey(:,1))-1,:) = [];
% elseif (isempty(edgey)) diam = norm(max(edgex) - min(edgex)); bound(min(edgex(:,2))+1:max(edgex(:,2))-1,:) = [];
% end

major = [stats(1).Centroid(1) - stats(1).MajorAxisLength*0.5 * cos(pi*stats(1).Orientation/180) ...
        stats(1).Centroid(2) + stats(1).MajorAxisLength*0.5 * sin(pi*stats(1).Orientation/180)];

tip = [];
for i = 1:length(veron)    
    dist = sqrt((veron(i,1) - major(1,2))^2 + (veron(i,2) - major(1,1))^2);
    if (dist < diam*0.25) tip = [tip; veron(i,:)]; end
end
%  
lin = [(1:size(tip,1))' (2:size(tip,1)+1)']; lin(end,2)=1;
%lin = [(1:size(veron,1))' (2:size(veron,1)+1)']; lin(end,2)=1;
% smoother = smooth(veron(:,1),veron(:,2),'sgolay');
% verons = [veron(:,1) smoother];
% [fitty, structy] = polyfit(veron(:,1),veron(:,2),5);
% 
curve = LineCurvature2D(tip,lin);
[valt,loct] = min(curve);
tipf = tip(loct,:);

posx = find(bound(:,1) == tipf(1));
posy = find(bound(:,2) == tipf(2));
inter = intersect(posx,posy);
boundn = circshift(bound,(-ceil(length(bound)*0.5)-inter));

check1 = []; check2 = []; val1 = []; val2 = []; posc1 = []; posc2 = [];

for i = 1:ceil(length(boundn)*0.5)    
    dist = pdist2(boundn(i,:),tipf);
    if (dist < diam) val1 = [val1 dist]; posc1 = [posc1 i]; check1 = [check1; boundn(i,:)]; end
end
for i = ceil(length(boundn)*0.5+1):length(boundn)    
    dist = pdist2(boundn(i,:),tipf);
    if (dist < diam) val2 = [val2 dist]; posc2 = [posc2 i]; check2 = [check2; boundn(i,:)]; end
end

check1s = circshift(check1,-1); 
check2s = circshift(check2,-1); 
len1 = sum(diag(pdist2(check1(1:end-1,:),check1s(1:end-1,:))));
len2 = sum(diag(pdist2(check2(1:end-1,:),check2s(1:end-1,:))));

change = movmean(diff(val1),3);
for k = 1:length(change)
    if(change(k) > 0) turnin = k; break;end
end
for k = length(change):-1:1
    if(change(k) > 0) turnout = k; break;end
end
tipb = [mean([check1(1,1) check2(end,1)]) mean([check1(1,2) check2(end,2)])];

T = drawline(L,check1(turnin,1),check1(turnin,2),check1(turnout,1),check1(turnout,2),1);
stats3 = regionprops(T,'FilledImage');
U = stats3(1).FilledImage;
W = bwmorph(U,'skel',inf);
Y = bwmorph(U,'thin',inf);

boundn(min(posc1)+turnin:min(posc1)+turnout,:) = [];
tipr = []; far = [];
for i = 1:length(boundn)    
    dist = pdist2(boundn(i,:),tipf);
    if (dist < 0.25*diam) tipr = [tipr; boundn(i,:)]; far = [far pdist2(boundn(i,:),tipb)]; end
end
[posn, farn] = max(far);
tipm = tipr(farn,:);

posx1 = find(boundn(:,1) == tipm(1));
posy1 = find(boundn(:,2) == tipm(2));
inter1 = intersect(posx1,posy1); 

total1 =[]; total2 =[];
for i = inter1:length(boundn)
    total1 = [total1; boundn(i,:)];
    if (boundn(i,2) == maxy) break;end
end
for i = inter1-1:-1:1    
    total2 = [total2; boundn(i,:)];
    if (boundn(i,2) == maxy) break;end
end

total1s = circshift(total1,-1); 
total2s = circshift(total2,-1); 
lenf1 = sum(diag(pdist2(total1(1:end-1,:),total1s(1:end-1,:))));
lenf2 = sum(diag(pdist2(total2(1:end-1,:),total2s(1:end-1,:))));

[fit1, struct1] = polyfit(total1(:,2),total1(:,1),7);
[fit2, struct2] = polyfit(total2(:,2),total2(:,1),7);
xfit = -10:2:maxy;
line1f = polyval(fit1,xfit);
line2f = polyval(fit2,xfit);

diamf = abs(line1f - line2f);

for l = 1:length(line1f)
    cline(l) = mean([line1f(l) line2f(l)]);
end

dy = gradient(cline);
dx = gradient(xfit);

% 
% for j = 1:length(check1)
%     d = abs(det([tipf-tipb,check1(j)-tipb]))/abs(tipf-tipb);
% end
%alphao = alphaShape(bound(:,1),bound(:,2),25);
%[bf,poser] = boundaryFacets(alphao);

figure
plot(boundn(:,2), boundn(:,1), 'r');
hold on
%plot(check1(:,2), check1(:,1), 'k*');
%plot(check2(:,2), check2(:,1), 'm*');
%plot(tipf(:,2), tipf(:,1), 'b*', 'LineWidth', 5);
%plot(tipb(:,2), tipb(:,1), 'b*', 'LineWidth', 5);
%plot(tipr(:,2), tipr(:,1), 'k*', 'LineWidth', 2);
%plot(tipm(:,2), tipm(:,1), 'r*', 'LineWidth', 2);

plot(total1(:,2), total1(:,1), 'k', 'LineWidth', 2);
plot(total2(:,2), total2(:,1), 'b', 'LineWidth', 2);

plot(xfit, line1f, 'k*', 'LineWidth', 2);
plot(xfit, line2f, 'b*', 'LineWidth', 2);
plot(xfit, cline, 'r*', 'LineWidth', 2);
quiver(xfit,cline,-dy,dx, 0, 'm')

%plot(ver(:,2), ver(:,1), 'b*', 'LineWidth', 2);
%plot(veron(:,2), veron(:,1), 'b*');
%hold on
%plot(major(:,2), major(:,1), 'r*');
%plot(tip(:,2), tip(:,1), 'k*', 'LineWidth', 5);
%ellipse_view(stats)
%plot(poser(:,2),poser(:,1), 'r*', 'LineWidth', 2)

%plot(xco,yco,'g', 'LineWidth', 2)


figure
subplot(3,3,1)
imshow(B);
subplot(3,3,2)
imshow(F);
subplot(3,3,3)
imshow(G);
subplot(3,3,4)
imshow(H);
subplot(3,3,5)
imshow(M);
subplot(3,3,6)
imshow(L); 
subplot(3,3,7)
imshow(U); 
subplot(3,3,8)
%imshow(W);
subplot(3,3,9)
%imshow(Y);

figure
subplot(1,2,1)
imshow(B)
subplot(1,2,2)
mask = cast(imcomplement(M), class(B)); 
img_masked = B .* repmat(mask, [1 1 3]); 
%mask2 = cast(P, class(B)); 
%img_masked = img_masked .* repmat(imcomplement(mask2), [1 1 3]);  
imshow(img_masked);

figure
plot(diamf)