clear all
close all

% Path to h5 file
path = '/home/gm/Documents/Scripts/MATLAB/Tip_results'; % Input folder path
fname = 'YC11'; % File name 
stp = 1;
smp = 210;
specific = [84 190]; % Frames to change

% Other Options
register = 1; % Register image
union = 1; % Take the union of the two image masks

[optimizer, metric] = imregconfig('multimodal');

pathf = [path '/' fname '/' fname];
BT1 = h5read([pathf '_back_proc.h5'],'/BT1');
BT2 = h5read([pathf '_back_proc.h5'],'/BT2');
M = h5read([pathf '_back_proc.h5'],'/M');
posfront = h5readatt([pathf '_back_proc.h5'],'/M','Crop');
Bedge = h5readatt([pathf '_back_proc.h5'],'/M','Bedge');
intensity1 = h5read([pathf '_back_proc.h5'],'/intensity1');
intensity2 = h5read([pathf '_back_proc.h5'],'/intensity2');
Bsum = h5read([pathf '_back_proc.h5'],'/Bsum');

YFP = h5read([pathf '_back.h5'],'/YFP');
CFP = h5read([pathf '_back.h5'],'/CFP');
    
for count = 1:length(specific)
    % Read image and add bleach correction
    disp(['Pre Processing:' num2str(count)]);
    
    A1 = YFP(:,:,count)';
    A2 = CFP(:,:,count)';
    
    B1 = imcrop(A1,posfront);
    B2 = imcrop(A2,posfront);
    
    if (register == 1) B2 = imregister(B2,B1,'translation',optimizer,metric); end
    
    Bu1 = uint8(single(B1).*255/4095);
    Bu2 = uint8(single(B2).*255/4095);
    
    level1 = graythresh(Bu1); Bl1 = imbinarize(Bu1,level1);
    level2 = graythresh(Bu2); Bl2 = imbinarize(Bu2,level2);
    
    B1 = B1.*uint16(Bl1);
    B2 = B2.*uint16(Bl2);
    
    % Orient image
    if (count==stp) type = find_orient(B1); end
    if (type == 1) B1 = imrotate(B1,-90); B2 = imrotate(BN2,-90);
    elseif (type == 3) B1 = imrotate(B1,90); B2 = imrotate(BN2,90);
    elseif (type == 4) B1 = imrotate(B1,180); B2 = imrotate(BN2,180);
    end
    
    B1 = B1(:,1:end-Bedge);
    B2 = B2(:,1:end-Bedge);
    
    % Union of images to ensure perfect overlap
    countn = specific(count);
    B = B1.*B2; B(B>0) = 1;
    if (union == 1)
        BT1(:,:,countn) = B1.*B;
        BT2(:,:,countn) = B2.*B;
    else
        BT1(:,:,countn) = B1;
        BT2(:,:,countn) = B2;
    end
        
    % Analysis
    Bsum(countn,1) = sum(B(:));
    Bsum(countn,2) = nnz(B1);
    Bsum(countn,3) = nnz(B2);
    
    Br1 = reshape(BT1(:,:,countn),[numel(BT1(:,:,countn)),1]);
    Br2 = reshape(BT2(:,:,countn),[numel(BT2(:,:,countn)),1]);
    
    intensity1(countn) = median(nonzeros(Br1));
    intensity2(countn) = median(nonzeros(Br2));
end

% Create ratio image and output values
M = BT1./BT2;
M(M==Inf) = 0;
M(isnan(M)) = 0;

name = [pathf '_back_proc.h5']; 
delete(name);
h5create(name,'/M',size(M));
h5create(name,'/BT1',size(BT1));
h5create(name,'/BT2',size(BT2));
h5create(name,'/intensity1',size(intensity1));
h5create(name,'/intensity2',size(intensity2));
h5create(name,'/Bsum',size(Bsum));

h5write(name,'/M',uint16(M));
h5writeatt(name,'/M','Crop',posfront);
h5writeatt(name,'/M','Bedge',Bedge);

h5write(name,'/BT1',uint16(BT1));
h5write(name,'/BT2',uint16(BT2));
h5write(name,'/intensity1',intensity1);
h5write(name,'/intensity2',intensity2);
h5write(name,'/Bsum',Bsum);

% Plotting signal results and decay
h = figure;
subplot(2,2,1)
plot(stp:smp,Bsum(stp:smp,1),'b');
hold on
plot(stp:smp,Bsum(stp:smp,2),'r');
plot(stp:smp,Bsum(stp:smp,3),'g');
axis([stp smp 0.5*max(Bsum(:,1)) max(Bsum(:,1))])
title('Total Number of Pixels');

% Scale each channel by maximum intensity for viewing
BT1max = max(BT1(:));
BT1 = BT1./BT1max;
BT2 = BT2./BT1max;
    
[Mmin Mmax Mmin_prc Mmax_prc] = channel_analysis(M,smp);
[B1min B1max B1min_prc B1max_prc] = channel_analysis(BT1,smp);
[B2min B2max B2min_prc B2max_prc] = channel_analysis(BT2,smp);

subplot(2,2,2)
hold on
plot(stp:smp,Mmax(stp:smp),'b*')
plot(stp:smp,Mmax_prc(stp:smp),'bd')
plot(stp:smp,Mmin(stp:smp),'r*')
plot(stp:smp,Mmin_prc(stp:smp),'rd')
grid on
axis([stp smp 0 max(Mmax)])
set(gca,'YMinorTick','on')
title('Percentiles of Ratio image');

subplot(2,2,3)
hold on
plot(stp:smp,B1max(stp:smp),'b*')
plot(stp:smp,B1max_prc(stp:smp),'bd')
plot(stp:smp,B1min(stp:smp),'r*')
plot(stp:smp,B1min_prc(stp:smp),'rd')
grid on
axis([stp smp 0 1])
set(gca,'YMinorTick','on');
title('Percentiles of YFP image');

subplot(2,2,4)
hold on
plot(stp:smp,B2max(stp:smp),'b*')
plot(stp:smp,B2max_prc(stp:smp),'bd')
plot(stp:smp,B2min(stp:smp),'r*')
plot(stp:smp,B2min_prc(stp:smp),'rd')
grid on
axis([stp smp 0 1])
set(gca,'YMinorTick','on');
title('Percentiles of CFP image');
hold off
    
savefig(h,[pathf '_back_proc2.fig']);
