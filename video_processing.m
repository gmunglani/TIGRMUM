function L = video_processing(pathf,fname,stp,smp,timestep,L,Cmin,Cmin_tmp,Cmax)

movie = [pathf '/' fname '_ratio.avi'];

V = VideoWriter(movie);
V.FrameRate = 50;
open(V);    

for count = 1:size(L,3)
    h = figure('visible', 'off');
    map = colormap(jet(255));
    map = vertcat([0 0 0],map);
    imshow(imgaussfilt(L(:,:,count),1.5),map); 
    hcb = colorbar;
    set(hcb,'FontSize',14)
    set(hcb,'YTick', [0 255])
    set(hcb,'YTickLabel', {num2str(Cmin_tmp),num2str(Cmax)})
    set(gca, 'CLim', [0,255]);
    disp(['Video Processing:' num2str((count+stp-1))]);
    txtstr = strcat('Time(s): ',num2str((count+stp-1)*timestep));
    text(10,10,txtstr,'color','white')
    drawnow();
    frame = getframe(gcf);
    writeVideo(V,frame);
    close(h);
end

close(V);
disp(['Cmax:' num2str(Cmax)]);
disp(['Cmin:' num2str(Cmin*Cmax)]);

