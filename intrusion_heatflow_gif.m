%% Intrusion heat flow, gif making
% William Davis 13/07/16

% Pascal and Olesen 2009
clc, clearvars, close all

load('T_mat.mat')

%// Create file name.
FileName = 'Diffusion2.gif';

%% Time stepping
for i = 1:nt+1 
    figure(1)
    contourf(X/1000,Y/1000,T_mat(:,:,i),'LineStyle','none','LevelStep',100) % 1600/9
    c = colorbar('southoutside','Limits',[0,1600]);
    c.Label.String = 'Temp. (\circC)';
    colormap(jet)
    ylabel('Depth, km')
    axis ij
    axis equal

    line([0,x_length/kmSca],[crust_thickness/kmSca,crust_thickness/kmSca],'Color','k',...
       'LineWidth',1)
    line([(x_length-diapir_length)/(2*kmSca),(x_length-diapir_length)/(2*kmSca)],...
       [crust_thickness/kmSca,y_length/kmSca],'Color','k',...
       'LineWidth',1)
    line([(x_length+diapir_length)/(2*kmSca),(x_length+diapir_length)/(2*kmSca)],...
       [crust_thickness/kmSca,y_length/kmSca],'Color','k',...
       'LineWidth',1)

    timestep = num2str(i,'%03.0f');
    ageMa = num2str(round(i*dtMa,1),'%04.1f');
    title([' FD, time step: ',timestep,'     Time(Ma): ',ageMa])
    xlabel ('Distance, km'), ylabel('Depth, km')

    % Making gif
    set(gca,'nextplot','replacechildren','visible','on')
    drawnow
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);

    % On the first loop, create the file. In subsequent loops, append.
    if i==1
        imwrite(imind,cm,FileName,'gif','DelayTime',0,'loopcount',inf);
    else
        imwrite(imind,cm,FileName,'gif','DelayTime',0,'writemode','append');
    end
end

%{
for k = 1:numel(ImageCell)

    if k ==1

%// For 1st image, start the 'LoopCount'.
        imwrite(ImageCell{k},FileName,'gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(ImageCell{k},FileName,'gif','WriteMode','append','DelayTime',1);
    end

end
%}