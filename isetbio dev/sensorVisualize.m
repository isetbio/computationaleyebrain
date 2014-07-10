function sensorVisualize(sensor, varargin)
%% Visualize the cone absorption image / time sequencing data
%    function sensorVisualize(sensor, [varargin])
%
%
%
%  (HJ) ISETBIO TEAM, 2014

%% Init
if notDefined('sensor'), error('sensor required'); end

%% Photons
v = sensorGet(sensor, 'volts');
if isempty(v), error('No absorption data found in sensor'); end

if ismatrix(v)
    vcAddAndSelectObject('sensor', sensor);
    sensorWindow;
end

if ndims(v) == 3 % time sequencing data
    v = v/max(v(:));
    tmp = sensor;
    rgb = zeros(size(v,1), size(v,2), 3, size(v,3));
    for ii = 1 : size(v,3)
        tmp = sensorSet(tmp,'volts',v(:,:,ii));
        rgb(:,:,:,ii) = sensorGet(tmp,'rgb','volts',1,0);
    end
    rgb = rgb/max(rgb(:));
    
    %% Make the running average   
    % This many milliseconds of temporal integration
    tInt = 50;
    tave = ones(1,1,1,tInt)/tInt;
    cones = convn(rgb,tave,'valid');
    
    %% Make a movie - turn this into a function
    
    figure, set(gcf, 'Color','white')
    cones = uint8(255 * cones);
    image(cones(:,:,:,1)); axis tight
    set(gca, 'nextplot','replacechildren', 'Visible','off');
    
    %# create AVI object
    nFrames = size(cones,4);
    vidObj = VideoWriter('coneMovie.avi');
    vidObj.Quality = 60;
    vidObj.FrameRate = 24;
    open(vidObj);
    
    %# create movie
    for k = 1 : nFrames
        image(cones(:,:,:,k));
        writeVideo(vidObj, getframe(gca));
    end
    close(gcf)
    
    %# save as AVI file, and open it using system video player
    close(vidObj);
end

end