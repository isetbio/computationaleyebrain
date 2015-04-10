%% t_displayPixelLayoutComparison
%    In this script, we simulate a horizontal line on different types of
%    displays. Then, we compare the scene radiance and optical images at
%    different viewing distance
%
%  (HJ) ISETBIO TEAM, 2015

%% Init
ieInit; % init a new isetbio session
vd = [0.25 1];
nDist = length(vd);
I = zeros(21, 30); I(11, :) = 1; % create image
doSub = true; % subpixel rendering

vcNewGraphWin; % new figure window

%% CRT Display
d = displayCreate('CRT-HP');
scene = sceneFromFile(I, 'rgb', [], d, [], doSub, [], [20 20]);
% vcAddObject(scene); sceneWindow;

width = sceneGet(scene, 'width');

oi = oiCreate('human');
for ii = 1 : nDist
    scene = sceneSet(scene, 'distance', vd(ii));
    scene = sceneSet(scene, 'h fov', atand(width / vd(ii)));
    oi = oiCompute(scene, oi);
    % vcAddObject(oi); oiWindow;
    
    subplot(nDist, 3, 1 + 3*(ii-1)); imshow(oiGet(oi, 'rgb image'));
    drawnow; title(sprintf('CRT-HP (%.2f m)', vd(ii)));
end

%% LCD Chevron Shape Display
d = displayCreate('Dell-Chevron');
scene = sceneFromFile(I, 'rgb', [], d, [], doSub, [], [20 20]);
% vcAddObject(scene); sceneWindow;

width = sceneGet(scene, 'width');

for ii = 1 : nDist
    scene = sceneSet(scene, 'distance', vd(ii));
    scene = sceneSet(scene, 'h fov', atand(width / vd(ii)));
    oi = oiCompute(scene, oi);
    % vcAddObject(oi); oiWindow;
    
    subplot(nDist, 3, 2 + 3*(ii-1)); imshow(oiGet(oi, 'rgb image'));
    drawnow; title(sprintf('Dell-Chevron (%.2f m)', vd(ii)));
end

%% OLCD Display
d = displayCreate('OLED-Sony');
scene = sceneFromFile(I, 'rgb', [], d, [], doSub, [], [20 20]);
% vcAddObject(scene); sceneWindow;

width = sceneGet(scene, 'width');

oi = oiCreate('human');
for ii = 1 : nDist
    scene = sceneSet(scene, 'distance', vd(ii));
    scene = sceneSet(scene, 'h fov', atand(width / vd(ii)));
    oi = oiCompute(scene, oi);
    % vcAddObject(oi); oiWindow;
    
    subplot(nDist, 3, 3 + 3*(ii-1)); imshow(oiGet(oi, 'rgb image'));
    drawnow; title(sprintf('OLED-Sony (%.2f m)', vd(ii)));
end