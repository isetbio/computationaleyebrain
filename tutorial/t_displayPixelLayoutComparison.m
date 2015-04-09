%% t_displayPixelLayoutComparison
%    In this script, we simulate a font on different types of displays.
%    Then, we compare the scene radiance and optical images at different
%    viewing distance
%
%  (HJ) ISETBIO TEAM, 2015

%% Init
ieInit; % init a new isetbio session
vd = [0.01 2];
nDist = length(vd);
font = fontCreate; % create a default font
vcNewGraphWin;

%% CRT Display
d = displayCreate('CRT-HP');
scene = sceneFromFont(font, d);
% vcAddObject(scene); sceneWindow;

oi = oiCreate('human');
for ii = 1 : nDist
    scene = sceneSet(scene, 'distance', vd(ii));
    oi = oiCompute(scene, oi);
    % vcAddObject(oi); oiWindow;
    
    subplot(nDist, 3, 1 + 3*(ii-1)); imshow(oiGet(oi, 'rgb image'));
end

%% LCD Chevron Shape Display
d = displayCreate('Dell-Chevron');
scene = sceneFromFont(font, d);
% vcAddObject(scene); sceneWindow;

for ii = 1 : nDist
    scene = sceneSet(scene, 'distance', vd(ii));
    oi = oiCompute(scene, oi);
    % vcAddObject(oi); oiWindow;
    
    subplot(nDist, 3, 2 + 3*(ii-1)); imshow(oiGet(oi, 'rgb image'));
end

%% OLCD Display
d = displayCreate('OLED-Sony');
scene = sceneFromFont(font, d);
% vcAddObject(scene); sceneWindow;

oi = oiCreate('human');
for ii = 1 : nDist
    scene = sceneSet(scene, 'distance', vd(ii));
    oi = oiCompute(scene, oi);
    % vcAddObject(oi); oiWindow;
    
    subplot(nDist, 3, 3 + 3*(ii-1)); imshow(oiGet(oi, 'rgb image'));
end