% s_csfScene

%% Create a Gabor patch
params.freq = 3;   % Per FOV
params.contrast = 1;
params.ph = pi;
params.ang = 0;
params.row = 128; params.col = 128;
params.GaborFlag = .2;

scene = sceneCreate('harmonic',params);
scene = sceneSet(scene,'fov',1);
vcAddObject(scene);
sceneWindow;

%% END