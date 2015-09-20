%% s_colorContours
%
%   Compute color discrimination contours.
%
%   This calculates threshold increment for a uniform patch on a solid
%   background. The noise arises from photon absorptions and second-site
%   additive noise.
%
%  (HJ) ISETBIO TEAM, 2014

%% Init parameters
ieInit;

%% Set up the stimulus parameters
ref     = [0 0 0];    % Background alone, no contrast
dirList = [0 15 25 30 40 42 45 47 50 55 65 90 115 135 150];  % color directions in degrees
% dirList = [0 30 45 60 90 135];
dirList = [dirList dirList - 180];    % Color directions made symmetric
sensorSz  = [45, 45]; % Number of cones

% Set the cone backgrounds maybe other parameters up here ...

%% Compute classification accuracy
%  Set up command
cmd = '[thresh, params] = ccThresh(ref, dirList(ii), params);';
params.ccParams.sensorSz = sensorSz;

% Cone densities for blank, L, M and S
params.ccParams.cone = coneCreate;
params.ccParams.cone.spatialDensity = [0 0.6 0.3 0.1];
d = displayCreate('OLED-Sony', 'wave', 400:10:700);
d = displaySet(d, 'gamma', 'linear');
params.ccParams.d = d;

cprintf('*Keywords', 'Color Discrimination Contour\n');
dataDir = '~/SimResults/ColorContour/Detection/';
for ii = 1 : length(dirList)
    % print progress
    fprintf('\t(%d deg L-M)\n', dirList(ii));
    
    % compute threshold
    eval(cmd);
    fname = fullfile(dataDir, sprintf('%d.mat', dirList(ii)));
    save(fname, 'thresh', 'params');
    
    % print results
    fprintf('\t\tThreshold:%.4f\n', thresh);
end

%% Plot color contour
threshPts = zeros(length(dirList), 3);

for ii = 1 : length(dirList)
    fName = fullfile(dataDir, sprintf('%d.mat', dirList(ii)));
    if ~exist(fName, 'file'), continue; end
    data = load(fName);
    curDir = dirList(ii); % current direction
    threshPts(ii, :) = ref + data.thresh*[cosd(curDir) sind(curDir) 0];
end

% fit ellipse and plot
figure; hold on;
plot(threshPts(:,1), threshPts(:,2), 'xr'); % plot points
[zg, ag, bg, alphag] = fitellipse(threshPts(:,1:2));
plotellipse(zg, ag, bg, alphag, 'b--')

axis equal; grid on;
xlabel('L contrast'); ylabel('M contrast');