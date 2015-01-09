%% s_colorContours
%   Compute color discrimination contours.
%
%   This calculates threshold increment for a uniform patch on a solid
%   background.  The noise arises from photon absorptions and second-site
%   additive noise.
%
%  (HJ) ISETBIO TEAM, 2014

%% Init parameters
s_initISET;

%% Set up the stimulus parameters
ref     = [0 0 0];    % Background alone, no contrast
dirList = [0 15 25 30 40 42 45 47 50 55 65 90 115 135 150];  % color directions in L-M
dirList = [dirList dirList - 180];    % Color directions made symmetric
cropSz  = 22;    % Number of cones is 2*cropSz + 1

% Set the cone backgrounds maybe other parameters up here ...

%% Compute classification accuracy
%  Set up proclus command
cmd = '[thresh, expData] = ccThresh(ref, dirList(jobindex), params);';
cmd = [cmd 'save(sprintf(''~/ccContour%d.mat'',jobindex));'];
params.ccParams.cropSz = cropSz;
params.ccParams.rgbDensities = [0 0.6 0.3 .1];  % Blank, L, M and S proportions

try % try use proclus
    sgerun2(cmd, 'colorContour', 1, 1:length(dirList));
catch % compute locally
    for ii = 1 : length(dirList)
        [thresh, expData] = ccThresh(ref, dirList(ii), params);
        save(sprintf('~/ccContour%d.mat', ii));
    end
end

%% Plot color contour
threshPts = zeros(length(dirList), 3);

for ii = 1 : length(dirList)
    fName = sprintf('./ccContour%d.mat', ii);
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