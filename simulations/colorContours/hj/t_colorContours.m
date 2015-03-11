%% t_colorContours
%   Compute color discrimination contours.
%
%   This calculates threshold increment for a uniform patch on a solid
%   background.  The noise arises from photon absorptions and second-site
%   additive noise.
%
%  (HJ) ISETBIO TEAM, 2015

%% Init parameters
ieInit;

%% Set up the stimulus parameters
ref     = [0 0 0];    % Background contrast, no contrast
dirList = [0 40 45 50 90 135];  % color directions in L-M
dirList = [dirList dirList - 180];    % Color directions made symmetric
cropSz  = 22;    % Number of cones is 2*cropSz + 1

% Set the cone backgrounds maybe other parameters up here ...

%% Compute classification accuracy
%  Set up parameters
params.ccParams.cropSz = cropSz; % Size of patch used for classification
params.ccParams.rgbDensities = [0 .6 .3 .1]; %Ratio of K, L, M and S cones

d = displayCreate('LCD-Apple');
d = displaySet(d, 'gamma', repmat(linspace(0,1,256)', [1 3]));
params.ccParams.d = d;

% Compute detection threshold for each direction
threshPts = zeros(length(dirList), 3);

for ii = 1 : length(dirList)
    fprintf('Computing color contour for degree: %d\n', dirList(ii));
    tic;
    [thresh, expData] = ccThresh(ref, dirList(ii), params);
    direction = [cosd(dirList(ii)) sind(dirList(ii)) 0];
    threshPts(ii, :) = ref + thresh * direction;
    toc;
end

%% Plot color contour
% fit ellipse and plot
figure; hold on;
plot(threshPts(:,1), threshPts(:,2), 'xr'); % plot points
threshPts = [threshPts; -threshPts]; % symmetric by origin
[zg, ag, bg, alphag] = fitellipse(threshPts(:,1:2));
plotellipse(zg, ag, bg, alphag, 'b--')

axis equal; grid on;
xlabel('L contrast'); ylabel('M contrast');