%% s_colorContours
%    Compute color discrimination contours
%
%  (HJ) ISETBIO TEAM, 2014

%% Init parameters
s_initISET;
ref     = [0 0 0];
dirList = 0:10:359;  % directions
cropSz  = 24;

%% Compute classification accuracy
%  Set up proclus command
cmd = '[thresh, expData] = ccThreshold(ref, dirList(jobindex), params);';
cmd = [cmd 'save(sprintf(''~/ccContour%d.mat'',jobindex));'];
params.ccParams.cropSz = cropSz;
try % try use proclus
    sgerun2(cmd, 'colorContour', 1, 1:length(simParams));
catch % compute locally
    for ii = 1 : length(dirList)
        [thresh, expData] = ccThreshold(ref, dirList(ii), params);
        save(sprintf('~/ccContour%d.mat', ii));
    end
end

%% Plot color contour
