function [jndDist, acc, err] = coPixelVisibilityThreshold(d, params)
%% function coPixelVisibilityThreshold(params)
%    Compute classification accuracy of white display pixels vs uniform
%    scene on certain display and viewing distance
%
%  Inputs:
%    d     : isetbio display structure, see displayCreate for more detail
%    params: parameter structure, could include:
%      tDist         - vector of viewing distances to be tested
%      threshold     - threshold of accuracy to be marked as 'jnd'
%      OTHER         - See coPixelVisibilityAcc for detail
%
%  Outputs:
%    jndDist: just noticeable distance
%    acc:     classification accuray
%    err:     classification standard deviation
%
%  See also:
%    coPixelVisibilityAcc, s_coPixelVisibility_Proclus
%
% (HJ) ISETBIO TEAM, 2014

%% Check inputs & init
%  check inputs
if notDefined('d'), error('display required'); end
if notDefined('params'), params = []; end

%  init parameters
if ischar(d), d = displayCreate(d); end

try thresh = params.threshold; catch, thresh = 0.8; end
try tDist  = params.tDist;     catch, tDist  = 0.7:0.02:0.8; end

%% Compute classification accuracy for each distance
%  init parameters
acc = zeros(length(tDist), 1);
err = zeros(length(tDist), 1);

%  loop over tDist and compute
for ii = 1 : length(tDist)
    [acc(ii), err(ii)] = coPixelVisibilityAcc(d, tDist(ii), params);
    fprintf('Distance: %.2f\t Accuracy:%.2f\n', tDist(ii), acc(ii));
end

%% Interpolate and estimate merely visible distance
try
    tRange = min(tDist):((max(tDist)-min(tDist))/100):max(tDist);
    interpolatedAcc = interp1(tDist, acc, tRange, 'linear');
    [~, ind] = min(abs(interpolatedAcc - thresh));
    jndDist = tRange(ind);
catch
    jndDist = 0;
end

end