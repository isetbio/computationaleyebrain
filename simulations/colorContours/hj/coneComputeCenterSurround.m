function coneResp = coneComputeCenterSurround(coneResp, wc, kc, ks, params)
%% function coneComputeSSNoise(coneResp, [wc], [kc], [ks], [params])
%    Compute second site noise based on the cone absorptions
%    This function implements basic center surround idea of retina ganglion
%    cell model used in Bradley et. al paper: Retina-V1 model of
%    detectability across the visual field.
%    Formula for kernel:
%      D = wc * Gc - (1 - wc) * Gs
%
%  Inputs:
%    coneResp - matrix of cone absorption with no second site noise
%    wc       - weights for center mechanism
%    kc       - standard deviation for center Gaussian filter
%    ks       - standard deviation for surround Gaussian filter
%    params   - parameter structure, include
%             .addNoise   - logical, indicate whether or not add noise
%             .fanoFactor - fanoFactor for noise
%
%  Outputs:
%    coneResp      - matrix of cone response with second site noise added
%
%  See also:
%    coneComputeSSNoise
%
%  HJ ISETBIO Team, 2015

%% Check inputs
if notDefined('coneResp'), error('Cone response input required'); end
if notDefined('wc'), wc = 0.53; end
if notDefined('kc'), kc = 1; end
if notDefined('ks'), ks = 10.1; end
if notDefined('params'), params = []; end

try addNoise = params.addNoise; catch, addNoise = false; end
try fanoFactor = params.fanoFactor; catch, fanoFactor = 4; end

%% Generate Gaussian filters for Center-Surround
fc = fspecial('gaussian', round(3*ks), kc); % center filter
fs = fspecial('gaussian', round(3*ks), ks); % surround filter

%% Filter cone response
coneResp = wc * imfilter(coneResp, fc) - (1-wc) * imfilter(coneResp, fs);

%% Add noise
if addNoise
    stdMat = sqrt(fanoFactor * abs(coneResp));
    coneResp = coneResp + stdMat.*randn(size(coneResp));
end