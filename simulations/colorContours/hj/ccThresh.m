function [threshold, params] = ccThresh(rContrast, direction, params)
%% Compute threshold in certain direction
%    function ccThreshold(rColor, direction, params)
%
%  Inputs:
%    rContrast - reference contrast in LMS
%    direction - angle from reference to test color in degrees (LM plane)
%    params    - parameter structure, could contain:
%                .rseed    - random seed for the experiment
%                .pCorrect - percent of correctness for threshold
%                .bgColor  - background color
%                .ccParams - parameter structure for ccAcc function
%  
%  Outputs:
%    threshold - threshold found, if failed, NaN is returned
%    params    - parameters used for the simulations, includes
%                .rseed    - random seed for the experiment
%                .pCorrect - percent of correctness for threshold
%                .bgColor  - background color
%                .ccParams - parameter structure for ccAcc function
%                .expData  - experiment data structure, include:
%                  .rContrast  - 1x3 vector, reference contrast
%                  .direction  - angle from reference to test color
%                  .bgColor    - 1x3 vector, background RGB
%                  .display    - name of display used for experiment
%                  .mContrast  - Nx3 vector, test / match color
%                  .acc        - Nx1 vector, classification accuracy
%                  .err        - Nx1 vector, standard deviation of .acc
%                  .pCorrect   - targe percent of correctness at threshold
%
%  Examples:
%    rContrast = [0 0 0];
%    direction = 45;
%    [thresh, params] = ccThreshold(rContrast, direction);
%
%  See also:
%    ccAcc
%
%  (HJ) ISETBIO TEAM, 2014

%% Init
%  Check inputs
if notDefined('rContrast'), error('reference contrast required'); end
if notDefined('direction'), error('direction required'); end
if notDefined('params'), params = []; end

% Init display
expData.rContrast = rContrast(:)';
if ~isfield(params, 'ccParams'), params.ccParams = []; end
if isfield(params.ccParams, 'd')
    d = params.ccParams.d;
else
    d = displayCreate('OLED-Sony');
    d = displaySet(d, 'gamma', 'linear');
    params.ccParams.d = d;
end
expData.display = d;

% Init background color
if isfield(params, 'bgColor'), bgColor = params.bgColor;
else bgColor = [0.5 0.5 0.5]; params.bgColor = bgColor; end
expData.bgColor = bgColor;

% Init reference color in RGB
rColor = coneContrast2RGB(d, rContrast, bgColor);

% Init direction
expData.diretion = direction;
direction = [cosd(direction) sind(direction) 0];

% Init mContrast, acc, err
expData.acc = []; expData.err = []; expData.mContrast = [];

% Init pCorrect
if isfield(params, 'pCorrect'), pCorrect = params.pCorrect;
else pCorrect = 0.8; params.pCorrect = pCorrect; end

%% Binary search
lDist = 0; lAcc = 0.5;
s = sensorCreate('human');
% Estimate upper bound
uDist = 0.1; uAcc = 0;
while uAcc < pCorrect
    mContrast = rContrast + uDist * direction;
    mColor = coneContrast2RGB(d, mContrast, bgColor, s);
    uAcc   = ccAcc(rColor, mColor, params.ccParams);
    uDist  = uDist + 0.1;
    assert(uDist < 0.3, 'Cannot classify even with obvious difference');
end

% Binary search
while true
    curDist = (lDist + uDist) / 2;
    mContrast = rContrast + curDist * direction;
    mColor = coneContrast2RGB(d, mContrast, bgColor, s);
    [curAcc, curErr, ccParams] = ccAcc(rColor, mColor, params.ccParams);
    if curAcc > pCorrect
        uDist = curDist;
        uAcc = curAcc;
    else
        lDist = curDist;
        lAcc = curAcc;
    end
    
    % print debug info
    fprintf('\t\tContrast: %f\t Accuracy:%f\n', curDist, curAcc);
    
    expData.acc = cat(1, expData.acc, curAcc);
    expData.err = cat(1, expData.err, curErr);
    expData.mContrast = cat(1, expData.mContrast, mContrast(:)');
    if uAcc < pCorrect + 0.05 && lAcc > pCorrect - 0.05
        break;
    end
    
    if uAcc < pCorrect + 0.01
        threshold = uDist;
        return;
    elseif lAcc > pCorrect - 0.01
        threshold = lDist;
        return;
    end
end

% Interpolate and find threshold
threshold = lDist + (uDist - lDist)*(pCorrect - lAcc)/(uAcc - lAcc);

% set expData and ccParams to params
params.expData = expData;
params.ccParams = ccParams;

end