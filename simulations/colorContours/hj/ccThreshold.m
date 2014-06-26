function [threshold, expData] = ccThreshold(rContrast, direction, params)
%% Compute threshold in certain direction
%    function ccThreshold(rColor, direction, params)
%
%  Inputs:
%    rContrast - reference contrast in LMS
%    direction - angle from reference to test color in degrees
%    params    - parameter structure, could contain:
%                .pCorrect - percent of correctness for threshold
%                .bgColor  - background color
%                .ccParams - parameter structure for ccAcc function
%  
%  Outputs:
%    threshold - threshold found, if failed, NAN is returned
%    expData   - structure containing experiment information, including:
%                .rContrast  - 1x3 vector, reference contrast
%                .bgColor    - 1x3 vector, background RGB
%                .display    - name of display used for experiment
%                .mContrast  - Nx3 vector, test / match color
%                .acc        - Nx1 vector, classification accuracy
%                .err        - Nx1 vector, standard deviation of .acc
%                .pCorrect   - targe percent of correctness at threshold
%
%  Examples:
%    rContrast = [0 0 0];
%    direction = 45;
%    [thresh, expData] = ccThreshold(rContrast, direction);
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
if isfield(params.ccParams, 'd'), d = params.ccParams.d;
else d = 'OLED-Sony'; end
expData.display = d;

% Init background color
if isfield(params, 'bgColor'), bgColor = params.bgColor;
else bgColor = [0.5 0.5 0.5]; end
expData.bgColor = bgColor;

% Init reference color in RGB
rColor = coneContrast2RGB(d, rContrast, bgColor);

% Init direction
direction = [cosd(direction) sind(direction) 0];

% Init mContrast, acc, err
expData.acc = []; expData.err = []; expData.mContrast = [];

% Init pCorrect
if isfield(params, 'pCorrect'), pCorrect = params.pCorrect;
else pCorrect = 0.8; end

%% Binary search
lDist = 0; lAcc = 0.5;

% Estimate upper bound
uDist = 0.1; uAcc = 0;
while uAcc < pCorrect
    mContrast = rContrast + uDist * direction;
    mColor = coneContrast2RGB(d, mContrast, bgColor);
    uAcc  = ccAcc(rColor, mColor, params.ccParams);
    uDist = uDist + 0.1;
    assert(uDist < 0.6, 'Cannot classify even with obvious difference');
end

% Binary search
while true
    curDist = (lDist + uDist) / 2;
    mContrast = rContrast + curDist * direction;
    mColor = coneContrast2RGB(d, mContrast, bgColor);
    [curAcc, curErr] = ccAcc(rColor, mColor, params.ccParams);
    if curAcc > pCorrect
        uDist = curDist;
        uAcc = curAcc;
    else
        lDist = curDist;
        lAcc = curAcc;
    end
    expData.acc = cat(1, expData.acc, curAcc);
    expData.err = cat(1, expData.err, curErr);
    expData.mContrast = cat(1, expData.mContrast, mContrast(:)');
    if uAcc < pCorrect + 0.05 && lAcc > pCorrect - 0.05
        break;
    end
    
    fprintf('curDist:%f\t curAcc:%f\n', curDist, curAcc);
end

% Interpolate and find threshold
threshold = lDist + (uDist - lDist)*(pCorrect - lAcc)/(uAcc - lAcc);

end