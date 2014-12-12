function [threshold, expData] = ccThresh(rContrast, direction, params)
%% Compute threshold in certain direction
%
%    [threshold, expData] = ccThresh(rColor, direction, params)
%
%  Inputs:
%    rContrast - reference contrast in LMS
%    direction - angle in LM space, L = cos(d)*s and M = sin(d)*s
%    params    - parameter structure, could contain:
%                .pCorrect - percent correct at threshold
%                .bgColor  - background RGB 
%                     ** Should be mean L,M,S (Stockman) Absorptions **
%                .ccParams - parameter structure for ccAcc (accuracy)
%  
%  Outputs:
%    threshold - threshold contrast, NAN on failure 
%    expData   - contains experiment information for replication:
%                .rContrast  - 1x3 vector, reference contrast
%                .bgColor    - 1x3 vector, background RGB
%                .display    - name of display used for experiment
%                .mContrast  - Nx3 vector, test / match color
%                .acc        - Nx1 vector, classification accuracy
%                .err        - Nx1 vector, standard deviation of .acc
%                .pCorrect   - targe percent of correctness at threshold
%
%  The ccParams structure contains information about the experiment.  See
%  the function ccAcc.
%
%  Examples:
%    rContrast = [0 0 0];
%    direction = 45;
%    [thresh, expData] = ccThresh(rContrast, direction);
%
%  See also: ccAcc, s_colorContours
%
%  (HJ) ISETBIO TEAM, 2014

%% Programming TODO
%  Philosophy now will be to make this all an LMS function and use
%  auxiliary rgb2lms type functions for people who want to start at a
%  calibrated display and use this.
%
%  Change direction to be a unit 3-vector in LMS space or a general
%  3-vector.
%  Make background (L,M,S) mean.
%  To simplify this maybe we should write
%         scene = sceneCreate('LMS image',lmsData)
%
%  Let's make sure expData has everything organized so we can easily call
%  ccThresh again and replicate the results.  Means 
%    add direction slot
%    put parameters inside of a params struct so expData.params is useful
%    A little bit of restructuring.
%

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

% Dist is viewing distance, and Acc is classification accuracy.
lDist = 0; lAcc = 0.5;  % Lower bound on the distance and accuracy
uDist = 0.2; uAcc = 0;  % Upper bound on the distance and accuracy

% Adjust viewing distance so that accuracy is at least equal to pCorrect
while uAcc < pCorrect
    mContrast = rContrast + uDist * direction;
    mColor = coneContrast2RGB(d, mContrast, bgColor);
    uAcc  = ccAcc(rColor, mColor, params.ccParams);
    uDist = uDist + 0.1;
    assert(uDist < 0.6, 'Cannot classify even with obvious difference');
end

% Binary search for viewing distance that corresponds to pCorrect
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
    
    % print debug info
    fprintf('curDist:%f\t curAcc:%f\n', curDist, curAcc);
    
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

% Interpolate and return threshold as an average.
threshold = lDist + (uDist - lDist)*(pCorrect - lAcc)/(uAcc - lAcc);

end