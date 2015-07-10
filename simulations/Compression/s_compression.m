% Script to analyze the visibility of compression artifacts
%
% The main script is compressionVisiblity.
%
% HJ/BW Vistasoft

%% Load Image List

% We might use the rdata command in isetbio
% remote.host = 'http://scarlet.stanford.edu/validation/SCIEN/ISETBIO/VESA';
% remote.directory = 'Full_1080p_ref_images';
% imgList = rdata('ls',remote,'.bmp');

%  Web directory for image data set
webdir = 'http://scarlet.stanford.edu/validation/SCIEN/ISETBIO/VESA';

%  Get list for full size images
imgList = lsScarlet([webdir '/Full_1080p_ref_images'], '.bmp');

%% Set display and viewing parameters

% Model the display
d = displayCreate('LCD-Apple');

% Adjust viewing distance to match pixel per degree
ppd = 60;                                  % 60 pixels per degree
vd  = ppd / displayGet(d, 'dpi') / tand(1) * 0.0254;  % Meter
d   = displaySet(d, 'viewing distance', vd);
invGamma = displayGet(d, 'inverse gamma'); % inverse gamma table
rgb2xyz  = displayGet(d, 'rgb2xyz');       % rgb to xyz tranformation matrix
wp = displayGet(d, 'white xyz');           % white point

% Eye position and spatial integration size
pSzDeg  = 0.08;               % spatial integration size in degree
% Represented here in terms of dots on the display
pSzDots = round(pSzDeg * displayGet(d, 'dots per deg')); 
fprintf('Spatial integration %.3f arc minutes (%.3f deg)\n',60*pSzDeg, pSzDeg);
conesPerDeg = 100;  
fprintf('This spans a region containing %.1f cones\n', (conesPerDeg*pSzDeg)^2);

% Number of trials
nTrial = 15;

% Initialize parameters that will summarize accuracy
acc = zeros(length(imgList), 16);           % in shape of n_images, n_algorithms
expAcc = nan(length(imgList), 16);
cieDeltaE  = zeros(length(imgList), 16);
scieDeltaE = zeros(length(imgList), 16);

%% Load experimental data 

%  The behavioral data are stored in the github folder
expData = importdata('2014_VESA.csv');
expData.textdata = expData.textdata(2:end, 1);

%% Compute accuracy for algorithms and images

for ii = 1 : length(imgList)
    % get image name but strip off the extension '.bmp'
    imgName = [webdir '/cropped_images/' imgList(ii).name(1:end-4)];
    
    % load cropped image - original
    imgO = im2double(imread([imgName '_x0.bmp']));
    imgSz = [size(imgO, 1), size(imgO,2)];

    % There are various types of compression labeled by x1 ~ x16
    for jj = 1 : 16
        % try loading compressed images
        try 
            imgC = im2double(imread(sprintf([imgName '_x%d.bmp'], jj)));
        catch
            warning('Image missing:%s_x%d.bmp', imgName, jj);
        end
        
        % print log
        fprintf('%s_x%d:', imgList(ii).name(1:end-4), jj);
        
        % randomize eye position
        % Can we use the eye movement code in ISETBIO for this?
        pos = floor(bsxfun(@times, rand(nTrial, 2), imgSz - pSzDots)) + 1;
        
        % compute accuracy for each trial
        curAcc = zeros(nTrial, 1);
        for kk = 1 : nTrial
            
            % Crop out a portion of the original and compressed image. The
            % cropped portion corresponds to the eye positions at pos.
            pO = imcrop(imgO, [pos(kk, [2 1]) pSzDots pSzDots]);
            pC = imcrop(imgC, [pos(kk, [2 1]) pSzDots pSzDots]);
            
            % compute accuracy for the cropped image, in display RGB.
            curAcc(ii) = compressionVisibility(pO, pC, d);
        end
        
        % compute probability sum
        curAcc(curAcc < 0.5) = 0.5;
        acc(ii, jj) = 1 - prod(2 - 2*curAcc)/2;
        
        % print log
        fprintf('%.2f\n', acc(ii, jj));
        
        % try load experiment data
        try
            expImg = sprintf('%s_x%d.bmp', imgList(ii).name(1:end-4), jj);
            indx = cellfun(@(s) ~isempty(strfind(s, expImg)), expData.textdata);
            expAcc(ii, jj) = expData.data(indx);
        catch
        end
        
        % compute DeltaEab for CIELab
        imgO_Linear = ieLUTLinear(imgO, invGamma);
        imgC_Linear = ieLUTLinear(imgC, invGamma);
        imgO_XYZ = imageLinearTransform(imgO_Linear, rgb2xyz);
        imgC_XYZ = imageLinearTransform(imgC_Linear, rgb2xyz);
        
        % Keep this out for a while.  We wil put it back after we move
        % deltaEab into ISETBIO
        % de = deltaEab(imgO_XYZ, imgC_XYZ, wp);
        % cieDeltaE(ii, jj) = median(de(:));
        
        % compue DeltaEab for SCIELab
        % params.sampPerDeg = ppd;
        % params.imageFormat = 'xyz10';
        % params.filters = []; 
        % params.deltaEversion = '1976';
        % de = scielab(imgO_XYZ, imgC_XYZ, wp, params);
        % scieDeltaE(ii, jj) = median(de(:));
    end
end

%% Visualize
vcNewGraphWin; 
% acc and expAcc are different different images and compression levels
% acc(img,clevel)
plot(acc(:), expAcc(:), 'o');
xlabel('Predicted Accuracy'); ylabel('Measurement Accuracy');

save result.mat acc expAcc cieDeltaE scieDeltaE

%% Analyze how the predicted and measured difficulties correlate

%  See rank correlation
indx = ~isnan(expAcc);
exp_diff = bsxfun(@minus, expAcc(indx), expAcc(indx)');
sim_diff = bsxfun(@minus, acc(indx), acc(indx)');
rank_corr = sum(exp_diff(:) .* sim_diff(:) >= 0) / numel(exp_diff);

% See rank correlation
indx = abs(exp_diff) > 0.05;
rank_corr_large = sum(exp_diff(indx) .* sim_diff(indx) >= 0)/sum(indx(:));

% Mean and variation accuracy of each image averaged across compression
% levels.
sim_imgMean = mean(acc, 2); sim_imgStd = std(acc, [], 2);
exp_imgMean = nanmean(expAcc, 2); exp_imgStd = nanstd(expAcc, [], 2);
indx = ~isnan(exp_imgMean);

vcNewGraphWin;
plot(sim_imgMean(:),exp_imgMean(:),'o')
grid on;
xlabel('Simulated'); ylabel('Experiment')
axis equal; identityLine;

vcNewGraphWin; hold on;
X = 1:sum(indx); X = X(:);
bar(X, [sim_imgMean(indx) exp_imgMean(indx)]);
errorbar([X-0.15 X+0.15], [sim_imgMean(indx) exp_imgMean(indx)], ...
    [sim_imgStd(indx) exp_imgStd(indx)], '.');
grid on;
ax = gca; ax.XTick = X; ax.XTickLabel = {imgList(indx).name};
ax.XTickLabelRotation = 45;
legend('ISETBIO Pred', 'Measurement Data');

% Mean and variation of each algorithm averaged across images
sim_algMean = mean(acc); sim_algStd = std(acc);
exp_algMean = nanmean(expAcc); exp_algStd = nanstd(expAcc);
indx = ~isnan(exp_algMean);

vcNewGraphWin;
title('Algorithm comparisons')
plot(sim_algMean(:),exp_algMean(:),'o')
grid on;
xlabel('Simulated'); ylabel('Experiment')
axis equal; identityLine;

vcNewGraphWin; hold on;
X = 1:sum(indx); X = X(:);
bar(X, [sim_algMean(indx)' exp_algMean(indx)']);
errorbar([X-0.15 X+0.15], [sim_algMean(indx)' exp_algMean(indx)'], ...
    [sim_algStd(indx)' exp_algStd(indx)'], '.');
grid on;

%% End