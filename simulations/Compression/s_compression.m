% Script to analyze the visibility of compression artifacts
%
% The main script is compressionVisiblity.
%
% HJ/BW Vistasoft

%% Load Image List
%  Web directory for image data set
webdir = 'http://scarlet.stanford.edu/validation/SCIEN/ISETBIO/VESA';

%  Get list for full size images
imgList = lsScarlet([webdir '/Full_1080p_ref_images'], '.bmp');

%% Set display and viewing parameters
% Model the display
d = displayCreate('LCD-Apple');

% Experiment conditions
% viewing distance
vd = 1;  % Meter
d = displaySet(d, 'viewing distance', vd);

% eye position and spatial integration size
pSzDeg  = 0.12; % spatial integration size in degree
pSzDots = round(pSzDeg * displayGet(d, 'dots per deg')); % convert to dots

nTrial = 10;

% Init for accuracy
acc = zeros(length(imgList), 16);

for ii = 1 : length(imgList)
    % get image name with out '.bmp'
    imgName = [webdir '/cropped_images/' imgList(ii).name(1:end-4)];
    
    % load cropped image - original
    imgO = im2double(imread([imgName '_x0.bmp']));
    imgSz = [size(imgO, 1), size(imgO,2)];

    % There are various types of compression labeled by x1 ~ x16
    for jj = 1 : 16
        % try loading compressed images
        try imgC = im2double(imread(sprintf([imgName '_x%d.bmp'], jj)));
        catch
            warning('Image missing:%s_x%d.bmp', imgName, jj);
        end
        
        % print log
        fprintf('%s_x%d:', imgList(ii).name(1:end-4), jj);
        
        % randomize eye position
        pos = floor(bsxfun(@times, rand(nTrial, 2), imgSz - pSzDots)) + 1;
        
        % compute accuracy for each trial
        curAcc = zeros(nTrial, 1);
        for kk = 1 : nTrial
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
    end
end

%% End