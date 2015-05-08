% Script to analyze the visibility of compression artifacts
%
% The main script is compressionVisiblity.
% This script is to teach BW what is going on.
%
% HJ/BW Vistasoft

%% How to download an image
webdir = 'http://scarlet.stanford.edu/validation/SCIEN/ISETBIO/VESA/img';

% The ones ending in x0 are uncompressed
% urlwrite(fullfile(webdir,'t_FineTextRendering14_x0.bmp'),'original.bmp');
imgO = im2double(imread(fullfile(webdir,'t_FineTextRendering14_x0.bmp')));
% vcNewGraphWin; imshow(imgO)

% There are various types of compression
% urlwrite(fullfile(webdir,'t_FineTextRendering14_x9.bmp'),'compressed.bmp');
imgC = im2double(imread(fullfile(webdir,'t_FineTextRendering14_x9.bmp')));
% vcNewGraphWin; imshow(imgC)
imgSz = [size(imgO, 1), size(imgO,2)];

%% Set display and viewing parameters
% Model the display
d = displayCreate('LCD-Apple');

% Experiment conditions
% viewing distance
vd = 1;  % Meter
d = displaySet(d, 'viewing distance', vd);

% eye position and spatial integration size
pSzDeg  = 0.15; % spatial integration size in degree
pSzDots = pSzDeg * displayGet(d, 'dots per deg'); % convert to pixels

nTrial = 5;

% upper left position of region of interest (eye pos)
pos = floor(bsxfun(@times, rand(nTrial, 2), imgSz - pSzDots));

%% Compute classification accuracy
acc = zeros(nTrial, 1);  % nfold accuracy
pRange = 1:pSzDots;      % Crop region

for ii = 1 : nTrial
    % crop image to the region of interest
    % Let's use imcrop
    pO = imgO(pos(ii, 1) + pRange, pos(ii,2)+pRange, :);
    pC = imgC(pos(ii, 1) + pRange, pos(ii,2)+pRange, :);
    
    % compute acc for the cropped image, in display RGB.
    acc(ii) = compressionVisibility(pO, pC, d);
    
    % print log
    fprintf('Trial: %d\t Accuracy: %.2f\n', ii, acc(ii));
end

%% End