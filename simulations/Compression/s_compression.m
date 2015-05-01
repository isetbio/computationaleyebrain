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


%% Set display and viewing parameters

% Model the display
d = displayCreate('LCD-Apple');

% viewing distance and other conditions could be set here.
sceneO = sceneFromFile(imgO,'rgb',[],d);
vcAddObject(sceneO); sceneWindow;

sceneC = sceneFromFile(imgC,'rgb',[],d);
vcAddObject(sceneC); sceneWindow;



%%  Model the human observer

% Choose a model of the optics

% Decide on eye movements

% render cone mosaics 

%% Run classifiers on them.


%% End