% runAllPresets
%
% Run color contours with all the presets

%% Clear out the junk.  Remember where you are.
%
% Sometime, but not always s_initISET clears debugger stop points.
% So I (DHB) commented it out pending a better understanding.
clear; close all;  %s_initISET
talkD = pwd;
saveFlag = 0;

%% Make sure we are in the right place
cd(fileparts(mfilename('fullpath'))); %#ok<MCCD>

%% Cell array of what to run
thePresets = {'NoSurroundNoSecondSiteNoise' 'NoSurroundWithSecondSiteNoise' ...
    'RandomSurroundNoSecondSiteNoise' 'RandomSurroundWithSecondSiteNoise' ...
    'SelectiveSurroundNoSecondSiteNoise' 'SelectiveSurroundWithSecondSiteNoise' ...
    'MacularPigmentVary'};

%thePresets = {'NoSurroundNoSecondSiteNoise'};

%% Run 'em
for p = 1:length(thePresets)
    colorContours(thePresets{p});
end


