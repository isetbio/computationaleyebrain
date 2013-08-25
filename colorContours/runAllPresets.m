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

%% Run 'em
runBasicNoSurround;
runBasicDetermSurround;
runBasicNoSurroundWithNoise;
runBasicDetermSurroundWithNoise;
runBasicRDrawSurround;
runMacularPigmentVary;
