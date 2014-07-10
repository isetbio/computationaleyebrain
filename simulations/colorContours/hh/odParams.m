function inertP = odParams(visualfield, wave) 
% Create optical density default parameters
%
%
% Inputs:
%   visualfield:  'fovea' or 'periphery' or specific subject
%   wave:          Wavelength list
%
%
% see also ieReadHumanQE.m
% 
%  (c) VISTA lab 2012 HH
%%
if notDefined('wave'), wave = (400:10:700)'; end
if notDefined('visualfield'), visualfield = 'fovea'; end

inertP.visfield = visualfield;
inertP.wave     = wave;

switch lower(visualfield)
    
    % stockman foveal params (at 2 deg)
    case {'f','fov','fovea','stf','stockmanfovea'}
        
        inertP.lensTransF = 1;
        inertP.macDensity = 0.28;
        inertP.LPOD    = 0.5;
        inertP.MPOD    = 0.5;
        inertP.SPOD    = 0.4;
        inertP.melPOD  = 0.5;
        
    % stockman peripheral params (at 10 deg)    
    case {'p','peri','periphery','stp','stockmanperi','stockmanperiphery'}
        
        inertP.lensTransF = 1;
        inertP.macDensity = 0;
        inertP.LPOD    = 0.38;
        inertP.MPOD    = 0.38;
        inertP.SPOD    = 0.3;
        inertP.melPOD  = 0.5;
        
    case {'s1f'}
        
        inertP.lensTransF = 0.7467;
        inertP.macDensity = 0.6910;
        inertP.LPOD    = 0.4964; 
        inertP.MPOD    = 0.2250;
        inertP.SPOD    = 0.1480;
        inertP.melPOD  = 0.3239;
        inertP.visfield = 'f';
        
    case {'s1p'}
        
        inertP.lensTransF = 0.7467;
        inertP.macDensity = 0;
        inertP.LPOD    = 0.4964 ./ 0.5 .* 0.38; 
        inertP.MPOD    = 0.2250 ./ 0.5 .* 0.38;
        inertP.SPOD    = 0.1480 ./ 0.4 .* 0.3;
        inertP.melPOD  = 0.3239;
        inertP.visfield = 'p';
        
    case {'s2f'}
        
        inertP.lensTransF = 0.7637;
        inertP.macDensity = 0.5216;
        inertP.LPOD    = 0.4841; 
        inertP.MPOD    = 0.2796;
        inertP.SPOD    = 0.2072;
        inertP.melPOD  = 0.3549;
        inertP.visfield = 'f';

    case {'s2p'}
        inertP.lensTransF = 0.7637;
        inertP.macDensity = 0;
        inertP.LPOD    = 0.4841 ./ 0.5 .* 0.38; 
        inertP.MPOD    = 0.2796 ./ 0.5 .* 0.38;
        inertP.SPOD    = 0.2072 ./ 0.4 .* 0.3;
        inertP.melPOD  = 0.3549;
        inertP.visfield = 'p';
end
        