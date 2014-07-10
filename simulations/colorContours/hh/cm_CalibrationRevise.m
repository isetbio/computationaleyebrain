function CalreviseMtx = cm_CalibrationRevise(currAxisPhotmeter, prevAxisPhotmeter, subjectname, POD, macularfactor, lensfactor, lambdashift)
%
% CalreviseMtx = cm_CalibrationRevise(currAxisPhotmeter,prevAxisPhotmeter,subjectname, subjectname, POD, macularfactor, lensfactor, lambdashift)
%
% A function of the code is to translocate 4x1 vector in one condition to
% another condition. 
%
% An aim of this function first was to revise an error of a calibration
% done by uncalibrated PR-650 because there were differences at a power of
% short andlong wavelength between pr650 and pr715.
%  
% Current version extends to have more functions
% 1. convert vectors with different a) photopigment optical densities.
%                                   b) macular pigments
%                                   c) parameterized lens pigments
%                                   d) peak sensitivities
%
%
% <Input>
%   currAxisPhotmeter ... a name of photometer for calculation as current
%   prevAxisPhotmeter ... a name of photometer for measurement
%   subjectname       ... char, if includes 'fov', calc foveal data.
%   POD               ... a) scalor, the ratio of photopigment optica densities compare to standard color observer
%                         b) 1x4 vector, each photopigment optical density 
%   macularfactor     ... scaler, macular pigment density
%   lensfactor        ... scaler, a factor of lens pigment density
%   lambdashift       ... shift of pigment function's peak sensitivities
%
% <Output>
%   CalreviseMtx      ... 4x4 matrix, calibration revised matrix
%                         v' = CalreviseMtx * v(4,1) 
%
% Example 1
%   CalreviseMtx = cm_CalibrationRevise
%
% Example 2
%   currAxisPhotmeter = 'pr715'; prevAxisPhotmeter = 'pr650';
%   CalreviseMtx = cm_CalibrationRevise(currAxisPhotmeter, prevAxisPhotmeter)
%
% Example 3
%   POD = [0.4 0.4 0.3 0.5];
%   CalreviseMtx = cm_CalibrationRevise([],[],'s1f',POD)
%
%
% see also cm_loadStairCaseData.m
%
% C) vista lab 2011, HH
%% prep params
if ~exist('currAxisPhotmeter','var') || isempty(currAxisPhotmeter)
    currAxisPhotmeter = 'pr715';
end

if ~exist('prevAxisPhotmeter','var') || isempty(prevAxisPhotmeter)
    prevAxisPhotmeter = 'pr715';
end

if ~exist('subjectname','var') || isempty(subjectname)
    subjectname = 'none';
end

if ~exist('POD','var') || isempty(POD)
    POD = 1;
end

if ~exist('lensfactor','var') || isempty(lensfactor)
    lensfactor = 1;
end

if ~exist('lambdashift','var') || isempty(lambdashift)
    lambdashift = 0;
end

if  ~isempty(strfind(subjectname,'f'))
    fovealflag = true;
else
    fovealflag = false;
end

if ~exist('macularfactor','var') || isempty(macularfactor)
    if fovealflag == true;
        macularfactor = 0.28;
    else
        macularfactor = 0;
    end
end

wls = cm_getDefaultWls;

%% calicurate pigment-isolate stimulu as a combination of 6 lights
%  for standard color observer at current photometers axis
[sensorIsolateAmp1P AbsSpectra densities] = cm_prepCalibrationMtx(fovealflag, prevAxisPhotmeter);

%% get matrix which comvert previous axis to current axis 
meanLEDAmps = 0.5 * ones(1,6);
 
% get mean led spd
meanLEDspd = cm_getledSPD(wls, currAxisPhotmeter);

% Photopigment optical densities
if length(POD) == 1
    PODdensities = densities * POD;
elseif length(POD) == 4
    PODdensities = POD;
end

% call macular function
tmp = macular(macularfactor, wls);

% calc transmittance based on macular and lens pigment densities
lensTransmit = cm_LensTransmittance(lensfactor,wls,'stockman2');
eyeTransmit = tmp.transmittance .* lensTransmit';

% calc cone fundamentals with different pigment densities
varLMSIfun = cm_variableLMSI_PODandLambda(AbsSpectra, PODdensities, lambdashift, eyeTransmit');

CalreviseMtx = cm_ledamp2percont(sensorIsolateAmp1P,varLMSIfun, meanLEDspd, meanLEDAmps) * 100;

fprintf('[%s]:Current  axis is  calibrated with %s. \n', mfilename, currAxisPhotmeter);
fprintf('[%s]:Previous axis was calibrated with %s. \n', mfilename, prevAxisPhotmeter);

