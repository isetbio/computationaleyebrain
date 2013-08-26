function [theVectors,oneConeEachClassStartIndices,isetLMSIsomerizations,isetLMSQuantalEfficiencies,isetLMSQuantalEfficiencyWavelengths] = ...
    getConeResponsesToStimRGB(cSensor,rgbVector,theParams,staticParams,runtimeParams,staticComputedValues)
% [theVectors,oneConeEachClassStartIndices,isetLMSIsomerizations,isetLMSQuantalEfficiencies,isetLMSQuantalEfficiencyWavelengths] = ...
%   getConeResponsesToStimRGB(cSensor,rgbVector,theParams,staticParams,runtimeParams,staticComputedValues)
%
% Get matrix of sensor values for use in classification.  Noise is added (or not) depending
% on the value of theParams.noiseType
%
% Eventually, we'll want to generalize this routine to handle more general stimulus
% and sensor values formats.
%
% Also returns some additional information that is primarily useful for debugging.
%
% 8/26/13  dhb  Pulled this out of main code, wrote so it will work for blank and test stimuli.

%% Filename unique to calling condition, to prevent stepping on ourselves from parfor loop
imageFileName = sprintf('tmp_%s_%d_%d_%d_%d_%d.png',theParams.OBSERVER_STATE,theParams.DO_TAFC_CLASSIFIER,round(100*theParams.macularPigmentDensityAdjust),...
    round(1000*theParams.cdAngle),round(1000*theParams.testLevel),now);

%% Create a scene with the background spectrum
%
% To patch into scene creation routines, we create an image
% on disk with spatially uniform RGB values.  These are then treated
% as frame buffer values, and will be raised to the specified
% gamma value(gamma uncorrected) by the scene creation routine.
%
% [**] Find a way to do this (and the corresponding version
% for the test) that does not involve writing an image to disk and
% that skips the gamma correction stuff.
ungammedRgbVector = rgbVector.^(1/staticParams.isetGammaValue);
tmp = ones(staticParams.scenePixels,staticParams.scenePixels,3);
[tmp,row,col] = RGB2XWFormat(tmp);
tmp = tmp*diag(ungammedRgbVector(:));
theImg = XW2RGBFormat(tmp,row,col);
clear tmp
imwrite(theImg,imageFileName,'png');
theScene = sceneFromFile(imageFileName,'rgb',[],[staticParams.monitorName '.mat'],staticComputedValues.wavelengthsNm);
theScene = sceneSet(theScene,'name','background');
theScene = sceneSet(theScene,'fov',staticParams.fieldOfViewDegrees);
unix(['rm ' imageFileName]);
%vcAddAndSelectObject(sceneB);
%sceneWindow;
        
%% Pass the background through the optics
theOiD = oiCompute(staticComputedValues.oiD, theScene);
vcAddAndSelectObject(theOiD);
%oiWindow;

%% Compute noise free background sensor image
theCSensorNF = sensorComputeNoiseFree(cSensor,theOiD);
%vcAddAndSelectObject(theCSensorNF);
%sensorWindow;

%% Get noise free cone isomerizations for background
%
% Store the max for each sensor class to pass back.
% We pick the max to avoid the edge effects in the sensor image.
for ii = 1:staticParams.nSensorClasses
    theSensorValsNF{ii} = sensorGet(theCSensorNF,'electrons',staticParams.isetSensorConeSlots(ii));
    isetLMSIsomerizations(ii) = round(max(theSensorValsNF{ii}));
end

% Get iset LMS quantal efficiences and wavelengths
temp = sensorGet(theCSensorNF,'spectralqe')';
isetLMSQuantalEfficiencies = temp(staticParams.isetSensorConeSlots,:);
isetLMSQuantalEfficiencyWavelengths = sensorGet(theCSensorNF,'wave');
clear temp

%% Get multivariate sample distribution of LMS
% responses out of the sensor objects.
%
% Have option of doing this noise free here, so we can add all
% our noise at an opponent site.  That isn't realistic, but
% I'm charging ahead just to try to get some intuitions about
% what opponency does.
if (theParams.noiseType ~= 0)
    theVoltImage = sensorComputeSamples(theCSensorNF,staticParams.nDrawsPerTestStimulus,theParams.noiseType);
else
    theVoltImageNF = sensorGet(theCSensorNF,'volts');
    theVoltImage = theVoltImageNF(:,:,ones(1,staticParams.nDrawsPerTestStimulus));
end

%% I'm sure there is a reason why we put volts back in and then
% take electrons out.  Brian probably explained it to me.  But
% now I forget.
% [** Brian, can you write a comment here?]
for k = 1:staticParams.nDrawsPerTestStimulus
    for ii = 1:staticParams.nSensorClasses
        theCSensorTemp = sensorSet(theCSensorNF,'volts',theVoltImage(:,:,k));
        theSensorVals{k,ii} = sensorGet(theCSensorTemp,'electrons',staticParams.isetSensorConeSlots(ii));
    end
end
clear theVoltImage

%% We want to control the integration area.
%
% Here we'll do this by specifying the fraction of
% the total mosaic to use.
%
% First figure out length of sample vector
for ii = 1:staticParams.nSensorClasses
    nUse(ii) = round(theParams.fractionUse*length(theSensorVals{1,ii}));
end
nUseAll = sum(nUse);

%% Now draw samples
theVectors = zeros(nUseAll,staticParams.nDrawsPerTestStimulus);
typeVec = zeros(nUseAll,1);
oneConeEachClassStartIndices = zeros(staticParams.nSensorClasses,1);
for k = 1:staticParams.nDrawsPerTestStimulus
    if (rem(k,10) == 0)
        if (~runtimeParams.SIM_QUIET)
            fprintf('\tGetting cone catches for draw %d of %d\n',k,staticParams.nDrawsPerTestStimulus);
        end
    end
    
    startIndex = 1;
    for ii = 1:staticParams.nSensorClasses
        % Pull out a set of randomly chosen responses for this sensor class
        % and tuck it into the response vector, for both background and test.
        % Also store the sensor class index for each stored response.
        endIndex = startIndex + nUse(ii) - 1;
        
        % Currently, this shufles the responses of each class before
        % pulling out the requisite number of responses, for some reason.
        % It shouldn't have any 
        %
        % [** Better will be to fix this up to be more realistic.]
        temp = Shuffle(theSensorVals{k,ii});
        theVectors(startIndex:endIndex,k) = temp(1:nUse(ii));
        if (k == 1)
            typeVec(startIndex:endIndex) = ii;
            oneConeEachClassStartIndices(ii) = startIndex;
        end
        
        startIndex = endIndex+1;
    end
end

