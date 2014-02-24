function [sensor, gain, offset] = coneAdapt(sensor, typeAdapt)
%% Cone adaptation
%
%   [sensor, gain, offset] = coneAdaptation(sensor, typeAdapt)
%
% Implement adaptation models to produce the cone volts. Cone absorption
% sampels should be computed and stored in sensorGet(sensor, volts).
%
% The mean level is determined as a reference so that increments and
% decrements are meaningful. This is needed so that on- and off-cells can
% work. When we have a cone absorption that is above the mean, this should
% cause an off-cell to reduce its firing rate and on-cell to increase.
% Conversely, cone absorptions below the mean cause an on-cell to decrease
% and an off-cell to increase.
%
% The discussion in Rodieck book - The first steps in seeing - p. 177 and
% thereabouts for ideas about cone signals. The physiological question is
% what is the meaning of the offset or zero mean. Rodieck describes the
% effect of light as setting the mean transmitter release. In the dark,
% there is a relatively large dynamic range. As the light is on steadily,
% the release decreases and the range available for another flash also
% decreases. If you darken from the mean, however, the rate can increase.
%
% The way in which the mean is set must depend on a combination of the
% photoisomerizations and the recycling rate, probably through some
% equilibrium equation.  The recycling rate depends on the isomerization
% rate and they set an equilibrium that is slower and slower as the mean
% background gets brighter.  This is how the value of the offset gets set.
%
% In addition to Rodieck's discussion, there are famous papers by Boynton
% (e.g. adaptation in cones) and others at SRI expressing such a model
% based on cone ERPs, and probably Rushton.
%
% The other issue is the total gain.  The cones can only modulate over a
% dynamic range around the mean.  So, we set the gain as well to keep the
% modulation within some range, such as +/- 50 mV of the mean.
%
% The third issue is whether all the cones are the same, or there is
% space-variant adaptation.  That is for the future.
%
% Inputs:
%  sensor:     ISETBio sensor with cone absorption computed and stored
%  typeAdapt:  The adaptation model
%   typeAdapt= 0 - no adaptation
%   typeAdapt= 1 - a single gain map for the whole cone array
%   typeAdapt= 2 - one gain map computed for each cone class
%   typeAdapt= 3 - cone by cone adaptation (NYI)
%
% Output:
%   sensor:   ISETBio sensor with cone adapted volts set
%
%
% ALERT:  The amount of adaptation is hard coded at this moment.  It should
% be a parameter somewhere.  Worse yet, the amount of adaptation differs
% between types 1 and 2. Fix this.
%
% In the future, we will look at the time series and have a time-varying
% adaptation function.
%
% Examples:
%   sensor = coneAdapt(sensor, 1);
%
% (c) Stanford VISTA lab, 2014

%% Why are there 4 adaptation terms?
%  When one is black, it shouldn't become adapated

%% Check inputs and Init
if notDefined('sensor'),      error('sensor is required'); end
if ~sensorCheckHuman(sensor), error('unsupported species'); end
if notDefined('typeAdapt'),   typeAdapt = 2; end

%% Compute cone adaptations
vSwing = pixelGet(sensorGet(sensor,'pixel'),'voltageSwing');
volts  = sensorGet(sensor, 'volts');

if isempty(volts), error('cone absorptions should be pre-computed'); end

switch typeAdapt
    case 0 % no adaptation 
        gain = 1;
    case 1
        % Method 1: same gain for all cone type
        
        % Set the gain so that the returned mean voltage is 1/2 of voltage
        % swing.  This adaptation has not dynamic, sigh.
        gainMap = (0.8*vSwing) / mean(volts(:));
        adaptedData = gainMap * volts;
        
    case 2 
        % Method 2 : different gain for each cone type
        % For human, the cone types in sensor are 1~4 for K, L, M, S
        gainMap = ones(4, 1);
        
        for ii = 2 : 4 % L,M,S and we don't need to compute for K
            v = sensorGet(sensor,'volts',ii);
            gainMap(ii) = (vSwing/2)/mean(v(:));
        end
        
        nSamples = size(volts, 3);
        gainMap = gainMap(sensorGet(sensor, 'cone type'));
        
        adaptedData = volts .* repmat(gainMap, [1 1 nSamples]);
    otherwise
        error('unknown adaptation type');
end

% Set the zero level as the median.  See discussion in the header
% It seems unreasonable to me, anyway.
offset      = median(adaptedData(:));
adaptedData = adaptedData - offset;

% Set adapted data back to sensor
sensor = sensorSet(sensor, 'volts', adaptedData);

end
