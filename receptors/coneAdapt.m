function [sensor, adaptedData] = coneAdapt(sensor, typeAdapt)
%% Cone adaptation
%
%   [sensor, adaptedData] = coneAdaptation(sensor, typeAdapt)
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
% the release decreases and the range available for another flash. If you
% darken from the mean, the rate can increase.
%
% The way to compute the mean must depend on a combination of the
% photoisomerizations and the recycling rate, probably through some
% equilibrium equation. The recycling rate depends on the isomerization
% rate and they set an equilibrium that is slower and slower as the mean
% when background gets brighter.
%
% In addition to Rodieck's discussion, there are famous papers by Boynton
% (e.g. adaptation in cones) expressing such a model based on cone ERPs.
%
% The other issue is the total gain.  The cones can only modulate over a
% dynamic range around the mean.  So, we set the gain as well to keep the
% modulation within some range, such as +/- 40 mV of the mean (deplorizing
% voltage).
%
% The third issue is whether all the cones are the same, or there is
% space-variant adaptation. We now here consider the difference between
% different cone types, but we think it's spatial invariant.
%
% Inputs:
%  sensor:     ISETBio sensor with cone absorption computed and stored
%  typeAdapt:  The adaptation model
%   typeAdapt= 0 - no adaptation
%   typeAdapt= 1 - a single gain map for the whole cone array
%   typeAdapt= 2 - one gain map computed for each cone class
%   typeAdapt= 3 - non-linear cone adapt for each type of cone
%   typeAdapt= 4 - cone by cone adaptation (NYI)
%
% Output:
%   sensor:        ISETBio sensor with cone adaptation parameters set. The
%                  parameters include gain and offset. You could retrieve
%                  the adaptation parameters and data by sensorGet
%                  function calls.
%   adaptedData:   Adapted voltage 3D matrix, would be in same size as
%                  volts image in sensor.
%
%
% ALERT:  The amount of adaptation is hard coded at this moment. Maybe we
%         should accept a third input parameter to set it.
%
% In the future, we will look at the time series and have a time-varying
% adaptation function. Generally, cones take 200ms and rods take 800ms for
% adaptation.
%
% Examples:
%   sensor = coneAdapt(sensor, 1);
%   adaptedData = sensorGet(sensor, 'adapted volts');
%   gain = sensorGet(sensor, 'adaptation gain');
%   offset = sensorGet(sensor, 'adaptation volts');
%
% (c) Stanford VISTA Lab, 2014

%% Check inputs and Init
if notDefined('sensor'),      error('sensor is required'); end
if ~sensorCheckHuman(sensor), error('unsupported species'); end
if notDefined('typeAdapt'),   typeAdapt = 2; end

%% Compute cone adaptations
volts  = double(sensorGet(sensor, 'volts'));

if isempty(volts), error('cone absorptions should be pre-computed'); end

vSwing = pixelGet(sensorGet(sensor,'pixel'),'voltageSwing');

switch typeAdapt
    case 0 % no adaptation
        gainMap = 1;
        offset  = 0;
    case 1
        % Use same gain for all cone type
        
        % Set the gain so that the max - min get scaled to 0.8 * vSwig mV
        gainMap = 0.8 * vSwing / (max(volts(:)) - min(volts(:)));
        adaptedData = gainMap * volts;
        
        % Set the zero level as the median. Actually, we could use mean,
        % but to avoid some extreme bright points, we use median
        offset  = median(adaptedData(:));
    case 2
        % Use different gains for each cone type
        % For human, the cone types in sensor are 1~4 for K, L, M, S
        gainMap = ones(4, 1);
        
        for ii = 2 : 4 % L,M,S and we don't need to compute for K
            v = sensorGet(sensor,'volts',ii);
            if ~isempty(v)
                gainMap(ii) = 0.8 * vSwing / (max(v) - min(v));
            end
        end
        
        nSamples = size(volts, 3);
        gainMap = gainMap(sensorGet(sensor, 'cone type'));
        
        adaptedData = volts .* repmat(gainMap, [1 1 nSamples]);
        % Set the zero level as the median
        offset      = median(adaptedData(:));
    case 3
        % Use non-linear cone adaptation for each cone type
        % Formula:
        %   R(Is|Ia) = R(p(Ia)Is|0) - R(p(Ia)Ia|0)
        %   R(Is|0)  = Is^n R_max / (Is^n + Kr^n)
        %   p(Ia)    = Kd / (Ia + Kd)
        % For human:
        %   n = 0.7; R_max = 40 mV;
        %   Kr = 166; Kd = 194; (units: trolands)
        % See Dawis (1982) paper for more details
        
        % Question: Does it make sense to use the same Kr and Kd for all
        % the three cone types?
        
        % From the idea of chromatic adaptation, we should scale the three
        % channels differently to make the adapted white point "white". If
        % we took median of cone absorptions as white point and as
        % adaptation offset, then we should be fine with chromatic
        % adaptation.
%         n = 0.7; R_max = 0.04;
        
        % Compute white point for L,M,S
%         adaptPoint  = zeros(4,1);
%         for ii = 2 : 4 % ii = 1 is for K, which is ignored
%             v = double(sensorGet(sensor, 'volts', ii));
%             adaptPoint(ii) = median(v);
%         end
%         
        % Convert unit for Kr and Kd
        % Kr = 166; Kd = 194; in units of trolands
        % Here, the conversion is just an approximation.
%         pixel = sensorGet(sensor, 'pixel');
%         photon2volts = pixelGet(pixel, 'conversion gain');         
%         Kr = 16.6 * photon2volts; Kd = 19.4 * photon2volts;
%         
%         coneType = sensorGet(sensor, 'cone type');
%         adaptPoint = adaptPoint(coneType);
%         adaptPoint = repmat(adaptPoint, [1 1 size(volts,3)]);
%         
%         % Compute equivalent gain by R(p(Ia)Is|0)
%         gainMap = volts .* Kd ./ (adaptPoint + Kd); % compute p(Ia)Is
%         gainMap = gainMap .^ (n-1) * R_max ./ (gainMap.^n + Kr.^n);
%         
%         % Compute offset map by R(p(Ia)Ia|0)
%         offset = adaptPoint .* Kd ./ (adaptPoint + Kd);
%         offset = offset.^n * R_max ./ (offset.^n + Kr.^n);
%         
%         % Compute adapted data
%         adaptedData = volts .* gainMap;
%         adaptedData(isnan(adaptedData)) = 0;
        
        % I'm not sure if this is right. Just throw out an error first
        error('NYI');
    otherwise
        error('unknown adaptation type');
end

% Compute adapted data
adaptedData = adaptedData - offset;

% Set adaptation parameters back to sensor
sensor = sensorSet(sensor, 'adaptation gain', gainMap);
sensor = sensorSet(sensor, 'adaptation offset', offset);

end
