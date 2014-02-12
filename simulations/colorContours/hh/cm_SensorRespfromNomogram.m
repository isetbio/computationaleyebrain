function Sensor = cm_SensorRespfromNomogram(peaksenstivity, modelparams, sensor_axial_density, wavelength, methods)
%
% Sensor = cm_SensorRespfromNomogram(peaksenstivity, modelparams, sensor_axial_density, wavelength, methods)
%
% Using nomogram, the code estimates 'Sensor' response.
%
% % Example - let's estimate ipRGC response function
% peaksenstivity       = 482;
% modelparams          = [1 0];
% sensor_axial_density = 0.5;
% wavelength           = cm_getDefaultWls;
% Sensor = cm_SensorRespfromNomogram(peaksenstivity, modelparams, sensor_axial_density, wavelength);
% figure, plot(wavelength,Sensor)
%
% see also cm_prepCalibrationMtx.m
%
% HH (c) Vista lab 2012. 
% 
%% prep
if ieNotDefined('peaksenstivity'), help cm_SensorRespfromNomogram; return ; end
if ieNotDefined('modelparams'), modelparams = [1 0]; end % peripheral condition
if ieNotDefined('wavelength'),  wavelength  = cm_getDefaultWls; end
if ieNotDefined('sensor_axial_density'), sensor_axial_density = 0.5; end % tentative density for ipRGC

if ieNotDefined('methods'),
    methods.Lenssource = 'stockman2';
    methods.LMSsource  = 'StockmanSharpe';
end

%% sensor responsivity using by A1 Based photopigment nomogram
%  with transmittance of cornea, lense and macular pigments

% Calculate the whole eye transmittance, combining macular, lens, and other
% factors
lensTransmitF  = modelparams(1);
maculardensity = modelparams(2); % macular density, which varies across indivisual

% Note that the macular density has a huge effect. Below we plot the
% melanopsin with a range of macular densities.
mt = macular(maculardensity, wavelength);
lt = cm_LensTransmittance(lensTransmitF, wavelength, methods.Lenssource);

% eye transmittance (et)
et = mt.transmittance .* lt';

sensorAbosorbance = PhotopigmentNomogram(wavelength',peaksenstivity,methods.LMSsource);
sensorAbsorbtance = AbsorbanceToAbsorbtance(sensorAbosorbance, wavelength, sensor_axial_density); % 3 x w

% nomogram combined with complete eye transmittance (et)
Sensor = sensorAbsorbtance' .* et;
Sensor = Sensor ./ max(Sensor);

return