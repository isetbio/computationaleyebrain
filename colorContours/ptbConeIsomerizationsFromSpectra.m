function [isoPerCone,pupilDiamMm,photoreceptors,irradianceWattsPerM2] = ptbConeIsomerizationsFromSpectralRadiance(spd_input,wls_input,pupilDiamMm,focalLengthMm,integrationTimeSec)
% [isoPerCone,pupilDiamMm,photoreceptors] = ptbConeIsomerizationsFromSpectralRadiance(spd_input,wls_input,pupilDiamMm,focalLengthMm,integrationTimeSec)
%
% Compute LMS human cone isomerizations from spectral radiance in Watts/[m2-sr-nm].
%
% This routine is set up for a quick commparison to isetbio calculations.  The underlying
% code is demonstrated and (sort of) documented in PTB routine IsomerizationsInEyeDemo.
%
% 8/4/13  dhb  Wrote it.

%% Set up PTB photoreceptors structure
% 
% We'll do the computations at the wavelength
% spacing passed in for the spectrum of interest.
whatCalc = 'LivingHumanFovea';
photoreceptors = DefaultPhotoreceptors(whatCalc);
photoreceptors.eyeLengthMM.source = num2str(focalLengthMm);
photoreceptors.nomogram.S = WlsToS(wls_input);
S = photoreceptors.nomogram.S;
photoreceptors = FillInPhotoreceptors(photoreceptors);

% Convert units to power per wlband rather than power per nm.
% Units of power per nm is the PTB way, for better or worse.
radianceWattsPerM2Sr = spd_input*S(2);
		
% Find pupil area, needed to get retinal irradiance, if not passed.
%
% In that case, pupil area based on the luminance of stimulus according
% to the algorithm specified in the photoreceptors structure.
load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);
theXYZ = T_xyz*radianceWattsPerM2Sr; theLuminance = theXYZ(2);
if (nargin < 3 || isempty(pupilDiamMm))
    [pupilDiamMm,pupilAreaMm2] = PupilDiameterFromLum(theLuminance,photoreceptors.pupilDiameter.source);
else
    pupilAreaMm2 = pi*((pupilDiamMm/2)^2);
end

% Convert radiance of source to retinal irradiance
irradianceWattsPerUm2 = RadianceToRetIrradiance(radianceWattsPerM2Sr,S, ...
    pupilAreaMm2,photoreceptors.eyeLengthMM.value);

% Pass back to calling routine in areal units of M2 and
% spectral units of 'per nm'.
irradianceWattsPerM2 = 1e12*irradianceWattsPerUm2/S(2);

%% Do the work in toolbox function
[isoPerConeSec,absPerConeSec,photoreceptors] = ...
    RetIrradianceToIsoRecSec(irradianceWattsPerUm2,S,photoreceptors);
isoPerCone = isoPerConeSec*integrationTimeSec;