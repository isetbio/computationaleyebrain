function [isoPerCone,pupilDiamMm,photoreceptors] = ptbConeIsomerizationsFromSpectralRadiance(spd_isetbio,wls_isetbio,integrationTimeSec)
% [isoPerCone,pupilDiamMm,photoreceptors] = ptbConeIsomerizationsFromSpectralRadiance(spd_isetbio,wls_isetbio,integrationTimeSec)
%
% Compute LMS human cone isomerizations from spectral radiance in Watts/[m2-sr-nm].
%
% This routine is set up for a quick commparison to isetbio calculations.  The underlying
% code is demonstrated and (sort of) documented in PTB routine IsomerizationsInEyeDemo.
%
% 8/4/13  dhb  Wrote it.

%% Set up PTB photoreceptors structure
whatCalc = 'LivingHumanFovea';
photoreceptors = DefaultPhotoreceptors(whatCalc);
photoreceptors.eyeLengthMM.source = 'LeGrand';
photoreceptors = FillInPhotoreceptors(photoreceptors);
S = photoreceptors.nomogram.S;

%% XYZ color matching functions
load T_xyz1931
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);
T_Y = T_xyz(2,:);

%% Convert isetbio radiance, which is in units of power/nm to PTB
% convention of power/wlband.
%
% Spline to same wavelength spaceing used here
S_isetbio = WlsToS(wls_isetbio);
spd_temp = SplineRaw(S_isetbio,spd_isetbio,S);
radianceWattsPerM2Sr = spd_temp/S(2);
		
% Find pupil area, needed to get retinal irradiance.  We compute
% pupil area based on the luminance of stimulus according to the
% algorithm specified in the photoreceptors structure.
theXYZ = T_xyz*radianceWattsPerM2Sr; theLuminance = theXYZ(2);
[pupilDiamMm,pupilAreaMm2] = PupilDiameterFromLum(theLuminance,photoreceptors.pupilDiameter.source);
photopicLuminanceCdM2 = T_Y*radianceWattsPerM2Sr;

% Convert radiance of source to retinal irradiance and convert to quantal units.
irradianceWattsPerUm2 = RadianceToRetIrradiance(radianceWattsPerM2Sr,S, ...
    pupilAreaMm2,photoreceptors.eyeLengthMM.value);

%% Do the work in toolbox function
[isoPerConeSec,absPerConeSec,photoreceptors] = ...
    RetIrradianceToIsoRecSec(irradianceWattsPerUm2,S,photoreceptors);
isoPerCone = isoPerConeSec*integrationTimeSec;