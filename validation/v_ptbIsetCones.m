%% v_ptbIsetCones
%
%
% HJ/BW ISETBIO Team, 2013

%% This is how PTB filles in the photoreceptor values
whatCalc = 'CIE2Deg';
focalLengthMm = 16.7;
wave = 400:700;
photoreceptors = DefaultPhotoreceptors(whatCalc);
photoreceptors.eyeLengthMM.source = num2str(focalLengthMm);
photoreceptors.nomogram.S = WlsToS(wave(:));
S = photoreceptors.nomogram.S;
% if (nargin > 5 && ~isempty(mPigmentAdjustment))
%     photoreceptors.macularPigmentDensity.adjustDen = mPigmentAdjustment;
% end
photoreceptors = FillInPhotoreceptors(photoreceptors);

%% These are the normalized cone fundamentals

vcNewGraphWin;
% They can be specified for calculating with photons
plot(wave,photoreceptors.quantalFundamentals)

% Or they can be specified for calculating with energy
plot(wave,photoreceptors.energyFundamentals)

%% Now, get the ISET stuff


%% END