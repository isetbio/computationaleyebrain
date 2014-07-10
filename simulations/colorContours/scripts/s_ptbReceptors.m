%% s_ptbReceptors
%

%% Example of how PTB creates photoreceptor absorbance functions

whatCalc = 'CIE2Deg';
photoreceptors = DefaultPhotoreceptors(whatCalc);

focalLengthMm = 17;
photoreceptors.eyeLengthMM.source = num2str(focalLengthMm);

wave = 400:1:770;
photoreceptors.nomogram.S = WlsToS(wave(:));
S = photoreceptors.nomogram.S;

% mPigmentAdjustment = 0;
% if (nargin > 5 && ~isempty(mPigmentAdjustment))
%     photoreceptors.macularPigmentDensity.adjustDen = mPigmentAdjustment;
% end

photoreceptors = FillInPhotoreceptors(photoreceptors);

vcNewGraphWin;
plot(wave,photoreceptors.energyFundamentals');
set(gca,'xtick',400:50:700);

%%  Now, read the Stockman fundamentals stored in ISETBIO

stockman = vcReadSpectra('stockmanEnergy',wave);
vcNewGraphWin; plot(wave,stockman);
set(gca,'xtick',400:50:700);

%% Compare PTB and ISETBIO directly

vcNewGraphWin;
for ii=1:3
   plot(stockman(:,ii),photoreceptors.energyFundamentals(ii,:)');
   hold on
end
axis equal
grid on

%% 
vcNewGraphWin;
plot(wave,stockman);
hold on;
plot(wave,photoreceptors.energyFundamentals');
title('PTB and ISETBIO Stockman Energy Fundamentals');
xlabel('wavelength (nm)')
ylabel('normalized sensitivity')

%% Macular data

% ISETBIO and PTB have different default macular densities.
% But the functions could be made the same
m = macularCreate;
m = macularSet(m,'density',.35);
macT= macularGet(m,'transmittance');
w = macularGet(m,'wave');

vcNewGraphWin;
plot(wave,photoreceptors.macularPigmentDensity.transmittance,'b-');
hold on
plot(w,macT,'r-');
legend('PTB','ISETBIO')
title('Macular pigment transmission')

%% Lens data

% ISETBIO and PTB have different default macular densities.
% But the functions could be made the same
lens = lensCreate;
lensT = lensGet(lens,'transmittance');
w = lensGet(lens,'wave');

vcNewGraphWin;
plot(wave,photoreceptors.lensDensity.transmittance,'b-');
hold on
plot(w,lensT,'r-');
legend('PTB','ISETBIO')
title('Lens transmission')

%% Remove lens and macular transmittance loss from the cones

% stockman = cones*lens*macular
% So we want
% cones alone = stockman/(lens*macular)

stockman = vcReadSpectra('stockmanEnergy',w);
cones = diag(1 ./ (lensT(:).*macT(:)))*stockman;

% scale the cones here
mx = max(cones);
cones = cones*diag(1./mx);

% Plot em up
vcNewGraphWin;
plot(w,cones);


%% End


