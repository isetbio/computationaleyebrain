%% s_ptbReceptors
%

%% Example of how PTB creates photoreceptor absorbance functions

whatCalc = 'CIE2Deg';
photoreceptors = DefaultPhotoreceptors(whatCalc);

focalLengthMm = 17;
photoreceptors.eyeLengthMM.source = num2str(focalLengthMm);

wave = 400:5:770;
photoreceptors.nomogram.S = WlsToS(wave(:));
S = photoreceptors.nomogram.S;

% mPigmentAdjustment = 0;
% if (nargin > 5 && ~isempty(mPigmentAdjustment))
%     photoreceptors.macularPigmentDensity.adjustDen = mPigmentAdjustment;
% end

photoreceptors = FillInPhotoreceptors(photoreceptors);

vcNewGraphWin;
plot(wave,photoreceptors.energyFundamentals');
set(gca,'xtick',400:50:800);

%%  Now, read the Stockman fundamentals stored in ISETBIO

stockman = vcReadSpectra('stockmanEnergy',wave);
vcNewGraphWin; plot(wave,stockman);
set(gca,'xtick',400:50:800);

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



%%
