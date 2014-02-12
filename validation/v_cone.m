%% v_cone
%
% Test cone, lens and macular function calls
%
% BW/HJ ISETBIO Team 2013

%%
cone = coneCreate;
wave = coneGet(cone,'wave');

%%
vcNewGraphWin([],'tall');
subplot(4,1,1)
plot(wave,coneGet(cone,'cone spectral absorptance'));
title('Cone spectral absorptance')

subplot(4,1,2)
lens = coneGet(cone,'lens');
plot(wave,lensGet(lens,'transmittance'))
title('Lens transmittance')

subplot(4,1,3)
macular = coneGet(cone,'macular');
plot(wave,macularGet(macular,'transmittance'))
title('Macular transmittance')

subplot(4,1,4)
plot(wave,coneGet(cone,'effective spectral absorptance'))
title('Cone-ocular absorptance')

%% Plot again, but change the macular pigment density to 0

m = coneGet(cone,'macular');
m = macularSet(m,'density',0);
cone = coneSet(cone,'macular',m);

%%
vcNewGraphWin([],'tall');
subplot(4,1,1)
plot(wave,coneGet(cone,'cone spectral absorptance'));
title('Cone spectral absorptance')

subplot(4,1,2)
lens = coneGet(cone,'lens');
plot(wave,lensGet(lens,'transmittance'))
title('Lens transmittance')

subplot(4,1,3)
macular = coneGet(cone,'macular');
plot(wave,macularGet(macular,'transmittance'))
title('Macular transmittance')

subplot(4,1,4)
plot(wave,coneGet(cone,'effective spectral absorptance'))
title('Cone-ocular absorptance')

%% End