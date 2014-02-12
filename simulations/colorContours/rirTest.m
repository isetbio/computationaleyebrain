sz = sceneGet(scene,'size');
% rect = [sz(2)/4,sz(1)/4,sz(2)/2,sz(1)/2];
rect = [1,1,sz(2)-1,sz(1)-1];


roiLocs = ieRoi2Locs(rect);
radiance = vcGetROIData(scene,roiLocs,'energy');
radiance = mean(radiance);

sz   = oiGet(oi,'size');
% rect = [sz(2)/4,sz(1)/4,sz(2)/2,sz(1)/2];
rect = [1,1,sz(2)-1,sz(1)-1];

roiLocs    = ieRoi2Locs(rect);
irradiance = vcGetROIData(oi,roiLocs,'energy');
irradiance = mean(irradiance);

vcNewGraphWin; plot(radiance(:) ./ irradiance(:));
