%
% Thinking about scene2display with DHB

xyz = sceneGet(scene,'xyz');

srgb = xyz2srgb(xyz);   % This could be rgb = xyz2display(xyz,d);
% This could also be coordinate system invariant, so that it might work
% from xyz or from stockman or from cones ....

vcNewGraphWin;
imagesc(srgb)
axis image
