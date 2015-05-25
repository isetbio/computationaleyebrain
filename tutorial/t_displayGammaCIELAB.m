%% t_displayGammaCIELAB
%
%  When a display matched in all ways has a gamma function of 1.8 or 2.2,
%  the same RGB produces two different XYZ values.  How different are they?
%
%  (JEF/BW) ISETBIO TEAM, 2015

%% Init variables

ieInit;         % init a new isetbio session


%% Set up the display and sample levels

d = displayCreate('LCD-Apple');
lrgb2xyz = displayGet(d,'lrgb2xyz');
wp = displayGet(d,'white point');
nLevels = displayGet(d,'nlevels');
l = (1:nLevels)/nLevels; l = l(:);

gTable = zeros(nLevels,3);
g1 = 1.8; gTable1 = repmat(l.^g1,1,3);
g2 = 2.4; gTable2 = repmat(l.^g2,1,3);

    
s = 1:50:nLevels;
[R,G,B] = meshgrid(s,s,s);
v = [R(:), G(:), B(:)];
%%
% v = [512 512 512];
% v = [128 512 802];
% v = [1023 1023 1023];
nSamp = size(v,1);
d = zeros(nSamp,1);
for vv = 1:nSamp
    
    lrgb = ieLUTDigital(v(vv,:),gTable1);
    xyz1 = lrgb*rgb2xyz;
    
    lrgb = ieLUTDigital(v(vv,:),gTable2);
    xyz2 = lrgb*rgb2xyz;
    
    d(vv) = deltaEab(xyz1,xyz2,wp);
end

%%
vcNewGraphWin
hist(d,50);
xlabel('\Delta E');
ylabel('Count');
title(sprintf('Comparison between \\gamma = %.1f and %.1f ',g1,g2))
median(d)

%% END