%% t_displayLUTComparison
%
%    Explore the consequences of having 8, 10 and 12 bit LUT control of the
%    display primaries.
%
%    The method is this
%
%      Load a display and set the dac size (8, 10, 12 bits)
%      Calculate the display spectral radiance for a small RGB patch
%      Calculate the spectral radiance of nearby RGB values
%      Compute the CIELAB (delta E) w.r.t. the display white point
%
%  (BW) ISETBIO TEAM, 2015

%% 
ieInit

%% Use the 10bit Apple display as an example

d10 = displayCreate('LCD-Apple');
% displayPlot(d10,'gamma');
fprintf('First display has a %i-bit dac size\n',displayGet(d10,'dac size'));
wp = displayGet(d10,'white point');  % XYZ of display white

% Row and column size of the scene data
r = 4; c = 4;
baseColor = [.5 .5 .5];   % Middle of range
% baseColor = [.2 .2 .2];   % Middle of range
% baseColor = [.8 .8 .8];   % Middle of range

%% Pick and rgb value and calculate the spectral radiance

% We will sample at delta intervals around these points
N = 20; [X,Y,Z] = sphere(N);
delta = [];
for radius = 1:2
    tmp = round(radius*[X(:),Y(:),Z(:)]);
    delta = [delta; unique(tmp,'rows')]; %#ok<AGROW>
end
dE = zeros(size(delta,1),2);

%% Set the base color for the display

nLevels = displayGet(d10,'nlevels');

rgb0 = round(baseColor*nLevels);
img = RGB2XWFormat(ones(r,c,3));
img = XW2RGBFormat(img*diag(rgb0),r,c);
scene = sceneFromFile(img,'rgb',[],d10);
xyz0 = mean(RGB2XWFormat(sceneGet(scene,'xyz')));
% vcAddObject(scene); sceneWindow;

% Calculate the delta E to nearby RGB values
for ii=1:size(delta,1)
    rgb1 = rgb0 + delta(ii,:);
    img = RGB2XWFormat(ones(r,c,3));
    img = XW2RGBFormat(img*diag(rgb1),r,c);
    scene = sceneFromFile(img,'rgb',[],d10);
    xyz1 = mean(RGB2XWFormat(sceneGet(scene,'xyz')));
    dE(ii,1) = norm(ieXYZ2LAB(xyz0,wp) - ieXYZ2LAB(xyz1,wp));
end

%% Resample the gamma curve by interpolation

newLevels = 2^8;

% This keeps the display about the same, but just gives finer DAC control
g = displayGet(d10,'gamma');
nLevels = displayGet(d10,'nlevels');
dvs = (0:(nLevels-1))/ nLevels;
newDVs = (0:newLevels - 1)/newLevels;

newG = zeros(newLevels,3);
for ii=1:3
    newG(:,ii) = interp1(dvs,g(:,ii),newDVs,'pchip');
end
% vcNewGraphWin; plot(dvs,g(:,2),'--',newDVs,newG(:,2),'.r')

d8 = displaySet(d10,'gamma',newG);
fprintf('Second display has a %i-bit dac size\n',displayGet(d8,'dac size'));

%% Pick and rgb value and calculate the spectral radiance for d10

nLevels = displayGet(d8,'nlevels');

% For the base color
rgb0 = round(baseColor*nLevels);
img = RGB2XWFormat(ones(r,c,3));
img = XW2RGBFormat(img*diag(rgb0),r,c);
scene = sceneFromFile(img,'rgb',[],d8);
xyz0 = mean(RGB2XWFormat(sceneGet(scene,'xyz')));
% vcAddObject(scene); sceneWindow;

% For nearby RGB values
for ii=1:size(delta,1)
    rgb1 = rgb0 + delta(ii,:);
    img = RGB2XWFormat(ones(r,c,3));
    img = XW2RGBFormat(img*diag(rgb1),r,c);
    scene = sceneFromFile(img,'rgb',[],d8);
    xyz1 = mean(RGB2XWFormat(sceneGet(scene,'xyz')));
    dE(ii,2) = norm(ieXYZ2LAB(xyz0,wp) - ieXYZ2LAB(xyz1,wp));
end

%% Make the comparison plot

vcNewGraphWin; hist(dE,20)
h = findobj(gca,'Type','patch');
set(h(1),'FaceColor',[.8 0 .2],'EdgeColor','k');
set(h(2),'FaceColor',[.2 0 .8],'EdgeColor','k');
xlabel('\Delta E'); ylabel('Count');
legend('10 bit','8 bit');
title(sprintf('Base RGB [%.2f %.2f %.2f]',baseColor(1), baseColor(2), baseColor(3)))
grid on

% set(gca,'xlim',[0 2.5])

%% Now compare the two gamma conditions

dGamma = d8;
nLevels = displayGet(dGamma,'nlevels');
gam = 1.8;

newG = zeros(nLevels,3);
for ii=1:3
    newG(:,ii) = ((0:(nLevels-1)) / nLevels) .^ gam;
end
% vcNewGraphWin; plot(newG(:,2),'.r')

dGamma = displaySet(dGamma,'gamma',newG);
fprintf('Third display has gamma of 1.8\n');
% 
dGamma = displaySet(dGamma,'name',sprintf('Gam = %.1f',gam));
vcAddObject(dGamma); displayWindow;

%% delta E Error distribution
dE = zeros(size(delta,1),2);

% For the base color
rgb0 = round(baseColor*nLevels);
img = RGB2XWFormat(ones(r,c,3));
img = XW2RGBFormat(img*diag(rgb0),r,c);
scene = sceneFromFile(img,'rgb',[],dGamma);
xyz0 = mean(RGB2XWFormat(sceneGet(scene,'xyz')));
% vcAddObject(scene); sceneWindow;

% For nearby RGB values
for ii=1:size(delta,1)
    rgb1 = rgb0 + delta(ii,:);
    img = RGB2XWFormat(ones(r,c,3));
    img = XW2RGBFormat(img*diag(rgb1),r,c);
    scene = sceneFromFile(img,'rgb',[],dGamma);
    xyz1 = mean(RGB2XWFormat(sceneGet(scene,'xyz')));
    dE(ii,1) = norm(ieXYZ2LAB(xyz0,wp) - ieXYZ2LAB(xyz1,wp));
end

%% Adjust the gamma

nLevels = displayGet(dGamma,'nlevels');
gam = 2.2;

newG = zeros(nLevels,3);
for ii=1:3
    newG(:,ii) = ((0:(nLevels-1)) / nLevels) .^ gam;
end
% vcNewGraphWin; plot(dvs,g(:,2),'--',newDVs,newG(:,2),'.r')

dGamma = displaySet(dGamma,'gamma',newG);
fprintf('Fourth display has gamma of 2.2\n');

dGamma = displaySet(dGamma,'name',sprintf('Gam = %.1f',gam));
vcAddObject(dGamma); displayWindow;

%% Error distribution

% For the base color
rgb0 = round(baseColor*nLevels);
img = RGB2XWFormat(ones(r,c,3));
img = XW2RGBFormat(img*diag(rgb0),r,c);
scene = sceneFromFile(img,'rgb',[],dGamma);
xyz0 = mean(RGB2XWFormat(sceneGet(scene,'xyz')));
% vcAddObject(scene); sceneWindow;

% For nearby RGB values
for ii=1:size(delta,1)
    rgb1 = rgb0 + delta(ii,:);
    img = RGB2XWFormat(ones(r,c,3));
    img = XW2RGBFormat(img*diag(rgb1),r,c);
    scene = sceneFromFile(img,'rgb',[],dGamma);
    xyz1 = mean(RGB2XWFormat(sceneGet(scene,'xyz')));
    dE(ii,2) = norm(ieXYZ2LAB(xyz0,wp) - ieXYZ2LAB(xyz1,wp));
end

%% Make the gamma plot comparison plot

vcNewGraphWin; hist(dE,20)
h = findobj(gca,'Type','patch');
set(h(1),'FaceColor',[.8 0 .2],'EdgeColor','k');
set(h(2),'FaceColor',[.2 0 .8],'EdgeColor','k');
xlabel('\Delta E'); ylabel('Count');
legend('1.8 gamma','2.2 gamma');
title(sprintf('Base RGB [%.2f %.2f %.2f]',baseColor(1), baseColor(2), baseColor(3)))
grid on


%% End