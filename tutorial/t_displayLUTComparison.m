%% t_displayLUTComparison
%
%    Explore the consequences of having 8 and 10 bit LUT control of the
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
n = which('t_displayLUTComparison.m');
chdir(fileparts(n))

%% Use the 10bit Apple display as an example

% Specific display
d10 = displayCreate('LCD-Apple');
% displayPlot(d10,'gamma');

% Store display white point for cielab calculation
wp = displayGet(d10,'white point');  % XYZ of display white

% Row and column size of the scene data.  Small for speed.  Larger won't
% change anything.
r = 2; c = 2;

% Sample spacing for the base colors.
samps = (0.2:0.3:0.8);
[U,V,W] = meshgrid(samps,samps,samps);
baseColors = [U(:),V(:),W(:)];
nBase = size(baseColors,1);

% Select the delta steps around the base colors We sample at discrete
% intervals around the base colors.  These discrete intervals are turned
% into fractional steps appropriate for the number of levels (10 bit or 8
% bit) below.

% This is enough sample directions to get all the discrete directions, 26
% for 1 step and more for two steps
N = 40; [X,Y,Z] = sphere(N);
delta = [];
nsteps = 2;
for rr = 1:nsteps  
    tmp = round(rr*[X(:),Y(:),Z(:)]);
    delta = [delta; unique(tmp,'rows')]; %#ok<AGROW>
end
fprintf('%i delta steps\n',size(delta,1));

%% Loop on the base color for the display

dE = zeros(size(delta,1),nBase,2);

wb = waitbar(0,'Base color fraction');
for bc = 1:nBase
    
    waitbar(bc/nBase,wb);

    % Get the XYZ of the base color
    % This could be a displayGet(d,'xyz',rgbValue) some day.
    nLevels = displayGet(d10,'nlevels');
    rgb0  = baseColors(bc,:);
    img   = RGB2XWFormat(ones(r,c,3));
    img   = XW2RGBFormat(img*diag(rgb0),r,c);
    scene = sceneFromFile(img,'rgb',[],d10);
    xyz0  = mean(RGB2XWFormat(sceneGet(scene,'xyz')));
    % vcAddObject(scene); sceneWindow;
    
    % Calculate the delta E to nearby RGB values
    for ii=1:size(delta,1)
        rgb1 = rgb0 + (delta(ii,:)/nLevels);  % Notice scaling of the step
        img  = RGB2XWFormat(ones(r,c,3));
        img  = XW2RGBFormat(img*diag(rgb1),r,c);
        scene = sceneFromFile(img,'rgb',[],d10);
        xyz1  = mean(RGB2XWFormat(sceneGet(scene,'xyz')));
        
        % We don't seem to have a deltaE(xyz1,xyz2,wp) function in ISETBIO
        % yet.
        dE(ii,bc,1) = norm(ieXYZ2LAB(xyz0,wp) - ieXYZ2LAB(xyz1,wp));
    end
    
    %% Resample the gamma to 256 levels by interpolation
    
    newLevels = 2^8;
    
    % This keeps the display about the same, but just gives finer DAC control
    g       = displayGet(d10,'gamma');
    nLevels = displayGet(d10,'nlevels');
    dvs     = (0:(nLevels-1))/ nLevels;
    newDVs  = (0:newLevels - 1)/newLevels;
    
    newG = zeros(newLevels,3);
    for ii=1:3
        newG(:,ii) = interp1(dvs,g(:,ii),newDVs,'pchip');
    end
    % vcNewGraphWin; plot(dvs,g(:,2),'--',newDVs,newG(:,2),'.r')
    
    % Create the equivalent 8 bit display
    d8 = displaySet(d10,'gamma',newG);
    
    %% Pick and rgb value and calculate the spectral radiance for d10
    
    nLevels = displayGet(d8,'nlevels');
    
    % For the base color and the new display
    scene = sceneFromFile(img,'rgb',[],d8);
    xyz0 = mean(RGB2XWFormat(sceneGet(scene,'xyz')));
    
    % For nearby RGB values
    for ii=1:size(delta,1)
        rgb1 = rgb0 + (delta(ii,:)/nLevels);
        img = RGB2XWFormat(ones(r,c,3));
        img = XW2RGBFormat(img*diag(rgb1),r,c);
        scene = sceneFromFile(img,'rgb',[],d8);
        xyz1 = mean(RGB2XWFormat(sceneGet(scene,'xyz')));
        dE(ii,bc,2) = norm(ieXYZ2LAB(xyz0,wp) - ieXYZ2LAB(xyz1,wp));
    end
    
end
delete(wb);

%% Make the comparison plots

vcNewGraphWin; 
xvalues = 0:.05:3;
subplot(2,1,1)
dE10 = dE(:,:,1); hist(dE10(:),xvalues);
set(gca,'xlim',[0 3],'ylim',[0 500]);
h = findobj(gca,'Type','patch');
set(h(1),'FaceColor',[.8 0 .2],'EdgeColor','k');
xlabel('\Delta E'); ylabel('Count');
grid on
title(sprintf('10-bit display, %i steps',nsteps))

subplot(2,1,2)
dE8 = dE(:,:,2);  hist(dE8(:),xvalues);
set(gca,'xlim',[0 3],'ylim',[0 500]);
h = findobj(gca,'Type','patch');
set(h(1),'FaceColor',[.2 .2 .8],'EdgeColor','k');
xlabel('\Delta E'); ylabel('Count');
grid on
title(sprintf('8-bit display, %i steps',nsteps))

%%
print -depsc Display1.eps

%% End