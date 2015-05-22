%% t_displayPixelLayoutSCIELAB
%
%    We simulate a letter on different types of displays and compare
%    the S-CIELAB error maps at different viewing distances.
%
%    This calculation is based on a simpler model than the full isetbio.
%    Here, we do not use the optics.  Rather, we use the display model and
%    compute XYZ for each point. We let S-CIELAB do the blurring in
%    opponent-colors space.
%
%  (HJ/BW) ISETBIO TEAM, 2015

%% Init variables

ieInit;         % init a new isetbio session

% Smaller field of view means subject is further away
fovList = (0.15:0.1:0.8);

% displayName = {'CRT-Dell', 'crt', 'Dell-Chevron', 'LCD-HP', 'OLED-Sony'};
displayName = {'Dell-Chevron', 'LCD-HP'};
nd = length(displayName);  %Number of displays

% Make the letter list, all upper and lower case letters
strLower = char(97:122);
strUpper = char(65:90);
letterList = cell(52,1);
jj=1;
for ii=1:26
    letterList{jj} = strLower(ii);
    letterList{jj+1} = strUpper(ii);
    jj = jj + 2;
end
% letterList  = {'g','X'};
nLetters = length(letterList);


%% Show the subpixel images, though they are clear in the scene

for ii = 1 : nd
    % create display and make images of the subpixels
    d = displayCreate(displayName{ii});
    
    % show sub-pixel image zoomed
    vcNewGraphWin;
    dmap = displayGet(d, 'dixel intensity map');
    dmap = dmap(:, :, 1:3); % use RGB primaries only
    imshow(dmap / max(dmap(:))); title(displayGet(d, 'name')); drawnow;
    
end

%% Font

% Below, from the dpi and fov we calculate the viewing distance
% We set up the mean error and the range of the delta E values
medianE = zeros(length(fovList),nLetters);
rng = zeros(length(fovList),nLetters,2);
wbar = waitbar(0,'Letters');
for ll=1:nLetters
    waitbar(ll/nLetters,wbar,'Letters');
    font = fontCreate(letterList{ll});             % letter 'g' or 'w'
    font = fontSet(font,'size',12);
    font.bitmap = 1 - font.bitmap; % flip data to be white letter on black
    padsz = [5 5]; padval = 0;     % padding size and val for sceneFromFont
    
    % Make scenes of the same font on the two displays
    % Match dpi and spd
    clear scene
    clear wp
    for ii=1:nd
        d = displayCreate(displayName{ii});
        
        % Make the SPDs the same
        if ii == 1,    spd = displayGet(d,'spd');
        elseif ii > 1, d = displaySet(d,'spd',spd);
        end
        
        wp{ii} = displayGet(d,'white point');
        d = displaySet(d,'dpi',fontGet(font,'dpi'));
        
        scene{ii} = sceneFromFont(font, d, [], [], padsz, padval);
        scene{ii} = sceneSet(scene{ii},'name',displayName{ii});
        % scene{ii} = sceneAdjustLuminance(scene{ii},meanLuminance);
        vcAddObject(scene{ii});
    end
    % sceneWindow;
    
    % Adjust field of view to set viewing distance
    for ff=1:length(fovList)
        for ii=1:2
            scene{ii} = sceneSet(scene{ii},'fov',fovList(ff));
        end
        
        % Now run scielab on the two scenes
        sampPerDeg = 1/sceneGet(scene{1},'degrees per sample');
        
        % Run S-CIELAB based on CIELAB 2000
        params.deltaEversion = '2000';
        params.sampPerDeg    = sampPerDeg;
        params.imageFormat   = 'xyz';
        params.filterSize    = sampPerDeg;
        params.filters = [];
        params.filterversion = 'distribution';
        errorImage = scielab(sceneGet(scene{1},'xyz'), sceneGet(scene{2},'xyz'), wp, params);
        
        % Crop the error image a bit to eliminate the padding.
        % This rect will depend on the character.
        rect = [30     1   184   266];  % Georgia g, 14
        errorImage = imcrop(errorImage,rect);
        
        medianE(ff,ll) = median(errorImage(:));
        rng(ff,ll,1) = prctile(errorImage(:),0.1);
        rng(ff,ll,2) = prctile(errorImage(:),0.9);
        
        % Show the error images, though this is not all that important
        %     vcNewGraphWin;
        %     imagesc(errorImage);
        %     colorbar;
        
    end
    
end
delete(wbar)

%% Plot the error as a function of viewing distance (meters)

dpi  = fontGet(font,'dpi');
cols = sceneGet(scene{1},'cols');
colsInches = cols/dpi;

% a = tan(colsInches/vd);
vd = colsInches ./ atan(fovList(:));
% 0.0254 m/inch
vd = 0.0254*vd * 100;  % In centimeters

medians = mean(medianE,2);
stdE = std(medianE,1,2);
p = errorbar(vd,medians,stdE);
set(p,'linewidth',2);
grid on; 
xlabel('Viewing distance (cm)')
ylabel('Median \Delta E_S')
set(gca,'ylim',[0 1])


%% END