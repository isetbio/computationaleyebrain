%% t_displayPixelLayoutFont
%
%    OBSOLETE.
%    see t_displayPixelLayoutSCIELAB
%
%    In this script, we simulate a letter on different types of
%    displays. Then, we compare the scene radiance and optical images at
%    different viewing distance
%
%  See also:  t_displayPixelLayoutLine, t_displayPixelLayoutSCIELAB
%
%  (HJ/BW) ISETBIO TEAM, 2015

%% Init
ieInit; % init a new isetbio session
vd = [0.25 1];
nDist = length(vd);

% Create a line
doSub = true; % subpixel rendering


%% Font on different display
font = fontCreate; % letter g
font.bitmap = 1 - font.bitmap; % flip it to be white letter on black
padsz = [5 5]; padval = 0; % padding size and val for sceneFromFont

displayName = {'CRT-Dell', 'crt', 'Dell-Chevron', 'LCD-HP', 'OLED-Sony'};
nd = length(displayName);

oi = oiCreate('human');  % This is the Marimont and Wandell model, I think

vcNewGraphWin;
for ii = 1 : nd
    % create display
    d = displayCreate(displayName{ii});
    
    % show sub-pixel image zoomed
    subplot(2 + nDist, nd, ii);
    dmap = displayGet(d, 'dixel intensity map');
    dmap = dmap(:, :, 1:3); % use RGB primaires only
    imshow(dmap / max(dmap(:))); title(displayGet(d, 'name')); drawnow;
    
    % create and show font scene on the display
    scene = sceneFromFont(font, d, [], [], padsz, padval);
    % vcAddObject(scene); sceneWindow;
    subplot(2 + nDist, nd, ii + nd);
    imshow(sceneGet(scene, 'rgb image')); drawnow;
    
    % compute optical image at two viewing distances
    width = sceneGet(scene, 'width');
    for jj = 1 : nDist
        % adjust scene field of view according to viewing distance
        scene = sceneSet(scene, 'distance', vd(jj));
        scene = sceneSet(scene, 'h fov', atand(width / vd(jj)));
        
        % compute optical image
        oi = oiCompute(oi, scene);
        
        % show optical image
        subplot(2+nDist, nd, ii + nd * (jj + 1));
        imshow(oiGet(oi, 'rgb image')); drawnow;
    end
end

%% Show scene radiance plot
% vcNewGraphWin;
% for ii = 1 : nd
%     % create display
%     d = displayCreate(displayName{ii});
%     
%     % show sub-pixel image
%     subplot(3, nd, ii);
%     dmap = displayGet(d, 'dixel intensity map');
%     dmap = dmap(:, :, 1:3); % use RGB primaires only
%     imshow(dmap / max(dmap(:))); title(displayGet(d, 'name')); drawnow;
%     
%     % create and show font scene
%     scene = sceneFromFont(font, d, [], [], padsz, padval);
%     subplot(3, nd, ii + nd);
%     imshow(sceneGet(scene, 'rgb image')); drawnow;
%     
%     % compute opitcal image at different viewing distance
%     [uData,hg] = plotScene(scene, 'radiance hline', [1 230]); close(hg);
%     subplot(3, nd, ii + 2*nd);
%     mesh(uData.pos,uData.wave,uData.data)
%     xlabel('Position (mm)'); ylabel('Wavelength (nm)');
%     zlabel('Radiance (q/s/nm/m^2)');
% end

%% Show OI irradiance plot
% sampleWave = [450 550 650];
% for jj = 1 : nDist
%     vcNewGraphWin;
%     for ii = 1 : nd
%         % create display
%         d = displayCreate(displayName{ii});
%         
%         % show sub-pixel image
%         subplot(6, nd, ii);
%         dmap = displayGet(d, 'dixel intensity map');
%         dmap = dmap(:, :, 1:3); % use RGB primaires only
%         imshow(dmap / max(dmap(:))); title(displayGet(d, 'name')); drawnow;
%         
%         % create and show font scene
%         scene = sceneFromFont(font, d, [], [], padsz, padval);
%         width = sceneGet(scene, 'width');
%         
%         % adjust scene field of view according to viewing distance
%         scene = sceneSet(scene, 'distance', vd(jj));
%         scene = sceneSet(scene, 'h fov', atand(width / vd(jj)));
%         
%         % compute optical image
%         oi = oiCompute(oi, scene);
%         
%         % show optical image
%         subplot(6, nd, ii + nd);
%         imshow(oiGet(oi, 'rgb image')); drawnow;
%         
%         % show optical image at different wavelength
%         for kk = 1 : length(sampleWave)
%             subplot(6, nd, ii + nd * (kk+1));
%             p = oiGet(oi, 'photons', sampleWave(kk));
%             imshow(p / max(p(:)));
%         end
%         
%         % show optical image
%         [uData,hg] = plotOI(oi, 'irradiance hline', [1 288]); close(hg);
%         subplot(6, nd, ii + 5*nd);
%         mesh(uData.pos,uData.wave,uData.data)
%         xlabel('Position (um)');
%         ylabel('Wavelength (nm)'); zlabel('Irradiance (q/s/nm/m^2)');
%     end
% end

%% END