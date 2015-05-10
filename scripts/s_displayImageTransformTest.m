%% s_displayImageTransformTest
%   This script is used to test built-in color transform and display color
%   converter.
%
%  (HJ) ISETBIO TEAM, 2015

% Init
ieInit;

% Load image and create display
d = displayCreate('OculusDK2_PR650');
I_srgb = im2double(imread('gretag-sRGB-350.jpg'));

% Compute linear rgb and xyz value
I_lrgb = srgb2lrgb(I_srgb);
I_xyz  = srgb2xyz(I_srgb);

% Now xyz is normalized xyz and we would like to assume some mean-luminance
% for the image on the display
% Actually, if we assume mean luminance equals to 1, we can skip this step
meanL = 50; % mean luminance level
I_XYZ = I_xyz * meanL;

% Convert to display linear rgb
I_dlrgb = imageLinearTransform(I_XYZ, inv(displayGet(d, 'rgb2xyz')));
I_dlrgb(I_dlrgb < 0) = 0;

% Apple display inverse gamma distortion
I_dac = ieLUTLinear(I_dlrgb, displayGet(d, 'inverse gamma table'));

% Now, let's go backward and compute srgb values from display dac values
scene = sceneFromFile((I_dac-1)/(displayGet(d, 'nlevels')-1),'rgb',[],d);
I_rsrgb = sceneGet(scene, 'rgb image');

% Show the original and reconstructed image
vcNewGraphWin;
subplot(1,2,1); imshow(I_srgb); title('Original Image');
subplot(1,2,2); imshow(I_rsrgb); title('Reconstructed Image');

vcNewGraphWin; hold on;
plot(I_srgb(:), I_rsrgb(:)); xlabel('Original'); ylabel('Reconstructed');
plot([0 1], [0 1], 'r--', 'LineWidth', 2);
grid on;
