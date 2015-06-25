%% s_colorblindVideoxy
%    This scripts create a video that compares performance of linear and
%    Brettel's dichromatic color transform in xy space
%
%  (HJ) ISETBIO TEAM, 2015

%% Init
ieInit;
load('Alpern_data.mat');
XYZ  = ieReadSpectra('XYZ', wave);
cbType = 'Tritan'; % color blind type
n = length(wave);

%% Plot visible region in xy
vcNewGraphWin; % new figure window
hold on;

xy = bsxfun(@rdivide, XYZ, sum(XYZ, 2));
plot([xy(:,1);xy(1,1)], [xy(:,2);xy(1,2)], '--', 'LineWidth', 2);

%% Plot Brettel transform
wp = ieXYZFromEnergy(1e-3*ones(1, n), wave);
cbXYZ = lms2xyz(xyz2lms(reshape(XYZ, [n 1 3]), cbType, 'Brettel', wp));
% cbXYZ = lms2xyz(xyz2lms(reshape(XYZ, [n 1 3]), 1, wp));
cbXYZ = squeeze(cbXYZ);

cbXYZ_b = cbXYZ;
cbxy_b = bsxfun(@rdivide, cbXYZ, sum(cbXYZ, 2));

indx = inpolygon(cbxy_b(:,1), cbxy_b(:,2), xy(:,1), xy(:,2));
plot(cbxy_b(indx,1), cbxy_b(indx,2), 'o')

%% Plot Linear transform
cbXYZ = lms2xyz(xyz2lms(reshape(XYZ, [n 1 3]), cbType, 'Linear'));
cbXYZ = squeeze(cbXYZ);

cbXYZ_l = cbXYZ;
cbxy_l = bsxfun(@rdivide, cbXYZ, sum(cbXYZ, 2));
indx = inpolygon(cbxy_l(:,1), cbxy_l(:,2), xy(:,1), xy(:,2));
plot(cbxy_l(indx,1), cbxy_l(indx,2), 'o')

xlabel('CIE-x'); ylabel('CIE-y'); grid on;

%% Non-negative Constraints
LMS  = squeeze(xyz2lms(reshape(XYZ, [n 1 3])));
dLMS = dColorTransform(LMS, 3);
cbXYZ = lms2xyz(reshape(dLMS, [n 1 3]));
cbXYZ = squeeze(cbXYZ);

cbXYZ_n = cbXYZ;
cbxy_n = bsxfun(@rdivide, cbXYZ, sum(cbXYZ, 2));

indx = inpolygon(cbxy_n(:,1), cbxy_n(:,2), xy(:,1), xy(:,2));
plot(cbxy_n(indx,1), cbxy_n(indx,2), 'o')

%% Comparison with Alpern Data in units of DeltaE
%  Assuming the luminance level is 50 cd/m^2 and equal energy gray
%  whitepoint with luminance 100
lum = 50;    % stimulus brightness
wplum = 100; % white point brightness
alpern_XYZ = bsxfun(@times, [data, 1-sum(data, 2)], lum./data(:,2));
cbXYZ_bs   = bsxfun(@times, cbXYZ_b, lum./cbXYZ_b(:,2));
cbXYZ_ns   = bsxfun(@times, cbXYZ_n, lum./cbXYZ_n(:,2));

wp = wp / wp(2) * wplum;
de = zeros(length(wave), 2); % first column for brettel, second for SET
de(:, 1) = deltaEab(alpern_XYZ, cbXYZ_bs, wp, '2000');
de(:, 2) = deltaEab(alpern_XYZ, cbXYZ_ns, wp, '2000');

%% Make a video
hfig = vcNewGraphWin; hold on; grid on;
plot([xy(:,1);xy(1,1)], [xy(:,2);xy(1,2)], '--', 'LineWidth', 2);

videoObj = VideoWriter('Tritanope.avi');
videoObj.FrameRate = 2;
open(videoObj);

% plot Alpern
lColor = [136 204 238; 221 204 119; 204 102 119]/255;
for ii = 1 : n
    plot(data(ii, 1), data(ii, 2), 'o', 'Color', lColor(1,:));
    if ii > 1
        plot(data(ii-1:ii, 1), data(ii-1:ii, 2), ...
            'LineWidth', 2, 'Color', lColor(1,:))
    end

    plot(cbxy_b(ii, 1), cbxy_b(ii,2), 'o', 'Color', lColor(2,:));

    plot(cbxy_l(ii, 1), cbxy_l(ii, 2), 'x', 'Color', lColor(3,:));
    % plot(cbxy_n(ii, 1), cbxy_n(ii, 2), 'x', 'Color', lColor(3,:));
    
    axis([0 0.8 0 0.9]);
    drawnow; writeVideo(videoObj, getframe(hfig));
end

close(videoObj);