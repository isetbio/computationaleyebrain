%% s_colorblindVideoxy
%    This scripts create a video that compares performance of linear and
%    Brettel's dichromatic color transform in xy space
%
%  (HJ) ISETBIO TEAM, 2015

%% Init
ieInit;
wave = 420:5:700;  % wavelength samples in nm
XYZ  = ieReadSpectra('XYZ', wave);
cbType = 'Tritan'; % color blind type
n = length(wave);

%% Plot visible region in xy
vcNewGraphWin; % new figure window
hold on;

xy = bsxfun(@rdivide, XYZ, sum(XYZ, 2));
plot([xy(:,1);xy(1,1)], [xy(:,2);xy(1,2)], '--');

%% Plot Brettel transform
wp = ieXYZFromEnergy(1e-3*ones(1, n), wave);
cbXYZ = lms2xyz(xyz2lms(reshape(XYZ, [n 1 3]), cbType, 'Brettel', wp));
cbXYZ = squeeze(cbXYZ);

cbxy_b = bsxfun(@rdivide, cbXYZ, sum(cbXYZ, 2));
plot(cbxy_b(:,1), cbxy_b(:,2), 'o')

%% Plot Linear transform
cbXYZ = lms2xyz(xyz2lms(reshape(XYZ, [n 1 3]), cbType, 'Linear'));
cbXYZ = squeeze(cbXYZ);

cbxy_l = bsxfun(@rdivide, cbXYZ, sum(cbXYZ, 2));
plot(cbxy_l(:,1), cbxy_l(:,2), 'o')

xlabel('CIE-x'); ylabel('CIE-y'); grid on;

%% Non-negative Constraints
LMS  = squeeze(xyz2lms(reshape(XYZ, [n 1 3])));
dLMS = dColorTransform(LMS, 3);
cbXYZ = lms2xyz(reshape(dLMS, [n 1 3]));
cbXYZ = squeeze(cbXYZ);

cbxy_n = bsxfun(@rdivide, cbXYZ, sum(cbXYZ, 2));
plot(cbxy_n(:,1), cbxy_n(:,2), 'o')

%% 
vcNewGraphWin; hold on;
for ii = 1 : n
    plot(xy(ii,1), xy(ii,2), '.b');
    plot(cbxy_b(ii, 1), cbxy_b(ii, 2), 'or');
    plot(cbxy_l(ii, 1), cbxy_l(ii, 2), 'xg');
    drawnow; WaitSecs(0.5);
end