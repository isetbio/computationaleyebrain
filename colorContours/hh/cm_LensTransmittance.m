function [lensTrans, wavelength, den_lens] = cm_LensTransmittance(lensTransmitttanceFactor,wavelength,source)

% Calculate lens tranmittance from
%
%   [lensTrans, wavelength, den_lens] =
%   cm_LensTransmittance(lensTransmitttanceFactor,wavelength,source)
%
%
% Copyright HH, Vista Lab, 2010
%
% TODO:  Convert to vcReadSpectra() format for the loads.  If necessary, we
% will put the data into ISET directory.  The decomposition of the Stockman
% for a two part model as suggested by Smith and Pokorny is a good idea
% that Hiroshi implemented.  We should write the code for that.
%
% Hiroshi to put more comments including references to all methods.

if ~exist('wavelength','var'), wavelength = cm_getDefaultWls; end
if ~exist('source','var'), source = 'stockman2'; end

switch source
    
    case {'stockmanF','fixed','Fixed','F'}
        % default stockman's lens pigment densitiy.
        den_lens_ssf = vcReadSpectra('lensDensity',wavelength);
        den_lens = SplineSrf(S_lens_ssf,den_lens_ssf,wavelength,1);
        
    case {'stockman1','stockman'};
        %  lensTransmitttanceFactor should be within 25% (0.75 to 1.25)
        %  Stockman et al. 1999. I'm still not sure how I will deal with
        %  it, that is, if I could multiply the function with a scale factor
        %  linearly or not.
        den_lens_ssf = vcReadSpectra('lensDensity',wavelength);
        den_lens = den_lens_ssf * lensTransmitttanceFactor;
        
    case {'Savage','Linear'} % Savage et al. 1993
        den_lens = lensTransmitttanceFactor * 10.^(5.543-(0.013439 .* wavelength));
        
    case {'twofactors','smith-pokorny'}
        % interpolated smith-pokorny's 1987 function, which fit well
        % stockman 1999, of course but not perfect. 
        % lensTramsmittanceFactor should be  0.76 to 1.56 (see below); 
        % k = 1.00 + 0.02 * (Age -32) in case of 20 < Age < 60
        load den_lens_SmithPokorny.mat
        den_lens = pchipIntLens1' * lensTransmitttanceFactor + pchipIntLens2';
        
    case {'two_stockman','stock2','stockman2'}
        % Decompose Stockman's (1999) lens function into two factors. First
        % truncated data from 460nm to the end, and interpolated linearly
        % (dens
        load('decompStockLens'); % please read comments in the mat file
        
        if ~exist('wave','var'), wave = nm; end
        
        % Pigment densities are calculated with two components
        % 1) shorter wavelength part (fixed)
        % 2) short to middle wavelength part parameterized with one scalar
        
        den_lens = SplineSrf(S_lens_s2f, dens_lens1, wave, 1) * lensTransmitttanceFactor ...
            + SplineSrf(S_lens_s2f, dens_lens2, wave,1);
        
        den_lens = den_lens(find(wavelength(1) == wave):find(wavelength(end) == wave));

                 
    otherwise
        error('Unknown method: %s\n',source)
end

lensTrans = 10.^(-den_lens)';

return

%% just see and compare
% n = 50; C = jet(n); figure, hold on
% Method = 'stockman1';
% Method = 'stockman2';
% for ii = 1:n
%     [lensTrans, wavelength, den_lens] = cm_LensTransmittance(ii*0.1, cm_getDefaultWls, Method);
%     plot(wavelength, lensTrans,'color',C(ii,:));
% end;
% xlim([400 650]);ylim([0 1]); grid on
% 
% % I prefer Stockman2 model so far.... 
