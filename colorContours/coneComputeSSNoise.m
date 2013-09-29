function coneResp = coneComputeSSNoise(coneRespSSNF, coneType, opts)
%% function coneComputeSSNoise(coneRespSSNF, coneType, [opts])
%    Compute second site noise based on the cone absorptions
%    Currently, this only handles RGB uniform scene, we need to think about
%    how to handle more general scenes
%    In computation, L and M cones are randomly wired up and each L cone
%    response is substracted by a portion of random M cone and the same
%    thing applies to M cones
%
%  Inputs:
%    coneRespSSNF  - matrix of cone absorption with second site noise
%    coneType      - matrix of same size as coneRespSSNF, coneType
%                    specifies the type of each cone by number 1 ~ 4 for
%                    (K, L, M, S) correspondingly
%    opts          - option structure, contain weights, surrounding type
%                    and fano factor information
%
%  Outputs:
%    coneResp      - matrix of cone response with second site noise added
%
%  See also:
%    ccAccuracy
%
%  BW/HJ/DHB ISETBIO Team, 2013

%% Check inputs and init
if nargin < 1, error('cone absoprtion matrix required'); end
if nargin < 2, error('cone type matrix required'); end
if nargin < 3, opts = []; end
if size(coneRespSSNF, 1) ~= size(coneType, 1) || ...
   size(coneRespSSNF, 2) ~= size(coneType, 2)
    error('cone response and cone type matrix size mismatch');
end

coneType = repmat(coneType, [1 1 size(coneRespSSNF,3)]);

if ~isfield(opts, 'weights'), opts.weights = 0.2; end
if ~isfield(opts, 'ssFanoFactor'), opts.ssFanoFactor = 1; end

%% Wire L and M cones
LConeResp = coneRespSSNF(coneType == 3);
MConeResp = coneRespSSNF(coneType == 4);

nLCones = length(LConeResp);
nMCones = length(MConeResp);

coneResp = coneRespSSNF;
coneResp(coneType == 3) = LConeResp - opts.weights * ...
                          MConeResp(randi(nMCones, [nLCones 1]));
coneResp(coneType == 4) = MConeResp - opts.weights * ...
                          LConeResp(randi(nLCones, [nMCones 1]));
%% Add second site noise
%  The noise is governed by the specified Fano factor, although we compute
%  the variance from the magnitude of the sample value, rather from a mean
%  over time. That is mostly a matter of convenience, and perhaps not
%  unreasonable biophysically in any case

if opts.ssFanoFactor > 0
    stdMat = sqrt(opts.ssFanoFactor*abs(coneResp));
    coneResp = coneResp + stdMat.*randn(size(coneResp));
end

end