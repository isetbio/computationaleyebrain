function coneCreate
% Default cone structure
%
%
% See the PTB routines,
%     DefaultPhotoreceptors
%     ptbConeIsomerizationsFromSpectra
%     FillInPhotoreceptors
%   
%  And check HH's routines.  They need to be integrated into isetbio
%  structures/calls.
%
% I am leaning to coneCreate, rodCreate, iprgcCreate rather than
% photoreceptorCreate.  At least, I will start at the lower level and then
% maybe create the top level one when the separate ones are created.
%
% There are many types of cones, also.  I am thinking they all have the
% same slots, but the parameters differ depending on fovea, periphery, and
% so forth.  I am guessing we stick macular pigment into the parameter
% field, and optical density, and so forth, and that when we do a
% coneGet(cone,'spectral qe') we combine the whole thing.
%
% 
end
