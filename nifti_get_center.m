% ---------------------------------------------------------------------
function center = nifti_get_center(vol)
% get center in World coords of an SPM volume

dx = norm(vol.mat*[1 0 0 0]' - vol.mat*[0 0 0 0]');
dy = norm(vol.mat*[0 1 0 0]' - vol.mat*[0 0 0 0]');
dz = norm(vol.mat*[0 0 1 0]' - vol.mat*[0 0 0 0]');
center = (0.5.*[dx dy dz].*(vol.dim-1))'; % this is what Jeninson does in rmsdiff.cc (but is wrong? for obliques!!)
%center = 0.5.*(vol.mat(1:3,1:3) * (vol.dim-1)'); % this is correct but we'll go with Jenkinson's code for now
return