function [BTH,LBL,neighs] = getNeighborhood(X,L,r,alpha,beta)
% getNeighborhood: extracts patches in a neighborhood from an array of
% 3D volumes.
% USAGE: [BTH,LBL] = getNeighborhood(X,L,r,alpha,beta)
% INPUTS:
% X:     Array of 3D intensity volumes. size(X) = [W H S N]: (W)idth, (H)eight,
%        number of (S)lices, and (N)umber of volumes.
% L:     Array of 3D label volumes. size(X) = [W H S N]: (W)idth, (H)eight, 
%        number of (S)lices, and (N)umber of volumes.
% r:     Neighborhood center: r = [r1 r2 r3]
% alpha: Patch radius. alpha>0.
% beta:  Neighborhood radius. beta>0.
% OUTPUTS:
% BTH:   Matrix of intensity patches. size(BTH) = [Q M]. Q: number of 
%        elements in a patch, M: number of patches
% LBL:   Vector of patch labels. size(BTH) = [1 M]. M: number of patches.
%        If L is empty, LBL will be empty.
% 
% WARNING: No boundary check!
%
% Created on Wed Oct  7 10:23:42 2015
% Mauricio Orbes Arteaga - GCPDS
% David Cardenas Pena - GCPDS

siz = size(X);
if numel(siz)==3
  N = 1;
else
  N = siz(4);
end

X = X(:);
if ~isempty(L)
  L = L(:);
end

[tmp1,tmp2,tmp3,tmp4] = ndgrid(r(1)+(-beta:beta),r(2)+(-beta:beta),r(3)+(-beta:beta),1:N);
neighs = sub2ind(siz,tmp1(:),tmp2(:),tmp3(:),tmp4(:));

[tmp1,tmp2,tmp3] = ndgrid(r(1)+(-alpha:alpha),r(2)+(-alpha:alpha),r(3)+(-alpha:alpha));
patches = sub2ind(siz(1:3),tmp1(:),tmp2(:),tmp3(:)) - sub2ind(siz(1:3),r(1),r(2),r(3));

indices = bsxfun(@plus,patches,neighs');

BTH = reshape(X(indices),size(indices));
if isempty(L)
  LBL = [];
else
  LBL = reshape(L(neighs'),size(neighs'));
end