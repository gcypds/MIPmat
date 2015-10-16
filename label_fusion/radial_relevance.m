function [mrho,srho,rho] = radial_relevance(A,rmax)
%USAGE: [mrho srho] = radial_relevance(A,rmax)
%INPUTS:
% A: Cell array projection matrices. Each cell holds a matrix of N rows.
%    Number of columns may change.
% rmax: Maximum radius.
%OUTPUTS:
% mrho: Mean radial relevance computed along the cells
% srho: Standard deviarion of the radial relevance computed along the cells
%
% Created on Fri Oct 16 15:26:46 2015
% David Cardenas Pena - GCPDS

sum2 = @(X)sum(X,2);
A = cellfun(@abs,A,'UniformOutput',false);
sumA = cellfun(sum2,A,'UniformOutput',false);
sumA = cell2mat(sumA);

r = -rmax:rmax;
[x(:,:,:,1),x(:,:,:,2),x(:,:,:,3)] = ndgrid(r,r,r);
r = max(abs(x),[],4);

R = unique(r);
rho = zeros(numel(A),numel(R));
for i=1:numel(R)
  ind = r(:) == R(i);
  rho(:,i) = sum(sumA(ind,:),1)'/sum(ind); %mean by radius
end

rho = bsxfun(@times,rho,1./sum(rho,2));
mrho = mean(rho,1);
srho = std(rho,[],1);
