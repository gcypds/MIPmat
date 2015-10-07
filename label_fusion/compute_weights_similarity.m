function W = compute_weights_similarity(X,Y)
%Similarity-based weight computation.
%USAGE:
% W = compute_weights_npbs(X,Y)  returns a matrix W containing the 
%     similarities between each pair of observations in the p-by-nx data
%     matrix X and p-by-n data matrix Y. Columns of X and Y correspond to
%     observations, rows correspond to variables. W is an nx-by-ny matrix, 
%     with the (i,j) entry equal to similarity between observation i in X
%     and observation j in Y, given by:
%     w(i,j) = exp(-(X(:,i)-Y(:,j))'*(X(:,i)-Y(:,j))/(2*h^2)).
%     The scale parameter h is automatically estimated.
%
% Created on Wed Oct  7 11:53:17 2015
% David Cardenas Pena - GCPDS


W = pdist2(X',Y');
h = median(pdist(X')); %scale parameter
W=exp(-(W.^2)/(2*h^2));