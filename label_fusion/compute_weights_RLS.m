
function W=compute_weights_RLS(X,Y,rho)
%Weight computation by solving a regularized least squares problem.
%USAGE:
% W = compute_weights_RLS(X,Y,rho)  returns a matrix W containing the 
%     solution to the problem:
%              min ||WX-Y||^2_2 + rho*||X||_1
%     given the observations in the p-by-nx data matrix X and p-by-n data 
%     matrix Y. Columns of X and Y correspond to observations, rows 
%     correspond to variables. W is an nx-by-ny matrix, with the (i,j) 
%     weight of the observation i in X to build the observation j in Y.
%     The regularization parameter rho is a non-negative scalar.
%
% Created on Wed Oct  7 14:52:55 2015
% Mauricio Orbes Arteaga - GCPDS

opts=[];
% Starting point
opts.init=2;        % starting from a zero point
opts.tFlag=5;       % run .maxIter iterations
opts.maxIter=100;   % maximum number of iterations
% normalization
opts.nFlag=0;       % without normalization
% regularization
opts.rFlag=1;       % the input parameter 'rho' is a ratio in (0, 1)
%opts.rsL2=0.01; 

[W, funVal]= nnLeastR(X, Y, rho, opts);