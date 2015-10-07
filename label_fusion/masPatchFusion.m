function L=masPatchFusion(imgqry,supdata,supdata_lbl,opts,schspels)
% Generate the label segmentation by Multi-Atlas Segmentation using
% Patch-wise label Fusion.
% USAGE:
% L = masPatchFusion(imgqry,supdata,supdata_lbl,opts,schspels)
% INPUTS:
% imgqry       - Intensity MRI wich is wanted to segment.
% supdata      - Is a 4-D Matrix size: [m x n x p x NS] wich holds a set of 
%                NS Intensity Atlases.
% supdata_lbl  - Is a 4-D Matrix size: [m x n x p x NS] wich holds a set of 
%                NS Labeled Atlases.
% opts         - Option strucuture with fields:
%                 alpha: patch radius (default 0)
%                 beta:  neighborhood radius (default 0)
%                 ss: Binary flag to perform patch selection by structure
%                  similarity (default false). 
%                 weighting_method: Method to compute the weights used in 
%                  patch label fusion:
%                  0 = All weights are set to 1 (default);
%                  1 = w are computed according to [1] (similarity measure)
%                  2 = w are computed according to [2] (sparse restriction)
%                 feat_method: Method to compute patch features:
%                  0 = Patch intensities are used as features (default).
%                  1 = Linear mapping of patch intensities is used as
%                      features. Mapping is automatically computed using
%                      Centered Kernel Alignment according to [3].
% schspels      - (optional) Holds the set of spels where the labels are
%                 estimated. If it is not specified, all elements in imgqry
%                 are labeled.
% 
% OUTPUT:   
% L             - Resulting segmentation on the specified
%                 spels.
% 
% REQUIREMENTS:
% SLEP package, available in http://yelab.net/software/SLEP/
% MLmat package, available in https://github.com/gcypds/MLmat
% 
%RELATED PAPERS:
% [1] Coupe, P., Manjon, J.V., Fonov, V., Pruessner, J., Robles, M., Collins
%    D.L., 2011. Patch-based segmentation using expert priors: application
%    to hippocampus and ventricle segmentation. Neuroimage 54, 904-954
% 
% [2] Tong, T.,Wolz, R., Hajnal, J.V., Rueckert, D., 2012. Segmentation of 
%    brain images via sparse patch representation. MICCAIWorkshop on 
%    Sparsity Techniques inMedical Imaging, Nice, France.
%
% [3] Cortes, C., Mohri, M., & Rostamizadeh, A. 2012. Algorithms for 
%   learning kernels based on centered alignment. The Journal of Machine 
%   Learning Research, 13(1), 795-828.
%
% Created on Wed Oct  7 11:18:05 2015
% Mauricio Orbes Arteaga - GCPDS
% David Cardenas Pena - GCPDS

if isempty(opts)
 alpha = 1;
 beta  = 2;
 ss  = false;
 weighting_method = 0;
 feat_method = 0;
else
 ss  = opts.ss;
 alpha = opts.alpha;
 beta = opts.beta;
 weighting_method = opts.weighting_method;
 feat_method = opts.feat_method;
end



if isempty(schspels)
    schspels=1:numel(imgqry);
end

[x,y,z]=ind2sub(size(imgqry),schspels);
L=zeros(size(imgqry));

for indc=1:numel(schspels)
  
  %Neighborhood contruction:
  Aph=getNeighborhood(imgqry,[],[x(indc) y(indc) z(indc)],alpha,0);
  [BTH,LBL]=getNeighborhood(supdata,supdata_lbl,[x(indc) y(indc) z(indc)],alpha,beta);    
  
  %Similarity structure selection
  if ss
    [BTH,LBL]=ssim(Aph,BTH,LBL,0.9);
  end
  
  %Feature computation
  switch feat_method
  %case 0: Patch intensities are the features.
  case 1
    BTH = zscore([BTH Aph]');
    Aph = BTH(end,:)';
    BTH = BTH(1:end-1,:)';
    A = kMetricLearningMahalanobis(BTH',pdist2(LBL',LBL')==0,LBL',0,false,false,[1e-4 1e-5]);
    BTH = A'*BTH;
    Aph = A'*Aph;
  end
    
  %Weights computation:
  switch weighting_method        
  case 0
    w = ones(1,size(BTH,2));
  case 1
    w = compute_weights_similarity(BTH,Aph);
  case 2
    w = compute_weights(BTH,Aph,0.01);
  end
  
  %Label fusion estimation:
  L(x(indc),y(indc),z(indc))=weightedVoting(w,LBL');
    
end