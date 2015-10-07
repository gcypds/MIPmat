function L=masPatchFusion(imgqry,supdata,supdata_lbl,band_method,NG,PTH,schspels)
% Generate the label segmentation by Multi-Atlas Segmentation using
% Patch-wise label Fusion.
% USAGE:
% L = Patch_fusion_MAS(imgqry,supdata,supdata_lbl,band_method,NG,PTH,schspels)
% INPUTS:
% imgqry       - Intensity MRI wich is wanted to segment.
% supdata      - Is a 4-D Matrix size: [m x n x p x NS] wich holds a set of 
%                NS Intensity Atlases.
% supdata_lbl  - Is a 4-D Matrix size: [m x n x p x NS] wich holds a set of 
%                NS Labeled Atlases.
% band_method  - Type of method used to compute the weights used in patch
%                label fusion methods:
%                  0 = All weights are set to 1 (same value);
%                  1 = w are computed according to [1] (similarity measure);
%                  2 = w are computed according to [2] (sparse restriction);
% NG            - Radius of Search neighborhood
% PTH           - Radius of Patch
% schspels      - (optional) Holds the set of spels where the labels are
%                 estimated. If it is not specified, all elements in imgqry
%                 are labeled.
% 
% OUTPUT:   
% L             - Resulting segmentation on the specified
%                 spels(schspels ).
% 
% REQUIREMENTS:
% SLEP package, available in http://yelab.net/software/SLEP/
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
% Created on Wed Oct  7 11:18:05 2015
% Mauricio Orbes Arteaga - GCPDS


if isempty(schspels)
    schspels=1:numel(imgqry);
end

[x,y,z]=ind2sub(size(imgqry),schspels);
L=zeros(size(imgqry));

for indc=1:numel(schspels)
  
  %Neighborhood contruction:
  Aph=getNeighborhood(imgqry,[],[x(indc) y(indc) z(indc)],PTH,0);
  [BTH,LBL]=getNeighborhood(supdata,supdata_lbl,[x(indc) y(indc) z(indc)],PTH,NG);    
  [BTH,LBL]=ssim(Aph,BTH,LBL,0.9);
    
  %Weights computation:
  switch band_method        
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