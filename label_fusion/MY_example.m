%Example MAS segmenation script
% Created on Wed Oct  7 15:04:41 2015
% Mauricio Orbes Arteaga - GCPDS

clear all
close all
clc

weighting_method=1; %similarity-based weighting
alpha=floor(3/2);   %Patch radius
beta=floor(5/2);    %Neighborhood radius

load sample_data/bw_db

%Segment subject Y4 using subjects Y5,L5 and Y6,L6:
imgqry = Y4;                       %Query image
supdata = cat(4,Y5,Y6);            %Intensity atlases
supdata_lbl = cat(4,L5,L6);        %Label atlases

mask = zeros(size(imgqry));
mask(80:100,80:100,80:100) = ones; %Mask to segment
indices = find(mask(:));           %Indices to segment

%Perform Multi Atlas Segmentation by Patch-wise Label Fusion
O = masPatchFusion(imgqry,supdata,supdata_lbl,weighting_method,beta,alpha,indices);

Accuracy = 100*sum(O(indices)==L4(indices))/numel(indices);