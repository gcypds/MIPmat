function [B,L]=ssim(X,Y,L,ths)
% ssim computes the structural similarity between two matrix datasets:
% USAGE:
% [Z,H] = ssim(X,Y,L,th): returns a n-by-mz matrix Z with the subset of 
%       columns from the n-by-my matrix Y (my>=mz) and an m-by-mz matrix H
%       from the m-by-my matrix L with a structual similarity between each 
%       observation in Y and the column vector X larger than the real-valued
%       threshold th. Columns of Y correspond to observations, rows 
%       correspond to variables.
%
% https://en.wikipedia.org/wiki/Structural_similarity  
% 
% Created on Wed Oct  7 10:53:17 2015
% Mauricio Orbes Arteaga - GCPDS
% David Cardenas Pena - GCPDS

if isempty(Y) 
  Y = X;
end

ui=mean(X);
uj=mean(Y);

term1= (2*ui*uj)./(ui^2+uj.^2);
         
si=std(X).^2;
ssj=std(Y).^2;
sisj=abs(1/(size(Y,1)-1)*sum(repmat((X-ui),[1,size(Y,2)]).*(Y-repmat(uj,[size(Y,1) 1]))));
         
term2 = 2*sisj./(si+ssj);
ss = term1.*term2;
         
ssc = ss>ths;
         
        if sum(ssc)==0;
           ssc=(ss>(max(ss)*ths));
        end
        
        B=Y(:,ssc);
        L=L(ssc);
         
end