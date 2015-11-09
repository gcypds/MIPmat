function L = weightedVoting(W,LBL)

lb=unique(LBL);
   
pooling = zeros(numel(lb),size(W,2));
for i=1:numel(lb)
  pooling(i,:)=sum(W(LBL==lb(i),:),1)./sum(W,1);
end
          
[~,ind] = max(pooling,[],1);
L=lb(ind);