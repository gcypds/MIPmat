function L = weightedVoting(W,LBL)

lb=unique(LBL);
   
pooling = zeros(1,numel(lb));
for i=1:numel(lb)
  pooling(i)=sum(W(LBL==lb(i)))./sum(W);
end
          
[~,ind] = max(pooling);
L=lb(ind);