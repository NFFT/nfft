function x=rand_subsample(X,M)

I=randperm(size(X,1));
x=X(I(1:M),:);