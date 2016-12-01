clear all; clc; close all
theta = 4;
L = 20;
r = 1;


dist = @(x1,y1,x2,y2) sqrt((y2-y1).^2 + (x2-x1).^2);

degree = [];
for iter = 1:400
%     iter
    N = poissrnd(theta*(L + 4*r)^2);
    X_init = unifrnd(-L/2 - 2*r,L/2 + 2*r,N,2);
    
    cond = (X_init(:,1) >= -L/2 & X_init(:,1) <= L/2 ) & (X_init(:,2) >= -L/2 & X_init(:,2) <= L/2 ); 
    X = X_init(cond,:);

    if iter == 1
        scatter(X(:,1),X(:,2),30,'.');
    end
    
    ind_init = 1:length(X_init);
    ind = 1:length(X);
    for i = 1:length(X)
        self_index = find((X_init(:,1) == X(i,1)) & (X_init(:,2) == X(i,2)) ); 
        temp = dist(X(i,1),X(i,2),X_init(ind_init(ind_init~=self_index),1),X_init(ind_init(ind_init~=self_index),2)); % distances between ith node and all other nodes
        degree = [degree; sum(temp < r)];
    end
end
figure;
[b,x]= histnorm(degree,11);
bar(x,b)
hold on
xx = 0:50;
d = (exp(-theta*pi*r^2)*(theta *pi*r^2).^xx)./(factorial(xx));
plot(xx,d,'-r')