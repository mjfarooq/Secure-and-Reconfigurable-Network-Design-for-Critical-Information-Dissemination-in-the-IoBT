clear all; clc; close all
lambda1 = 3;
lambda2 = 1;
L = 20;
r1 = 1;
r2 = 2;

dist = @(x1,y1,x2,y2) sqrt((y2-y1).^2 + (x2-x1).^2);

degree = [];
for iter = 1:200
    N1 = poissrnd(lambda1*(L + 4*r1)^2); 
    N2 = poissrnd(lambda2*(L + 4*r2)^2);

    X_init = unifrnd(-L/2 - 2*r1,L/2 + 2*r1,N1,2);
    cond_x = (X_init(:,1) >= -L/2 & X_init(:,1) <= L/2 ) & (X_init(:,2) >= -L/2 & X_init(:,2) <= L/2 );
    
    Y_init = unifrnd(-L/2 - 2*r2,L/2 + 2*r2,N2,2);
    cond_y = (Y_init(:,1) >= -L/2 & Y_init(:,1) <= L/2 ) & (Y_init(:,2) >= -L/2 & Y_init(:,2) <= L/2 );
    
    X = X_init(cond_x,:); % coordinates for layer 1
    Y = Y_init(cond_y,:); % coordinates for layer 2
    
    if iter == 1
        scatter(X(:,1), X(:,2),30,'.');
        hold on
        scatter(Y(:,1), Y(:,2),80,'r.');
    end
    
    ind_init = 1:length(X_init);
    ind_init2 = 1:length(Y_init);
    ind = 1:length(X);
    ind2 = 1:length(Y);
    
    for i = 1:length(X)
        self_index = find((X_init(:,1) == X(i,1)) & (X_init(:,2) == X(i,2)) ); 
        temp = dist(X(i,1),X(i,2),X_init(ind_init(ind_init~=self_index),1),X_init(ind_init(ind_init~=self_index),2)); % distances between ith node and all other nodes
        temp1 = dist(X(i,1),X(i,2),Y_init(:,1),Y_init(:,2));
        degree = [degree; sum(temp < r1)+sum(temp1 < r1)];
    end
    
    for j = 1:length(Y)
        self_index = find((Y_init(:,1) == Y(j,1)) & (Y_init(:,2) == Y(j,2)) );
        temp2 = dist(Y(j,1),Y(j,2),Y_init(ind_init2(ind_init2~=self_index),1),Y_init(ind_init2(ind_init2~=self_index),2)); % distances between ith node and all other nodes
        %temp3 = dist(Y(j,1),Y(j,2),Y(ind2(ind2~=j),1),Y(ind2(ind2~=j),2));
        temp3 = dist(Y(j,1),Y(j,2),X_init(:,1),X_init(:,2));
        degree = [degree; sum(temp2 < r2) + sum(temp3 < r2)];
    end
    
end
figure
[b,x]= histnorm(degree,40);
bar(x,b)
hold on
xx = 0:80;
d = (exp(-lambda1*pi*r1^2)*(lambda1 *pi*r1^2).^xx)./(factorial(xx));
d1 = (exp(-lambda2*pi*r2^2)*(lambda2 *pi*r2^2).^xx)./(factorial(xx));
d2 = exp(-(lambda1*pi*r1^2 + lambda2*pi*r2^2) )*((lambda1 *pi*r1^2 + lambda2*pi*r2^2 ).^xx)./(factorial(xx));
da = exp(-(lambda1+lambda2)*pi*r1^2 )*(((lambda1+lambda2) *pi*r1^2 ).^xx)./(factorial(xx));
db = exp(-(lambda1+lambda2)*pi*r2^2 )*(((lambda1+lambda2) *pi*r2^2 ).^xx)./(factorial(xx));
d3 = (lambda1/(lambda1 + lambda2))*da + (lambda2/(lambda1 + lambda2))*db;
%plot(xx,d,'-r')
%plot(xx,d1,'-blue')
%plot(xx,d2,'-g')
plot(xx,d3,'-r')
