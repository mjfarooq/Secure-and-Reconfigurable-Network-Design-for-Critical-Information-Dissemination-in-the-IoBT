

clc; clear all; close all; format long
lambda = 1*[15 5]; % Intensity of nodes of each type
r = 0.3*[1 2]; % Transmission range of nodes in each network layer
M = length(lambda); % Number of layers
Lambda = sum(lambda);

L = 10;
N1 = poissrnd(lambda(1)*L^2); 
N2 = poissrnd(lambda(2)*L^2);
X = unifrnd(-L/2,L/2,N1,2);
Y = unifrnd(-L/2,L/2,N2,2);
scatter(X(:,1), X(:,2),30,'.');
hold on
circles(0,0,r(1),'edgecolor','blue','facecolor', 'none');
scatter(Y(:,1), Y(:,2),80,'r.');
circles(1,1,r(2),'edgecolor','r','facecolor', 'none');

Mean_degree = sum(pi*lambda.*r.^2);

l1 = lambda(1)*pi*r(1)^2;
l2 = lambda(2)*pi*r(2)^2;
k = 0:80;
P1 = exp(-l1).*(l1.^k).*(1./factorial(k));
P2 = exp(-l2).*(l2.^k).*(1./factorial(k));
P_bar = (1/Lambda)* (   lambda(1)*exp(-Lambda*pi*r(1)^2) * (Lambda*pi*r(1)^2).^k .* (1./ factorial(k))  +  lambda(2)*exp(-Lambda*pi*r(2)^2) * (Lambda*pi*r(2)^2).^k .* (1./ factorial(k))  );

attack_V = 0.001:0.01:1;
for i = 1:length(attack_V)
    attack_probability = attack_V(i)
    mu = 1 - attack_probability;
    % t = 0.3;`
    % S = ( exp(-l1) * mu* t/l1   ) * k.^2 .* l1.^k .* (1./factorial(k)) .* (1./(1 + mu*t*k));
    % sum(S)
    
    a =@(t) 1 + 1/(mu*t);
    b =@(t) 2 + 1/(mu*t);
    h1 = @(t) real(   (exp(-l1)/l1)* -(-l1)^(-1/(mu*t)) * (gammainc(-l1,a(t),'lower')*gamma(a(t))  -  gammainc(-l1,b(t),'lower')*gamma(b(t))) ) - t;
    h2 = @(t) real(   (exp(-l2)/l2)* -(-l2)^(-1/(mu*t)) * (gammainc(-l2,a(t),'lower')*gamma(a(t))  -  gammainc(-l2,b(t),'lower')*gamma(b(t))) ) - t;
    h_bar =@(t) real(  (lambda(1)/Lambda) *(exp(-Lambda*pi*r(1)^2))* -(-Lambda*pi*r(1)^2)^(-1/(mu*t)) * (gammainc(-Lambda*pi*r(1)^2,a(t),'lower')*gamma(a(t))  -  gammainc(-Lambda*pi*r(1)^2,b(t),'lower')*gamma(b(t)))  +  (lambda(2)/Lambda) *(exp(-Lambda*pi*r(2)^2))* -(-Lambda*pi*r(2)^2)^(-1/(mu*t)) * (gammainc(-Lambda*pi*r(2)^2,a(t),'lower')*gamma(a(t))  -  gammainc(-Lambda*pi*r(2)^2,b(t),'lower')*gamma(b(t))) ) - t;
    %theta1(i) = fsolve(h1,1);
    theta1_approx(i) = max(0, 1 - (1/(mu*lambda(1)*pi*r(1)^2))  );
    %theta2(i) = fsolve(h2,1);
    theta2_approx(i) = max(0, 1 - (1/(mu*lambda(2)*pi*r(2)^2))  );    
    %theta_bar(i) = fsolve(h_bar,1);
    theta_bar_approx(i) = max(0, 1 - (1/(mu*sum(lambda.*pi.*r.^2)))  );
    
    rho_k1 = mu.*k.*theta1_approx(i) ./ (1 + mu.*k.*theta1_approx(i));
    rho_k2 = mu.*k.*theta2_approx(i) ./ (1 + mu.*k.*theta2_approx(i));
    rho_kbar = mu.*k.*theta_bar_approx(i) ./ (1 + mu.*k.*theta_bar_approx(i));
    
    rho1(i) = sum(rho_k1.*P1);
    rho2(i) = sum(rho_k2.*P2);
    rho_bar(i) = sum(rho_kbar.*P_bar);
    
end

cutoff_1 = 1 - 1/(1+lambda(1)*pi*r(1)^2);
cutoff_2 = 1- 1/(1+lambda(2)*pi*r(2)^2);
cutoff_bar = 1 - 1/(1+lambda(1)*pi*r(1)^2 + lambda(2)*pi*r(2)^2);

figure;
%plot(attack_V,theta1,'-blue');
plot(attack_V,theta1_approx,':blue', 'Linewidth',1.2);
hold on;
%plot(attack_V,theta2,'-red');
plot(attack_V,theta2_approx,':red', 'Linewidth',1.2);
%plot(attack_V,theta_bar,'--g');
plot(attack_V,theta_bar_approx,':m', 'Linewidth',1.2);

figure;
aa = plot(attack_V, rho1, '-blue', 'Linewidth',1.2)
hold on;
plot(attack_V, theta1_approx, ':blue', 'Linewidth',1.2)
bb = plot(attack_V, rho2, '-r', 'Linewidth',1.2)
plot(attack_V, theta2_approx, ':r', 'Linewidth',1.2)
cc = plot(attack_V, rho_bar, '-m', 'Linewidth',1.2)
plot(attack_V, theta_bar_approx, ':m', 'Linewidth',1.2)
%legend([aa bb cc], 'Layer 1', 'Layer 2', 'Combined');
xlabel('Probability of attack')
ylabel('Information prevalence')
line([cutoff_1 cutoff_1],[0 1])
line([cutoff_2 cutoff_2],[0 1])
line([cutoff_bar cutoff_bar],[0 1])
grid on;

figure;
aa2 = plot(attack_V, rho1*lambda(1), '-b', 'Linewidth',1.2)
hold on;
bb2 = plot(attack_V, rho2*lambda(2), '-r', 'Linewidth',1.2)
cc2 = plot(attack_V, rho_bar*Lambda, '-m', 'Linewidth',1.2)
%legend([aa2 bb2 cc2], 'Layer 1', 'Layer 2', 'Combined');
xlabel('Probability of attack')
ylabel('Number of Informed Nodes')
line([cutoff_1 cutoff_1],[0 1])
line([cutoff_2 cutoff_2],[0 1])
line([cutoff_bar cutoff_bar],[0 1])
grid on;


figure;
plot(P_bar)

format short




