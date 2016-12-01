

clc; clear all; close all;
lambda = 3*[2 1]; % Intensity of nodes of each type
r = 1*[0.5 1]; % Transmission range of nodes in each network layer
M = length(lambda); % Number of layers
Lambda = sum(lambda);

Mean_degree = sum(pi*lambda.*r.^2);

l1 = lambda(1)*pi*r(1)^2;
l2 = lambda(2)*pi*r(2)^2;
k = 0:200;
P1 = exp(-l1).*(l1.^k).*(1./factorial(k));
P2 = exp(-l2).*(l2.^k).*(1./factorial(k));
P_bar = (1/Lambda)* (   lambda(1)*exp(-lambda(1)*pi*r(1)^2) * (lambda(1)*pi*r(1)^2).^k .* (1./ factorial(k))  +  lambda(2)*exp(-lambda(2)*pi*r(2)^2) * (lambda(2)*pi*r(2)^2).^k .* (1./ factorial(k))  );

attack_V = 1e-6:0.01:1;
for i = 1:length(attack_V)
    attack_probability = attack_V(i);
    mu = 1 - attack_probability;
    % t = 0.3;`
    % S = ( exp(-l1) * mu* t/l1   ) * k.^2 .* l1.^k .* (1./factorial(k)) .* (1./(1 + mu*t*k));
    % sum(S)
    
    a =@(t) 1 + 1/(mu*t);
    b =@(t) 2 + 1/(mu*t);
    h1 = @(t) real(   (exp(-l1)/l1)* -(-l1)^(-1/(mu*t)) * (gammainc(-l1,a(t))*gamma(a(t))  -  gammainc(-l1,b(t))*gamma(b(t))) ) - t;
    h2 = @(t) real(   (exp(-l2)/l2)* -(-l2)^(-1/(mu*t)) * (gammainc(-l2,a(t))*gamma(a(t))  -  gammainc(-l2,b(t))*gamma(b(t))) ) - t;
    h_bar =@(t) real(  (lambda(1)/Lambda) *(exp(-Lambda*pi*r(1)^2)/Lambda*pi*r(1)^2)* -(-Lambda*pi*r(1)^2)^(-1/(mu*t)) * (gammainc(-Lambda*pi*r(1)^2,a(t),'lower')*gamma(a(t))  -  gammainc(-Lambda*pi*r(1)^2,b(t),'lower')*gamma(b(t)))  +  (lambda(2)/Lambda) *(exp(-Lambda*pi*r(2)^2)/Lambda*pi*r(2)^2)* -(-Lambda*pi*r(2)^2)^(-1/(mu*t)) * (gammainc(-Lambda*pi*r(2)^2,a(t),'lower')*gamma(a(t))  -  gammainc(-Lambda*pi*r(2)^2,b(t),'lower')*gamma(b(t))) ) - t;
    theta1 = fsolve(h1,1);
    theta2 = fsolve(h2,1);
    theta_bar(i) = fsolve(h_bar,1)
    theta_bar_approx(i) = max(0, 1 - (1/(mu*sum(lambda.*pi.*r.^2)))  );
    
    rho_k1 = mu.*k.*theta1 ./ (1 + mu.*k.*theta1);
    rho_k2 = mu.*k.*theta2 ./ (1 + mu.*k.*theta2);
    rho_kbar = mu.*k.*theta_bar(i) ./ (1 + mu.*k.*theta_bar(i));
    rho_kbar_approx = mu.*k.*theta_bar_approx(i) ./ (1 + mu.*k.*theta_bar_approx(i));
    
    rho1(i) = sum(rho_k1.*P1);
    rho2(i) = sum(rho_k2.*P2);
    rho_bar(i) = sum(rho_kbar.*P_bar);
    rho_approx(i) = sum(rho_kbar_approx.*P_bar);
    
end
plot(attack_V,theta_bar,'--g')
hold on
plot(attack_V,theta_bar_approx,'-m')

figure;
plot(k, rho_k1, '--')
hold on
plot(k, rho_k2, '-r')
plot(k, rho_kbar, '-g')
plot(k, rho_kbar_approx, '-m')



figure
plot(attack_V, rho1, '--')
hold on
plot(attack_V, rho2, '-r')
plot(attack_V, rho_bar, '-g')
plot(attack_V, rho_approx, '-m')
axis([0 1 0 1])

% g = [];
% for t = 0:0.01:1
%     S = ( exp(-l1) * mu* t/l1   ) * k.^2 .* l1.^k .* (1./factorial(k)) .* (1./(1 + mu*t*k));

%     g = [g sum(S) - t];
% end
% plot(g)

