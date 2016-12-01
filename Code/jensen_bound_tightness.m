clc; clear all; close all;
l = 25;
r = 0.2;
k = 0:100;
P = exp(-l*pi*r^2).*(l*pi*r^2).^k./factorial(k);

m = l*pi*r^2; % mean of degree distribution


ind = 1;
for a = 0.01:0.01:1;
    
    e = 1;
    x = 1;
    tol = 1e-3;
    while e > tol && x ~= 0
        e = abs(x - (1/(l*pi*r^2))*sum( (a.*x.*k.^2).*P./(1 + a.*x.*k) ))
        x = x - 0.0001;
    end
    
    x_opt(ind) = x
    
    x_jensen(ind) = max(0, 1 - 1/(a*l*pi*r^2));
    ind = ind+1;
end

pl1 = plot(0.01:0.01:1,x_opt,'-b')
hold on
pl2 = plot(0.01:0.01:1,x_jensen,'--r')
axis([0 1 0 1])
grid on
legend([pl1 pl2],'Exact Solution', 'Approximation')
xlabel('Information spreading rate, \alpha')
ylabel('\Theta(\alpha)')

a_V = 0.1:0.1:0.9;
for i = 1:length(x_opt)
    rho_opt(i,:) = (k.*a_V(i).*x_opt(i))./(1+k.*a_V(i).*x_opt(i)) ;
    rho_jensen(i,:) = (k.*a_V(i).*x_jensen(i))./(1+k.*a_V(i).*x_jensen(i));
    loose_approximation(i) = max(0, 1 - 1/(a_V(i)*l*pi*r^2));
    c1 = (a_V(i)*x_jensen(i))/((1+m*a_V(i)*x_jensen(i))^2); % first coefficient of taylor series
    c2 = -2*((a_V(i)*x_jensen(i))^2)/((1+m*a_V(i)*x_jensen(i))^3); % second coefficient of taylor series
    taylor_approximation(i) = m*a_V(i)*x_jensen(i)/(1+m*a_V(i)*x_jensen(i)) - (((a_V(i)*x_jensen(i))^2)/((1 + m*a_V(i)*x_jensen(i))^3))*(m) + (((a_V(i)*x_jensen(i))^3)/((1 + m*a_V(i)*x_jensen(i))^4))*(m + 3*m^2);
end

figure;
plot(mean(rho_opt,2),'-b')
hold on
plot(mean(rho_jensen,2),'--r')
plot(loose_approximation, ':r')
plot(taylor_approximation, '-.r')
axis([1 9 0 1])


%% Taylor expansion of rho_k and the approximate expectation

for i = 1:length(a_V)
    c1 = (a_V(i)*x_jensen(i))/((1+a_V(i)*x_jensen(i))^2); % first coefficient of taylor series
    c2 = -2*((a_V(i)*x_jensen(i))^2)/((1+a_V(i)*x_jensen(i))^3); % second coefficient of taylor series
    rho_taylor(i,:) = (m*a_V(i)*x_jensen(i))/(1+ m*a_V(i)*x_jensen(i)) + c1*(a_V(i) - m) - 0.5*c2*(a_V(i) - m)^2; % Taylor series expansion at mean value
end


at = 1;
for k = 0:100
    g(k+1) = k*at/(1+k*at);
    g_approx(k+1) = m*at/(1+m*at) + ((at)/((1+m*at)^2))*(k - m) + 0.5*((-2*at^2)/((1 + m*at)^3))*(k - m)^2  + (1/6)*((6*at^3)/((1 + m*at)^4))*(k - m)^3;
end
 


figure
plot(g,'-b')
hold on
plot(g_approx,'--r')
plot(P)
axis([0 50 0 1])

a = 0.05;
k = 0:100;
ind_x = 1;
for x = 0:0.1:1;
    fp(ind_x) = (1/(l*pi*r^2))*sum( (a.*k.^2).*P./(1 + a.*x.*k) );
    ind_x = ind_x + 1;
end
plot(fp)
hold on;
plot(ones(1,length(0:0.1:1)))




