
clc; clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R1_min = 0.1; % minimum range is 100 meter 
R1_max = 2; % max range is 500 meter 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R2_min = 0.01; % minimimum range is 10 meter
R2_max = 0.8; % max range is 100 meters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L_min = 1;%0.1; % minimum deployment 0.02 per km2
L_max = 15;%10; % maximum deployment 0.1 per km^2 
p_min = 0;
p_max = 0.4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T1 = 0.6;
T2 = 0.6;
Tc = 0.8;

w1 = 100;
w2 = 50;
c = 100;
options = optimoptions('fmincon','Algorithm','sqp');

eta = 4; % Path-loss exponent

ind = 0;
delta_V = 0.1:0.1:0.9;
for delta = 0.1:0.1:0.9 % threat level
    ind = ind+1;
    alpha = 1 - delta;
  
    % Variables are x = [p  lambda r1 r2]
    fun = @(x) (w1 * x(1)*x(2)  +  w2 * (1 - x(1)) * x(2) + c * (  x(1)*x(2)*x(3)^eta  + x(2) * x(4)^eta) );
    lb = [p_min, L_min, R1_min,R2_min];
    ub = [p_max, L_max, R1_max,R2_max];
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    x0 = [0.1, 5, 0.5, 0.1];
    
    [x_opt(ind,:) fval(ind) flag(ind)] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,@(x) const(x,alpha,T1,T2,Tc),options)
end

w1 * x_opt(:,1).*x_opt(:,2)  +  w2 * (1 - x_opt(:,1)) .* x_opt(:,2)
c * (  x_opt(:,1).*x_opt(:,2).*x_opt(:,3).^eta  + x_opt(:,2) .* x_opt(:,4).^eta)

flag

figure;
pp1 = plot(delta_V,1000*x_opt(:,3),'-*r', 'Linewidth', 1.2); % Radius
hold on
pp2 = plot(delta_V,1000*x_opt(:,4),'--ob', 'Linewidth', 1.2); % Radius
grid on
xlabel('Threat level, $\delta$','Interpreter','latex')
ylabel('Transmission Range (m)')
legend([pp1,pp2],'r_1','r_2')
title('$T_1 = 0.6, T_2 = 0.6, T_{c} = 0.8$','Interpreter','latex')

figure;
pp3 = plot(delta_V,x_opt(:,1).*x_opt(:,2),'-*r', 'Linewidth', 1.2) % Intensity
hold on
pp4 = plot(delta_V,x_opt(:,2),'--ob', 'Linewidth', 1.2) % Intensity
grid on
xlabel('Threat level, $\delta$','Interpreter','latex')
ylabel('Device Intensity (devices/km^2)')
legend([pp3,pp4],'\lambda_1', '\lambda_2')
title('$T_1 = 0.6, T_2 = 0.6, T_{c} = 0.8$','Interpreter','latex')

figure
plot(delta_V,fval, '-xb' ,'Linewidth', 1.2) % Plot of the Cost function
grid on
xlabel('Threat level, $\delta$','Interpreter','latex')
ylabel('Cost Function')
title('$T_1 = 0.6, T_2 = 0.6, T_{c} = 0.8$','Interpreter','latex')
