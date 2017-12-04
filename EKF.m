% IMPLEMENTATION OF EXTENDED KALMAN FILTER
clear all;clc;close all;warning off;format bank;
rho_0 = 3.4e-3;
g = 32.2;
k_rho = 22000;
P_0 = diag([500 2*10^4 2.5*10^5]);
u_0 = [10^5;-6000;2000];
R_t = [0 0 0;0 2 0;0 0 0];
Q_t = 100;
H_t = [1 0 0];
tf=20;
dt=0.1;
t = 0.1:dt:tf;
%calculate the Jacobian matrix
% syms x y z real
% G = jacobian([y*dt+x;y+dt*(-g+rho_0*exp(-x/k_rho)*y^2/(2*z));z],[x;y;z]);

figure
j=1;
%%
%true value of the state itself- sample from the given distribution
x_t=[normrnd(10^5,sqrt(500));normrnd(-6000,sqrt(2*10^4));normrnd(2000,sqrt(2.5*10^5))];
tic
for time=t
    % find the mean prediction
    if time==0.1
    s_u(:,j) = gmeanfunc(u_0,dt)
    else
    s_u(:,j) = gmeanfunc(s_u(:,j-1),dt);
    end
    % find the covariance prediction
    G_t(:,:,j) = jacobbi(s_u(1,j),s_u(2,j),s_u(3,j));
    if time==0.1
    P_t(:,:,j) = G_t(:,:,j)*P_0*G_t(:,:,j)' + dt^2*R_t;
    else
    P_t(:,:,j) = G_t(:,:,j)*P_t(:,:,j-1)*G_t(:,:,j)' + dt^2*R_t; 
    end
    % find the kalman gain matrix
    K_t(:,:,j) = P_t(:,:,j)*H_t'*(H_t*P_t(:,:,j)*H_t'+ Q_t)^-1;
    
    %dynamics truth simulation
    x_t(:,j+1)= truthfunc(x_t(:,j),dt);
    % measurement value
    zm(j) = x_t(1,j+1)+normrnd(0,sqrt(100));
    
    % kalman update step
    s_u(:,j) = s_u(:,j)+K_t(:,:,j)*(zm(j)-s_u(1,j));
    P_t(:,:,j) = (eye(3)-K_t(:,:,j)*H_t)*P_t(:,:,j);
    if time==5
    s_u(:,j)
    P_t(:,:,j)
    end
    
    j=j+1;
end
toc
% plot(0.1:dt:tf,s_u(2,:),0.1:dt:tf,x_t(2,2:end));
plot(0.1:dt:tf,zm,'-b',0.1:dt:tf,s_u(1,:),'-r')
legend('measured state','filtered state')
xlabel('Time (sec)');
ylabel('x1 (feet)');

function snext = gmeanfunc(s,dt)
% euler integration method
rho_0 = 3.4e-3;g = 32.2;
k_rho = 22000;
snext=zeros(3,1);
snext(1,1) = s(1) + s(2)*dt;
snext(2,1) = s(2) + dt*(-g+rho_0*exp(-s(1)/k_rho)*s(2)^2/(2*s(3)));
snext(3,1) = s(3);
end

function snextt = truthfunc(s,dt)
rho_0 = 3.4e-3;g = 32.2;
k_rho = 22000;
snextt=zeros(3,1);
snextt(1,1)=s(1) + s(2)*dt;
snextt(2,1)=s(2) + dt*(-g+rho_0*exp(-s(1)/k_rho)*s(2)^2/(2*s(3)) + normrnd(0,sqrt(2)));
snextt(3,1)=s(3);
end

function G_t = jacobbi(x,y,z,dt)
rho_0 = 3.4e-3;
g = 32.2;
k_rho = 22000;
G_t = [                                      1,                               1/10,                                    0;
     -(17*y^2*exp(-x/22000))/(2200000000*z), (17*y*exp(-x/22000))/(50000*z) + 1, -(17*y^2*exp(-x/22000))/(100000*z^2);
                                  0,                                  0,                                    1];
end