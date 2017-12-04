% Implementation of UKF
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
n=size(u_0,1);
% find lambda parameter
alpha = 10^-2;kappa = 0;beta=2;
l=alpha^2*(n+kappa)-n;
% weights of the UKF points
w_m = [l/(l+n),0.5/(l+n)+zeros(1,2*n)];
w_c = [l/(l+n)+(1-alpha^2+beta),0.5/(l+n)+zeros(1,2*n)];
j=1;
x_t=[normrnd(10^5,sqrt(500));normrnd(-6000,sqrt(2*10^4));normrnd(2000,sqrt(2.5*10^5))];
tic
for time=t
    % find the sigma points
    if time == 0.1
    X = [u_0,u_0*ones(1,n)+sqrt(n+l)*chol(P_0)',u_0*ones(1,n)-sqrt(n+l)*chol(P_0)'];
    else
    X = [s_u(:,j-1),s_u(:,j-1)*ones(1,n)+sqrt(n+l)*chol(P_t(:,:,j-1),'lower'),s_u(:,j-1)*ones(1,n)-sqrt(n+l)*chol(P_t(:,:,j-1),'lower')];
    end
    % propagate the points into the dynamics and find mean and covariance
    [Xp(:,:,j),u_b(:,j),P_b(:,:,j)] = propxsi(1,X,n,dt,w_m,w_c);
   
    % measurements
    Xpp=[u_b(:,j),u_b(:,j)*ones(1,n)+sqrt(n+l)*chol(P_b(:,:,j),'lower'),u_b(:,j)*ones(1,n)-sqrt(n+l)*chol(P_b(:,:,j),'lower')]; %resample
    % propagate the measurement equation and find mean and covariance
    [Zpp,z_u(j),P_z(j)] = propxsi(2,Xpp,n,dt,w_m,w_c);
    
    P_xy=zeros(3,1);
    for i=1:2*n+1
        P_xy = P_xy+w_c(i)*(Xpp(:,i)-u_b(:,j))*(Zpp(i)-z_u(j))';
    end
    % find the kalman gain
    K = P_xy*P_z(j)^-1;
    
    %dynamics truth simulation
    x_t(:,j+1)= truthfunc(x_t(:,j),dt);
    
    % measurement value
    zm(j) = x_t(1,j+1)+normrnd(0,sqrt(100));
    
    % do the kalman update
    s_u(:,j) = u_b(:,j)+K*(zm(j)-z_u(j));
    P_t(:,:,j) = P_b(:,:,j)-K*P_xy';
    j=j+1;
end
toc
% plot(t,s_u(2,:),'-r',t,x_t(2,2:end),'-b')
% legend('filtered','truth','measured')
plot(t,zm,'-b',t,s_u(1,:),'-r',t,u_b(1,:),'-g')
% plot(t,u_b(1,:))
legend('measured state','filtered state','Dynamics')
xlabel('Time (sec)');
ylabel('x1 (feet)');

function [Xp,s_u,P_t] = propxsi(val,X,n,dt,w_m,w_c)
R_t = [0 0 0;0 2 0;0 0 0];
Q_t = 100;
if val==1
    % function to propogate the points through dynamics
    % and find the predicted mean and covariance.
    Xp=[];s_u=zeros(n,1);P_t=zeros(n,n);
    % propogate the points
    for i=1:2*n+1
        Xp = [Xp,gmeanfunc(X(:,i),dt)];
    end
    % find the mean of the propogated points
    for i=1:2*n+1
        s_u = s_u + w_m(i)*Xp(:,i);
    end
    % find covariance of the propogated points
    for i=1:2*n+1
        P_t = P_t+w_c(i)*(Xp(:,i)-s_u)*(Xp(:,i)-s_u)';
    end
    % add the noise term
    P_t = P_t+R_t.*dt^2;
end
if val==2
    % function to propogate the points through dynamics
    % and find the predicted mean and covariance.
    Xp=[]; s_u=0; P_t=0;
    % propogate the points
    for i=1:2*n+1
        Xp = [Xp,hmeanfunc(X(:,i))];
    end
    % find the mean of the propogated points
    for i=1:2*n+1
        s_u = s_u + w_m(i)*Xp(:,i);
    end
    % find covariance of the propogated points
    for i=1:2*n+1
        P_t = P_t+w_c(i)*(Xp(:,i)-s_u)*(Xp(:,i)-s_u)';
    end
    P_t = P_t + dt^2.*Q_t;
end

end

function snext = gmeanfunc(s,dt)
% function to propagate the dynamics
rho_0 = 3.4e-3;g = 32.2;
k_rho = 22000;
snext=zeros(3,1);
snext(1,1) = s(1) + s(2)*dt;
snext(2,1) = s(2) + dt*(-g+rho_0*exp(-s(1)/k_rho)*s(2)^2/(2*s(3)));
snext(3,1) = s(3);
end

function y = hmeanfunc(s)
y=s(1);
end

function snextt = truthfunc(s,dt)
rho_0 = 3.4e-3;g = 32.2;
k_rho = 22000;
snextt=zeros(3,1);
snextt(1,1)=s(1) + s(2)*dt;
snextt(2,1)=s(2) + dt*(-g+rho_0*exp(-s(1)/k_rho)*s(2)^2/(2*s(3)) + normrnd(0,sqrt(2)));
snextt(3,1)=s(3);
end