% Implementation of UKF
clear all;clc;close all;warning off;
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
w_m = [l/(l+n) 0.5/(l+n)+zeros(1,2*n)];
w_c = [l/(l+n)+(1-alpha^2+beta) 0.5/(l+n)+zeros(1,2*n)];
j=1;
x_t=[normrnd(10^5,sqrt(500));normrnd(-6000,sqrt(2*10^4));normrnd(2000,sqrt(2.5*10^5))];
tic
for time =t
    % find the sigma points
    X=[u_0,u_0*ones(1,n)+sqrt(n+l)*chol(P_0)',u_0*ones(1,n)-sqrt(n+l)*chol(P_0)'];
    % propagate the points into the dynamics and find mean and covariance
    [Xp(:,:,j),s_u(:,j),P_t(:,:,j)] = propxsi(1,X,n,dt,w_m,w_c);
    % add noise term
    P_t(:,:,j) = P_t(:,:,j)+R_t.*dt^2;
    % measurements
    Xpp=[s_u(:,j),s_u(:,j)*ones(1,n)+sqrt(n+l)*chol(P_t(:,:,j))',s_u(:,j)*ones(1,n)-sqrt(n+l)*chol(P_t(:,:,j))']; %resample
    [Ypp,s_uy(j),P_ty(j)] = propxsi(2,Xpp,n,dt,w_m,w_c);
    P_ty(j) = P_ty(j)+Q_t;
    
    P_xy=zeros(3,1);
    for i=1:2*n+1
        P_xy = P_xy+w_c(i)*(Xpp(:,i)-s_u(:,j))*(Ypp(i)-s_uy(j))';
    end
    % find the kalman gain
    K = P_xy*P_ty(j)^-1;
    
    %dynamics truth simulation
    x_t(:,j+1)= truthfunc(x_t(:,j),dt);
    % measurement value
    zm(j) = x_t(1,j+1)+normrnd(0,sqrt(100));
    
    % do the kalman update
    u_upd(:,j) = s_u(:,j)+K*(zm(j)-s_uy(:,j));
    P_upd(:,:,j) = P_t(:,:,j)-K*P_ty(j)*K';
    j=j+1;
end
toc
plot(t,u_upd(1,:),'-r',t,zm,'-b')
legend('measured state','filtered state')
xlabel('Time (sec)');
ylabel('x1 (feet)');

function [Xsi,s_u,P_t] = propxsi(val,X,n,dt,w_m,w_c)
if val==1
    % function to propogate the points through dynamics
    % and find the predicted mean and covariance.
    Xsi=[];s_u=zeros(n,1);P_t=zeros(n,n);
    % propogate the points
    for i=1:2*n+1
        Xsi = [Xsi,gmeanfunc(X(:,i),dt)];
    end
    % find the mean of the propogated points
    for i=1:2*n+1
        s_u = s_u + w_m(i)*Xsi(:,i);
    end
    % find covariance of the propogated points
    for i=1:2*n+1
        P_t = P_t+w_c(i)*(Xsi(:,i)-s_u)*(Xsi(:,i)-s_u)';
    end
end
if val==2
    % function to propogate the points through dynamics
    % and find the predicted mean and covariance.
    Xsi=[]; s_u=0; P_t=0;
    % propogate the points
    for i=1:2*n+1
        Xsi = [Xsi,hmeanfunc(X(:,i),dt)];
    end
    % find the mean of the propogated points
    for i=1:2*n+1
        s_u = s_u + w_m(i)*Xsi(:,i);
    end
    % find covariance of the propogated points
    for i=1:2*n+1
        P_t = P_t+w_c(i)*(Xsi(:,i)-s_u)*(Xsi(:,i)-s_u)';
    end
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

function y = hmeanfunc(s,dt)
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