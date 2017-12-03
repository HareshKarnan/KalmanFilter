clear all;
close all;
clc
dt=0.1;
[t,y]=ode45(@dynamics,0.1:0.1:20,[10^5;-6000;2000],[],dt) 
plot(t,y(:,1))