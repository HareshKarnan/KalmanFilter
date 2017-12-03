function snextt = dynamics(t,s,dt)
rho_0 = 3.4e-3;g = 32.2;
k_rho = 22000;
snextt=zeros(3,1);
snextt(1,1)=s(2);
snextt(2,1)=-g+rho_0*exp(-s(1)/k_rho)*s(2)^2/(2*s(3));
snextt(3,1)=0;
end