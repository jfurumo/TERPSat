%John Furumo
%UMD TERPSat
%Orbit Decay Model
%11/12/2020

clear;
clf;
clc;
close all;

format long;

%% NOTES
%simple exponential model for atmospheric density
%circular orbit assumption
%add elliptical orbits - drag paradox

%cannot used drag acceleration to lower orbit velocity over time because
%this translates to higher radius using vis-viva equation

%drag causes orbit radius to decrease and velocity to increase...

%% CONSTANTS
mu = 3.986E+05; %Earth gravitational parameter (km^3/s^2)
R = 6378; %Earth mean equatorial radius (km)
rho_0 = 1.226; %Earth surface atmosphere density (kg/m^3)
%h_s=7524/1000; %Earth atmosphere scale height (km)
h_s=8.5; %Earth atmosphere scale height (km)
cD = 2; %cubesat drag coefficient

%% CUBESAT PARAMETERS
%preallocate memory
A = zeros(3,5);
beta = zeros(3,5);

%1U
l = 0.1; %length (m)
w = 0.1; %width (m)
h = 0.1; %height (m)
A(1,1) = l*w; %small face area (m^2)
A(1,2) = l*h; %large face area (m^2)
A(1,3) = sqrt(l^2+w^2)*h; %2-face (small) projected area (m^2)
A(1,4) = sqrt(l^2+h^2)*w; %2-face (large) projected area (m^2)
A(1,5) = NaN; %3-face projected area (m^2)
m(1) = 1.33; %mass (kg)
for i=1:5
    beta(1,i) = m(1)/(A(1,i)*cD); %ballistic coefficient (kg/m^2)
end

%2U
l = 0.1; %length (m)
w = 0.1; %width (m)
h = 0.2; %height (m)
A(2,1) = l*w; %small face area (m^2)
A(2,2) = l*h; %large face area (m^2)
A(2,3) = sqrt(l^2+w^2)*h; %2-face (small) projected area (m^2)
A(2,4) = sqrt(l^2+h^2)*w; %2-face (large) projected area (m^2)
A(2,5) = NaN; %3-face projected area (m^2)
m(2) = 2.66; %mass (kg)
for i=1:5
    beta(2,i) = m(2)/(A(2,i)*cD); %ballistic coefficient (kg/m^2)
end


%3U
l = 0.1; %length (m)
w = 0.1; %width (m)
h = 0.3; %height (m)
A(3,1) = l*w; %small face area (m^2)
A(3,2) = l*h; %large face area (m^2)
A(3,3) = sqrt(l^2+w^2)*h; %2-face (small) projected area (m^2)
A(3,4) = sqrt(l^2+h^2)*w; %2-face (large) projected area (m^2)
A(3,5) = NaN; %3-face projected area (m^2)
m(3) = 4; %mass (kg)
for i=1:5
    beta(3,i) = m(3)/(A(3,i)*cD); %ballistic coefficient (kg/m^2)
end

%% ORBIT DECAY MODEL
% n = 101; %number of data points
% alt = linspace(500,0,n); %range of orbit altitudes (km)
% rho = rho_0*exp(-(alt)./h_s)
%a = R + alt; %range of orbit semi-major axes (km)
%P=2*pi*sqrt(a.^3/mu); %orbit period (s)

% figure('name','density as function of altitude')
% plot(alt,rho)
% figure('name','altitude as function of density')
% plot(rho,alt)

%%
beta = beta(1,1)
n = 1; %dataset index
alt = 300; %initial orbit altitude (km)
r(1) = alt+R;
%rho(n) = rho_0*exp(-(alt)/h_s) %atmospheric density at initial altitude (kg/m^3)
%v = zeros(1,n)';
%v(1) = sqrt(mu/(r(1))); %initial orbit velocity (km/s)
%a_D(n) = (0.5*rho(n)*v(n)^2)/beta
%t_step = 30*86400; %time step (s)

%t = 0;

for n=2:10000
%while alt > 0
    %while v(n) > 0
%         %Change in energy  due to drag method
%         rho(n-1) = rho_0*exp(-(r(n-1)-R)/h_s);
%         a_D(n-1) = (0.5*rho(n-1)*v(n-1)^2)/beta;
%         d(n-1) = v(n-1)*t_step;
%         epsilon(n) = (v(n-1)^2)/2-(mu/r(n-1)) - a_D(n-1)*d(n-1);
%         r(n) = mu/(v(n-1)^2);
%         v(n) = sqrt(2*(epsilon(n)+mu/r(n)));
        
%         %change in velocity due to drag method
%         rho(n) = rho_0*exp(-(alt)/h_s) %atmospheric density at initial altitude (kg/m^3)
%         %v(t-1) = sqrt(mu/(R+alt(t-1))); %initial orbit velocity (km/s)
%         a_D(n-1) = (0.5*rho(n-1)*v(n-1)^2)/beta
%         v(n) = v(n-1) - a_D(n-1)*t_step
%         r(n) = mu/(v(n)^2)
%         alt = r(n)-R
%         t = t+t_step %advance simulation by one day
%         n = n+1 %increment dataset index
    %end
    
    
end
t_step = 30*86400; %time step (s)
beta = beta(1,1)
n = 1 %dataset index
alt_0 = 300 %initial orbit altitude (km)
r(n) = alt_0+R
P(n)=2*pi*sqrt(r(n)^3/mu) %orbit period (s)
rho(n) = rho_0*exp(-(r(n)-R)/h_s)
dP_dt(n) = (-3*pi*r(n)*rho(n))/beta

n = n+1
P(n) = P(n-1)+dP_dt(n-1)*t_step
r(n) = ((P(n)^2*mu)/(4*pi^2))^(1/3)
rho(n) = rho_0*exp(-(r(n)-R)/h_s)
dP_dt(n) = (-3*pi*(r(n)+R)*rho(n))/beta

while r > 180+R
    n = n+1;
    P(n) = P(n-1)+dP_dt(n-1)*t_step
    r(n) = ((P(n)^2*mu)/(4*pi^2))^(1/3)
    rho(n) = rho_0*exp(-(r(n)-R)/h_s)
    dP_dt(n) = (-3*pi*(r(n)+R)*rho(n))/beta
end    


%plot(t,r)
%atmospheric density as function of altitude
%rho=rho_0*exp(-(r-r0)/h_s); %kg/m^3

%a_drag = (0.5*rho*v^2)/beta