%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Solving 1-D Euler system of equations with 5th order
%          Weighted Essentially Non-Oscilaroty (MOL-WENO5-LF)
%
%        dq_i/dt + df_i/dx = 0, for x \in [a,b] and i =1,. ..,D
%
%           coded by Manuel A. Diaz, manuel.ade'at'gmail.com 
%            Institute of Applied Mechanics, NTU, 2012.08.25
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code solves the Sod's shock tube problem (IC=1)
%
% t=0                                 t=tEnd
% Density                             Density
%   ****************|                 *********\
%                   |                           \
%                   |                            \
%                   |                             ****|
%                   |                                 |
%                   |                                 ****|
%                   ***************                       ***********
%
% coded by Manuel A. Diaz, 2012.12.27. Last modif: 20.06.2018.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs:
% [1] Shu, Chi-Wang. "Essentially non-oscillatory and weighted essentially 
%     non-oscillatory schemes for hyperbolic conservation laws." Advanced 
%     numerical approximation of nonlinear hyperbolic equations. Springer, 
%     Berlin, Heidelberg, 1998. 325-432.
% [2] Jiang, Guang-Shan, and Cheng-chin Wu. "A high-order WENO finite
%     difference scheme for the equations of ideal magnetohydrodynamics."
%     Journal of Computational Physics 150.2 (1999): 561-594.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: 
% 1. A fully conservative finite volume implementation of the method of
% lines (MOL) using WENO5 associated with SSP-RK33 time integration method. 
% 2. Sharpenning of contact discontinuities is NOT implemented here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; %close all; clc;
global gamma

%% Parameters
tFinal  = 0.2;  % Desired output time 
CFL     = 0.55;	% CFL number
nE      = 200;  % Number of cells/Elements
gamma   = 1.4;  % Ratio of specific heats for ideal di-atomic gas
IC      = 07;	% 10 IC cases are available
plotFig = true; % Plot figures every 10 time steps
scheme  = 6;    % 1: Component-wise FD, flux splitting: WENO5 from Ref.[1]
                % 2: Component-wise FD, flux splitting: WENO5 from Ref.[2]
                % 3: Characteristic-wise FD, flux splitting: WENO5
                % 4: Component-wise FV: WENO5
                % 5: Primitive-wise FV: WENO5
                % 6: Characteristic-wise FV: WENO5
                
% Discretize spatial domain
a=0; b=1; dx=(b-a)/nE; nx=nE+1; x=linspace(a,b,nx); bl=1:3; br=nx-2:nx;

% Set IC
[r0,u0,p0,tFinalIC,~]=Euler_IC1d(x,IC); tFinal=min(tFinal,tFinalIC);
E0 = p0./((gamma-1))+0.5*r0.*u0.^2;  % Total Energy density
a0 = sqrt(gamma*p0./r0);	% Speed of sound
q0=[r0; r0.*u0; E0];        % vec. of conserved properties

% Exact solution (needs to be improved!)
[xe,re,ue,pe,ee,te,Me,se] = ...
    EulerExact(r0(1),u0(1),p0(1),r0(nx),u0(nx),p0(nx),tFinal);

% Discretize time domain
lambda0=max(abs(u0)+a0); dt0=CFL*dx/lambda0;  % using the system's largest eigenvalue

% Set reconstruction scheme
switch scheme
    case 1, WENO5LF1d = @FD_compWise_WENO5LF1d; method_name = 'Component-wise FD';
    case 2, WENO5LF1d = @FD_compWise_WENO5LFv2; method_name = 'Component-wise FD';
    case 3, WENO5LF1d = @FD_charWise_WENO5LF1d; method_name = 'Characteristic-wise FD';
    case 4, WENO5LF1d = @FV_compWise_WENO5LF1d; method_name = 'Component-wise FV';
    case 5, WENO5LF1d = @FV_primWise_WENO5LF1d; method_name = 'Primitive-wise FV';
    case 6, WENO5LF1d = @FV_charWise_WENO5LF1d; method_name = 'Characteristic-wise FV';
end

% Prepare figure
figure(1); set(gcf,'position',[100,100,800,550]);

%% Solver Loop
% Load initial condition
q=q0; it=0; dt=dt0; t=0; lambda=lambda0;

tic
while t<tFinal
    % Iteration current time
    if t+dt>tFinal; dt=tFinal-t; end; t=t+dt;
    
    % RK Initial step
    qo = q;
    
    % 1st stage
    L=WENO5LF1d(lambda,q,dx);     q = qo-dt*L;
    q(:,bl)=qo(:,bl); q(:,br)=qo(:,br); % Neumann BCs
    
    % 2nd Stage
    L=WENO5LF1d(lambda,q,dx);     q = 0.75*qo+0.25*(q-dt*L);
    q(:,bl)=qo(:,bl); q(:,br)=qo(:,br); % Neumann BCs
    
    % 3rd stage
    L=WENO5LF1d(lambda,q,dx);     q = (qo+2*(q-dt*L))/3;
    q(:,bl)=qo(:,bl); q(:,br)=qo(:,br); % Neumann BCs
    
    % compute primary properties
    r=q(1,:); u=q(2,:)./r; E=q(3,:); p=(gamma-1)*(E-0.5*r.*u.^2);
    a=sqrt(gamma*p./r); if min(p)<0; error('negative pressure found!'); end
    
    % Compute dt for next time step
    lambda=max(abs(u)+a); dt=CFL*dx/lambda; 
    
    % Update iteration counter
    it=it+1;
    
    % Plot figure
    if plotFig && rem(it,10) == 0
        subplot(2,2,1); plot(x,r,'.b');
        subplot(2,2,2); plot(x,u,'.m');
        subplot(2,2,3); plot(x,p,'.k');
        subplot(2,2,4); plot(x,E,'.r');
        drawnow
    end
end
cputime=toc; disp(cputime);

%% Post Process
% Solution fields
r=q(1,:); u=q(2,:)./r; E=q(3,:); p=(gamma-1)*(E-0.5*r.*u.^2);

% Calculation of flow parameters
a = sqrt(gamma*p./r); M = u./a; % Mach number [-]
p_ref = 101325;         % Reference air pressure (N/m^2)
r_ref= 1.225;           % Reference air density (kg/m^3)
s_ref = 1/(gamma-1)*(log(p/p_ref)+gamma*log(r_ref./r)); 
                        % Entropy w.r.t reference condition
s = log(p./r.^gamma);	% Dimensionless Entropy
Q = r.*u;               % Mass Flow rate per unit area
e = p./((gamma-1)*r);	% internal Energy

% Final plot
s1=subplot(2,3,1); plot(x,r,'or',xe,re,'k'); xlabel('x(m)'); ylabel('Density (kg/m^3)');
s2=subplot(2,3,2); plot(x,u,'or',xe,ue,'k'); xlabel('x(m)'); ylabel('Velocity (m/s)');
s3=subplot(2,3,3); plot(x,p,'or',xe,pe,'k'); xlabel('x(m)'); ylabel('Pressure (Pa)');
s4=subplot(2,3,4); plot(x,s,'or',xe,se,'k'); xlabel('x(m)'); ylabel('Entropy/R gas');
s5=subplot(2,3,5); plot(x,M,'or',xe,Me,'k'); xlabel('x(m)'); ylabel('Mach number');
s6=subplot(2,3,6); plot(x,e,'or',xe,ee,'k'); xlabel('x(m)'); ylabel('Internal Energy (kg/m^2s)');
title(s2,[method_name,' WENO-LF Euler Solver']);

% save 
save([method_name,'.mat'],'x','r','xe','re');