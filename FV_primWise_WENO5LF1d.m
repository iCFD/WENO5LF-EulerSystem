function res = FV_primWise_WENO5LF1d(a,q,dx)
% *************************************************************************
%
%    Component-wise Finite Volume-1d for the Euler Equations 
%
% Based on:
% ---------
%   Shu, Chi-Wang. "Essentially non-oscillatory and weighted essentially 
%   non-oscillatory schemes for hyperbolic conservation laws." Advanced 
%   numerical approximation of nonlinear hyperbolic equations. Springer, 
%   Berlin, Heidelberg, 1998. 325-432. 
%
% coded by Manuel Diaz, 02.10.2012, NTU Taiwan.
%       last updated on 2018.06.20, NHRI Taiwan.
% *************************************************************************
%
% Domain cells (I{i}) reference:
%
%                |           |   u(i)    |           |
%                |  u(i-1)   |___________|           |
%                |___________|           |   u(i+1)  |
%                |           |           |___________|
%             ...|-----0-----|-----0-----|-----0-----|...
%                |    i-1    |     i     |    i+1    |
%                |-         +|-         +|-         +|
%              i-3/2       i-1/2       i+1/2       i+3/2
%
% ENO stencils (S{r}) reference:
%
%
%                               |___________S2__________|
%                               |                       |
%                       |___________S1__________|       |
%                       |                       |       |
%               |___________S0__________|       |       |
%             ..|---o---|---o---|---o---|---o---|---o---|...
%               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
%                                      -|
%                                     i+1/2
%
%
%               |___________S0__________|
%               |                       |
%               |       |___________S1__________|
%               |       |                       |
%               |       |       |___________S2__________|
%             ..|---o---|---o---|---o---|---o---|---o---|...
%               | I{i-2}| I{i-1}|  I{i} | I{i+1}| I{i+2}|
%                               |+
%                             i-1/2
%
% WENO stencil: S{i} = [ I{i-2},...,I{i+2} ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Careful!: by using circshift over our domain, we are implicitly creating a
% favorable code that automatically includes periodical boundary conditions.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global gamma; nx=size(q,2); R=3; I=R:(nx-R); % R: stencil size 

% Compute primitive variables at solution points
w(1,:) = q(1,:);
w(2,:) = q(2,:)./q(1,:);
w(3,:) = (gamma-1)*( q(3,:) - 0.5*(q(2,:).^2)./q(1,:));

%% Right Flux
% Choose the positive fluxes, 'v', to compute the left cell boundary flux:
% $u_{i+1/2}^{-}$
vmm = w(:,I-2);
vm  = w(:,I-1);
v   = w(:, I );
vp  = w(:,I+1);
vpp = w(:,I+2);

% Polynomials
p0n = (2*vmm - 7*vm + 11*v)/6;
p1n = ( -vm  + 5*v  + 2*vp)/6;
p2n = (2*v   + 5*vp - vpp )/6;

% Smooth Indicators (Beta factors)
B0n = 13/12*(vmm-2*vm+v  ).^2 + 1/4*(vmm-4*vm+3*v).^2; 
B1n = 13/12*(vm -2*v +vp ).^2 + 1/4*(vm-vp).^2;
B2n = 13/12*(v  -2*vp+vpp).^2 + 1/4*(3*v-4*vp+vpp).^2;

% Constants
d0n = 1/10; d1n = 6/10; d2n = 3/10; epsilon = 1e-6;

% Alpha weights 
alpha0n = d0n./(epsilon + B0n).^2;
alpha1n = d1n./(epsilon + B1n).^2;
alpha2n = d2n./(epsilon + B2n).^2;
alphasumn = alpha0n + alpha1n + alpha2n;

% ENO stencils weigths
w0n = alpha0n./alphasumn;
w1n = alpha1n./alphasumn;
w2n = alpha2n./alphasumn;

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
hn = w0n.*p0n + w1n.*p1n + w2n.*p2n;

%% Left Flux 
% Choose the negative fluxes, 'u', to compute the left cell boundary flux:
% $u_{i-1/2}^{+}$ 
umm = w(:,I-1);
um  = w(:, I );
u   = w(:,I+1);
up  = w(:,I+2);
upp = w(:,I+3);

% Polynomials
p0p = ( -umm + 5*um + 2*u  )/6;
p1p = ( 2*um + 5*u  - up   )/6;
p2p = (11*u  - 7*up + 2*upp)/6;

% Smooth Indicators (Beta factors)
B0p = 13/12*(umm-2*um+u  ).^2 + 1/4*(umm-4*um+3*u).^2; 
B1p = 13/12*(um -2*u +up ).^2 + 1/4*(um-up).^2;
B2p = 13/12*(u  -2*up+upp).^2 + 1/4*(3*u -4*up+upp).^2;

% Constants
d0p = 3/10; d1p = 6/10; d2p = 1/10; epsilon = 1e-6;

% Alpha weights 
alpha0p = d0p./(epsilon + B0p).^2;
alpha1p = d1p./(epsilon + B1p).^2;
alpha2p = d2p./(epsilon + B2p).^2;
alphasump = alpha0p + alpha1p + alpha2p;

% ENO stencils weigths
w0p = alpha0p./alphasump;
w1p = alpha1p./alphasump;
w2p = alpha2p./alphasump;

% Numerical Flux at cell boundary, $u_{i-1/2}^{+}$;
hp = w0p.*p0p + w1p.*p1p + w2p.*p2p;

% Compute conservative variables at faces
qn(1,:) = hn(1,:);
qn(2,:) = hn(2,:).*hn(1,:);
qn(3,:) = hn(3,:)./(gamma-1) + 0.5*hn(1,:).*hn(2,:).^2;

qp(1,:) = hp(1,:);
qp(2,:) = hp(2,:).*hp(1,:);
qp(3,:) = hp(3,:)./(gamma-1) + 0.5*hp(1,:).*hp(2,:).^2;

%% Compute finite volume residual term, df/dx.
LF=zeros(size(q)); res=zeros(size(q));

LF(:,I) = 0.5*(F(qn)+F(qp)-abs(a).*(qp-qn)); % Lax friedrichs flux
res(:,I) = (LF(:,I)-LF(:,I-1))/dx; % L = -df(q)/dx.

% Flux contribution of the LEFT MOST FACE: left face of cell j=1.
res(:,3) = res(:,3)-LF(:,3)/dx;
 
% Flux contribution of the RIGHT MOST FACE: right face of cell j=nx-1.
res(:,nx-2)=res(:,nx-2)+LF(:,nx-2)/dx;
end

% Compute flux vector
function flux = F(q)
    global gamma
    
    % primary properties
    rho=q(1,:); u=q(2,:)./rho; E=q(3,:); p=(gamma-1)*(E-0.5*rho.*u.^2);
    
    % flux vector of conserved properties
    flux=[rho.*u; rho.*u.^2+p; u.*(E+p)];
end