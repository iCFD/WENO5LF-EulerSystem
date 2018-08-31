function res = FD_charWise_WENO5LF1d(a,q,dx)
% *************************************************************************
%
%    Characteristic-wise Finite Difference-1d for the Euler Equations 
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
% NOTE: the reconstruction is performed using characteristic decomposition
% NOTE: Roe averages are assumed for the properties at the cell interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global gamma; gamma1=gamma-1; N=size(q,2); R=3; I=R:(N-R); % R: stencil size 

%% Compute eigenvectors at the cell interfaces 
evr = zeros(3,3,N);
evl = zeros(3,3,N);
hn = zeros(size(q)); 
hp = zeros(size(q));

for i=I, is=i+(-2:3);
    % Compute properties at cell interfaces using Roe avegares
    r_sqrtl = sqrt(q(1,i-1));
    r_sqrtr = sqrt(q(1, i ));
    pl = gamma1*(q(3,i-1) - 0.5*(q(2,i-1)^2)/q(1,i-1));
    pr = gamma1*(q(3, i ) - 0.5*(q(2, i )^2)/q(1, i ));
    r_sq2 = r_sqrtl + r_sqrtr;
    u = (q(2,i-1)/r_sqrtl + q(2,i)/r_sqrtr)/r_sq2;
    H = (((q(3,i-1)+pl)/r_sqrtl + (q(3,i)+pr)/r_sqrtr))/r_sq2;
    c2 = gamma1*(H - 0.5*u^2);
    c = sqrt(c2);

    % Construct matrix of right eigenvectors
    %      _                    _ 
    %     |                      |
    %     |   1      1       1   |
    %     |                      |
    % R = |  u-c     u      u+c  |
    %     |                      |
    %     |  H-uc   u^2/2   H+uc |
    %     |_                    _|
    
    evr(:,:,i) = [...
          1  ,  1  ,  1   ;...
         u-c ,  u  , u+c  ;...
        H-u*c,u^2/2,H+u*c];

    % Construct matrix of left eigenvectors
    %                          _                                       _ 
    %                         |                                         |
    %                         |  uc/(gamma-1)+u^2/2  -c/(gamma-1)-u   1 |
    %                         |                                         |
    % R^{-1}=(gamma-1)/(2c^2)*|  2(H-u^2)             2u             -2 |
    %                         |                                         |
    %                         | -uc/(gamma-1)+u^2/2   c/(gamma-1)-u   1 |
    %                         |_                                       _|

    evl(:,:,i) = gamma1/(2*c^2)*[...
         c*u/gamma1+u^2/2,-(c/gamma1+u), 1 ;...
              2*(H-u^2)  ,    2*u      ,-2 ;...
        -c*u/gamma1+u^2/2, c/gamma1-u  , 1];

    % Project q and F to local characteristic fields
	w=evl(:,:,i)*q(:,is); g=evl(:,:,i)*F(q(:,is));

    % Reconstruct Right Flux:
    % Using the positive fluxes from the LF-splitting for $u_{i+1/2}^{+}$
    vmm = 0.5*(g(:,1)+a*w(:,1));
    vm  = 0.5*(g(:,2)+a*w(:,2));
    v   = 0.5*(g(:,3)+a*w(:,3));
    vp  = 0.5*(g(:,4)+a*w(:,4));
    vpp = 0.5*(g(:,5)+a*w(:,5));

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

    % Numerical Flux at cell boundary, $u_{i+1/2}^{+}$;
    hn(:,i) = evr(:,:,i)*(w0n.*p0n + w1n.*p1n + w2n.*p2n);

    % Reconstruct Left Flux 
    % Using the negative fluxes from the LF-splitting for $u_{i+1/2}^{-}$
    umm = 0.5*(g(:,2)-a*w(:,2));
    um  = 0.5*(g(:,3)-a*w(:,3));
    u   = 0.5*(g(:,4)-a*w(:,4));
    up  = 0.5*(g(:,5)-a*w(:,5));
    upp = 0.5*(g(:,6)-a*w(:,6));

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

    % Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
    hp(:,i) = evr(:,:,i)*(w0p.*p0p + w1p.*p1p + w2p.*p2p);
end

%% Compute finite volume residual term, df/dx.
res = (hp-circshift(hp,[0,1])+hn-circshift(hn,[0,1]))/dx;

end

% Compute flux vector
function flux = F(q)
    global gamma
    
    % primary properties
    rho=q(1,:); u=q(2,:)./rho; E=q(3,:); p=(gamma-1)*(E-0.5*rho.*u.^2);
    
    % flux vector of conserved properties
    flux=[rho.*u; rho.*u.^2+p; u.*(E+p)];
end