function res = FD_compWise_WENO5LFv2(a,q,dx)
% *************************************************************************
%
%    Component-wise Finite Volume-1d for the Euler Equations 
%
% Based on:
% ---------
% [1] Shu, Chi-Wang. "Essentially non-oscillatory and weighted essentially 
%     non-oscillatory schemes for hyperbolic conservation laws." Advanced 
%     numerical approximation of nonlinear hyperbolic equations. Springer, 
%     Berlin, Heidelberg, 1998. 325-432.
% [2] Jiang, Guang-Shan, and Cheng-chin Wu. "A high-order WENO finite
%     difference scheme for the equations of ideal magnetohydrodynamics."
%     Journal of Computational Physics 150.2 (1999): 561-594.  
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
% NOTE: 
% This formulation follows exactly the formulation provided by Jiang in [2]
% is a special reformulation of the method suitable for developing the WENO
% reconstructions in characterstic and primitive froms.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The reconstruction is performed component-wise;
% NOTE: equations = components, no characteristic reconstruction is used!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E = size(q,1);
N = size(q,2);
R=3; I=R:(N-R); % R: stencil size 
epweno=1E-40;

qr=zeros(size(q));
ql=zeros(size(q));
LF=zeros(size(q)); 
res=zeros(size(q));

% Loop over all components
for e=1:E

    % Compute and store the differences of the cell averages
    for i=2:N 
        dqmh(i) = q(e,i)-q(e,i-1); % dp_{i-1/2}
    end 

    % the reconstruction
    for idx=1:2

        % idx=1: construct hn
        % idx=2: construct hp

        im=(-1)^(idx+1);
        i1=im; in1=-im; in2=-2*im;

        for i=R:N-R+1

            t1=im*(dqmh(i+in2)-dqmh(i+in1));
            t2=im*(dqmh(i+in1)-dqmh(  i  ));
            t3=im*(dqmh(  i  )-dqmh(i+i1 ));

            IS1=13.*t1^2+3.*(   dqmh(i+in2)-3.*dqmh(i+in1))^2;
            IS2=13.*t2^2+3.*(   dqmh(i+in1)+   dqmh(  i  ))^2;
            IS3=13.*t3^2+3.*(3.*dqmh(  i  )-   dqmh(i+i1 ))^2;

            IS1=(epweno+IS1)^2;
            IS2=(epweno+IS2)^2;
            IS3=(epweno+IS3)^2;
            s1 =IS2*IS3;
            s2 =6.*IS1*IS3;
            s3 =3.*IS1*IS2;
            t0 =1./(s1+s2+s3);
            s1 =s1*t0;
            s3 =s3*t0;

            h(idx,i) = (s1*(t2-t1)+(0.5*s3-0.25)*(t3-t2))/3 ...
                     +(-q(e,i-2)+7*(q(e,i-1)+q(e,i))-q(e,i+1))/12;
        end 
    end 

   qr(e,R-1:N-R  )=h(1,R:N-R+1);
   ql(e,R  :N-R+1)=h(2,R:N-R+1);

end 

%% Compute finite volume residual term, df/dx.
LF(:,I) = 0.5*(F(qr(:,I))+F(ql(:,I+1))-abs(a).*(ql(:,I+1)-qr(:,I))); % Lax friedrichs flux
res(:,I) = (LF(:,I)-LF(:,I-1))/dx; % L = -df(q)/dx.

% Flux contribution of the LEFT MOST FACE: left face of cell j=1.
res(:,3) = res(:,3)-LF(:,3)/dx;
 
% Flux contribution of the RIGHT MOST FACE: right face of cell j=nx-1.
res(:,N-2)=res(:,N-2)+LF(:,N-2)/dx;
end

% Compute flux vector
function flux = F(q)
    global gamma
    
    % primary properties
    rho=q(1,:); u=q(2,:)./rho; E=q(3,:); p=(gamma-1)*(E-0.5*rho.*u.^2);
    
    % flux vector of conserved properties
    flux=[rho.*u; rho.*u.^2+p; u.*(E+p)];
end