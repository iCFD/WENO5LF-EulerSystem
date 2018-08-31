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

% Set plotting defaults
set(0,'defaultTextInterpreter','latex')
set(0,'DefaultTextFontName','Times',...
'DefaultTextFontSize',14,...
'DefaultAxesFontName','Times',...
'DefaultAxesFontSize',14,...
'DefaultLineLineWidth',1.5,...
'DefaultAxesBox','on',...
'defaultAxesLineWidth',1.5,...
'DefaultFigureColor','w',...
'DefaultLineMarkerSize',7.75)

% Compare
load('Component-wise FD.mat'); r_FD_1=r;
load('Component-wise FV.mat'); r_FV_1=r;
load('Primitive-wise FV.mat'); r_FV_2=r;
load('Characteristic-wise FD.mat'); r_FD_2=r;
load('Characteristic-wise FV.mat'); r_FV_3=r;

% Plot Full Figure
figure(1)
plot(xe,re,'-k'); hold on 
plot(x,r_FD_1,'-');
plot(x,r_FD_2,'-');
plot(x,r_FV_1,'-');
plot(x,r_FV_2,'-');
plot(x,r_FV_3,'-'); hold off
title('The Gas Dynamics Euler Equations (EE)'); 
xlabel('$x$','interpreter','latex'); ylabel('$\varrho(x)$','interpreter','latex');
legend({'Exact','FD Component-wise','FD Characteristic-wise','FV Component-wise','FV Primitive-wise','FV Characteristic-wise'},...
    'interpreter','latex','location','southwest'); 
legend boxoff
fig=gcf; fig.PaperUnits='inches'; fig.PaperPosition=[0 0 6 4];
print('-depsc','SodShockTube1d.eps');

% Plot subfigures 
figure(2)
subplot(131)
plot(xe,re,'-k'); hold on 
plot(x,r_FD_1,'-x');
plot(x,r_FD_2,'-o');
plot(x,r_FV_1,'-+');
plot(x,r_FV_2,'-s');
plot(x,r_FV_3,'-^'); hold off
title('Expansion Fan'); 
xlabel('$x$','interpreter','latex'); ylabel('$\varrho(x)$','interpreter','latex');
legend({'Exact','FD-comp','FD-char','FV-comp','FV-prim','FV-char'},...
    'interpreter','latex'); 
legend boxoff
xlim([0.5,0.65]);

subplot(132)
plot(xe,re,'-k'); hold on 
plot(x,r_FD_1,'-x');
plot(x,r_FD_2,'-o');
plot(x,r_FV_1,'-+');
plot(x,r_FV_2,'-s');
plot(x,r_FV_3,'-^'); hold off
title('Contact'); 
xlabel('$x$','interpreter','latex'); ylabel('$\varrho(x)$','interpreter','latex');
legend({'Exact','FD-comp','FD-char','FV-comp','FV-prim','FV-char'},...
    'interpreter','latex'); 
legend boxoff
xlim([0.7,0.8]);

subplot(133)
plot(xe,re,'-k'); hold on 
plot(x,r_FD_1,'-x');
plot(x,r_FD_2,'-o');
plot(x,r_FV_1,'-+');
plot(x,r_FV_2,'-s');
plot(x,r_FV_3,'-^'); hold off
title('Shock'); 
xlabel('$x$','interpreter','latex'); ylabel('$\varrho(x)$','interpreter','latex');
legend({'Exact','FD-comp','FD-char','FV-comp','FV-prim','FV-char'},...
    'interpreter','latex'); 
legend boxoff
xlim([0.85,0.9]);
fig=gcf; fig.PaperUnits='inches'; fig.PaperPosition=[0 0 12 6];
print('-depsc','SodShockTube1d_basicStructures.eps');