%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ina Vollmer
% 2018
% This script implements the singular value decomposition to find the true
% components of XANES spectra collected at different X-ray emission energies
% to find site-selective structural differences.
% Formulas are implemented according to:  
% G. Smolentsev, G. Guilera, M. Tromp, S. Pascarelli, A.V. Soldatov, 
% Local structure of reaction intermediates probed by time-resolved x-ray 
% absorption near edge structure spectroscopy, 
% The Journal of Chemical Physics, 130 (2009) 174508.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%an .xlsx file has to be in the same folder as this script containing the
%spectra
files = dir('*.xlsx');
filename = files(3,1).name;
Y1 = xlsread(filename,1,'C12:J296');
X = xlsread(filename,1,'B12:B296');

%singular value decomposition using built in function of matlab according to
%formula (1) rows (i) correspond to the energy while colums (j) correspond
%to the different energies at which the spectra were taken in our case
%k refers to number of true components to describe the spectra
%S, s_ik := absorption coefficient for the component number k
%L, l_km := elements of the diagonal matrix with the diagonal elements 
%           sorted in decreased order
Y = Y1';
Y(4,:)=[];
[S,L,V] = svd(Y);
%visualize the components of S
% figure(1)
% plot(X,S(:,1),X,S(:,2),X,S(:,3),X,S(:,4),X,S(:,5),X,S(:,6),X,S(:,7),X,S(:,8))
% make S smaller
% first n components
np = 2;
e = size(Y);
ns = e(1);
S1 = V(1:np,:);
x0 = -ones(ns,np);
fun = @(x)Y - x*S1;
options = optimoptions('fsolve','Display','off');
[x,fval,exitflag,output] = fsolve(fun,x0,options);
% constrain the concentrations so that they add up to 1 according to:
% P. Glatzel, L. Jacquamet, U. Bergmann, F.M.F. de Groot, S.P. Cramer, 
% Site-Selective EXAFS in Mixed-Valence Compounds Using High-Resolution 
% Fluorescence Detection:? A Study of Iron in Prussian Blue, Inorganic 
% Chemistry, 41 (2002) 3121-3127.
% number of unknowns
u = np^2+ns*np;
c0 = 10*ones(u,1);
A = [];
b = [];
Aeq = [];
beq = [];
lb = zeros(u,1);
ub = ones(u,1);
lb(1:np^2,1) = -Inf;
ub(1:np^2,1) = Inf;
nonlcon = [];
c = lsqnonlin(@(c)fun1(c,x),c0,lb,ub)
figure(2)
m = [c(1),c(3);c(2),c(4)];
for i = 1:(2*ns)
    if i<=ns
    cf(i,1) = c(np^2+i);
    else i>=ns
    cf(i-ns,2) = c(np^2+i);
    end
end
Yf = cf*m^-1*S1;
% % spectra of the components
F = m^-1*S1;
hold on
p = plot(X',F(2,:),'--',X',F(1,:),'--','LineWidth',1);
p(1).Color = 'red';
p(2).Color = 'green';
%figure(3)
plot(X',Y(5,:),'red',X',Y(7,:),'green')
%plot(X',Yf(5,:),'red',':',X',Yf(7,:),':','green','LineWidth',1)
xlabel('\bf \it energy \rm / eV','FontSize',16,'FontName','Calibri')
ylabel('\bf \it normalized absorption \rm / \mu (E)','FontSize',16,'FontName','Calibri')
lgd = legend('component2','component1','carb','He pretreat');
%,'carb fit','He pretreat fit'
lgd.FontName = 'Calibri';
lgd.FontSize =16;
legend('Location','southeast')
legend('boxoff')

axis([7100 7165 0 Inf])
ax = gca;
c = ax.LineWidth;
ax.LineWidth = 2;
ax.TickDir = 'out';
ax.FontSize = 14;
ax.FontName = 'Calibri';
    function F = fun1(c,x)
    F(1) = c(5) -(x(1,1)*c(1)+x(1,2)*c(2));
    F(2) = 1-c(5) -(x(1,1)*c(3)+x(1,2)*c(4)); 
    F(3) = c(6) -(x(2,1)*c(1)+x(2,2)*c(2));
    F(4) = 1-c(6) -(x(2,1)*c(3)+x(2,2)*c(4));  
    F(5) = c(7) -(x(3,1)*c(1)+x(3,2)*c(2));
    F(6) = 1-c(7) -(x(3,1)*c(3)+x(3,2)*c(4));
    F(7) = c(8) -(x(4,1)*c(1)+x(4,2)*c(2));
    F(8) = 1-c(8) -(x(4,1)*c(3)+x(4,2)*c(4));
    F(9) = c(9) -(x(5,1)*c(1)+x(5,2)*c(2));
    F(10) = 1-c(9) -(x(5,1)*c(3)+x(5,2)*c(4));
    F(11) = c(10) -(x(6,1)*c(1)+x(6,2)*c(2));
    F(12) = 1-c(10) -(x(6,1)*c(3)+x(6,2)*c(4));
    F(13) = c(11) -(x(7,1)*c(1)+x(7,2)*c(2));
    F(14) = 1-c(11) -(x(7,1)*c(3)+x(7,2)*c(4));
    F(15) = 1- c(12)-c(5);
    F(16) = 1- c(13)-c(6);
    F(17) = 1- c(14)-c(7);
    F(18) = 1- c(15)-c(8);
    F(19) = 1- c(16)-c(9);
    F(20) = 1- c(17)-c(10);
    F(21) = 1- (c(18)+c(11));
    F(22) = 1-c(11);
    end



