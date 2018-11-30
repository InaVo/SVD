%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Ina Vollmer
% 2018
% This script implements the singular value decomposition to find the true
% components of XANES spectra collected at different X-ray emission energies
% to find site-selective structural differences.
% Most formulas are implemented according to:  
% P. Glatzel, L. Jacquamet, U. Bergmann, F.M.F. de Groot, S.P. Cramer, 
% Site-Selective EXAFS in Mixed-Valence Compounds Using High-Resolution 
% Fluorescence Detection:? A Study of Iron in Prussian Blue, 
% Inorganic Chemistry, 41 (2002) 3121-3127.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%an .xlsx file has to be in the same folder as this script containing the
%spectra
files = dir('*.xlsx');
filename = files(2,1).name;
%read in spectra
Y1 = xlsread(filename,1,'C12:J296');
Y = Y1';%transpose of read in matrix to make compatible with notation from paper.
% read in vector of energies
E1 = xlsread(filename,1,'B12:B296');
E = E1';%transpose of read in vector to make compatible with notation from paper.
% singular value decomposition using built in function of matlab
% rows (i) correspond to the spectra colums (j) correspond to the energies 
[V,L,S] = svd(Y);
% total number of spectra
e = size(Y);
ns = e(1);
% select only the spectra of the first np components
np = 2; % first np components
Sp = S(1:np,:);
%create initial guess for solver
x0 = -ones(ns,np);
fun = @(x)Y - x*Sp;
options = optimoptions('fsolve','Display','off');
[x,fval,exitflag,output] = fsolve(fun,x0,options);
% Find real concentrations
c0 = 10*ones(20,1); % create initial guesses
%create constraints for solver
A = [];
b = [];
Aeq = [];
beq = [];
%lower and upper bounds for the concentrations
lb = zeros(20,1);
ub = ones(20,1);
lb(1:4,1) = -Inf;
ub(1:4,1) = Inf;
nonlcon = [];
% calling solver to find real concentrations
c_real = lsqnonlin(@(c)fun1(c,x),c0,lb,ub);%function fun1 is defined at the end of script
% The first np^2 components of c_real contain the entries of T
T = [c_real(1),c_real(3);c_real(2),c_real(4)];
% creating matrix with real concentrations
C_real = [c_real(5),c_real(13);c_real(6),c_real(14);c_real(7),c_real(15);c_real(8),c_real(16);c_real(9),c_real(17);c_real(10),c_real(18);c_real(11),c_real(19);c_real(12),c_real(20)];
% calcuating fitted spectra
Yf = C_real*T^-1*Sp;
% % spectra of the components
S_real = T^-1*Sp;

%plotting the fit together with original spectra
figure(2)
hold on
p = plot(E,Yf(6,:),'--',E,Yf(8,:),'--','LineWidth',1); %fit
p(1).Color = 'red';
p(2).Color = 'green';

plot(E,Y(6,:),'red',E,Y(8,:),'green') % original spectra
%style settings:
xlabel('\bf \it energy \rm / eV','FontSize',16,'FontName','Calibri')
ylabel('\bf \it normalized absorption \rm / \mu (E)','FontSize',16,'FontName','Calibri')
lgd = legend('fit CH_4 flow 6.058 eV','fit He pretreat','carb','He pretreat');
lgd.FontName = 'Calibri';
lgd.FontSize =16;
legend('Location','southeast')
legend('boxoff')

axis([7100 7165 0 Inf])
ax = gca;
c_real = ax.LineWidth;
ax.LineWidth = 2;
ax.TickDir = 'out';
ax.FontSize = 14;
ax.FontName = 'Calibri';

%plotting the real components
figure(3)
hold on
p = plot(E,S_real(2,:),'--',E,S_real(1,:),'--','LineWidth',1);
p(1).Color = 'red';
p(2).Color = 'green';
xlabel('\bf \it energy \rm / eV','FontSize',16,'FontName','Calibri')
ylabel('\bf \it normalized absorption \rm / \mu (E)','FontSize',16,'FontName','Calibri')
lgd = legend('component2','component1');
lgd.FontName = 'Calibri';
lgd.FontSize =16;
legend('Location','southeast')
legend('boxoff')

axis([7100 7165 0 Inf])
ax = gca;
c_real = ax.LineWidth;
ax.LineWidth = 2;
ax.TickDir = 'out';
ax.FontSize = 14;
ax.FontName = 'Calibri';

% creating the functions for lsqnonlin solver
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
    F(15) = c(12) -(x(8,1)*c(1)+x(8,2)*c(2));
    F(16) = 1-c(12) -(x(8,1)*c(3)+x(8,2)*c(4));
    F(17) = 1- c(13)-c(5);
    F(18) = 1- c(14)-c(6);
    F(19) = 1- c(15)-c(7);
    F(20) = 1- c(16)-c(8);
    F(21) = 1- c(17)-c(9);
    F(22) = 1- c(18)-c(10);
    F(23) = 1- (c(19)+c(11));
    F(24) = 1- (c(20)+c(12));
    end



