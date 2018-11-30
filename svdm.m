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
filename = files(1,1).name;
Y = xlsread(filename,1,'C12:J296');
X = xlsread(filename,1,'B12:B296');
%X = H(:,1);
%Y = H(:,2:3);
%singular value decomposition using built in function of matlab according to
%formula (1) rows (i) correspond to the energy while colums (j) correspond
%to the different energies at which the spectra were taken in our case
%k refers to number of true components to describe the spectra
%S, s_ik := absorption coefficient for the component number k
%L, l_km := elements of the diagonal matrix with the diagonal elements 
%           sorted in decreased order
[S,L,V] = svd(Y);
%W = L*V:=concentration dependence for the component k
%w_kj concentration dependence of spectra j on component k
W = L*V;
%visualize the components of S
% figure(1)
% plot(X,S(:,1),X,S(:,2),X,S(:,3),X,S(:,4),X,S(:,5),X,S(:,6),X,S(:,7),X,S(:,8))
% make S smaller
% first n components
n = 2;
S1 = S(:,1:n);
x0 = -ones(n,8);
fun = @(x)Y - S1*x;
options = optimoptions('fsolve','Display','off');
[x,fval,exitflag,output] = fsolve(fun,x0,options);
% constrain the concentrations so that they add up to 1 according to:
% P. Glatzel, L. Jacquamet, U. Bergmann, F.M.F. de Groot, S.P. Cramer, 
% Site-Selective EXAFS in Mixed-Valence Compounds Using High-Resolution 
% Fluorescence Detection:? A Study of Iron in Prussian Blue, Inorganic 
% Chemistry, 41 (2002) 3121-3127.
d0 = -ones(1,n);
un = ones(1,8);
fun1 = @(d)un - d*x;
options = optimoptions('fsolve','Display','off');
[d,fval1,exitflag1,output1] = fsolve(fun1,d0,options);
figure(2)
plot(X,Y(:,6))
% for i = 1:10
%     m4 = 0+0.1*i;
%     m2 = -m4*c(8,2)/c(8,1);
%     m1 = d(1)+ m4*c(8,2)/c(8,1);
%     m3 = d(2) - m4;
% m = [m1,m2;m3,m4];
% cf = c*m
% end
b6 = 0.1;
    m4 = -0.4;
    m2 = b6/x(1,6)-m4*x(2,6)/x(1,6);
    m1 = d(1)+ m2;
    m3 = d(2) - m4;
m = [m1,m2;m3,m4];
cf = m^-1*x;
Yf = S1*m*cf;
% spectra of the components
F = S1*m^-1;
hold on
plot(X,Yf(:,6),X,F(:,1))
plot(X,F(:,2))
%plot(X,Y(:,8))
%     m4 = 0.883/c(6,1)*1/(1-c(6,1)*c(8,2)/c(6,1)*c(8,1));
%     m2 = -m4*c(8,2)/c(8,1);
%     m1 = d(1)+ m4*c(8,2)/c(8,1);
%     m3 = d(2) - m4;
% m = [m1,m2;m3,m4];
% cf = c*m
