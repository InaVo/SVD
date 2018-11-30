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
S1 = S(:,1:2);
x0 = -ones(2,8);
fun = @(x)Y - S1*x;
options = optimoptions('fsolve','Display','off');
[x,fval,exitflag,output] = fsolve(fun,x0,options);
% constrain the concentrations so that they add up to 1 according to:
% P. Glatzel, L. Jacquamet, U. Bergmann, F.M.F. de Groot, S.P. Cramer, 
% Site-Selective EXAFS in Mixed-Valence Compounds Using High-Resolution 
% Fluorescence Detection:? A Study of Iron in Prussian Blue, Inorganic 
% Chemistry, 41 (2002) 3121-3127.
d0 = -ones(2,1);
un = ones(8,1);
fun1 = @(d)un - x'*d;
options = optimoptions('fsolve','Display','off');
[d,fval1,exitflag1,output1] = fsolve(fun1,d0,options);
c = x;
for i = 1:5
    m4 = 0+0.2*i;
    m2 = -m4*c(2,8)/c(1,8);
    m1 = d(1)+ m4*c(2,8)/c(1,8);
end
m = [m1,m2;m3,m4];
cf = c*m;
Yfit = S1*cf;
figure(2)
plot(X,Yfit(:,6),X,Y(:,6))