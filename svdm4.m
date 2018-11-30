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
%X = H(:,1);
%Y = H(:,2:3);
%singular value decomposition using built in function of matlab according to
%formula (1) rows (i) correspond to the energy while colums (j) correspond
%to the different energies at which the spectra were taken in our case
%k refers to number of true components to describe the spectra
%S, s_ik := absorption coefficient for the component number k
%L, l_km := elements of the diagonal matrix with the diagonal elements 
%           sorted in decreased order
Y = Y1';
%Y(4,:)=[];
[S,L,V] = svd(Y);
%visualize the components of S
% figure(1)
% plot(X,S(:,1),X,S(:,2),X,S(:,3),X,S(:,4),X,S(:,5),X,S(:,6),X,S(:,7),X,S(:,8))
% make S smaller
% first n components
np = 3;
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
c0 = ones(u,1);
A = [];
b = [];
Aeq = [];
beq = [];
lb = zeros(u,1);
ub = ones(u,1);
lb(1:np^2,1) = -Inf;
ub(1:np^2,1) = Inf;
nonlcon = [];
c = lsqnonlin(@(c)fun1(c,x,np,ns),c0,lb,ub)

for i = 1:np^2
    if i <=np
        m(i,1)=c(i);
    elseif i>np && i <=2*np
        m(i-np,2) = c(i);
    elseif i>2*np
        m(i-2*np,3) = c(i);
    end
end       
        

for i = 1:(np*ns)
    if i<=ns
    cf(i,1) = c(np^2+i);
    elseif i>ns && i<=2*ns
    cf(i-ns,2) = c(np^2+i);
    else i>2*ns
    cf(i-2*ns,3) = c(np^2+i);
    end
end
Yf = cf*m^-1*S1;
% % spectra of the components
F = m^-1*S1;

figure(1)
hold on

p = plot(X',F(2,:),'--',X',F(1,:),'--',X',F(3,:),'--','LineWidth',1);
p(1).Color = 'red';
p(2).Color = 'green';
p(3).Color = 'blue';

plot(X',Y(6,:),'red',X',Y(8,:),'green')
%plot(X',Yf(5,:),'red',':',X',Yf(7,:),':','green','LineWidth',1)
xlabel('\bf \it energy \rm / eV','FontSize',16,'FontName','Calibri')
ylabel('\bf \it normalized absorption \rm / \mu (E)','FontSize',16,'FontName','Calibri')
lgd = legend('component2','component1','component3','carb','He pretreat');
%,'carb fit','He pretreat fit'
lgd.FontName = 'Calibri';
lgd.FontSize =16;
legend('Location','southeast')
legend('boxoff')

axis([7100 7165 0 1.4])
ax = gca;
c = ax.LineWidth;
ax.LineWidth = 2;
ax.TickDir = 'out';
ax.FontSize = 14;
ax.FontName = 'Calibri';
figure(2)
hold on
plot(X',Yf(6,:),':',X',Yf(8,:),':','LineWidth',1.5)
plot(X',Y(6,:),'red',X',Y(8,:),'green')
xlabel('\bf \it energy \rm / eV','FontSize',16,'FontName','Calibri')
ylabel('\bf \it normalized absorption \rm / \mu (E)','FontSize',16,'FontName','Calibri')
    function F = fun1(c,x,np,ns)
    y = np^2+1;
    l = 1;
    for i = 1:np*ns
        if mod(i,3) == 2
            F(i) = 1-c(2*ns+y)-c(y) -(x(l,1)*c(4)+x(l,2)*c(5)+x(l,3)*c(6));
        elseif mod(i,3) == 0
            F(i) = 1-c(y+ns)-c(y) -(x(l,1)*c(7)+x(l,2)*c(8)+x(l,3)*c(9));
            y = y +1;
            l = l+1;
        else
            F(i) = c(y) -(x(l,1)*c(1)+x(l,2)*c(2)+x(l,3)*c(3));
        end
    end
    for i = 1:ns
    F(np*ns+i) = 1- c(ns+np^2+i)-c(np^2+i)-c(2*ns+np^2+i);
    end
    F(np*ns+ns+1) = c(np^2+ns)-1;
    F(np*ns+ns+2) = c(np^2+ns+ns);
%     F(np*ns+ns+3) = c(np^2+ns*np);
%     F
    end



