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
filename = files(3,1).name;

E1 = xlsread(filename,1,'B30:B296');% reading in vector of energies
E = E1';%transpose to make compatible with notation

Y1 = xlsread(filename,1,'C30:J296');% Reading in spectra
Y2 = Y1';%transpose to make compatible with notation

e = size(Y2);% number of experimental spectra
ns1 = e(1);
%--------------------------------------------------------------------------
% In this section, we have to define, which spectra to use for the fitting
% The code will automatically create the right amount of concentrations and
% functions to solve. Only the exta constraints in function fun1 at the end 
% of the script and the variables nox and ncarb need to be adjusted.
nox = 7;
ncarb = 5;
% If three components should be used, the function fun1 and the variable np 
% as well as the plotting have to be adjusted
np = 2;

% Optional: select a few spectra for SVD randomly starting from hend. The random
% selection seems to improve the noise level
% the first spectra is always the first one
%{
r = 1;
r1 = zeros(ns1);
hend = ns1/2-5;
Y(1,:)=Y2(r,:);
for i = 2:hend
    cond = find(r1==r,1);
    while isempty(cond) == 0 || r <=hend
        r = random('unid',ns1);
        cond = find(r1==r,1);
    end
    if isempty(cond) == 1
    Y(i,:)=Y2(r,:);
    r1(i) = r;
    r = random('unid',ns1);
    end
end
xlswrite('randgen.xlsx',r1(:,1),1,'J1:J28')
Y = Y2(5:7,:);
Y(1,:) = Y2(8,:);
Y(2,:) = Y2(5,:);
Y(3,:) = Y2(7,:);
Y(4,:) = Y2(6,:);
%}
% instead of randomly choosing, we can also define, which spectra to use
% directly:
Y = Y2;
Y(4,:)=[];
e = size(Y);
ns = e(1);  

%singular value decomposition using built in function of matlab
%(i) correspond to the spectra, columns(j) correspond to energies
[V,L,S] = svd(Y);
% first np components chosen according to the highest singular values of L
% the singular values are ordered by magnitude on the diagonal of L
Sp = S(1:np,:);

x0 = ones(ns,np);%creating initial guesses
fun = @(x)Y - x*Sp;% two-norm to be minimized by fsolve
options = optimoptions('fsolve','Display','off'); % some options
[x,fval,exitflag,output] = fsolve(fun,x0,options);
%--------------------------------------------------------------------------
% find real concentrations
u = np^2+ns*np;% number of unknowns
c0 = 10*ones(u,1);% create initial guesses
A = zeros(e(2),u);
A(:,1:2)=Sp';
A(:,2) = -A(:,2);
A(:,3:4)=Sp';
A(:,2) = -A(:,3);%some constraints are not used
b = 1*ones(e(2),1);
Aeq = [];
beq = [];
nonlcon = [];% no nonlinear constraints

lb = zeros(u,1);%upper and lower bounds for concentrations
ub = ones(u,1);
lb(1:np^2,1) = -Inf;% the first np^2 components of C_real are the components 
ub(1:np^2,1) = Inf;% of T and do not need bounds 

%finding real concentrations
[C,d] = fun1(x,np,ns);
options = optimoptions('lsqlin','Algorithm','interior-point');
[c,resnorm, residual,exitflag1] = lsqlin(C,d,A,b,Aeq,beq,lb,ub,c0,options);
c1 = c;
%--------------------------------------------------------------------------
% create invertible matrix T from solver results
for i = 1:np^2
    if i <=np
        T(i,1)=c(i);
    elseif i>np && i <=2*np
        T(i-np,2) = c(i);
    elseif i>2*np
        T(i-2*np,3) = c(i);
    end
end 
%--------------------------------------------------------------------------
% create matrix with real concentrations from solver results
for i = 1:(np*ns)
    if i<=ns
    cf(i,1) = c(np^2+i);
    elseif i>ns<=2*ns
    cf(i-ns,2) = c(np^2+i);
    else i>=2*ns;
    cf(i-2*ns,3) = c(np^2+i);
    end
end
% -------------------------------------------------------------------------
% calculating the fitted spectra
Yf = x*Sp; % alternatively cf*T^-1*Sp
% % spectra of the real components found
S_real = T^-1*Sp;
%--------------------------------------------------------------------------
%plotting
figure(1)
hold on
% plotting the real components:
p = plot(E,S_real(1,:),'--',E,S_real(2,:),'--','LineWidth',1);
p(1).Color = 'red';
p(2).Color = 'blue';
% plotting all spectra:
for i = 1:ns
plot(E,Y(i,:))
end
% plotting the fits obtained for the He treated and most carburized spectra
p1 = plot(E,Yf(nox,:),':',E,Yf(ncarb,:),':','LineWidth',1);
p1(1).Color = 'red';
p1(2).Color = 'blue';
% Style settings for the plotting:
xlabel('\bf \it energy \rm / eV','FontSize',16,'FontName','Calibri')
ylabel('\bf \it normalized absorption \rm / \mu (E)','FontSize',16,'FontName','Calibri')
% creating the labels:
labels = ['component 1   ';'component 2   '];
for i = 1:ns
    if i<10
    labels(i+2,:) = ['spectra ',int2str(i),'     '];
    else
        labels(i+2,:) = ['spectra ',int2str(i),'    '];
    end
end
labels(ns+3,:) = ['fit oxidized  '];
labels(ns+4,:) = ['fit carburized'];
lgd = legend(labels);
% more style settings:
lgd.FontName = 'Calibri';
lgd.FontSize =16;
legend('Location','northwest')
legend('boxoff')
axis([7090 7180 -0.01 1.4])
ax = gca;
c = ax.LineWidth;
ax.LineWidth = 2;
ax.TickDir = 'out';
ax.FontSize = 14;
ax.FontName = 'Calibri';
%--------------------------------------------------------------------------
% Optional: fitting all read in spectra with found components S_real
%{
x10 = -ones(ns1,np);
fun = @(x1)Y2 - x1*S_real;
options = optimoptions('fsolve','Display','off');
[x1,fval,exitflag,output] = fsolve(fun,x10,options);
%}
%--------------------------------------------------------------------------
% Optional: plotting concentrations obtained using the results of fitting of all spectra that 
% are read in using the found real components S_real
%{
figure(2)
plot(x1,'o')
xlabel('spectra')
ylabel('composition')
%}

%--------------------------------------------------------------------------
% Optional: plotting the result of fitting of all spectra that are read in using the 
% found real components S_real
%{
figure(4)
hold on
Yfitall = x1*S_real;
plot(E,Yfitall(nox,:),'red')
plot(E,Y2(nox,:),'green')
plot(E,Yfitall(ncarb,:),'red')
plot(E,Y2(ncarb,:),'green')
%}

%--------------------------------------------------------------------------
% Optional: plotting all components found initially
%{
% figure(5)
% hold on
% for i = 1:ns
% plot(E,S(i,:))
% end
% legend('all components found by SVD')
%}
%--------------------------------------------------------------------------
% definition of functions for lsqnonlin
    function [C,d] = fun1(x,np,ns)
    C = zeros(np*ns+ns+4,np^2+ns*np);
    y = np^2+1;
    l = 1;
    for i = 1:np*ns
        if mod(i,2) == 0
            C(i,y) = 1;
            C(i,3)= x(l,1);
            C(i,4) = x(l,2);
            d(i,1) = 1;
%         elseif mod(i,3) == 0 %only used for three components
%             F(i) = 1-c(y+ns)-c(y) -(x(l,1)*c(7)+x(l,2)*c(8)+x(l,3)*c(9));
            y = y +1;
            l = l+1;
        elseif mod(i,2) == 1
            C(i,y) = 1;
            C(i,1) = -x(l,1);
            C(i,2) = -x(l,2);
            d(i,1) = 0;
        end
    end
    for i = 1:ns
        C(np*ns+i,ns+np^2+i) = 1;
        C(np*ns+i,np^2+i) = 1;
        d(np*ns+i,1) = 1;
    end
    %additional constraints
    C(np*ns+ns+1,np^2+4) =  1;
    d(np*ns+ns+1,1) = 0.850;
    
    C(np*ns+ns+2,np^2+6) = 1;
    d(np*ns+ns+2,1) = 0.905;
    
    C(np*ns+ns+3,np^2+5) = 1;
    d(np*ns+ns+3,1) = 0.74;
    
    C(np*ns+ns+4,np^2+7) = 1;
    d(np*ns+ns+4,1) = 1;
    end



