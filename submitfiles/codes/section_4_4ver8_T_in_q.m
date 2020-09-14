%This script calls the function PGD_TransientHeat_TX which implements the
%solution of the transient heat problem 4.1 in order to generate all the
%figures of section 4.4
%
%Copyright (c) 2013, Francisco Chinesta (Ecole Centrale de Nantes), Roland Keunings (Universite catholique de Louvain), Adrien Leygue(CNRS)
%Author: Adrien Leygue.
%All rights reserved.
%See License file for more details

close all
clear
%domain size & meshes
Lx = 20;
Nx = 101;
Lt = 20;
% Nt = 1001;
x = linspace(0,Lx,Nx)';

% [T_in,Nt] = input_T();
load('20200303_T_in4.mat');
Nt = size(T_in,1);
t = linspace(0,Lt,Nt)';

%PGD parameters
%Fixed point tolerance
epsilon = 1e-8;
%PGD enrichment tolerance set to 0 to enforce the number of enrichments
epsilon_tilde = 0;
%Maximum dumber of enrichments
Max_terms = 10;
%Maximum number of iterations in the fixed point loop
Max_fp_iter = 500;


%Source term in separated form
Fx = zeros(size(x));
Ft = zeros(size(t));

%thermal diffusivity
k = 98.8;

%heat flux
prompt = {'Heat flux q:'};
dlgtitle = 'Input heat flux';
dims = [1 35];
% definput = {'Heat flux q:'};
q1 = inputdlg(prompt,dlgtitle,dims);%,definput);
q1 = str2num(q1{1});
% q1 = -55000;%input('q1?');
% q1 = q1*10^-6;
q = q1*ones(Nt,1);

%Reference solution Eq. (4.35)
% U_ex = reference_solution_4_35(t,x,k,200);

%PGD solution
[X,T] = PGD_TransientHeat_TXver8_FEM(t,T_in,x,Max_terms,Max_fp_iter,epsilon,epsilon_tilde,Ft,Fx,k,q);

%Normalization of the computed enrichment terms
Xnorm = bsxfun(@rdivide,X,sqrt(trapz(x,X.^2)));
Tnorm = bsxfun(@rdivide,T,sqrt(trapz(t,T.^2)));

%Plot X_i
%Figure 4.1
figure(1);
subplot(3,2,1);
set(gca,'fontsize',14);
handles=plot(x,Xnorm(:,1:4));
legend(handles,cellfun(@(in) ['X_' num2str(in) '(x)'],num2cell(1:4),'uniformoutput',false))
xlabel('$x$','interpreter','latex')
ylabel('$\frac{X_i(x)}{\|X_i(x)\| \,\,\, \,\,}$','interpreter','latex','fontsize',18)

%Plot T_i
%Figure 4.2
% figure;
subplot(3,2,2);
set(gca,'fontsize',14);
handles=plot(t,Tnorm(:,1:4));
legend(handles,cellfun(@(in) ['T_' num2str(in) '(t)'],num2cell(1:4),'uniformoutput',false))
xlabel('$t$','interpreter','latex')
ylabel('$\frac{T_i(t)}{\|T_i(t)\| \,\,\, \,\,}$','interpreter','latex','fontsize',18)

%Computation of the error berween the reference and PGD solution
% E_N = zeros(1,Max_terms);
for j=1:(Max_terms)
    titlename = sprintf('term = %d',j-1);
    %reconstruction of the PGD solution
    U_pgd = X(:,1:j)*T(:,1:j)';
    % E_N(j) = trapz(t,trapz(x,(U_ex - U_pgd).^2));
    figure(2);
    subplot(2,5,j)
    mesh(U_pgd(:,2:end))
    view(45,45)
    xlabel('T')
    ylabel('X')
    xlim([0 Nt+4])
    ylim([0 Nx+4])
%     zlim([0 inf])
    title(titlename)

    %reconstruction of the PGD solution
    
    U_pgd2(:,:,j) = X(:,j)*T(:,j)';
    % E_N(j) = trapz(t,trapz(x,(U_ex - U_pgd).^2));
    figure(3);
    subplot(2,5,j)
    mesh(U_pgd2(:,2:end,j))
    view(45,45)
    xlabel('T')
    ylabel('X')
    xlim([0 Nt+4])
    ylim([0 Nx+4])
%     zlim([0 inf])
    title(titlename)

end

figure(1);
subplot(3,2,3);
plot(U_pgd(end,2:end));
title('U-pgd(x:end)')
xlabel('t')

% figure;
subplot(3,2,4);
plot(U_pgd(1,2:end));
title('U-pgd(x:1)')
xlabel('t')

subplot(3,2,5)
plot(U_pgd(2:end,1));
title('U-pgd(t:1)')
xlabel('x')

subplot(3,2,6)
plot(U_pgd(2:end,end));
title('U-pgd(t:end)')
xlabel('x')

figure;
[tm,xm] = meshgrid(t(2:end),x);
mesh(tm,xm,U_pgd(:,2:end))
set(gca,'fontsize',18);
xlabel('T(s)')
ylabel('X(mm)')
xlim([0 Lt+1])
ylim([0 Lx+1])
% title('U-pgd')
view(45,45)

ind = [1 2 3 10];
for j=1:4
figure(6);

titlename = sprintf('N = %d',ind(j));
subplot(2,2,j)
mesh(tm,xm,U_pgd2(:,2:end,ind(j)))
set(gca,'fontsize',12);
xlabel('T(s)')
ylabel('X(mm)')
xlim([0 Lt+1])
ylim([0 Lx+1])
% title('U-pgd')
view(45,45)
title(titlename)
end

% figure;
% for i=2:10
%     plot(U_pgd(:,1+100*(i-1)))
%     hold on
% end
% legend('2','3','4','5','6','7','8','9','10')
% hold off

makegraph;
