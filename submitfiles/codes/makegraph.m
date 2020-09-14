
% load('U_pgd.mat')
load('20200303_temp_plate.mat')
t = linspace(0,20,1001);
t = t(2:end);
n = 4;%input('number 4 or 5');
Tmaxmin = [172.5158 169.7726; 171.033 167.544];
Tmax = Tmaxmin(1,n-3) + 273.15;
Tmin = Tmaxmin(2,n-3) + 273.15;

T = T';
T_thin = U_pgd(1,:)*(Tmax-Tmin)+Tmin;
T_pgd = U_pgd(end,:)*(Tmax-Tmin)+Tmin;
T_thick = T*(Tmax-Tmin)+Tmin;
figure;
plot(t,T_thin,'linewidth',1.5,'LineStyle','-.')
hold on
plot(t,T_pgd,'linewidth',1.5,'LineStyle','--')
plot(t,T_thick,'linewidth',1.5,'LineStyle','-')
hold off

xlabel('Time [s]')
ylabel('Temperature [K]')
set(gca,'fontsize',18)
ylim([443 446])
yticks([443 444 445 446])
legend('Plate surface','Calculated by PGD','Mold surface')

str = [sprintf('U(X=L)_q='), sprintf(num2str(q1))];print(gcf,'-dpng','-r600',[str,'.png']);

% str = [sprintf('pca'), sprintf('%d',j), sprintf('-2')];print(gcf,'-dpng','-r600',[str,'.png']);