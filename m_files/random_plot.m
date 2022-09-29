%% z-t diagram of density
figure()
set(gcf,'paperpositionmode','auto','renderer','zbuffer')
% set(gcf,'WindowState','maximized')
L=zmax-zmin;
% if ~exist('hidefraction')
  hidefraction=0.0;
% end
zrange = (z>zmin+hidefraction*L & z<zmax-hidefraction*L);
trange = (pt(1:end-1)>0& pt(1:end-1)<300);
pp=(squeeze(densitymatrix(1,zrange,trange)));
ss=size(pp);
pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
pz=zcorn(zrange); pz = [pz (pz(end)+dz(length(pz)))]; pz = (5.2e7 - pz)/1e3;
dt_timestepsEfield = abs(dt*(timestepsEfield(1)-timestepsEfield(2)));
pt = pt(trange); pt = [pt(1)-dt_timestepsEfield pt];
surf(pt,pz,pp)
clear pp
set(gca,'fontname','times','fontsize',14)
view(2)
axis([min(pt) max(pt) min(pz) max(pz)])
% axis([min(pt2) max(pt2) 4.2e7 max(pz)])
shading flat
hh=colorbar;
set(hh,'fontname','times','fontsize',14)
hh.Label.String = 'n [m^{-3}]';
hh.Label.FontSize = 18;
colormap(jet);
% caxis([-2 1])
xlabel('t  [s]','fontname','times','fontsize',18)
ylabel('z  [km]','fontname','times','fontsize',18)
grid off

%% Animation of density + Efield
load outp/density.mat
load outp/Efield.mat
pt=[timestepsEfield timestepsEfield(end)+dt]*dt;
potentialmatrix =[zeros(length(timestepsEfield),1) ...
               -cumsum(Efieldmatrix.*(ones(length(timestepsEfield),1)*dz),2)];

z = 5.2e7 - z + 400e3;
pt=[timestepsEfield timestepsEfield(end)+dt]*dt;
ii = 1;

figure(6)
t = tiledlayout(1,2);

for i_t = 1:25:numel(pt)-1
  nexttile(1)
  plot(densitymatrix(ii,(1:end),1),z(1:end)/1000,'linewidth',2)
  hold on
  plot(densitymatrix(ii,(1:end),i_t),z(1:end)/1000,'linewidth',2)
  hold off
%   plot((densitymatrix(ii,(1:end),i_t)-densitymatrix(ii,(1:end),1))./densitymatrix(ii,(1:end),1),z(1:end)/1000,'linewidth',2)
%   plot((densitymatrix(ii,(1:end),i_t)-densitymatrix(ii,(1:end),1)),z(1:end)/1000,'linewidth',2)
  set(gca,'fontname','times','fontsize',14)
  ylabel('z [km]','fontname','times','fontsize',18)
  xlabel('n [m^{-3}]','fontname','times','fontsize',18)
  legend('Density','location','northwest')
  ylim([400 8000])
%   xlim([0 0.1])
%   xlim([0 3e6])
%   xlim([0 1e5])
  grid on
  title(['t =',num2str(pt(i_t),'%.3f'),' s'],'FontSize',16);
  drawnow;
  
  nexttile(2)
%   plot(potentialmatrix(1,(2:end)).*1000,z(1:end)/1000,'linewidth',2)
%   hold on
%   plot(potentialmatrix(i_t,(2:end)).*1000,z(1:end)/1000,'linewidth',2)
%   hold off
  plot(Efieldmatrix(1,(1:end)).*1000,z(1:end)/1000,'linewidth',2)
  hold on
  plot(Efieldmatrix(i_t,(1:end)).*1000,z(1:end)/1000,'linewidth',2)
  hold off
  set(gca,'fontname','times','fontsize',14)
  ylabel('z [km]','fontname','times','fontsize',18)
  xlabel('E [V/km]','fontname','times','fontsize',18)
  legend('Efield','location','northwest')
  ylim([4000 8000])
%   xlim([-2 0.5])
  grid on

end