% make videos with the data from ketchup
%

print_fig = 0;
save_input = 0;
plot_distrib = 1;
load outp/g.mat

%% --- Efield, over altitude + z-t diagram --- %
load outp/Efield.mat
load outp/current.mat
clear legend_str

clear E_integral
% for ii = 2:length(z)
%   E_integral(:,ii) = trapz(z(1:ii),Efieldmatrix(:,1:ii),2);
% end
% for ii = length(zcorn):-1:1
%   E_integral(:,ii) = trapz(zcorn(end:-1:ii),Efieldmatrix(:,end:-1:ii),2);
% end
pt=[timestepsEfield timestepsEfield(end)+dt]*dt;
E_integral(:,:) = - cumsum(Efieldmatrix.*dz,2,'reverse'); % in (V)

colours ='brkcm';
figure()
set(gcf,'paperpositionmode','auto')
figPHigh = [1014, 226, 1000, 800];
set(gcf,'position',figPHigh)
set(gcf,'WindowState','maximized')

doMovie = 0;

if doMovie
  try
    vidObj = VideoWriter('Integral_Efield_over_iono.avi');
  %   vidObj.Quality = 100;
    open(vidObj);
    doMovie = 1;
  catch
    doMovie = 0;
    disp(['Failed to initialize VideoWriter Object for movie: ',movieOut])
  end
end

t = tiledlayout(1,3);

% z-t diagram initialization
nexttile(2,[1 2])
L=zmax-zmin;
hidefraction=0.00;
zrange = (z>zmin+hidefraction*L & z<zmax-hidefraction*L);
trange = (pt(1:end-1)>0& pt(1:end-1)<100);
pp=(Efieldmatrix(trange,zrange).')*1000;
ss=size(pp);
pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
pz=zcorn(zrange); pz = [pz (pz(end)+dz(length(pz)))];
pt=[timestepsEfield timestepsEfield(end)+dt]*dt; %redondant (already done in "time to plot" section)
dt_timestepsEfield = abs(dt*(timestepsEfield(1)-timestepsEfield(2)));
pt2 = pt(trange); pt2 = [pt2(1)-dt_timestepsEfield pt2];

% pp = pp([end:-1:1],:);
% pz = pz([end:-1:1]);
pz = 5.2e7 - pz + 400e3;

surf(pt2,pz/1000,pp)
clear pp
set(gca,'fontname','times','fontsize',14)
view(2)
% axis([min(pt2) max(pt2) 4.75e7 5.195e7])
axis([min(pt2) max(pt2) 400 5000])
shading flat
hh=colorbar;
set(hh,'fontname','times','fontsize',14)
hh.Label.String = 'E [V/km]';
hh.Label.FontSize = 18;
colormap(jet);
caxis([-0.04 0.04]);
xlabel('t  [s]','fontname','times','fontsize',18)
% ylabel('z  [m]','fontname','times','fontsize',18)
yticklabels('')
grid off
X = xline(0);

z = 5.2e7 - z + 400e3;

for i_t = 1:numel(pt)-1
  % plot of E over space
  nexttile(1)
  
  plot((Efieldmatrix(1,(1:end)) - Efieldmatrix(i_t,(1:end))).*1000,z(1:end)/1000,'linewidth',2)
  set(gca,'fontname','times','fontsize',14)
  ylabel('z [km]','fontname','times','fontsize',18)
  xlabel('\DeltaE [V/km]','fontname','times','fontsize',18)
  legend('Efield','location','northwest')
  ylim([400 3000])
  xlim([-1 1])
  grid on
   
%   plot(E_integral(i_t,(1:end)),z(1:end)/1000,'linewidth',2)
%   set(gca,'fontname','times','fontsize',14)
%   ylabel('z [km]','fontname','times','fontsize',18)
%   xlabel('E [V/km]','fontname','times','fontsize',18)
%   legend('Efield','location','northwest')
%   ylim([400 5000])
%   grid on

%   plot(Efieldmatrix(i_t,(1:end)).*1000,z(1:end)/1000,'linewidth',2)
%   set(gca,'fontname','times','fontsize',14)
%   ylabel('z [km]','fontname','times','fontsize',18)
%   xlabel('E [V/km]','fontname','times','fontsize',18)
%   legend('Efield','location','northwest')
%   ylim([400 2000])
%   xlim([-1 10])
%   grid on

  
  nexttile(2,[1 2])
  delete(X)
  X = xline(pt(i_t));
  
  title(['t =',num2str(pt(i_t),'%.3f'),' s'],'FontSize',16);

  drawnow
  if doMovie == 1
    writeVideo(vidObj,getframe(gcf));
  end
end
close(vidObj);
output = [];

%% --- Efield, z-t diagram --- %

load outp/Efield.mat
pt=[timestepsEfield timestepsEfield(end)+dt]*dt;

figure()
set(gcf,'paperpositionmode','auto')
% set(gcf,'WindowState','maximized')
figPHigh = [1014, 226, 875, 648];
set(gcf,'position',figPHigh)


L=zmax-zmin;
hidefraction=0.01;
zrange = (z>zmin+hidefraction*L & z<zmax-hidefraction*L);
trange = (pt(1:end-1)>0& pt(1:end-1)<100);
pp=(Efieldmatrix(trange,zrange).').*1000;
ss=size(pp);
pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
pz=zcorn(zrange); pz = [pz (pz(end)+dz(length(pz)))];
pt=[timestepsEfield timestepsEfield(end)+dt]*dt; %redondant (already done in "time to plot" section)
dt_timestepsEfield = abs(dt*(timestepsEfield(1)-timestepsEfield(2)));
pt2 = pt(trange); pt2 = [pt2(1)-dt_timestepsEfield pt2];
surf(pt2,pz,pp)
clear pp
set(gca,'fontname','times','fontsize',14)
view(2)
axis([min(pt2) max(pt2) 4.75e7 max(pz)])
shading flat
hh=colorbar;
set(hh,'fontname','times','fontsize',14)
hh.Label.String = 'E [mV/m]';
hh.Label.FontSize = 18;
colormap(jet);
caxis([-0.04 0.04]);
xlabel('t  [s]','fontname','times','fontsize',18)
ylabel('z  [m]','fontname','times','fontsize',18)
grid off

doMovie = 0;
try
  vidObj = VideoWriter('zt_diagram.avi');
  open(vidObj);
  doMovie = 1;
catch
  doMovie = 0;
  disp(['Failed to initialize VideoWriter Object for movie: ',movieOut])
end

X = xline(0)
for i_t = 1:numel(pt)-1
  delete(X)
  X = xline(pt(i_t));
  
  drawnow
  if doMovie == 1
    writeVideo(vidObj,getframe(gcf));
  end
end
close(vidObj);
output = [];

