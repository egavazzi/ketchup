% ketchup_b6plot plots some things the b6 version of ketchup has
% written into its output files.
%
% HG 2013-02-25

print_fig = 0;
save_input = 0;
plot_distrib = 1;
load outp/g.mat

%% Calculating epsilon_r
disp('---------------------------')
mass_e = 9.10938215e-31;
mass_i = 1.672621637e-27;
epsilon_0 = 8.85418781762039e-12;
e_charge = 1.602176487e-19;

% This is to get the variables of interest from the file inputb6
% dt = finding_var('dt')
% const_a = finding_var('const_a')
% nbr_iterations = finding_var('Niter')

dt = 2e-5
const_a = 100
nbr_iterations = 800000

n0R = 1.0e9;
n0L = 3.0e5;
t = nbr_iterations*dt

omega_p2 = (n0R*e_charge^2/mass_e/epsilon_0) ...
          + (n0R*e_charge^2/mass_i/epsilon_0) ...
          + (n0L*e_charge^2/mass_e/epsilon_0) ...
          + (n0L*e_charge^2/mass_i/epsilon_0);

omega_p = sqrt(omega_p2);
epsilon_r = (const_a*omega_p*dt)^2


new_period = 2*pi*(const_a*dt)
%2*pi./omega_p*sqrt(epsilon_r) %same result

%% Hard-coded output-directory
if print_fig
  FigDir = sprintf('Gunell_et_al_2013/%.2e_%.0fs',epsilon_r,t)
%   FigDir = sprintf('Figures-%s',datestr(now,'yyyymmdd'))
  S1 = (fileparts(which('ketchup_b6plot.m')));
  results_dir = fullfile(S1,'../../',FigDir);
  mkdir(results_dir)
  mkdir(fullfile(results_dir,'fzvzmu'))
  mkdir(fullfile(results_dir,'fzvz'))
end
%% Saving the inputs in the output-directory
% Allow for future analysis
if save_input
  copyfile('inputb6.m',results_dir)

  fileID = fopen(fullfile(results_dir,'epsilon.txt'),'w');
  fprintf(fileID,'Epsilon_r = %d \n',epsilon_r);
  fprintf(fileID,'dt = %d \n',dt);
  fprintf(fileID,'a = %d \n',const_a);
  fclose(fileID);
end
%% Setting time steps to plot
t_to_plot = [40]; % CHOOSE TIME TO PLOT (in s)

load outp/Efield.mat
pt=[timestepsEfield timestepsEfield(end)+dt]*dt;
[~,index_t_to_plot] = min(abs(pt(:)-t_to_plot));

%% --- Efield and Vp --- %
load outp/Efield.mat
load outp/current.mat
clear legend_str

potentialmatrix =[zeros(length(timestepsEfield),1) ...
               -cumsum(Efieldmatrix.*(ones(length(timestepsEfield),1)*dz),2)];

colours ='brkcm';
figure()
set(gcf,'paperpositionmode','auto')
set(gcf,'WindowState','maximized')


subplot(2,1,1)
plot(z,Efieldmatrix(index_t_to_plot,:),'linewidth',2)
set(gca,'fontname','times','fontsize',14)
ylabel('E (V/km)','fontname','times','fontsize',18)
xline(z(end),'--')
for i = 1:length(index_t_to_plot)
  legend_str(i) = {[num2str(timestepsEfield(index_t_to_plot(i))*dt) 's']};
end
legend(legend_str,'location','northwest')
xlim([4e7 5.5e7])
% ylim([0 2e-4])
grid on
xticklabels('')
xticks(0:0.5e7:5.5e7)


subplot(2,1,2)
plot(zcorn,potentialmatrix(index_t_to_plot,:),'linewidth',2)
set(gca,'fontname','times','fontsize',14)
xlabel('z (m)','fontname','times','fontsize',18)
ylabel('Potential (V)','fontname','times','fontsize',18)
xline(z(end),'--')
for i = 1:length(index_t_to_plot)
  legend_str(i) = {[num2str(timestepsEfield(index_t_to_plot(i))*dt) 's']};
end
legend(legend_str,'location','northwest')
sgtitle(['\epsilon_r = ',num2str(epsilon_r,3)])
grid on
xlim([0 5.5e7])
xticks(0:0.5e7:5.5e7)


if print_fig
  print('-dpng','-painters',fullfile(results_dir,'E_and_Vp_over_z.png'))
  print('-depsc2','-painters',fullfile(results_dir,'E_and_Vp_over_z.eps'))
end
%% --- Density --- %
load outp/density.mat

figure();clf
set(gcf,'paperpositionmode','auto')
set(gcf,'WindowState','maximized')
sn = size(densitymatrix);
if length(sn)==3
  for jj=1:sn(1)
    subplot(sn(1),1,jj)
    hold on
    for ii=index_t_to_plot(1:1:length(index_t_to_plot))
        semilogy(z,densitymatrix(jj,:,ii),'linewidth',2)
    end
  end

  xlabel('z','fontname','times','fontsize',18)
  for jj=1:sn(1)
    subplot(sn(1),1,jj)
    set(gca,'fontname','times','fontsize',14)
    ylabel('n','fontname','times','fontsize',18)
    title(['Species ' num2str(jj)],'fontname','times','fontsize',16)
    axis([4e7 5.5e7 -inf inf],'autoy')
    xline(z(end),'--')
    for i = 1:length(index_t_to_plot)
      legend_str(i) = {[num2str(timestepsdensity(index_t_to_plot(i))*dt) 's']};
    end
    legend(legend_str,'location','northwest')
  end
else
  plot(z,densitymatrix)
  set(gca,'yscale','log','fontname','times','fontsize',14)
  xlabel('z','fontname','times','fontsize',18)
  ylabel('n','fontname','times','fontsize',18)
  grid on
end
sgtitle(['\epsilon_r = ',num2str(epsilon_r,3)])

if print_fig
  print('-dpng','-painters',fullfile(results_dir,'Densities_over_z.png'))
  print('-depsc2','-painters',fullfile(results_dir,'Densities_over_z.eps'))
end
%% --- Current --- %
load outp/current.mat

figure(3)
set(gcf,'paperpositionmode','auto')
set(gcf,'WindowState','maximized')
plot(z,currentmatrix(index_t_to_plot,:),'linewidth',2)
set(gca,'fontname','times','fontsize',14)
ylabel('i','fontname','times','fontsize',18)
xlabel('z','fontname','times','fontsize',18)
% axis([z(1) z(end) -inf inf],'autoy')
xline(z(end),'--')
for i = 1:length(index_t_to_plot)
  legend_str(i) = {[num2str(timestepscurrent(index_t_to_plot(i))*dt) 's']};
end
legend(legend_str,'location','northwest')
title(['\epsilon_r = ',num2str(epsilon_r,3)])

if print_fig
  print('-dpng','-painters',fullfile(results_dir,'Current_over_z.png'))
  print('-depsc','-painters',fullfile(results_dir,'Current_over_z.eps'))
end
%% --- Efield(z,t) --- %
load outp/Efield.mat
pt=[timestepsEfield timestepsEfield(end)+dt]*dt;

figure()
set(gcf,'paperpositionmode','auto','renderer','zbuffer')
% set(gcf,'WindowState','maximized')
L=zmax-zmin;
% if ~exist('hidefraction')
  hidefraction=0.0;
% end
zrange = (z>zmin+hidefraction*L & z<zmax-hidefraction*L);
trange = (pt(1:end-1)>0& pt(1:end-1)<300);
pp=(Efieldmatrix(trange,zrange).').*1000;
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
hh.Label.String = 'E [V/km]';
hh.Label.FontSize = 18;
colormap(jet);
caxis([-2 1])
xlabel('t  [s]','fontname','times','fontsize',18)
ylabel('z  [km]','fontname','times','fontsize',18)
grid off
% title(['\epsilon_r = ',num2str(epsilon_r,3)])
% 
% 
% ylim([0 1000])
% caxis([-0.01 0.01])

if print_fig
print('-dpng','-painters',fullfile(results_dir,'Efield_over_z&t.png'))
%print('-depsc','-painters',fullfile(results_dir,'Efield_over_z&t.eps'))
end
%% --- Distribution function --- %
if plot_distrib
  dd=dir('outp');
  timestepsdistr=[];
  for ii = 1:length(dd)
    if length(dd(ii).name)>=17
      if strcmp(dd(ii).name(1:6),'fzvzmu') & strcmp(dd(ii).name(14:17),'.mat')
        timestepsdistr=[timestepsdistr str2num(dd(ii).name(7:13))];
      end 
    end
  end

  load outp/Bfield.mat

  if ~exist('drawdistr')
    drawdistr = logical(1);
  end


  if drawdistr
  %   figure(5);clf
  %   set(gcf,'paperpositionmode','auto')
  %   finitemass=[];
  %   for ii=1:Nspecies
  %     if ~isnan(particle(ii).mass) & ~isinf(particle(ii).mass)
  %       finitemass=[finitemass ii];
  %     end
  %   end


    for ii=5%:10:length(timestepsdistr)
      distrfile=['outp/fzvzmu' num2str(timestepsdistr(ii),'%0.7i') '.mat'];
      load(distrfile)

      figure(1000+ii)
      set(gcf,'WindowState','maximized')
      ketchup_b6fplot(particle,fzvzmustruct,Nz,zcorn,z,dz,Nspecies)
      sgtitle({['\epsilon_r = ',num2str(epsilon_r,3)];['t = ',num2str(timestepsdistr(ii)*dt),'s']})
      if print_fig
        print(sprintf('-f10%02d',ii),'-dpng','-painters', ...
              fullfile(results_dir,'fzvzmu',sprintf('testfzvzmu_10%02d.png',ii)))
      end
  %     print(sprintf('-f10%02d',ii),'-depsc','-painters', ...
  %           fullfile(results_dir,'fzvzmu',sprintf('testfzvzmu_10%02d.eps',ii)))       

      figure(2000+ii)
      set(gcf,'WindowState','maximized')
      ketchup_b6fzvzplot(particle,fzvzmustruct,Nz,zcorn,z,dz,Nspecies)
      sgtitle({['\epsilon_r = ',num2str(epsilon_r,3)];['t = ',num2str(timestepsdistr(ii)*dt),'s']})
      if print_fig
        print(sprintf('-f20%02d',ii),'-dpng','-painters', ...
              fullfile(results_dir,'fzvz',sprintf('testfzvz_20%02d.png',ii)))
      end
  %     print(sprintf('-f20%02d',ii),'-depsc','-painters', ...
  %           fullfile(results_dir,'fzvz',sprintf('testfzvz_20%02d.eps',ii)))

  %     figure(5)
  %     kTz=ketchupTz(particle,fzvzmustruct,Nz,zcorn,z,dz,Nspecies);
  %     kTp=ketchupTp(particle,fzvzmustruct,Nz,zcorn,z,dz,Nspecies,B);
  %     stz=size(kTz);
  %     for jj=1:stz(1)
  %       subplot(stz(1),3,(jj-1)*3+1)
  %       hold on
  %       plot(z,kTz(jj,:),colours(mod(ii-1,5)+1))
  %       hold off
  %       title(['species ' num2str(finitemass(jj))], ...
  %             'fontname','times','fontsize',18)
  %       subplot(stz(1),3,(jj-1)*3+2)
  %       hold on
  %       plot(z,kTp(jj,:),colours(mod(ii-1,5)+1))
  %       hold off
  %       subplot(stz(1),3,(jj-1)*3+3)
  %       hold on
  %       plot(z,(kTz(jj,:)+2*kTp(jj,:))/3,colours(mod(ii-1,5)+1))
  %       hold off
  %       title(['species ' num2str(finitemass(jj))], ...
  %             'fontname','times','fontsize',18)
  %     end
    end

  %   figure(5)
  %   for jj=1:stz(1)
  %     subplot(stz(1),3,(jj-1)*3+1)
  %     set(gca,'fontname','times','fontsize',14)
  %     ylabel('k_{B}T_{z}/e','fontname','times','fontsize',18)
  %     xlabel('z','fontname','times','fontsize',18)
  %     title(['species ' num2str(finitemass(jj))], ...
  %           'fontname','times','fontsize',18)
  %     grid on
  %     aa=axis;axis([min(z) max(z) aa(3:4)]);
  %     subplot(stz(1),3,(jj-1)*3+2)
  %     grid on
  %     aa=axis;axis([min(z) max(z) aa(3:4)]);
  %     set(gca,'fontname','times','fontsize',14)
  %     ylabel('k_{B}T_{\perp}/e','fontname','times','fontsize',18)
  %     xlabel('z','fontname','times','fontsize',18)
  %     title(['species ' num2str(finitemass(jj))], ...
  %           'fontname','times','fontsize',18)
  %     subplot(stz(1),3,(jj-1)*3+3)
  %     grid on
  %     aa=axis;axis([min(z) max(z) aa(3:4)]);
  %     set(gca,'fontname','times','fontsize',14)
  %     ylabel('k_{B}(T_{z} + 2T_{\perp})/(3e)','fontname','times','fontsize',18)
  %     xlabel('z','fontname','times','fontsize',18)
  %     title(['species ' num2str(finitemass(jj))], ...
  %           'fontname','times','fontsize',18)
  %   end
  end
end
%% --- Potential --- %
figure(6)
set(gcf,'WindowState','maximized')
plot(zcorn,potentialmatrix(end,:),'k','linewidth',2)
set(gca,'fontname','times','fontsize',14)
xlabel('z','fontname','times','fontsize',18)
ylabel('V_{p}','fontname','times','fontsize',18)
xline(z(end),'--')
title(['\epsilon_r = ',num2str(epsilon_r,3)])

if print_fig
  print('-dpng','-painters',fullfile(results_dir,'Vp_over_z.png'))
  print('-depsc','-painters',fullfile(results_dir,'Vp_over_z.eps'))
end
% % Characteristic frequencies
% omegap2=0;
% for ii=1:Nspecies
%   if ~isnan(particle(ii).mass)
%     omegap2 = omegap2 + elc^2 * particle(ii).n0 /(particle(ii).mass*eps0);
%   end
% end
% omegap0 = sqrt(omegap2);
% epsilon_r = max((omegap0*const_a*dt)^2,1);
% nlast=densitymatrix(:,:,end);
% masses = zeros(Nspecies,1);
% for ii=1:Nspecies
%   masses(ii)=particle(ii).mass;
% end
% mfin=~isnan(masses);
% fp = sqrt(sum(elc^2*nlast(mfin,:)./ ...
%               (eps0*epsilon_r*masses(mfin)*ones(1,Nz)),1))/(2*pi);
% fpe = fp; % approximately at least
% fpi = fpe*sqrt(me/mp);
% fge = elc*B.'/(me*2*pi);
% flh = fpi./sqrt(1+fpe.^2./fge.^2);
% fuh = sqrt(fpe.^2 + fge.^2);
% 
% % Debye length
% lambdaD=sqrt(eps0*epsilon_r*kTz./nlast(mfin,:)/elc);
% if Nspecies>=3
%   lambdaDe=1./sqrt(sum(lambdaD([1 3],:).^(-2)));
% else
%   lambdaDe=lambdaD(1,:);
% end
% lambdaDtot=1./sqrt(sum(lambdaD.^(-2)));
% 
% % $$$ B=B(:).';dB=dB(:).';
% % $$$ lambdaDplus=0.5*(dB./B)*eps0*epsilon_r.*kTz./nlast(mfin,:)/elc + ...
% % $$$         sqrt(0.25*((dB./B)*eps0*epsilon_r.*kTz./nlast(mfin,:)/elc).^2 ...
% % $$$              + eps0*epsilon_r*kTz./nlast(mfin,:)/elc);
% % $$$ lambdaDminus=0.5*(dB./B)*eps0*epsilon_r.*kTz./nlast(mfin,:)/elc - ...
% % $$$         sqrt(0.25*((dB./B)*eps0*epsilon_r.*kTz./nlast(mfin,:)/elc).^2 ...
% % $$$              + eps0*epsilon_r*kTz./nlast(mfin,:)/elc);
% 
% if ~exist('drawfft')
%   drawfft = logical(1);
% end
% 
% if drawfft
%   figure(7);clf
%   set(gcf,'paperpositionmode','auto','renderer','zbuffer')
%   Nfft = 256;
%   if length(timestepsEfield)<Nfft
%     Nfft=2^(floor(log(length(timestepsEfield))/log(2)));
%   end
%   Powermatrix = [];
%   fs=1/(dump_period_fields*dt);
%   for ii = 1:Nz
%     E=Efieldmatrix(:,ii);
%     [Pxx,Fscale] = pwelch(E-mean(E),hanning(Nfft), ...
%                      round(0.65*Nfft),Nfft,fs);
%     Powermatrix = cat(2,Powermatrix,Pxx);
%   end
%   pp = Powermatrix;
%   ss = size(pp);
%   pp = [[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
%   maxpp=max(max(pp));
%   pz = zcorn;
%   pf = [Fscale; 2*Fscale(end)-Fscale(end-1)];
%   surf(pz,pf,log10(pp))
%   clear pp
%   set(gca,'fontname','times','fontsize',14)
%   view(2)
%   axis([min(pz) max(pz) min(pf) max(pf)])
%   shading flat
%   hh=colorbar;
%   set(hh,'fontname','times','fontsize',14)
%   xlabel('z','fontname','times','fontsize',18)
%   ylabel('f','fontname','times','fontsize',18)
%   grid off
%  
%   % draw a plasma frequency curve  
%   hold on
%   plot3(z,fp,ones(size(fp))*log10(maxpp*2),'w')
%   plot3(z,fpi,ones(size(fpi))*log10(maxpp*2),'k')
% % $$$   plot3(z,flh,ones(size(flh))*log10(maxpp*2),'k')
% % $$$   plot3(z,fge,ones(size(fge))*log10(maxpp*2),'g')
% % $$$   plot3(z,fuh,ones(size(fuh))*log10(maxpp*2),'r')
%   hold off
% end
