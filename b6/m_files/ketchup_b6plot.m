% ketchup_b6plot plots some things the b6 version of ketchup has
% written into its output files.
%
% HG 2013-02-25

inputb6
load outp/g.mat

% construct xi-vectors
dxi=1/Nz;
xicorn=dxi*[0:Nz];
xi=0.5*(xicorn(1:end-1) + xicorn(2:end));


load outp/Efield.mat

potentialmatrix=[zeros(length(timestepsEfield),1) ...
               -cumsum(Efieldmatrix.*(ones(length(timestepsEfield),1)*dz),2)];

colours='brkcm';

figure(1)
set(gcf,'paperpositionmode','auto')
subplot(2,1,1)
plot(z,Efieldmatrix.')
set(gca,'fontname','times','fontsize',14)
ylabel('E','fontname','times','fontsize',18)
subplot(2,1,2)
plot(zcorn,potentialmatrix.')
set(gca,'fontname','times','fontsize',14)
xlabel('z','fontname','times','fontsize',18)
ylabel('Potential','fontname','times','fontsize',18)


% --- Density --- %
load outp/density.mat

figure(2);clf
set(gcf,'paperpositionmode','auto')
sn=size(densitymatrix);
if length(sn)==3
  for jj=1:sn(1)
    subplot(sn(1),1,jj)
    hold on
    for ii=1:sn(3)
      plot(z,densitymatrix(jj,:,ii),colours(mod(ii-1,5)+1))
    end
  end

  xlabel('z','fontname','times','fontsize',18)
  for jj=1:sn(1)
    subplot(sn(1),1,jj)
    set(gca,'fontname','times','fontsize',14)
    ylabel('n','fontname','times','fontsize',18)
    title(['Species ' num2str(jj)],'fontname','times','fontsize',16)
  end
else
  plot(z,densitymatrix)
  set(gca,'yscale','log','fontname','times','fontsize',14)
  xlabel('z','fontname','times','fontsize',18)
  ylabel('n','fontname','times','fontsize',18)
  grid on
end


% --- Current --- %
load outp/current.mat

figure(3)
set(gcf,'paperpositionmode','auto')
plot(z,currentmatrix.')
set(gca,'fontname','times','fontsize',14)
ylabel('i','fontname','times','fontsize',18)
xlabel('z','fontname','times','fontsize',18)


figure(4)
set(gcf,'paperpositionmode','auto','renderer','zbuffer')
L=zmax-zmin;
if ~exist('hidefraction')
  hidefraction=0.01;
end
zrange = (z>zmin+hidefraction*L & z<zmax-hidefraction*L);
pp=Efieldmatrix(:,zrange).';
ss=size(pp);
pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
pz=zcorn(zrange);pz=[pz (pz(end)+dz(length(pz)))];
pt=[timestepsEfield timestepsEfield(end)+dt]*dt;
surf(pt,pz,pp)
clear pp
set(gca,'fontname','times','fontsize',14)
view(2)
axis([min(pt) max(pt) min(pz) max(pz)])
shading flat
hh=colorbar;
set(hh,'fontname','times','fontsize',14)
xlabel('t  [s]','fontname','times','fontsize',18)
ylabel('z  [m]','fontname','times','fontsize',18)
grid off
% $$$ cax=caxis;Ncmap=64;Nneg=floor(Ncmap*(-cax(1))/diff(cax));
% $$$ cax(2)=cax(1)*(1-Ncmap/Nneg);caxis(cax)
% $$$ %colormap([summer(Nneg); jet(Ncmap-Nneg)])
% $$$ %colormap([summer(Nneg); fliplr(autumn(Ncmap-Nneg))])
% $$$ colormap([summer(Nneg); wshift(2,autumn(Ncmap-Nneg),[0 1])])


% --- distribution function --- %
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
  figure(5);clf
  set(gcf,'paperpositionmode','auto')
  finitemass=[];
  for ii=1:Nspecies
    if ~isnan(particle(ii).mass) & ~isinf(particle(ii).mass)
      finitemass=[finitemass ii];
    end
  end


  for ii=1:length(timestepsdistr)
    distrfile=['outp/fzvzmu' num2str(timestepsdistr(ii),'%0.7i') '.mat'];
    load(distrfile)
    figure(1000+ii)
    ketchup_b6fplot(particle,fzvzmustruct,Nz,zcorn,z,dz,Nspecies)
    figure(2000+ii)
    ketchup_b6fzvzplot(particle,fzvzmustruct,Nz,zcorn,z,dz,Nspecies)
    figure(5)
    kTz=ketchupTz(particle,fzvzmustruct,Nz,zcorn,z,dz,Nspecies);
    kTp=ketchupTp(particle,fzvzmustruct,Nz,zcorn,z,dz,Nspecies,B);
    stz=size(kTz);
    for jj=1:stz(1)
      subplot(stz(1),3,(jj-1)*3+1)
      hold on
      plot(z,kTz(jj,:),colours(mod(ii-1,5)+1))
      hold off
      title(['species ' num2str(finitemass(jj))], ...
            'fontname','times','fontsize',18)
      subplot(stz(1),3,(jj-1)*3+2)
      hold on
      plot(z,kTp(jj,:),colours(mod(ii-1,5)+1))
      hold off
      subplot(stz(1),3,(jj-1)*3+3)
      hold on
      plot(z,(kTz(jj,:)+2*kTp(jj,:))/3,colours(mod(ii-1,5)+1))
      hold off
      title(['species ' num2str(finitemass(jj))], ...
            'fontname','times','fontsize',18)
    end
  end
  
  figure(5)
  for jj=1:stz(1)
    subplot(stz(1),3,(jj-1)*3+1)
    set(gca,'fontname','times','fontsize',14)
    ylabel('k_{B}T_{z}/e','fontname','times','fontsize',18)
    xlabel('z','fontname','times','fontsize',18)
    title(['species ' num2str(finitemass(jj))], ...
          'fontname','times','fontsize',18)
    grid on
    aa=axis;axis([min(z) max(z) aa(3:4)]);
    subplot(stz(1),3,(jj-1)*3+2)
    grid on
    aa=axis;axis([min(z) max(z) aa(3:4)]);
    set(gca,'fontname','times','fontsize',14)
    ylabel('k_{B}T_{\perp}/e','fontname','times','fontsize',18)
    xlabel('z','fontname','times','fontsize',18)
    title(['species ' num2str(finitemass(jj))], ...
          'fontname','times','fontsize',18)
    subplot(stz(1),3,(jj-1)*3+3)
    grid on
    aa=axis;axis([min(z) max(z) aa(3:4)]);
    set(gca,'fontname','times','fontsize',14)
    ylabel('k_{B}(T_{z} + 2T_{\perp})/(3e)','fontname','times','fontsize',18)
    xlabel('z','fontname','times','fontsize',18)
    title(['species ' num2str(finitemass(jj))], ...
          'fontname','times','fontsize',18)
  end
end

figure(6)
plot(zcorn,potentialmatrix(end,:),'k','linewidth',2)
set(gca,'fontname','times','fontsize',14)
xlabel('z','fontname','times','fontsize',18)
ylabel('V_{p}','fontname','times','fontsize',18)



% Characteristic frequencies
omegap2=0;
for ii=1:Nspecies
  if ~isnan(particle(ii).mass)
    omegap2 = omegap2 + elc^2 * particle(ii).n0 /(particle(ii).mass*eps0);
  end
end
omegap0 = sqrt(omegap2);
epsilon_r = max((omegap0*const_a*dt)^2,1);
nlast=densitymatrix(:,:,end);
masses = zeros(Nspecies,1);
for ii=1:Nspecies
  masses(ii)=particle(ii).mass;
end
mfin=~isnan(masses);
fp = sqrt(sum(elc^2*nlast(mfin,:)./ ...
              (eps0*epsilon_r*masses(mfin)*ones(1,Nz)),1))/(2*pi);
fpe = fp; % approximately at least
fpi = fpe*sqrt(me/mp);
fge = elc*B.'/(me*2*pi);
flh = fpi./sqrt(1+fpe.^2./fge.^2);
fuh = sqrt(fpe.^2 + fge.^2);

% Debye length
lambdaD=sqrt(eps0*epsilon_r*kTz./nlast(mfin,:)/elc);
if Nspecies>=3
  lambdaDe=1./sqrt(sum(lambdaD([1 3],:).^(-2)));
else
  lambdaDe=lambdaD(1,:);
end
lambdaDtot=1./sqrt(sum(lambdaD.^(-2)));

% $$$ B=B(:).';dB=dB(:).';
% $$$ lambdaDplus=0.5*(dB./B)*eps0*epsilon_r.*kTz./nlast(mfin,:)/elc + ...
% $$$         sqrt(0.25*((dB./B)*eps0*epsilon_r.*kTz./nlast(mfin,:)/elc).^2 ...
% $$$              + eps0*epsilon_r*kTz./nlast(mfin,:)/elc);
% $$$ lambdaDminus=0.5*(dB./B)*eps0*epsilon_r.*kTz./nlast(mfin,:)/elc - ...
% $$$         sqrt(0.25*((dB./B)*eps0*epsilon_r.*kTz./nlast(mfin,:)/elc).^2 ...
% $$$              + eps0*epsilon_r*kTz./nlast(mfin,:)/elc);

if ~exist('drawfft')
  drawfft = logical(1);
end

if drawfft
  figure(7);clf
  set(gcf,'paperpositionmode','auto','renderer','zbuffer')
  Nfft = 256;
  if length(timestepsEfield)<Nfft
    Nfft=2^(floor(log(length(timestepsEfield))/log(2)));
  end
  Powermatrix = [];
  fs=1/(dump_period_fields*dt);
  for ii = 1:Nz
    E=Efieldmatrix(:,ii);
    [Pxx,Fscale] = pwelch(E-mean(E),hanning(Nfft), ...
                     round(0.65*Nfft),Nfft,fs);
    Powermatrix = cat(2,Powermatrix,Pxx);
  end
  pp = Powermatrix;
  ss = size(pp);
  pp = [[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
  maxpp=max(max(pp));
  pz = zcorn;
  pf = [Fscale; 2*Fscale(end)-Fscale(end-1)];
  surf(pz,pf,log10(pp))
  clear pp
  set(gca,'fontname','times','fontsize',14)
  view(2)
  axis([min(pz) max(pz) min(pf) max(pf)])
  shading flat
  hh=colorbar;
  set(hh,'fontname','times','fontsize',14)
  xlabel('z','fontname','times','fontsize',18)
  ylabel('f','fontname','times','fontsize',18)
  grid off
 
  % draw a plasma frequency curve  
  hold on
  plot3(z,fp,ones(size(fp))*log10(maxpp*2),'w')
  plot3(z,fpi,ones(size(fpi))*log10(maxpp*2),'k')
% $$$   plot3(z,flh,ones(size(flh))*log10(maxpp*2),'k')
% $$$   plot3(z,fge,ones(size(fge))*log10(maxpp*2),'g')
% $$$   plot3(z,fuh,ones(size(fuh))*log10(maxpp*2),'r')
  hold off
end
