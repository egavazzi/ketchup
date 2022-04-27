% ketchup_b6fplot plots distribution functions saved by the ketchup b3
% program. 
%
% ketchup_b6fplot(particle,fzvzmustruct,Nz,zcorn,z,dz,Nspecies)
%
% HG 2013-02-25

function ketchup_b6fplot(particle,fzvzmustruct,Nz,zcorn,z,dz,Nspecies)
    
set(gcf,'paperpositionmode','auto','renderer','zbuffer')

c1=jet(64);c2=jet(256);c3=jet(1024);
c=[c1(1:8,:);c2(32:64,:);c3(257:1024,:)];
colormap(c)

zzlist = [1 round([1 2 3 4]*Nz/4)];
finitemass=[];
for ii=1:Nspecies
  if ~isnan(particle(ii).mass) & ~isinf(particle(ii).mass)
    finitemass=[finitemass ii];
  end
end

% $$$ for jj=1:3
for jj=1:5
  for ii=1:length(finitemass)
    subplot(length(finitemass),5,(ii-1)*5+jj)
    zz=zzlist(jj);
    pp=fzvzmustruct(finitemass(ii)).f(:,:,zz);
    if sum(sum(pp))<1e-301           % Avoid Matlab error
      pp=ones(size(pp))*NaN;         % when displaying
    end                              % small numbers.
    ss=size(pp);
    pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
    vzplot=particle(finitemass(ii)).vzcorn + ...
           fzvzmustruct(finitemass(ii)).ivzoffset(zz) * ...
           particle(finitemass(ii)).dvz;
    muplot=particle(finitemass(ii)).mucorn;
    surf(muplot,vzplot,pp)
    view(2);grid off;shading flat
    axis([min(particle(finitemass(ii)).mucorn) ...
          max(particle(finitemass(ii)).mucorn) ...
          (min(particle(finitemass(ii)).vzcorn) + ...
           fzvzmustruct(finitemass(ii)).ivzoffset(zz) * ...
           particle(finitemass(ii)).dvz) ...
          (max(particle(finitemass(ii)).vzcorn) + ...
           fzvzmustruct(finitemass(ii)).ivzoffset(zz) * ...
           particle(finitemass(ii)).dvz)])
    set(gca,'fontname','times','fontsize',14)
    if ii==1
      title(['z=' num2str(z(zz)) 'm'],'fontname','times','fontsize',18)
    end
    
    if ii==length(finitemass)
      xlabel('\mu   [Am^{2}]','fontname','times','fontsize',18)
    end
    if jj==1
      ylabel('v_{z}  [m/s]','fontname','times','fontsize',18)
    end
  end
end



