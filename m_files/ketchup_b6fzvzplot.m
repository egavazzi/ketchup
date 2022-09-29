% ketchup_b6fzvzplot plots distribution functions saved by the ketchup b3
% program. This one integrates over mu to plot f(z,vz).
%
% ketchup_b6fzvzplot(particle,distribution,Nz,zcorn,z,dz,Nspecies)
%
% HG 2010-01-22 (b1 version)
%
% HG 2010-10-26 (b2 version)
%
% HG 2011-01-05 (b3 version)
%
% HG 2013-02-25 (b6 version)

function ketchup_b6fzvzplot(particle,distribution,Nz,zcorn,z,dz,Nspecies)
    
set(gcf,'paperpositionmode','auto','renderer','zbuffer')

c4=[];
for nn=1:4
  c1=jet(16*4^nn);
  c4=[c4;c1((1+(nn-1)*2*4^nn):nn*2*4^nn,:)];
end
c4=[c4;c1(nn*2*4^nn+1:end,:)];

c3=[];
for nn=1:3
  c1=jet(16*4^nn);
  c3=[c3;c1((1+(nn-1)*2*4^nn):nn*2*4^nn,:)];
end
c3=[c3;c1(nn*2*4^nn+1:end,:)];

c2=[];
for nn=1:2
  c1=jet(16*4^nn);
  c2=[c2;c1((1+(nn-1)*2*4^nn):nn*2*4^nn,:)];
end
c2=[c2;c1(nn*2*4^nn+1:end,:)];

colormap(c3)



finitemass=[];
for ii=1:Nspecies
  if ~isnan(particle(ii).mass) & ~isinf(particle(ii).mass)
    finitemass=[finitemass ii];
  end
end

for ii=1:length(finitemass)
  if ndims(distribution(finitemass(ii)).f)>2
    f=[];
    for jj=1:Nz
      f(:,jj)=distribution(finitemass(ii)).f(:,:,jj) * ...
              particle(finitemass(ii)).dmu.';
    end
  else
    f=distribution(finitemass(ii)).f;
  end
  subplot(length(finitemass),1,ii)
  blockends=[0 find(diff(distribution(finitemass(ii)).ivzoffset)~=0) Nz];
  hold on
  for kk=2:length(blockends)
    pp=f(:,blockends(kk-1)+1:blockends(kk));
    ss=size(pp);
    pp=[[pp zeros(ss(1),1)];zeros(1,ss(2)+1)];
    zplot=zcorn(blockends(kk-1)+1:blockends(kk)+1);
    vzplot=particle(finitemass(ii)).vzcorn + ...
           distribution(finitemass(ii)).ivzoffset(blockends(kk)) * ...
           particle(finitemass(ii)).dvz;
    surf(zplot,vzplot,pp)
%    surf(zplot,vzplot,log10(pp))
  end

  view(2);grid off;shading flat
  axis([min(zcorn) max(zcorn) ...
        (min(particle(finitemass(ii)).vzcorn) + ...
         min(distribution(finitemass(ii)).ivzoffset) * ...
         particle(finitemass(ii)).dvz) ...
        (max(particle(finitemass(ii)).vzcorn) + ...
         max(distribution(finitemass(ii)).ivzoffset) * ...
         particle(finitemass(ii)).dvz)])

  colorbar
  set(gca,'fontname','times','fontsize',14)
  xlabel('z','fontname','times','fontsize',18)
  ylabel('v_{z}','fontname','times','fontsize',18)
  % Decreased title font size.  HG 2015-05-08
  title(['species ' num2str(finitemass(ii))], ...
        'fontname','times','fontsize',16)
  hold off

end
