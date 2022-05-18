% ketchup_b2transfplot plots coordinate transformation functions for
% ketchup. 
%
% ketchup_b2transfplot(UHPpole,UHPresid,Nz,zmin,zmax) uses the upper half 
%         plane poles supplied in UHPpole and the corresponding residues
%         in UHPresid. 
%
% ketchup_b2transfplot('filename.m',Nz,zmin,zmax) evaluates the commands
%         in filename.m and if that produces usable UHPpole and UHPresid
%         variables it continues to plot them. 
%
% ketchup_b2transfplot('filename.xxx',Nz,zmin,zmax) where .xxx is any
%         extension other than .m will cause ketchup_b2transfplot to try
%         to read UHPpole and UHPresid from that file.

function ketchup_b2transfplot(varargin)

if isstr(varargin{1})
  filename=varargin{1};
  if strcmp(filename(end),'m')
    eval(filename(1:end-2))
  else
    aa=load(filename);
    if isstruct(aa)
      UHPpole=aa.UHPpole;
      UHPresid=aa.UHPresid;
    else
      UHPpole  = aa(:,1) + 1i*aa(:,2);
      UHPresid = aa(:,3) + 1i*aa(:,4);
    end
  end
  Nz = varargin{2};
  zmin = varargin{3};
  zmax = varargin{4};
else
  UHPpole = varargin{1};
  UHPresid = varargin{2};
  Nz = varargin{3};
  zmin = varargin{4};
  zmax = varargin{5};
end


figure(1)
plot(UHPpole,'k*')
set(gca,'fontname','times','fontsize',14)
xlabel('\Re\{b\}','fontname','times','fontsize',18)
ylabel('\Im\{b\}','fontname','times','fontsize',18)


% xi represents cell centres, and xihalf cell borders. g and ghalf and z
% and zhalf are defined in the same way.
dxi = 1/Nz;
xi  = [0:Nz-1]*dxi + 0.5*dxi;
xihalf = dxi*[0:Nz];

[g,gp] = ketchup_b2transform(UHPpole,UHPresid,xi,zmin,zmax);
ghalf = ketchup_b2transform(UHPpole,UHPresid,xihalf,zmin,zmax);

dg=diff(ghalf);


figure(2)
subplot(2,1,1)
plot(xi,g,'k')
set(gca,'fontname','times','fontsize',14)
xlabel('\xi','fontname','times','fontsize',18)
ylabel('z=g(\xi)','fontname','times','fontsize',18)
subplot(2,1,2)
plot(xi,gp,'k')
set(gca,'fontname','times','fontsize',14)
xlabel('\xi','fontname','times','fontsize',18)
ylabel('g''(\xi)','fontname','times','fontsize',18)

figure(3)
subplot(2,1,1)
semilogy(xi,dg,'k')
set(gca,'fontname','times','fontsize',14)
xlabel('\xi','fontname','times','fontsize',18)
ylabel('\Delta z','fontname','times','fontsize',18)
subplot(2,1,2)
plot(g,dg,'k')
set(gca,'fontname','times','fontsize',14)
xlabel('z','fontname','times','fontsize',18)
ylabel('\Delta z','fontname','times','fontsize',18)
