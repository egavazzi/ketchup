%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [resid]=idpole2res(pole)
n=length(pole);
resid=1:n;
a=zeros(n);
for ii=1:n
  for j=1:n
    a(ii,j)=pole(j)-pole(ii);
  end
  a(ii,ii)=1;
end
resid=1./prod(a);