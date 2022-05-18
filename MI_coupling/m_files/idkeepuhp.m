%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% keep the residues and poles in the Upper Half Plane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[pole2 resid2]=idkeepuhp(pole1,resid1)
function [pole2,resid2]=idkeepuhp(pole1,resid1)

n=length(pole1);
pole2  = pole1(imag(pole1)>0);
resid2 = resid1(imag(pole1)>0);

if length(pole2) < n/2
  pole2(n/2)=1e30
end
