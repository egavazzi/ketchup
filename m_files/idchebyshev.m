function root=idchebyshev(n,r,vc)
polynom=zeros(2*n+1,1);
polynom(2*n+1)=1;
t0=[1];t1=[1 0];tii=[2 0 -1];
if n<1,
  tii=t0;
elseif n<2
  tii=t1;
else
  for ii=2:n
    tii=2*[t1 0]-[0 0 t0];
    t0=t1;t1=tii;
  end
end
polynom=polynom + r^2*conv(tii,tii).';
root=roots(polynom)*vc;
%figure(1);grid;plot(real(root),imag(root),'o');
%v=-3*vc:0.01:3*vc;
%figure(2);grid;plot(v,1./polyval(polynom,v/vc));
