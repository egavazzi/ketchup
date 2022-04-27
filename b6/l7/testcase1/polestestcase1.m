% This is an input file for ketchup_b2transfplot.
% It is a modified version of an old simple pole expansion file.

ALLpole =[];
ALLresid=[];
n_vec    =[3 0 0];                        % # of terms for Maxwell expansion
nb_vec   =[0 6 1];                        % # of terms for Chebyshev 
vtbut_vec=[0.05 0.25 1];                   % Cut velocity in cs units
vdbut_vec=[0.0  0.0 1.0];                    % Chebyshev center velocity
ripple   =[0.1 0.1 0.1];                    % ripple
op2      =[0.03 0.35 0.69];                  % densities
vt       =[0.4 1.0 1.0];                   % thermal velocity cs==1
vd       =[0.0 0.0 0.0];                   % drift velocity cs==1

op2=op2/sum(op2);			%use total plasma frequency

for jj=1:3
  n= n_vec(jj);  
  nb=nb_vec(jj);
  vtbut=vtbut_vec(jj);
  vdbut=vdbut_vec(jj);
  
  pol=0*(1:2*n+1);
  pol(2*n+1)=1;
  for ii=1:n                             % Poles for Maxwellian
    pol(2*(n-ii)+1)=1.0/(2^ii*prod(1:ii));
  end
  sol=roots(pol);
  res=idpole2res(sol);
 
  z0=idchebyshev(nb,ripple(jj),vtbut); 
  z0=z0+vdbut;

  allpole =[sol;z0];                    % poles from expansion and mask
  allresid=idpole2res(allpole);
  [uhppole,uhpresid]=idkeepuhp(allpole,allresid);

  allresid=allresid/(2i*pi*sum(uhpresid)); % normalise zeroth moment
  for kk=1:length(allpole)                 % return to variables with dimension
    ALLpole =[ALLpole allpole(kk)*vt(jj)+vd(jj)];
    ALLresid=[ALLresid op2(jj)*allresid(kk)];
  end
  
end

[UHPpole,UHPresid]=idkeepuhp(ALLpole,ALLresid);
