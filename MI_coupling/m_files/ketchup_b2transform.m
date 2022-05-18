% ketchup_b2transform computes the z-vector and g' given a xi vector as input
%
% [z,gp] = ketchup_b2transform(UHPpole,UHPresid,xi,zmin,zmax)

function [g,gp] = ketchup_b2transform(UHPpole,UHPresid,xi,zmin,zmax)
  
  ALLpole=[UHPpole(:).' UHPpole(:)'];
  ALLresid=[UHPresid(:).' UHPresid(:)'];

% g0 is the value at xi=0 and that serves as an offset, making ghalf start
% exactly at zmin. 
% intg is the integral from xi=0 to xi=1. This is used to normalise the
% transformation so that ghalf ends exactly at zmax.
  g0 = real(sum(ALLresid.'.*log( -ALLpole.')));
  intg = real(sum(ALLresid.'.*log(ones(size(ALLpole)).' - ALLpole.')))-g0;

  gp = sum(ALLresid.'*ones(size(xi))./ ...
           (ones(size(ALLpole)).'*xi - ALLpole.'*ones(size(xi))))* ...
       (zmax-zmin)/intg;
  g  = zmin + (-g0 + sum(ALLresid.'*ones(size(xi)).* ...
                         log(ones(size(ALLpole)).'*(xi) - ...
                             ALLpole.'*ones(size(xi)))))* ...
       (zmax-zmin)/intg;

  gp = real(gp); g=real(g);
