% ketchup_b2transfsave saves UHPpole and UHPresid in a file that can be
% read by ketchup.
% 
% ketchup_b2transfsave(UHPpole,UHPresid,outfile)


function ketchup_b2transfsave(UHPpole,UHPresid,outfile)

  if nargin<3
    outfile='transfb2.dat'
  end
  
  aa=[real(UHPpole(:)) imag(UHPpole(:)) real(UHPresid(:)) imag(UHPresid(:))];
  
  eval(['save -ascii -double ' outfile ' aa'])

