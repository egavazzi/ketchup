% ketchup_b6fmuconv converts output files to large .mat files
%
% The output files treated by this program are the full particle
% distribution files. The program reads those and stores one .mat file
% per species. Thus this program need only be used for very large
% simulation runs. For smaller files, use ketchup_b6conv as usual.
%
% HG 2013-08-20

ccc=pwd;
[a0,a1]=system('uname -n'); [a0,a2]=system('ps |grep -i matlab');
dummy=[a1 a2 datestr(now)];
dummy=dummy(double(dummy)~=10); dummy=dummy(double(dummy)~=32);

if ~exist('absolutelynomessages')
  absolutelynomessages=logical(0);
end

% ------- input data ------- %
% The easiest way to get the general parameters
inputb6;

% If the fields_per_file parameter wasn't defined in the input file,
% default to 1 field per file.
if ~exist('fields_per_file')
  fields_per_file = 1;
end

% Now read the species specific data
particle=struct();
fid=fopen('inputb6.m');
for ii=1:Nspecies
  theline=fgetl(fid);
  if length(theline)<5,
    theline=[theline '     '];
  end
  while ~(strcmp(theline(1:5),'%SPEC') | strcmp(theline(1:5),'%spec'))
    theline=fgetl(fid);
    if length(theline)<5,
      theline=[theline '     '];
    end
  end
  while ~(strcmp(theline(1:4),'%END') | strcmp(theline(1:4),'%end'))
    eval(theline)
    theline=fgetl(fid);
    if length(theline)<4,
      theline=[theline '    '];
    end
  end
  particle(ii).Nvz=Nvz;
  particle(ii).vzmin=vzmin;
  particle(ii).vzmax=vzmax;
  particle(ii).Nmu=Nmu;
  particle(ii).mumin=mumin;
  particle(ii).mumax=mumax;
  particle(ii).muexp=muexp;
  particle(ii).mass=mass;
  particle(ii).charge=charge;
  particle(ii).n0=n0;
  particle(ii).vz0=vz0;
  particle(ii).kTz=kTz;
  particle(ii).kTp=kTp;
  particle(ii).n0L=n0L;
  particle(ii).vz0L=vz0L;
  particle(ii).kTzL=kTzL;
  particle(ii).kTpL=kTpL;
  particle(ii).n0R=n0R;
  particle(ii).vz0R=vz0R;
  particle(ii).kTzR=kTzR;
  particle(ii).kTpR=kTpR;
end
fclose(fid);

% Find the number of processors used from the job.bat file, if there is one. 
if exist([pwd '/job.bat'])
  fid=fopen('job.bat');
  theline=fgetl(fid);
  while isempty(strfind(theline,'Nprocs')) & ...
        (isempty(strfind(theline,'PBS')) | ...
         isempty(strfind(theline,'-l')) | ...
         isempty(strfind(theline,'nodes')) )
    theline=fgetl(fid);
    if ~ischar(theline), break, end
  end
  fclose(fid);
  if ~isempty(strfind(theline,'Nprocs'))
    eval([theline(strfind(theline,'Nprocs'):end) ';'])
    Nprocsref=Nprocs;
    clear Nprocs
  elseif ~isempty(strfind(theline,'nodes'))
    eval([theline(strfind(theline,'nodes'):end) ';'])
    Nprocs=nodes;
    Nprocsref=Nprocs;
    clear Nprocs
  else
    Nprocsref=NaN;
  end
else
  Nprocsref=NaN;
end

% construct xi-vectors
dxi=1/Nz;
xicorn=dxi*[0:Nz];
xi=0.5*(xicorn(1:end-1) + xicorn(2:end));
% compute transformation, i.e. z-vectors and g'
aa=load(transffilename);
UHPpole  = aa(:,1) + 1i*aa(:,2);
UHPresid = aa(:,3) + 1i*aa(:,4);
[z,gp] = ketchup_b2transform(UHPpole,UHPresid,xi,zmin,zmax);
zcorn = ketchup_b2transform(UHPpole,UHPresid,xicorn,zmin,zmax);
dz    = diff(zcorn);

for ii=1:Nspecies
  particle(ii).dvz=(particle(ii).vzmax-particle(ii).vzmin)/particle(ii).Nvz;
  particle(ii).vzcorn=particle(ii).vz0 + ...
      particle(ii).vzmin + particle(ii).dvz*[0:particle(ii).Nvz];
  particle(ii).vz = ...
      0.5*(particle(ii).vzcorn(1:end-1)+particle(ii).vzcorn(2:end));

  vmu=[1:particle(ii).Nmu];
  particle(ii).mu = particle(ii).mumin + ...
      0.5*((vmu.^particle(ii).muexp+(vmu-1).^particle(ii).muexp) / ...
           particle(ii).Nmu^particle(ii).muexp) * ...
      (particle(ii).mumax-particle(ii).mumin);
  particle(ii).dmu = ((vmu.^particle(ii).muexp - ...
                       (vmu-1).^particle(ii).muexp) / ...
                      particle(ii).Nmu^particle(ii).muexp) * ...
      (particle(ii).mumax-particle(ii).mumin);
  particle(ii).mucorn = particle(ii).mu-0.5*particle(ii).dmu; 
  particle(ii).mucorn = [particle(ii).mucorn particle(ii).mu(end) + ...
                      0.5*particle(ii).dmu(end)];
end

cd outp

% --- distribution function f(z,vz,mu) --- %
% one file per timestep and species. For compatibility with other
% programs the distribution function is put in a structured array of all
% species, but only the species in the variable thisspecies contains a
% three-dimensional array for f(z,vz,mu).

% Prevent two processes from performing simultaneous conversions
if exist([pwd '/datfiles/lock.fzvzmu'])
  if ~absolutelynomessages
    disp('Another process is already working on this directory.')
    disp('If this is not the case, remove the file')
    disp([pwd '/datfiles/lock.fzvzmu'])
  end
else
  all_is_fine=logical(0);
  try
    dlmwrite('datfiles/lock.fzvzmu',dummy,'')
    fid=fopen('datfiles/lock.fzvzmu','r');
    lock=textscan(fid,'%s');
    fclose(fid);
    if strcmp(dummy,lock{1})
      all_is_fine=logical(1);
    end
  catch
    all_is_fine=logical(0);
  end

  if all_is_fine
    try
      dd=dir('datfiles/fzvzmu');
      fzvzmustruct=struct();
      for ii=1:Nspecies
        fzvzmustruct(ii).timestep=0;
      end

      % pick out the files containing species 1 and process 0. 
      % Then start the loop.
      for speciesnumber=Nspecies:-1:1
        if ~isnan(particle(speciesnumber).mass) & ...
              ~isinf(particle(speciesnumber).mass)
          break
        end
      end
      startfiles=[];Nprocs=0;
      for ii=3:length(dd)
        if length(dd(ii).name)>=20
          if strcmp(dd(ii).name(1:6),'fzvzmu')
            if strcmp(dd(ii).name(15:16),num2str(speciesnumber,'%0.2d')) & ...
                  strcmp(dd(ii).name(22:end),'.ketchup.dat')
              procid = str2num(dd(ii).name(18:21));
              Nprocs = max(Nprocs,procid);
              if strcmp(dd(ii).name(17:21),'p0000')
                startfiles = [startfiles ii];
              end
            end
          end
        end
      end
      Nprocs = Nprocs + 1; % Numbers from 0 to Nprocs-1

      % To prevent attempts to process files that are being written we wait
      % twenty seconds if there are less than two time steps to process.
      if length(startfiles)<2 & length(startfiles)>0 & ~(Nprocs<Nprocsref)
        pause(20)
      end

      if length(startfiles)>0 & ~(Nprocs<Nprocsref)
        for ii = startfiles
          % If not all files of the largest non-infinite mass species
          % number exist, that is an error
          infilesexistnot = logical(zeros(1,Nprocs));
          for jj = 0:Nprocs-1
            infilesexistnot(jj+1) = ...
                ~exist(['datfiles/fzvzmu/' dd(ii).name(1:14) ...
                        num2str(speciesnumber,'%0.2d') ...
                        'p' num2str(jj,'%0.4d') '.ketchup.dat']);
          end        
          if sum(infilesexistnot)>0
            error(['fzvzmu: infilesexistnot=',num2str(infilesexistnot)])
          end
        
          thistimestep=str2num(dd(ii).name(7:13));
          for thisspecies = 1:Nspecies
            if ~isnan(particle(thisspecies).mass) & ...
                  ~isinf(particle(thisspecies).mass)
              fzvzmustruct(thisspecies).timestep=thistimestep;
              fzvzmustruct(thisspecies).f = ...
                  zeros(particle(thisspecies).Nvz,particle(thisspecies).Nmu,Nz);
              ivzoffset=[];fcounter=1;

              for jj = 0:Nprocs-1
                infile = ['datfiles/fzvzmu/' dd(ii).name(1:14) ...
                          num2str(thisspecies,'%0.2d') ...
                          'p' num2str(jj,'%0.4d') '.ketchup.dat'];
                fid=fopen(infile,'r');
                if fid<0
                  error(['Error reading file ' infile])
                end

                for kk = 1:Nz
                  instruct = textscan(fid,'%*s%*s%f',1);
                  if isempty(instruct{1}), break, end
                  ivzoffset=[ivzoffset instruct{1}];
                  instruct = textscan(fid,'%f');
                  fzvzmustruct(thisspecies).f(:,:,fcounter) = ...
                      reshape(instruct{1}, [particle(thisspecies).Nmu ...
                                      particle(thisspecies).Nvz]).';
                  fcounter = fcounter + 1;
                end
      
                fclose(fid);
              end
              fzvzmustruct(thisspecies).ivzoffset = ivzoffset;
              % save this species
              outfile = ['fzvzmu' num2str(thistimestep,'%0.7i') ...
                         's' num2str(thisspecies,'%0.2i') '.mat'];
              save(outfile,'-v7.3', ...
                   'fzvzmustruct','thistimestep','thisspecies', ...
                   'particle','Nz','dz','zcorn','z','dt','Niter', ...
                   'dump_period_fields','dump_period_distr','dump_start', ...
                   'shift_test_period','zmin','zmax','Nspecies','voltage')
              % reset this species in the structure
              fzvzmustruct(thisspecies).ivzoffset=[];
              fzvzmustruct(thisspecies).f=[];
              fzvzmustruct(thisspecies).timestep=0;
            end
          end

          delete([infile(1:30) '*p*.ketchup.dat']);

          if ~absolutelynomessages
            disp(['fzvzmu: timestep ' num2str(thistimestep) ' done!'])
          end
        end
      end
      clear fzvzmustruct
      delete('datfiles/lock.fzvzmu')
    catch
      felfelfel=lasterror;
      disp(felfelfel.message)
      delete('datfiles/lock.fzvzmu')
    end % end try
  end
end % end if exist('datfiles/lock.fzvzmu')



cd(ccc)
