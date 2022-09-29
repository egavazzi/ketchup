U_final =  3000;
nbre_iteration = 320000; 

% Creating a matrix [nbre_iteration x 2], with the index of iteration in
% the first column, and the corresponding voltage in the second column. The
% voltage is calculated as a gradual increase from 0V to U_final

iteration =  [1:1:nbre_iteration].';
voltage = 3000/nbre_iteration .* iteration;

fileName = fopen('volt.dat','w');
fprintf(fileName,'%f %f\n', [iteration,voltage].');
fclose(fileName);