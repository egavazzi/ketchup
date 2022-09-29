function variable = finding_var(variable_of_interest)
% Find the value of the variable in the file inputb6.m
% Note that the syntax in the inputb6.m MUST be in the form 
%                 "var = somenumbers;"
% with a space before and after the '=' sign.
% 
% Calling:
%   variable = finding_var(variable_of_interest)
% 
% Input:
%   variable_of_interest - must be a char and match name in inputb6
%
% Output:
%   variable - value written in inputb6 of variable_of_interest


filetext = fileread('inputb6.m');

expr = append('[^\n]*',variable_of_interest,' =[^;]*');
matches = regexp(filetext,expr,'match');

expr2 = '[^\s]*\d.*';
matches2 = regexp(matches{1},expr2,'match');

variable = str2num(matches2{1});
end