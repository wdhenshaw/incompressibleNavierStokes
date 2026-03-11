% run regression tests in the check folder
% 
% runChecks -replace=[0|1]
% 
% runChecks              : run all cases and compared to saved check files
% runChecks -replace=1   : run all cases and save new check files (do this if you add a new case)
%
function runChecks( varargin )

  clearvars -except varargin; 

  addpath(genpath(pwd)); % allow matlab to find files in subfolders


  % par.replace=0;  % set to 1 to save new check files
  checkDir = './check';
  % par.echo=0; 

  % --- read command line args ---
  for i = 1 : nargin
    line = varargin{i};
    % par = assignCommandLineOption( line, par, par.echo );
  end
  replace=0; % par.replace;

  eval(sprintf('check -replace=%d',replace));
  % eval('check/check');
  % check

return
end