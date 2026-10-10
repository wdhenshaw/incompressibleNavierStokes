% ------------------------------------------------
% Look for command line options of the form
%     -name=value
% and set
%     par.name = value
% If par has a member called name (or if addNamesToPar==1).
% echo          : optional argument, echo=0 means do not echo when a variable is set
% addNamesToPar : optional argument,
%                 0 (default) = only set par.name if par already has a member called name,
%                 1 = add name to par if it is not already a member
%
%  Examples:
%    par.nu=1e-2
%    par.N=10;
%  To set an array value:
%    par.myArray(2)=5
%
function par = assignCommandLineOption( line, par, echo, addNamesToPar )

  if( nargin<3 ) echo=1; end
  if( nargin<4 ) addNamesToPar=0; end
  % look for a substring
  %  -name= or -name =
  n1 = strfind(line,'-');
  n2 = strfind(line,'=');

  if( length(n1)>0 && length(n2)>0 )
  	name = line(n1(1)+1:n2(1)-1);
  	value = line(n2(1)+1:end);

		if ~isnan(str2double(value))
		  % fprintf('The string value=[%s] is a valid number.\n',value);
		else
		  % fprintf('The string value=[%s] is not a valid number.\n',value);
		  value = sprintf("'%s'",value);
		end

    % check for an array assigment:
    m1=strfind(name,'(');
    m2=strfind(name,')');

    % fprintf('line=[%s] name=[%s], value=[%s], array? [m1,m2]=[(,)] = [%d,%d]\n',line,name,value,m1,m2)
    isArray=0;
    if( length(m1)>0 && length(m2)>0 )
      isArray=1;
      index = int32(str2double(name(m1(1)+1:m2(1)-1)));
      name = name(1:m1(1)-1);
      % fprintf('Array asssigment found... name=[%s] index=%d\n',name,index);
    end

    if( isfield(par,name) || addNamesToPar==1 )
     	% fprintf('par.%s exists\n',name);
      % lambdaMax=par.lambdaMax;
      if( isArray==0 )
        cmd = sprintf('par.%s=%s;',name,value);
      else
        cmd = sprintf('par.%s(%d)=%s;',name,index,value);
      end
      % fprintf('cmd=[%s]\n',cmd);
      eval(cmd);
      if( echo==1 && strcmp(name,'debug')==0 )
        if( isArray==0 )
         fprintf('AssignCommandLineOption: Setting: par.%s=%s\n',name,value);
        else
         fprintf('AssignCommandLineOption: Setting: par.%s(%d)=%s\n',name,index,value);
        end
      end
   end

   % pause
  end

return
end
