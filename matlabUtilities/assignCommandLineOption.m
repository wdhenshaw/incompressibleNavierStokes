% ------------------------------------------------
% Look for command line options of the form
%     -name=value
% and set 
%     par.name = value 
% IF par has a member called name.
% echo : optional argument, echo=0 means do not echo when a variable is set
%
function par = assignCommandLineOption( line, par, echo )

	if( nargin<3 ) echo=1; end
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
		% fprintf('line=[%s] name=[%s], value=[%s]\n',line,name,value)

		if( ~isfield(par,name) && echo )
			warning(sprintf('Parameter "%s" not an existing parameter!', name));
		end
		% fprintf('par.%s exists\n',name);
		% lambdaMax=par.lambdaMax;
		cmd = sprintf('par.%s=%s;',name,value);
		% fprintf('cmd=[%s]\n',cmd);
		eval(cmd);
		if( echo==1 && strcmp(name,'debug')==0 ) fprintf('AssignCommandLineOption: Setting: par.%s=%s\n',name,value); end

		% pause
	end

	return
end
