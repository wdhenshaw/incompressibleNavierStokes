
%
%  Run a grid convergence test 
% 
% Examples:
%  runGridConvergence -ts=ab2 -tf=.25 -ms=trig -numResolutions=3
%  runGridConvergence -ts=im2 -tf=1.0 -ms=trig -numResolutions=3
%
%  runGridConvergence -ts=im2 -tf=1.0 -ms=none -knownSolution=TaylorGreen -numResolutions=3
% 
%  runGridConvergence -ts=ab2 -ms=poly -numResolutions=2
%
function runGridConvergence( varargin )


 clearvars -except varargin; 


 par.numResolutions=3;
 par.N0=10;
 par.ts = 'ab2'; 
 par.ms='trig'; % 'poly' 
 par.knownSolution='none';
 par.tf=0.25;
 par.idebug=0; 
 par.map = 'Cartesian';           % 'Cartesian', 'Rectangle', 'Annulus', 'TFI', ...
 par.bcs='dddd'; 
 par.motion='none';
 par.dtMax=1e8; 

 par.gravity = 0;
 par.gamma = 0;
 

 par.echo = 0;

 % --- read command line args ---
  for i = 1 : nargin
    line = varargin{i};
    par = assignCommandLineOption( line, par, par.echo );

  end

  numComp = 4; % [p,u,v,div]
  maxErr = zeros(numComp,par.numResolutions);

  for ires=1:par.numResolutions

    Nx= par.N0*2^(ires-1); 


    cmd = sprintf('ins -ts=%s -tzScale=1 -tf=%g -ms=%s -knownSolution=%s -idebug=%d -nu=0.1 -bcs=%s -N0=%d -map=%s -motion=%s -dtMax=%g -plotOption=-1 -gravity=%g -gamma=%g -computeErrors=1 -ampfs=1e-4 -icfs=cos -ya=-1 -yb=0 ;',...
              par.ts,par.tf,par.ms,par.knownSolution, par.idebug,par.bcs, Nx, par.map, par.motion,par.dtMax, par.gravity, par.gamma);

	% Add necessary flags for some known solutions
	% if( strcmp(par.knownSolution, 'GravityCapillaryWave') )
	% 	cmd = sprintf('%s\b ', cmd);
	% end

    % cmd = sprintf('ins -ts=ab2 -tzScale=1 -tf=.1 -ms=poly -idebug=1 -nu=0.1 -degreex=2 -degreet=2 -bc1=noSlipWall -bc2=dirichlet  -bc3=noSlipWall -bc4=dirichlet -plotOption=0 -N0=%d;',Nx);
    fprintf('Running [%s]\n',cmd);
    eval(cmd); 

    % output results from the run are found here:
    Nv(ires) = Nx;
    dtv(ires) = ans.par.dt;
    fprintf(' ires=%d: Nx=%d maxErr=[%9.2e,%9.2e,%9.2e,%9.2e]\n',ires,Nx,ans.maxErr(1),ans.maxErr(2),ans.maxErr(3),ans.maxErr(4));
    maxErr(1:4,ires) = ans.maxErr(1:4);

    cpuv(ires) = ans.cpu;
    % pause

    par.dtMax= par.dtMax/2; % reduce dtMax for the implicit method
 
  end 

  fprintf('  N     dt      p-err  ratio   u-err  ratio   v-err  ratio    div   ratio    cpu(s)  ratio\n');
  fprintf(' -------------------------------------------------------------------------------------- \n');
  for ires=1:par.numResolutions
    if ires==1 
      fprintf('%4d  %8.2e %8.2e       %8.2e       %8.2e       %8.2e        %8.2e\n',Nv(ires),dtv(ires),maxErr(1,ires),maxErr(2,ires),maxErr(3,ires),maxErr(4,ires),cpuv(ires));
    else
      fprintf('%4d  %8.2e %8.2e %4.1f  %8.2e %4.1f  %8.2e %4.1f  %8.2e %4.1f   %8.2e %4.1f\n',...
          Nv(ires),dtv(ires),...
          maxErr(1,ires),maxErr(1,ires-1)/maxErr(1,ires),...
          maxErr(2,ires),maxErr(2,ires-1)/maxErr(2,ires),...
          maxErr(3,ires),maxErr(3,ires-1)/maxErr(3,ires),...
          maxErr(4,ires),maxErr(4,ires-1)/maxErr(4,ires),...
          cpuv(ires),cpuv(ires)/cpuv(ires-1));
    end

  end

  %
  % Output a LaTeX table of results
  %
  insPar = ans.par; 

  extra='';
  if( ~strcmp(insPar.ms,'none') )     extra = strcat(extra,sprintf('MS%s'),insPar.ms);          end 
  if( ~strcmp(insPar.motion,'none') ) extra = strcat(extra,sprintf('Motion%s'),insPar.motion);  end 

  name=sprintf('insTS%sBC%sMap%s%s',insPar.ts,insPar.bcLabel,insPar.map,extra);

  tableDir = 'doc/tables';

  latexFileName = sprintf('%s/%s.tex',tableDir,name);

  output = fopen(latexFileName,'w');


  %  Output results as a LaTeX table
  fprintf(output,"\n%% ------------ INS Table for LaTeX from runGridConvergence.m  -------------------------------\n");

  fprintf(output,"\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|} \\hline\n"); 
  fprintf(output," \\multicolumn{10}{|c|}{INS: ts=%s, BC=%s, Map=%s, Motion=%s, MS=%s } \\\\ \\hline \n",insPar.ts,insPar.bcLabel,insPar.map,insPar.motion,insPar.ms);
  fprintf(output,"     N  &   $\\Delta t$ & $E_p$   &  r  &   $E_u$   &   r    &   $E_v$   & r   &  $|\\grad\\cdot\\uv|$  & r   \\\\ \\hline \n");
  
  ires=1;
  fprintf(output," %3d   &  %8.2e  &  %9.2e  &        &  %9.2e  &        &  %9.2e  &        &  %9.2e  &        \\\\ \n",...
          Nv(ires),dtv(ires),...
          maxErr(1,ires),...
          maxErr(2,ires),...
          maxErr(3,ires),...
          maxErr(4,ires) ... 
           );  
  for ires=2:par.numResolutions
    fprintf(output," %3d   &  %8.2e  &  %9.2e  & %3.1f  &  %9.2e  & %3.1f  &  %9.2e  & %3.1f  &  %9.2e  & %3.1f  \\\\ \n",...
          Nv(ires),dtv(ires),...
          maxErr(1,ires),maxErr(1,ires-1)/maxErr(1,ires),...
          maxErr(2,ires),maxErr(2,ires-1)/maxErr(2,ires),...
          maxErr(3,ires),maxErr(3,ires-1)/maxErr(3,ires),...
          maxErr(4,ires),maxErr(4,ires-1)/maxErr(4,ires) ...
           );
  end 

  fprintf(output," \\hline \n");
  fprintf(output,"\\end{tabular}\n");  


  fclose(output);
  if( 1==1 || par.idebug >0 ) fprintf('Wrote table of ins results to file=[%s]\n',latexFileName); end

  return

end

