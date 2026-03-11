
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

    cmd = sprintf('ins -ts=%s -tzScale=1 -tf=%g -ms=%s -knownSolution=%s -idebug=%d -nu=0.1 -bcs=nnnn -N0=%d -plotOption=-1;',...
              par.ts,par.tf,par.ms,par.knownSolution, par.idebug,Nx);
    % cmd = sprintf('ins -ts=ab2 -tzScale=1 -tf=.1 -ms=poly -idebug=1 -nu=0.1 -degreex=2 -degreet=2 -bc1=noSlipWall -bc2=dirichlet  -bc3=noSlipWall -bc4=dirichlet -plotOption=0 -N0=%d;',Nx);
    % fprintf('Running [%s]\n',cmd);
    eval(cmd); 

    % output results from the run are found here:
    Nv(ires) = Nx;
    % fprintf(' ires=%d: maxErr=[%9.2e,%9.2e,%9.2e,%9.2e]\n',ires,ans.maxErr(1),ans.maxErr(2),ans.maxErr(3),ans.maxErr(4));
    maxErr(1:4,ires) = ans.maxErr(1:4);

    cpuv(ires) = ans.cpu;
    % pause

  end 

  fprintf('  N     p-err  ratio   u-err  ratio   v-err ratio    div  ratio    cpu(s)  ratio\n');
  fprintf(' --------------------------------------------------------------------------------- \n');
  for ires=1:par.numResolutions
    if ires==1 
      fprintf('%4d  %8.2e       %8.2e       %8.2e       %8.2e        %8.2e\n',Nv(ires),maxErr(1,ires),maxErr(2,ires),maxErr(3,ires),maxErr(4,ires),cpuv(ires));
    else
      fprintf('%4d  %8.2e %4.1f  %8.2e %4.1f  %8.2e %4.1f  %8.2e %4.1f   %8.2e %4.1f\n',...
          Nv(ires),...
          maxErr(1,ires),maxErr(1,ires-1)/maxErr(1,ires),...
          maxErr(2,ires),maxErr(2,ires-1)/maxErr(2,ires),...
          maxErr(3,ires),maxErr(3,ires-1)/maxErr(3,ires),...
          maxErr(4,ires),maxErr(4,ires-1)/maxErr(4,ires),...
          cpuv(ires),cpuv(ires)/cpuv(ires-1));
    end

  end

return
end

% fprintf('\\bigskip\n');
% fprintf('%% ....................................... (a) .......................................................\n');
% fprintf('\\noindent(a) FE results\n');
% 
% % FE - poly 
% fprintf('\\begin{lstlisting}[frame=single,caption={insPP FE+CD2 poly(2,1)}]\n');
% insPP -ts=fe -tf=.5 -ms=poly -idebug=1 -nu=0.1 -degreex=2 -degreet=1 -bc1=dirichlet -bc2=dirichlet -bc3=dirichlet -bc4=dirichlet -N0=10 -numResolutions=1
% fprintf('\\end{lstlisting}\n');
% 
% fprintf('\\bigskip\n');
% % FE - trig
% fprintf('\\begin{lstlisting}[frame=single,caption={insPP FE+CD2 trig}]\n');
% insPP -ts=fe -tf=.25 -ms=trig -idebug=1 -nu=0.1 -bc1=dirichlet -bc2=dirichlet -bc3=dirichlet -bc4=dirichlet -N0=10 -numResolutions=3
% fprintf('\\end{lstlisting}\n');
% 
% 
% 
% fprintf('\\bigskip\n');
% fprintf('%% ....................................... (b) .......................................................\n');
% fprintf('\\noindent(b)  AB2 results.\n');
% 
% fprintf('\\begin{lstlisting}[frame=single,caption={insPP AB2 + poly(2,2)}]\n');
% insPP -ts=ab2 -tf=.5 -ms=poly -idebug=1 -nu=0.1 -degreex=2 -degreet=2 -bc1=noSlipWall -bc2=dirichlet  -bc3=noSlipWall -bc4=dirichlet -N0=10 -numResolutions=1
% fprintf('\\end{lstlisting}\n');
% 
% fprintf('\\begin{lstlisting}[frame=single,caption={insPP AB2 + trig}]\n');
% insPP -ts=ab2 -tf=.25 -ms=trig -idebug=1 -nu=0.1 -bc1=noSlipWall -bc2=dirichlet -bc3=noSlipWall -bc4=dirichlet -N0=10 -numResolutions=3
% fprintf('\\end{lstlisting}\n');
% 
% 
% fprintf('\\bigskip\n');
% fprintf('%% ....................................... (c) .......................................................\n');
% fprintf('\\noindent(c) IMEX results. \n');
% 
% fprintf('\\begin{lstlisting}[frame=single,caption={insPP IMEX + poly(2,1)}]\n');
% insPP -ts=im2 -tf=.5 -ms=poly -idebug=1 -nu=0.1 -degreex=2 -degreet=1 -bc1=noSlipWall -bc2=dirichlet  -bc3=noSlipWall -bc4=dirichlet -N0=10 -numResolutions=1
% fprintf('\\end{lstlisting}\n');
% 
% 
% fprintf('\\bigskip\n');
% fprintf('\\begin{lstlisting}[frame=single,caption={insPP IMEX + trig}]\n');
% insPP -ts=im2 -tf=.75 -ms=trig -idebug=1 -nu=0.1 -bc1=noSlipWall -bc2=dirichlet -bc3=noSlipWall -bc4=dirichlet -N0=10 -numResolutions=3
% fprintf('\\end{lstlisting}\n');
% 