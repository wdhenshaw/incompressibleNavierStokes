% 
% regression tests for ins.m
%
% check -replace=[0|1]
% 
% check              : run all cases and compared to saved check files
% check -replace=1   : run all cases and save new check files (do this if you add a new case)
%
function check( varargin )


 clearvars -except varargin; 


 par.replace=0;  % set to 1 to save new check files
 checkDir = './check'; 
 par.echo=0;

 % --- read command line args ---
  for i = 1 : nargin
    line = varargin{i};
     par = assignCommandLineOption( line, par, par.echo );
  end


  numChecks=0;

  % NOTE: DO NOT PUT A SEMI-COLON AT THE END OF THE CMD 
  numChecks=numChecks+1;
  checkFileName{numChecks}='insAB2PolyTZ.check';
  cmd{numChecks}='ins -ts=ab2 -tzScale=1 -tf=.1 -ms=poly -knownSolution=none -nu=0.1 -degreex=2 -degreet=2 -bcs=nnnn -N0=10 -idebug=0 -plotOption=-1';


  numChecks=numChecks+1;
  checkFileName{numChecks}='insPC2PolyTZ.check';
  cmd{numChecks}='ins -ts=pc2 -tzScale=1 -tf=.1 -ms=poly -knownSolution=none -nu=0.1 -degreex=2 -degreet=2 -bcs=nnnn -N0=10 -idebug=0 -plotOption=-1';

  % noslip and slip walls:
  numChecks=numChecks+1;
  checkFileName{numChecks}='insPC2PolyTZnsns.check';
  cmd{numChecks}='ins -ts=pc2 -tzScale=1 -tf=.1 -ms=poly -knownSolution=none -nu=0.1 -degreex=2 -degreet=2 -bcs=nsns -N0=10 -idebug=0 -plotOption=-1';

  % inflow and outflow 
  numChecks=numChecks+1;
  checkFileName{numChecks}='insPC2PolyTZionn.check';
  cmd{numChecks}='ins -ts=pc2 -tzScale=1 -tf=.1 -ms=poly -knownSolution=none -nu=0.1 -degreex=2 -degreet=2 -bcs=ionn -N0=10 -idebug=0 -plotOption=-1';


  numChecks=numChecks+1;
  checkFileName{numChecks}='insIM2PolyTZ.check';
  cmd{numChecks}='ins -ts=im2 -tzScale=1 -tf=.1 -ms=poly -knownSolution=none -nu=0.1 -degreex=2 -degreet=2 -bcs=nnnn -N0=10 -idebug=0 -plotOption=-1';

  % noslip and slip walls:
  numChecks=numChecks+1;
  checkFileName{numChecks}='insIM2PolyTZsnsn.check';
  cmd{numChecks}='ins -ts=im2 -tzScale=1 -tf=.1 -ms=poly -knownSolution=none -nu=0.1 -degreex=2 -degreet=2 -bcs=snsn -N0=10 -idebug=0 -plotOption=-1';


  numChecks=numChecks+1;
  checkFileName{numChecks}='insPC2TrigTZ.check';
  cmd{numChecks}='ins -ts=pc2 -tzScale=1 -tf=.1 -ms=trig -knownSolution=none -nu=0.1 -bcs=nnnn -N0=10 -idebug=0 -plotOption=-1';

  numChecks=numChecks+1;
  checkFileName{numChecks}='insIM2TrigTZ.check';
  cmd{numChecks}='ins -ts=im2 -tzScale=1 -tf=.1 -ms=trig -knownSolution=none -nu=0.1 -bcs=nnnn -N0=10 -idebug=0 -plotOption=-1';


  numChecks=numChecks+1;
  checkFileName{numChecks}='insIM2TaylorGreen.check';
  cmd{numChecks}='ins -ts=im2 -tf=.2 -ms=none -knownSolution=TaylorGreen -nu=0.1  -bcs=nnnn -N0=20 -idebug=0 -plotOption=-1';


  % AB2 -- map=TFI 
  numChecks=numChecks+1;
  checkFileName{numChecks}='insAB2TrigTFI.check';
  cmd{numChecks}='ins -ts=ab2 -tf=.2 -ms=trig -knownSolution=none -nu=0.1  -bcs=nnnn -N0=20 -idebug=0 -map=TFI -plotOption=-1';

  % IM2 -- map=TFI 
  numChecks=numChecks+1;
  checkFileName{numChecks}='insIM2TrigTFI.check';
  cmd{numChecks}='ins -ts=im2 -tf=.2 -ms=trig -knownSolution=none -nu=0.1  -bcs=nnnn -N0=20 -idebug=0 -map=TFI -plotOption=-1';


  % IM2 ROTATE
  numChecks=numChecks+1;
  checkFileName{numChecks}='insIM2TrigRotate.check';
  cmd{numChecks}='ins -ts=im2 -tf=.2 -ms=trig -knownSolution=none -nu=0.1  -bcs=nnnn -N0=10 -idebug=0 -motion=rotate -plotOption=-1';  

  % PC2 ROTATE
  numChecks=numChecks+1;
  checkFileName{numChecks}='insPC2TrigRotate.check';
  cmd{numChecks}='ins -ts=pc2 -tf=.2 -ms=trig -knownSolution=none -nu=0.1  -bcs=nnnn -N0=10 -idebug=0 -motion=rotate -plotOption=-1';  

  % PC2 with traction and slip wall BC 
  numChecks=numChecks+1;
  checkFileName{numChecks}='insPC2PolyRotatedSquareTraction.check';
  cmd{numChecks}='ins -ts=pc2 -tf=1 -ms=poly -degreex=2 -degreet=2 -knownSolution=none -nu=0.05 -cdv=0 -bcs=stsn -idebug=0 -map=rotatedSquare -N0=5 -plotOption=-1';

  % IM2 + traction BC + combined implicit
  numChecks=numChecks+1;
  checkFileName{numChecks}='insIM2PolyRotatedSquareTraction.check';
  cmd{numChecks}='ins -ts=im2 -tf=1 -ms=poly -degreex=2 -degreet=1 -knownSolution=none -nu=0.05 -cdv=0 -bcs=stnn -idebug=0 -map=rotatedSquare -N0=5 -plotOption=-1';


   numFailed=0; 
  for icheck=1:numChecks

    oldCheckFile = sprintf('%s/%s',checkDir,checkFileName{icheck});
    newCheckFile = 'ins.check';

    checkFileExists=isfile(oldCheckFile);
    % if( checkFileExists==0 )
    %   fprintf('File=[%s] does not exist.\n',oldCheckFile);
    % end
    if( checkFileExists==0 || par.replace==1 )
      newCheckFile = oldCheckFile;
    end
    
    cmd{icheck}=sprintf('%s -checkFileName=%s;',cmd{icheck},newCheckFile);

    fprintf('Running cmd=[%s]\n',cmd{icheck});
    eval(cmd{icheck});

    if( checkFileExists )
      % ------ DIFF OLD AND NEW CHECK FILES -----
      fprintf('diff check files: < new file,  > old file\n');
      rt(icheck) = system(sprintf('diff %s %s',newCheckFile,oldCheckFile));
      if( rt(icheck)==0 )
        fprintf('Test %d: %s : success.\n',icheck,oldCheckFile);
      else
        fprintf('Test %d: %s : failed. File [%s] did not match [%s]\n',icheck,oldCheckFile,newCheckFile,oldCheckFile);
        numFailed=numFailed+1;
      end
    else
      fprintf('Test %2d: %-30s : new file created.\n',icheck,checkFileName{icheck});
      rt(icheck)=12345;
    end

  end

  % ----- print a SUMMARY of results ----
  fprintf('\n ----- check SUMMARY ----\n');
  for icheck=1:numChecks
    if( rt(icheck)==0 )
      fprintf('Test %2d: %-40s : success.\n',icheck,checkFileName{icheck});
    elseif( rt(icheck)==12345 )
      fprintf('Test %2d: %-40s : new file created.\n',icheck,checkFileName{icheck});
    else
      fprintf('Test %2d: %-40s : failed.\n',icheck,checkFileName{icheck});
    end

  end
  if( numFailed==0 )
    fprintf('===== SUCCESS. All regression tests passed. ======\n');
  else
    fprintf('===== ERROR: %d regression tests FAILED. ======\n',numFailed);

  % ----- CHANGE TO CHECK DIRECTORY -----
 %%  cd ./check;

  return
end
