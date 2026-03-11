%
% Save a check file, used for regression tesing
% 
function writeCheckFile( par )

  % -------- write a check file for regression testing -----------------------
  checkFile = fopen(par.checkFileName,'w');
  fprintf(checkFile,'--------------------- INS Check File --------------------------\n');
  fprintf(checkFile,'ms=%s, ts=%s, bcs=%s, knownSolution=%s\n',par.ms,par.ts,par.bcs,par.knownSolution);
  fprintf(checkFile,'nu=%9.2e, cfl=%9.2e, tFinal=%9.2e, cdv=%9.2e\n',par.nu,par.cfl,par.tFinal,par.cdv);
  fprintf(checkFile,'ad=%d, ad21=%9.2e, ad22=%9.2e\n',par.ad,par.ad21,par.ad22);

  pErrMax = par.maxErr(1);
  uErrMax = par.maxErr(2);
  vErrMax = par.maxErr(3);
  divMax  = par.maxErr(4);
  fprintf(checkFile,'Nx=%d Ny=%d Nt=%d dt=%8.2e\n',par.Nx,par.Ny,par.Nt,par.dt);
  fprintf(checkFile,'uNorm=%8.2e, vNorm=%8.2e, pNorm=%8.2e\n',par.uNorm(1),par.uNorm(2),par.uNorm(3));
  fprintf(checkFile,'maxDivU=%8.2e, maxGradU=%8.2e\n',par.maxDivU,par.maxGradU);
  if( par.computeErrors )
    fprintf(checkFile,'max-Err(p,u,v)=(%8.2e,%8.2e,%8.2e)\n',pErrMax,uErrMax,vErrMax);
  end


  % fprintf(checkFile,'nd=%d;\n',par.nd);
  % fprintf(checkFile,'Np=%d;\n',par.Np);
  % fprintf(checkFile,'numToDeflate=%d;\n',par.numToDeflate);
  % fprintf(checkFile,'numFreq=%d;\n',par.numFreq);
  % for freq=1:par.numFreq
  %   fprintf(checkFile,'frequency(%d)=%14.8e;\n',freq,par.frequencyArraySaved(freq));
  % end
  % fprintf(checkFile,'damp=%g, filterTimeDerivative=%d, useFilterWeights=%d, viFactor=%g \n',...
  %       par.damp,par.filterTimeDerivative,par.useFilterWeights,par.viFactor);
  % fprintf(checkFile,'isign=%g : exp( isign*I*omegat), phase=%g : for cos( omega*t + phase*pi ) forcing \n',...
  %       par.isign,par.phase);
  % fprintf(checkFile,'upwind=%d (WaveHoltz), upwindh=%d (Helmholtz), upwindScaleFactor=%g, implicitUpwind=%d\n',...
  %       par.upwind,par.upwindh,par.upwindScaleFactor,par.implicitUpwind);   
  % fprintf(checkFile,'useUpwindCoeff=%d : 1= compute UPW coeff arrays uwx(:,:), uwy(:,:) using getUpwindCoefficients\n', par.useUpwindCoeff); 
  % if( par.useUpwindCoeff && par.superGrid )
  %   fprintf(checkFile,'useUpwindCoeff=1 and superGrid=1 : turn off dissipation in the interior\n'); 
  % end       
  % fprintf(checkFile,'superGrid=%d, superGridWidth=%.2g,  superGridh=%i (HelmHoltz)\n',par.superGrid,par.superGridWidth,par.superGridh);            
  % fprintf(checkFile,'Nf=%d, Np=%d, adjustOmega=%d, useDiscreteBetaFunction=%d\n',par.numFreq,par.Np,par.adjustOmega,par.useDiscreteBetaFunction);
  % fprintf(checkFile,'implicit=%d, tFinal=%8.4f, minStepsPerPeriod=%d, Nt=%d, adjustForcingForImplicitTimeStepping=%d\n',...
  %      par.implicit,par.tFinal,par.minStepsPerPeriod,par.Nt,par.adjustForcingForImplicitTimeStepping);
  % fprintf(checkFile,'deflateWaveHoltz=%d, agmres=%d, abicgstab=%d, deflateForcing=%d\n',par.deflateWaveHoltz,par.agmres,par.abicgstab,par.deflateForcing);  

  % fprintf(checkFile,'ax=%14.8e; bx=%14.8e; ay=%14.8e; by=%14.8e;\n',par.ax,par.bx,par.ay,par.by);
  % fprintf(checkFile,'damp=%14.8e;\n',par.damp);
  % fprintf(checkFile,'Nv=[%d,%d];\n',par.Nv(1),par.Nv(2));
  % fprintf(checkFile,'dxv=[%14.8e,%14.8e];\n',par.dxv(1),par.dxv(2));
  % fprintf(checkFile,'ppw=%14.8e;\n',par.ppw);
  % fprintf(checkFile,'restart=%14.8e; %% GMRES restart\n',par.restart);
  % fprintf(checkFile,'maxit=%14.8e; %% GMRES maxit\n',par.maxit);
  % fprintf(checkFile,'bcLabel=''%s'';\n',par.bcLabel);
  % fprintf(checkFile,'useDiscreteBetaFunction=%d\n',par.useDiscreteBetaFunction);
  

  % if( par.computeEigenmodes==0 )
  %   % ------- WAVEHOLTZ ------
  %   fprintf(checkFile,'totalWaveSolves=%d;\n',par.totalWaveSolves);
  %   fprintf(checkFile,'applyPreconditioner=%d;\n',par.applyPreconditioner);
   
  %   fprintf(checkFile,'numPC=%d;\n',par.numPC);
  %   for( ipc=1:par.numPC )
  %     fprintf(checkFile,'pc{%d}.damp=%14.8e;\n',ipc,par.pc{ipc}.damp);
  %     fprintf(checkFile,'pc{%d}.its=%d;\n',ipc,par.pc{ipc}.its);
  %   end  

  %   for( freq=1:par.numFreq )
  %     fprintf(checkFile,"freq=%d: max || H v -f ||/omega^2  = %8.2e (residual in WaveHoltz solution v)\n",freq,par.maxResidual(freq));
  %     if( par.computeEigenmodes==0 )
  %       fprintf(checkFile,"freq=%d: max || v - Helmholtz ||   = %8.2e (error in WaveHoltz - Helmholtz)\n",freq,par.errHelmholtz(freq));
  %     end
  %   end

  %   if par.useFixPoint
  %     fprintf(checkFile,"FPI:    average-CR=%4.2f, ACR(deflate=%d)=%4.2f, its=%d\n",par.aveCR,par.numToDeflate,par.rateDeflated,par.numFixPointIterations);
  %   end
  %   if par.useGMRES
  %     fprintf(checkFile,"%s: average-CR=%4.2f, FPI-ACR=%4.2f, its=%d, ||res||_h=%9.2e [wave-solves=%d]\n", ...
  %          upper(par.krylovType),par.aveKrylovCR,par.rateDeflated,par.numKrylovIterations,par.resKrylovMax,par.totalWaveSolves);
  %   end

  % else
  %   % ----- EIGENWAVE -------


  %  fprintf(checkFile,'numEigsToCompute=%d, numComputed=%d, numWaveSolves=%d, waveSolves-per-eig=%.3g\n',...
  %      par.numEigsToCompute,par.numComputed,numberOfMatrixVectorProducts,numberOfMatrixVectorProducts/max(1,par.numComputed));

  %   for( j=1:par.numComputed )

  %     [minDist,lamTrueIndex] = min(abs(par.lambdaRQ(j)-par.lambdaDiscrete));

  %     par.eigNumber(j) = lamTrueIndex;      % used to mark deflated eigs

  %     lamTrue = par.lambdaDiscrete(lamTrueIndex);
  %     relErr = abs( par.lambdaRQ(j) - lamTrue )/max(1,lamTrue);

  %     betaTrue = getWaveHoltzIterationEigenvalue( par.lambdaRQ(j),par );
  %     if( par.filterTimeDerivative==0 )
  %       fprintf(checkFile,'j=%3d: [beta,betaRQ]=[%9.2e,%9.2e], eig=%3d, par.lambdaRQ= %9.6g, lamTrue=%9.6g, rel-err-lam=%9.2e\n',...
  %             j,par.lamBeta(j),betaTrue,lamTrueIndex,par.lambdaRQ(j),lamTrue,relErr);
  %     else
  %       fprintf(checkFile,'j=%3d: par.lambdaRQ=[%10.6g,%10.6g]\n',j,real(par.lambdaRQ(j)),imag(par.lambdaRQ(j)));
  %     end
  %   end
  %   for( j=1:par.numComputed )
  %     fprintf(checkFile,"j=%2d : rel-residual  || A phi - lam^2 phi || / lam^2 = %9.3e (max-norm)\n",j,par.maxEigRes(j));
  %   end

  % end

  fclose(checkFile);
  if( par.idebug>0 ) fprintf('wrote checkFile=[%s]\n',par.checkFileName); end

end