#!/bin/bash

############ INPUTS ####################
# source /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64-slc5/setup.sh
# source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.34.05/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh

## conditionally append to the JobID global variable
## NOTE: make sure to have the appropriate condition
append_jobID () {
  if [ $(echo "${1} != 0" | bc -l) -eq 1 ]; then
    JobID=${JobID}_${2}_${1}
  fi
}

## either do this before or here! building the executables won't work else-wise!
export CHIC_BINNING=1 # use chi_c1 or chi_c2 binning? ## exporting this, such that make can read it

Cdir=$PWD

cd ..
basedir=$PWD
cd macros

# input arguments
for nState in 6;do    #1,2,3,Upsi(1S,2S,3S); 4=Jpsi, 5=PsiPrime, 6=chic1 and chic2
  for FidCuts in 32; do #defines the set of cuts to be used, see macros/polFit/effsAndCuts.h
    cd $Cdir

    COPY_AND_COMPILE=0

    rapMin=1     #takes bins, not actual values
    rapMax=1     #if you only want to process 1 y bin, rapMax = rapMin
    ptMin=1      #takes bins, not acutal values
    ptMax=1      # chic1 pt binning
    if [ ${CHIC_BINNING} -eq 2 ]; then
      ptMax=4
    fi

    Plotting=1   #plotting macro: 1 = plot all, 2 = plot mass, 3 = plot lifetime
    #plotting macro: 4 = plot lifetimeSR1, 5 = plot lifetimeSR2, 6 = plot lifetimeLSB, 7 = plot lifetimeRSB, 8 = plot lifetimeFullRegion

    PlottingJpsi=1   #plotting macro: 1 = plot all, 2 = plot mass, 3 = plot lifetimeSBs, 4= lifetimeSR, 5= PlotMassRap, 6= plot lifetime, 7=pedagogical

    PlottingDataDists=1 #0...all, 1...1D plots, 2...2D plots

    runChiMassFitOnly=false
    incChic0InMassFit=true # include the Chic0 in the mass(-lifetime) fit of the chic
    includeBkgInMassFit=true # if false the chic mass model is simply the sum of the chic1 and the chic2 mass
    correctCtau=false   #correct pseudo-proper lifetime
    rejectCowboys=true
    rejectSeagulls=false # not checked if cowboys are rejected!
    RequestTrigger=true
    MC=false
    MCclosure=false # run MCclosure, CAUTION with setting what (and what not) to execute (not everything works!)
    drawRapPt2D=false  #draw Rap-Pt 2D map of Psi
    FixRegionsToInclusiveFit=false
    rapFixTo=1
    ptFixTo=1
    doFractionUncer=true #chic. false for MCclosure
    useRefittedChic=true # use the refitted mass for chic or use M_chic - M_jpsi + M_jpsi_pdg, NOTE: this will be set to false in bkgHistos_leptonBased.C if onia::KinParticleChi == false
    subtractNP=false #default == false
    cutDimuon10Gev=true # apply a 10 GeV cut in dimuon pt (at event selection stage)
    singleMuAccShift=5.6 # overall shift that is applied to the single muon acceptance pT cuts for FidCuts 31 or the uniform acceptance cut for FidCuts 32

    #PlotFitPar:::
    AddInclusiveResult=false #Inclusive defined by rapFixTo, ptFixTo

    doCtauUncer=true   #redundant for chic
    PolLSB=false       #measure polarization of the left sideband
    PolRSB=false       #measure polarization of the right sideband
    PolNP=false      #measure polarization of the non prompt events
    forceBinning=true  #set binning of Psi1S consistently to non prompt binning and Psi2S consistently to background binning #redundtant for chic
    folding=true       #folding is applied to all background histograms
    normApproach=false #normalization #used before binning algorithm #redundtant
    ctauScen=3         #0:default(1s:2.5,2s:2.0), 1:(1s:3.5,2s:3.0), 2:(1s:1.5,2s:1.0), 3:100mm 1S and 2S, #redundant for chic
    FracLSB=-1         #-1:defalut, 0, 100
    scaleFracBg=false  # not used in chic
    fitMassPR=false
    fitMassNP=false

    ## signal region definitions (default is -1 -> see commonVar.h to check used values, between 0 and 1 for adjusted ranges)
    chic1MassLow=-1
    chic1MassHigh=-1
    chic2MassLow=-1
    chic2MassHigh=-1
    nSigPR=-1
    nSigNP=-1

    injectParams=false # change the values of fit values to other user-defined ones (defined in the $altFileName file)
    # file in which fractions to be used in runBkgHistos_new when no fit values are present
    altFileName=defineFitVars.in
    fillBkgRandom=false # fill the bkg histograms (cosThetaPhi and pTRapMass) with randomly drawn events (particle gun MC for non-empty bkg histos)

    DataID=Psi$[nState-3]S_ctauScen0_FracLSB-1_16Mar2013
    polDataPath=${basedir}/Psi/Data/${DataID}

    #Define JobID
    # JobID=chic_11April2016_nonRefit_useRef_${useRefittedChic}_rejCBs # fChi1MassLow = 0.1
    # JobID=jpsi_13May2016_MCclosure_rejCow_${rejectCowboys}_rejSea_${rejectSeagulls}
    # JobID=chic_23May2016_MCclosure_test
    JobID=chic_17July2016_rejCow_pt10Cut_chic1Binning
    # JobID=chic_21June2016_rejCow_pt10Cut_chic2Binning_corrLib # added for comparison with old library
    # JobID=chic_15June2016_rejCow_pt10Cut_chic2Binning
    # JobID=chic_30March2016_ML30 # fChi1MassLow = 0.3

    ## safety feature to have unique job ids for different settings of these
    # append_jobID ${FracLSB} fLSB
    # append_jobID ${chic1MassLow} c1massL
    # append_jobID ${chic2MassLow} c2massL
    # append_jobID ${chic1MassHigh} c1massH
    # append_jobID ${chic2MassHigh} c2massH
    # append_jobID ${nSigPR} nSigPR
    # append_jobID ${nSigNP} nSigNP

    append_jobID ${singleMuAccShift} accCut

    ################ EXECUTABLES #################

    #following flags decide if the step is executed (1) or not (0):
    #IMPORTANT: for MC set execute_runWorkspace, execute_MassFit and execute_runLifetimeFit to 0
    execute_runChiData=0			           		#independent of rapMin, rapMax, ptMin, ptMax
    execute_runWorkspace=0	    					#independent of rapMin, rapMax, ptMin, ptMax
    execute_runMassFit=0				    	    #can be executed for different pt and y bins
    execute_runLifetimeFit=1  				    #can be executed for different pt and y bins
    execute_runPlotJpsiMassLifetime=0    			#can be executed for different pt and y bins
    execute_PlotJpsiFitPar=0              			#can be executed for different pt and y bins
    execute_runChiMassLifetimeFit=1		  	    	#can be executed for different pt and y bins
    execute_runDefineRegionsAndFractions=1			#can be executed for different pt and y bins
    execute_runPlotMassLifetime=0   				#can be executed for different pt and y bins
    execute_PlotFitPar=0              				#can be executed for different pt and y bins
    execute_runPlotDataDistributions=0		 		#This step only has to be executed once for each set of cuts (indep. of FracLSB and nSigma)
    execute_runBkgHistos=1          				#can be executed for different pt and y bins
    execute_PlotCosThetaPhiBG=0 		 			#This step only has to be executed once for each set of cuts (indep. of FracLSB and nSigma)
    execute_PlotMassRapPtBG=0 		 			#This step only has to be executed once for each set of cuts (indep. of FracLSB and nSigma)
    execute_PlotCosThetaPhiDistribution=0 			#This step only has to be executed once for each set of cuts (indep. of FracLSB and nSigma)

    #################################
    #PsiRelics:::
    execute_runData=0			           #independent of rapMin, rapMax, ptMin, ptMax


    # input files
    # In case of more input Files: define inputTreeX and adapt the line starting with inputTrees, at the moment up to 4 files implemented
    if [ ${nState} -eq 4 ]
    then
      inputTree1=/scratch/knuenz/Polarization/RootInput/Psi/TTree_Onia2MuMu_v30_PromptRecoAB_10May2012_Jpsi.root
      if [ ${MC} = 'true' ]
      then
        inputTree1=/scratch/ikratsch/Polarization/Jpsi/InputFiles/TTree_Psi1S_Gun_Pt9p5_70p5_19Dec2012.root
      fi
      if [ ${MCclosure} = 'true' ]; then
        inputTree1=/afs/hephy.at/data/ikraetschmer01/ikraetschmer/TnP2012/InputFiles/validation/onia2MuMu_tree_validation.root
      fi
    fi

    if [ ${nState} -eq 5 ]
    then
      inputTree1=/scratch/knuenz/Polarization/RootInput/Psi/TTree_Onia2MuMu_v30_PromptRecoAB_10May2012_Psi.root
      if [ ${MC} = 'true' ]
      then
        inputTree1=/scratch/ikratsch/Polarization/Jpsi/InputFiles/TTree_Psi2S_Gun_Pt6p5_50p5_19Dec2012.root
      fi
    fi


    if [ ${nState} -eq 6 ]
    then
      inputTree1=/scratch/knuenz/Polarization/RootInput/ChicPol/chic_rootuple_subFeb2014.root
      if [ ${MC} = 'true' ]
      then
        inputTree1=/scratch/knuenz/Polarization/RootInput/ChicPol/chic_rootuple_MC_15M_sel.root
      fi
      if [ ${MCclosure} = 'true' ]; then
        inputTree1=/afs/hephy.at/data/ikraetschmer01/ikraetschmer/Polarization/2012/InputFiles/ChiC/tuple-mc-chic-rec.root
      fi
    fi




    # Make directories
    CutDir=${Cdir}/DataFiles/SetOfCuts${FidCuts}_${JobID}

    WorkDir=${CutDir}
    mkdir -p ${CutDir}
    mkdir -p ${WorkDir}

    mkdir -p DataFiles
    mkdir -p ${WorkDir}/tmpFiles/backupWorkSpace
    mkdir -p ${WorkDir}/Figures
    mkdir -p ${WorkDir}/PDF
    mkdir -p ${WorkDir}/Fit


    if [ ${COPY_AND_COMPILE} -eq 1 ]
    then

      cp ../interface/commonVar.h ${WorkDir}/commonVar.h
      # Copy files to directory
      cp Makefile ${WorkDir}/Makefile
      cp ../interface/rootIncludes.inc ${WorkDir}/rootIncludes.inc
      cp ../interface/clarg_parsing.h  ${WorkDir}/clarg_parsing.h

      cp runData.cc ${WorkDir}/runData.cc
      cp runChiData.cc ${WorkDir}/runChiData.cc
      cp PolData.C ${WorkDir}/PolData.C
      cp PolData.h ${WorkDir}/PolData.h
      cp PolChiData.C ${WorkDir}/PolChiData.C
      cp PolChiData.h ${WorkDir}/PolChiData.h
      cp ../interface/effsAndCuts_chi.h ${WorkDir}/effsAndCuts.h

      cp runMC_helperStructs.h ${WorkDir}/runMC_helperStructs.h
      cp runMC.cc ${WorkDir}/runMC.cc

      cp runWorkspace.cc ${WorkDir}/runWorkspace.cc
      cp createWorkspace.C ${WorkDir}/createWorkspace.C

      cp runMassFit.cc ${WorkDir}/runMassFit.cc
      cp massFit.cc ${WorkDir}/massFit.cc

      cp runChiMassLifetimeFit.cc ${WorkDir}/runChiMassLifetimeFit.cc
      cp chiMassLifetimeFit.cc ${WorkDir}/chiMassLifetimeFit.cc
      cp runDefineRegionsAndFractions.cc ${WorkDir}/runDefineRegionsAndFractions.cc
      cp DefineRegionsAndFractions.cc ${WorkDir}/DefineRegionsAndFractions.cc



      cp runLifetimeFit.cc ${WorkDir}/runLifetimeFit.cc
      cp lifetimeFit.cc ${WorkDir}/lifetimeFit.cc
      cp ../interface/calculatePar.cc ${WorkDir}/calculatePar.cc
      cp ../interface/RooUtils.h ${WorkDir}/RooUtils.h

      cp runPlotMassLifetime.cc ${WorkDir}/runPlotMassLifetime.cc
      cp PlotMassLifetime.cc ${WorkDir}/PlotMassLifetime.cc

      cp runPlotJpsiMassLifetime.cc ${WorkDir}/runPlotJpsiMassLifetime.cc
      cp PlotJpsiMassLifetime.cc ${WorkDir}/PlotJpsiMassLifetime.cc

      cp runPlotDataDistributions.cc ${WorkDir}/runPlotDataDistributions.cc
      cp PlotDataDistributions.cc ${WorkDir}/PlotDataDistributions.cc

      cp PlotFitPar.cc ${WorkDir}/PlotFitPar.cc
      cp PlotJpsiFitPar.cc ${WorkDir}/PlotJpsiFitPar.cc

      cp runBkgHistos.cc ${WorkDir}/runBkgHistos.cc
      cp bkgHistos.C ${WorkDir}/bkgHistos.C
      cp bkgHistos_leptonBased.C ${WorkDir}/bkgHistos_chi.C
      cp calcPol.C ${WorkDir}/calcPol.C
      cp bkgHistos_MCclosure.C  ${WorkDir}/bkgHistos_MCclosure.C

      cp runBkgHistos_new.cc ${WorkDir}/runBkgHistos_new.cc
      cp BkgHistoProducer*.h ${WorkDir}/
      cp bkgHistos_helper.h ${WorkDir}/bkgHistos_helper.h

      cp PlotCosThetaPhiBG.cc ${WorkDir}/PlotCosThetaPhiBG.cc
      cp PlotMassRapPtBG.cc ${WorkDir}/PlotMassRapPtBG.cc
      cp PlotCosThetaPhiDistribution.cc ${WorkDir}/PlotCosThetaPhiDistribution.cc

      #cp ../latex/Mass_fitParameter.tex ${WorkDir}/Mass_fitParameter.tex
      #cp ../latex/Lifetime_fitParameter.tex ${WorkDir}/Lifetime_fitParameter.tex
      #cp ../latex/myStyle.tex ${WorkDir}/myStyle.tex
      #cp ../latex/evaluateCtau.tex ${WorkDir}/evaluateCtau.tex
      #cp ../latex/NumEvents.tex ${WorkDir}/NumEvents.tex

      cp ../latex/cosThetaPhi_chic1_BG.tex        ${WorkDir}/cosThetaPhi_chic1_BG.tex
      cp ../latex/cosThetaPhi_chic2_BG.tex        ${WorkDir}/cosThetaPhi_chic2_BG.tex
      #cp ../latex/cosThetaPhi_BG.tex        ${WorkDir}/cosThetaPhi_BG.tex
      #cp ../latex/cosThetaPhi_BG_highct.tex ${WorkDir}/cosThetaPhi_BG_highct.tex
      #cp ../latex/cosThetaPhi_NPBG.tex      ${WorkDir}/cosThetaPhi_NPBG.tex
      #cp ../latex/cosThetaPhi_TBG.tex       ${WorkDir}/cosThetaPhi_TBG.tex
      #cp ../latex/cosThetaPhi.tex           ${WorkDir}/cosThetaPhi.tex
      cp ../latex/massRapPt_chic1_BG.tex        ${WorkDir}/massRapPt_chic1_BG.tex
      cp ../latex/massRapPt_chic2_BG.tex        ${WorkDir}/massRapPt_chic2_BG.tex
      cp ../latex/cosThetaPhi_chic1.tex           ${WorkDir}/cosThetaPhi_chic1.tex
      cp ../latex/cosThetaPhi_chic2.tex           ${WorkDir}/cosThetaPhi_chic2.tex
      #cp ../latex/MassLifetime.tex       ${WorkDir}/MassLifetime.tex

      cd ${WorkDir}

      make -k

    fi

    if [ ${injectParams} = "true" ]; then
      cp ${basedir}/macros/${altFileName} ${WorkDir}
    fi
    cd ${WorkDir}

    inputTrees="inputTree=${inputTree1} inputTree=${inputTree2} inputTree=${inputTree3} inputTree=${inputTree4}"
    if [ ${execute_runData} -eq 1 ]
    then
      ./runData ${inputTrees} rejectCowboys=${rejectCowboys} FidCuts=${FidCuts} nState=${nState} MC=${MC} RequestTrigger=${RequestTrigger}
    fi

    if [ ${execute_runChiData} -eq 1 ]
    then
      ## TODO: the chic MC (closure) can be handled by the normal runChiData!!
      # if [ ${MCclosure} = 'true' ]; then
      #   ./runMC ${inputTrees} rejectCowboys=${rejectCowboys} FidCuts=${FidCuts} nState=${nState} RequestTrigger=${requestTrigger} rejectSeagulls=${rejectSeagulls}
      # else
        ./runChiData ${inputTrees} rejectCowboys=${rejectCowboys} FidCuts=${FidCuts} nState=${nState} MC=${MC} RequestTrigger=${RequestTrigger} correctCtau=${correctCtau} useRefittedChic=${useRefittedChic} cutDimuon10Gev=${cutDimuon10Gev} muAccShift=${singleMuAccShift} rejectSeagulls=${rejectSeagulls}
      # fi
    fi

    if [ ${execute_runWorkspace} -eq 1 ]
    then
      ./runWorkspace nState=${nState} correctCtau=${correctCtau} drawRapPt2D=${drawRapPt2D} useRefittedChic=${useRefittedChic} mcClosure=${MCclosure}
    fi

    if [ ${execute_runMassFit} -eq 1 ]
    then
      cp runMassFit runMassFit_rap${rapMin}_pt${ptMin}
      ./runMassFit_rap${rapMin}_pt${ptMin} rapMin=${rapMin} rapMax=${rapMax} ptMin=${ptMin} ptMax=${ptMax} nState=${nState} fitMassPR=${fitMassPR} fitMassNP=${fitMassNP} MC=${MC}
      rm runMassFit_rap${rapMin}_pt${ptMin}
    fi

    if [ ${execute_runLifetimeFit} -eq 1 ]
    then
      cp runLifetimeFit runLifetimeFit_rap${rapMin}_pt${ptMin}
      ./runLifetimeFit_rap${rapMin}_pt${ptMin} rapMin=${rapMin} rapMax=${rapMax} ptMin=${ptMin} ptMax=${ptMax} nState=${nState} MC=${MC}
      rm runLifetimeFit_rap${rapMin}_pt${ptMin}
    fi

    if [ ${execute_runPlotJpsiMassLifetime} -eq 1 ]
    then
      cp runPlotJpsiMassLifetime runPlotJpsiMassLifetime_$[nState-3]S_rap${rapMin}_pt${ptMin}
      ./runPlotJpsiMassLifetime_$[nState-3]S_rap${rapMin}_pt${ptMin} rapMin=${rapMin} rapMax=${rapMax} ptMin=${ptMin} ptMax=${ptMax} nState=${nState} PlottingJpsi=${PlottingJpsi}
      rm runPlotJpsiMassLifetime_$[nState-3]S_rap${rapMin}_pt${ptMin}
    fi

    if [ ${execute_PlotJpsiFitPar} -eq 1 ]
    then
      ./PlotJpsiFitPar nState=${nState} doCtauUncer=${doCtauUncer} rapMin=${rapMin} rapMax=${rapMax} ptMin=${ptMin} ptMax=${ptMax} rapFixTo=${rapFixTo} ptFixTo=${ptFixTo} AddInclusiveResult=${AddInclusiveResult}
    fi

    if [ ${execute_runChiMassLifetimeFit} -eq 1 ]
    then
      cp runChiMassLifetimeFit runChiMassLifetimeFit_rap${rapMin}_pt${ptMin}
      ./runChiMassLifetimeFit_rap${rapMin}_pt${ptMin} runChiMassFitOnly=${runChiMassFitOnly} rapMin=${rapMin} rapMax=${rapMax} ptMin=${ptMin} ptMax=${ptMax} nState=${nState} MC=${MC} MCclosure=${MCclosure} useChic0=${incChic0InMassFit} useBkgMassFit=${includeBkgInMassFit}
      rm runChiMassLifetimeFit_rap${rapMin}_pt${ptMin}
    fi

    if [ ${execute_runDefineRegionsAndFractions} -eq 1 ]
    then
      rootfile=fit_Chi_rap${rapMin}_pt${ptMin}.root
      #cp tmpFiles/backupWorkSpace/ws_MassLifetimeFit_Chi_rap${rapMin}_pt${ptMin}.root tmpFiles/backupWorkSpace/ws_DefineRegionsAndFractions_Chi_rap${rapMin}_pt${ptMin}.root
      cp runDefineRegionsAndFractions runDefineRegionsAndFractions_rap${rapMin}_pt${ptMin}
      ./runDefineRegionsAndFractions_rap${rapMin}_pt${ptMin} runChiMassFitOnly=${runChiMassFitOnly} doFractionUncer=${doFractionUncer} rapMin=${rapMin} rapMax=${rapMax} ptMin=${ptMin} ptMax=${ptMax} nState=${nState} FixRegionsToInclusiveFit=${FixRegionsToInclusiveFit} rapFixTo=${rapFixTo} ptFixTo=${ptFixTo} chic1MassLow=${chic1MassLow} chic1MassHigh=${chic1MassHigh} chic2MassLow=${chic2MassLow} chic2MassHigh=${chic2MassHigh} nSigPR=${nSigPR} nSigNP=${nSigNP}
      rm runDefineRegionsAndFractions_rap${rapMin}_pt${ptMin}
    fi

    if [ ${execute_runPlotDataDistributions} -eq 1 ]
    then
      cp runPlotDataDistributions runPlotDataDistributions_$[nState-3]S_rap${rapMin}_pt${ptMin}
      mkdir Figures/PlotDataDistributions
      ./runPlotDataDistributions_$[nState-3]S_rap${rapMin}_pt${ptMin} rapMin=${rapMin} rapMax=${rapMax} ptMin=${ptMin} ptMax=${ptMax} nState=${nState} PlottingDataDists=${PlottingDataDists}
      rm runPlotDataDistributions_$[nState-3]S_rap${rapMin}_pt${ptMin}
    fi

    if [ ${execute_runPlotMassLifetime} -eq 1 ]
    then
      cp runPlotMassLifetime runPlotMassLifetime_$[nState-3]S_rap${rapMin}_pt${ptMin}
      ./runPlotMassLifetime_$[nState-3]S_rap${rapMin}_pt${ptMin} rapMin=${rapMin} rapMax=${rapMax} ptMin=${ptMin} ptMax=${ptMax} nState=${nState} Plotting=${Plotting}
      rm runPlotMassLifetime_$[nState-3]S_rap${rapMin}_pt${ptMin}
    fi

    if [ ${execute_PlotFitPar} -eq 1 ]
    then
      ./PlotFitPar nState=${nState} doCtauUncer=${doCtauUncer} rapMin=${rapMin} rapMax=${rapMax} ptMin=${ptMin} ptMax=${ptMax} rapFixTo=${rapFixTo} ptFixTo=${ptFixTo} AddInclusiveResult=${AddInclusiveResult}
      #pdflatex Lifetime_fitParameter.tex
      #pdflatex Mass_fitParameter.tex
      #pdflatex evaluateCtau.tex
      #pdflatex evaluateCtau.tex
      #pdflatex NumEvents.tex
      #pdflatex NumEvents.tex
      #mv Lifetime_fitParameter.pdf PDF/Lifetime_fitParameter.pdf
      #mv Mass_fitParameter.pdf PDF/Mass_fitParameter.pdf
      #mv evaluateCtau.pdf PDF/evaluateCtau.pdf
      #mv NumEvents.pdf PDF/NumEvents.pdf
    fi

    if [ ${execute_runBkgHistos} -eq 1 ]
    then
      if [ ${MCclosure} = 'true' ]; then
        cp runBkgHistos_new runBkgHistos_$[nState-3]S_rap${rapMin}_pt${ptMin}
      else
        cp runBkgHistos runBkgHistos_$[nState-3]S_rap${rapMin}_pt${ptMin}
      fi
      ./runBkgHistos_$[nState-3]S_rap${rapMin}_pt${ptMin} rapMin=${rapMin} rapMax=${rapMax} ptMin=${ptMin} ptMax=${ptMax} nState=${nState} MC=${MC} doCtauUncer=${doCtauUncer} PolLSB=${PolLSB} PolRSB=${PolRSB} PolNP=${PolNP} ctauScen=${ctauScen} FracLSB=${FracLSB} forceBinning=${forceBinning} folding=${folding} normApproach=${normApproach} scaleFracBg=${scaleFracBg} polDataPath=${polDataPath} subtractNP=${subtractNP} useRefittedChic=${useRefittedChic} mcClosure=${MCclosure} injectionFile=${altFileName} randomFill=${fillBkgRandom} injectParams=${injectParams}
      rm runBkgHistos_$[nState-3]S_rap${rapMin}_pt${ptMin}
    fi

    if [ ${execute_PlotCosThetaPhiBG} -eq 1 ]
    then
      ./PlotCosThetaPhiBG nState=${nState}
      pdflatex cosThetaPhi_chic1_BG.tex
      pdflatex cosThetaPhi_chic2_BG.tex
      mv cosThetaPhi_chic1_BG.pdf PDF/cosThetaPhi_chic1_BG.pdf
      mv cosThetaPhi_chic2_BG.pdf PDF/cosThetaPhi_chic2_BG.pdf
    fi

    if [ ${execute_PlotMassRapPtBG} -eq 1 ]
    then
      ./PlotMassRapPtBG nState=${nState}
      pdflatex massRapPt_chic1_BG.tex
      pdflatex massRapPt_chic2_BG.tex
      mv massRapPt_chic1_BG.pdf PDF/massRapPt_chic1_BG.pdf
      mv massRapPt_chic2_BG.pdf PDF/massRapPt_chic2_BG.pdf
    fi

    if [ ${execute_PlotCosThetaPhiDistribution} -eq 1 ]
    then
      ./PlotCosThetaPhiDistribution nState=${nState} DataPath=${WorkDir}
      pdflatex cosThetaPhi_chic1.tex
      pdflatex cosThetaPhi_chic2.tex
      mv cosThetaPhi_chic1.pdf PDF/cosThetaPhi_chic1.pdf
      mv cosThetaPhi_chic2.pdf PDF/cosThetaPhi_chic2.pdf
    fi

    #rm runData
    #rm runChiData
    #rm runWorkspace
    #rm runMassFit
    #rm runLifetimeFit
    #rm runPlotMassLifetime
    #rm runBkgHistos
    #rm PlotFitPar
    #rm PlotCosThetaPhiBG
    #rm PlotCosThetaPhiDistribution
    #rm *.tex
    #rm *.aux
    #rm *.log
    #rm *.so
    #rm *.d
    #rm *.nav
    #rm *.out
    #rm *.snm
    #rm *.toc

  done
done
