#!/bin/sh

########## INPUTS ##########

nState=6	#6...chic1, 7...chic2

JobID=ToyMC_Test_20160530

#nGenerations=50
nGenerations=25

rapBinMin=1
rapBinMax=1
ptBinMin=5
ptBinMax=5

polScenSig=3
polScenBkg=3
frameSig=3
frameBkg=3

nEff=1050 # 105 MCtruth Jpsi
UseMCeff=false
nDileptonEff=1
UseMCDileptoneff=true
nRhoFactor=1

FidCuts=11

nSample=10000		#default:10000, includes burn-in
nSkipGen=0

#GENERATION SETTINGS
ConstEvents=1000
UseConstEv=true #if false, the number of events is taken from ToyMC.h

UseDifferingEff=false #if false, the next five lines do not matter
nEffRec=1050 #1101
UseMCReceff=false
nDileptonEffRec=1
UseMCDileptonReceff=true
nRecRhoFactor=1

gen=true
rec=true
fit=true
plot=false
deletePseudoData=false

MPValgo=3 		#1...mean,2...gauss,3...gauss-loop with chi2<2
nSigma=1
NewAccCalc=false

########################################

useAmapApproach=false       #if false, the next two lines do not matter, but must be given a valid value
nAmap=32104                 #frame/state/sigma/ID ( ID= 2 digits )
nDenominatorAmap=1 		    #the number here corresponds to the same notation as nEff

homedir=$PWD
cd ${homedir}
cd ..
cd ..
basedir=$PWD
cd macros/polFit
storagedir=/afs/cern.ch/work/k/knuenz/storage/ChicPol #please define the directory storagedir in the file macros/polFit/storagedir

ScenDir=Sig_frame${frameSig}scen${polScenSig}_Bkg_frame${frameBkg}scen${polScenBkg}

mkdir -p ${storagedir}/${JobID}/${ScenDir}

cp ${basedir}/macros/polFit/polGenRecFitPlot.cc ${storagedir}/${JobID}/${ScenDir}/polGenRecFitPlot.cc
cp ${basedir}/macros/polFit/polRapPtPlot.cc ${storagedir}/${JobID}/${ScenDir}/polRapPtPlot.cc
cp ${basedir}/macros/polFit/PlotFinalResults.cc ${storagedir}/${JobID}/${ScenDir}/PlotFinalResults.cc
cp ${basedir}/macros/polFit/Makefile ${storagedir}/${JobID}/${ScenDir}/Makefile
cp ${basedir}/macros/polFit/polGen.C ${storagedir}/${JobID}/${ScenDir}/polGen.C
cp ${basedir}/macros/polFit/polRec.C ${storagedir}/${JobID}/${ScenDir}/polRec.C
cp ${basedir}/macros/polFit/polFit.C ${storagedir}/${JobID}/${ScenDir}/polFit.C
cp ${basedir}/macros/polFit/polPlot.C ${storagedir}/${JobID}/${ScenDir}/polPlot.C

cp ../../interface/rootIncludes.inc ${storagedir}/${JobID}/${ScenDir}/rootIncludes.inc
cp ../../interface/commonVar.h ${storagedir}/${JobID}/${ScenDir}/commonVar.h
cp ../../macros/polFit/ToyMC.h ${storagedir}/${JobID}/${ScenDir}/ToyMC.h
cp ../../interface/effsAndCuts_chi.h ${storagedir}/${JobID}/${ScenDir}/effsAndCuts.h

cd ${storagedir}/${JobID}/${ScenDir}
cp ${basedir}/macros/polFit/runToyMC.sh .

touch polGenRecFitPlot.cc
make

rap_=${rapBinMin}
while [ $rap_ -le ${rapBinMax} ]
do
  pT_=${ptBinMin}
  while [ $pT_ -le ${ptBinMax} ]
  do

    nGen_=${nSkipGen}
    nGen_=$[nGen_+1]
    nMaxGen=$[nGenerations+nSkipGen]
    while [ $nGen_ -le $nMaxGen ]
    do

      plot=${plot}

      if [ ${nGen_} -eq 1 ]
      then
        plot=true
      fi


      cp polGenRecFitPlot polGenRecFitPlot_rap${rap_}_pt${pT_}_Gen${nGen_}
      #./polGenRecFitPlot_rap${rap_}_pt${pT_}_Gen${nGen_} ThisGen=${nGen_} JobID=${JobID} storagedir=${storagedir} basedir=${basedir} nGeneration=${nGenerations} polScenSig=${polScenSig} frameSig=${frameSig} polScenBkg=${polScenBkg} frameBkg=${frameBkg} rapBinMin=${rap_} rapBinMax=${rap_} ptBinMin=${pT_} ptBinMax=${pT_} nEff=${nEff} nDiEff=${nDileptonEff} nRecEff=${nEffRec} nRecDiEff=${nDileptonEffRec} FidCuts=${FidCuts} nSample=${nSample} ConstEvents=${ConstEvents} nSkipGen=${nSkipGen} UseConstEv=${UseConstEv} gen=${gen} rec=${rec} fit=${fit} plot=${plot} UseDifferingEff=${UseDifferingEff} UseMCeff=${UseMCeff} UseMCReceff=${UseMCReceff} UseMCDileptoneff=${UseMCDileptoneff} UseMCDileptonReceff=${UseMCDileptonReceff} nRhoFactor=${nRhoFactor} nRecRhoFactor=${nRecRhoFactor} MPValgo=${MPValgo} nSigma=${nSigma} nState=${nState} NewAccCalc=${NewAccCalc} deletePseudoData=${deletePseudoData} useAmapApproach=${useAmapApproach} nAmap=${nAmap} nDenominatorAmap=${nDenominatorAmap}
      ./polGenRecFitPlot_rap${rap_}_pt${pT_}_Gen${nGen_} ${nGen_}ThisGen ${JobID}=JobID ${storagedir}=storagedir ${basedir}=basedir ${nGenerations}nGenerations ${polScenSig}polScenSig ${frameSig}frameSig ${polScenBkg}polScenBkg ${frameBkg}frameBkg ${rap_}rapBinMin ${rap_}rapBinMax ${pT_}ptBinMin ${pT_}ptBinMax ${nEff}nEff ${nDileptonEff}nDiEff ${nEffRec}nRecEff ${nDileptonEffRec}nRecDiEff ${FidCuts}FidCuts ${nSample}nSample ${ConstEvents}ConstEvents ${nSkipGen}nSkipGen UseConstEv=${UseConstEv} gen=${gen} rec=${rec} fit=${fit} plot=${plot} UseDifferingEff=${UseDifferingEff} UseMCeff=${UseMCeff} UseMCReceff=${UseMCReceff} UseMCDileptoneff=${UseMCDileptoneff} UseMCDileptonReceff=${UseMCDileptonReceff}  ${nRhoFactor}nRhoFactor ${nRecRhoFactor}nRecRhoFactor ${MPValgo}MPValgo ${nSigma}nSigma ${nState}nState NewAccCalc=${NewAccCalc} deletePseudoData=${deletePseudoData} useAmapApproach=${useAmapApproach} ${nAmap}nAmap ${nDenominatorAmap}nDenominatorAmap
      rm polGenRecFitPlot_rap${rap_}_pt${pT_}_Gen${nGen_}

      nGen_=$[nGen_+1]
    done
    pT_=$[pT_+1]
  done
  rap_=$[rap_+1]
done


#rm polGen.C
#rm polRec.C
#rm polFit.C
#rm polPlot.C
