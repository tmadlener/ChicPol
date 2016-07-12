#!/bin/bash
#SBATCH -J DataFits
#SBATCH -D /afs/hephy.at/user/t/tmadlener/CMSSW_5_3_11/src/ChicPol/macros/polFit/
#SBATCH -o /afs/hephy.at/work/t/tmadlener/logfiles/slurm_arrays/runDataFits_%A_%a.out

echo "-------------------- "$(date)" --------------------"
## takes pt and rap bin as inputs!!!
ptBin=${1}
rapBin=${2}
nSample=${3}
nGenerations=${4}

JobID=${5}
storagedir=${6} # start from within storage directory
basedir=${7}
FidCuts=${8}
TreeID=${9}
datadir=${10}
nSkipGen=${11}
nState=${12}

muAccShift=0
if [ $# -gt 12 ]; then
  muAccShift=${13}
fi

# echo $TreeID
# echo $datadir
# echo $nSkipGen
# echo $nState

## fit steering
MPValgo=3         # 1...mean, 2...gauss, 3...gauss-loop with chi2<2

## other setup parameters
polScenSig=3
polScenBkg=3
frameSig=1
frameBkg=1

## efficiency setup (see runDataFits.sh for some more information on the different possibilities)
nEff=100003  # tmadlener, 24.06.2016: 100003 uses hybrid efficiencies for seagulls only, 100001 uses parametrized
# efficiencies (but for seagulls & cowboys combined!)
# WARNING: MAKE SURE TO CHANGE THE PREPROCESSOR FLAG USE_TF1_EFFICIENCIES in polFit.C accordingly!
# 100004 is MC truth (have not checked if this is the same as some other option)

DataType=DATA # MC for MC closure, DATA for data, concerns single muon efficiencies only!
UseMCeff=false
nDileptonEff=1
UseMCDileptoneff=true
nRhoFactor=1

## setup the different flags for running the polarization framework
StatVarTotBGfraction=0 # apply statistical fluctuations on f_background
StatVarTotBGmodel=0    # apply statistical fluctuations on Bg model
StatVarRho=0           # apply statistical fluctuations on rho factor
StatVarEff=0           # apply statistical fluctuations on the single muon efficiencies

useAmapApproach=false
nAmap=1
nDenominatorAmap=1
cutDeltaREllDpt=false

useRefittedChic=true

ConstEvents=15000
UseConstEv=true
gen=false
rec=false
fit=true
plot=false
NewAccCalc=false

## from some older batch system???? (leave this at zero for slurm job submission)
useBatch=0

## have to cd in to the directory for every array job! (Am not sure why this is not taken from the call site)
cd ${storagedir}/${JobID}

## get the current array job id and setup the files appropriately and start the fit
iFit=$SLURM_ARRAY_TASK_ID
# iFit=1
jobFileID=_rap_${rapBin}_pt_${ptBin}_Fit_${iFit}

execbl="polGenRecFitPlot"${jobFileID}

cp polGenRecFitPlot ${execbl}

./${execbl} ${iFit}ThisGen ${JobID}=JobID ${storagedir}=storagedir ${basedir}=basedir ${nGenerations}=nGenerations ${polScenSig}polScenSig ${frameSig}frameSig ${polScenBkg}polScenBkg ${frameBkg}frameBkg ${rapBin}_rapBinMin ${rapBin}rapBinMax ${ptBin}ptBinMin ${ptBin}ptBinMax ${nEff}nEff ${nDileptonEff}nDiEff ${FidCuts}FidCuts ${nSample}nSample ${ConstEvents}ConstEvents ${nSkipGen}nSkipGen UseConstEv=${UseConstEv} gen=${gen} rec=${rec} fit=${fit} plot=${plot} ${TreeID}=TreeID ${datadir}=realdatadir UseMCeff=${UseMCeff} UseMCDileptoneff=${UseMCDileptoneff} ${nRhoFactor}nRhoFactor ${MPValgo}MPValgo NewAccCalc=${NewAccCalc} useAmapApproach=${useAmapApproach} nAmap=${nAmap} nDenominatorAmap=${nDenominatorAmap} useBatch=${useBatch} StatVarTotBGfraction=${StatVarTotBGfraction} StatVarTotBGmodel=${StatVarTotBGmodel} StatVarRho=${StatVarRho} StatVarEff=${StatVarEff} ${nState}=nState dataType=${DataType} muAccShift=${muAccShift}

exitcode=$?

if [ ${exitcode} -eq 0 ]; then
  rm polGenRecFitPlot${jobFileID} ## do some cleanup
else
  exit ${exitcode}
fi

echo "========================= "$(date)" ========================="
