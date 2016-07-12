#!/bin/bash

## new script that makes better use of the slurm batch system available on the heplx
#################################################################################
# ## NOTE: only the bare minimum setup is done here and then passed further on  #
# ##       All other setup is done in batchRunDataFits.sh                       #
# ##       Also the setup of the target directory is done by another script that#
# ##       checks if the source files are there (but not their timesptamps) so  #
# ##       make sure that the most recent source files are there (or none)      #
#################################################################################

## ONLY ONE OF THESE FLAGS CAN BE SET AT A TIME (can be changed via cl arg)
FITTING=1 ## execute this first and when all batch jobs have finished execute merging
MERGING=0 ## although in principle this should check if all jobs are completed, proper working is not guaranteed!

## read arguments to alter the upper two values
if [ $# -gt 0 ]; then
  if [ ${1} == "fit" ]; then
    FITTING=1
    MERGING=0
  fi
  if [ ${1} == "merge" ]; then
    MERGING=1
    FITTING=0
  fi
fi

export CHIC_BINNING=1

# directory in which all this will take place and where the results will be stored
storagedir=/afs/hephy.at/work/t/tmadlener/ChiPol/results
# storagedir=/afs/hephy.at/work/t/tmadlener/ChiPol/testSlurmJobs ## testing directory
# directory from which the data files will be obtained
## NOTE: this makes this variable relative to the basedir! (This is important in the setup script!)
datadir_start=${basedir}/macros/DataFiles

## JobID and DataID setup. DataID is the name of the folder inside DataFiles (FidCuts is then computed from there)
## JPSI MC
# nState=4
# DataID=SetOfCuts11_jpsi_13May2016_MCclosure_rejCow_true_rejSea_false
# JobID=jpsi_19June2016_mcClosure_newEff
# TreeID=Psi1S

# for chicMassLow in 0.01; do # 0.05 0.15 0.2; do
# for nSig in 3.5; do
# for polReg in LSB RSB NP; do

for muAccShift in 0.25 0.5 -0.25 -0.5; do

## CHIC DATA
# nState=6
nState=$[$CHIC_BINNING+5] ## use the CHIC_BINNING defined in the shell to determine nState
# DataID=SetOfCuts11_chic_21June2016_rejCow_pt10Cut_chic$[$nState-5]Binning_corrLib_fLSB_${FracLSB}
# JobID=chic$[$nState-5]_21June2016_rejCow_cutDiMu10_chic$[$nState-5]Binning_corrLib_fLSB_${FracLSB}
# DataID=SetOfCuts11_chic_21June2016_rejCow_pt10Cut_chic$[$nState-5]Binning_corrLib
# JobID=chic$[$nState-5]_21June2016_rejCow_cutDiMu10_chic$[$nState-5]Binning_corrLib_mcTruthEff
# DataID=SetOfCuts11_chic_30June2016_rejCow_pt10Cut_chic$[$nState-5]Binning_c1massL_${chicMassLow} #_nSigPR_${nSig}_nSigNP_${nSig}
# JobID=chic$[$nState-5]_30June2016_rejCow_cutDiMu10_chic$[$nState-5]Binning_c1massL_${chicMassLow} #_nSigPR_${nSig}_nSigNP_${nSig}
# DataID=SetOfCuts11_chic_30June2016_rejCow_pt10Cut_chic$[$nState-5]Binning_pol${polReg}_c1massL_0.05
# JobID=chic$[$nState-5]_30June2016_rejCow_cutDiMu10_chic$[$nState-5]Binning_pol${polReg}_c1massL_0.05

DataID=SetOfCuts31_chic_11July2016_muAccShift_${muAccShift}
JobID=chic1_11July2016_muAccShift_${muAccShift}

# TreeID=chic$[$nState-5] # switch to Psi1S for MC closure

## CHIC MC
# nState=6
# DataID=SetOfCuts11_chic_21June2016_MCclosure_rejCow_pt10Cut_chic1Binning_corrLib
# JobID=chic$[$nState-5]_21June2016_MCclosure_rejCow_chic1Binning_corrLib
# DataID=SetOfCuts11_chic_06July2016_pGunMC_chic1Binning_acceptAll
# JobID=chic${CHIC_BINNING}_06July2016_pGunMC_chic${CHIC_BINNING}Binning_acceptAll_PtRapMassSubtr
# DataID=SetOfCuts11_chic_07July2016_pGunMC_chic${CHIC_BINNING}
# JobID=chic1_07July2016_pGunMC
# DataID=SetOfCuts11_chic_07July2016_officialMC_chic${CHIC_BINNING}
# JobID=chic1_07July2016_officialMC
TreeID=chic$[$nState-5]

datadir=${datadir_start}/${DataID}

## some directory setup
homedir=$(pwd)
cd ../../
basedir=$(pwd)
cd macros/polFit

## rap and pt ranges
rapBinMin=1
rapBinMax=1
ptBinMin=1
ptBinMax=5 # set to 4 for when chic2Binning is used

if [ ${nState} -eq 7 ]; then
  ptBinMax=4
fi

## number of samples and number of fits setup
nSample=20000
nFits=30
nSkipGen=0

## Get FidCuts form DataID
FidCuts=$(echo $DataID | awk -F'_' '{print $1}')
FidCuts=${FidCuts:9} # strip the SetOfCuts from the previous output

batchIDFile=BatchJobIDs.log

#####################################
# STARTING THE BATCH JOBS           #
#####################################
if [ ${FITTING} -eq 1 ]; then
  ## Do the setup of the directory only once
  bash setupPolFit.sh ${basedir} ${storagedir} ${JobID} ${datadir}

  if [ $? -ne 0 ]; then
    echo "ERROR DURING SETUP"
    exit 111
  fi

  cd ${storagedir}/${JobID}
  cp ${homedir}/batchRunDataFits.sh .

  ## start a slurm array job for each pt and rap bin that is wanted
  iRap=${rapBinMin}
  iPT=${ptBinMin}
  while [ ${iRap} -le ${rapBinMax} ]; do
    while [ ${iPT} -le ${ptBinMax} ]; do
      ## submit a job and save the output message to a variable
      subtext=$(sbatch --array=$[$nSkipGen+1]-${nFits} ${storagedir}/${JobID}/batchRunDataFits.sh ${iPT} ${iRap} ${nSample} ${nFits} ${JobID} ${storagedir} ${basedir} ${FidCuts} ${TreeID} ${basedir}${datadir}/tmpFiles ${nSkipGen} ${muAccShift})
      echo ${subtext}

      batchJobID=${subtext##* }
      echo ${iRap} ${iPT} ${batchJobID} >> ${batchIDFile} ## store all the batch job IDs in a file for easier checking afterwards

      ## testing without the batch system
      # bash batchRunDataFits.sh ${iPT} ${iRap} ${nSample} ${nFits} ${JobID} ${storagedir} ${basedir} ${FidCuts} ${TreeID} ${basedir}${datadir}/tmpFiles ${nSkipGen} ${nState} ${muAccShift}

      # echo "Submitted batch job for "$[$nFits-$nSkipGen]" fits for (pt,rap) bin ("${iPT}","${iRap}"). The job id is: "${batchJobID}

      iPT=$[$iPT+1]
    done
    iRap=$[$iRap+1]
  done
fi

######################################
# MERGING THE OUTPUT FILES           #
######################################
if [ ${MERGING} -eq 1 ]; then
  cd ${storagedir}/${JobID}

  ## read in the batch job ids one by one and check if they have succesfully finished
  while read line; do
    iPT=$(echo $line | awk '{print $2}')
    iRap=$(echo $line | awk '{print $1}')

    batchJobID=$(echo $line | awk '{print $3}') ## get batch job id
    ## the following only works if sacct returns only one unique status (check that first)
    if [ $(sacct -j ${batchJobID} -n --format=state | sort -u | wc -l) -eq 1 ]; then
      jobState=$(sacct -j ${batchJobID} -n --format=state | sort -u) ## if all wen't well this is only one entry
      exitCode=$(sacct -j ${batchJobID} -n --format=exitcode | uniq) ## as is this

      echo ${batchJobID}" status is "${jobState}" with exitCode "${exitCode}
      if [ ${jobState} == "COMPLETED" ] && [ ${exitCode} == "0:0" ]; then
        echo "merging for (rap, pt) bin ("${iRap}","${iPT}")"
        ## use hadd to add all the files that are present. Use -f flag for forced recration of output file
        hadd -v 1 -f results_${TreeID}_rap${iRap}_pT${iPT}.root results_Fit*_${TreeID}_rap${iRap}_pT${iPT}.root
        if [ $? -eq 0 ]; then
          echo "merge successful"
          sed -i '/'${batchJobID}'/d' ${batchIDFile}
          echo ${iRap} ${iPT} ${batchJobID} >> ${batchIDFile}.merged
        fi
      fi
    else
      echo ${batchJobID}" status is ambiguous! Run sacct -j "${batchJobID}" to see what is happening"
    fi
  done < ${batchIDFile}

  [ -s ${batchIDFile} ] || rm ${batchIDFile} # delete file if empty
fi

cd ${homedir} # return back to where we started

done # chic1MassLow
