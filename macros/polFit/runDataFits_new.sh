#!/bin/bash

## new script that makes better use of the slurm batch system available on the heplx
#################################################################################
# ## NOTE: only the bare minimum setup is done here and then passed further on  #
# ##       All other setup is done in batchRunDataFits.sh                       #
# ##       Also the setup of the target directory is done by another script     #
# ##       checks if the source files are there (but not their timesptamps) so  #
# ##       make sure that the most recent source files are there (or none)      #
#################################################################################

## ONLY ONE OF THESE FLAGS CAN BE SET AT A TIME
FITTING=1 ## execute this first and when all batch jobs have finished execute merging
MERGING=0 ## although in principle this should check if all jobs are completed, proper working is not guaranteed!

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

FracLSB=0

## CHIC DATA
nState=7
# DataID=SetOfCuts11_chic_21June2016_rejCow_pt10Cut_chic$[$nState-5]Binning_corrLib_fLSB_${FracLSB}
# JobID=chic$[$nState-5]_21June2016_rejCow_cutDiMu10_chic$[$nState-5]Binning_corrLib_fLSB_${FracLSB}
DataID=SetOfCuts11_chic_21June2016_rejCow_pt10Cut_chic$[$nState-5]Binning_corrLib
JobID=chic$[$nState-5]_21June2016_rejCow_cutDiMu10_chic$[$nState-5]Binning_corrLib_statVarEff
TreeID=chic$[$nState-5] # switch to Psi1S for MC closure

## CHIC MC (STILL TODO: RUN)
# nState=6
# DataID=SetOfCuts11_chic_21June2016_MCclosure_rejCow_cutDiMu10_chic1Binning_corrLib
# JobID=chic_21June2016_MCclosure_rejCow_chic1Binning_corrLib
# TreeID=chic1

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
ptBinMax=5

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
      subtext=$(sbatch --array=$[$nSkipGen+1]-${nFits} ${storagedir}/${JobID}/batchRunDataFits.sh ${iPT} ${iRap} ${nSample} ${nFits} ${JobID} ${storagedir} ${basedir} ${FidCuts} ${TreeID} ${basedir}${datadir}/tmpFiles ${nSkipGen})
      echo ${subtext}

      batchJobID=${subtext##* }
      echo ${iRap} ${iPT} ${batchJobID} >> ${batchIDFile} ## store all the batch job IDs in a file for easier checking afterwards

      ## testing without the batch system
      # bash batchRunDataFits.sh ${iPT} ${iRap} ${nSample} ${nFits} ${JobID} ${storagedir} ${basedir} ${FidCuts} ${TreeID} ${basedir}${datadir}/tmpFiles ${nSkipGen} ${nState}

      echo "Submitted batch job for "$[$nFits-$nSkipGen]" fits for (pt,rap) bin ("${iPT}","${iRap}"). The job id is: "${batchJobID}

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
    jobState=$(sacct -j ${batchJobID} -n --format=state | uniq) ## if all wen't well this is only one entry
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

  done < ${batchIDFile}

  if [ -s ${batchIDFile} ]; then
    rm ${batchIDFile}
  fi
fi
