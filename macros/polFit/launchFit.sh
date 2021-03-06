#!/bin/bash

subJob=1
batchType='bsub -q 1nd' #'qsub -q cmsq'

workPath=$PWD
nSigma=3.00
#originFile=runDataFitsJob
originFile=runToyMCJob

#jobType=noRhoFactor
jobType=ToyMC
batchDir=${workPath}/batchJob/${jobType}
mkdir -p ${batchDir}
mkdir -p ${batchDir}/logfile

cp ${originFile}.sh ${batchDir}/${originFile}.sh
cd ${batchDir}

for nState in 4;do

##############
##### 1S #####
##############
if [ $nState -eq 4 ]
then
for rapBin in 1 2; do 
for ptBin in 1 2 3 4 5 6 7 8 9 10 11 12; do 
for nSkipGen in 0; do
#for nSkipGen in 0 1 2 3 4 5 6 7 8 9; do
JobID=Psi$[nState-3]S_${nSigma}Sigma_12Dec2012
cp ${originFile}.sh ${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
chmod 755 ${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
sed "s:HOMEDIR:${workPath}:g"   -i ${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
sed "s:NState:${nState}:g"      -i ${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
sed "s:NSkipGen:${nSkipGen}:g"  -i ${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
sed "s:RapBinMin:${rapBin}:g"   -i ${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
sed "s:RapBinMax:${rapBin}:g"   -i ${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
sed "s:PtBinMin:${ptBin}:g"     -i ${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
sed "s:PtBinMax:${ptBin}:g"     -i ${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh

echo "#!/bin/bash 
source /afs/cern.ch/user/z/zhlinl/rootset.sh
cd ${batchDir}
./${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh > logfile/log_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}SkipGen
rm ${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
rm job_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
" > job_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
chmod 755 job_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh

if [ $subJob -eq 1 ]
then
${batchType} job_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh  ## bsub: -q 1nh 1nh80 1nd 2nd 1nw
fi
done
done
done
fi

##############
##### 2S #####
##############
if [ $nState -eq 5 ]
then
for rapBin in 1 2 3; do 
for ptBin in 2 3 4 5 6; do 
for nSkipGen in 0; do
#for nSkipGen in 0 1 2 3 4 5 6 7 8 9; do
cp ${originFile}.sh ${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
chmod 755 ${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
sed "s:HOMEDIR:${workPath}:g"   -i ${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
sed "s:NState:${nState}:g"      -i ${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
sed "s:NSkipGen:${nSkipGen}:g"  -i ${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
sed "s:RapBinMin:${rapBin}:g"   -i ${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
sed "s:RapBinMax:${rapBin}:g"   -i ${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
sed "s:PtBinMin:${ptBin}:g"     -i ${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
sed "s:PtBinMax:${ptBin}:g"     -i ${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh

echo "#!/bin/bash 
source /afs/cern.ch/user/z/zhlinl/rootset.sh
cd ${batchDir}
./${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh > logfile/log_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}SkipGen
rm ${originFile}_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
rm job_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
" > job_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh
chmod 755 job_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh

if [ $subJob -eq 1 ]
then
${batchType} job_Psi$[nState-3]S_${rapBin}rap_${ptBin}pt_${nSkipGen}.sh  ## bsub: -q 1nh 1nh80 1nd 2nd 1nw
fi
done
done
done
fi

done
