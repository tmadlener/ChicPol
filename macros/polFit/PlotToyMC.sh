#!/bin/sh

homedir=$PWD
cd ..
cd ..
basedir=$PWD
cd macros/polFit
storagedir=/afs/cern.ch/work/k/knuenz/storage/ChicPol #please define the directory storagedir in the file macros/polFit/storagedir

########## INPUTS ##########

for nState in 6; do

cp ../../interface/rootIncludes.inc               rootIncludes.inc
cp ../../interface/commonVar.h    commonVar.h
cp ../../interface/ToyMC.h        ToyMC.h
cp ../../interface/effsAndCuts_chi.h  effsAndCuts.h
touch polRapPtPlot.cc
make

for JobID in ToyMC_Test_20160530; do

echo ${JobID}


if [ $nState -eq 6 ]
then
ptBinMin=5
ptBinMax=5
fi
if [ $nState -eq 7 ]
then
ptBinMin=5
ptBinMax=5
fi

frameSig=3
for polScenSig in 3;do

frameBkg=3
for polScenBkg in 3;do

nGenerations=25

MPValgo=3 		#1...mean,2...gauss,3...gauss-loop with chi2<2
additionalName=MPV${MPValgo}

############################

TreeID=Chic$[nState-5]

cd ${basedir}/macros/polFit

rapBinMin=1 #don't change
rapBinMax=1 #don't change

ScenDir=Sig_frame${frameSig}scen${polScenSig}_Bkg_frame${frameBkg}scen${polScenBkg}

mkdir ${basedir}/macros/polFit/FiguresToyMC
mkdir ${basedir}/macros/polFit/FiguresToyMC/${JobID}
mkdir ${basedir}/macros/polFit/FiguresToyMC/${JobID}/${ScenDir}

cd ${storagedir}/${JobID}
mkdir ${ScenDir}

cp ${basedir}/macros/polFit/polRapPtPlot .

./polRapPtPlot ${ptBinMin}ptBinMin ${ptBinMax}ptBinMax ${rapBinMin}rapBinMin ${rapBinMax}rapBinMax ${frameSig}frameSig ${polScenSig}polScen ${MPValgo}MPValgo ${nGenerations}nGenerations ${ScenDir}=dirstruct ${nState}nState

mv ${ScenDir}/TGraphResults_${TreeID}_temp.root ${ScenDir}/TGraphResults_${TreeID}.root 

cp ${basedir}/latex/PullSummaryResults.tex ${ScenDir}/PullSummaryResults_${ScenDir}.tex
cp ${basedir}/latex/ParameterSummaryResults.tex ${ScenDir}/ParameterSummaryResults_${ScenDir}.tex
cp ${basedir}/latex/ToyResults.tex ${ScenDir}/ToyResults_${ScenDir}.tex

pdflatex ToyNumericalResults_${ScenDir}.tex
mv ToyNumericalResults_${ScenDir}.pdf ${basedir}/macros/polFit/FiguresToyMC/${JobID}/${ScenDir}/ToyNumericalResults_${ScenDir}_${additionalName}.pdf
rm *.aux
rm *.log

cd ${ScenDir}
pdflatex PullSummaryResults_${ScenDir}.tex
pdflatex ParameterSummaryResults_${ScenDir}.tex

rap_=${rapBinMin}
while [ $rap_ -le ${rapBinMax} ]
do
pT_=${ptBinMin}
while [ $pT_ -le ${ptBinMax} ]
do

pdflatex "\newcommand\rappt{rap${rap_}pt${pT_}}\input{ToyResults_${ScenDir}.tex}"
mv ToyResults_${ScenDir}.pdf ${basedir}/macros/polFit/FiguresToyMC/${JobID}/${ScenDir}/ToyResults_${ScenDir}_rap${rap_}pt${pT_}_${additionalName}.pdf

pT_=$[pT_+1]
done
rap_=$[rap_+1]
done

mv PullSummaryResults_${ScenDir}.pdf ${basedir}/macros/polFit/FiguresToyMC/${JobID}/${ScenDir}/PullSummaryResults_${ScenDir}_${additionalName}.pdf
mv ParameterSummaryResults_${ScenDir}.pdf ${basedir}/macros/polFit/FiguresToyMC/${JobID}/${ScenDir}/ParameterSummaryResults_${ScenDir}_${additionalName}.pdf

rm *.aux
rm *.log
rm *.tex

cd ..
rm polRapPtPlot

done
done
done

cd ${basedir}/macros/polFit
rm polRapPtPlot
done
