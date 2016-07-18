#include "muonTests.C"

#include <vector>
#include <string>

#include "/afs/hephy.at/user/t/tmadlener/snippets/vector_stuff.h"

/** return the pointer to the last element of a c-style array. */
template<typename T, size_t N>
T * end(T (&ra)[N]) {
    return ra + N;
}

/** ugly global vector, but at least in this way it can easily be changed, and it is in anyway only used once!
 * Leave this empty if you don't want to check for any triggers.
 */
const char* caTriggerNames[] = {"HLT_Dimuon8_Jpsi_v3", "HLT_Dimuon8_Jpsi_v4", "HLT_Dimuon8_Jpsi_v5",
                                "HLT_Dimuon8_Jpsi_v6", "HLT_Dimuon8_Jpsi_v7",
                                "HLT_Dimuon10_Jpsi_v3", "HLT_Dimuon10_Jpsi_v4", "HLT_Dimuon10_Jpsi_v5",
                                "HLT_Dimuon10_Jpsi_v6"};
const std::vector<std::string> triggerNames(caTriggerNames, end(caTriggerNames));
// const std::vector<std::string> triggerNames;

/** another global vector, storing all the possible tree-names to loop over them in a somewhat automated fashion. */
const char* caTreeNames[] = {"rootuple/chicTree", "data", "selectedData"};
std::vector<std::string> treeNames(caTreeNames, end(caTreeNames));

/** main. */
int main(int argc, char* const argv[])
{
  std::string ifn = argv[1];
  TFile* inFile = TFile::Open(ifn.c_str());
  TTree* inTree = NULL;

  for (size_t i = 0; i < treeNames.size(); ++i) {
    inTree = static_cast<TTree*>(inFile->Get(treeNames[i].c_str()));
    if (inTree) {
      std::cout << "Found TTree \'" << treeNames[i] << "\' in file \'" << ifn << "\'." << std::endl;
      break;
    }
  }
  if (!inTree) {
    std::cerr << "Cannot find TTree in file: \'" << ifn << "\'" << std::endl;
    return -1;
  }

  muonTests t(triggerNames, inTree);
  t.Loop();

  return 0;
}
