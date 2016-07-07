#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

#include "../interface/rootIncludes.inc" // hacky solution!

#include "../macros/bkgHistos_helper.h"
// #include "snippets/string_stuff.h"
#include "snippets/vector_stuff.h"

std::vector<std::string> readLines(const std::string& fname)
{
  std::ifstream ifs(fname);
  std::vector<std::string> lines;
  std::string line;
  while (std::getline(ifs, line)) {
    if (!startsWith(line, "#")) lines.push_back(line);
  }
  ifs.close();

  return lines;
}

int main(int argc, char* argv[])
{
  if (argc < 3) {
    std::cerr << "need a .root file and two files of names (1 for variables, 1 for pdf) to print!" << std::endl;
    return -1;
  }

  std::ofstream outfileStream;
  bool writeToFile = false;
  if (argc == 5) {
    outfileStream.open(argv[4], std::ios::out);
    writeToFile = true;
  }

  TFile* rf = TFile::Open(argv[1], "READ");
  RooWorkspace* ws = static_cast<RooWorkspace*>(rf->Get("ws_masslifetime"));

  // extract pt from root file name
  const std::vector<std::string> fns = splitString(argv[1], '/');
  const char* ptBin = fns.back().substr(fns.back().length()-6, 1).c_str(); // NOTE: very dirty!!

  std::stringstream snapName;
  snapName << "m_snapshot_rap1_pt" << ptBin;
  ws->loadSnapshot(snapName.str().c_str());

  std::vector<std::string> printVars = readLines(argv[2]);

  std::ostream & ofs = writeToFile ? outfileStream : std::cout;
  for (size_t i = 0; i < printVars.size(); ++i) {
    ofs << printVars[i] << " = " << getVarVal(ws, printVars[i]) << " +/- " << getVarError(ws, printVars[i]) << std::endl;
  }

  std::vector<std::string> printPdfs = readLines(argv[3]);
  for (size_t i = 0; i < printPdfs.size(); ++i) {
    RooAbsPdf* pdf = getPdf(ws, printPdfs[i]);
    pdf->Print();
    ofs << printPdfs[i] << ": ";
    pdf->printStream(ofs, 4, RooPrintable::StyleOption::kInline);
    ofs << std::endl;
  }

  ofs.flush();

  return 0;
}
