#include <string>
#include <vector>
#include <iostream>
#include <map>
#include <locale>
#include <algorithm>

/** map containting all the possible graph objects to plot. */
const std::map<std::string, std::string> createLabelMap()
{
  std::map<std::string, std::string> graphObjects;
  graphObjects["ltilde"] = "#tilde{#lambda}";
  graphObjects["lth"] = "#lambda_{#theta}";
  graphObjects["tph"] = "#lambda_{#phi}";
  graphObjects["ltp"] = "#lambda_{#theta#phi}";

  return graphObjects;
}
/*const*/ std::map<std::string, std::string> labelMap = createLabelMap();

/** map containing for different frame options. */
const std::map<std::string, std::vector<std::string> > createFrameMap()
{
  std::map<std::string, std::vector<std::string> > frameMap;
  frameMap["hx"] = { std::string("HX") };
  frameMap["px"] = { std::string("PX") };
  frameMap["cs"] = { std::string("CS") };
  frameMap["all"] = { std::string("HX"), std::string("PX"), std::string("CS") };

  return frameMap;
}
/*const*/ std::map<std::string, std::vector<std::string> > frameMap = createFrameMap();

/** put together the complete names of the graphs stored in the files.
 *NOTE: this can probably made constexpr in c++11.
 */
const std::vector<std::string> getGraphNames(const std::string& obj, const std::vector<std::string>& frames,
                                             const std::string& rapBin, const std::string& delim = "_")
{
  std::vector<std::string> completeNames;
  // for (const auto& frame : frames) { // no range-based for loops pre c++11
  for (size_t i = 0; i < frames.size(); ++i) {
    const std::string& frame = frames[i];
    completeNames.push_back(obj + delim + frame + delim + rapBin);
  }
  return completeNames;
}

/** helper struct to return the command line arguments in an ordered manner.
 * COULDDO: implement some sort of caching for the argument processing.
 */
struct cl_args {
  explicit cl_args(const std::vector<std::string>& argv, const std::string& rapBin);
  std::vector<std::string> ifn;
  std::string ofn;
  std::vector<std::string> graphs;
  std::string yAxisTitle;
  std::vector<std::string> legendEntries;
  std::vector<std::string> frames;
  std::string lambda;

private:
  const std::string getFileEnding(const std::string& arg)
  {
    size_t posLastDot = arg.find_last_of(".");
    if (posLastDot >= arg.length()) return std::string("");
    // std::cout << arg.substr(posLastDot, arg.length()) << std::endl;
    return arg.substr(posLastDot, arg.length());
  }

  bool isRootFile(const std::string& arg)
  {
    return getFileEnding(arg) == std::string(".ROOT");
  }

  bool isPdfFile(const std::string& arg)
  {
    return getFileEnding(arg) == std::string(".PDF");
  }

  bool isValidFrame(const std::string& arg)
  {
    return (arg == "CS" || arg == "HX" || arg == "PX" || arg == "ALL"); // hit one of the frames
  }

  bool isValidLambda(const std::string& arg)
  {
    return (arg == "LTILDE" || arg == "LTH" || arg == "LPH" || arg == "LTP"); // hit one of the lambdas
  }
};

cl_args::cl_args(const std::vector<std::string>& argv, const std::string& rapBin)
{
  for (size_t i = 0; i < argv.size(); ++i) {
    const std::string& arg = argv[i];
    std::string upper, lower;
    std::transform(arg.begin(), arg.end(), std::back_inserter(upper), ::toupper);
    std::transform(arg.begin(), arg.end(), std::back_inserter(lower), ::tolower);

    if (isRootFile(upper)) {
      ifn.push_back(arg);
    } else if (isPdfFile(upper)) {
      ofn.assign(arg);
    } else if (isValidFrame(upper)) {
      frames = frameMap[lower];
    } else if (isValidLambda(upper)) {
      lambda = lower;
      yAxisTitle = labelMap[lower];
    } else {
      legendEntries.push_back(arg);
    }
  }

  graphs = getGraphNames(lambda, frames, rapBin);

  // do some sort of basic checks and early catching
  if (legendEntries.size() < ifn.size()) {
    for (size_t i = legendEntries.size(); i < ifn.size(); ++i) {
      legendEntries.push_back(std::string(""));
    }
  }

}
