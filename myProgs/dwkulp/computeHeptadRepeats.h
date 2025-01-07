// computeHeptadRepeats.h
#ifndef COMPUTEHEPTADREPEATS_H
#define COMPUTEHEPTADREPEATS_H

#include <string>
#include <vector>
#include "System.h"
#include "CartesianPoint.h"
#include "Helanal.h"
#include "Line.h"

struct Options {
    Options(){
        required.push_back("pdbFile");
        optional.push_back("pymol");
        optional.push_back("plot");
    }
    std::string pdbFile;
    bool pymol;
    bool plot;
    std::vector<std::string> required;
    std::vector<std::string> optional;
};

Options setupOptions(int argc, char* argv[]);
void computeHeptadRepeats(const Options& opt);

#endif // COMPUTEHEPTADREPEATS_H