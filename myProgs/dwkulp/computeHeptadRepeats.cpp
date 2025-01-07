// computeHeptadRepeats.cpp
#include "computeHeptadRepeats.h"
#include "OptionParser.h"
#include "PyMolVisualization.h"
#include "AtomSelection.h"
#include "MslTools.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;
using namespace MSL;

Options setupOptions(int argc, char* argv[]) {
    Options opt;

    OptionParser OP;
    OP.setRequired(opt.required);
    OP.setAllowed(opt.optional);
    OP.autoExtendOptions();
    OP.readArgv(argc, argv);

    if (OP.countOptions() == 0) {
        cout << "Usage:" << endl;
        cout << "  computeHeptadRepeats --pdbFile <PDB_FILE> [options]" << endl;
        cout << endl;
        cout << "Required options:" << endl;
        cout << "  --pdbFile <PDB_FILE>" << endl;
        cout << endl;
        cout << "Optional options:" << endl;
        cout << "  --pymol" << endl;
        cout << endl;
        exit(0);
    }

    opt.pdbFile = OP.getString("pdbFile");
    if (OP.fail()) {
        cerr << "ERROR: pdbFile not specified.\n";
        exit(1);
    }

    opt.pymol = OP.getBool("pymol");
    if (OP.fail()) {
        opt.pymol = false;
    }
    opt.plot = OP.getBool("plot");
    if (OP.fail()) {
        opt.plot = false;
    }

    return opt;
}

void computeHeptadRepeats(const Options& opt) {
    System sys;
    sys.readPdb(opt.pdbFile);

    vector<AtomPointerVector> chains;
    for (uint c = 0; c < sys.chainSize(); c++) {
        chains.push_back(sys.getChain(c).getAtomPointers());
    }

    if (chains.size() < 2) {
        cerr << "ERROR: At least two chains are required.\n";
        exit(1);
    }
    // Create sel1 and sel2 as CA atoms from chains
    AtomSelection as1(chains[0]);
    AtomPointerVector sel1 = as1.select("name CA");

    AtomSelection as2(chains[1]);
    AtomPointerVector sel2 = as2.select("name CA");

    // Find two sets of 7 that are close to each other
    map<string, double> distances;

    ofstream fout;
    int windowCount = 0;
    PyMolVisualization py;
    if (opt.pymol){
        fout.open("helanal.py");
    }

    for (uint i = 0; i < sel1.size() - 6; i += 7) {
        double minDistance = numeric_limits<double>::max();
        uint bestStart = 0;

        // Loop over all 7-residue windows in the second chain
        // To find the best match for the current window in the first chain
        for (uint j = 0; j < sel2.size() - 6; j++) {
            double currentDistance = 0.0;
            for (uint k = 0; k < 7; k++) {
                currentDistance += sel1(i + k).distance(sel2(j + k));
            }
            if (currentDistance < minDistance) {
                minDistance = currentDistance;
                bestStart = j;
            }
        }
        // Compute helanal for the two sets of 7 residues
        Helanal h;
        CartesianPoint helanalAxisA(0.0, 0.0, 0.0);
        CartesianPoint helanalCenterA(0.0, 0.0, 0.0);
        CartesianPoint helanalAxisB(0.0, 0.0, 0.0);
        CartesianPoint helanalCenterB(0.0, 0.0, 0.0);

        for (uint n = 0; n < 7; n++) {
            h.update(sel1(i + n).getCoor(), sel1(i + n + 1).getCoor(), sel1(i + n + 2).getCoor(), sel1(i + n + 3).getCoor());
            helanalAxisA += h.getAxis();
            helanalCenterA += h.getCenter();

            h.update(sel2(bestStart + n).getCoor(), sel2(bestStart + n + 1).getCoor(), sel2(bestStart + n + 2).getCoor(), sel2(bestStart + n + 3).getCoor());
            helanalAxisB += h.getAxis();
            helanalCenterB += h.getCenter();
        }

        CartesianPoint centerA = helanalCenterA / 7;
        CartesianPoint centerB = helanalCenterB / 7;
        CartesianPoint dir1 = helanalAxisA.getUnit();
        CartesianPoint dir2 = helanalAxisB.getUnit();

        Line helixA(centerA, dir1);
        Line helixB(centerB, dir2);

        CartesianPoint bundleAxis = (dir1 + dir2) / 2;
        CartesianPoint bundleMidpoint = (centerA + centerB) / 2;
        Line bundleZ(bundleMidpoint, bundleAxis);

        // Compute distances from each residue to the bundle axis
        for (uint n = 0; n < 7; n++) {
            CartesianPoint projectionA = bundleZ.projection(sel1(i + n).getCoor());
            double distanceA = projectionA.distance(sel1(i + n).getCoor());
            string posIdA = sel1(i + n).getPositionId();
            distances[posIdA] = distanceA;

            CartesianPoint projectionB = bundleZ.projection(sel2(bestStart + n).getCoor());
            double distanceB = projectionB.distance(sel2(bestStart + n).getCoor());
            string posIdB = sel2(bestStart + n).getPositionId();
            distances[posIdB] = distanceB;
        }


        if (opt.pymol) {
            windowCount++;
            
            py.createArrow(centerA, dir1, MslTools::stringf("%d-HelA", windowCount));
            py.createArrow(centerB, dir2, MslTools::stringf("%d-HelB", windowCount));
            py.createArrow(bundleMidpoint, bundleAxis, MslTools::stringf("%d-BundleZ",windowCount));
            CartesianPoint mutA = helixA.pointOfMinDistanceToLine(helixB);
            CartesianPoint mutB = helixB.pointOfMinDistanceToLine(helixA);
            CartesianPoint mutMid = (mutA + mutB) / 2;
            Line toMut(bundleMidpoint, mutMid - bundleMidpoint);
            CartesianPoint c = toMut.getCenter();
            CartesianPoint d = toMut.getDirection();
            py.createArrow(c, d, MslTools::stringf("%d-toMut",windowCount));

            fout << py << endl;
        }

    }
    fout.close();

    string datafile =MslTools::stringf("%s-output.csv", MslTools::getFileName(opt.pdbFile).c_str());
    ofstream outFile(datafile.c_str());
    outFile << "Chain,ResNum,ResId,Distance" << endl;
    for (const auto& entry : distances) {
        string posId = entry.first;
        double distance = entry.second;
        // get a Position from sys using posId
        Position &pos = sys.getPosition(posId);
        string chainId = pos.getChainId();
        int resNum = pos.getResidueNumber();
        string resId = pos.getResidueName();
        outFile << chainId << "," << resNum << "," << resId << "," << distance << endl;
    }

    outFile.close();

    // Iterate over positions, find 14 residue windows that have two heptads. define heptad when residue 2 and 5 are closeset to the helix axis out of the 7 residues.
    // use distances object to look up the distances for each position
    fprintf(stdout, "%-25s, %4s seq: %s\n", "Structure", "ResN", "gAbcDefgAbcDef");
    for (uint i = 0; i < sys.getChain(0).positionSize();i++){
        // find 14 residue windows that have two heptads. define heptad when residue 2 and 5 are closeset to the helix axis out of the 7 residues.
        int j = i;
        if (j+14>=sys.getChain(0).positionSize()) break;
    
        bool isHeptad1 = (
            distances[sys.getChain(0).getPosition(j + 1).getPositionId()] < distances[sys.getChain(0).getPosition(j).getPositionId()] &&
            distances[sys.getChain(0).getPosition(j + 1).getPositionId()] < distances[sys.getChain(0).getPosition(j + 2).getPositionId()] &&
            distances[sys.getChain(0).getPosition(j + 4).getPositionId()] < distances[sys.getChain(0).getPosition(j + 3).getPositionId()] &&
            distances[sys.getChain(0).getPosition(j + 4).getPositionId()] < distances[sys.getChain(0).getPosition(j + 5).getPositionId()] &&
            distances[sys.getChain(0).getPosition(j + 1).getPositionId()] < 5.6 &&
            distances[sys.getChain(0).getPosition(j + 4).getPositionId()] < 5.6 &&
            sys.getChain(0).getPosition(j + 1).getResidueName() != "ARG" &&
            sys.getChain(0).getPosition(j + 1).getResidueName() != "LYS" &&
            sys.getChain(0).getPosition(j + 1).getResidueName() != "HIS" &&
            sys.getChain(0).getPosition(j + 1).getResidueName() != "GLU" &&
            sys.getChain(0).getPosition(j + 4).getResidueName() != "ARG" &&
            sys.getChain(0).getPosition(j + 4).getResidueName() != "LYS" &&
            sys.getChain(0).getPosition(j + 4).getResidueName() != "HIS" &&
            sys.getChain(0).getPosition(j + 4).getResidueName() != "GLU");

        bool isHeptad2 = (
            distances[sys.getChain(0).getPosition(j + 8).getPositionId()] < distances[sys.getChain(0).getPosition(j + 7).getPositionId()] &&
            distances[sys.getChain(0).getPosition(j + 8).getPositionId()] < distances[sys.getChain(0).getPosition(j + 9).getPositionId()] &&
            distances[sys.getChain(0).getPosition(j + 11).getPositionId()] < distances[sys.getChain(0).getPosition(j + 10).getPositionId()] &&
            distances[sys.getChain(0).getPosition(j + 11).getPositionId()] < distances[sys.getChain(0).getPosition(j + 12).getPositionId()] &&
            distances[sys.getChain(0).getPosition(j + 8).getPositionId()] < 5.6 &&
            distances[sys.getChain(0).getPosition(j + 11).getPositionId()] < 5.6 &&
            sys.getChain(0).getPosition(j + 8).getResidueName() != "ARG" &&
            sys.getChain(0).getPosition(j + 8).getResidueName() != "LYS" &&
            sys.getChain(0).getPosition(j + 8).getResidueName() != "HIS" &&
            sys.getChain(0).getPosition(j + 8).getResidueName() != "GLU" &&
            sys.getChain(0).getPosition(j + 11).getResidueName() != "ARG" &&
            sys.getChain(0).getPosition(j + 11).getResidueName() != "LYS" &&
            sys.getChain(0).getPosition(j + 11).getResidueName() != "HIS" &&
            sys.getChain(0).getPosition(j + 11).getResidueName() != "GLU");

        if (isHeptad1 && isHeptad2) {
            //  get 14 amino acid sequence starting at i
            string seq = "";
            for (uint k = 0; k < 14; k++) {
                seq += MslTools::getOneLetterCode(sys.getChain(0).getPosition(j + k).getResidueName());
            }
            //cout << "Found two heptads in a 14-residue window starting at position " << sys.getChain(0).getPosition(j).getResidueNumber() << " seq: "<<seq<<endl;
            // formated print
            fprintf(stdout, "%-25s, %4d seq: %s\n", MslTools::getFileName(opt.pdbFile).c_str(), sys.getChain(0).getPosition(j).getResidueNumber(), seq.c_str());
        }
    }




    // Plot using ggplot2 distance vs position, color by amino acid type in R
    // if opt.plot
    if (opt.plot) {

        cout << "Plotting using ggplot2..." << endl;
        // Write R script to plot
        string plotscript = MslTools::stringf("%s-plot.r", MslTools::getFileName(opt.pdbFile).c_str());
        string helix_analysis_out = MslTools::stringf("helix_analysis_%s.pdf", MslTools::getFileName(opt.pdbFile).c_str());
        ofstream rScript(plotscript);
        rScript << "library(ggplot2)" << endl;
        rScript << "data <- read.csv(\"" << datafile << "\", header=T,sep=\",\")" << endl;
        rScript << "ggplot(data, aes(x=ResNum, y=Distance, color=ResId, group=Chain)) + geom_point() + geom_line() + theme_minimal() + facet_wrap(~Chain)" << endl;
        rScript << "ggsave(\""<<helix_analysis_out<< "\")" << endl;
        
                
        string helix_analysis_out_cat = MslTools::stringf("helix_analysis_category_%s.pdf", MslTools::getFileName(opt.pdbFile).c_str());
        rScript << "# Define categories with uppercase three-letter codes" << endl;
        rScript << "hydrophobic <- c(\"ALA\", \"VAL\", \"ILE\", \"LEU\", \"MET\", \"MSE\", \"PHE\", \"TRP\", \"TYR\")" << endl;
        rScript << "hydrophilic <- c(\"SER\", \"THR\", \"ASN\", \"GLN\", \"CYS\", \"GLY\", \"PRO\")" << endl;
        rScript << "charged <- c(\"ASP\", \"GLU\", \"LYS\", \"ARG\", \"HIS\")" << endl;
        rScript << "" << endl;
        rScript << "# Assign category" << endl;
        rScript << "data$category <- with(data, ifelse(ResId %in% hydrophobic, \"Hydrophobic\"," << endl;
        rScript << "                            ifelse(ResId %in% hydrophilic, \"Hydrophilic\", \"Charged\")))" << endl;
        rScript << "" << endl;
        rScript << "data.sel <- data[data$ResNum < 140,];" << endl;
        rScript << "" << endl;
        rScript << "ggplot(data.sel, aes(x=ResNum, y=Distance, color=category, group=Chain)) + geom_point() + geom_line() + scale_color_manual(values = c(\"Hydrophobic\" = \"blue\", \"Hydrophilic\" = \"green\", \"Charged\" = \"red\")) + theme_minimal() + facet_wrap(~Chain,nrow=2)" << endl;
        rScript << "ggsave(\""<<helix_analysis_out_cat<< "\")" << endl;
        rScript.close();
        // Run R script
        int retVal = system(MslTools::stringf("Rscript %s",plotscript.c_str()).c_str());
        
    }


 
}

int main(int argc, char* argv[]) {
    Options opt = setupOptions(argc, argv);
    computeHeptadRepeats(opt);
    return 0;
}