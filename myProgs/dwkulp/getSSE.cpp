#include <iostream>
#include <cstdlib>
#include <queue>
#include <algorithm>
#include <map>

#include "System.h"
#include "OptionParser.h"
#include "MslOut.h"
#include "SysEnv.h"
#include "SasaCalculator.h"
#include "PhiPsiStatistics.h"
#include "MonteCarloManager.h"
#include "Quench.h"
#include "PDBTopology.h"
#include "AtomSelection.h"
#include "PyMolVisualization.h"
#include "RandomNumberGenerator.h"
#include "getSSE.h"

using namespace std;
using namespace MSL;

static MslOut MSLOUT("getSSE");
static SysEnv SYSENV;


int main(int argc, char *argv[]) {

  // Parse commandline options
  Options opt = setupOptions(argc, argv);

  // Read PDB structure
  System sys;
  sys.readPdb(opt.pdb);

  // Read MSA file if provided
  vector<map<string, int>> aaDiversity;
  if (!opt.msa.empty()) {
    aaDiversity = computeAADiversity(sys, opt.msa, opt.pdb);
  }

  SasaCalculator scalc(sys.getAtomPointers());
  scalc.calcSasa();

  map<string, int> sse_map;
  vector<pair<string, pair<int, int>>> sse_pos;
  stringstream allSSE;
  stringstream previousSSE;
  string lastSSE = "";

  // Compute SASA for each residue and store per SSE  
  for (uint p = 1; p < sys.positionSize() - 1; p++) {
    if (!sys.getPosition(p).atomExists("CA")) continue;
    // if (!sys.getPosition(p).getAtom("CA").getSelectionFlag("resurface")) continue;

    double phi = PhiPsiStatistics::getPhi(sys.getPosition(p - 1).getCurrentIdentity(), sys.getPosition(p).getCurrentIdentity());
    double psi = PhiPsiStatistics::getPsi(sys.getPosition(p).getCurrentIdentity(), sys.getPosition(p + 1).getCurrentIdentity());

    string sse = "";
    if (phi < -35 && phi > -90 && psi > -70 && psi < 0) {
      sse = "H";
    } else {
      bool isSheet = checkBetaSheet(sys, p);
      sse = isSheet ? "E" : "O";
    }

    // SSE Boundary
    if (sse != lastSSE) {
      if (previousSSE.str().length() >= 6) {
        sse_map[lastSSE]++;
        int p_first = p - previousSSE.str().length();

        pair<int, int> pInt = std::make_pair(p_first, p);
        pair<string, pair<int, int>> pStrInts = std::make_pair(lastSSE, pInt);
        sse_pos.push_back(pStrInts);

        fprintf(stdout, "Secondary Structure Element is: %s, residues: %d-%d\n", lastSSE.c_str(), sys.getPosition(p_first).getResidueNumber(), sys.getPosition(p).getResidueNumber());
      }
      previousSSE.str("");
    }

    lastSSE = sse;
    previousSSE << sse;
    allSSE << sse;
  }

  if (previousSSE.str().length() >= 6) {
    sse_map[lastSSE]++;
    int p_first = sys.positionSize() - 1 - previousSSE.str().length();

    pair<int, int> pInt = std::make_pair(p_first, sys.positionSize() - 1);
    pair<string, pair<int, int>> pStrInts = std::make_pair(lastSSE, pInt);
    sse_pos.push_back(pStrInts);

    fprintf(stdout, "Secondary Structure Element is: %s, residues: %d-%d\n", lastSSE.c_str(), sys.getPosition(p_first).getResidueNumber(), sys.getPosition(sys.positionSize() - 1).getResidueNumber());
  }

  // Need to add a - for last Amino Acid
  allSSE << "-";

  fprintf(stdout, "FULL SSE: -%s\n", allSSE.str().c_str());
  map<string, int>::iterator it;
  for (it = sse_map.begin(); it != sse_map.end(); it++) {
    fprintf(stdout, "SSE %s %d\n", it->first.c_str(), it->second);
  }

  // Store SASA data about each SSE if option is set
  vector<sasaData> sasa_pos;

  if (opt.sasa) {
    for (uint i = 0; i < sse_pos.size(); i++) { // for each SSE
      double sse_sasa = 0.0;
      int exposed_positions = 0;
      int total_positions = 0;
      for (int sse1 = sse_pos[i].second.first; sse1 <= sse_pos[i].second.second; sse1++) { // for each position in SSE
        total_positions += 1;

        double normSasa = scalc.getResidueSasa(sys.getPosition(sse1).getPositionId()) / refSasa[MslTools::getOneLetterCode(sys.getPosition(sse1).getResidueName())];
        if (normSasa > 1.0) {
          normSasa = 1.0;
        }
        sse_sasa += normSasa;

        if (normSasa > 0.3) {
          exposed_positions += 1;
        }
        fprintf(stdout, "SSE SASA %2d %8s %6.2f %6.2f\n", i, sys.getPosition(sse1).getPositionId().c_str(), normSasa, sse_sasa);
      }
      sasaData sd;
      sd.total_pos = total_positions;
      sd.exposed_pos = exposed_positions;
      sd.total_nSasa = sse_sasa;
      sasa_pos.push_back(sd);
    }
  }

  // Check for interactions between SSE X and Y
  map<int, vector<int>> interactingSSEs;
  vector<int> listSSEs;
  for (uint i = 0; i < sse_pos.size(); i++) {
    if (sse_pos[i].first == "O" && !opt.include_loops) continue; // skip loops

    for (uint j = i + 1; j < sse_pos.size(); j++) {
      if (sse_pos[j].first == "O" && !opt.include_loops) continue; // skip loops

      int numContactingResidues = 0;
      // Find neighbors between every pair of positions of these two SSEs
      for (int sse1 = sse_pos[i].second.first; sse1 <= sse_pos[i].second.second; sse1++) {
        for (int sse2 = sse_pos[j].second.first; sse2 <= sse_pos[j].second.second; sse2++) {
          bool neighbor = sys.getPosition(sse1).getCurrentIdentity().isNeighbor(opt.neighbor_dist, sys.getPosition(sse2).getCurrentIdentity());
          if (neighbor) {
            numContactingResidues++;
          }
        }
      }

      if (numContactingResidues > opt.num_contacting) {
        fprintf(stdout, "INTERACTING SSEs: %s-%s %d %d %d %d\n", sse_pos[i].first.c_str(), sse_pos[j].first.c_str(), sse_pos[i].second.first, sse_pos[i].second.second, sse_pos[j].second.first, sse_pos[j].second.second);

        map<int, vector<int>>::iterator aIt = interactingSSEs.find(i);
        if (aIt != interactingSSEs.end()) {
          aIt->second.push_back(j); // might add twice?
        } else {
          vector<int> tmp;
          tmp.push_back(j);
          interactingSSEs[i] = tmp;
        }
        aIt = interactingSSEs.find(j);
        if (aIt != interactingSSEs.end()) {
          aIt->second.push_back(i);
        } else {
          vector<int> tmp;
          tmp.push_back(i);
          interactingSSEs[j] = tmp;
        }
        if (std::find(listSSEs.begin(), listSSEs.end(), i) == listSSEs.end()) {
          listSSEs.push_back(i);
        }
        if (std::find(listSSEs.begin(), listSSEs.end(), j) == listSSEs.end()) {
          listSSEs.push_back(j);
        }
      }
    }
  }

  map<int, vector<int>>::iterator aIt;
  for (aIt = interactingSSEs.begin(); aIt != interactingSSEs.end(); aIt++) {
    fprintf(stdout, "%d[%d-%d]: ", aIt->first, sse_pos[aIt->first].second.first, sse_pos[aIt->first].second.second);
    for (int b = 0; b < aIt->second.size(); b++) {
      fprintf(stdout, " %d[%d-%d]", aIt->second[b], sse_pos[aIt->second[b]].second.first, sse_pos[aIt->second[b]].second.second);
    }
    fprintf(stdout, "\n");
  }

  // Write out all pairs of interacting SSEs that are exposed
  if (opt.sasa) {
    map<int, vector<int>>::iterator aIt;
    for (aIt = interactingSSEs.begin(); aIt != interactingSSEs.end(); aIt++) {
      for (int b = 0; b < aIt->second.size(); b++) {
        if (sasa_pos[aIt->first].exposed_pos >= opt.exposed_pos_threshold && sasa_pos[aIt->second[b]].exposed_pos >= opt.exposed_pos_threshold) {
            int p1 = sse_pos[aIt->first].second.first;
            int p2 = sse_pos[aIt->first].second.second;
            int p3 = sse_pos[aIt->second[b]].second.first;
            int p4 = sse_pos[aIt->second[b]].second.second;
            fprintf(stdout, "EXPOSED SSEs: %3d[%4d-%-4d or %8s-%8s] ; %2d, %6.2f + ", aIt->first, p1, p2, sys.getPosition(p1).getPositionId().c_str(), sys.getPosition(p2).getPositionId().c_str(), sasa_pos[aIt->first].exposed_pos, sasa_pos[aIt->first].total_nSasa);
            fprintf(stdout, "%3d[%4d-%-4d or %8s-%-8s] ; %2d, %6.2f \n", aIt->second[b], p3, p4, sys.getPosition(p3).getPositionId().c_str(), sys.getPosition(p4).getPositionId().c_str(), sasa_pos[aIt->second[b]].exposed_pos, sasa_pos[aIt->second[b]].total_nSasa);

            // Write PDB file for the pair of SSEs
            AtomContainer ats;
            for (int pos = p1; pos <= p2; pos++) {
                ats.addAtoms(sys.getPosition(pos).getAtomPointers());
            }
            for (int pos = p3; pos <= p4; pos++) {
                ats.addAtoms(sys.getPosition(pos).getAtomPointers());
            }
            stringstream ss;
            ss << MslTools::getFileName(opt.pdb) << "_SASA_SSEs_" << aIt->first << "_" << aIt->second[b] << ".pdb";
            ats.writePdb(ss.str());
        }
      }
    }
  }
  if (opt.msa != "") {
    for (auto &interaction : interactingSSEs) {
      int sse1 = interaction.first;
      for (int sse2 : interaction.second) {
        double avgFreqSSE1 = 0.0;
        double avgFreqSSE2 = 0.0;
        int countSSE1 = 0;
        int countSSE2 = 0;

        for (int pos = sse_pos[sse1].second.first; pos <= sse_pos[sse1].second.second; pos++) {
          string aa = MslTools::getOneLetterCode(sys.getPosition(pos).getResidueName());
          avgFreqSSE1 += ((aaDiversity[pos][aa] / sys.positionSize())*100);
          countSSE1++;
        }
        for (int pos = sse_pos[sse2].second.first; pos <= sse_pos[sse2].second.second; pos++) {
          string aa = MslTools::getOneLetterCode(sys.getPosition(pos).getResidueName());
          avgFreqSSE2 += ((aaDiversity[pos][aa] / sys.positionSize()) * 100);
          countSSE2++;
        }

        avgFreqSSE1 /= countSSE1;
        avgFreqSSE2 /= countSSE2;

        if (avgFreqSSE1 > opt.aa_freq_threshold && avgFreqSSE2 > opt.aa_freq_threshold) {
            fprintf(stdout, "INTERACTING SSEs with high AA frequency: %-2s-%-2s %4d-%-4d %4d-%-4d\n", 
                sse_pos[sse1].first.c_str(), sse_pos[sse2].first.c_str(), 
                sse_pos[sse1].second.first, sse_pos[sse1].second.second, 
                sse_pos[sse2].second.first, sse_pos[sse2].second.second);

            // Write PDB file for the pair of SSEs
            AtomContainer ats;
            for (int pos = sse_pos[sse1].second.first; pos <= sse_pos[sse1].second.second; pos++) {
                ats.addAtoms(sys.getPosition(pos).getAtomPointers());
            }
            for (int pos = sse_pos[sse2].second.first; pos <= sse_pos[sse2].second.second; pos++) {
                ats.addAtoms(sys.getPosition(pos).getAtomPointers());
            }
            stringstream ss;
            ss << MslTools::getFileName(opt.pdb) << "_MSA_SSEs_" << sse1 << "_" << sse2 << ".pdb";
            ats.writePdb(ss.str());
        }
      }
    }
  }

  if (opt.msa != "" && opt.sasa) {
    for (auto &interaction : interactingSSEs) {
      int sse1 = interaction.first;
      for (int sse2 : interaction.second) {
        double avgFreqSSE1 = 0.0;
        double avgFreqSSE2 = 0.0;
        int countSSE1 = 0;
        int countSSE2 = 0;

        for (int pos = sse_pos[sse1].second.first; pos <= sse_pos[sse1].second.second; pos++) {
          string aa = MslTools::getOneLetterCode(sys.getPosition(pos).getResidueName());
          avgFreqSSE1 += ((aaDiversity[pos][aa] / sys.positionSize()) * 100);
          countSSE1++;
        }
        for (int pos = sse_pos[sse2].second.first; pos <= sse_pos[sse2].second.second; pos++) {
          string aa = MslTools::getOneLetterCode(sys.getPosition(pos).getResidueName());
          avgFreqSSE2 += ((aaDiversity[pos][aa] / sys.positionSize()) * 100);
          countSSE2++;
        }

        avgFreqSSE1 /= countSSE1;
        avgFreqSSE2 /= countSSE2;

        if (avgFreqSSE1 > opt.aa_freq_threshold && avgFreqSSE2 > opt.aa_freq_threshold &&
            sasa_pos[sse1].exposed_pos > opt.exposed_pos_threshold && sasa_pos[sse2].exposed_pos > opt.exposed_pos_threshold) {
            fprintf(stdout, "INTERACTING SSEs with high AA frequency and exposure: %-2s-%-2s %4d-%-4d %4d-%-4d\n", 
                    sse_pos[sse1].first.c_str(), sse_pos[sse2].first.c_str(), 
                    sse_pos[sse1].second.first, sse_pos[sse1].second.second, 
                    sse_pos[sse2].second.first, sse_pos[sse2].second.second);

            // Write PDB file for the pair of SSEs
            AtomContainer ats;
            for (int pos = sse_pos[sse1].second.first; pos <= sse_pos[sse1].second.second; pos++) {
                ats.addAtoms(sys.getPosition(pos).getAtomPointers());
            }
            for (int pos = sse_pos[sse2].second.first; pos <= sse_pos[sse2].second.second; pos++) {
                ats.addAtoms(sys.getPosition(pos).getAtomPointers());
            }
            stringstream ss;
            ss << MslTools::getFileName(opt.pdb) << "_MSA_SASA_SSEs_" << sse1 << "_" << sse2 << ".pdb";
            ats.writePdb(ss.str());
        }
      }
    }
  }

   
  if (opt.pick_random_SSEs != 0) {
    map<string, bool> tripletUsedAlready;
    RandomNumberGenerator rng;
    for (uint n = 0; n < opt.pick_random_SSEs; n++) {
      int k = 0;
      int randomIndex1 = -1;
      int randomIndex2 = -1;
      int randomIndex3 = -1;

      while (k < 1000) {
        randomIndex1 = listSSEs[rng.getRandomInt(listSSEs.size() - 1)];
        randomIndex2 = interactingSSEs[randomIndex1][rng.getRandomInt(interactingSSEs[randomIndex1].size() - 1)];
        randomIndex3 = interactingSSEs[randomIndex2][rng.getRandomInt(interactingSSEs[randomIndex2].size() - 1)];
        vector<int> v = {randomIndex1, randomIndex2, randomIndex3};
        sort(v.begin(), v.end());
        stringstream key;
        key << v[0] << "-" << v[1] << "-" << v[2];
        if (randomIndex1 != randomIndex2 &&
          randomIndex2 != randomIndex3 &&
          randomIndex1 != randomIndex3 &&
          tripletUsedAlready.find(key.str()) == tripletUsedAlready.end()) {
          // fprintf(stdout,"RND: %d,%d,%d\n",randomIndex1,randomIndex2,randomIndex3);
          break;
        }
        k++;
      }
      if (k == 1000) {
        fprintf(stdout, "ERROR no random triplet of SSEs found after 1000 tries\n");
      } else {
        vector<int> v = {randomIndex1, randomIndex2, randomIndex3};
        sort(v.begin(), v.end());
        stringstream key;
        key << v[0] << "-" << v[1] << "-" << v[2];
        tripletUsedAlready[key.str()] = true;
        fprintf(stdout, "%d Random triplet: %d,%d,%d\n", k, randomIndex1, randomIndex2, randomIndex3);

        // Get PDB file of this triplet.
        AtomContainer ats;
        for (int sse1 = sse_pos[randomIndex1].second.first - 1; sse1 < sse_pos[randomIndex1].second.second; sse1++) {
          ats.addAtoms(sys.getPosition(sse1).getAtomPointers());
        }
        for (int sse2 = sse_pos[randomIndex2].second.first - 1; sse2 < sse_pos[randomIndex2].second.second; sse2++) {
          ats.addAtoms(sys.getPosition(sse2).getAtomPointers());
        }
        for (int sse3 = sse_pos[randomIndex3].second.first - 1; sse3 < sse_pos[randomIndex3].second.second; sse3++) {
          ats.addAtoms(sys.getPosition(sse3).getAtomPointers());
        }
        stringstream ss;
        ss << MslTools::getFileName(opt.pdb) << "_SSEs" << n << ".pdb";
        ats.writePdb(ss.str());
      }
    }
  }
}



bool checkBetaSheet(System &_sys, int _pos) {

  Position &pos = _sys.getPosition(_pos);
  if (!(pos.atomExists("N") && pos.atomExists("O"))) return false;

  Atom &posN = pos.getAtom("N");
  Atom &posO = pos.getAtom("O");

  for (uint i = 0; i < _sys.positionSize(); i++) {
    if (abs(((int)i) - _pos) <= 2) continue;

    Position &pos2 = _sys.getPosition(i);
    if (!(pos2.atomExists("N") && pos2.atomExists("O"))) continue;
    Atom &pos2N = pos2.getAtom("N");
    Atom &pos2O = pos2.getAtom("O");

    int hbonds = 0;
    // Anti-parallel
    // pos == pos2
    double dist = 3.5;
    if (posN.distance(pos2O) < dist) {
      hbonds++;
    }
    if (posO.distance(pos2N) < dist) {
      hbonds++;
    }

    if (hbonds == 2) { return true; }

    // Anti-parallel
    // Check pos-1 to i+1 AND pos+1 to i-1
    hbonds = 0;
    if (_pos > 0 && _sys.getPosition(_pos - 1).atomExists("N") && _sys.getPosition(_pos - 1).atomExists("O")) {

      if (i < _sys.positionSize() - 1 && _sys.getPosition(i + 1).atomExists("N") && _sys.getPosition(i + 1).atomExists("O")) {
        if (_sys.getPosition(_pos - 1).getAtom("N").distance(_sys.getPosition(i + 1).getAtom("O")) < dist ||
            _sys.getPosition(_pos - 1).getAtom("O").distance(_sys.getPosition(i + 1).getAtom("N")) < dist) {
          hbonds++;
        }
      }
    }

    if (_pos < _sys.positionSize() - 1 && _sys.getPosition(_pos + 1).atomExists("N") && _sys.getPosition(_pos + 1).atomExists("O")) {

      if (i > 0 && _sys.getPosition(i - 1).atomExists("N") && _sys.getPosition(i - 1).atomExists("O")) {
        if (_sys.getPosition(_pos + 1).getAtom("N").distance(_sys.getPosition(i - 1).getAtom("O")) < dist ||
            _sys.getPosition(_pos + 1).getAtom("O").distance(_sys.getPosition(i - 1).getAtom("N")) < dist) {
          hbonds++;
        }
      }
    }
    if (hbonds == 2) { return true; }

    // Parallel conditions
    // check pos to i-1 AND pos to i+1
    hbonds = 0;
    if (i > 0 && _sys.getPosition(i - 1).atomExists("N") && _sys.getPosition(i - 1).atomExists("O")) {

      if (posN.distance(_sys.getPosition(i - 1).getAtom("O")) < dist || posO.distance(_sys.getPosition(i - 1).getAtom("N")) < dist) {
        hbonds++;
      }
    }

    if (i < _sys.positionSize() - 1 && _sys.getPosition(i + 1).atomExists("N") && _sys.getPosition(i + 1).atomExists("O")) {

      if (posN.distance(_sys.getPosition(i + 1).getAtom("O")) < dist || posO.distance(_sys.getPosition(i + 1).getAtom("N")) < dist) {
        hbonds++;
      }
    }

    if (hbonds == 2) { return true; }

    // check pos-1 to i AND pos+1 to i
    if (_pos > 0 && _sys.getPosition(_pos - 1).atomExists("N") && _sys.getPosition(_pos - 1).atomExists("O")) {

      if (_sys.getPosition(_pos - 1).getAtom("N").distance(_sys.getPosition(i).getAtom("O")) < dist ||
          _sys.getPosition(_pos - 1).getAtom("O").distance(_sys.getPosition(i).getAtom("N")) < dist) {
        hbonds++;
      }
    }

    if (_pos < _sys.positionSize() - 1 && _sys.getPosition(_pos + 1).atomExists("N") && _sys.getPosition(_pos + 1).atomExists("O")) {

      if (_sys.getPosition(_pos + 1).getAtom("N").distance(_sys.getPosition(i).getAtom("O")) < dist ||
          _sys.getPosition(_pos + 1).getAtom("O").distance(_sys.getPosition(i).getAtom("N")) < dist) {
        hbonds++;
      }
    }

    if (hbonds == 2) { return true; }
  }

  return false;
}


Options setupOptions(int theArgc, char * theArgv[]){
  Options opt;

  OptionParser OP;

  OP.setRequired(opt.required);
  OP.setAllowed(opt.optional);
  OP.setDefaultArguments(opt.defaultArgs); // a pdb file value can be given as a default argument without the --pdbfile option
  OP.autoExtendOptions(); // if you give option "solvat" it will be autocompleted to "solvationfile"
  OP.readArgv(theArgc, theArgv);

  if (OP.countOptions() == 0){
    cout << "Usage:" << endl;
    cout << "  getSSEs --pdb <PDB_FILE> [options]" << endl;
    cout << endl;
    cout << "Required options:" << endl;
    for (uint i = 0; i < opt.required.size(); i++){
      cout << "  --" << opt.required[i] << endl;
    }
    cout << endl;
    cout << "Optional options:" << endl;
    for (uint i = 0; i < opt.optional.size(); i++){
      cout << "  --" << opt.optional[i] << endl;
    }
    cout << endl;
    exit(0);
  }

  opt.pdb = OP.getString("pdb");
  if (OP.fail()){
    cerr << "ERROR 1111 pdb not specified.\n";
    exit(1111);
  }

  opt.pick_random_SSEs = OP.getInt("pick_random_SSEs");
  if (OP.fail()){
    opt.pick_random_SSEs = 0;
  }

  opt.sasa  = OP.getBool("sasa");
  if (OP.fail()){
    opt.sasa = false;
  }

  opt.include_loops = OP.getBool("loops");
  if (OP.fail()){
    opt.include_loops = false;
  }

  opt.neighbor_dist = OP.getDouble("neighbor_dist");
  if (OP.fail()){
    opt.neighbor_dist = 6.0;
  }

  opt.num_contacting = OP.getInt("num_contacting");
  if (OP.fail()){
    opt.num_contacting = 6;
  }

  opt.msa = OP.getString("msa");
  if (OP.fail()){
    opt.msa = "";
  }
  opt.aa_freq_threshold = OP.getDouble("aa_freq_threshold");
  if (OP.fail()){
    opt.aa_freq_threshold = 0.0;
  }

  opt.exposed_pos_threshold = OP.getInt("exposed_pos_threshold");
  if (OP.fail()){
    opt.exposed_pos_threshold = 0;
  }

  MSLOUT.stream() << "Options:\n"<<OP<<endl;

  refSasa["G"] = 83.91;
  refSasa["A"] = 116.40;
  refSasa["S"] = 125.68;
  refSasa["C"] = 141.48;
  refSasa["P"] = 144.80;
  refSasa["T"] = 148.06;
  refSasa["D"] = 155.37;
  refSasa["V"] = 162.24;
  refSasa["N"] = 168.87;
  refSasa["E"] = 187.16;
  refSasa["Q"] = 189.17;
  refSasa["I"] = 189.95;
  refSasa["L"] = 197.99;
  refSasa["H"] = 198.51;
  refSasa["K"] = 207.49;
  refSasa["M"] = 210.55;
  refSasa["F"] = 223.29;
  refSasa["Y"] = 238.30;
  refSasa["R"] = 249.26;
  refSasa["W"] = 265.42;

  return opt;
}


vector<map<string, int> > computeAADiversity(System &sys, const string &msaFilePath, const string &pdbFileName) {
    ifstream msaFile(msaFilePath);
    if (!msaFile.is_open()) {
      cerr << "ERROR: Unable to open MSA file: " << msaFilePath << endl;
      exit(1);
    }

    // Read MSA sequences
    vector<string> msaSequences;
    string line, sequence, pdbSequenceMSA;
    int sequenceFound = 0;
    while (getline(msaFile, line)) {
      if (line.empty()) continue; // Skip empty lines
      if (line[0] == '>') {
        if (!sequence.empty()) {
          msaSequences.push_back(sequence);

          if (sequenceFound == 1) {
            pdbSequenceMSA = sequence;
            sequenceFound = 2; // Anything above 1 will flag as found
          }

          sequence.clear();
        }
        if (line.substr(1) == pdbFileName) {
          sequenceFound = 1;
        }
      } else {
        sequence += line;
      }
    }
    if (!sequence.empty()) {
      msaSequences.push_back(sequence);

      if (sequenceFound == 1) {
        pdbSequenceMSA = sequence;
        sequenceFound = 2; // Anything above 1 will flag as found
      }
    }
    msaFile.close();

    if (sequenceFound == 0) {
      cerr << "ERROR: No sequence name matching PDB file: " << pdbFileName << endl;
      exit(1);
    }

    // Check if the sequence matches the system
    stringstream tmp;
    for (size_t i = 0; i < sys.positionSize(); i++) {
      tmp << MslTools::getOneLetterCode(sys.getPosition(i).getResidueName());
    }
    string pdbSequence = tmp.str();

    // Remove dashes from pdbSequenceMSA
    string pdbSequenceMSA_copy = pdbSequenceMSA;
    pdbSequenceMSA_copy.erase(remove(pdbSequenceMSA_copy.begin(), pdbSequenceMSA_copy.end(), '-'), pdbSequenceMSA_copy.end());

    // Check if the sequence matches the system
    if (pdbSequenceMSA_copy != pdbSequence) {
      cerr << "ERROR: Sequence in MSA does not match the sequence in the PDB file." << endl;
      cerr << "PDB sequence: " << pdbSequence << endl;
      cerr << "MSA sequence: " << pdbSequenceMSA_copy << endl;
      exit(1);
    }

    // Create a mapping of indices from pdbSequenceMSA to pdbSequence
    map<int, int> pdbToMSAIndex;
    int pdbIndex = 0;
    for (size_t i = 0; i < pdbSequenceMSA.size(); ++i) {
      if (pdbSequenceMSA[i] != '-') {
        pdbToMSAIndex[pdbIndex] = i;
        pdbIndex++;
      }
    }

    // Compute amino acid diversity at each position using pdbToMSAIndex
    vector<map<string, int>> aaDiversity(sys.positionSize());
    for (const auto& seq : msaSequences) {
      for (size_t i = 0; i < sys.positionSize(); i++) {
        if (pdbToMSAIndex.find(i) != pdbToMSAIndex.end()) {
          int pdbIndex = pdbToMSAIndex[i];
          aaDiversity[pdbIndex][seq.substr(i,1)]++;
        }
      }
    }

    return aaDiversity;
  }