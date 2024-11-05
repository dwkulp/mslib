#include <iostream>
#include <cstdlib>
#include <queue>
#include <algorithm>
#include <map>
#include <utility>

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
#include "Frame.h"
#include "barrelDecider_v7.h"
#include "Line.h"
#include "PDBTopology.h"

using namespace std;
using namespace MSL;


// MslOut for nice output
static MslOut MSLOUT("betaDecider");
static SysEnv SYSENV;


int main(int argc, char *argv[]) {

  // Parse commandline options
  Options opt = setupOptions(argc,argv);

  // Read PDB structure
  System sys;
  sys.readPdb(opt.pdb);

  // Count the number of strand residues
  int strandCount = 0;

  // Store list of the position for start, end of each strand
  vector< pair<int,int> > strands;

  // Loop over all positions
  for (uint i = 0; i < sys.positionSize();i++){

    // Check if this position is in a beta sheet
    if (checkBetaSheet(sys,i)){

      // Print out only in debug mode
      MSLOUT.debug() << " SHEET: "<<sys.getPosition(i).getPositionId()<<endl;


      // Increment Strand Count
      strandCount++;


    } else {

      // This residue is not a strand, so could be end of strand

      // Add ending strand to list of strands if more than 4 residues in a row were beta
      if (strandCount > 4){
        strands.push_back(pair<int,int>(i-strandCount,i-1));
      }

      // Reset strand count
      strandCount = 0;
    }

  }


  // Output strand defintions to screen and make a PyMOL script containing selections for each strand
  PyMolVisualization pymol;

  cout << "Number of strands detected: "<<strands.size()<<endl;

  // Print out each strand
  for (uint i = 0; i < strands.size();i++){

    // Give strand a name
    string strandName = MslTools::stringf("Strand-%03d",i+1);

    // Print out strand defition
    cout << strandName <<" is between "<<sys.getPosition(strands[i].first).getPositionId()<< " and "<<sys.getPosition(strands[i].second).getPositionId()<<endl;

    // Save for pymol
    stringstream strandSel;
    strandSel << "chain "<<sys.getPosition(strands[i].first).getChainId()<< " and resi "<<sys.getPosition(strands[i].first).getResidueNumber()<<"-"<<sys.getPosition(strands[i].second).getResidueNumber();
    string strandSelStr = strandSel.str();
    string strandColor = "blue";
    pymol.createSelection(strandName,strandSelStr, strandColor);
  }


  // Write out PyMOL Script in PyMol use "run PDB-STRANDS.py"
  ofstream fout;
  fout.open(MslTools::stringf("%s-STRANDS.py", MslTools::getFileName(opt.pdb).c_str()));
  fout << pymol;


  // Create an adjacency list
  vector<vector<int> > adList;
  for (uint i = 0; i < strands.size();i++){
    
    vector<int> neighbors_of_strand_i;
    for (uint j = 0; j < strands.size();j++){
      if (i == j) continue;
      if (neighboringStrand(sys,strands[i],strands[j])){
        neighbors_of_strand_i.push_back(j);
      }
    }
    adList.push_back(neighbors_of_strand_i);
  }

  //Read in graphs
  Graph g(adList.size()+1);
  for (uint i = 0; i < adList.size();i++){
    //cout << "i: " << i+1 << endl;
    for (uint j = 0; j < adList[i].size();j++){
      //cout << "j: " << adList[i][j]+1 << endl;
      //cout << i+1 << ": " << adList[i][j]+1 << endl;
      g.addEdge(i+1,adList[i][j]+1);
    }
  }

  vector<int> g_path;

  if(g.isBarrel(g_path)) 
    cout << "Graph contains barrel\n"; 
  else
    cout << "Graph doesn't contain barrel\n";

  
  for (uint i = 0; i < g_path.size();i++) {
    cout << "GLOBAL PATH: " << g_path[i] << endl;
  }

  // get axis
  int min_strand = 0;
  min_strand = getShortestStrand(sys, strands);
  cout << "shortest strand " << min_strand << endl;

  pair<AtomPointerVector,AtomPointerVector> bTopBot= getEndAtoms(sys,strands,g_path);

  Frame top;
  top.setName(MslTools::stringf("Top"));
  top.computeFrameFromPCA(bTopBot.first);
  fout << top.toString();

  Frame bot;
  bot.setName(MslTools::stringf("Bottom"));
  bot.computeFrameFromPCA(bTopBot.second);
  fout << bot.toString();

  Line axis((top.getCenter()+bot.getCenter())/2,top.getCenter()-bot.getCenter());
  axis.setName("Baxis");
  axis.setOutputFormat("pymol");
  fout << axis.toString();
  fout.close();

  int outwardHydrophobics = 0;
  int allStrandResidues = 0;
  for (uint i = 0; i < g_path.size(); i++) {

    int start1 = strands[g_path[i]-1].first;
    int end1   = strands[g_path[i]-1].second;
    for (uint j = start1; j <= end1;j++){

      Atom *CB = NULL;
      if (sys.getPosition(j).getResidueName() == "GLY"){
	CB = PDBTopology::getPseudoCbeta(sys.getPosition(j).getCurrentIdentity());
      } else {
	CB =&sys.getPosition(j).getAtom("CB");
      }
      double dist = CartesianGeometry::distance(CB->getCoor(),axis.projection(CB->getCoor()));
      cout << "DIST TO AXIS: "<<sys.getPosition(j).getCurrentIdentity().getIdentityId()<<" "<<dist<<endl;
      if (dist > 6 &&
	  (CB->getResidueName() == "PHE" ||
	  CB->getResidueName() == "LEU" ||
	  CB->getResidueName() == "ILE" ||
	  CB->getResidueName() == "VAL" ||
	  CB->getResidueName() == "TYR" ||
	  CB->getResidueName() == "MET" ||
	  CB->getResidueName() == "TRP" ||
	   CB->getResidueName() == "ALA")){
	outwardHydrophobics++;
     
      }
      allStrandResidues++;
    }

  }
  double freq = (double)outwardHydrophobics/(double)allStrandResidues*100;
  cout << "OUTWARD HYDROPHOBICS: "<<outwardHydrophobics<<" out of "<<allStrandResidues<<" percent: "<<freq<<endl;

  // top
  //AtomPointerVector apv_top = getTopAtoms(sys, strands, min_strand);
  //cout << apv_top.getGeometricCenter() << endl;

  // bottom
  //AtomPointerVector apv_bot = getBotAtoms(sys, strands, min_strand);
  //cout << apv_bot.getGeometricCenter() << endl;

}


// Function to help setup input options
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
    cout << endl;
    cout << "betaDecider --pdb PDB\n";

    cout << "\nprogram options: "<<endl;
    for (uint i = 0; i < opt.required.size();i++){
      cout <<"R  --"<<opt.required[i]<<"  "<<endl;
    }
    cout <<endl;
    for (uint i = 0; i < opt.optional.size();i++){
      cout <<"O  --"<<opt.optional[i]<<"  "<<endl;
    }
    cout << endl;
    exit(0);
  }

  opt.pdb = OP.getString("pdb");
  if (OP.fail()){
    cerr << "ERROR 1111 pdb not specified.\n";
    exit(1111);
  }

  return opt;
}
pair<AtomPointerVector,AtomPointerVector> getEndAtoms(System &_sys, vector< pair<int,int> > &_strands,vector<int> &_path){

    AtomPointerVector side1ats;
    AtomPointerVector side2ats;
    for (uint i = 0; i < _path.size(); i++) {

      int side1 = _strands[_path[i]-1].first;
      int side2 = _strands[_path[i]-1].second;
      cout <<"I0: "<<i<<" "<<_path[i]<<" "<<_path[i]<<_strands[_path[i]-1].first<<" "<<_strands[_path[i]-1].second<<endl;
      cout<< "I1: "<<i<<_sys.getPosition(side1).getAtom("CA").getPositionId()<<endl;
      cout<< "I2: "<<i<<_sys.getPosition(side2).getAtom("CA").getPositionId()<<endl;

      if (i == 0){
	side1ats.push_back(&_sys.getPosition(side1).getAtom("CA"));
	side2ats.push_back(&_sys.getPosition(side2).getAtom("CA"));
      }  else {
      
	int side1last = (uint)side1ats.size()-1;
	int side2last = (uint)side2ats.size()-1;
	double dist1 = side1ats[side1last]->distance(_sys.getPosition(_strands[_path[i]-1].first).getAtom("CA")); 
	double dist2 = side2ats[side2last]->distance(_sys.getPosition(_strands[_path[i]-1].first).getAtom("CA")); 	

	if (dist1 < dist2){
	  side1ats.push_back(&_sys.getPosition(_strands[_path[i]-1].first).getAtom("CA"));
	  side2ats.push_back(&_sys.getPosition(_strands[_path[i]-1].second).getAtom("CA"));
	} else {
	  side1ats.push_back(&_sys.getPosition(_strands[_path[i]-1].second).getAtom("CA"));
	  side2ats.push_back(&_sys.getPosition(_strands[_path[i]-1].first).getAtom("CA"));
	}
      }
    }

    return (pair<AtomPointerVector,AtomPointerVector>(side1ats,side2ats));

}

int getShortestStrand(System &_sys, vector< pair<int,int> > &_strands){
  
  int min = 0;
  int min_strand = 0;
  int length = 0;

  for (uint i = 0; i < _strands.size(); i++) {
    int one = _sys.getPosition(_strands[i].first).getIndexInSystem();
    int two = _sys.getPosition(_strands[i].second).getIndexInSystem();

    for (uint j = one; j < two; j++ ) {
        // cout << _sys.getPosition(j).getAtom("CA") << endl;
        length++;
    }

    if(min == 0){
        min = length;
        min_strand = i;
    }

    else if(min>length){
        min = length;
        min_strand = i;
    }

    cout << "min " << min << " length " << length << endl;

    length = 0;

  }
  return min_strand;
}

AtomPointerVector getTopAtoms(System &_sys, vector< pair<int,int> > &_strands, int min){

    AtomPointerVector top_atoms;

    int first = _sys.getPosition(_strands[min].first).getIndexInSystem();

    for (uint i = 0; i < _strands.size(); i++) {

        int one = _sys.getPosition(_strands[i].first).getIndexInSystem();
        int two = _sys.getPosition(_strands[i].second).getIndexInSystem();

        double dist_min = 0;
        int pos_min = 0;

        for (uint j = one; j < two; j++ ) {

            double dist = _sys.getPosition(first).getAtom("CA").distance(_sys.getPosition(j).getAtom("CA"));

            //cout << "distance " << dist << endl;

            if(dist_min == 0) {

                dist_min = dist;
                pos_min = j;

                //cout << "step 1" << endl;
            }

            else if(dist_min > dist){
                dist_min = dist;
                pos_min = j;

                //cout << "step 2" << endl;
            }

        }

        cout << dist_min << endl;
        cout << _sys.getPosition(pos_min).getAtom("CA") << endl;
        top_atoms.push_back(&_sys.getPosition(pos_min).getAtom("CA"));

    }
    return top_atoms;
}

AtomPointerVector getBotAtoms(System &_sys, vector< pair<int,int> > &_strands, int min){

    AtomPointerVector bot_atoms;

    int second = _sys.getPosition(_strands[min].second).getIndexInSystem();

    for (uint i = 0; i < _strands.size(); i++) {

        int one = _sys.getPosition(_strands[i].first).getIndexInSystem();
        int two = _sys.getPosition(_strands[i].second).getIndexInSystem();

        double dist_min = 0;
        int pos_min = 0;

        for (uint j = one; j < two; j++ ) {

            double dist = _sys.getPosition(second).getAtom("CA").distance(_sys.getPosition(j).getAtom("CA"));

            //cout << "distance " << dist << endl;

            if(dist_min == 0) {

                dist_min = dist;
                pos_min = j;

                //cout << "step 1" << endl;
            }

            else if(dist_min > dist){
                dist_min = dist;
                pos_min = j;

                //cout << "step 2" << endl;
            }

        }

        cout << dist_min << endl;
        cout << _sys.getPosition(pos_min).getAtom("CA") << endl;
        bot_atoms.push_back(&_sys.getPosition(pos_min).getAtom("CA"));

    }
    return bot_atoms;
}

bool checkBetaSheet(System &_sys, int _pos){
  Position &pos = _sys.getPosition(_pos);
  if (!pos.atomExists("CA")) return false;
  //cout << "Pos: "<<pos.getPositionId()<<endl;

  if (pos.getPhi() < -35 &&
        pos.getPhi() > -90 &&
        pos.getPsi() > -70 &&
        pos.getPsi() < 0){
       //HELIX
       return false;
    } else {
      double dist_NO = 3.4;
      int hb_num = 2;

      if (!(pos.atomExists("N") && pos.atomExists("O"))) return false;

      Atom &posN = pos.getAtom("N");
      Atom &posO = pos.getAtom("O");

      for (uint i = 0; i< _sys.positionSize();i++){
        if (abs(((int)i)-_pos) <= hb_num) continue;

        Position &pos2 = _sys.getPosition(i);
        if (!(pos2.atomExists("N") && pos2.atomExists("O"))) continue;
        Atom &pos2N = pos2.getAtom("N");
        Atom &pos2O = pos2.getAtom("O");

        int hbonds = 0;
        // Anti-parallel
        // pos == pos2
        if (posN.distance(pos2O) < dist_NO){
            hbonds++;
        }
        if (posO.distance(pos2N) < dist_NO){
            hbonds++;
        }

        if (hbonds == hb_num)  { return true;}
    
        // Anti-parallel
        // Check pos-1 to i+1 AND pos+1 to i-1
        hbonds = 0;
        if (_pos > 0 && _sys.getPosition(_pos-1).atomExists("N") && _sys.getPosition(_pos-1).atomExists("O")){

            if (i < _sys.positionSize()-1 && _sys.getPosition(i+1).atomExists("N") && _sys.getPosition(i+1).atomExists("O")){
            if (_sys.getPosition(_pos-1).getAtom("N").distance(_sys.getPosition(i+1).getAtom("O")) < dist_NO ||
                _sys.getPosition(_pos-1).getAtom("O").distance(_sys.getPosition(i+1).getAtom("N")) < dist_NO){
              hbonds++;
            }
            }
        }

      if (_pos < _sys.positionSize()-1 && _sys.getPosition(_pos+1).atomExists("N") && _sys.getPosition(_pos+1).atomExists("O")){

          if (i > 0 && _sys.getPosition(i-1).atomExists("N") && _sys.getPosition(i-1).atomExists("O")){
          if (_sys.getPosition(_pos+1).getAtom("N").distance(_sys.getPosition(i-1).getAtom("O")) < dist_NO ||
              _sys.getPosition(_pos+1).getAtom("O").distance(_sys.getPosition(i-1).getAtom("N")) < dist_NO){
            hbonds++;
          }
          }
      }
      if (hbonds == hb_num) {return true;}


      // Parallel conditions
      // check pos to i-1 AND pos to i+1
      hbonds= 0;
      if (i > 0 && _sys.getPosition(i-1).atomExists("N") && _sys.getPosition(i-1).atomExists("O")){
      
          if (posN.distance(_sys.getPosition(i-1).getAtom("O")) < dist_NO || posO.distance(_sys.getPosition(i-1).getAtom("N")) < dist_NO){
          hbonds++;
          }
      }

      if (i < _sys.positionSize()-1 && _sys.getPosition(i+1).atomExists("N") && _sys.getPosition(i+1).atomExists("O")){
      
          if (posN.distance(_sys.getPosition(i+1).getAtom("O")) < dist_NO || posO.distance(_sys.getPosition(i+1).getAtom("N")) < dist_NO){
          hbonds++;
          }
      }
   
      if (hbonds == hb_num ) {return true;}

      // check pos-1 to i AND pos+1 to i
      if (_pos > 0 && _sys.getPosition(_pos-1).atomExists("N") && _sys.getPosition(_pos-1).atomExists("O")){
      
          if (_sys.getPosition(_pos-1).getAtom("N").distance(_sys.getPosition(i).getAtom("O")) < dist_NO ||
          _sys.getPosition(_pos-1).getAtom("O").distance(_sys.getPosition(i).getAtom("N")) < dist_NO){
          hbonds++;
          }
      }

      if (_pos < _sys.positionSize()-1 && _sys.getPosition(_pos+1).atomExists("N") && _sys.getPosition(_pos+1).atomExists("O")){
      
          if (_sys.getPosition(_pos+1).getAtom("N").distance(_sys.getPosition(i).getAtom("O")) < dist_NO || 
          _sys.getPosition(_pos+1).getAtom("O").distance(_sys.getPosition(i).getAtom("N")) < dist_NO){
          hbonds++;
          }
      }

      if (hbonds == hb_num) {return true;}
      }

      return false;
    }
}


bool neighboringStrand(System &_sys, pair<int,int> &_s1, pair<int,int> &_s2){

  double nb_dist = 3.4;
  int hb_num = 2;
  int hbonds=0;
  for (uint i = _s1.first; i <= _s1.second;i++){
    for (uint j = _s2.first; j <= _s2.second;j++){
      if (_sys.getPosition(i).atomExists("N") && _sys.getPosition(j).atomExists("O")){
  if (_sys.getPosition(i).getAtom("N").distance(_sys.getPosition(j).getAtom("O")) < nb_dist){
    hbonds++;
  }
      }
      if (_sys.getPosition(i).atomExists("O") && _sys.getPosition(j).atomExists("N")){
  if (_sys.getPosition(i).getAtom("O").distance(_sys.getPosition(j).getAtom("N")) < nb_dist){
    hbonds++;
  }
      }
    }
  }

  if (hbonds > hb_num-1){
    return true;
  }

  return false;

}

Graph::Graph(int V) 
{ 
    this->V = V; 
    adj = new list<int>[V]; 
} 
  
void Graph::addEdge(int v, int w) 
{ 
    adj[v].push_back(w); // Add w to v’s list. 
} 

bool Graph::seeRoot(int v, int root){
    list<int>::iterator i;
    for(i = adj[v].begin(); i != adj[v].end(); ++i) {
        if(*i == root){return true;}
    }
    return false;
}

// Recursive function to check whether root connects back to root
bool Graph::rootConnect(int last, int v, int root, bool visited[], vector<int> path, vector<int> &_globalPath){

  // check visit
  if(visited[v]){
    path.pop_back();
    return false;
  }

  visited[v] = true;
  
  // connect check
  if(seeRoot(v,root) && path.size()>4){
    _globalPath = path;
    string road = "";
    for (int p=0; p<path.size(); p++){
      road += to_string(path[p]) + " ";
    }
    cout << "path: " << road << endl;
    cout << "strands number in barrel: " << path.size() << endl;
    return true;
  }

  else{
    // Recur for all the vertices adjacent to this vertex 
    list<int>::iterator i; 
    for(i = adj[v].begin(); i != adj[v].end(); ++i) {

        if(*i != last){
            path.push_back(*i);
            if(rootConnect(v,*i,root,visited,path,_globalPath)){return true;}
        }
        
    }
    path.pop_back();
    return false;
  }
}

// Whether a barrel exist from current root
bool Graph::rootBarrel(int root, vector<int> &_globalPath){
  // reset visited
  bool *visited = new bool[V];
  vector<int> path;

  for(int i = 0; i < V; i++) 
  { 
    visited[i] = false;
  } 

  // Recur for all the vertices adjacent to this vertex 
  list<int>::iterator i; 

  path.push_back(root);
  for(i = adj[root].begin(); i != adj[root].end(); ++i){
    path.push_back(*i);
    if(rootConnect(root,*i, root, visited, path,_globalPath)){return true;}
  }
  return false;
}

// Return true if is_barrel
bool Graph::isBarrel(vector<int> &_globalPath){
  for(int i = 0; i < V; i++){
    if(rootBarrel(i,_globalPath)){return true;}
  }
  return false;
}
