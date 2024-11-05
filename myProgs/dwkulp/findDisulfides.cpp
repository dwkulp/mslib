#include "MslTools.h"
#include "OptionParser.h"
#include "release.h"
#include "System.h"
#include "Chain.h"
#include "Residue.h"
#include "PDBFormat.h"
#include "PDBWriter.h"
#include "AtomContainer.h"
#include "Transforms.h"
#include "Timer.h"

#include "MslOut.h"
#include "findDisulfides.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;
using namespace MSL;
using namespace MslTools;

// MslOut 
static MslOut MSLOUT("findDisulfides");

//void readDisulfPdb(string _pdbfile,vector<pair<AtomContainer *, AtomContainer *> > &_container);
void readDisulfPdb(string _pdbfile,vector<pair<AtomContainer *, AtomContainer *> > &_container, bool _binary);

int main(int argc, char *argv[]) {

    Options opt = setupOptions(argc, argv);
    Transforms tm;

    System sys;
    sys.readPdb(opt.pdb);

    Timer t;
    double start = t.getWallTime();
    vector<pair<AtomContainer *, AtomContainer *> > disulfs;
    readDisulfPdb(opt.disulfPdb,disulfs,opt.binaryDB);
    MSLOUT.stream() << "Done reading disulfide database: "<<disulfs.size()<<" took "<<(t.getWallTime() - start)<<" seconds."<<endl;


    // Flush the stdout buffer
    fflush(stdout);
    
    PDBWriter pout;
    int cys_index = 1;

    ofstream fastaFile;
    if (opt.fasta) {
      fastaFile.open(MslTools::stringf("DS_%s.fasta", MslTools::getFileName(opt.pdb).c_str()));
    }
	  
    // Search all pairs of residues
    if (opt.specific_positions.size() == 0){

      // Instead try all pairs of positions on the first chain, when Ca are < 6 Angstroms from each other and not next to one another
      for (uint p1 = 0; p1 < sys.positionSize();p1++){
	Position &pos1 = sys.getPosition(p1);
	if (!pos1.atomExists("CA")) { continue; }

	for (uint p2 = p1+2; p2 < sys.positionSize();p2++){
	  Position &pos2 = sys.getPosition(p2);
	  if (!pos2.atomExists("CA")) { continue; }

	  if (pos1.getAtom("CA").distance(pos2.getAtom("CA")) < 6){
	    opt.specific_positions.push_back(pos1.getPositionId());
	    opt.specific_positions.push_back(pos2.getPositionId());
	  }	  

	}
      }
    }

    

    // Search specific disulfides
    for (uint i = 0 ; i < opt.specific_positions.size();i+=2){
      if (!sys.positionExists(opt.specific_positions[i])){
	cerr << "ERROR 2222 Position: "<<opt.specific_positions[i]<< " does not exist in "<<opt.pdb<<endl;
	exit(2222);
      }
      if (!sys.positionExists(opt.specific_positions[i+1])){
	cerr << "ERROR 2222 Position: "<<opt.specific_positions[i+1]<< " does not exist in "<<opt.pdb<<endl;
	exit(2222);
      }
	
      Position &p1 = sys.getPosition(opt.specific_positions[i]);
      if (! (p1.atomExists("N") && p1.atomExists("CA") && p1.atomExists("C")) ){
	cerr << "WARNING position: "<<p1<<" is missing an atom: N,CA or C "<<endl;
	continue;
      }
      AtomPointerVector pos_bb;
      pos_bb.push_back(&p1.getAtom("N"));
      pos_bb.push_back(&p1.getAtom("CA"));
      pos_bb.push_back(&p1.getAtom("C"));

      Position &p2 = sys.getPosition(opt.specific_positions[i+1]);
      if (! (p2.atomExists("N") && p2.atomExists("CA") && p2.atomExists("C")) ){
	cerr << "WARNING position: "<<p2<<" is missing an atom: N,CA or C "<<endl;
	continue;
      }

      pos_bb.push_back(&p2.getAtom("N"));
      pos_bb.push_back(&p2.getAtom("CA"));
      pos_bb.push_back(&p2.getAtom("C"));


      // Skip local disulfides
      if (opt.skip_local_disulfides != -1 && p1.getChainId() == p2.getChainId() && abs(p1.getResidueNumber() - p2.getResidueNumber()) < opt.skip_local_disulfides) continue;

      // Only consider inter-chain
      if (opt.only_inter_chain && p1.getChainId() == p2.getChainId()) continue;

      double lowestRMSD =MslTools::doubleMax;
      pair<AtomPointerVector,AtomPointerVector>  bestDisulfide;
      int numDisulfs = 0;
      for (uint d = 0; d < disulfs.size();d++){
	try {
	  string id   = "";

	  if (disulfs[d].first->atomSize() == 0 || disulfs[d].second->atomSize() == 0){
	    MSLOUT.stream() << "Disulfide [ "<<d<<" ] has 0 atoms."<<endl;
	    continue;
	  }
	  AtomPointerVector cys1_bb;
	  id   = MslTools::getAtomId(disulfs[d].first->getAtom(0).getChainId(),
				     disulfs[d].first->getAtom(0).getResidueNumber(),
				     disulfs[d].first->getAtom(0).getResidueIcode(),
				     "N");
	  if (!disulfs[d].first->atomExists(id)){
	    MSLOUT.stream() << "Disulfide missing atom: "<<id<<endl;
	    continue;
	  }
	  cys1_bb.push_back( &disulfs[d].first->getAtom(id));

	  id   = MslTools::getAtomId(disulfs[d].first->getAtom(0).getChainId(),
				     disulfs[d].first->getAtom(0).getResidueNumber(),
				     disulfs[d].first->getAtom(0).getResidueIcode(),
				     "CA");
	  if (!disulfs[d].first->atomExists(id)){
	    MSLOUT.stream() << "Disulfide missing atom: "<<id<<endl;
	    continue;
	  }
	  cys1_bb.push_back( &disulfs[d].first->getAtom(id));

	  id   = MslTools::getAtomId(disulfs[d].first->getAtom(0).getChainId(),
				     disulfs[d].first->getAtom(0).getResidueNumber(),
				     disulfs[d].first->getAtom(0).getResidueIcode(),
				     "C");
	  if (!disulfs[d].first->atomExists(id)){
	    MSLOUT.stream() << "Disulfide missing atom: "<<id<<endl;
	    continue;
	  }
	  cys1_bb.push_back( &disulfs[d].first->getAtom(id));

	  id   = MslTools::getAtomId(disulfs[d].first->getAtom(0).getChainId(),
				     disulfs[d].first->getAtom(0).getResidueNumber(),
				     disulfs[d].first->getAtom(0).getResidueIcode(),
				     "SG");
	  if (!disulfs[d].first->atomExists(id)){
	    MSLOUT.stream() << "Disulfide missing atom: "<<id<<endl;
	    continue;
	  }
	  Atom &SG1 = disulfs[d].first->getAtom(id);


	  AtomPointerVector cys2_bb;
	  id   = MslTools::getAtomId(disulfs[d].second->getAtom(0).getChainId(),
				     disulfs[d].second->getAtom(0).getResidueNumber(),
				     disulfs[d].second->getAtom(0).getResidueIcode(),
				     "N");
	  if (!disulfs[d].second->atomExists(id)){
	    MSLOUT.stream() << "Disulfide missing atom: "<<id<<endl;
	    continue;
	  }
	  cys2_bb.push_back( &disulfs[d].second->getAtom(id));

	  id   = MslTools::getAtomId(disulfs[d].second->getAtom(0).getChainId(),
				     disulfs[d].second->getAtom(0).getResidueNumber(),
				     disulfs[d].second->getAtom(0).getResidueIcode(),
				     "CA");
	  if (!disulfs[d].second->atomExists(id)){
	    MSLOUT.stream() << "Disulfide missing atom: "<<id<<endl;
	    continue;
	  }
	  cys2_bb.push_back( &disulfs[d].second->getAtom(id));

	  id   = MslTools::getAtomId(disulfs[d].second->getAtom(0).getChainId(),
				     disulfs[d].second->getAtom(0).getResidueNumber(),
				     disulfs[d].second->getAtom(0).getResidueIcode(),
				     "C");
	  if (!disulfs[d].second->atomExists(id)){
	    MSLOUT.stream() << "Disulfide missing atom: "<<id<<endl;
	    continue;
	  }
	  cys2_bb.push_back( &disulfs[d].second->getAtom(id));


	  
	  id   = MslTools::getAtomId(disulfs[d].second->getAtom(0).getChainId(),
				     disulfs[d].second->getAtom(0).getResidueNumber(),
				     disulfs[d].second->getAtom(0).getResidueIcode(),
				     "SG");
	  if (!disulfs[d].second->atomExists(id)){
	    MSLOUT.stream() << "Disulfide missing atom: "<<id<<endl;
	    continue;
	  }
	  Atom &SG2 = disulfs[d].second->getAtom(id);
	  
	  double sulfur_dist = SG1.distance(SG2);
	  if (sulfur_dist > 2.3 || sulfur_dist < 1.8){
	    //MSLOUT.stream() << "OH NO! Disulfide not formed. "<<sulfur_dist<<endl;
	    continue;
	  }
	  


	  AtomPointerVector cys1_cys2_ats = disulfs[d].first->getAtomPointers() + disulfs[d].second->getAtomPointers();
	  AtomContainer cys_fwd_bb;
	  cys_fwd_bb.addAtoms(cys1_bb);
	  cys_fwd_bb.addAtoms(cys2_bb);

	  AtomContainer cys_rev_bb;
	  cys_rev_bb.addAtoms(cys2_bb);
	  cys_rev_bb.addAtoms(cys1_bb);


	  // Alignment time.
	  double rmsd1 = MslTools::doubleMax;
	  if (!tm.rmsdAlignment(cys_fwd_bb.getAtomPointers(), pos_bb)){
	    MSLOUT.stream() << "ERROR BB alignment1"<<endl;
	    continue;
	  } else {
	    rmsd1 = cys_fwd_bb.getAtomPointers().rmsd(pos_bb);
	  }

	  if (rmsd1 < opt.tol){
	    
	    MSLOUT.stream() << "Found a disulfide1: "<<rmsd1<<" "<<sulfur_dist<<endl;

	    AtomPointerVector cys1_cys2_bb = cys1_bb+cys2_bb;
	    if (!tm.rmsdAlignment(cys1_cys2_bb, pos_bb,cys1_cys2_ats)){
	      MSLOUT.stream() << "ERROR FULL alignment1"<<endl;
	      continue;
	    } 

	    if (opt.writeOutAll){
	      string fname = MslTools::stringf("%s_cys_cys_%04d_%04d_%06d.pdb", MslTools::getFileName(opt.pdb).c_str(),cys_index++);
	      pout.open(fname);
	      pout.write(cys1_cys2_ats);
	      pout.close();
	    }

	    if (rmsd1 < lowestRMSD){
	      lowestRMSD = rmsd1;
	      bestDisulfide.first  = disulfs[d].first->getAtomPointers();
	      bestDisulfide.second = disulfs[d].second->getAtomPointers();
	    }

	    numDisulfs++;
	    continue;
	  }

	  double rmsd2 = MslTools::doubleMax;
	  if (!tm.rmsdAlignment(cys_rev_bb.getAtomPointers(), pos_bb)){
	    MSLOUT.stream() << "ERROR alignment1"<<endl;
	    continue;
	  } else {
	    rmsd2 = cys_rev_bb.getAtomPointers().rmsd(pos_bb);
	  }

	  if (rmsd2 < opt.tol){

	    MSLOUT.stream() << "Found a disulfide2: "<<rmsd2<<" "<<sulfur_dist<<endl;

	    AtomPointerVector cys2_cys1_bb = cys2_bb+cys1_bb;
	    if (!tm.rmsdAlignment(cys2_cys1_bb, pos_bb,cys1_cys2_ats)){
	      MSLOUT.stream() << "ERROR FULL alignment1"<<endl;
	      continue;
	    } 
	    if (opt.writeOutAll){
	      string fname = MslTools::stringf("%s_cys_cys_%04d_%04d_%06d.pdb", MslTools::getFileName(opt.pdb).c_str(),p1.getResidueNumber(),p2.getResidueNumber(),cys_index++);
	      pout.open(fname);
	      pout.write(cys1_cys2_ats);
	      pout.close();
	    }

	    if (rmsd2 < lowestRMSD){
	      lowestRMSD = rmsd2;
	      bestDisulfide.first  = disulfs[d].second->getAtomPointers();
	      bestDisulfide.second = disulfs[d].first->getAtomPointers();
	    }

	    numDisulfs++;
	    continue;
	  }	  

	  //MSLOUT.fprintf(stdout,"RMSDs: %6.2f  %6.2f\n",rmsd1,rmsd2);
	} catch (...){
	  cout << "EXCEPTION CAUGHT, try to continue: "<<d<<endl;
	  continue;
	}
      }
      string nativeDisulfFlag = "";
      if (p1.getResidueName() == "CYS" && p2.getResidueName() == "CYS"){
	nativeDisulfFlag = "****";
      }

      if (numDisulfs > 0){

	AtomContainer disulf1(bestDisulfide.first);
	AtomContainer disulf2(bestDisulfide.second);


	string cb1   = MslTools::getAtomId(disulf1.getAtom(0).getChainId(),
					   disulf1.getAtom(0).getResidueNumber(),
					   disulf1.getAtom(0).getResidueIcode(),
					   "CB");

	string sg1   = MslTools::getAtomId(disulf1.getAtom(0).getChainId(),
					   disulf1.getAtom(0).getResidueNumber(),
					   disulf1.getAtom(0).getResidueIcode(),
					   "SG");

	string cb2   = MslTools::getAtomId(disulf2.getAtom(0).getChainId(),
					   disulf2.getAtom(0).getResidueNumber(),
					   disulf2.getAtom(0).getResidueIcode(),
					   "CB");

	string sg2   = MslTools::getAtomId(disulf2.getAtom(0).getChainId(),
					   disulf2.getAtom(0).getResidueNumber(),
					   disulf2.getAtom(0).getResidueIcode(),
					   "SG");

	double sg1_sg2         = disulf1.getAtom(sg1).distance(disulf2.getAtom(sg2));
	double cb1_sg1_sg2     = disulf1.getAtom(cb1).angle(disulf1.getAtom(sg1),disulf2.getAtom(sg2));
	double cb2_sg2_sg1     = disulf2.getAtom(cb2).angle(disulf2.getAtom(sg2),disulf1.getAtom(sg1));
	double cb1_sg1_sg2_cb2 = disulf1.getAtom(cb1).dihedral(disulf1.getAtom(sg1),disulf2.getAtom(sg2),disulf2.getAtom(cb2));

	MSLOUT.fprintf(stdout, "DATA: %12s %12s , %6d disulfides, best RMSD: %6.2f, Geometry: %7.2f %7.2f %7.2f %7.2f, %s\n",
		       p1.getCurrentIdentity().getIdentityId().c_str(),
		       p2.getCurrentIdentity().getIdentityId().c_str(),
		       numDisulfs, 
		       lowestRMSD,
		       sg1_sg2,
		       cb1_sg1_sg2,
		       cb2_sg2_sg1,
		       cb1_sg1_sg2_cb2,
		       nativeDisulfFlag.c_str());

	if (opt.fasta){

	  fastaFile << MslTools::stringf(">%s%s_DS_%04d\n",MslTools::getFileName(opt.pdb).c_str(),opt.appendName.c_str(), i);
	  for (uint c = 0; c < sys.chainSize();c++){
	    for (uint r = 0; r < sys.getChain(c).positionSize();r++){
	      Residue &res = sys.getChain(c).getResidue(r);
	      string aa = MslTools::getOneLetterCode(res.getResidueName());
	      
	      if (res.getPositionId() == sys.getPosition(opt.specific_positions[i]).getPositionId() || res.getPositionId() == sys.getPosition(opt.specific_positions[i+1]).getPositionId()) {
		aa = "C";
	      }
	      if (aa != "X")
		fastaFile << MslTools::stringf("%1s",aa.c_str());
	    }
	  }
	  fastaFile << endl;
	  
      }

      // Make model of lowest.
      if (opt.modelBest){

	try{
	  if (numDisulfs == 0 || lowestRMSD == MslTools::doubleMax){
	    MSLOUT.stream() << "NO DISULFIDE HAS BEEN FOUND WITH MATCHING GEOMETRY, try increasing the --tol arguement or a larger --disulfPdb database\n";
	    continue;
	  }

	  MSLOUT.stream() << "BEST DISFULIDE HAS RMSD: "<<lowestRMSD<<endl;


	  System newSys(sys);
	  Position &new_p1 = newSys.getPosition(opt.specific_positions[i]);
	  new_p1.removeIdentity("CYS"); // in case CYS is already there
	  new_p1.addIdentity(bestDisulfide.first,"CYS");
	  new_p1.setActiveIdentity("CYS");
	  new_p1.getCurrentIdentity().getAtom("N").setCoor(p1.getCurrentIdentity().getAtom("N").getCoor());
	  new_p1.getCurrentIdentity().getAtom("C").setCoor(p1.getCurrentIdentity().getAtom("C").getCoor());
	  new_p1.getCurrentIdentity().getAtom("O").setCoor(p1.getCurrentIdentity().getAtom("O").getCoor());

	  Position &new_p2 = newSys.getPosition(opt.specific_positions[i+1]);
	  new_p2.removeIdentity("CYS"); // in case CYS is already there
	  new_p2.addIdentity(bestDisulfide.second,"CYS");
	  new_p2.setActiveIdentity("CYS");
	  new_p2.getCurrentIdentity().getAtom("N").setCoor(p2.getCurrentIdentity().getAtom("N").getCoor());
	  new_p2.getCurrentIdentity().getAtom("C").setCoor(p2.getCurrentIdentity().getAtom("C").getCoor());
	  new_p2.getCurrentIdentity().getAtom("O").setCoor(p2.getCurrentIdentity().getAtom("O").getCoor());


	  string fname = MslTools::stringf("%s_bestDisulf_%04d_%04d.pdb", MslTools::getFileName(opt.pdb).c_str(),p1.getResidueNumber(),p2.getResidueNumber());

	  newSys.writePdb(fname);


	} catch (...) {
	  cout << "ERROR 3433 exception caught with best model: "<<i<<" "<<lowestRMSD<<" "<<opt.specific_positions[i]<<" "<<bestDisulfide.first.size()<<" "<<bestDisulfide.second.size()<<endl;
	  exit(3433);
	}

      }
      
      }

      // Flush the stdout buffer
      fflush(stdout);
      
    }
      
    if (opt.fasta) {
      fastaFile.close();
    }
    MSLOUT.stream() << "Done."<<endl;
}

Options setupOptions(int theArgc, char * theArgv[]){
    // Create the options
    Options opt;

    // Parse the options
    OptionParser OP;
    OP.setRequired(opt.required);	
    OP.setDefaultArguments(opt.defaultArgs); // the default argument is the --configfile option
    OP.readArgv(theArgc, theArgv);

    if (OP.countOptions() == 0){
	cout << "Usage: getTripletCaMeasurements " << endl;
	cout << endl;
	cout << "\n";
	cout << "pdblist PDB\n";
	cout << endl;
	exit(0);
    }

    opt.pdb = OP.getString("pdb");
    if (OP.fail()){
	cerr << "ERROR 1111 no pdb specified."<<endl;
	exit(1111);
    }
    opt.disulfPdb = OP.getString("disulfPdb");
    if (OP.fail()){
	cerr << "ERROR 1111 no disulfPdb specified."<<endl;
	exit(1111);
    }
    opt.specific_positions = OP.getStringVectorJoinAll("positions");
    opt.tol = OP.getDouble("tol");
    if (OP.fail()){
      opt.tol = 0.5;
    }
    opt.modelBest = OP.getBool("modelBest");
    if (OP.fail()){
      opt.modelBest = false;
    }
    opt.writeOutAll = OP.getBool("writeOutAll");
    if (OP.fail()){
      opt.writeOutAll = false;
    }

    opt.skip_local_disulfides = OP.getInt("skip_local_disulfides");
    if (OP.fail()){
      opt.skip_local_disulfides = -1;
    }
    opt.only_inter_chain = OP.getBool("only_inter_chain");
    if (OP.fail()){
      opt.only_inter_chain = false;
    }
    opt.fasta = OP.getBool("fasta");
    opt.appendName = OP.getString("appendName");
    if (OP.fail()){
      opt.appendName="";
    }
      
    opt.binaryDB = OP.getBool("binaryDB");
    if (OP.fail()){
      opt.binaryDB = false;
    }
    return opt;
}

void readDisulfPdb(string _pdbfile,vector<pair<AtomContainer *, AtomContainer *> > &_container, bool _binary){
  

     if (!_binary){

	ifstream fs;
	fs.open(_pdbfile.c_str());
	if (fs.fail()) {
	  return;
	}

	string currentResidue = "";
	vector<PDBFormat::AtomData> currentResidueAtoms;
	int cysPair  = 0;
	while (true) {
		string line;
		getline(fs,line);

		if (fs.fail()) {
			// no more lines to read from file, quit the while loop
			break;
		}
		string header = "";
		if (line.size() >= PDBFormat::S_RECORD_NAME + PDBFormat::L_RECORD_NAME) {
		  header 	= line.substr(PDBFormat::S_RECORD_NAME, PDBFormat::L_RECORD_NAME);
		}

		if (header == "ATOM  " || header == "HETATM"){
		  PDBFormat::AtomData atom = PDBFormat::parseAtomLine(line,true);
		  
		  // Residue description string
		  stringstream resDescription;
		  resDescription << atom.D_CHAIN_ID <<":"<<atom.D_RES_NAME<<":"<<atom.D_RES_SEQ<<":"<<atom.D_I_CODE;
		  
		  // New residue means its time to decide which atoms will get added to "atoms"
		  if (currentResidue != "" && currentResidue != resDescription.str()){

		    AtomContainer *atoms = new AtomContainer();
		    vector<PDBFormat::AtomData>::iterator it;
		    for (it = currentResidueAtoms.begin();it != currentResidueAtoms.end();it++){
		      Atom a;
		      // atom name, residue name, residue icode, chain id, coor, element
		      a.setName(it->D_ATOM_NAME);
		      a.setResidueName(it->D_RES_NAME);
		      a.setResidueIcode(it->D_I_CODE);
		      a.setResidueNumber(it->D_RES_SEQ);
		      a.setChainId(it->D_CHAIN_ID);
		      a.setCoor(it->D_X,it->D_Y, it->D_Z);
		      a.setElement(it->D_ELEMENT_SYMBOL);
		      a.setTempFactor(it->D_TEMP_FACT);
		      a.setSegID(it->D_SEG_ID);

		      atoms->addAtom(a);
		    }

		    // Add atoms to pair.first
		    if (cysPair == 0){

		      _container.push_back(pair<AtomContainer *, AtomContainer *>(NULL,NULL));
		      _container.back().first = atoms;
		      cysPair++;
		      // Add atoms to pair.second
		    } else if (cysPair == 1){
		      _container.back().second = atoms;
		      //MSLOUT.fprintf(stdout,"Added new disulfide [ %10d ]\n", _container.size());
		      cysPair = 0;
		    } 


		    currentResidueAtoms.clear();
		  } // END NEW RESIDUE

		  currentResidueAtoms.push_back(atom);
		  currentResidue = resDescription.str();

		} // ATOM FIELD



	}
     } else {

       Timer t;
       double start = t.getWallTime();

       // READ A BINARY FILE
       AtomPointerVector ats;
       ats.load_checkpoint(_pdbfile);
       cout << "LOADED CHECKPOINT: "<<(t.getWallTime() - start)<<endl;
       fflush(stdout);

       AtomContainer *currentResidueAtoms = new AtomContainer();
       string currentResidue = "";
       int cysPair = 0;
       for (uint i = 0; i < ats.size();i++){
	 
	 
	 stringstream resDescription;
	 resDescription << ats(i).getPositionId();
		  
	 // New residue means its time to decide which atoms will get added to "atoms"
	 if (currentResidue != "" && currentResidue != resDescription.str()){
	   
	   // Add atoms to pair.first
	   if (cysPair == 0){

	     _container.push_back(pair<AtomContainer *, AtomContainer *>(NULL,NULL));
	     _container.back().first = currentResidueAtoms;
	     cysPair++;
	     // Add atoms to pair.second
	   } else if (cysPair == 1){
	     _container.back().second = currentResidueAtoms;
	     //MSLOUT.fprintf(stdout,"Added new disulfide [ %10d ]\n", _container.size());
	     cysPair = 0;
	   } 	   

	   currentResidueAtoms = new AtomContainer();
	 } // END CURRENT RESIDUE

	 currentResidueAtoms->addAtom(ats(i));
	 currentResidue = resDescription.str();


       } // For i in ats

       ats.deletePointers();

     }

}

//void readDisulfPdb(string _pdbfile,vector<pair<AtomContainer *, AtomContainer *> > &_container){
//  
//
//	ifstream fs;
//	fs.open(_pdbfile.c_str());
//	if (fs.fail()) {
//	  return;
//	}
//
//	string currentResidue = "";
//	vector<PDBFormat::AtomData> currentResidueAtoms;
//	int cysPair  = 0;
//	while (true) {
//		string line;
//		getline(fs,line);
//
//		if (fs.fail()) {
//			// no more lines to read from file, quit the while loop
//			break;
//		}
//		string header = "";
//		if (line.size() >= PDBFormat::S_RECORD_NAME + PDBFormat::L_RECORD_NAME) {
//		  header 	= line.substr(PDBFormat::S_RECORD_NAME, PDBFormat::L_RECORD_NAME);
//		}
//
//		if (header == "ATOM  " || header == "HETATM"){
//		  PDBFormat::AtomData atom = PDBFormat::parseAtomLine(line,true);
//		  
//		  // Residue description string
//		  stringstream resDescription;
//		  resDescription << atom.D_CHAIN_ID <<":"<<atom.D_RES_NAME<<":"<<atom.D_RES_SEQ<<":"<<atom.D_I_CODE;
//		  
//		  // New residue means its time to decide which atoms will get added to "atoms"
//		  if (currentResidue != "" && currentResidue != resDescription.str()){
//
//		    AtomContainer *atoms = new AtomContainer();
//		    vector<PDBFormat::AtomData>::iterator it;
//		    for (it = currentResidueAtoms.begin();it != currentResidueAtoms.end();it++){
//		      Atom a;
//		      // atom name, residue name, residue icode, chain id, coor, element
//		      a.setName(it->D_ATOM_NAME);
//		      a.setResidueName(it->D_RES_NAME);
//		      a.setResidueIcode(it->D_I_CODE);
//		      a.setResidueNumber(it->D_RES_SEQ);
//		      a.setChainId(it->D_CHAIN_ID);
//		      a.setCoor(it->D_X,it->D_Y, it->D_Z);
//		      a.setElement(it->D_ELEMENT_SYMBOL);
//		      a.setTempFactor(it->D_TEMP_FACT);
//
//		      atoms->addAtom(a);
//		    }
//
//		    // Add atoms to pair.first
//		    if (cysPair == 0){
//
//		      _container.push_back(pair<AtomContainer *, AtomContainer *>(NULL,NULL));
//		      _container.back().first = atoms;
//		      cysPair++;
//		      // Add atoms to pair.second
//		    } else if (cysPair == 1){
//		      _container.back().second = atoms;
//		      //MSLOUT.fprintf(stdout,"Added new disulfide [ %10d ]\n", _container.size());
//		      cysPair = 0;
//		    } 
//
//
//		    currentResidueAtoms.clear();
//		  } // END NEW RESIDUE
//
//		  currentResidueAtoms.push_back(atom);
//		  currentResidue = resDescription.str();
//
//		} // ATOM FIELD
//
//
//
//	}
//}
