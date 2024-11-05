#include <string>
#include <vector>
using namespace std;


class Graph 
{ 
    int V;    // No. of vertices 
    list<int> *adj;    // Pointer to an array containing adjacency lists 
    bool rootConnect(int last,int v, int root, bool visited[], vector<int> path); // Recursive function to check whether root connects back to root
    bool rootBarrel(int root); // Whether a barrel exist from current root
    bool seeRoot(int v, int root);


public: 
    Graph(int V);   // Constructor 
    void addEdge(int v, int w);   // to add an edge to graph 
    bool isBarrel(); // Return true is is_barrel
    
};

// Define Input Options and Store them.
struct Options {

    // Set up options here...
    Options(){

        // Input pdb
        required.push_back("pdb");

        defaultArgs.push_back("configfile");

    }

    // Storage for the vales of each optional
    string pdb;
    bool debug;

    // Storage for different types of options
    vector<string> required;
    vector<string> optional;
    vector<string> defaultArgs;


};


// Helper function to clean up main.
Options setupOptions(int theArgc, char * theArgv[]);
int getShortestStrand(System &_sys, vector< pair<int,int> > &_strands);
AtomPointerVector getTopAtoms(System &_sys, vector< pair<int,int> > &_strands, int _min_strand);
AtomPointerVector getBotAtoms(System &_sys, vector< pair<int,int> > &_strands, int _min_strand);
bool checkBetaSheet(System &_sys, int _pos);
bool neighboringStrand(System &_sys, pair<int,int> &_s1, pair<int,int> &_s2);
