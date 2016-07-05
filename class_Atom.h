#include <string>
using namespace std;
class Atom
{
private:

protected:

public:

	static const int spin_dim = 3;

	string species;             //Species of a atom,refer to the POSCAR file

	int atom_label;//Label of the atom

	double position[3]; //Atom position in the xyz space

	double spin[spin_dim];

	double On_Site;
	//	int tag;
};



