#include "iqtree.h"

void printCopyright(ostream &out);
string copyright()
{
	stringstream stream;
	printCopyright(stream);
	return stream.str();
}
