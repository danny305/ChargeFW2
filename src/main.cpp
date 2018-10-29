#include "formats/SDF.h"
#include "structures/MoleculeSet.h"
#include "PeriodicTable.h"

int main() {
    SDF reader;
    MoleculeSet m = reader.read_file("/home/krab1k/Tmp/HOH.sdf");
    m.info();
    return 0;
}