#include "kraken2_headers.h"
#include "mmscanner.h"

using namespace std;
using namespace kraken2;

int main() {
  string seq = "ACGATCGACGACG";
  MinimizerScanner scanner(10, 5, 0);
  scanner.LoadSequence(seq);
  uint64_t *mmp;
  while ((mmp = scanner.NextMinimizer()) != nullptr)
    printf("%016x\n", *mmp);
  return 0;
}
