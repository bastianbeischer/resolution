#include "Analysis.hh"

int main(int argc, char* argv[])
{
  Analysis ana;

  if (argc > 1)
    ana.processFile(argv[1]);

  return 0;
}
