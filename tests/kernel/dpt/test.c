#include "nfft3.h"
#include "dpt.h"

#include <stdlib.h>
#include <stdio.h>

int main()
{
  dpt_set dpts = dpt_init(1,5,0U);
  dpt_finalize(dpts);
  return EXIT_SUCCESS;
}