from dolfin import *
code = \
"""
#include "dolfin.h"

namespace dolfin{

void allow_nonzero_allocation(GenericMatrix* A)
{
  #ifdef HAS_PETSC
  if (A && has_type<PETScMatrix>(*A))
  {
    PETScMatrix& petsc_A = A->down_cast<PETScMatrix>();
    // what?
    MatSetOption(petsc_A.mat(), MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  }
  #endif
}
}
"""
compiled_module = compile_extension_module(code=code)
