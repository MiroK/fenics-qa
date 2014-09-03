#!/usr/bin/env python
from ufl import *
set_level(DEBUG)
A = Constant(element)*dx(mesh)
