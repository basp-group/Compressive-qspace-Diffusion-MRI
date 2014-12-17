CO4dMRI: Convex Optimisation for diffusion MRI imaging
                
  ----------------------------------------------------------------

DESCRIPTION
  This is a Matlab implementation of the optimisation algorithms for
  diffusion MRI imaging presented in:
    [1] A. Auria, A. Daducci, J.-P. Thiran and Y. Wiaux.  
    "Structured sparsity for spatially coherent fibre orientation estimation in diffusion MRI", submitted to NeuroImage, 2014.

AUTHORS
  A. Auria (http://people.epfl.ch/anna.auria)
  A. Daducci (http://people.epfl.ch/alessandro.daducci)
  J.-P. Thiran (http://people.epfl.ch/jean-philippe.thiran)
  Y. Wiaux (http://people.epfl.ch/yves.wiaux)

EXPERIMENTS
  To test a reconstruction using the algorithm L2L0_NL,
  run the script Demo_L2L0nl.m. 

  
REFERENCES
  When referencing this code, please cite our related paper:
    [1] A. Auria, A. Daducci, J.-P. Thiran and Y. Wiaux.  
    "Structured sparsity for spatially coherent fibre orientation estimation in diffusion MRI", submitted 
    to NeuroImage, 2014.

INSTALLATION 
  To run the main solver, one of the functions of SPGL1 Toolbox (Spectral Projected Gradient for L1, https://www.math.ucdavis.edu/~mpf/spgl1/) is used. For the sake of completeness, it is already included in subfolder "spgl1-1.7". 
Otherwise, the code should run as is, without any additional installation.

DOWNLOAD
  https://github.com/basp-group/co-dmri

SUPPORT
  If you have any questions or comments, feel free to contact Anna
  Auria at: anna {DOT} auria {AT} epfl {DOT} ch.

NOTES
  The code is not optimized and is given for educational purpose.
  
LICENSE
  co4dmri: Convex optimization for optical interferometry
  Copyright (C) 2014 Anna Auria, Alessandro Daducci, Jean-Philippe
  Thiran and Yves Wiaux.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details (LICENSE.txt).

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301, USA.
