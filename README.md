# COULCC-25 v1.0
!-----------------------------------------------------------------------

  Copyright (c) 1983, 2007, 2012 A. R Barnett and I. J. Thompson
  Permission is hereby granted, free of charge,
  to any person obtaining a copy of this software and associated documentation
  files (the "Software"), to deal in the Software without restriction, including without
  limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
  and/or sell copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:


  The above copyright notice and this permission notice shall be included in all copies or
  substantial portions of the Software.


  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
  INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
  PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
  OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
  OTHER DEALINGS IN THE SOFTWARE.

!-----------------------------------------------------------------------

   COMPLEX COULOMB WAVEFUNCTION PROGRAM USING STEED'S METHOD

!-----------------------------------------------------------------------

  A CONTINUED-FRACTION ALGORITHM FOR COULOMB

  FUNCTIONS OF COMPLEX ORDER WITH COMPLEX ARGUMENTS

!-----------------------------------------------------------------------

  REF. IN COMP. PHYS. COMMUN. 36 (1985) 363

!-----------------------------------------------------------------------

  A. R. Barnett           Manchester  March   1981

  modified I.J. Thompson  Daresbury, Sept. 1983 for Complex Functions

  modified A.R. flores    WashU, April 2025 for Modern Compilers

!-----------------------------------------------------------------------

  original program  RCWFN       in    CPC  8 (1974) 377-395

  original program  RCWFF       in    CPC 11 (1976) 141-142

  original program  COULFG      in    CPC 27 (1982) 147-166

  description of real algorithm in    CPC 21 (1981) 297-314

  description of complex algorithm    JCP XX (1985) YYY-ZZZ

  this version written up       in    CPC XX (1985) YYY-ZZZ

# Building and testing
## fpm
### GCC (`gfortran`)
```
fpm test --compiler gfortran --profile release --flag "-cpp"
```
## cmake
