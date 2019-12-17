#! /bin/bash

if ! python -c 'import numpydoc'; then easy_install --user numpydoc; fi
if ! python -c 'import sphinx'; then easy_install --user sphinx; fi

outdir=source/api
sphinx-apidoc -H "Comprehensive API reference" -M -f -o $outdir ../Corrfunc/ ../Corrfunc/tests.py ../Corrfunc/call_correlation_functions.py ../Corrfunc/call_correlation_functions_mocks.py

# Fix the blank sub-modules in the Corrfunc file
for docfile in $outdir/Corrfunc.rst
do
    # Delete three lines following the "submodules"
    sed -e '/Submodules/{N;N;d;}' $docfile > xx
    mv xx  $docfile
done


# Fix the duplicate entries for the various pair-counters
# (e.g., Corrfunc.theory.DD *and* Corrfunc.theory.DD.DD)
for docfile in $outdir/Corrfunc.mocks.rst $outdir/Corrfunc.theory.rst
do
    # Delete ALL lines following this "submodule" line in the theory/mocks
    # auto-generated documentation
    sed  -n '/Submodules/q;p' $docfile > xx
    mv xx  $docfile
done
