# Evaluation of alignment libraries

This is a selfcontained benchmark to test which of the many existing C/C++ libraries for sequence alignment work best for our use case: alinigng short reads with high error rates (in particular ancient DNA damage at the read ends) to genome subsequences that span the read length plus some extra intervals outside, of some hundred base pairs max in each direction of the read ends.

Compile:

```
cd spear-map-paper/eval-align-libs
cmake -B build && cmake --build build
```
