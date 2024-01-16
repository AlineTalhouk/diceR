## Test environments
* local R installation, R 4.3.2
* ubuntu 22.04 (on GitHub Actions), R 4.3.2, devel
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

Suggests or Enhances not in mainstream repositories:
  mixedClust
  
I updated to skip certain unit tests if Suggested packages are unavailable. In the case of mixedClust, it is archived.
