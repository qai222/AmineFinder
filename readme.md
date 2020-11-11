# AmineFinder
This was developed to automatically find valid amines from
a long IUPAC name in the Cambridge Structure Database, e.g. 
```
IUPAC_name = "catena-(1,8-diaminooctane bis(bis(μ4-phosphato)-(μ2-hydroxo)-tri-tin))"
amine_name = "1,8-diaminooctane"
```
This code heavily relies on the `OPSIN` package:
> Lowe, Daniel M., et al. "Chemical name to structure: OPSIN, an open source solution." (2011): 739-753.

## install
1. download [opsin](https://github.com/dan2097/opsin)
   ```
   wget https://github.com/dan2097/opsin/releases/download/2.5.0/opsin-2.5.0-jar-with-dependencies.jar -O opsin-2.5.0.jar
   ```
2. install packages:
   ```
   pip install python-Levenshtein pyyaml    
   ```
   or using the [requirements.txt](requirements.txt)

## usage
See [test.py](test.py) for examples. In general the `MoleculeFinder` class can
be used to identify other types of molecule so long as a set of keywords can be defined. 
See the doc string in [mfinder.py](mfinder.py) for more info.

