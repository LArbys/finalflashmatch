# Isolate Flash Match Tuning sample

We take in a list of runs and files and produce one file with the necessary ingredients for the flash matching code.

This involves bringing together the following items:

### From LArLite

* vertex: selected vertex information from (stage2) vertex file
* ophits/opflash: optical flash information, opreco and opflash, from (stage1) opreco file
* tracks: tracker output from 3D reco
* hits: gaushits from (stage1) reco2d file


## Building

Can be a little tricky. Build dllee_unified basically. Then go into LLCVProcessor and build the core. Go into app and InterTool and build there.

These are LArCV modules, so go into LArCV top directory and type make to copy symbols into larlitecv.so library.