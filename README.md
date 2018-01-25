# ATACseq-workflow

In house script for processing ATACseq data, AFTER peakcalling

To do:
* Neaten up code to improve useability in the future
* Change normalisation strategy - Currently uses CQN for coverage data, and resampling for tracks. Replace both with deeptools normalisation?
* Optimise various stages

General notes:
* Blacklist file can easily be modified, and has been altered for Col1a1 spikes
