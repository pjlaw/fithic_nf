# fithic_nf

This is effectively a fork of [FitHiC2](https://github.com/ay-lab/fithic), converted to a nextflow pipeline, using the guidelines described in the [Nature Protocols manuscript](https://www.nature.com/articles/s41596-019-0273-0).

Main changes:
 - Uses [JuicerTools](https://github.com/aidenlab/juicer/wiki/Download) to dump the interactions, and perform the normalisation
 - Modified the CombineNearbyInteraction.py script to pre-cluster neighbouring bins, vastly improving the run time for high resolution data

#### #TODO
 - make Docker/singularity of fithic rather than conda
 - Update nextflow to DSL2
