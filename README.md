# HybridDIA
HybridDIA is a new pethod thatcombines simultaneous targeted and discovery (DIA) proteomics methods. HybridDIA is based on an Application Programming Interface, which enables the ability to combine trgeted and discovery proteomics by intercalating DIA scans with accurate triggering of predefined peptide targets.

The pipeline to analyze HybridDIA files is comprised of three steps:
1. Extract MSx Scans from the raw file.
2. Upload to [Skyline](https://skyline.ms/project/home/software/Skyline/begin.view) to extract the Fragment Intensities of the Heavy (Internal Standard or IS) and Light (Endogeneous or ENDO) peptides.
3. Injection Time based intensity correction and data visualization.

You can find a complete guide in: [Guide]("anamdv/HybridDIA/Guide")
