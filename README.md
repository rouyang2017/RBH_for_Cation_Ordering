- Recommendation-based Basin-Hopping (RBH) for thermodynamic global optimization of cation ordering in multi-cation perovskites
Reference: Y. Zhang, Z. Li, Z.-K. Han, R. Ouyang, Global optimization of cation ordering in perovskites by recommendation-based basin-hopping, under review (2024).

- Modification to the traditional BH:   
a) A recommender system based on an on-the-fly machine learning force field (MLFF) as implemented in VASP is used to guide the Monte Carlo move.   
b) A scheme that combines DFT and on-the-fly MLFF is used to ensure fast and accurate minimization toward local minima.  

- Three examples of using the RBH for cation ordering optimization are provided in the folders BSCFbulk, BSCFsurface, and CSCFsurface.

- Running RBH   
1. Setting your parameters in the program RBH.f90, and compile it: ifort RBH.f90 -o RBH
2. Places the files INCAR_0, INCAR_1, INCAR_2, INCAR_3, POSCAR, KPOINTS, POTCAR, and the code RBH in the same working directory
3. Put the command ./RBH in your job submission script and submit your job.   
