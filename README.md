CalaTK
==

The CalaTK project develops innovative methods and tools for longitudinal atlases with a focus on neurodevelopment.  
The computational toolbox is developed with the objective to analyze the neural developmental patterns observed in macaque structural and diffusion  on LDDMM, growth model LDDMM, LDDMM with geodesic shooting and/or initial momentum, or geometric metamorphostensor magnetic resonance (MR) images.

A number of algorithms are available including registration and atlas building basedis with LDDMM.

Unlike existing atlas­-building methods, we explicitly use longitudinal (or temporal) information, both for the structural atlas as well as for the diffusion tensor atlas.
This will be achieved directly within the registration framework by modeling expected changes in image intensity for the structural images (to handle contrast inversion at the early stage of brain development) and by using subject-­specific growth models. 
The proposed atlas­-building strategy is specifically tailored for the construction of longitudinal atlases. 
The longitudinal approach is expected to significantly improve estimation accuracy.

CalaTK is implemented as a C++ library and makes extensive use of ITK.  Algorithms are provided as command-line applications and as extensions in 3D Slicer, an open-source medical imaging application.  Algorithm configuration is performed with a system based on JSON files that is human-readable, machine-editable, and easily archived for reproducible analysis.
