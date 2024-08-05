# Project Overview #
**Vegetation Interaction with Sand Movement in Coastal Dynamics**

Coastal environments are dynamic systems where natural processes move sediment to build landforms. The primary sources of energy to move sediment are wind and water driven by current, waves, tides, storm events, temperature and density stratification. Coastal vegetation plays a pivotal role in initiating planform changes, such as dune and bar formation. These morphological alterations arise from various dynamic processes, including aerodynamics and hydrodynamics. Aerodynamics explains the mechanisms through which wind transports sand particles, with vegetation acting as an obstruction that reduces wind velocity, leading to the deposition of transported sand and the subsequent formation of dunes. Conversely, hydrodynamics describes the processes governing sediment transport, the interaction with aquatic vegetation, and the deposition of sediments, ultimately leading to bar formation. Although dune and bar formation processes may appear similar, the mechanisms and factors involved in sediment and sand transport are markedly different, necessitating the use of distinct models to simulate these dynamics. We employed Aeolis to model the aerodynamic processes of dune formation, while the PyDeltaRCM morphodynamic model was used to simulate the hydrodynamic processes involved in bar formation.




## Overview of AeoLiS ## 
**Model:** Aeolis is an open source, 2D depth average aeolian sediment transport model
**Developer:** Bas M. Hoonhout, Sierd de Vries
AeoLiS is a process-based model designed for simulating aeolian sediment transport, particularly in environments where supply-limiting factors are significant, such as coastal areas. It can simulate various surface configurations, including moisture, shells, vegetation, and non-erodible elements. The bed surface configuration influences aeolian sediment transport by altering both the sediment transport capacity and sediment availability.
**Reference:** https://doi.org/10.1002/2015JF003692 



## Overview of PyDeltaRCM ##
**Model:** PyDeltaRCM is an open source, 2D depth average morphodynamic model. 
**Developer:** Mariela Perignon (Python version of the Matlab deltaRCM model by Man Liang). 
PyDeltaRCM is a parcel-based cellular flux routing and sediment transport model that can simulate the planform evaluation due to sediment movement and vegetation model. It is a  reduced-complexity model (RCM), different from the process-based morphodynamic models based on detailed computational fluid dynamics by employing stochastic parcel-based cellular routing schemes for water and sediment transport. The model requires different parameters and variables as input including topography, flow discharge, bed roughness and produces outputs including depth-averaged flow field, water surface elevation and evolving bed topography. 
**Reference:** https://doi.org/10.21105/joss.03398  




