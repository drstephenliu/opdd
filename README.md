# $O$PDD
A fast converging algorithm to solve the nonlinear problem in the following generalized form:

$$\left( \mathrm{OPDD} \right): \quad \min_{\mathbf{\pi}} \frac{1}{2} \Vert \mathbf{y} - \mathbf{B} \mathrm{exp} \left( -\mathbf{M} \mathbf{\pi} \right) \Vert_{\mathbf{\zeta}}^{2} + \sum_{\kappa} \beta_\kappa \Vert \pi_\kappa \Vert_{\mathrm{Huber}} \quad ,$$

here, with a particular emphasis on projection-domain decomposition problem in multi-energy radiography.

## Reference
The optimization framework of OPDD algorithm is a derivation from our paper below:
- Liu, S. Z., Tivnan, M., Osgood, G. M., Siewerdsen, J. H., Stayman, J. W., & Zbijewski, W. (2022) "Model-based three-material decomposition in dual-energy CT using the volume conservation constraint," *Phys. Med. Biol.*, **67**(14), 145006. DOI: https://doi.org/10.1088/1361-6560/ac7a8b

## Contact
- Wojciech Zbijewski: wzbijewski@jhu.edu
- Stephen Z. Liu: szliu@jhmi.edu
