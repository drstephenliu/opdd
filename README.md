# OPDD
A fast converging algorithm to solve the nonlinear problem in the following generalized form:

$$\left( \mathrm{OPDD} \right): \quad \min_{\mathbf{\pi}} \frac{1}{2} \Vert \mathbf{y} - \mathbf{B} \mathrm{exp} \left( -\mathbf{M} \mathbf{\pi} \right) \Vert_{\mathbf{\zeta}}^{2} + \mathbf{\beta}^T \mathcal{R}_{\mathrm{Huber}} \left( \mathbf{\pi} \right) \quad ,$$

with smoothness regularization, optimal line search and momentum. Here, we show a particular emphasis on the projection-domain decomposition problem in multi-energy x-ray radiography/CT/CBCT imaging.

## Reference
The optimization framework for OPDD is a derivation from our paper below:
- Liu, S. Z., Tivnan, M., Osgood, G. M., Siewerdsen, J. H., Stayman, J. W., & Zbijewski, W. (2022) "Model-based three-material decomposition in dual-energy CT using the volume conservation constraint," *Phys. Med. Biol.*, **67**(14), 145006. DOI: https://doi.org/10.1088/1361-6560/ac7a8b

## Contact
- Stephen Z. Liu: szliu@jhmi.edu
- Wojciech Zbijewski: wzbijewski@jhu.edu

## Description
We here present the definition of each matrix/vector in the context of multi-energy imaging. However, keep in mind that all terms in OPDD can be customized for general purposes.

(Under construction...)
