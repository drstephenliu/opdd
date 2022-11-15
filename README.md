# OPDD
A fast converging algorithm to solve the nonlinear problem in the following generalized form:

$$\min_{\mathbf{\pi}} \frac{1}{2} \Vert \mathbf{y} - \mathbf{B} \mathrm{exp} \left( -\mathbf{M} \mathbf{\pi} \right) \Vert_{\mathbf{z}}^{2} + \sum_{\kappa} \beta_\kappa \Vert \pi_\kappa \Vert_{\mathrm{Huber}}$$
