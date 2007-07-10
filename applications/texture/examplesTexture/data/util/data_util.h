/**
 * @defgroup data_util Texture: Utility functions for preprocessing data.
 * @ingroup texture_examples
 * @author Matthias Schmalz
 * @{
 */

/**
 * Calculates the presentation of an angle in @f$[-\pi, \pi)@f$.
 */
double normalise_phi(double angle);

/**
 * Calculates the presentation of an elevation angle in @f$[0, \pi]@f$.
 * Raises an assertion failure if this is not possible.
 */
double normalise_theta(double angle);

/**
 * Converts the degree presentation of an angle to the radiant presentation.
 */
double deg2rad(double angle);

/**
 * Calculates the mirrored version of an azimuth angle.
 * @f$\phi \mapsto \phi + PI@f$.
 */
double invert_phi(double phi);

/**
 * Calculates the mirrored version of an elevation angle.
 * @f$\theta \mapsto \pi - \theta@f$.
 */
double invert_theta(double theta);

/**
 * Applies ::normalise_phi to each member of h_phi and ::normalise_theta
 * to each member of ::h_theta.
 * @par N1 is the length of h_phi and h_theta.
 */
void normalise_h(int N1, double *h_phi, double *h_theta);

/**
 * Applies ::normalise_phi and ::normalise_theta to the angles in r.
 * @par N2 is the number of points in r.
 * @par r stores the points. If the i-th point is @f$( \phi, \theta ) @f$, 
 * then
 * @f$r(2 i) = \phi@f$ and @f$r(2 i + 1) = \theta@f$.
 */
void normalise_r(int N2, double *r);

/**
 * Mirrores the nodes from r using ::invert_phi and ::invert_theta.
 * @par N2 is the number of points in r after mirroring.
 * @par r stores the points. If the i-th point is @f$(\phi, \theta)@f$, then
 * @f$r(2 i) = \phi@f$ and @f$r(2 i + 1) = \theta@f$.
 */
void expand_r(int N2, double *r);

/**
 * Mirrores the pole figures from h_phi and h_theta using ::invert_phi and 
 * ::invert_theta.
 * @par N1 is the number of points in h_phi / h_theta after mirroring.
 * @par h_phi and h_theta store the points. 
 */
void expand_h(int N1, double *h_phi, double *h_theta);
/**
 * @}
 */
