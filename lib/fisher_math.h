#ifndef TWK_FISHER_MATH_H_
#define TWK_FISHER_MATH_H_


/* Log gamma function
 * \log{\Gamma(z)}
 * AS245, 2nd algorithm, http://lib.stat.cmu.edu/apstat/245
 */
double kf_lgamma(double z);

/* complementary error function
 * \frac{2}{\sqrt{\pi}} \int_x^{\infty} e^{-t^2} dt
 * AS66, 2nd algorithm, http://lib.stat.cmu.edu/apstat/66
 */
double kf_erfc(double x);

/* The following computes regularized incomplete gamma functions.
 * Formulas are taken from Wiki, with additional input from Numerical
 * Recipes in C (for modified Lentz's algorithm) and AS245
 * (http://lib.stat.cmu.edu/apstat/245).
 *
 * A good online calculator is available at:
 *
 *   http://www.danielsoper.com/statcalc/calc23.aspx
 *
 * It calculates upper incomplete gamma function, which equals
 * kf_gammaq(s,z)*tgamma(s).
 */

double kf_gammap(double s, double z);
double kf_gammaq(double s, double z);

/* Regularized incomplete beta function. The method is taken from
 * Numerical Recipe in C, 2nd edition, section 6.4. The following web
 * page calculates the incomplete beta function, which equals
 * kf_betai(a,b,x) * gamma(a) * gamma(b) / gamma(a+b):
 *
 *   http://www.danielsoper.com/statcalc/calc36.aspx
 */
double kf_betai(double a, double b, double x);

/*
 *    n11  n12  | n1_
 *    n21  n22  | n2_
 *   -----------+----
 *    n_1  n_2  | n
 */
double kt_fisher_exact(int n11, int n12, int n21, int n22, double *_left, double *_right, double *two);

inline const double chi_squared(const double& n11, const double& n12, const double& n21, const double& n22){
	const double marginA_R = n11 + n12;
	const double marginB_R = n21 + n22;
	const double marginA_B = n11 + n21;
	const double marginB_B = n12 + n22;
	const double n_total = n11 + n12 + n21 + n22;

	// Expected
	return(pow(n11 - (marginA_R*marginA_B / n_total), 2) / (marginA_R*marginA_B / n_total) +
	       pow(n12 - (marginA_R*marginB_B / n_total), 2) / (marginA_R*marginB_B / n_total) +
	       pow(n21 - (marginB_R*marginA_B / n_total), 2) / (marginB_R*marginA_B / n_total) +
	       pow(n22 - (marginB_R*marginB_B / n_total), 2) / (marginB_R*marginB_B / n_total));
}

#endif /* TWK_FISHER_MATH_H_ */
