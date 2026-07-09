# Sums of powers via central finite differences and Newton's formula

## Abstract

In this manuscript, we derive closed formulas for multifold sums of powers of integers by combining the central Newton interpolation formula with hockey-stick identities for binomial coefficients. We further provide Wolfram Mathematica programs for the efficient verification of the derived identities.


## Related projects

- [Newton's interpolation formula and sums of powers (2025)](https://github.com/kolosovpetro/NewtonsInterpolationFormulaAndSumsOfPowers)
- [Sums of powers via central finite differences and Newton's formula (2025)](https://github.com/kolosovpetro/SumsOfPowersViaCentralFiniteDifferencesAndNewtonFormula)
- [Sums of powers via backward finite differences and Newton's formula (2026)](https://github.com/kolosovpetro/SumsOfPowersViaBackwardFiniteDifferencesAndNewtonFormula)
- [Sums of powers of integers: A complete framework for closed formulas (2026)](https://github.com/kolosovpetro/SumsOfPowersACompleteFrameworkForClosedForms)

## OEIS

### Central factorial numbers

- https://oeis.org/A036969 — Triangle read by rows: T(n,k) = T(n-1,k-1) + k^2*T(n-1,k), 1 < k <= n, T(n,1) = 1.
- https://oeis.org/A269945 — Triangle read by rows. Stirling set numbers of order 2, T(n, n) = 1, T(n, k) = 0 if k < 0 or k > n, otherwise T(n, k) = T(n-1, k-1) + k^2*T(n-1, k), for 0 <= k <= n.
- https://oeis.org/A008957 — Triangle of central factorial numbers T(2*n,2*n-2*k), k >= 0, n >= 1 (in Riordan's notation).
- https://oeis.org/A395862 — Triangle read by rows: T(n,k) = numerator(CF(n,k)), where CF(n,k) = (1/k!) * Sum_{j=0..k} (-1)^j * binomial(k,j) * (k/2-j)^n.
- https://oeis.org/A370703 — Triangle read by rows: T(n, k) = denominator([x^k] n! [t^n] (t/2 + sqrt(1 + (t/2)^2))^(2*x)).
- https://oeis.org/A395860 — Triangle read by rows: T(n,k) = numerator(CF(n,k)), where CF(n,k) = (1/k!) * Sum_{j=0..k} (-1)^j * binomial(k,j) * (1+k/2-j)^n.
- https://oeis.org/A395861 — Triangle read by rows: T(n,k) = denominator(CF(n,k)), where CF(n,k) = (1/k!) * Sum_{j=0..k} (-1)^j * binomial(k,j) * (1+k/2-j)^n.
- https://oeis.org/A394466 — Triangle read by rows: T(n,k) = numerator(CF(n,k)), where CF(n,k) = (1/k!) * Sum_{j=0..k} (-1)^j * binomial(k,j) * (2+k/2-j)^n.
- https://oeis.org/A395314 — Triangle read by rows: T(n,k) = denominator(CF(n,k)), where CF(n,k) = (1/k!) * Sum_{j=0..k} (-1)^j * binomial(k,j) * (2+k/2-j)^n.
- https://oeis.org/A396559 — Triangle read by rows: T(n,k) = numerator(CF(n,k)), where CF(n,k) = (1/k!) * Sum_{j=0..k} (-1)^j * binomial(k,j) * (3+k/2-j)^n.
- https://oeis.org/A269945 — Triangle read by rows. Stirling set numbers of order 2, T(n, n) = 1, T(n, k) = 0 if k < 0 or k > n, otherwise T(n, k) = T(n-1, k-1) + k^2*T(n-1, k), for 0 <= k <= n.
- https://oeis.org/A394692 — Triangle read by rows: T(n,k) = (1/(2k)!) * Sum_{j=0..2k} (-1)^j * binomial(2k,j) * (1+k-j)^(2n).
- https://oeis.org/A395456 — Triangle read by rows: T(n,k) = (1/(2k)!) * Sum_{j=0..2k} (-1)^j * binomial(2k,j) * (2+k-j)^(2n).
- https://oeis.org/A395457 — Triangle read by rows: T(n,k) = (1/(2k)!) * Sum_{j=0..2k} (-1)^j * binomial(2k,j) * (3+k-j)^(2n).

### Central finite differences

- https://oeis.org/A387597 — Triangle read by rows: T(n,k) = Sum_{j=0..2k} (-1)^j * binomial(2k,j) * (0+k-j)^(2n).
- https://oeis.org/A392337 — Triangle read by rows: T(n,k) = Sum_{j=0..2k} (-1)^j * binomial(2k,j) * (1+k-j)^(2n).
- https://oeis.org/A390029 — Triangle read by rows: T(n,k) = Sum_{j=0..2k} (-1)^j * binomial(2k,j) * (2+k-j)^(2n).

## Metadata

- **Initial release date:** January 3, 2026.
- **MSC2010:** 05A19, 05A10, 41A15, 11B83.
- **Keywords:** Sums of powers, Newton's interpolation formula, Finite differences, Binomial coefficients, Faulhaber's formula, Bernoulli numbers, Bernoulli polynomials, Interpolation, Approximation, Discrete convolution, Combinatorics, Polynomial identities, Central factorial numbers, Stirling numbers, Eulerian numbers, Worpitzky identity, Pascal's triangle, OEIS.
- **License:** This work is licensed under a [CC BY 4.0 License](https://creativecommons.org/licenses/by/4.0/).
- **DOI:** https://doi.org/10.5281/zenodo.18096789
- **Web Version:** https://kolosovpetro.github.io/sums-of-powers-central-differences/
- **Sources:** https://github.com/kolosovpetro/SumsOfPowersViaCentralFiniteDifferencesAndNewtonFormula
- **ORCID:** https://orcid.org/0000-0002-6544-8880
- **Email:** kolosovp94@gmail.com

## References

- Knuth, D. E. (1993). Johann Faulhaber and sums of powers. Mathematics of Computation, 61(203), 277–294. https://arxiv.org/abs/math/9207222
- Butzer, P. L., Schmidt, K., Stark, E. L., & Vogt, L. (1989). Central factorial numbers; their main properties and some applications. Numerical Functional Analysis and Optimization, 10(5–6), 419–488. https://doi.org/10.1080/01630568908816313
- Steffensen, J. F. (1933). On the definition of the central factorial. Journal of the Institute of Actuaries (1886–1994), 64(2), 165–168. https://www.jstor.org/stable/41137516
- Steffensen, J. F. (1927). Interpolation. Williams & Wilkins. https://www.amazon.com/-/de/Interpolation-Second-Dover-Books-Mathematics-ebook/dp/B00GHQVON8
- Riordan, J. (1968). Combinatorial identities. Wiley, New York. https://www.amazon.com/-/de/Combinatorial-Identities-Probability-Mathematical-Statistics/dp/0471722758
- Graham, R. L., Knuth, D. E., & Patashnik, O. (1994). Concrete Mathematics: A Foundation for Computer Science (2nd ed.). Addison-Wesley Publishing Company, Inc. https://archive.org/details/concrete-mathematics
- Knuth, D. E. (1992). Two notes on notation. https://arxiv.org/abs/math/9205211
- Kolosov, P. (2026). Sums of powers of integers: A complete framework for closed formulas. https://doi.org/10.5281/zenodo.20548019
- Kolosov, P. (2025). Newton's interpolation formula and sums of powers. https://doi.org/10.5281/zenodo.18040979
- Kolosov, P. (2026). Sums of powers via backward finite differences and Newton's formula. https://doi.org/10.5281/zenodo.18118011