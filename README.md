# Sums of Powers via Central Finite Differences and Newton's Formula

## Abstract

In this manuscript, we derive closed formulas for multifold sums of powers of integers by combining the central Newton interpolation formula with hockey-stick identities for binomial coefficients. We further provide Wolfram Mathematica programs for the efficient verification of the derived identities.

## Related projects

- [Newton's interpolation formula and sums of powers (2025)](https://github.com/kolosovpetro/NewtonsInterpolationFormulaAndSumsOfPowers)
- [Sums of Powers via Central Finite Differences and Newton's Formula (2025)](https://github.com/kolosovpetro/SumsOfPowersViaCentralFiniteDifferencesAndNewtonFormula)
- [Sums of powers via backward finite differences and Newton's formula (2026)](https://github.com/kolosovpetro/SumsOfPowersViaBackwardFiniteDifferencesAndNewtonFormula)
- [Sums of powers of integers: A complete framework for closed formulas (2026)](https://github.com/kolosovpetro/SumsOfPowersACompleteFrameworkForClosedForms)

## References

1. Knuth, D. E. (1993). Johann Faulhaber and sums of powers. *Mathematics of Computation, 61*(203), 277–294. https://arxiv.org/abs/math/9207222
2. Steffensen, J. F. (1933). On the Definition of the Central Factorial. *Journal of the Institute of Actuaries (1886–1994), 64*(2), 165–168. https://www.jstor.org/stable/41137516
3. Butzer, P. L., Schmidt, K., Stark, E. L., & Vogt, L. (1989). Central Factorial Numbers: Their Main Properties and Some Applications. *Numerical Functional Analysis and Optimization, 10*(5–6), 419–488. https://doi.org/10.1080/01630568908816313
4. Steffensen, J. F. (1927). *Interpolation*. Williams & Wilkins. https://www.amazon.com/-/de/Interpolation-Second-Dover-Books-Mathematics-ebook/dp/B00GHQVON8
5. Riordan, J. (1968). *Combinatorial Identities* (Vol. 217). Wiley, New York. https://www.amazon.com/-/de/Combinatorial-Identities-Probability-Mathematical-Statistics/dp/0471722758
6. Graham, R. L., Knuth, D. E., & Patashnik, O. (1994). *Concrete Mathematics: A Foundation for Computer Science* (2nd ed.). Addison-Wesley. https://archive.org/details/concrete-mathematics
7. Knuth, D. E. (1992). Two Notes on Notation. arXiv:math/9205211. https://arxiv.org/abs/math/9205211
8. Kolosov, P. (2026, June). *Sums of Powers of Integers: A Complete Framework for Closed Formulas*. Zenodo. https://doi.org/10.5281/zenodo.20548019
9. Kolosov, P. (2025, December). *Newton's Interpolation Formula and Sums of Powers*. Zenodo. https://doi.org/10.5281/zenodo.18040979
10. Kolosov, P. (2026, January). *Sums of Powers via Backward Finite Differences and Newton's Formula*. Zenodo. https://doi.org/10.5281/zenodo.18118011
11. Sloane, N. J. A., et al. (2003). *The On-Line Encyclopedia of Integer Sequences*. OEIS. https://oeis.org/

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
- **Keywords:** Sums of powers, power sums, multifold sums, central finite differences, Newton interpolation formula, binomial sums, central factorial numbers, Stirling numbers, Eulerian numbers, Faulhaber's formula.
- **License:** This work is licensed under a [CC BY 4.0 License](https://creativecommons.org/licenses/by/4.0/).
- **DOI:** https://doi.org/10.5281/zenodo.18096789
- **Web Version:** https://kolosovpetro.github.io/sums-of-powers-central-differences/
- **Sources:** https://github.com/kolosovpetro/SumsOfPowersViaCentralFiniteDifferencesAndNewtonFormula
- **ORCID:** https://orcid.org/0000-0002-6544-8880
- **Email:** kolosovp94@gmail.com
