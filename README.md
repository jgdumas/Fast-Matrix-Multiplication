--------------------------------------------------------------------------------
# Fast-Matrix-Multiplication
### Accuracy of Fast matrix multiplication Algorithms via 2x2 recursion
--------------------------------------------------------------------------------

**Authors**:  Jean-Guillaume Dumas, Clément Pernet, Alexandre Sedoglavic
- [ J-G. Dumas, C. Pernet, A. Sedoglavic; Strassen's algorithm is not optimally accurate, ISSAC 2024, Raleigh, NC USA, pp. 254-263.](https://hal.science/hal-04441653)


**About**:
This is a Fork of the
[Complex-Matrix-Multiplication](https://github.com/zhen06/Complex-Matrix-Multiplication) repository,
created for the numerical experiments of the paper:
- [Dai, Z., Lim, LH; Numerical stability and tensor nuclear norm. Numer. Math. 155, 345--376 (2023)](https://link.springer.com/article/10.1007/s00211-023-01377-5).

- We add more accurate variants of fast matrix multiplication via 2x2 recursion



--------------------------------------------------------------------------------
**`FMM-matlab-benchmarks`, Accuracy benchmarks with Matlab**:
| Program | Growth Factor | Description |
| :---    |     :---:     |        ---: |
| `strassenw.m` 	| 17.8530 | original Winograd's algorithm|
| `strassen.m` 		| 14.8284 | original Strassen's algorithm|
| `DPS_evenpow.m` 	| 12.2034 | using only powers of 2|
| `DPS_smallrat.m` 	| 12.2034 | small rational coefficients|
| `DPS_intermediate.m` 	| 12.0695 | fast rational |
| `DPS_integral.m` 	| 12.0662  | large rational coefficients|
| `DPS.m` 		| 12.0660 | most accurate, with sqrt(3) |
|  |  |  |



**Faster Algorithms**,
Faster variants, using an alternative basis from:
- [G. Beniamini, N. Cheng, O. Holtz, E. Karstadt, O. Schwartz; Sparsifying the Operators of Fast Matrix Multiplication Algorithms](https://arxiv.org/abs/2008.03759).

| Program | Complexity bound | Description |
| :---    |     :---:     |        ---: |
| `Strassen_aternative.m` | $5n^{\log_2(7)}$ | Strassen's algorithm, sparsified |
| `Winograd_aternative.m` | $5n^{\log_2(7)}$ | Winograd's algorithm, sparsified |
| `DPS_aternative.m` | $5n^{\log_2(7)}$ | DPS's algorithm, sparsified |
|  |  |  |

All "aternative" variants use left/right and inverse change of basis (`*CoB*` files) together with an inner sparse multiplication (`*mul*` files).

**Benchmarks**,
Accuracy comparison with symbolic matrix multiplication:
- `accuracy_2x2_real.m`: comparing algorithms accuracy
- `accuracy_alternative_real.m`: comparing alternative change of basis
- `gallery_alternative_real.m`: change of basis on large condition number



--------------------------------------------------------------------------------
**`FMM-plinopt-codegen`, Progam generation via the PLinOpt library**:

Command-line scritps generating optimized matlab programs.

Examples:
`DPS.plo`,
`DPS_evenpow.plo`,
`DPS_smallrat.plo`,
`DPS_intermediate.plo`,
`DPS_integral.plo`,
`DPS_CoB.plo`.
- [PLinOpt library](https://github.com/jgdumas/plinopt): Routines handling linear, bilinear & trilinear programs.



--------------------------------------------------------------------------------
**`FMM-maple-proofs`, Proofs of accuracy bounds in Maple**:

- `Minimal2Norm`: Proof of the Frobenius norm minimal point of Fast Matrix Multiplication in Strassen's orbit
- `LowerBoundGamma`: Proof of a lower bound on the growth factor of Fast Matrix Multiplication in Strassen's orbit
