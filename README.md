--------------------------------------------------------------------------------
# Fast-Matrix-Multiplication
### Accuracy of Fast matrix multiplication Algorithms
--------------------------------------------------------------------------------

**Authors**:  Jean-Guillaume Dumas, Clément Pernet, Alexandre Sedoglavic
- [ J-G. Dumas, C. Pernet, A. Sedoglavic; Strassen's algorithm is not optimally accurate, ISSAC 2024, Raleigh, NC USA, pp. 254-263.](https://hal.science/hal-04441653)


**About**:
This is a Fork of the
[Complex-Matrix-Multiplication](https://github.com/zhen06/Complex-Matrix-Multiplication) repository,
created for the numerical experiments of the paper:
- [Dai, Z., Lim, LH; Numerical stability and tensor nuclear norm. Numer. Math. 155, 345--376 (2023)](https://link.springer.com/article/10.1007/s00211-023-01377-5).
- We add more accurate variants of fast matrix multiplication, especially 2x2 recursion

**Requirements**:
- [Matlab](https://fr.mathworks.com/products/matlab.html)
- [PLinOpt](https://github.com/jgdumas/plinopt)


**Installation**:
- See [`auto-docker.run`](https://github.com/jgdumas/Fast-Matrix-Multiplication/blob/main/auto-docker.run)

--------------------------------------------------------------------------------
**`FMM-matlab-benchmarks`, Accuracy benchmarks with Matlab**:
| Program | Growth Factor | Normalized GF | Description |
| :---    |     :---:     |        ---: |        ---: |
| `strassenw.m` 	| 17.8530 | 2.23163| original Winograd's algorithm|
| `strassen.m` 		| 14.8284 | 1.85355| original Strassen's algorithm|
| `DPS_evenpow.m` 	| 12.2034 | 1.52543| using only powers of 2|
| `DPS_smallrat.m` 	| 12.2034 | 1.52543| small rational coefficients|
| `DPS_intermediate.m` 	| 12.0695 | 1.50869| fast rational |
| `DPS_integral.m` 	| 12.0662 | 1.50828| large rational coefficients|
| `DPS.m` 		| 12.0660 | 1.50825| most accurate, with sqrt(3) |
| `FMM_3_3_6.m`		| 395.1 | 7.31667| fast 3x3x6 |
| `FMM_3_6_3.m`		| 395.1 | 7.31667| fast 3x6x3 |
| `FMM_6_3_3.m`		| 395.1 | 7.31667| fast 6x3x3 |
| `FMMa_3_3_6.m`	| 104.1 | 1.92778| fast and accurate 3x3x6 |
| `FMMa_3_6_3.m`	| 104.1 | 1.92778| fast and accurate 3x6x3 |
| `FMMa_6_3_3.m`	| 104.1 | 1.92778| fast and accurate 6x3x3 |
|  |  |  |  |



**Faster Algorithms**,
Faster variants, using an alternative basis from:
- [G. Beniamini, N. Cheng, O. Holtz, E. Karstadt, O. Schwartz; Sparsifying the Operators of Fast Matrix Multiplication Algorithms](https://arxiv.org/abs/2008.03759).

| Program | Complexity bound | Description |
| :---    |     :---:     |        ---: |
| `Strassen_aternative.m` | $5n^{\log_2(7)}$ | Strassen's algorithm, sparsified |
| `Winograd_aternative.m` | $5n^{\log_2(7)}$ | Winograd's algorithm, sparsified |
| `DPS_aternative.m` | $5n^{\log_2(7)}$ | DPS's accurate algorithm, sparsified |
| `FMMa336_alternative.m` | $O(n^{\log_{54}(40^3)})$ | Accurate & sparsified <3;3;6> algorithm |
| `FMMa363_alternative.m` | $O(n^{\log_{54}(40^3)})$ | Accurate & sparsified <3;6;3> algorithm |
| `FMMa633_alternative.m` | $O(n^{\log_{54}(40^3)})$ | Accurate & sparsified <6;3;3> algorithm |
|  |  |  |

All "aternative" variants use left/right and inverse change of basis (`*CoB*` files) together with an inner sparse multiplication (`*mul*` files).

**Benchmarks**,
Accuracy comparison with symbolic matrix multiplication:
- `accuracy_2x2_real.m`: comparing algorithms accuracy
- `accuracy_alternative_real.m`: comparing alternative change of basis
- `gallery_alternative_real.m`: change of basis on large condition number
- `accuracy_3x3x6_real.m`: comparing algorithms accuracy for (3x3).(3x6) matrix ratios
- `accuracy_3x6x3_real.m`: comparing algorithms accuracy for (3x6).(6x3) matrix ratios
- `accuracy_6x3x3_real.m`: comparing algorithms accuracy for (6x3).(3x3) matrix ratios
- `accuracy_54.m`: comparing algorithms accuracy <3;3;6>, <3;6;3> and <6;3;3> matrix ratios



--------------------------------------------------------------------------------
**`FMM-plinopt-codegen`, Progam generation via the PLinOpt library**:

- [fmm.univ-lille.fr](https://fmm.univ-lille.fr/): catalogue of fast matrix multiplication algorithms.
- [PLinOpt library](https://github.com/jgdumas/plinopt): Routines handling linear, bilinear & trilinear programs.


Command-line scripts generating optimized Matlab programs from L,R,P matrices of fast algorithms:
- `sms2matlab.sh`: L.sms R.sms P.sms filename


Examples:
- `sms2matlab.sh data/Lo.sms data/Ro.sms data/Po.sms data/DPS -p 3 1013`
- `sms2matlab.sh data/Ld.sms data/Rd.sms data/Pd.sms data/DPS_evenpow`
- `sms2matlab.sh data/Lr.sms data/Rr.sms data/Pr.sms data/DPS_smallrat`
- `sms2matlab.sh data/Lj.sms data/Rj.sms data/Pj.sms data/DPS_intermediate`
- `sms2matlab.sh data/Li.sms data/Ri.sms data/Pi.sms data/DPS_integral`
- `sms2matlab.sh -a -p 3 1013 data/Lo.sms data/Ro.sms data/Po.sms data/DPS`
- `sms2matlab.sh -n data/La.sms data/Ra.sms data/Pa.sms data/DPS_mul`
- `sms2matlab.sh -p 3 1013 -c data/CLo.sms data/CRo.sms data/CPo.sms data/DPS`
- `336.plo`

Tools:
- `replacer`: regex replacement of variable names in straight-line programs
- `MM.rpl`: Matlab Matrix Multiplication generator from Straight-Line Programs
- `CoB.rpl`: Matlab change of basis generator from Straight-Line Program

--------------------------------------------------------------------------------
**`FMM-maple-proofs`, Proofs of accuracy bounds in Maple**:

- `Minimal2Norm`: Proof of the Frobenius norm minimal point of Fast Matrix Multiplication in Strassen's orbit
- `LowerBoundGamma`: Proof of a lower bound on the growth factor of Fast Matrix Multiplication in Strassen's orbit

--------------------------------------------------------------------------------
