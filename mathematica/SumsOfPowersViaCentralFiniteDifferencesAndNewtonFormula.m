(* ::Package:: *)

BeginPackage["CentralDifferences`"]

(*BEGIN: Definitions *)
CentralDifferenceInZero::usage=""
CentralDifference::usage=""
T::usage=""
OddPowerIdentity1::usage=""
KnuthOddPowerIdentity::usage=""
TestIdentity::usage=""
CentralFactorial::usage=""
CentralFactorialNumber::usage=""
RiordanPowerIdentity::usage=""
RiordanPowerIdentity1::usage=""
RiordanPowerIdentity2::usage=""
RiordanPowerIdentity3::usage=""
FallingFactorial::usage=""
NewtonsFormulaForCentralDifferences::usage=""
NewtonsFormulaForCentralDifferencesShifted::usage=""
(*END: Definitions *)

(* =========================================================================DOCS END=================================================================== *)

(*BEGIN: Define 0^x = 1 for all x *)
Begin["`Private`"]
Unprotect[Power];
Power[0|0., 0|0.] = 1;
Protect[Power];
(*END: Define 0^x = 1 for all x *)

(* =========================================================================FUNCTIONS BEGIN=========================================================== *)

(*BEGIN: Definitions *)
CentralDifferenceInZero[n_, k_] := Sum[(-1)^(j) * Binomial[k, j] * (k/2 - j)^n, {j, 0, k}];
T[n_, k_] :=  1 / (k!) *  CentralDifferenceInZero[n,k];
CentralDifference[x_, n_, k_] := Sum[(-1)^j * Binomial[k, j] * (x + k/2 -j)^n, {j, 0, k}];
OddPowerIdentity1[n_, m_] := Sum[Binomial[n+k-1, 2k-1] * 1/(2k) * CentralDifference[0, 2m, 2k], {k, 1, m}];
KnuthOddPowerIdentity[n_, m_] := Sum[(2k-1)! * T[2m, 2k] * Binomial[n+k-1, 2k-1], {k, 1, m}];
TestIdentity[n_, m_, t_] := Sum[Binomial[n+k-1, 2k-1] * 1/(2k) * CentralDifference[t, 2m, 2k], {k, 1, m}];

CentralFactorial[x_, k_] := 0 /; k<0;
CentralFactorial[x_, k_] := 1 /; k==0;
CentralFactorial[x_, k_] := x * Product[(x + k/2 -j), {j, 1, k-1}] /; k>0;

CentralFactorialNumber[n_, k_] := 1/k! * CentralDifferenceInZero[n, k];
RiordanPowerIdentity[x_, m_] := Sum[T[m,k] * CentralFactorial[x, k], {k, 1, m}];
RiordanPowerIdentity1[x_, m_] := Sum[1/k! * CentralDifference[0, m, k] * CentralFactorial[x, k], {k, 1, m}];
RiordanPowerIdentity2[x_, m_, t_] := Sum[1/k! * CentralDifference[t, m, k] * CentralFactorial[x, k], {k, 1, m}];
RiordanPowerIdentity3[n_, m_] := Sum[T[2m, 2k] * CentralFactorial[n, 2k], {k, 1, m}];
FallingFactorial[x_, n_] := Product[x-k, {k, 0, n-1}];
NewtonsFormulaForCentralDifferences[n_, m_] := Sum[CentralFactorial[n, j] * 1/j! * CentralDifference[0, m, j], {j, 0, m}];
NewtonsFormulaForCentralDifferencesShifted[x_, h_, n_] := Sum[ 1/k! * CentralDifference[h, n, k] * CentralFactorial[x, k], {k, 0, n}];
(*END: Definitions *)

End[ ]
EndPackage[ ]



