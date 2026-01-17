(* ::Package:: *)

BeginPackage["CentralDifferences`"]

(*BEGIN: Definitions *)
CentralDifference::usage=""
CentralFactorial::usage=""
FallingFactorial::usage=""
CentralFactorialNumber2ndKind::usage=""
RiordanPowerIdentity::usage=""
NewtonsFormulaForCentralDifferencesShifted::usage=""
MultifoldSumOfPowersRecurrence::usage=""

ValidateCentralFactorialsInTermsOfFalling::usage=""

ValidateBinomialFormOfCentralFactorials::usage=""

NewtonsFormulaForPowersInZero::usage=""
ValidateNewtonsFormulaForPowersInZero::usage=""

OrdinarySumsOfOddPowersInCentralDifferences::usage=""
ValidateOrdinarySumsOfOddPowersInCentralDifferences::usage=""

MultifoldSumsOfOddPowersInCentralDifferences::usage=""
ValidateMultifoldSumsOfOddPowersInCentralDifferences::usage=""

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
MultifoldSumOfPowersRecurrence[r_, n_, m_]:= 0;
MultifoldSumOfPowersRecurrence[r_, n_, m_]:= n^m /; r==0;
MultifoldSumOfPowersRecurrence[r_, n_, m_]:= Sum[MultifoldSumOfPowersRecurrence[r-1, k, m], {k, 1, n}] /; r>0;

CentralDifference[x_, n_, k_] := Sum[(-1)^j * Binomial[k, j] * (x + k/2 -j)^n, {j, 0, k}];

CentralFactorial[x_, k_] := 0 /; k<0;
CentralFactorial[x_, k_] := 1 /; k==0;
CentralFactorial[x_, k_] := x * Product[(x + k/2 -j), {j, 1, k-1}] /; k>0;

FallingFactorial[x_, n_] := Product[x-k, {k, 0, n-1}];

CentralFactorialNumber2ndKind[n_, k_] := 1/k! * CentralDifference[0, n, k];

RiordanPowerIdentity[x_, m_] := Sum[CentralFactorialNumber2ndKind[m,k] * CentralFactorial[x, k], {k, 1, m}];

NewtonsFormulaForCentralDifferencesShifted[x_, h_, n_] := Sum[ 1/k! * CentralDifference[h, n, k] * CentralFactorial[x, k], {k, 0, n}];

ValidateCentralFactorialsInTermsOfFalling[max_] := Table[CentralFactorial[n,k] - n * FallingFactorial[n + k/2 -1, k-1], {n, -max, max}, {k, 1, max}] //Flatten

ValidateBinomialFormOfCentralFactorials[max_] := Table[CentralFactorial[n,k]/ k! - (n/k) * Binomial[n+k/2-1, k-1], {n, -max, max}, {k, 1, max}] //Flatten

NewtonsFormulaForPowersInZero[n_, m_] := Sum[CentralFactorial[n, j] * 1/j! * CentralDifference[0, m, j], {j, 0, m}];

ValidateNewtonsFormulaForPowersInZero[max_] := Table[n^m - NewtonsFormulaForPowersInZero[n, m], {n, 1, max}, {m, 1, max}] //Flatten

OrdinarySumsOfOddPowersInCentralDifferences[n_, m_] := Sum[1/(2k) * Binomial[n+k, 2k] * CentralDifference[0, 2m, 2k], {k, 1, m}];

ValidateOrdinarySumsOfOddPowersInCentralDifferences[max_] := Table[MultifoldSumOfPowersRecurrence[1, n, 2m-1] - OrdinarySumsOfPowersInCentralDifferences[n, m], {n,1, max}, {m, 1, max}] //Flatten

MultifoldSumsOfOddPowersInCentralDifferences[r_, n_, m_] := Sum[1/(2k) * Binomial[n+k-1+r, 2k-1+r] * CentralDifference[0, 2m, 2k], {k, 1, m}];
ValidateMultifoldSumsOfOddPowersInCentralDifferences[max_] := Table[MultifoldSumOfPowersRecurrence[r, n, 2m-1] - MultifoldSumsOfOddPowersInCentralDifferences[r, n, m], {n,1, max}, {m, 1, max}, {r, 0, max}] //Flatten

(*END: Definitions *)

End[ ]
EndPackage[ ]



