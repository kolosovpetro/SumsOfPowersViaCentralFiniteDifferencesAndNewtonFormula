(* ::Package:: *)

BeginPackage["SumsOfPowersViaCentralFiniteDifferencesAndNewtonFormula`"]

(*BEGIN: Definitions *)
CentralFactorial::usage=""
FallingFactorial::usage=""
CentralDifference::usage=""
NewtonSeriesForPowerInZero::usage=""
NewtonSeriesForPower::usage=""
MultifoldSumOfPowersRecurrence::usage=""
MultifoldSumsOfPowers::usage=""
MultifoldSumsOfPowersOddPowers::usage=""
EvenPower::usage=""
SumsOfEvenPower::usage=""
CentralFactorialNumber::usage=""
NewtonSeriesForPowerInCentralFactorials::usage=""
OddPowersNewtonSeries::usage=""
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
CentralFactorial[n_, k_]:= n * Product[n+ k/2 - j, {j, 1, k-1}];
FallingFactorial[n_, k_] := Product[n-j, {j, 0, k-1}];
CentralDifference[x_, n_, k_] := Sum[(-1)^j * Binomial[k, j] * (x + k/2 -j)^n, {j, 0, k}];
CentralFactorialNumber[n_, k_]:= 1/k! * CentralDifference[0, n, k];
NewtonSeriesForPowerInZero[n_, m_] := Sum[1/(k+1) * Binomial[n+(k+1)/2-1, k] * CentralDifference[0, m+1, k+1], {k, 0, m+1}];
NewtonSeriesForPower[t_, n_, m_] := Sum[1/(k+1) * Binomial[n-t+(k+1)/2-1, k] * CentralDifference[t, m+1, k+1], {k, 0, m+1}];
NewtonSeriesForPowerInCentralFactorials[a_, x_, n_] := Sum[CentralFactorial[x-a, k] * 1/k! * CentralDifference[a, n, k], {k, 0, n}];
MultifoldSumOfPowersRecurrence[r_, n_, m_]:= 0;
MultifoldSumOfPowersRecurrence[r_, n_, m_]:= n^m /; r==0;
MultifoldSumOfPowersRecurrence[r_, n_, m_]:= Sum[MultifoldSumOfPowersRecurrence[r-1, k, m], {k, 1, n}] /; r>0;
MultifoldSumsOfPowers[r_, n_, m_]:= Sum[1/(k+1) * Binomial[n+(k+1)/2 - 1 + r, k+r] * CentralDifference[0, m+1, k+1], {k, 0, m+1}];
MultifoldSumsOfPowersOddPowers[r_, n_, m_]:= Sum[1/(k+1) * Binomial[n+(k+1)/2 - 1 + r, k+r] * CentralDifference[0, 2m+2, k+1], {k, 0, 2m}];
EvenPower[n_, m_] := Sum[1/(2k+1) * Binomial[n+(2k+1)/2-1, 2k] * CentralDifference[0, 2m+1, 2k+1], {k, 0, m}];
SumsOfEvenPower[r_, n_, m_] := Sum[1/(2k+1) * Binomial[n+(2k+1)/2-1+r, 2k+r] * CentralDifference[0, 2m+1, 2k+1], {k, 0, m}];
OddPowersNewtonSeries[t_, n_, m_] := Sum[1/(2k) * Binomial[n+k-t-1, 2k-1] * CentralDifference[t, 2m, 2k], {k, 1, m}];
(*END: Definitions *)

End[ ]
EndPackage[ ]






