(* ::Package:: *)

(* ::Input:: *)
(*tester[x_]:= x^2*)


(* ::Title:: *)
(*Analytical Solution (just run it!)*)


(* ::Section:: *)
(*MATCHING AT THE Z POLE*)


sin2thetamZ=0.23126;
UmatchTEMP=1/2 ({
 {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2gVu},
 {1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2gVd},
 {0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2gVu},
 {0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2gVd},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 2gVd},
 {0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2gVe},
 {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2gVe},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2gVe},
 {-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2gAu},
 {-1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2gAd},
 {0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2gAu},
 {0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2gAd},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 1, 0, 0, 2gAd},
 {0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2gAe},
 {0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, 0, 0, 2gAe},
 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 2gAe}
});
Umatch=UmatchTEMP/.gVu->1/2-4/3 swsquared/.gAu->-(1/2)/.gVd->-(1/2)+2/3 swsquared/.gAd->1/2/.gVe->-(1/2)+2swsquared/.gAe->1/2/.swsquared->sin2thetamZ;


(* ::Section:: *)
(*ANOMALOUS DIMENSIONS - general form*)


(* ::Subsubsection:: *)
(*Yukawa Interactions*)


OneGenerationBlock=({
 {(Subscript[\[Lambda], u]^2+Subscript[\[Lambda], d]^2)/2, -Subscript[\[Lambda], u]^2/2, -Subscript[\[Lambda], d]^2/2, 0, 0, (Subscript[\[Lambda], u]^2+Subscript[\[Lambda], d]^2)/2},
 {-Subscript[\[Lambda], u]^2, Subscript[\[Lambda], u]^2, 0, 0, 0, Subscript[\[Lambda], u]^2},
 {-Subscript[\[Lambda], d]^2, 0, Subscript[\[Lambda], d]^2, 0, 0, Subscript[\[Lambda], d]^2},
 {0, 0, 0, Subscript[\[Lambda], e]^2/2, -(Subscript[\[Lambda], e]^2/2), Subscript[\[Lambda], e]^2/2},
 {0, 0, 0, -Subscript[\[Lambda], e]^2, Subscript[\[Lambda], e]^2, Subscript[\[Lambda], e]^2},
 {3(Subscript[\[Lambda], u]^2-Subscript[\[Lambda], d]^2), -3 Subscript[\[Lambda], u]^2, 3Subscript[\[Lambda], d]^2, -Subscript[\[Lambda], e]^2, Subscript[\[Lambda], e]^2, 3(Subscript[\[Lambda], u]^2+Subscript[\[Lambda], d]^2)+Subscript[\[Lambda], e]^2}
});
SM1=OneGenerationBlock;
SM2=OneGenerationBlock/.Subscript[\[Lambda], u]->Subscript[\[Lambda], c]/.Subscript[\[Lambda], d]->Subscript[\[Lambda], s]/.Subscript[\[Lambda], e]->Subscript[\[Lambda], \[Mu]];
SM3=OneGenerationBlock/.Subscript[\[Lambda], u]->Subscript[\[Lambda], t]/.Subscript[\[Lambda], d]->Subscript[\[Lambda], b]/.Subscript[\[Lambda], e]->Subscript[\[Lambda], \[Tau]];
ZeroBlock={0,0,0,0,0};
First5Rows=Table[Flatten[Append[Append[Append[Table[SM1[[j,i]],{i,1,5}],ZeroBlock],ZeroBlock],SM1[[j,6]]]],{j,1,5}];
Second5Rows=Table[Flatten[Append[Append[Append[ZeroBlock,Table[SM2[[j,i]],{i,1,5}]],ZeroBlock],SM2[[j,6]]]],{j,1,5}];
Third5Rows=Table[Flatten[Append[Append[Append[ZeroBlock,ZeroBlock],Table[SM3[[j,i]],{i,1,5}]],SM3[[j,6]]]],{j,1,5}];
LastRow=Flatten[Append[Append[Append[Table[SM1[[6,i]],{i,1,5}],Table[SM2[[6,i]],{i,1,5}]],Table[SM3[[6,i]],{i,1,5}]],SM1[[6,6]]+SM2[[6,6]]+SM3[[6,6]]]];
GammaYukawaSYMB=1/(8Pi^2) Append[Flatten[Join[{First5Rows,Second5Rows,Third5Rows}],1],LastRow];


(* ::Subsubsection:: *)
(*Hypercharge Interactions*)


ArrayOfHypercharges={Subscript[y, q],Subscript[y, u],Subscript[y, d],Subscript[y, l],Subscript[y, e],Subscript[y, q],Subscript[y, u],Subscript[y, d],Subscript[y, l],Subscript[y, e],Subscript[y, q],Subscript[y, u],Subscript[y, d],Subscript[y, l],Subscript[y, e],Subscript[y, H]};
ArrayOfColorNumber={3,3,3,1,1,3,3,3,1,1,3,3,3,1,1,1};
ArrayOfWeakNumber={2,1,1,2,1,2,1,1,2,1,2,1,1,2,1,1};
GammaHyperchargeSYMB=4/3 gprime^2/(16Pi^2) Table[ArrayOfColorNumber[[j]]ArrayOfWeakNumber[[j]] ArrayOfHypercharges[[i]]ArrayOfHypercharges[[j]],{i,1,16},{j,1,16}];


(* ::Subsubsection:: *)
(*Fermion Masses Interactions*)


ArrayOfCouplingsToTheZ={gVu,gVd,gVu,gVd,gVd,gVe,gVe,gVe,gAu,gAd,gAu,gAd,gAd,gAe,gAe,gAe};
ArrayOfColorNumberBELOW={3,3,3,3,3,1,1,1,3,3,3,3,3,1,1,1};
ArrayOfFermionMasses=4 GF/Sqrt[2] {0,0,0,0,0,0,0,0,Subscript[m, u]^2,Subscript[m, d]^2,Subscript[m, c]^2,Subscript[m, s]^2,Subscript[m, b]^2,Subscript[m, e]^2,Subscript[m, \[Mu]]^2,Subscript[m, \[Tau]]^2};
GammaFermionMassesSYMB=1/(2Pi^2) Table[ArrayOfColorNumberBELOW[[j]]ArrayOfFermionMasses[[j]] ArrayOfCouplingsToTheZ[[i]]ArrayOfCouplingsToTheZ[[j]],{i,1,16},{j,1,16}];


(* ::Subsubsection:: *)
(*Electromagnetic Interactions*)


ArrayOfElectricCharges={Subscript[Q, u],Subscript[Q, d],Subscript[Q, u],Subscript[Q, d],Subscript[Q, d],Subscript[Q, e],Subscript[Q, e],Subscript[Q, e],Subscript[Q, u],Subscript[Q, d],Subscript[Q, u],Subscript[Q, d],Subscript[Q, d],Subscript[Q, e],Subscript[Q, e],Subscript[Q, e]};
ArrayOfColorNumberBELOW={3,3,3,3,3,1,1,1,3,3,3,3,3,1,1,1};
ArrayOfVectorCouplings={1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0};
GammaElectromagneticSYMB=8/3 el^2/(16Pi^2) Table[ArrayOfVectorCouplings[[i]]ArrayOfVectorCouplings[[j]]ArrayOfColorNumberBELOW[[j]] ArrayOfElectricCharges[[i]]ArrayOfElectricCharges[[j]],{i,1,16},{j,1,16}];


(* ::Section:: *)
(*ANOMALOUS DIMENSIONS - numerical*)


(* ::Subsubsection:: *)
(*Input parameter*)


HiggsVEV=246.21971;
g1ATmtop=Sqrt[5/3]0.35830;
g2ATmtop=0.64779;


(* ::Subsubsection:: *)
(*Anomalous Dimensions Matrices*)


GammaYukawaSYMB2=GammaYukawaSYMB/.Subscript[\[Lambda], u]->0/.Subscript[\[Lambda], d]->0/.Subscript[\[Lambda], e]->0/.Subscript[\[Lambda], c]->0/.Subscript[\[Lambda], s]->0/.Subscript[\[Lambda], \[Mu]]->0/.Subscript[\[Lambda], \[Tau]]->lambdatau/.Subscript[\[Lambda], b]->lambdabottom/.Subscript[\[Lambda], t]->lambdatop;
GammaHyperchargeSYMB2=GammaHyperchargeSYMB/.Subscript[y, q]->1/6/.Subscript[y, u]->2/3/.Subscript[y, d]->-(1/3)/.Subscript[y, l]->-(1/2)/.Subscript[y, e]->-1/.Subscript[y, H]->1/2;
GammaSM=GammaYukawaSYMB2+GammaHyperchargeSYMB2;
GammaSMfinal=GammaSM/.gprime->0.357003/.lambdatop->0.973644/.lambdabottom->0.0164834/.lambdatau->0.00994705;
GammaFermionMassesSYMB2=GammaFermionMassesSYMB/.Subscript[m, u]->0/.Subscript[m, d]->0/.Subscript[m, e]->0/.Subscript[m, c]->0/.Subscript[m, s]->0/.Subscript[m, \[Mu]]->0/.gVu->1/2-4/3 swsquared/.gAu->-(1/2)/.gVd->-(1/2)+2/3 swsquared/.gAd->1/2/.gVe->-(1/2)+2swsquared/.gAe->1/2/.Subscript[m, \[Tau]]->1.777/.GF->1/(Sqrt[2]HiggsVEV^2)/.swsquared->sin2thetamZ;
GammaElectromagneticSYMB2=GammaElectromagneticSYMB/.Subscript[Q, u]->2/3/.Subscript[Q, d]->-(1/3)/.Subscript[Q, e]->-1;
GammaEMSM=GammaFermionMassesSYMB2+GammaElectromagneticSYMB2;
GammaEMSMfinal=GammaEMSM/.el->0.313521/.Subscript[m, b]->2.86/.Subscript[m, \[Tau]]->1.777;


(* ::Section:: *)
(*Evolution Matrix*)


tnuclear=Log[1/91.1875];
UEMSM=MatrixExp[GammaEMSMfinal tnuclear];
EvolutionSM[t_]:=MatrixExp[-GammaSMfinal t];
Ufull[t_]:=UEMSM.Umatch.EvolutionSM[t];


(* ::Section:: *)
(*Function for the benchmarks*)


ArrayNuclearVectorAnalytical[\[CapitalLambda]_]:=ArrayNuclearAnalytical[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,\[CapitalLambda]];
ArrayNuclearAxialAnalytical[\[CapitalLambda]_]:=ArrayNuclearAnalytical[-1,1,1,-1,1,-1,1,1,-1,1,-1,1,1,-1,1,0,\[CapitalLambda]];
ArrayNuclearVectorQuarksAnalytical[\[CapitalLambda]_]:=ArrayNuclearAnalytical[1,1,1,0,0,1,1,1,0,0,1,1,1,0,0,0,\[CapitalLambda]];
ArrayNuclearAxialQuarksAnalytical[\[CapitalLambda]_]:=ArrayNuclearAnalytical[-1,1,1,0,0,-1,1,1,0,0,-1,1,1,0,0,0,\[CapitalLambda]];
ArrayNuclearVectorLeptonsAnalytical[\[CapitalLambda]_]:=ArrayNuclearAnalytical[0,0,0,1,1,0,0,0,1,1,0,0,0,1,1,0,\[CapitalLambda]];
ArrayNuclearAxialLeptonsAnalytical[\[CapitalLambda]_]:=ArrayNuclearAnalytical[0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1,0,\[CapitalLambda]];
ArrayNuclearVectorThirdAnalytical[\[CapitalLambda]_]:=ArrayNuclearAnalytical[0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,0,\[CapitalLambda]];
ArrayNuclearAxialThirdAnalytical[\[CapitalLambda]_]:=ArrayNuclearAnalytical[0,0,0,0,0,0,0,0,0,0,-1,1,1,-1,1,0,\[CapitalLambda]];


(* ::Title:: *)
(*Other stuff...*)


(* ::Text:: *)
(*Functions for calculating running of the coefficients*)


(* ::Text:: *)
(*Calculating the yakawa-driven running in the SM\[Chi]*)


(* ::Input:: *)
(*Calc\[Gamma]\[Lambda][]:=*)
(*( *)
(*Clear[\[Lambda]u, \[Lambda]d, \[Lambda]c, \[Lambda]s, \[Lambda]t, \[Lambda]b, \[Lambda]e, \[Lambda]\[Mu], \[Lambda]\[Tau]];*)
(**)
(*\[Gamma]\[Lambda] = ConstantArray[0,{16,16}];*)
(**)
(*\[Gamma]\[Lambda][[1,1]] = (\[Lambda]u^2 + \[Lambda]d^2)/2;*)
(*\[Gamma]\[Lambda][[1,2]] = -\[Lambda]u^2/2;*)
(*\[Gamma]\[Lambda][[1,3]] = -\[Lambda]d^2/2;*)
(*\[Gamma]\[Lambda][[1,16]] = (\[Lambda]u^2 + \[Lambda]d^2)/2;*)
(**)
(*\[Gamma]\[Lambda][[2,1]] = -\[Lambda]u^2;*)
(*\[Gamma]\[Lambda][[2,2]] = \[Lambda]u^2;*)
(*\[Gamma]\[Lambda][[2,16]] = \[Lambda]u^2;*)
(**)
(*\[Gamma]\[Lambda][[3,1]] = -\[Lambda]d^2;*)
(*\[Gamma]\[Lambda][[3,3]] = \[Lambda]d^2;*)
(*\[Gamma]\[Lambda][[3,16]] = \[Lambda]d^2;*)
(**)
(*\[Gamma]\[Lambda][[4,4]] = \[Lambda]e^2/2;*)
(*\[Gamma]\[Lambda][[4,5]] = -\[Lambda]e^2/2;*)
(*\[Gamma]\[Lambda][[4,16]] = \[Lambda]e^2/2;*)
(**)
(*\[Gamma]\[Lambda][[5,4]] = - \[Lambda]e^2;*)
(*\[Gamma]\[Lambda][[5,5]] = \[Lambda]e^2;*)
(*\[Gamma]\[Lambda][[5,16]] = \[Lambda]e^2;*)
(**)
(*\[Gamma]\[Lambda][[6,6]] = (\[Lambda]c^2 + \[Lambda]s^2)/2;*)
(*\[Gamma]\[Lambda][[6, 7]] = -\[Lambda]c^2/2;*)
(*\[Gamma]\[Lambda][[6,8]] = -\[Lambda]s^2/2;*)
(*\[Gamma]\[Lambda][[6, 16]] = (\[Lambda]c^2 + \[Lambda]s^2)/2;*)
(**)
(*\[Gamma]\[Lambda][[7, 6]] = -\[Lambda]c^2;*)
(*\[Gamma]\[Lambda][[7,7]] = \[Lambda]c^2;*)
(*\[Gamma]\[Lambda][[7,8]] = \[Lambda]c^2;*)
(**)
(*\[Gamma]\[Lambda][[8, 6]] = -\[Lambda]s^2;*)
(*\[Gamma]\[Lambda][[8, 8]] = \[Lambda]s^2;*)
(*\[Gamma]\[Lambda][[8,16]] = \[Lambda]s^2;*)
(**)
(*\[Gamma]\[Lambda][[9,9]] = \[Lambda]\[Mu]^2/2;*)
(*\[Gamma]\[Lambda][[9,10]] = -\[Lambda]\[Mu]^2/2;*)
(*\[Gamma]\[Lambda][[9,16]] = \[Lambda]\[Mu]^2/2;*)
(**)
(*\[Gamma]\[Lambda][[10,9]] = -\[Lambda]\[Mu]^2;*)
(*\[Gamma]\[Lambda][[10,10]] = \[Lambda]\[Mu]^2;*)
(*\[Gamma]\[Lambda][[10,16]] = \[Lambda]\[Mu]^2;*)
(**)
(*\[Gamma]\[Lambda][[11,11]] = (\[Lambda]t^2 + \[Lambda]b^2)/2;*)
(*\[Gamma]\[Lambda][[11,12]] = -\[Lambda]t^2/2;*)
(*\[Gamma]\[Lambda][[11,13]] = -\[Lambda]b^2/2;*)
(*\[Gamma]\[Lambda][[11,16]] = (\[Lambda]t^2 + \[Lambda]b^2)/2;*)
(**)
(*\[Gamma]\[Lambda][[12,11]] = -\[Lambda]t^2;*)
(*\[Gamma]\[Lambda][[12,12]] = \[Lambda]t^2;*)
(*\[Gamma]\[Lambda][[12,16]] = \[Lambda]t^2;*)
(**)
(*\[Gamma]\[Lambda][[13,11]] = -\[Lambda]b^2;*)
(*\[Gamma]\[Lambda][[13,13]] = \[Lambda]b^2;*)
(*\[Gamma]\[Lambda][[13,16]] = \[Lambda]b^2;*)
(**)
(*\[Gamma]\[Lambda][[14,14]] = \[Lambda]\[Tau]^2/2;*)
(*\[Gamma]\[Lambda][[14,15]] = -\[Lambda]\[Tau]^2/2;*)
(*\[Gamma]\[Lambda][[14,16]] = \[Lambda]\[Tau]^2/2;*)
(**)
(*\[Gamma]\[Lambda][[15,14]] = -\[Lambda]\[Tau]^2;*)
(*\[Gamma]\[Lambda][[15,15]] = \[Lambda]\[Tau]^2;*)
(*\[Gamma]\[Lambda][[15,16]] = \[Lambda]\[Tau]^2;*)
(**)
(*\[Gamma]\[Lambda][[16,1]] = 3(\[Lambda]u^2 - \[Lambda]d^2);*)
(*\[Gamma]\[Lambda][[16,2]] = -3\[Lambda]u^2;*)
(*\[Gamma]\[Lambda][[16,3]] = 3\[Lambda]d^2;*)
(*\[Gamma]\[Lambda][[16,4]] = -\[Lambda]e^2;*)
(*\[Gamma]\[Lambda][[16,5]] = \[Lambda]e^2;*)
(**)
(*\[Gamma]\[Lambda][[16,6]] = 3(\[Lambda]c^2 - \[Lambda]s^2);*)
(*\[Gamma]\[Lambda][[16,7]] = -3\[Lambda]c^2;*)
(*\[Gamma]\[Lambda][[16,8]] = 3\[Lambda]s^2;*)
(*\[Gamma]\[Lambda][[16,9]] = -\[Lambda]\[Mu]^2;*)
(*\[Gamma]\[Lambda][[16,10]] = \[Lambda]\[Mu]^2;*)
(**)
(*\[Gamma]\[Lambda][[16,11]] = 3(\[Lambda]t^2 - \[Lambda]b^2);*)
(*\[Gamma]\[Lambda][[16,12]] = -3\[Lambda]t^2;*)
(*\[Gamma]\[Lambda][[16,13]] = 3\[Lambda]b^2;*)
(*\[Gamma]\[Lambda][[16,14]] = -\[Lambda]\[Tau]^2;*)
(*\[Gamma]\[Lambda][[16,15]] = \[Lambda]\[Tau]^2;*)
(**)
(*\[Gamma]\[Lambda][[16,16]] = 3*(\[Lambda]u^2 + \[Lambda]d^2 + \[Lambda]c^2 + \[Lambda]s^2 + \[Lambda]b^2 + \[Lambda]t^2);*)
(*\[Gamma]\[Lambda][[16,16]] += \[Lambda]e^2 + \[Lambda]\[Mu]^2 + \[Lambda]\[Tau]^2;*)
(**)
(*\[Gamma]\[Lambda] *= 1/(8*Pi^2);*)
(*Return[\[Gamma]\[Lambda]];*)
(*)*)
(*Calc\[Gamma]\[Lambda]::Usage = *)
(*"Returns the yukawa contribution to the anomalous dimension matrix in the unbroken SM\[Chi].";*)


(* ::Text:: *)
(*Calculating the hypercharge-driven running in the SM\[Chi]*)


(* ::Input:: *)
(*Calc\[Gamma]Y[]:=*)
(*( *)
(*Clear[yq, yu, yd, yl, ye,yH];*)
(**)
(*\[Gamma]Y = ConstantArray[1,{16,16}];*)
(**)
(*\[Gamma]Y[[{1, 6, 11}, All]] *= yq;*)
(*\[Gamma]Y[[{2, 7, 12}, All]] *= yu;*)
(*\[Gamma]Y[[{3, 8, 13}, All]] *= yd;*)
(*\[Gamma]Y[[{4, 9, 14}, All]] *= yl;*)
(*\[Gamma]Y[[{5, 10, 15}, All]] *= ye;*)
(*\[Gamma]Y[[16, All]] *= yH;*)
(**)
(*\[Gamma]Y[[All,{1, 6, 11}]] *= 6*yq;*)
(*\[Gamma]Y[[All,{2, 7, 12}]] *= 3*yu;*)
(*\[Gamma]Y[[All,{3, 8, 13}]] *= 3*yd;*)
(*\[Gamma]Y[[All,{4, 9, 14}]] *= 2*yl;*)
(*\[Gamma]Y[[All,{5, 10, 15}]] *= ye;*)
(*\[Gamma]Y[[ All,16]] *= yH;*)
(**)
(*\[Gamma]Y *= 4*gp^2/(3*16*Pi^2);*)
(**)
(*Return[\[Gamma]Y];*)
(*)*)
(*Calc\[Gamma]Y::Usage = *)
(*"Returns the hypercharge contribution to the anomalous dimension matrix in the unbroken SM\[Chi].";*)
(**)


(* ::Text:: *)
(*Calculating the full anomalous dimension matrix in the SM\[Chi]*)


(* ::Input:: *)
(*Calc\[Gamma]SM\[Chi][]:=*)
(*( *)
(*Calculate the different contributions separately;*)
(**)
(*\[Gamma]SM\[Chi]= Calc\[Gamma]Y[] + Calc\[Gamma]\[Lambda][];*)
(**)
(*Return[\[Gamma]SM\[Chi]];*)
(**)
(*)*)
(*Calc\[Gamma]SM\[Chi]::Usage = *)
(*"Returns the full anomalous dimension matrix for the unbroken SM\[Chi] (including numerical values.";*)
(**)


(* ::Text:: *)
(*Calculating the EM-driven running in the EMSM\[Chi]*)


(* ::Input:: *)
(*Calc\[Gamma]EM[]:=*)
(*( *)
(*Clear[Qu,Qd,Qe];*)
(*(*Qu = 2/3;*)
(*Qd = 2/3;*)
(*Qe = 1;*)*)
(*\[Gamma]EM = ConstantArray[0,{16,16}];*)
(*\[Gamma]EM[[1;;8,1;;8]] += 1;*)
(*\[Gamma]EM[[{1,3},1;;8]] *= Qu;*)
(*\[Gamma]EM[[1;;8,{1,3}]] *= Qu;*)
(**)
(*\[Gamma]EM[[{2,4,5},1;;8]] *= Qd;*)
(*\[Gamma]EM[[1;;8,{2,4,5}]] *= Qd;*)
(**)
(*\[Gamma]EM[[1;;8,1;;5]] *= 3;*)
(**)
(*\[Gamma]EM[[6;;8,1;;8]] *= Qe;*)
(*\[Gamma]EM[[1;;8, 6;;8]] *= Qe;*)
(**)
(*\[Gamma]EM *= 8 *e^2/(3 *16 *Pi^2);*)
(**)
(*Return [\[Gamma]EM];*)
(*)*)
(*Calc\[Gamma]EM::Usage = *)
(*"Returns the EM interaction contribution to the  anomalous dimension matrix in the EMSM\[Chi].";*)
(**)


(* ::Text:: *)
(*Calculating the mass-driven running in the EMSM\[Chi]*)


(* ::Input:: *)
(*Calc\[Gamma]m[]:=*)
(*( *)
(*Clear[mu, md, mc, ms, mb, me , m\[Mu], m\[Tau], gAu, gAd, gAe, gVu, gVd, gVe];*)
(**)
(*\[Gamma]m = ConstantArray[0,{16,16}];*)
(*\[Gamma]m[[All, 9;;16]] += 1;*)
(**)
(*\[Gamma]m[[All,9]] *= 3 mu^2 gAu;*)
(*\[Gamma]m[[All,10]] *= 3 md^2 gAd;*)
(*\[Gamma]m[[All,11]] *= 3 mc^2 gAu;*)
(*\[Gamma]m[[All,12]] *= 3 ms^2 gAd;*)
(*\[Gamma]m[[All,13]] *= 3 mb^2 gAd;*)
(**)
(*\[Gamma]m[[All,14]] *= me^2 gAe;*)
(*\[Gamma]m[[All,15]] *= m\[Mu]^2 gAe;*)
(*\[Gamma]m[[All,16]] *= m\[Tau]^2 gAe;*)
(**)
(*\[Gamma]m[[{1, 3}, All]] *= gVu;*)
(*\[Gamma]m[[{2,4,5}, All]] *= gVd;*)
(*\[Gamma]m[[{6,7,8}, All]] *= gVe;*)
(**)
(*\[Gamma]m[[{9, 11}, All]] *= gAu;*)
(*\[Gamma]m[[{10,12,13}, All]] *= gAd;*)
(*\[Gamma]m[[{14,15,16}, All]] *= gAe;*)
(**)
(*\[Gamma]m *= Sqrt[2]*GF/(Pi^2);*)
(**)
(*Return [\[Gamma]m];*)
(*)*)
(*Calc\[Gamma]m::Usage = *)
(*"Returns the fermion mass contribution to the anomalous dimension matrix in the EMSM\[Chi].";*)
(**)


(* ::Text:: *)
(*Calculating the full anomalous dimension matrix in the EMSM\[Chi]*)


(* ::Input:: *)
(*Calc\[Gamma]EMSM\[Chi][]:=*)
(*( *)
(*Calculate the different contributions separately;*)
(**)
(*\[Gamma] = Calc\[Gamma]EM[] + Calc\[Gamma]m[];*)
(**)
(*Return[\[Gamma]];*)
(**)
(*)*)
(*Calc\[Gamma]EMSM\[Chi]::Usage = *)
(*"Returns the full anomalous dimension matrix for the EMSM\[Chi] (including numerical values.";*)
(**)


(* ::Text:: *)
(*Calculate the EWSB matching matrix*)


(* ::Input:: *)
(*CalcUmatch[]:=*)
(*( *)
(*Clear[ gAu, gAd, gAe, gVu, gVd, gVe];*)
(**)
(*Umatch = ConstantArray[0,{16,16}];*)
(**)
(*Umatch[[1,{1,2}]] = 1;*)
(*Umatch[[2,{1,3}]] = 1;*)
(*Umatch[[3, {6, 7}]] = 1;*)
(*Umatch[[4, {6, 8}]] = 1;*)
(*Umatch[[5, {11,13}]] = 1;*)
(*Umatch[[6, {4, 5}]] = 1;*)
(*Umatch[[7, {9, 10}]] = 1;*)
(*Umatch[[8, {14, 15}]] = 1;*)
(**)
(*Umatch[[9, 1]] = -1;*)
(*Umatch[[9, 2]] = 1;*)
(*Umatch[[10,1]] = -1;*)
(*Umatch[[10,3]] = 1;*)
(*Umatch[[11, 6]] = -1;*)
(*Umatch[[11, 7]] = 1;*)
(*Umatch[[12, 6]] = -1;*)
(*Umatch[[12,8]] = 1;*)
(*Umatch[[13, 11]] = -1;*)
(*Umatch[[13, 13]] = 1;*)
(*Umatch[[14, 4]] = -1;*)
(*Umatch[[14, 5]] = 1;*)
(*Umatch[[15, 9]] = -1;*)
(*Umatch[[15, 10]] = 1;*)
(*Umatch[[16, 14]] = -1;*)
(*Umatch[[16, 15]] = 1;*)
(**)
(*Umatch[[{1, 3}, 16]] = 2*gVu;*)
(*Umatch[[{2, 4, 5}, 16]] = 2*gVd;*)
(*Umatch[[{6, 7, 8}, 16]] = 2*gVe;*)
(**)
(*Umatch[[{9, 11}, 16]] = 2*gAu;*)
(*Umatch[[{10, 12, 13}, 16]] = 2*gAd;*)
(*Umatch[[{14, 15, 16}, 16]] = 2*gAe;*)
(**)
(*Umatch *= 1/2;*)
(**)
(*Return[Umatch];*)
(*)*)
(*CalcUmatch::Usage = *)
(*"Returns the matching matrix to be applied at the EWSB scale (\[Mu] ~ mZ) to go between the SM\[Chi] and EMSM\[Chi].";*)


(* ::Text:: *)
(*Set values of Standard Model parameters*)


(* ::Input:: *)
(*SetSMparameters[]:=*)
(*( *)
(*Set the numerical values of masses and couplings - PDG July 2014;*)
(**)
(*mu=2.3*^-3; (*in GeV*)*)
(*md=4.8*^-3;*)
(*ms = 95*^-3;*)
(*mc = 1.275;*)
(*mb = 4.66;*)
(*mt = 173.34;*)
(*me = 0.511*^-3;*)
(*m\[Mu] = 105.7*^-3;*)
(*m\[Tau] = 1776.82*^-3;*)
(**)
(*Qu = 2/3;*)
(*Qd = -1/3;*)
(*Qe = -1;*)
(**)
(*GF = 1.16637*^-5 ;(*in GeV^-2*)*)
(*\[Alpha]EM = 7.2974*^-3;*)
(*e = Sqrt[4*Pi*\[Alpha]EM];*)
(**)
(*sw2= 0.2312; *)
(* gAu =-1/2;*)
(*gAd = 1/2;*)
(*gAe = 1/2;*)
(*gVu = 1/2-sw2*4.0/3;*)
(*gVd = -1/2 + sw2*2.0/3;*)
(*gVe = -1/2 + sw2*2.0;*)
(**)
(*(*gVe = -0.03783;*)
(*gVu = 0.25;*)
(*gVd = -0.33;*)
(*gAe = -0.501;*)
(*gAu = 0.5;*)
(*gAd = -0.523;*)*)
(**)
(*yq = 1/6;*)
(*yu  = 2/3;*)
(*yd = -1/3;*)
(*yl = -1/2;*)
(*ye = -1;*)
(*yH = 1/2;*)
(**)
(*vH = 246;(*in GeV*)*)
(*\[Lambda]u=mu*Sqrt[2]/vH; *)
(*\[Lambda]d=md*Sqrt[2]/vH;*)
(*\[Lambda]s = ms*Sqrt[2]/vH;*)
(*\[Lambda]c = mc*Sqrt[2]/vH;*)
(*\[Lambda]b = mb* Sqrt[2]/vH;*)
(*\[Lambda]t = mt*Sqrt[2]/vH;*)
(*\[Lambda]e = me*Sqrt[2]/vH;*)
(*\[Lambda]\[Mu] = m\[Mu]*Sqrt[2]/vH;*)
(*\[Lambda]\[Tau] = m\[Tau]*Sqrt[2]/vH;*)
(**)
(*gp = 0.357;*)
(*Return[0];*)
(*)*)
(*SetSMparameters::Usage = *)
(*"Sets the numerical values of the SM parameters. Returns 0.";*)


(* ::Text:: *)
(*Initialise the running matrices and parameters*)


(* ::Input:: *)
(*InitRunning[]:=*)
(*( *)
(*mZ = 91.2;*)
(*\[Gamma]SM\[Chi] = Calc\[Gamma]SM\[Chi][];*)
(*\[Gamma]EMSM\[Chi] = Calc\[Gamma]EMSM\[Chi][];*)
(*Umatch = CalcUmatch[];*)
(*(*invUmatch = Inverse[Umatch];*)*)
(*SetSMparameters[];*)
(*)*)
(*InitRunning::Usage = *)
(*"Initialise anomalous dimension matrix and SM parameter values for use in calculating running.";*)
(**)


(* ::Text:: *)
(*Calculate the running of the couplings between two energies*)


(* ::Input:: *)
(*CalcRunning[Ci_ ,E1_, E2_]:=*)
(*( *)
(*Check whether E1, E2 are above or below EWSB-scale;*)
(*mZ = 91.2;*)
(**)
(*(*At the minute we're assuming that E1 > E2*)*)
(*Cf = Which[ *)
(*(E1 > mZ && E2 > mZ), MatrixExp[-\[Gamma]SM\[Chi]*Log[E1/E2]].Ci,*)
(*(E1 > mZ && E2  < mZ),MatrixExp[-\[Gamma]EMSM\[Chi]*Log[mZ/E2]].Umatch.MatrixExp[-\[Gamma]SM\[Chi]*Log[E1/mZ]].Ci,*)
(*(E1 < mZ && E2 < mZ), MatrixExp[-\[Gamma]EMSM\[Chi]*Log[E1/E2]].Ci];*)
(**)
(*Return[Cf];*)
(*)*)
(*CalcRunning::Usage = *)
(*"Run the coupling vector Ci from scale E1 to E2 (with E1 > E2). E1 and E2 may be either above or below the EWSB scale.";*)


(* ::Input:: *)
(**)


(* ::Text:: *)
(*Functions for manipulating coupling vectors*)


(* ::Text:: *)
(*Initialise vectors of couplings*)


initCouplings[]:=
( 
Return[ConstantArray[0,16]];
)
InitCouplings::Usage = 
"Returns a vector with 16 entries, all zero (suitable for storing EFT couplings).";



(* ::Text:: *)
(*Setting and getting values in coupling vectors*)


(* ::Input:: *)
(*OpInd[op_]:=*)
(*( *)
(*Nop = Length[op];*)
(*op2 = op;*)
(*If [Nop==0, *)
(*(Nop= 1;*)
(*op2 = ConstantArray[0,1];*)
(*op2[[1]] = op;)*)
(*];*)
(*indout = ConstantArray[0,Nop];*)
(*For[i=1,i<Nop+1,i++,*)
(*indout[[i]]  = Switch[op2[[i]], *)
(*"uL", 1,*)
(*"dL", 1,*)
(*"uR", 2,*)
(*"dR", 3,*)
(*"eL", 4,*)
(*"eR", 5,*)
(*"cL", 6,*)
(*"sL", 6,*)
(*"cR", 7,*)
(*"sR", 8,*)
(*"\[Mu]L", 9,*)
(*"\[Mu]R", 10,*)
(*"tL", 11,*)
(*"bL", 11,*)
(*"tR", 12,*)
(*"bR", 13,*)
(*"\[Tau]L", 14,*)
(*"\[Tau]R", 15,*)
(*"H", 16,*)
(*"uV", 1, *)
(*"dV", 2, *)
(*"cV", 3,*)
(*"sV", 4,*)
(*"bV", 5,*)
(*"eV", 6,*)
(*"\[Mu]V", 7,*)
(*"\[Tau]V", 8,*)
(*"uA", 9, *)
(*"dA", 10, *)
(*"cA", 11,*)
(*"sA", 12,*)
(*"bA", 13,*)
(*"eA", 14,*)
(*"\[Mu]A", 15,*)
(*"\[Tau]A", 16*)
(*];*)
(*];*)
(*If[Nop == 1,*)
(*Return[indout[[1]]];,*)
(*Return[indout];]*)
(*)*)
(*OpInd::Usage = *)
(*"Returns the coupling index corresponding to the operator or list of operators 'op'.";*)


(* ::Text:: *)
(*Generate coupling vectors with preset operator structures*)


setBenchmark[benchmark_]:=
(
Switch[benchmark,
"Higgs",
Return[{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0}],
"UniversalVector", 
Return[{1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.0}],
"UniversalAxial",
Return[{-1.0,1.0,1.0,-1.0,1.0,-1.0,1.0,1.0,-1.0,1.0,-1.0,1.0,1.0,-1.0,1.0,0.0}],
"QuarksVector",
Return[{1.0,1.0,1.0,0.0,0.0,1.0,1.0,1.0,0.0,0.0,1.0,1.0,1.0,0.0,0.0,0.0}],
"QuarksAxial",
Return[{-1.0,1.0,1.0,0.0,0.0,-1.0,1.0,1.0,0.0,0.0,-1.0,1.0,1.0,0.0,0.0,0.0}],
"LeptonsVector",
Return[{0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,1.0,1.0,0.0,0.0,0.0,1.0,1.0,0.0}],
"LeptonsAxial",
Return[{0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,-1.0,1.0,0.0,0.0,0.0,-1.0,1.0,0.0}],
"ThirdVector",
Return[{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0,0.0}],
"ThirdAxial",
Return[{0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,1.0,1.0,-1.0,1.0,0.0}]];

(*If no match was made, return zeros...*)
Print["Error in runDM.setBenchmark: Benchmark <<" <> benchmark <>">> not found..."];
Print["Options are: 'Higgs', 'UniversalVector', 'UniversalAxial', 'QuarksVector', 'QuarksAxial', 'LeptonsVector', 'LeptonsAxial', 'ThirdVector', 'ThirdAxial'..."];
Print["Returning empty coupling vector..."];
couplings = ConstantArray[0,16];
Return[couplings];
)



(* ::Text:: *)
(*Displaying the contents of a coupling vector*)


(* ::Input:: *)
(*DisplayCouplings[C_, E_, bilinear_]:=*)
(*( *)
(*mZ = 91.2;*)
(*DMcurr = "    \[Chi]" <> bilinear <> "\[Chi]";*)
(*If[E< mZ,*)
(*( *)
(*Print["EMSM\[Chi] operator couplings:"]; *)
(*Print[DMcurr <> "uVu: ", C[[1]]];*)
(*Print[DMcurr <> "dVd: ", C[[2]]];*)
(*Print[DMcurr <> "cVc: ", C[[3]]];*)
(*Print[DMcurr <> "sVs: ", C[[4]]];*)
(*Print[DMcurr <> "bVb: ", C[[5]]];*)
(*Print[DMcurr <> "eVe: ", C[[6]]];*)
(*Print[DMcurr <> "\[Mu]V\[Mu]: ", C[[7]]];*)
(*Print[DMcurr <> "\[Tau]V\[Tau]: ", C[[8]]];*)
(*Print[DMcurr <> "uAu: ", C[[9]]];*)
(*Print[DMcurr <> "dAd: ", C[[10]]];*)
(*Print[DMcurr <> "cAc: ", C[[11]]];*)
(*Print[DMcurr <> "sAs: ", C[[12]]];*)
(*Print[DMcurr <> "bAb: ", C[[13]]];*)
(*Print[DMcurr <> "eAe: ", C[[14]]];*)
(*Print[DMcurr <> "\[Mu]A\[Mu]: ", C[[15]]];*)
(*Print[DMcurr <> "\[Tau]A\[Tau]: ", C[[16]]];*)
(**)
(*)]*)
(*)*)
(**)


evolutionMat[E1_, E2_]:=
( 
mZ = 91.1875;

(*At the minute we're assuming that E1 > E2*)
Emat = Which[ 
(E1 > mZ && E2 > mZ), MatrixExp[-GammaSMfinal*Log[E1/E2]],
(E1 > mZ && E2  < mZ),MatrixExp[-GammaEMSMfinal*Log[mZ/E2]].Umatch.MatrixExp[-GammaSMfinal*Log[E1/mZ]],
(E1 < mZ && E2 < mZ), MatrixExp[-GammaEMSMfinal*Log[E1/E2]].Umatch];

Return[Emat];
)
evolutionMat::Usage = 
"Run the coupling vector Ci from scale E1 to E2 (with E1 > E2). E1 and E2 may be either above or below the EWSB scale.";


evolveCouplings[c_, E1_, E2_]:=
(
	Return[evolutionMat[E1, E2].c];
)
evolveCouplings::Usage = 
"   "


lightqCouplings[c_, E1_, E2_]:=
(
	cf = evolveCouplings[c, E1, E2];
	Return[cf[[{1,2,4,9,10,12}]]];
)
lightqCouplings::Usage = 
"   "
