(* ::Package:: *)

ComputeDI[samplepoint_]:=Module[{SubReduList,ReductionRules,ExtraRelations,tpoint,diffop,kinematicsT,loopMomentumScalarProductsT,mandelstamScalarProductsT,dhxbAll,ReducedDiffrential,NULLS,dhxbMatrices,dhxbMatricesSaveFormat,reductionRules,PureIntegralN,diffOp,eps,tPoint,repl,PureIntegralN2,Basis,NullInts,BasisChangMat,BasisChangMatInv,diffInts,diffIntsIBREDUCED,NonReduced,NullInts2,epsrepl,ZeroInts,ZeroIntsRepl,PureIntegralN3,PureIntegralN4,SubLaporta,NonZeros,DeleteList},
SubReduList=Get[ToString[StringForm["``/KIRA/splitres/phasespace_``.m",Directory[],IntegerString[samplepoint,10,5]]]];
tPoint=Thread[Rule[invariants,SubReduList[[2,1]]]];

reductionRules=ReductionRules=GetReductionTableFromList[SubReduList[[2,2]]];

PureIntegralN=Expand[PureIntegrals/.tPoint/.root[n_]:>PowerMod[n,1/2,$P],Modulus->$P];
diffOp=getDifferentialOperator[tPoint];
epsrepl={(ep->(4-d)/2)/.tPoint};

ExtraRelations={};
(*ExtraRelations=GetExtraRelationReplacements[reductionRules,tPoint];*)
PureIntegralN2=Expand[PureIntegralN/.reductionRules,Modulus->$P];
PureIntegralN3=Expand[PureIntegralN2/.ExtraRelations,Modulus->$P];

NonZeros=Cases[reductionRules,_hxb,Infinity]//DeleteDuplicates;
ZeroInts=Complement[Ureduce,NonZeros];
ZeroIntsRepl=Thread[Rule[ZeroInts,0]];
PureIntegralN4=Expand[PureIntegralN3/.ExtraRelations/.ZeroIntsRepl,Modulus->$P];
SubLaporta=Cases[PureIntegralN4,_hxb,Infinity]//DeleteDuplicates;
PureIntegralN4=Expand[PureIntegralN4/.epsrepl,Modulus->$P];

BasisChangMat=Table[CoefficientArrays[PureIntegralN4[[pos]],SubLaporta],{pos,1,Length[PureIntegralN]}];
BasisChangMat=Normal[BasisChangMat[[All,2]]];
BasisChangMatInv=Inverse[BasisChangMat,Modulus->$P];

loopMomentumScalarProductsT=Expand[loopMomentumScalarProducts/.tPoint,Modulus->$P];
mandelstamScalarProductsT=Expand[mandelstamScalarProducts/.tPoint,Modulus->$P];
DeleteList=Delete[{1,2,3,4,5,6,7,8,9,10,11},Position[tPoint[[All,1]],s16]];

diffInts=Table[
Expand[Expand[(applyDifferentialOperator[Collect[Coefficient[diffOp,ds[invariantsnodUpperCase[[var]]]],del[___]],Expand[PureIntegrals[[intpos]]/.tPoint[[DeleteCases[DeleteList,var+1]]]/.epsrepl,Modulus->$P]]/.tPoint/.root'[x_]:>1/(2 root[x]))/.root[x_]:>PowerMod[x,1/2,$P],Modulus->$P]/.loopMomentumScalarProductsT/.mandelstamScalarProductsT,Modulus->$P],{var,1,Length[invariantsnodUpperCase]},{intpos,1,Length[PureIntegrals]}];

diffIntsIBREDUCED=Expand[diffInts/.reductionRules/.ExtraRelations,Modulus->$P];
NonReduced=Cases[diffIntsIBREDUCED,_hxb,Infinity]//DeleteDuplicates;
NullInts2=Complement[NonReduced,SubLaporta];
diffIntsIBREDUCED=diffIntsIBREDUCED/.Thread[Rule[NullInts2,0]];


dhxbMatrices=Table[Normal[CoefficientArrays[diffIntsIBREDUCED[[i,j]],SubLaporta]/.{0}->{0,ConstantArray[0,Length[PureIntegralN]]}][[2]],{i,1,Length[invariants]-1},{j,1,Length[SubLaporta]}];

For[i=1,i<=Length[invariantsnod],i++,dhxbMatrices[[i]]=Expand[dhxbMatrices[[i]] . BasisChangMatInv,Modulus->$P];];


dhxbMatricesSaveFormat=Table[SparseArray[Flatten[dhxbMatrices[[i]]]],{i,1,Length[invariantsnod]}];
Return[{tPoint[[All,2]],dhxbMatricesSaveFormat}];
];
powerReconstruct[pars_,FitVecIn_]:=Module[{pars2,AnsatzPowersNumi,FitVec,AnsatzPowersDeno,MatNumiPart,onevc,parsCutOne,reconstructedFit,AnsatzPowersDenoOne,
fit,solrule,MatDenoPart,Solset,FinalMat,sol3,Ansatz,MaxPower},
	MaxPower=9;
AnsatzPowersNumi=FrobeniusSolve[{1,1},MaxPower];
	AnsatzPowersDeno=FrobeniusSolve[{1,1},MaxPower];
	pars2=pars/.{x_Integer:>{x,1}};
	(*MatNumiPart=Parallelize@Outer[m31ExpDot,pars2,AnsatzPowersNumi,1];*)

	MatNumiPart=Table[Table[Times@@PowerMod[pars2[[i]],AnsatzPowersNumi[[q]],$P],{q,1,Length[AnsatzPowersNumi]}],{i,1,Length[pars2]}];
	(*Table[Times@@PowerMod[pars2[[1]],AnsatzPowersNumi[[q]],$P],{q,1,Length[AnsatzPowersNumi]}]//Echo;*)
	
	FitVec=Expand[-FitVecIn,Modulus->$P];
	onevc=Table[{1},MaxPower+1];
	parsCutOne=MapThread[Join,{pars2,Transpose[{FitVec}]}];
	AnsatzPowersDenoOne=MapThread[Join,{AnsatzPowersDeno,onevc}];

	(*MatDenoPart=Parallelize@Outer[m31ExpDot,parsCutOne,AnsatzPowersDenoOne,1];*)
	
	MatDenoPart=Table[Table[Times@@PowerMod[parsCutOne[[i]],AnsatzPowersDenoOne[[j]],$P],{j,1,Length[AnsatzPowersDenoOne]}],{i,1,Length[parsCutOne]}];
	
	FinalMat=Expand[MapThread[Join,{MatNumiPart,MatDenoPart}],Modulus->$P];

	
	sol3=NullSpace[FinalMat,Modulus->$P];
	Ansatz=Sum[fac[n]t^n,{n,0,MaxPower}]/Sum[fac[n+MaxPower+1]t^n,{n,0,MaxPower}];
	Solset=Cases[Ansatz,fac[___],Infinity];
	solrule=Thread[Rule[Solset,sol3[[1]]]];
	fit=Together[Ansatz/.solrule]//ExpandDenominator//ExpandNumerator;
	Return[Expand[fit,Modulus->$P]]
]
powerReconstruct[pars_,FitVecIn_,MaxPower_]:=Module[{pars2,AnsatzPowersNumi,FitVec,AnsatzPowersDeno,MatNumiPart,onevc,parsCutOne,reconstructedFit,AnsatzPowersDenoOne,
fit,solrule,MatDenoPart,Solset,FinalMat,sol3,Ansatz},
AnsatzPowersNumi=FrobeniusSolve[{1,1},MaxPower];
	AnsatzPowersDeno=FrobeniusSolve[{1,1},MaxPower];
	pars2=pars/.{x_Integer:>{x,1}};
	(*MatNumiPart=Parallelize@Outer[m31ExpDot,pars2,AnsatzPowersNumi,1];*)

	MatNumiPart=Table[Table[Times@@PowerMod[pars2[[i]],AnsatzPowersNumi[[q]],$P],{q,1,Length[AnsatzPowersNumi]}],{i,1,Length[pars2]}];
	(*Table[Times@@PowerMod[pars2[[1]],AnsatzPowersNumi[[q]],$P],{q,1,Length[AnsatzPowersNumi]}]//Echo;*)
	
	FitVec=Expand[-FitVecIn,Modulus->$P];
	onevc=Table[{1},MaxPower+1];
	parsCutOne=MapThread[Join,{pars2,Transpose[{FitVec}]}];
	AnsatzPowersDenoOne=MapThread[Join,{AnsatzPowersDeno,onevc}];

	(*MatDenoPart=Parallelize@Outer[m31ExpDot,parsCutOne,AnsatzPowersDenoOne,1];*)
	
	MatDenoPart=Table[Table[Times@@PowerMod[parsCutOne[[i]],AnsatzPowersDenoOne[[j]],$P],{j,1,Length[AnsatzPowersDenoOne]}],{i,1,Length[parsCutOne]}];
	
	FinalMat=Expand[MapThread[Join,{MatNumiPart,MatDenoPart}],Modulus->$P];

	
	sol3=NullSpace[FinalMat,Modulus->$P];
	Ansatz=Sum[fac[n]t^n,{n,0,MaxPower}]/Sum[fac[n+MaxPower+1]t^n,{n,0,MaxPower}];
	Solset=Cases[Ansatz,fac[___],Infinity];
	solrule=Thread[Rule[Solset,sol3[[1]]]];
	fit=Together[Ansatz/.solrule]//ExpandDenominator//ExpandNumerator;
	Return[Expand[fit,Modulus->$P]]
]


GenerateAnsatz[vars_,degree_]:=Module[{Ansatz,AnsatzPowers},
	If[degree==0,1*fac[1]];
	AnsatzPowers=FrobeniusSolve[Table[1,{i,1,Length[vars]}],degree];
	Ansatz=Sum[fac[j]Product[vars[[i]]^AnsatzPowers[[j,i]],{i,1,Length[vars]}],{j,1,Length[AnsatzPowers]}];
	Return[Ansatz];
]

Get[FileNameJoin[{Directory[],"MonomialEvaluator.m"}]];
FFLib`$InstallDir="/home/dschiebe/Documents/packages/fflib-master";
Get[FFLib`$InstallDir<>"/src/fflib.m"];

ReconstructEntry[fitvecIN_,pars_,vars_,exp1_,exp2_]:=
Module[{variables,fitvec,flattpos,PowerMatrixUsed,AnsatzNumi,AnsatzDeno,varsA,AnsatzPowersNumi,AnsatzPowersDeno,AnsatzLengthN,AnsatzLengthD,AnsatzLength,Ansatz,Solset,AnsatzTerms,SampleSize,parsCut,dvec,nusedvars,dropl,MatNumiPart,onevc,parsCutOne,AnsatzPowersDenoOne,MatDenoPart,FinalMat,sol3,solrule,fit,reconstructedFit,usedSQRTS},
	If[exp1==0 && exp2==0,
	Return[0];
	];
	fitvec=fitvecIN;
	If[ Length[DeleteDuplicates[fitvec]]==1,Return[DeleteDuplicates[fitvec][[1]]]];
	variables=invariants[[2;;]];
	varsA=vars;
	AnsatzPowersNumi=FrobeniusSolve[Table[1,{i,1,Length[varsA]}],exp1];
	AnsatzPowersDeno=FrobeniusSolve[Table[1,{i,1,Length[varsA]}],exp2];
	AnsatzLengthN=Length[AnsatzPowersNumi];
	AnsatzLengthD=Length[AnsatzPowersDeno];
	AnsatzLength=AnsatzLengthN+AnsatzLengthD;
	
	SampleSize=Length[fitvec];


	AnsatzNumi=GenerateAnsatz[varsA,exp1];
	AnsatzDeno=GenerateAnsatz[varsA,exp2];
	AnsatzDeno=AnsatzDeno/.fac[n_]:>fac[n+AnsatzLengthN];
	Ansatz=AnsatzNumi/AnsatzDeno;
	Solset=Cases[Ansatz,fac[___],Infinity];
	
	fitvec=Drop[fitvec,SampleSize-(AnsatzLength)];
	parsCut=Drop[pars,SampleSize-(AnsatzLength)];
	
	
	nusedvars=Complement[variables,varsA];
	dropl=Table[Position[variables,nusedvars[[i]]],{i,1,Length[nusedvars]}][[All,All,1]];
		
	parsCut=Transpose[Delete[Transpose[parsCut],dropl]];
	
	MatNumiPart=Parallelize@Outer[m31ExpDot,parsCut,AnsatzPowersNumi,1];
	fitvec=Expand[-fitvec,Modulus->$P];
	onevc=Table[{1},AnsatzLengthD];
	parsCutOne=MapThread[Join,{parsCut,Transpose[{fitvec}]}];
	AnsatzPowersDenoOne=MapThread[Join,{AnsatzPowersDeno,onevc}];
	MatDenoPart=Parallelize@Outer[m31ExpDot,parsCutOne,AnsatzPowersDenoOne,1];
	FinalMat=MapThread[Join,{MatNumiPart,MatDenoPart}];
	sol3=NullSpace[FinalMat,Modulus->$P];
	solrule=Thread[Rule[Solset,sol3[[1]]]];
	fit=Together[Ansatz/.solrule,Modulus->$P];
	reconstructedFit=(Numerator[fit]/.n_Integer:>RationalReconstruct[n])/(Denominator[fit]/.n_Integer:>RationalReconstruct[n]);
	Return[reconstructedFit];
]
RationalReconstruct[a_]:=With[{v=LatticeReduce[{{a,1},{$P,0}}][[1]]},v[[1]]/v[[2]]];
dLog[letter_,variables_]:=Module[{varpos,dLogVals},
dLogVals=Table[D[Log[letter],variables[[varpos]]],{varpos,1,Length[variables]}];
Return[dLogVals];
]
getEvenDlogs[letters_,RationalFuncs_] := 
 Module[{ent,variables, Ansatz, eqns, solveSet, solution, vars},
 variables=Variables[letters];
  Ansatz = Sum[c[i] dLog[letters[[i]],variables], {i, 1, Length[letters]}];
  solveSet = 
   Sort[Cases[Ansatz, c[n_], Infinity] // DeleteDuplicates];
  eqns = Table[Ansatz[[i]] == RationalFuncs[[i]], {i, 1, Length[variables]}];
  eqns = Subtract @@@ DeleteCases[eqns, True];
  vars = Complement[Variables[eqns], solveSet];
  eqns = eqns // Together // Numerator;
  eqns = CoefficientRules[eqns, vars];
  eqns = Flatten[eqns][[All, 2]] == 0;
  solution = Together[Solve[eqns, solveSet]]//Echo;
  Return[{letters, solution[[1, All, 2]]}];
  ]
FancyDeleteDuplicates[delset_] := Module[{size, i, j, delset2, factor},
  delset2 = delset;
  size = Length[delset2];
  For[i = 1, i <= size, i++,
   For[j = i, j <= size, j++,
     factor = FullSimplify[delset2[[i]]/delset2[[j]]];
     If[NumberQ[factor],
      delset2[[j]] = delset2[[i]];
      ];
     ];
   ];
  Return[DeleteDuplicates[delset2]];
  ]
findDenominatorFactors[entry_]:=Module[{expr},
expr=Denominator[Factor[entry]];
If[expr==1,Return[{}]];
If[Length[expr]==0,
Return[{expr}]
];
If[expr[[0]]==Plus,
Return[{expr}]
];
Return[Table[expr[[i]],{i,1,Length[expr]}]];
]
