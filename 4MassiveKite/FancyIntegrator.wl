(* ::Package:: *)

SetDirectory["/home/dschiebe/Documents/packages/trager"];
FileNameJoin["Trager.wl"]//Get;
kernels=LaunchKernels[24];
QuietEcho[ParallelNeeds[FileNameJoin["Trager`"]]];
Print["HELLO FROM FANCY INTEGRATE"]
RationalReconstruct[a_]:=With[{v=LatticeReduce[{{a,1},{$P,0}}][[1]]},v[[1]]/v[[2]]];
SetDirectory["/home/dschiebe/Documents/phd_files/Reconstruct Pentagon"];
Get[FileNameJoin["MonomialEvaluator.m"]];
FFLib`$InstallDir="/home/dschiebe/Documents/packages/fflib-master";
Get[FFLib`$InstallDir<>"/src/fflib.m"];




SampleIntagrateLine[Integrand_,qring_,rootset_,Intvar_,var_,kinp0_,kinp1_]:=Module[{roots,IntFuncRoot,tval},
(*kinPoint = #->RandomInteger[{-10000,10000}]&/@var;*)
tval=RandomInteger[{0,1000000}];
kinPoint=Thread[Rule[var,kinp0+tval*kinp1]];
shift={Intvar->Intvar-23};
shiftBack={Intvar->Intvar+23};
TragerLogArg=Quiet[TragerIntegrate[Integrand /. kinPoint/.shift, qring /. kinPoint/.shift,rootset , Intvar, 3]/.Log[x_]:>x/.shiftBack//Expand//Together];

NUMISCALE=CoefficientRules[Numerator[TragerLogArg],Join[{Intvar},rootset]][[All,2]][[1]];
DENSCALE=CoefficientRules[Denominator[TragerLogArg],Join[{Intvar},rootset]][[All,2]][[1]];
NumArg=Numerator[TragerLogArg]/NUMISCALE;
DenArg=Denominator[TragerLogArg]/DENSCALE;
CoefNumi=CoefficientRules[NumArg,Join[{Intvar},rootset]][[All,2]];
CoefDeno=CoefficientRules[DenArg,Join[{Intvar},rootset]][[All,2]];
pars=kinPoint[[All,2]];
Return[{CoefNumi,CoefDeno,tval}];

]


SampleIntagrate[Integrand_,qring_,rootset_,Intvar_,var_,numForm_,denForm_]:=Module[{roots,IntFuncRoot,tval},
kinPoint = #->RandomInteger[{0,1000000}]&/@var;
shiftval=RandomInteger[{-10000,10000}];
shift={Intvar->Intvar-shiftval};
shiftBack={Intvar->Intvar+shiftval};
TragerLogArg=Quiet[TragerIntegrate[Integrand /. kinPoint/.shift//ExpandDenominator//ExpandNumerator, qring /. kinPoint/.shift//Expand,rootset , Intvar, 3]/.Log[x_]:>x/.shiftBack//Expand//Together];
NUMISCALE=CoefficientRules[Numerator[TragerLogArg],Join[{Intvar},rootset]][[All,2]][[1]];
DENSCALE=CoefficientRules[Denominator[TragerLogArg],Join[{Intvar},rootset]][[All,2]][[1]];
NumArg=Numerator[TragerLogArg]/NUMISCALE;
DenArg=Denominator[TragerLogArg]/DENSCALE;
CoefNumi=CoefficientRules[NumArg,Join[{Intvar},rootset]][[All,2]];
CoefDeno=CoefficientRules[DenArg,Join[{Intvar},rootset]][[All,2]];
If[Length[CoefNumi]=!=numForm || Length[CoefDeno]=!=denForm,
Return[SampleIntagrate[Integrand,qring,rootset,Intvar,var,numForm,denForm]];
];
pars=kinPoint[[All,2]];
Return[{CoefNumi,CoefDeno,pars}];

]


GenerateAnsatz[vars_,degree_]:=Module[{Ansatz,AnsatzPowers},
If[degree==0,1*fac[1]];
AnsatzPowers=FrobeniusSolve[Table[1,{i,1,Length[vars]}],degree];
Ansatz=Sum[fac[j]Product[vars[[i]]^AnsatzPowers[[j,i]],{i,1,Length[vars]}],{j,1,Length[AnsatzPowers]}];
Return[Ansatz];
]


GetPowers[Integrand_,qring_,rootset_,Intvar_,var_]:=Module[{Ansatz,AnsatzPowers,Pars,FitVec},
kinp0=Table[RandomInteger[10000],{i,1,Length[var]}];
kinp1=Table[RandomInteger[10000],{i,1,Length[var]}];
Nums=ParallelTable[SampleIntagrateLine[Integrand,qring,rootset,Intvar,var,kinp0,kinp1],{i,1,22}];
NumiList=Nums[[All,1]];
DenoList=Nums[[All,2]];
Pars=Nums[[All,3]];
PowersNum={};
For[Pos=1,Pos<=Length[NumiList[[1]]],Pos++,
	FitVec=NumiList[[All,Pos]];
(*	AnsatzNumi=Sum[fac[n]t^n,{n,0,5}];
	AnsatzLengthN=Length[AnsatzNumi];
	AnsatzDeno=Sum[fac[n]t^n,{n,0,5}];
	AnsatzDeno=AnsatzDeno/.fac[n_]:>fac[n+AnsatzLengthN];
	Ansatz=AnsatzNumi/AnsatzDeno;
	eqn=Table[(Ansatz/.Thread[Rule[t,Pars[[i]]]])==FitVec[[i]],{i,1,18}];
	solset=Cases[Ansatz,fac[n_],Infinity];
	sol=Quiet[Solve[eqn,solset][[1]]]//EchoTiming;
	fit=Together[Ansatz/.sol];
    Print[fit//ExpandNumerator//ExpandDenominator];*)
    
    fit=powerReconstruct[Pars,FitVec]//EchoTiming;
    fit//ExpandNumerator//ExpandDenominator//Echo;
    
	maxNumi=Exponent[Numerator[fit],t];
	maxDeno=Exponent[Denominator[fit],t];
	PowersNum=Join[PowersNum,{{maxNumi,maxDeno}}];
];
PowersDen={};
For[Pos=1,Pos<=Length[DenoList[[1]]],Pos++,
	FitVec=DenoList[[All,Pos]];
(*    AnsatzNumi=Sum[fac[n]t^n,{n,0,5}];
	AnsatzLengthN=Length[AnsatzNumi];
	AnsatzDeno=Sum[fac[n]t^n,{n,0,5}];
	AnsatzDeno=AnsatzDeno/.fac[n_]:>fac[n+AnsatzLengthN];
	Ansatz=AnsatzNumi/AnsatzDeno;
	eqn=Table[(Ansatz/.Thread[Rule[t,Pars[[i]]]])==FitVec[[i]],{i,1,18}];
	solset=Cases[Ansatz,fac[n_],Infinity];
	sol=Quiet[Solve[eqn,solset][[1]]]//EchoTiming;
	fit=Together[Ansatz/.sol]//EchoTiming;
	Print[fit//ExpandNumerator//ExpandDenominator];*)
	
	fit=powerReconstruct[Pars,FitVec]//EchoTiming;
    fit//ExpandNumerator//ExpandDenominator//Echo;
	
	maxNumi=Exponent[Numerator[fit],t];
	maxDeno=Exponent[Denominator[fit],t];
	PowersDen=Join[PowersDen,{{maxNumi,maxDeno}}];
];
Return[{PowersNum,PowersDen}]
]


ansatzTermsN[vars_,exp_]:=Module[{},
If[vars==0,Return[0]];
If[exp==0,Return[0]];
Return[Binomial[vars+exp-1,exp]];
]


FancyIntegrate[IntFunc_,Intvar_,var_]:=Module[{roots,IntFuncRoot,rootset},
If[IntFunc==0,
Return[0];
];
roots=Cases[IntFunc,root[x_],Infinity]//DeleteDuplicates;
rootset=Table[r[i],{i,1,Length[roots]}];
rootrepls=Thread[Rule[roots,rootset]];
roots=roots/.root[x_]:>x;
qring=Table[r[i]^2-roots[[i]],{i,1,Length[roots]}];
Integrand=IntFunc/.rootrepls;
kinPoint = #->RandomInteger[{0,1000000}]&/@var;
shift={Intvar->Intvar-23};
shiftBack={Intvar->Intvar+23};
TragerLog=Quiet[TragerIntegrate[Integrand /. kinPoint/.shift//ExpandDenominator//ExpandNumerator//Echo, qring /. kinPoint/.shift//Expand//Echo,rootset//Echo , Intvar, 3]/.shiftBack]//Echo;
TragerLogArg=Cases[{TragerLog},Log[___],Infinity][[1]]/.Log[x_]:>x//Expand//Together//Echo;

(*Print[TragerLog];
Print[TragerLogArg];
*)
Difference=((D[(TragerLog/.Reverse/@rootrepls)/.root[x_]:>Sqrt[x]/. kinPoint,Intvar])-(IntFunc/.root[x_]:>Sqrt[x]/. kinPoint))//FullSimplify//Echo;

If[Difference =!=0,
	Return[0]
];



NUMISCALE=CoefficientRules[Numerator[TragerLogArg],Join[{Intvar},rootset]][[All,2]][[1]];
DENSCALE=CoefficientRules[Denominator[TragerLogArg],Join[{Intvar},rootset]][[All,2]][[1]];
NumArg=Numerator[TragerLogArg]/NUMISCALE;
DenArg=Denominator[TragerLogArg]/DENSCALE;

CoefPowersNumi=CoefficientRules[NumArg,Join[{Intvar},rootset]][[All,1]];
CoefPowersDeno=CoefficientRules[DenArg,Join[{Intvar},rootset]][[All,1]];

CoefNumi=CoefficientRules[NumArg,Join[{Intvar},rootset]][[All,2]];
CoefDeno=CoefficientRules[DenArg,Join[{Intvar},rootset]][[All,2]];

numForm=Length[CoefNumi];
denForm=Length[CoefDeno];

Powers=GetPowers[Integrand,qring,rootset,Intvar,var];

(*Print[Powers];*)

PowersNum=Powers[[1]]//Echo;
(*PowersNum=PowersNum[[All,All,1]];*)
PowersDen=Powers[[2]]//Echo;
(*PowersDen=PowersDen[[All,All,1]];*)

AnsatzTermsNum=Table[ansatzTermsN[Length[var],PowersNum[[k,1]]]+ansatzTermsN[Length[var],PowersNum[[k,2]]],{k,1,Length[PowersNum]}];
AnsatzTermsDen=Table[ansatzTermsN[Length[var],PowersDen[[k,1]]]+ansatzTermsN[Length[var],PowersDen[[k,2]]],{k,1,Length[PowersDen]}];
MaxTermsNeeded=Max[Join[AnsatzTermsNum,AnsatzTermsDen]]//Echo;
If[MaxTermsNeeded>200,
Print[MaxTermsNeeded];
];


Nums=ParallelTable[SampleIntagrate[Integrand,qring,rootset,Intvar,var,numForm,denForm],{i,1,MaxTermsNeeded}];

NumiList=Nums[[All,1]]//Echo;
DenoList=Nums[[All,2]]//Echo;
Pars=Nums[[All,3]];
var//Echo;
FunctionsNumi={};
For[Pos=1,Pos<=Length[NumiList[[1]]],Pos++,
	FitVec=NumiList[[All,Pos]]//Echo;
	Pars//Echo;
	Ansatz1=GenerateAnsatz[var,PowersNum[[Pos,1]]]//Echo;
	Ansatz2=GenerateAnsatz[var,PowersNum[[Pos,2]]]/.fac[n_]:>fac[n+Length[Ansatz1]]//Echo;
	Ansatz=Ansatz1/Ansatz2;
	eqn=Table[(Ansatz/.Thread[Rule[var,Pars[[i]]]])==FitVec[[i]],{i,1,MaxTermsNeeded}];
	solset=Cases[Ansatz,fac[n_],Infinity];
	sol=Quiet[Solve[eqn,solset][[1]]]//Echo;
	FunctionsNumi=Join[FunctionsNumi,{Together[Ansatz/.sol]}];
];
FunctionsDeno={};
For[Pos=1,Pos<=Length[DenoList[[1]]],Pos++,
FitVec=DenoList[[All,Pos]];
Ansatz1=GenerateAnsatz[var,PowersDen[[Pos,1]]];
Ansatz2=GenerateAnsatz[var,PowersDen[[Pos,2]]]/.fac[n_]:>fac[n+Length[Ansatz1]];
Ansatz=Ansatz1/Ansatz2;
eqn=Table[(Ansatz/.Thread[Rule[var,Pars[[i]]]])==FitVec[[i]],{i,1,MaxTermsNeeded}];
solset=Cases[Ansatz,fac[n_],Infinity];
sol=Quiet[Solve[eqn,solset][[1]]];
FunctionsDeno=Join[FunctionsDeno,{Together[Ansatz/.sol]}];
];
newvars=Join[{Intvar},rootset];
NumeratorRec=Sum[FunctionsNumi[[i]]*Times@@Power[newvars,CoefPowersNumi[[i]]],{i,1,Length[FunctionsNumi]}];
DenominatorRec=Sum[FunctionsDeno[[i]]*Times@@Power[newvars,CoefPowersDeno[[i]]],{i,1,Length[FunctionsDeno]}];

NumeratorRec=Collect[Expand[NumeratorRec/.Reverse/@rootrepls],root[___],Together];
DenominatorRec=Collect[Expand[DenominatorRec/.Reverse/@rootrepls],root[___],Together];
FullArg=Together[NumeratorRec/DenominatorRec];
FullNumerator=Collect[Numerator[FullArg],root[___]];
FullDenominator=Collect[Denominator[FullArg],root[___]];

Return[Log[FullNumerator/FullDenominator]];
]


powerReconstruct[pars_,FitVecIn_]:=Module[{pars2,AnsatzPowersNumi,FitVec,AnsatzPowersDeno,MatNumiPart,onevc,parsCutOne,reconstructedFit,AnsatzPowersDenoOne,fit,solrule,MatDenoPart,Solset,FinalMat,sol3,Ansatz},
	AnsatzPowersNumi=FrobeniusSolve[{1,1},10];
	AnsatzPowersDeno=FrobeniusSolve[{1,1},10];
	pars2=pars/.{x_Integer:>{x,1}};
	(*MatNumiPart=Parallelize@Outer[m31ExpDot,pars2,AnsatzPowersNumi,1];*)

	MatNumiPart=ParallelTable[Table[Times@@PowerMod[pars2[[i]],AnsatzPowersNumi[[q]],$P],{q,1,Length[AnsatzPowersNumi]}],{i,1,Length[pars2]}];
	(*Table[Times@@PowerMod[pars2[[1]],AnsatzPowersNumi[[q]],$P],{q,1,Length[AnsatzPowersNumi]}]//Echo;*)
	
	FitVec=Expand[-FitVecIn,Modulus->$P];
	onevc=Table[{1},11];
	parsCutOne=MapThread[Join,{pars2,Transpose[{FitVec}]}];
	AnsatzPowersDenoOne=MapThread[Join,{AnsatzPowersDeno,onevc}];

	(*MatDenoPart=Parallelize@Outer[m31ExpDot,parsCutOne,AnsatzPowersDenoOne,1];*)
	
	MatDenoPart=ParallelTable[Table[Times@@PowerMod[parsCutOne[[i]],AnsatzPowersDenoOne[[j]],$P],{j,1,Length[AnsatzPowersDenoOne]}],{i,1,Length[parsCutOne]}];
	
	FinalMat=Expand[MapThread[Join,{MatNumiPart,MatDenoPart}],Modulus->$P];

	
	sol3=NullSpace[FinalMat,Modulus->$P];
	Ansatz=Sum[fac[n]t^n,{n,0,10}]/Sum[fac[n+11]t^n,{n,0,10}];
	Solset=Cases[Ansatz,fac[___],Infinity];
	solrule=Thread[Rule[Solset,sol3[[1]]]];
	fit=Together[Ansatz/.solrule]//ExpandDenominator//ExpandNumerator;
     reconstructedFit=(Numerator[fit]/.n_Integer:>RationalReconstruct[n])/(Denominator[fit]/.n_Integer:>RationalReconstruct[n]);
	Return[reconstructedFit]
]
