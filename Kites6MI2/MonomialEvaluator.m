(* ::Package:: *)

(* ::Subsubsection:: *)
(*Fast eval tech*)


m31Pow = With[{p = 2^31 - 1}, Compile[{{x, _Integer}, {n, _Integer}}, Module[{xn = 1}, Do[xn = Mod[xn*x, p], {n}]; xn]]];


m31ExpDot = With[{p=2^31-1}, 
    Compile[{{val, _Integer, 1}, {exp, _Integer, 1}}, 
        Module[{res = 1}, Do[res = Mod[res*m31Pow[val[[i]], exp[[i]]], p], {i,Length@exp}]; res], 
        CompilationOptions -> {"InlineExternalDefinitions" -> True}, RuntimeOptions -> "Speed"
    ]
 ];



