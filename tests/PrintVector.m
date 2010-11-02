(* ::Package:: *)

BeginPackage["PrintVector`"];

PrintVector::usage = "Print a vector."

Begin["`Private`"];

FormatRealNumber[x_]:="K(" <> ToString[NumberForm[x,NumberSigns->{"-","+"}]] <> ")";
FormatNumber[x_]:=FormatRealNumber[x]/;Im[x]==0;
FormatNumber[x_]:=FormatRealNumber[Re[x]]<>" + "<>FormatRealNumber[Im[x]] <> " * I";
TypeString[x_]:="R"/;Max[Abs[Im[x]]]==0;
TypeString[x_]:="C";

PrintVector[x_,name_]:=Module[{M=Length[x], s = "static const " <>TypeString[x] <>" "<> name <> "[] = \n{\n"},
For[j=1,j<=M,j++,s=s<>"  " <>FormatNumber[x[[j]]]<>",\n"];s=s<>"};";Print[s]];

End[];

EndPackage[];
