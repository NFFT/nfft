(* ::Package:: *)

BeginPackage["PrintVector`"];

FormatVector::usage = "Print a vector in C format."
FormatVectorRaw::usage = "Print a vector in raw format."
FormatInteger::usage = "Print a vector in raw format."
FormatIntegerVector::usage = "Print a vector in raw format."
FormatIntegerRaw::usage = "Print a vector in raw format."
FormatIntegerVectorRaw::usage = "Print a vector in raw format."

Begin["`Private`"];

FormatRealNumber[x_]:="K(" <> ToString[CForm[x]] <> ")";
FormatNumber[x_]:="\""<>x<>"\""/;StringQ[x];
FormatNumber[x_]:=FormatRealNumber[x]/;Abs[Im[x]]==0;
FormatNumber[x_]:=FormatRealNumber[Re[x]]<>" + "<>FormatRealNumber[Im[x]] <> " * I";
TypeString[x_]:="char*"/;StringQ[x[[1]]];
TypeString[x_]:="R"/;Max[Abs[Im[x]]]==0;
TypeString[x_]:="C";

FormatVector[x_,name_,type_:TypeString,formatter_:None]:=Module[{M=Length[x], s = "static const " <> If[StringQ[type],type,TypeString[x]] <>" "<> name <> "[] = \n{\n"},
For[j=1,j<=M,j++,s=s<>"  " <> If[TrueQ[formatter==None],FormatNumber[x[[j]]],formatter[x[[j]]]] <>",\n"];s=s<>"};";Return[s]];

FormatRealNumberRaw[x_]:=ToString[CForm[x]];
FormatNumberRaw[x_]:=FormatRealNumberRaw[x]/;Abs[Im[x]]==0;
FormatNumberRaw[x_]:=FormatRealNumberRaw[Re[x]]<>" "<>FormatRealNumberRaw[Im[x]];

FormatVectorRaw[x_]:=Module[{M=Length[x], s = ""},
For[j=1,j<=M,j++,s=s<>FormatNumberRaw[x[[j]]]<>"\n"];Return[s]];

FormatIntegerInternal[x_]:=ToString[NumberForm[IntegerPart[x],NumberSigns->{"-",""}]];
FormatIntegerVector[x_,name_]:=Module[{M=Length[x], s = "static const INT "<> name <> "[] = \n{\n"},
For[j=1,j<=M,j++,s=s<>"  " <>FormatIntegerInternal[x[[j]]]<>",\n"];s=s<>"};";Return[s]];
FormatInteger[x_,name_]:=Module[{M=Length[x], s = "static const INT "<> name <> " = "<>FormatIntegerInternal[x]<>";\n"},Return[s]];
FormatIntegerVectorRaw[x_]:=Module[{M=Length[x], s = ""},
For[j=1,j<=M,j++,s=s<>FormatIntegerInternal[x[[j]]]<>"\n"];Return[s]];
FormatIntegerRaw[x_]:=Module[{M=Length[x], s = FormatIntegerInternal[x]<>"\n"},Return[s]];

End[];

EndPackage[];












