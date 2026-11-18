(* Run with wolframscript -file Model/check_sarah.wls /path/to/SARAH.m /path/to/models-parent *)
args = Rest[$ScriptCommandLine];
If[Length[args] != 2,
  Print["Usage: check_sarah.wls SARAH.m MODELS_PARENT"]; Exit[2]
];
If[!FileExistsQ[args[[1]]] ||
   !FileExistsQ[FileNameJoin[{args[[2]], "ZeeTypeII", "ZeeTypeII.m"}]],
  Print["Missing SARAH.m or ZeeTypeII model snapshot"]; Exit[2]
];
Get[args[[1]]];
SARAH[InputDirectories] = Prepend[SARAH[InputDirectories], args[[2]]];
Print["SARAH_MODEL_PARENT=", args[[2]]];
Start["ZeeTypeII"];
CheckModel[];
Print["CHECKMODEL_EXECUTED: examine warnings, rotations and vertices before claiming validation"];
