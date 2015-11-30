(* Content-type: application/vnd.wolfram.cdf.text *)

(*** Wolfram CDF File ***)
(* http://www.wolfram.com/cdf *)

(* CreatedBy='Mathematica 10.0' *)

(*************************************************************************)
(*                                                                       *)
(*  The Mathematica License under which this file was created prohibits  *)
(*  restricting third parties in receipt of this file from republishing  *)
(*  or redistributing it by any means, including but not limited to      *)
(*  rights management or terms of use, without the express consent of    *)
(*  Wolfram Research, Inc. For additional information concerning CDF     *)
(*  licensing and redistribution see:                                    *)
(*                                                                       *)
(*        www.wolfram.com/cdf/adopting-cdf/licensing-options.html        *)
(*                                                                       *)
(*************************************************************************)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[      1064,         20]
NotebookDataLength[      6215,        222]
NotebookOptionsPosition[      6437,        207]
NotebookOutlinePosition[      6780,        222]
CellTagsIndexPosition[      6737,        219]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell[BoxData[
 RowBox[{"TrigToExp", "[", 
  RowBox[{"ArcCosh", "[", "x", "]"}], "]"}]], "Input"],

Cell[BoxData[
 RowBox[{"Log", "[", 
  RowBox[{"x", "+", 
   RowBox[{
    SqrtBox[
     RowBox[{
      RowBox[{"-", "1"}], "+", "x"}]], " ", 
    SqrtBox[
     RowBox[{"1", "+", "x"}]]}]}], "]"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"FullSimplify", "[", 
  RowBox[{
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"Exp", "[", 
       RowBox[{"2", 
        RowBox[{"ArcCosh", "[", "x", "]"}]}], "]"}], "-", "1"}], ")"}], "-", 
    RowBox[{"2", " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        RowBox[{"(", 
         RowBox[{"x", "-", "1"}], ")"}], 
        RowBox[{"(", 
         RowBox[{"x", "+", "1"}], ")"}]}], "+", 
       RowBox[{"x", 
        SqrtBox[
         RowBox[{
          RowBox[{"(", 
           RowBox[{"x", "-", "1"}], ")"}], 
          RowBox[{"(", 
           RowBox[{"x", "+", "1"}], ")"}]}]]}]}], ")"}]}]}], ",", 
   RowBox[{"Assumptions", "\[Rule]", 
    RowBox[{"{", 
     RowBox[{"x", ">", "1"}], "}"}]}]}], "]"}]], "Input"],

Cell[BoxData["0"], "Output"]
}, Open  ]],

Cell[BoxData[{
 RowBox[{
  RowBox[{
   SubscriptBox["U", "val"], "=", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"21", "/", "5"}], ",", 
     RowBox[{"27", "/", "10"}], ",", 
     RowBox[{"51", "/", "10"}], ",", 
     RowBox[{"19", "/", "5"}]}], "}"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   SubscriptBox["\[Tau]", "val"], "=", 
   RowBox[{"1", "/", "11"}]}], ";"}]}], "Input"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  SubscriptBox["\[Lambda]", "val"], "=", 
  RowBox[{"ArcCosh", "[", 
   RowBox[{"Exp", "[", 
    RowBox[{
     SubscriptBox["U", "val"], 
     RowBox[{
      SubscriptBox["\[Tau]", "val"], "/", "2"}]}], "]"}], 
   "]"}]}], "\[IndentingNewLine]", 
 RowBox[{"N", "[", "%", "]"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"ArcCosh", "[", 
    SuperscriptBox["\[ExponentialE]", 
     RowBox[{"21", "/", "110"}]], "]"}], ",", 
   RowBox[{"ArcCosh", "[", 
    SuperscriptBox["\[ExponentialE]", 
     RowBox[{"27", "/", "220"}]], "]"}], ",", 
   RowBox[{"ArcCosh", "[", 
    SuperscriptBox["\[ExponentialE]", 
     RowBox[{"51", "/", "220"}]], "]"}], ",", 
   RowBox[{"ArcCosh", "[", 
    SuperscriptBox["\[ExponentialE]", 
     RowBox[{"19", "/", "110"}]], "]"}]}], "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
  "0.6377500503954737`", ",", "0.5056270073617484`", ",", 
   "0.7074957248939457`", ",", "0.6048110364152367`"}], "}"}]], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["fn", "base"], "=", 
   RowBox[{
    RowBox[{"NotebookDirectory", "[", "]"}], "<>", 
    RowBox[{"FileBaseName", "[", 
     RowBox[{"NotebookFileName", "[", "]"}], "]"}]}]}], ";"}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"save", " ", "reference", " ", "values", " ", "to", " ", "disk"}], 
   " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{"Export", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["fn", "base"], "<>", "\"\<_expVu.dat\>\""}], ",", 
      RowBox[{"N", "[", 
       RowBox[{
        RowBox[{"Flatten", "[", 
         RowBox[{"{", 
          RowBox[{
           RowBox[{"Exp", "[", 
            RowBox[{"-", 
             SubscriptBox["\[Lambda]", "val"]}], "]"}], ",", 
           RowBox[{"Exp", "[", 
            SubscriptBox["\[Lambda]", "val"], "]"}]}], "}"}], "]"}], ",", 
        RowBox[{"2", "MachinePrecision"}]}], "]"}], ",", "\"\<Real64\>\""}], 
     "]"}], ";"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"Export", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["fn", "base"], "<>", "\"\<_expVd.dat\>\""}], ",", 
      RowBox[{"N", "[", 
       RowBox[{
        RowBox[{"Flatten", "[", 
         RowBox[{"{", 
          RowBox[{
           RowBox[{"Exp", "[", 
            SubscriptBox["\[Lambda]", "val"], "]"}], ",", 
           RowBox[{"Exp", "[", 
            RowBox[{"-", 
             SubscriptBox["\[Lambda]", "val"]}], "]"}]}], "}"}], "]"}], ",", 
        RowBox[{"2", "MachinePrecision"}]}], "]"}], ",", "\"\<Real64\>\""}], 
     "]"}], ";"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"Export", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["fn", "base"], "<>", "\"\<_delta.dat\>\""}], ",", 
      RowBox[{"N", "[", 
       RowBox[{
        RowBox[{"Flatten", "[", 
         RowBox[{"{", 
          RowBox[{
           RowBox[{
            RowBox[{"Exp", "[", 
             RowBox[{"2", 
              SubscriptBox["\[Lambda]", "val"]}], "]"}], "-", "1"}], ",", 
           RowBox[{
            RowBox[{"Exp", "[", 
             RowBox[{
              RowBox[{"-", "2"}], 
              SubscriptBox["\[Lambda]", "val"]}], "]"}], "-", "1"}]}], "}"}], 
         "]"}], ",", 
        RowBox[{"2", "MachinePrecision"}]}], "]"}], ",", "\"\<Real64\>\""}], 
     "]"}], ";"}]}]}]], "Input"]
},
WindowSize->{1358, 764},
WindowMargins->{{279, Automatic}, {86, Automatic}},
FrontEndVersion->"10.0 for Microsoft Windows (64-bit) (July 1, 2014)",
StyleDefinitions->"Default.nb"
]
(* End of Notebook Content *)

(* Internal cache information *)
(*CellTagsOutline
CellTagsIndex->{}
*)
(*CellTagsIndex
CellTagsIndex->{}
*)
(*NotebookFileOutline
Notebook[{
Cell[CellGroupData[{
Cell[1486, 35, 96, 2, 31, "Input"],
Cell[1585, 39, 206, 8, 36, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[1828, 52, 762, 26, 39, "Input"],
Cell[2593, 80, 28, 0, 31, "Output"]
}, Open  ]],
Cell[2636, 83, 396, 13, 52, "Input"],
Cell[CellGroupData[{
Cell[3057, 100, 313, 10, 52, "Input"],
Cell[3373, 112, 509, 14, 33, "Output"],
Cell[3885, 128, 172, 4, 31, "Output"]
}, Open  ]],
Cell[4072, 135, 242, 7, 31, "Input"],
Cell[4317, 144, 2116, 61, 92, "Input"]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature vv0ZZxkOOIw4UCgw9rY1Zqst *)
