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
NotebookDataLength[     39323,       1110]
NotebookOptionsPosition[     36840,       1004]
NotebookOutlinePosition[     37184,       1019]
CellTagsIndexPosition[     37141,       1016]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["Rectangular lattice", "Subsection"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"lattice", " ", "dimensions"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["n", "x"], "=", "4"}], ";"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{
     SubscriptBox["n", "y"], "=", "6"}], ";"}]}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"lattice", " ", "sites"}], ",", " ", 
    RowBox[{
     RowBox[{"in", " ", "column"}], "-", 
     RowBox[{"major", " ", 
      RowBox[{"(", "Fortran", ")"}], " ", "order"}]}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["latt", "sites"], "=", 
    RowBox[{"Flatten", "[", 
     RowBox[{
      RowBox[{"Transpose", "[", 
       RowBox[{"Outer", "[", 
        RowBox[{
         RowBox[{
          RowBox[{"{", "##", "}"}], "&"}], ",", 
         RowBox[{"Range", "[", 
          RowBox[{"0", ",", 
           RowBox[{
            SubscriptBox["n", "x"], "-", "1"}]}], "]"}], ",", 
         RowBox[{"Range", "[", 
          RowBox[{"0", ",", 
           RowBox[{
            SubscriptBox["n", "y"], "-", "1"}]}], "]"}]}], "]"}], "]"}], ",", 
      "1"}], "]"}]}], "\[IndentingNewLine]", 
   RowBox[{"Dimensions", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"3", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"3", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"3", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"3", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"3", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "5"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "5"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "5"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"3", ",", "5"}], "}"}]}], "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"24", ",", "2"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"total", " ", "number", " ", "of", " ", "lattice", " ", "sites"}], 
   " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   SubscriptBox["n", "sites"], "=", 
   RowBox[{"Length", "[", 
    SubscriptBox["latt", "sites"], "]"}]}]}]], "Input"],

Cell[BoxData["24"], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"nearest", " ", 
    RowBox[{"(", "direct", ")"}], " ", "neighbors"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["neigh", "nearest"], "=", 
    RowBox[{"Outer", "[", 
     RowBox[{
      RowBox[{
       RowBox[{"If", "[", 
        RowBox[{
         RowBox[{
          RowBox[{"Norm", "[", 
           RowBox[{"{", 
            RowBox[{
             RowBox[{"Mod", "[", 
              RowBox[{
               RowBox[{"First", "[", 
                RowBox[{"#1", "-", "#2"}], "]"}], ",", 
               SubscriptBox["n", "x"], ",", 
               RowBox[{"-", "1"}]}], "]"}], ",", 
             RowBox[{"Mod", "[", 
              RowBox[{
               RowBox[{"Last", "[", 
                RowBox[{"#1", "-", "#2"}], "]"}], ",", 
               SubscriptBox["n", "y"], ",", 
               RowBox[{"-", "1"}]}], "]"}]}], "}"}], "]"}], "\[Equal]", "1"}],
          ",", "1", ",", "0"}], "]"}], "&"}], ",", 
      SubscriptBox["latt", "sites"], ",", 
      SubscriptBox["latt", "sites"], ",", "1"}], "]"}]}], ";"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"check", ":", " ", "symmetric"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["neigh", "nearest"], ";"}], "\[IndentingNewLine]", 
   RowBox[{"Norm", "[", 
    RowBox[{"%", "-", 
     RowBox[{"Transpose", "[", "%", "]"}]}], "]"}]}]}]], "Input"],

Cell[BoxData["0"], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"check", ":", " ", 
    RowBox[{"every", " ", "site", " ", "has", " ", "4", " ", "neighbors"}]}], 
   " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{"Norm", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{"Total", "/@", 
      SubscriptBox["neigh", "nearest"]}], ")"}], "-", "4"}], "]"}]}]], "Input"],

Cell[BoxData["0"], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"next", "-", 
    RowBox[{"nearest", " ", 
     RowBox[{"(", "diagonal", ")"}], " ", "neighbors"}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["neigh", "nextnearest"], "=", 
    RowBox[{"Outer", "[", 
     RowBox[{
      RowBox[{
       RowBox[{"If", "[", 
        RowBox[{
         RowBox[{
          RowBox[{"Norm", "[", 
           RowBox[{"{", 
            RowBox[{
             RowBox[{"Mod", "[", 
              RowBox[{
               RowBox[{"First", "[", 
                RowBox[{"#1", "-", "#2"}], "]"}], ",", 
               SubscriptBox["n", "x"], ",", 
               RowBox[{"-", "1"}]}], "]"}], ",", 
             RowBox[{"Mod", "[", 
              RowBox[{
               RowBox[{"Last", "[", 
                RowBox[{"#1", "-", "#2"}], "]"}], ",", 
               SubscriptBox["n", "y"], ",", 
               RowBox[{"-", "1"}]}], "]"}]}], "}"}], "]"}], "\[Equal]", 
          SqrtBox["2"]}], ",", "1", ",", "0"}], "]"}], "&"}], ",", 
      SubscriptBox["latt", "sites"], ",", 
      SubscriptBox["latt", "sites"], ",", "1"}], "]"}]}], ";"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"check", ":", " ", "symmetric"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["neigh", "nextnearest"], ";"}], "\[IndentingNewLine]", 
   RowBox[{"Norm", "[", 
    RowBox[{"%", "-", 
     RowBox[{"Transpose", "[", "%", "]"}]}], "]"}]}]}]], "Input"],

Cell[BoxData["0"], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"check", ":", " ", 
    RowBox[{
     RowBox[{"every", " ", "site", " ", "has", " ", "4", " ", "next"}], "-", 
     RowBox[{"nearest", " ", "neighbors"}]}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{"Norm", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{"Total", "/@", 
      SubscriptBox["neigh", "nextnearest"]}], ")"}], "-", "4"}], 
   "]"}]}]], "Input"],

Cell[BoxData["0"], "Output"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Kinetic energy operator", "Subsection"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"time", " ", "step"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[Tau]", "val"], "=", 
    FractionBox["1", "8"]}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"t", " ", "hopping", " ", "parameter"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["t", "val"], "=", "1"}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    SuperscriptBox["t", "\[Prime]",
     MultilineFunction->None], " ", "hopping", " ", "parameter"}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["tp", "val"], "=", 
    FractionBox["1", "5"]}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"chemical", " ", "potential"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[Mu]", "val"], "=", 
    RowBox[{"-", 
     FractionBox["3", "17"]}]}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["k", "val"], "=", 
   RowBox[{"-", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{
       SubscriptBox["t", "val"], 
       SubscriptBox["neigh", "nearest"]}], "+", 
      RowBox[{
       SubscriptBox["tp", "val"], 
       SubscriptBox["neigh", "nextnearest"]}], "+", 
      RowBox[{
       SubscriptBox["\[Mu]", "val"], 
       RowBox[{"IdentityMatrix", "[", 
        SubscriptBox["n", "sites"], "]"}]}]}], ")"}]}]}], ";"}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   SuperscriptBox["\[ExponentialE]", 
    RowBox[{
     RowBox[{"-", "\[Tau]"}], " ", "k"}]], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["exp", "\[Tau]k"], "=", 
    RowBox[{"MatrixExp", "[", 
     RowBox[{"N", "[", 
      RowBox[{
       RowBox[{
        RowBox[{"-", 
         SubscriptBox["\[Tau]", "val"]}], 
        SubscriptBox["k", "val"]}], ",", 
       RowBox[{"6", "MachinePrecision"}]}], "]"}], "]"}]}], ";"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"N", "[", 
    RowBox[{
     SubscriptBox["exp", "\[Tau]k"], "\[LeftDoubleBracket]", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{"1", ",", "2", ",", 
        RowBox[{"-", "1"}]}], "}"}], ",", 
      RowBox[{"{", 
       RowBox[{"1", ",", "2", ",", 
        RowBox[{"-", "1"}]}], "}"}]}], "\[RightDoubleBracket]"}], "]"}], "//",
    "MatrixForm"}]}]], "Input"],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"1.0119815705757973`", "0.13216793496262425`", "0.041472525678860526`"},
     {"0.13216793496262425`", "1.0119815705757973`", "0.008505393159145882`"},
     {"0.041472525678860526`", "0.008505393159145882`", "1.0119815705757973`"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Hubbard-Stratonovich field", "Subsection"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"SeedRandom", "[", "42", "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"Block", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"\[Beta]", "=", "4"}], ",", "L"}], "}"}], ",", 
    RowBox[{
     RowBox[{"L", "=", 
      FractionBox["\[Beta]", 
       SubscriptBox["\[Tau]", "val"]]}], ";", 
     RowBox[{
      SubscriptBox["s", "val"], "=", 
      RowBox[{
       RowBox[{"2", 
        RowBox[{"RandomInteger", "[", 
         RowBox[{
          RowBox[{"{", 
           RowBox[{"0", ",", "1"}], "}"}], ",", 
          RowBox[{"{", 
           RowBox[{
            SubscriptBox["n", "sites"], ",", "L"}], "}"}]}], "]"}]}], "-", 
       "1"}]}]}]}], "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{"Dimensions", "[", 
  SubscriptBox["s", "val"], "]"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"24", ",", "32"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "examples", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["s", "val"], "\[LeftDoubleBracket]", 
    RowBox[{";;", ",", "1"}], "\[RightDoubleBracket]"}], 
   "\[IndentingNewLine]", 
   RowBox[{
    SubscriptBox["s", "val"], "\[LeftDoubleBracket]", 
    RowBox[{";;", ",", "2"}], "\[RightDoubleBracket]"}], 
   "\[IndentingNewLine]", 
   RowBox[{
    SubscriptBox["s", "val"], "\[LeftDoubleBracket]", 
    RowBox[{";;", ",", 
     RowBox[{"-", "1"}]}], "\[RightDoubleBracket]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}]}], "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}]}], "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}]}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"Flatten", "[", 
    RowBox[{"Transpose", "[", 
     SubscriptBox["s", "val"], "]"}], "]"}], "/.", 
   RowBox[{"{", 
    RowBox[{"1", "\[Rule]", "0"}], "}"}]}], "/.", 
  RowBox[{"{", 
   RowBox[{
    RowBox[{"-", "1"}], "\[Rule]", "1"}], "}"}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
  "1", ",", "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", 
   ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", 
   "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", 
   ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", 
   ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", 
   "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", 
   ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", 
   "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", 
   ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", 
   "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", 
   ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", 
   "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", 
   ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", 
   "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", 
   ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", 
   "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", 
   ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", 
   "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", 
   ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", 
   "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", 
   ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", 
   "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", 
   ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", 
   "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", 
   ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", 
   "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", 
   ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", 
   "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", 
   ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", 
   ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", 
   "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", 
   ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", 
   "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", 
   ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", 
   "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", 
   ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", 
   "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", 
   ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", 
   "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", 
   ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", 
   "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", 
   ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", 
   ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", 
   "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", 
   ",", "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", 
   "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", 
   ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", 
   ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", 
   "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", 
   ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", 
   "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", 
   ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", 
   "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", 
   ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", 
   "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", 
   ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", 
   ",", "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", 
   "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", ",", "1", 
   ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", 
   "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", 
   ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", 
   "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", 
   ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", 
   ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", 
   "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", 
   ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", 
   "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", 
   ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", 
   ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", 
   "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", 
   ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", 
   "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", 
   ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", 
   "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", 
   ",", "0", ",", "0", ",", "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", 
   "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", 
   ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", 
   "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", 
   ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", 
   "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", 
   ",", "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", 
   "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", "1", 
   ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", 
   "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", 
   ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", 
   "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", 
   ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", 
   "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", 
   ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", 
   "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", 
   ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", 
   "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", 
   ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", 
   "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", 
   ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", 
   "1", ",", "1", ",", "1"}], "}"}]], "Output"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Calculate time flow map", "Subsection"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
   "time", " ", "flow", " ", "map", " ", "generated", " ", "by", " ", "the", 
    " ", "Hubbard", " ", "Hamiltonian"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"HubbardTimeFlowMap", "[", 
    RowBox[{"\[Lambda]s_", ",", "exp\[Tau]k_"}], "]"}], ":=", 
   RowBox[{"Fold", "[", 
    RowBox[{
     RowBox[{
      RowBox[{
       RowBox[{"DiagonalMatrix", "[", 
        RowBox[{"Exp", "[", 
         RowBox[{"-", "#2"}], "]"}], "]"}], ".", "exp\[Tau]k", ".", "#1"}], 
      "&"}], ",", 
     RowBox[{"IdentityMatrix", "[", 
      RowBox[{"Length", "[", "exp\[Tau]k", "]"}], "]"}], ",", 
     RowBox[{"Transpose", "[", "\[Lambda]s", "]"}]}], "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["\[Lambda]", "val"], "=", 
   RowBox[{"3", "/", "4"}]}], ";"}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{
   SubscriptBox["A", "val"], "=", 
   RowBox[{"HubbardTimeFlowMap", "[", 
    RowBox[{
     RowBox[{
      SubscriptBox["\[Lambda]", "val"], 
      SubscriptBox["s", "val"]}], ",", 
     SubscriptBox["exp", "\[Tau]k"]}], "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{"Dimensions", "[", "%", "]"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"24", ",", "24"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["A", "val"], "\[LeftDoubleBracket]", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"1", ",", "2", ",", 
       RowBox[{"-", "1"}]}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"1", ",", "2", ",", 
       RowBox[{"-", "1"}]}], "}"}]}], "\[RightDoubleBracket]"}], "//", 
   "MatrixForm"}]}]], "Input"],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"3.015038897696885742974114592393091781469937294075130549756527248625540\
467374055602272226379385778772246113828240927677372320486`94.22238864282612*^\
11", "1.7682456718719777407818750859017402712053820665405724418030263089102253\
11374456390617460611328377907473471855624650242765799349`94.222388642826*^11",
       "3.57758094854373779430800408415395608919797230988843691410722857861867\
6847689262627473256503138475745850292904128101028969222121`94.22238864282598*^\
11"},
     {"9.345601128736616512026679216946133335882639576193476782863219399725329\
2256094287186555342045012766467821310334331787344897271`94.22238864282612*^\
10", "5.4815746237462729454170774568255878142197725147529278235886281752513508\
2821631786227179651109013053284859431126880208315944548`94.22238864282603*^\
10", "1.1100054877530853998054009228524303927676199745472657217544649154221790\
50642829266100142783163738195599123620603244168763579248`94.222388642826*^11"},
     {"1.119245551835973684052856709367887255922722530004000677320942096289488\
53008948083512475566903417218650561437141766959313139083`94.22238864282612*^\
11", "6.5653688026449125612073083525091977076824947110277121515331781170696134\
7952061200675557201496754850455942192326227552918452408`94.22238864282606*^\
10", "1.3350110591653227404893859342981156314910290887443570183272046302046558\
56544737799601948347130696309490062147281960635369270676`94.22238864282602*^\
11"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"N", "[", 
   RowBox[{"SingularValueList", "[", 
    SubscriptBox["A", "val"], "]"}], "]"}], "\[IndentingNewLine]", 
  RowBox[{"(*", " ", 
   RowBox[{"condition", " ", "number", " ", "is", " ", "very", " ", "large"}],
    " ", "*)"}]}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"Max", "[", "%", "]"}], "/", 
  RowBox[{"Min", "[", "%", "]"}]}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
  "3.0219852591785747`*^12", ",", "1.2237855757265606`*^10", ",", 
   "2.202554945363917`*^9", ",", "2.134291865167249`*^7", ",", 
   "1.1815250111148562`*^6", ",", "50255.61368382324`", ",", 
   "34207.30159771067`", ",", "4304.743791696316`", ",", 
   "1205.3281069375075`", ",", "574.3660127467067`", ",", 
   "203.20691911940733`", ",", "0.693703265742814`", ",", 
   "0.13058300313349894`", ",", "0.05602241011245954`", ",", 
   "0.027916180927763927`", ",", "0.0054449504441229`", ",", 
   "0.00045142057669584593`", ",", "0.0000612440059404503`", ",", 
   "8.274342098879295`*^-6", ",", "6.052496179507718`*^-6", ",", 
   "8.081079912124847`*^-7", ",", "4.075829477183962`*^-8", ",", 
   "7.762355499426161`*^-9", ",", "1.3077972709530828`*^-10"}], 
  "}"}]], "Output"],

Cell[BoxData["2.310744429812309`*^22"], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"after", " ", "adding", " ", "identity"}], ",", " ", 
    RowBox[{
    "smallest", " ", "singular", " ", "value", " ", "now", " ", "of", " ", 
     "order", " ", "1"}]}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"N", "[", 
    RowBox[{"SingularValueList", "[", 
     RowBox[{
      RowBox[{"IdentityMatrix", "[", 
       SubscriptBox["n", "sites"], "]"}], "+", 
      SubscriptBox["A", "val"]}], "]"}], "]"}], "\[IndentingNewLine]", 
   RowBox[{"(*", " ", 
    RowBox[{"condition", " ", "number"}], " ", "*)"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"Max", "[", "%", "]"}], "/", 
    RowBox[{"Min", "[", "%", "]"}]}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
  "3.0219852591794004`*^12", ",", "1.2237855757942867`*^10", ",", 
   "2.202554946010504`*^9", ",", "2.134291923751629`*^7", ",", 
   "1.1815255359533175`*^6", ",", "50256.11258534778`", ",", 
   "34207.62756073204`", ",", "4305.034280750286`", ",", 
   "1205.6385457832491`", ",", "574.8351775059595`", ",", 
   "203.74121107639343`", ",", "1.4622799323631732`", ",", 
   "1.0702541386926026`", ",", "1.0108416116579262`", ",", 
   "0.9913765962944453`", ",", "0.976464685623109`", ",", 
   "0.9553051610991544`", ",", "0.9514354128807613`", ",", 
   "0.9278755910350819`", ",", "0.894045284513373`", ",", 
   "0.8416572477215143`", ",", "0.7370206157114255`", ",", 
   "0.4884540946209657`", ",", "0.46563665277800215`"}], "}"}]], "Output"],

Cell[BoxData["6.490007264570231`*^12"], "Output"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Calculate Green\[CloseCurlyQuote]s function matrix", "Subsection"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["G", "val"], "=", 
   RowBox[{"Inverse", "[", 
    RowBox[{
     RowBox[{"IdentityMatrix", "[", 
      SubscriptBox["n", "sites"], "]"}], "+", 
     SubscriptBox["A", "val"]}], "]"}]}], ";"}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["G", "val"], "\[LeftDoubleBracket]", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"1", ",", "2", ",", 
       RowBox[{"-", "1"}]}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"1", ",", "2", ",", 
       RowBox[{"-", "1"}]}], "}"}]}], "\[RightDoubleBracket]"}], "//", 
   "MatrixForm"}]}]], "Input"],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"0.135228012837714371306533817161518567386381307878782428265512669592294\
0158847564181139618612935563226735045705134249`94.22238864282596", 
      RowBox[{
      "-", "0.2595300980988047268763380467432946248863589008568433420114746602\
383501885546882722334078464354384545700878010451632`94.22238864282596"}], 
      RowBox[{
      "-", "0.0075386027571279047961002006547927538214322010744746893531910064\
604378332155713694397937816273547502248210534760642`94.22238864282596"}]},
     {
      RowBox[{
      "-", "0.1825505988157490184914307351151524752916914732166034034056828702\
969584040061996533897177663603503771267433821647227`94.22238864282596"}], 
      "0.573962992593075532589664907463859458408460507635772689346034717456952\
8957933218630352254422734629285773435583501434`94.22238864282596", 
      "0.236692732071564859204586656547447250962779546763332924873267923254261\
1865109674395606642820043480869232050308485246`94.22238864282596"},
     {
      RowBox[{
      "-", "0.0063326087197388427108166044005499072391022655416723804203038264\
154402383326858332540513868328999558565017888979785`94.22238864282596"}], 
      "0.071738349652883096636998526353118202759404432860279150549755485681094\
7527765993267496572274552796526982281583066381`94.22238864282596", 
      "0.596155363575908051862747372576897010364489450095341948946604204469414\
2086832752016610243657316596677553530815395609`94.22238864282596"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"N", "[", 
  RowBox[{"SingularValueList", "[", 
   SubscriptBox["G", "val"], "]"}], "]"}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
  "2.147597260726728`", ",", "2.04727529365065`", ",", "1.3568141496757562`", 
   ",", "1.1881321080607836`", ",", "1.1185115757803006`", ",", 
   "1.0777306889649512`", ",", "1.051043493296297`", ",", 
   "1.0467859284350778`", ",", "1.0241025760822804`", ",", 
   "1.0086984136379527`", ",", "0.9892746682240906`", ",", 
   "0.9343575173851479`", ",", "0.6838635871750711`", ",", 
   "0.004908187178808154`", ",", "0.0017396290956630481`", ",", 
   "0.000829435989333225`", ",", "0.00023228618746927118`", ",", 
   "0.000029233246246750238`", ",", "0.000019898077040912053`", ",", 
   "8.463634255633306`*^-7", ",", "4.68539466823364`*^-8", ",", 
   "4.540181854765094`*^-10", ",", "8.171366126381733`*^-11", ",", 
   "3.3090829843145667`*^-13"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"largest", " ", "and", " ", "smallest", " ", "entry"}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"N", "[", 
    RowBox[{"Max", "[", 
     RowBox[{"Flatten", "[", 
      SubscriptBox["G", "val"], "]"}], "]"}], "]"}], "\[IndentingNewLine]", 
   RowBox[{"N", "[", 
    RowBox[{"Min", "[", 
     RowBox[{"Flatten", "[", 
      SubscriptBox["G", "val"], "]"}], "]"}], "]"}], "\[IndentingNewLine]", 
   RowBox[{"N", "[", 
    RowBox[{"Min", "[", 
     RowBox[{"Abs", "[", 
      RowBox[{"Flatten", "[", 
       SubscriptBox["G", "val"], "]"}], "]"}], "]"}], "]"}]}]}]], "Input"],

Cell[BoxData["0.9572111941297435`"], "Output"],

Cell[BoxData[
 RowBox[{"-", "1.2595958880218643`"}]], "Output"],

Cell[BoxData["0.00012575754287463709`"], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "determinant", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   SubscriptBox["detG", "val"], "=", 
   RowBox[{"Det", "[", 
    SubscriptBox["G", "val"], "]"}]}]}]], "Input"],

Cell[BoxData["2.\
859555739564023419808947168893190028129622025977029277045794873347031934054556\
663851419774620834688218097108346313982096765823`94.22238864282596*^-66"], \
"Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"save", " ", "as", " ", "reference", " ", "to", " ", "disk"}], " ",
    "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{"Export", "[", 
     RowBox[{
      RowBox[{
       RowBox[{"NotebookDirectory", "[", "]"}], "<>", 
       RowBox[{"FileBaseName", "[", 
        RowBox[{"NotebookFileName", "[", "]"}], "]"}], "<>", 
       "\"\<_G.dat\>\""}], ",", 
      RowBox[{"Flatten", "[", 
       RowBox[{"Transpose", "[", 
        RowBox[{"N", "[", 
         SubscriptBox["G", "val"], "]"}], "]"}], "]"}], ",", 
      "\"\<Real64\>\""}], "]"}], ";"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"Export", "[", 
     RowBox[{
      RowBox[{
       RowBox[{"NotebookDirectory", "[", "]"}], "<>", 
       RowBox[{"FileBaseName", "[", 
        RowBox[{"NotebookFileName", "[", "]"}], "]"}], "<>", 
       "\"\<_detG.dat\>\""}], ",", 
      RowBox[{"N", "[", 
       SubscriptBox["detG", "val"], "]"}], ",", "\"\<Real64\>\""}], "]"}], 
    ";"}]}]}]], "Input"]
}, Open  ]]
},
WindowSize->{1424, 867},
WindowMargins->{{247, Automatic}, {103, Automatic}},
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
Cell[1486, 35, 41, 0, 43, "Subsection"],
Cell[1530, 37, 320, 10, 72, "Input"],
Cell[CellGroupData[{
Cell[1875, 51, 947, 29, 72, "Input"],
Cell[2825, 82, 1469, 50, 52, "Output"],
Cell[4297, 134, 74, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[4408, 141, 295, 8, 52, "Input"],
Cell[4706, 151, 29, 0, 31, "Output"]
}, Open  ]],
Cell[4750, 154, 1122, 32, 52, "Input"],
Cell[CellGroupData[{
Cell[5897, 190, 336, 10, 72, "Input"],
Cell[6236, 202, 28, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[6301, 207, 358, 10, 52, "Input"],
Cell[6662, 219, 28, 0, 31, "Output"]
}, Open  ]],
Cell[6705, 222, 1167, 33, 59, "Input"],
Cell[CellGroupData[{
Cell[7897, 259, 340, 10, 72, "Input"],
Cell[8240, 271, 28, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[8305, 276, 426, 13, 52, "Input"],
Cell[8734, 291, 28, 0, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[8811, 297, 45, 0, 43, "Subsection"],
Cell[8859, 299, 228, 7, 67, "Input"],
Cell[9090, 308, 221, 7, 52, "Input"],
Cell[9314, 317, 309, 10, 67, "Input"],
Cell[9626, 329, 258, 8, 67, "Input"],
Cell[9887, 339, 486, 16, 31, "Input"],
Cell[10376, 357, 513, 16, 52, "Input"],
Cell[CellGroupData[{
Cell[10914, 377, 480, 14, 52, "Input"],
Cell[11397, 393, 789, 18, 71, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[12235, 417, 48, 0, 43, "Subsection"],
Cell[CellGroupData[{
Cell[12308, 421, 823, 26, 90, "Input"],
Cell[13134, 449, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[13246, 456, 567, 15, 92, "Input"],
Cell[13816, 473, 521, 14, 31, "Output"],
Cell[14340, 489, 521, 14, 31, "Output"],
Cell[14864, 505, 559, 16, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[15460, 526, 302, 10, 31, "Input"],
Cell[15765, 538, 8144, 105, 292, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[23958, 649, 45, 0, 43, "Subsection"],
Cell[24006, 651, 733, 19, 52, "Input"],
Cell[24742, 672, 124, 4, 31, "Input"],
Cell[CellGroupData[{
Cell[24891, 680, 347, 10, 52, "Input"],
Cell[25241, 692, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[25353, 699, 436, 13, 52, "Input"],
Cell[25792, 714, 1988, 35, 80, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[27817, 754, 392, 10, 72, "Input"],
Cell[28212, 766, 816, 15, 55, "Output"],
Cell[29031, 783, 49, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[29117, 788, 729, 19, 92, "Input"],
Cell[29849, 809, 782, 14, 55, "Output"],
Cell[30634, 825, 49, 0, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[30732, 831, 72, 0, 43, "Subsection"],
Cell[30807, 833, 253, 8, 31, "Input"],
Cell[CellGroupData[{
Cell[31085, 845, 436, 13, 52, "Input"],
Cell[31524, 860, 1985, 39, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[33546, 904, 123, 3, 31, "Input"],
Cell[33672, 909, 807, 14, 55, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[34516, 928, 651, 18, 92, "Input"],
Cell[35170, 948, 46, 0, 31, "Output"],
Cell[35219, 950, 63, 1, 31, "Output"],
Cell[35285, 953, 50, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[35372, 958, 220, 6, 52, "Input"],
Cell[35595, 966, 184, 3, 31, "Output"]
}, Open  ]],
Cell[35794, 972, 1030, 29, 72, "Input"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature qv06hrFYXGp0FAKTUWdFJ8YV *)
