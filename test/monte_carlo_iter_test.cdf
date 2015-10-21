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
NotebookDataLength[     55851,       1625]
NotebookOptionsPosition[     52786,       1498]
NotebookOutlinePosition[     53129,       1513]
CellTagsIndexPosition[     53086,       1510]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["General parameters", "Subsection"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"Coulomb", " ", "coupling", " ", "constant"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["U", "val"], "=", 
     RowBox[{"9", "/", "2"}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"N", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData["4.5`"], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"inverse", " ", "temperature"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[Beta]", "val"], "=", "2"}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"imaginary", "-", 
    RowBox[{"time", " ", "step"}]}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[Tau]", "val"], "=", 
    FractionBox["1", "8"]}], ";"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{
   SubscriptBox["\[Lambda]", "val"], "=", 
   RowBox[{"ArcCosh", "[", 
    RowBox[{"Exp", "[", 
     RowBox[{
      SubscriptBox["\[Tau]", "val"], 
      RowBox[{
       SubscriptBox["U", "val"], "/", "2"}]}], "]"}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"N", "[", "%", "]"}]}], "Input"],

Cell[BoxData["0.7856003600934582`"], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"number", " ", "of", " ", "time", " ", "steps"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   SubscriptBox["L", "val"], "=", 
   FractionBox[
    SubscriptBox["\[Beta]", "val"], 
    SubscriptBox["\[Tau]", "val"]]}]}]], "Input"],

Cell[BoxData["16"], "Output"]
}, Open  ]]
}, Open  ]],

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
   RowBox[{"direct", " ", "neighbors"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["latt", "neigh"], "=", 
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
    SubscriptBox["latt", "neigh"], ";"}], "\[IndentingNewLine]", 
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
      SubscriptBox["latt", "neigh"]}], ")"}], "-", "4"}], "]"}]}]], "Input"],

Cell[BoxData["0"], "Output"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Kinetic energy operator", "Subsection"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"symbolic", " ", 
    SuperscriptBox["\[ExponentialE]", 
     RowBox[{
      RowBox[{"-", "\[Tau]"}], " ", "k"}]], " ", "and", " ", 
    SuperscriptBox["\[ExponentialE]", 
     RowBox[{"\[Tau]", " ", "k"}]]}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["exp", "\[Tau]k"], "=", 
     RowBox[{"FullSimplify", "[", 
      RowBox[{"MatrixExp", "[", 
       RowBox[{
        RowBox[{"-", 
         SubscriptBox["\[Tau]", "val"]}], " ", 
        RowBox[{"(", 
         RowBox[{"-", 
          SubscriptBox["latt", "neigh"]}], ")"}]}], "]"}], "]"}]}], ";"}], 
   "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{
     SubscriptBox["invexp", "\[Tau]k"], "=", 
     RowBox[{"FullSimplify", "[", 
      RowBox[{"MatrixExp", "[", 
       RowBox[{
        SubscriptBox["\[Tau]", "val"], " ", 
        RowBox[{"(", 
         RowBox[{"-", 
          SubscriptBox["latt", "neigh"]}], ")"}]}], "]"}], "]"}]}], 
    ";"}]}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["exp", "\[Tau]k"], "\[LeftDoubleBracket]", 
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
     {
      RowBox[{
       FractionBox["1", "12"], " ", 
       RowBox[{"(", 
        RowBox[{"1", "+", 
         RowBox[{"6", " ", 
          RowBox[{"Cosh", "[", 
           FractionBox["1", "8"], "]"}]}], "+", 
         RowBox[{"2", " ", 
          RowBox[{"Cosh", "[", 
           FractionBox["1", "4"], "]"}]}], "+", 
         RowBox[{"2", " ", 
          RowBox[{"Cosh", "[", 
           FractionBox["3", "8"], "]"}]}], "+", 
         RowBox[{"Cosh", "[", 
          FractionBox["1", "2"], "]"}]}], ")"}]}], 
      RowBox[{
       FractionBox["1", "12"], " ", 
       RowBox[{"(", 
        RowBox[{
         RowBox[{"2", " ", 
          RowBox[{"(", 
           RowBox[{
            RowBox[{"Sinh", "[", 
             FractionBox["1", "8"], "]"}], "+", 
            RowBox[{"Sinh", "[", 
             FractionBox["3", "8"], "]"}]}], ")"}]}], "+", 
         RowBox[{"Sinh", "[", 
          FractionBox["1", "2"], "]"}]}], ")"}]}], 
      RowBox[{
       FractionBox["1", "12"], " ", 
       RowBox[{"(", 
        RowBox[{
         RowBox[{"-", "1"}], "-", 
         RowBox[{"Cosh", "[", 
          FractionBox["1", "8"], "]"}], "+", 
         RowBox[{"Cosh", "[", 
          FractionBox["3", "8"], "]"}], "+", 
         RowBox[{"Cosh", "[", 
          FractionBox["1", "2"], "]"}]}], ")"}]}]},
     {
      RowBox[{
       FractionBox["1", "12"], " ", 
       RowBox[{"(", 
        RowBox[{
         RowBox[{"2", " ", 
          RowBox[{"(", 
           RowBox[{
            RowBox[{"Sinh", "[", 
             FractionBox["1", "8"], "]"}], "+", 
            RowBox[{"Sinh", "[", 
             FractionBox["3", "8"], "]"}]}], ")"}]}], "+", 
         RowBox[{"Sinh", "[", 
          FractionBox["1", "2"], "]"}]}], ")"}]}], 
      RowBox[{
       FractionBox["1", "12"], " ", 
       RowBox[{"(", 
        RowBox[{"1", "+", 
         RowBox[{"6", " ", 
          RowBox[{"Cosh", "[", 
           FractionBox["1", "8"], "]"}]}], "+", 
         RowBox[{"2", " ", 
          RowBox[{"Cosh", "[", 
           FractionBox["1", "4"], "]"}]}], "+", 
         RowBox[{"2", " ", 
          RowBox[{"Cosh", "[", 
           FractionBox["3", "8"], "]"}]}], "+", 
         RowBox[{"Cosh", "[", 
          FractionBox["1", "2"], "]"}]}], ")"}]}], 
      FractionBox[
       RowBox[{
        SuperscriptBox[
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "+", 
           SuperscriptBox["\[ExponentialE]", 
            RowBox[{"1", "/", "4"}]]}], ")"}], "3"], " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", 
          SuperscriptBox["\[ExponentialE]", 
           RowBox[{"1", "/", "8"}]], "+", 
          SuperscriptBox["\[ExponentialE]", 
           RowBox[{"1", "/", "4"}]]}], ")"}]}], 
       RowBox[{"24", " ", 
        SqrtBox["\[ExponentialE]"]}]]},
     {
      RowBox[{
       FractionBox["1", "12"], " ", 
       RowBox[{"(", 
        RowBox[{
         RowBox[{"-", "1"}], "-", 
         RowBox[{"Cosh", "[", 
          FractionBox["1", "8"], "]"}], "+", 
         RowBox[{"Cosh", "[", 
          FractionBox["3", "8"], "]"}], "+", 
         RowBox[{"Cosh", "[", 
          FractionBox["1", "2"], "]"}]}], ")"}]}], 
      FractionBox[
       RowBox[{
        SuperscriptBox[
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "+", 
           SuperscriptBox["\[ExponentialE]", 
            RowBox[{"1", "/", "4"}]]}], ")"}], "3"], " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", 
          SuperscriptBox["\[ExponentialE]", 
           RowBox[{"1", "/", "8"}]], "+", 
          SuperscriptBox["\[ExponentialE]", 
           RowBox[{"1", "/", "4"}]]}], ")"}]}], 
       RowBox[{"24", " ", 
        SqrtBox["\[ExponentialE]"]}]], 
      RowBox[{
       FractionBox["1", "12"], " ", 
       RowBox[{"(", 
        RowBox[{"1", "+", 
         RowBox[{"6", " ", 
          RowBox[{"Cosh", "[", 
           FractionBox["1", "8"], "]"}]}], "+", 
         RowBox[{"2", " ", 
          RowBox[{"Cosh", "[", 
           FractionBox["1", "4"], "]"}]}], "+", 
         RowBox[{"2", " ", 
          RowBox[{"Cosh", "[", 
           FractionBox["3", "8"], "]"}]}], "+", 
         RowBox[{"Cosh", "[", 
          FractionBox["1", "2"], "]"}]}], ")"}]}]}
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
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"numerical", " ", "values"}], " ", "*)"}], "\[IndentingNewLine]", 
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
     {"1.0316390771107469`", "0.1282874159836006`", "0.01591196950710468`"},
     {"0.1282874159836006`", "1.0316390771107469`", "0.0019787011723065764`"},
     {"0.01591196950710468`", "0.0019787011723065764`", "1.0316390771107469`"}
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
 RowBox[{
  RowBox[{"(*", " ", "check", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{"Norm", "[", 
   RowBox[{
    RowBox[{"FullSimplify", "[", 
     RowBox[{
      SubscriptBox["exp", "\[Tau]k"], ".", 
      SubscriptBox["invexp", "\[Tau]k"]}], "]"}], "-", 
    RowBox[{"IdentityMatrix", "[", 
     SubscriptBox["n", "sites"], "]"}]}], "]"}]}]], "Input"],

Cell[BoxData["0"], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"compare", " ", "with", " ", "numerical", " ", "evaluation"}], " ",
    "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"Norm", "[", 
    RowBox[{
     SubscriptBox["exp", "\[Tau]k"], "-", 
     RowBox[{"MatrixExp", "[", 
      RowBox[{"N", "[", 
       RowBox[{
        RowBox[{"-", 
         SubscriptBox["\[Tau]", "val"]}], 
        RowBox[{"(", 
         RowBox[{"-", 
          SubscriptBox["latt", "neigh"]}], ")"}]}], "]"}], "]"}]}], "]"}], 
   "\[IndentingNewLine]", 
   RowBox[{"Norm", "[", 
    RowBox[{
     SubscriptBox["invexp", "\[Tau]k"], "-", 
     RowBox[{"MatrixExp", "[", 
      RowBox[{"N", "[", 
       RowBox[{
        SubscriptBox["\[Tau]", "val"], 
        RowBox[{"(", 
         RowBox[{"-", 
          SubscriptBox["latt", "neigh"]}], ")"}]}], "]"}], "]"}]}], 
    "]"}]}]}]], "Input"],

Cell[BoxData["9.59646170837076`*^-16"], "Output"],

Cell[BoxData["9.599818764028589`*^-16"], "Output"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Initial Hubbard-Stratonovich field", "Subsection"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"SeedRandom", "[", "42", "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   SubscriptBox["s", "0"], "=", 
   RowBox[{
    RowBox[{"2", 
     RowBox[{"RandomInteger", "[", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{"0", ",", "1"}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         SubscriptBox["n", "sites"], ",", 
         SubscriptBox["L", "val"]}], "}"}]}], "]"}]}], "-", "1"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"Dimensions", "[", 
  SubscriptBox["s", "0"], "]"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"24", ",", "16"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   SubscriptBox["s", "0"], "\[LeftDoubleBracket]", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{"1", ",", "2", ",", 
      RowBox[{"-", "1"}]}], "}"}], ",", 
    RowBox[{"{", 
     RowBox[{"1", ",", "2", ",", 
      RowBox[{"-", "1"}]}], "}"}]}], "\[RightDoubleBracket]"}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "1"}], ",", "1", ",", 
     RowBox[{"-", "1"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "1"}], ",", 
     RowBox[{"-", "1"}], ",", 
     RowBox[{"-", "1"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "1"}], ",", "1", ",", 
     RowBox[{"-", "1"}]}], "}"}]}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Flatten", "[", 
  RowBox[{
   RowBox[{
    RowBox[{"Transpose", "[", 
     SubscriptBox["s", "0"], "]"}], "/.", 
    RowBox[{"{", 
     RowBox[{"1", "\[Rule]", "0"}], "}"}]}], "/.", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "1"}], "\[Rule]", "1"}], "}"}]}], "]"}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
  "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", 
   ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", ",", 
   "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", 
   ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", 
   "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", 
   ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", ",", "0", ",", "0", ",", 
   "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", 
   ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", 
   ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", 
   ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", 
   "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", 
   ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", 
   "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", 
   ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", 
   "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", 
   ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", 
   "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", 
   ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", 
   "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", 
   ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", 
   ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", 
   ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", 
   "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", 
   ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", 
   "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", 
   ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", 
   "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", 
   ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", 
   "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", 
   ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", 
   ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", 
   "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", 
   ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", 
   "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", 
   ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", 
   ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", 
   ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", 
   "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", 
   ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", 
   ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", 
   "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", 
   ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", 
   "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", 
   ",", "1"}], "}"}]], "Output"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["DQMC step", "Subsection"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
   "time", " ", "flow", " ", "map", " ", "generated", " ", "by", " ", "the", 
    " ", "Hubbard", " ", "Hamiltonian"}], " ", "*)"}], "\n", 
  RowBox[{
   RowBox[{"HubbardTimeFlowMap", "[", 
    RowBox[{"expK_", ",", "expV_"}], "]"}], ":=", 
   RowBox[{"Fold", "[", 
    RowBox[{
     RowBox[{
      RowBox[{
       RowBox[{"DiagonalMatrix", "[", "#2", "]"}], ".", "expK", ".", "#1"}], 
      "&"}], ",", 
     RowBox[{"IdentityMatrix", "[", 
      RowBox[{"Length", "[", "expK", "]"}], "]"}], ",", 
     RowBox[{"Transpose", "[", "expV", "]"}]}], "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{
     RowBox[{"initialize", " ", "spin"}], "-", 
     RowBox[{"up", " ", "or", " ", "spin"}], "-", 
     RowBox[{"down", " ", 
      RowBox[{"Green", "'"}], "s", " ", "function"}]}], ",", " ", 
    RowBox[{"depending", " ", "on", " ", "the", " ", 
     RowBox[{"(", 
      RowBox[{"\[PlusMinus]", "1"}], ")"}], " ", "prefactor", " ", 
     RowBox[{"of", " ", "'"}], 
     RowBox[{"\[Lambda]", "'"}]}]}], " ", "*)"}], "\n", 
  RowBox[{
   RowBox[{"InitializeGreensFunction", "[", 
    RowBox[{"expK_", ",", "expV_"}], "]"}], ":=", 
   RowBox[{"Inverse", "[", 
    RowBox[{
     RowBox[{"IdentityMatrix", "[", 
      RowBox[{"Length", "[", "expK", "]"}], "]"}], "+", 
     RowBox[{"HubbardTimeFlowMap", "[", 
      RowBox[{"expK", ",", "expV"}], "]"}]}], "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"InitializeUpDownGreensFunctions", "[", 
   RowBox[{"expK_", ",", "\[Lambda]s_"}], "]"}], ":=", 
  RowBox[{"{", "\n", "\t", 
   RowBox[{
    RowBox[{"InitializeGreensFunction", "[", 
     RowBox[{"expK", ",", 
      RowBox[{"Exp", "[", 
       RowBox[{"-", "\[Lambda]s"}], "]"}]}], "]"}], ",", "\n", "\t", 
    RowBox[{"InitializeGreensFunction", "[", 
     RowBox[{"expK", ",", 
      RowBox[{"Exp", "[", 
       RowBox[{"+", "\[Lambda]s"}], "]"}]}], "]"}]}], "}"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"HubbardDQMCStep", "[", 
   RowBox[{
   "expK_", ",", "invexpK_", ",", "\[Lambda]_", ",", "s0_", ",", 
    "siteorder_List", ",", "rnd_List"}], "]"}], ":=", 
  RowBox[{"Module", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"s", "=", "s0"}], ",", "Gu", ",", "Gd", ",", "nsites", ",", "L",
       ",", "du", ",", "dd", ",", "cu", ",", "cd", ",", "l", ",", "i", ",", 
      "j"}], "}"}], ",", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{
      RowBox[{"{", 
       RowBox[{"nsites", ",", "L"}], "}"}], "=", 
      RowBox[{"Dimensions", "[", "s", "]"}]}], ";", "\[IndentingNewLine]", 
     RowBox[{"(*", " ", 
      RowBox[{"compute", " ", "the", " ", "initial", " ", 
       RowBox[{"Green", "'"}], "s", " ", "functions"}], " ", "*)"}], 
     "\[IndentingNewLine]", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{"Gu", ",", "Gd"}], "}"}], "=", 
      RowBox[{"InitializeUpDownGreensFunctions", "[", 
       RowBox[{"expK", ",", 
        RowBox[{"\[Lambda]", " ", "s"}]}], "]"}]}], ";", 
     "\[IndentingNewLine]", 
     RowBox[{"(*", " ", 
      RowBox[{"iterate", " ", "over", " ", "time", " ", "slices"}], " ", 
      "*)"}], "\[IndentingNewLine]", 
     RowBox[{"For", "[", 
      RowBox[{
       RowBox[{"l", "=", "1"}], ",", 
       RowBox[{"l", "\[LessEqual]", "L"}], ",", 
       RowBox[{"l", "++"}], ",", "\[IndentingNewLine]", 
       RowBox[{"(*", " ", 
        RowBox[{"Eq", " ", 
         RowBox[{"(", "16", ")"}]}], " ", "*)"}], "\[IndentingNewLine]", 
       RowBox[{
        RowBox[{"Gu", "=", 
         RowBox[{
          RowBox[{"DiagonalMatrix", "[", 
           RowBox[{"Exp", "[", 
            RowBox[{
             RowBox[{"-", "\[Lambda]"}], " ", 
             RowBox[{"s", "\[LeftDoubleBracket]", 
              RowBox[{";;", ",", "l"}], "\[RightDoubleBracket]"}]}], "]"}], 
           "]"}], ".", "expK", ".", "Gu", ".", "invexpK", ".", 
          RowBox[{"DiagonalMatrix", "[", 
           RowBox[{"Exp", "[", 
            RowBox[{
             RowBox[{"+", "\[Lambda]"}], " ", 
             RowBox[{"s", "\[LeftDoubleBracket]", 
              RowBox[{";;", ",", "l"}], "\[RightDoubleBracket]"}]}], "]"}], 
           "]"}]}]}], ";", "\[IndentingNewLine]", 
        RowBox[{"Gd", "=", 
         RowBox[{
          RowBox[{"DiagonalMatrix", "[", 
           RowBox[{"Exp", "[", 
            RowBox[{
             RowBox[{"+", "\[Lambda]"}], " ", 
             RowBox[{"s", "\[LeftDoubleBracket]", 
              RowBox[{";;", ",", "l"}], "\[RightDoubleBracket]"}]}], "]"}], 
           "]"}], ".", "expK", ".", "Gd", ".", "invexpK", ".", 
          RowBox[{"DiagonalMatrix", "[", 
           RowBox[{"Exp", "[", 
            RowBox[{
             RowBox[{"-", "\[Lambda]"}], " ", 
             RowBox[{"s", "\[LeftDoubleBracket]", 
              RowBox[{";;", ",", "l"}], "\[RightDoubleBracket]"}]}], "]"}], 
           "]"}]}]}], ";", "\[IndentingNewLine]", 
        RowBox[{"(*", " ", 
         RowBox[{"iterate", " ", "over", " ", "lattice", " ", "sites"}], " ", 
         "*)"}], "\[IndentingNewLine]", 
        RowBox[{"For", "[", 
         RowBox[{
          RowBox[{"j", "=", "1"}], ",", 
          RowBox[{"j", "\[LessEqual]", "nsites"}], ",", 
          RowBox[{"j", "++"}], ",", "\[IndentingNewLine]", 
          RowBox[{
           RowBox[{"i", "=", 
            RowBox[{"siteorder", "\[LeftDoubleBracket]", 
             RowBox[{"l", ",", "j"}], "\[RightDoubleBracket]"}]}], ";", 
           "\[IndentingNewLine]", 
           RowBox[{"(*", " ", 
            RowBox[{"Eq", " ", 
             RowBox[{"(", "13", ")"}]}], " ", "*)"}], "\[IndentingNewLine]", 
           RowBox[{"(*", " ", 
            RowBox[{"suggest", " ", "flipping", " ", "s", 
             RowBox[{"(", 
              RowBox[{"i", ",", "l"}], ")"}]}], " ", "*)"}], 
           "\[IndentingNewLine]", 
           RowBox[{"du", "=", 
            RowBox[{"1", "+", 
             RowBox[{
              RowBox[{"(", 
               RowBox[{"1", "-", 
                RowBox[{"Gu", "\[LeftDoubleBracket]", 
                 RowBox[{"i", ",", "i"}], "\[RightDoubleBracket]"}]}], ")"}], 
              RowBox[{"(", 
               RowBox[{
                RowBox[{"Exp", "[", 
                 RowBox[{
                  RowBox[{"+", "2"}], "\[Lambda]", " ", 
                  RowBox[{"s", "\[LeftDoubleBracket]", 
                   RowBox[{"i", ",", "l"}], "\[RightDoubleBracket]"}]}], 
                 "]"}], "-", "1"}], ")"}]}]}]}], ";", "\[IndentingNewLine]", 
           RowBox[{"dd", "=", 
            RowBox[{"1", "+", 
             RowBox[{
              RowBox[{"(", 
               RowBox[{"1", "-", 
                RowBox[{"Gd", "\[LeftDoubleBracket]", 
                 RowBox[{"i", ",", "i"}], "\[RightDoubleBracket]"}]}], ")"}], 
              RowBox[{"(", 
               RowBox[{
                RowBox[{"Exp", "[", 
                 RowBox[{
                  RowBox[{"-", "2"}], "\[Lambda]", " ", 
                  RowBox[{"s", "\[LeftDoubleBracket]", 
                   RowBox[{"i", ",", "l"}], "\[RightDoubleBracket]"}]}], 
                 "]"}], "-", "1"}], ")"}]}]}]}], ";", "\[IndentingNewLine]", 
           RowBox[{"If", "[", 
            RowBox[{"(*", 
             RowBox[{"RandomReal", "[", "]"}], "*)"}], 
            RowBox[{
             RowBox[{
              RowBox[{"rnd", "\[LeftDoubleBracket]", 
               RowBox[{"l", ",", "j"}], "\[RightDoubleBracket]"}], "<", 
              RowBox[{"du", " ", "dd"}]}], ",", "\[IndentingNewLine]", 
             RowBox[{"(*", " ", 
              RowBox[{"Eq", " ", 
               RowBox[{"(", "15", ")"}]}], " ", "*)"}], "\[IndentingNewLine]", 
             RowBox[{
              RowBox[{"cu", "=", 
               RowBox[{
                RowBox[{"(", 
                 RowBox[{
                  RowBox[{"Exp", "[", 
                   RowBox[{
                    RowBox[{"+", "2"}], "\[Lambda]", " ", 
                    RowBox[{"s", "\[LeftDoubleBracket]", 
                    RowBox[{"i", ",", "l"}], "\[RightDoubleBracket]"}]}], 
                   "]"}], "-", "1"}], ")"}], 
                RowBox[{"(", 
                 RowBox[{
                  RowBox[{"UnitVector", "[", 
                   RowBox[{"nsites", ",", "i"}], "]"}], "-", 
                  RowBox[{
                  "Gu", "\[LeftDoubleBracket]", "i", 
                   "\[RightDoubleBracket]"}]}], ")"}]}]}], ";", 
              "\[IndentingNewLine]", 
              RowBox[{"cd", "=", 
               RowBox[{
                RowBox[{"(", 
                 RowBox[{
                  RowBox[{"Exp", "[", 
                   RowBox[{
                    RowBox[{"-", "2"}], "\[Lambda]", " ", 
                    RowBox[{"s", "\[LeftDoubleBracket]", 
                    RowBox[{"i", ",", "l"}], "\[RightDoubleBracket]"}]}], 
                   "]"}], "-", "1"}], ")"}], 
                RowBox[{"(", 
                 RowBox[{
                  RowBox[{"UnitVector", "[", 
                   RowBox[{"nsites", ",", "i"}], "]"}], "-", 
                  RowBox[{
                  "Gd", "\[LeftDoubleBracket]", "i", 
                   "\[RightDoubleBracket]"}]}], ")"}]}]}], ";", 
              "\[IndentingNewLine]", 
              RowBox[{"Gu", "-=", 
               RowBox[{"KroneckerProduct", "[", 
                RowBox[{
                 RowBox[{
                  RowBox[{"Gu", "\[LeftDoubleBracket]", 
                   RowBox[{";;", ",", "i"}], "\[RightDoubleBracket]"}], "/", 
                  RowBox[{"(", 
                   RowBox[{"1", "+", 
                    RowBox[{
                    "cu", "\[LeftDoubleBracket]", "i", 
                    "\[RightDoubleBracket]"}]}], ")"}]}], ",", "cu"}], 
                "]"}]}], ";", "\[IndentingNewLine]", 
              RowBox[{"Gd", "-=", 
               RowBox[{"KroneckerProduct", "[", 
                RowBox[{
                 RowBox[{
                  RowBox[{"Gd", "\[LeftDoubleBracket]", 
                   RowBox[{";;", ",", "i"}], "\[RightDoubleBracket]"}], "/", 
                  RowBox[{"(", 
                   RowBox[{"1", "+", 
                    RowBox[{
                    "cd", "\[LeftDoubleBracket]", "i", 
                    "\[RightDoubleBracket]"}]}], ")"}]}], ",", "cd"}], 
                "]"}]}], ";", "\[IndentingNewLine]", 
              RowBox[{"(*", " ", 
               RowBox[{"actually", " ", "flip", " ", "spin"}], " ", "*)"}], 
              "\[IndentingNewLine]", 
              RowBox[{
               RowBox[{"s", "\[LeftDoubleBracket]", 
                RowBox[{"i", ",", "l"}], "\[RightDoubleBracket]"}], "=", 
               RowBox[{"-", 
                RowBox[{"s", "\[LeftDoubleBracket]", 
                 RowBox[{"i", ",", "l"}], "\[RightDoubleBracket]"}]}]}]}]}], 
            "]"}]}]}], "]"}]}]}], "]"}], ";", "\[IndentingNewLine]", 
     RowBox[{"(*", " ", 
      RowBox[{
       RowBox[{"return", " ", "updated", " ", "Hubbard"}], "-", 
       RowBox[{"Stratonovich", " ", "field", " ", "and", " ", "spin"}], "-", 
       RowBox[{"up", " ", "and", " ", "spin"}], "-", 
       RowBox[{"down", " ", 
        RowBox[{"Green", "'"}], "s", " ", "functions"}]}], " ", "*)"}], 
     "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"s", ",", "Gu", ",", "Gd"}], "}"}]}]}], "]"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["fn", "base"], "=", 
   RowBox[{
    RowBox[{"NotebookDirectory", "[", "]"}], "<>", 
    RowBox[{"FileBaseName", "[", 
     RowBox[{"NotebookFileName", "[", "]"}], "]"}]}]}], ";"}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"pre", "-", 
     RowBox[{"computed", " ", "random", " ", "numbers"}]}], ",", " ", 
    RowBox[{
    "to", " ", "enable", " ", "same", " ", "execution", " ", "path", " ", 
     "as", " ", "C", " ", "implementation"}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["rnd", "list"], "=", 
     RowBox[{"Partition", "[", 
      RowBox[{
       RowBox[{"Import", "[", 
        RowBox[{
         RowBox[{
          SubscriptBox["fn", "base"], "<>", "\"\<_rnd.dat\>\""}], ",", 
         "\"\<Real64\>\""}], "]"}], ",", 
       SubscriptBox["n", "sites"]}], "]"}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"Dimensions", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"16", ",", "24"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   SubscriptBox["rnd", "list"], "\[LeftDoubleBracket]", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{"1", ",", "2", ",", 
      RowBox[{"-", "1"}]}], "}"}], ",", 
    RowBox[{"{", 
     RowBox[{"1", ",", "2", ",", 
      RowBox[{"-", "1"}]}], "}"}]}], "\[RightDoubleBracket]"}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{
    "0.8915577426527812`", ",", "0.27670149969199476`", ",", 
     "0.31599748821824397`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
    "0.38364571450571405`", ",", "0.5322988842984337`", ",", 
     "0.16584223657106564`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
    "0.8150237712183559`", ",", "0.4622413020415062`", ",", 
     "0.04117686139563334`"}], "}"}]}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"pre", "-", 
    RowBox[{"determined", " ", "site", " ", "update", " ", "ordering"}]}], 
   " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["site", "order"], "=", 
     RowBox[{"Partition", "[", 
      RowBox[{
       RowBox[{
        RowBox[{"Import", "[", 
         RowBox[{
          RowBox[{
           SubscriptBox["fn", "base"], "<>", "\"\<_siteorder.dat\>\""}], ",", 
          "\"\<Integer32\>\""}], "]"}], "+", "1"}], ",", 
       SubscriptBox["n", "sites"]}], "]"}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"Dimensions", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"16", ",", "24"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "check", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{"Total", "[", 
   RowBox[{
    RowBox[{
     RowBox[{"Norm", "[", 
      RowBox[{
       RowBox[{"Sort", "[", "#", "]"}], "-", 
       RowBox[{"Range", "[", 
        SubscriptBox["n", "sites"], "]"}]}], "]"}], "&"}], "/@", 
    SubscriptBox["site", "order"]}], "]"}]}]], "Input"],

Cell[BoxData["0"], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"{", 
    RowBox[{
     SubscriptBox["s", "1"], ",", 
     SubscriptBox["G", 
      RowBox[{"u", ",", "1"}]], ",", 
     SubscriptBox["G", 
      RowBox[{"d", ",", "1"}]]}], "}"}], "=", 
   RowBox[{"Block", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"$MinPrecision", "=", 
       RowBox[{"4", "MachinePrecision"}]}], "}"}], ",", 
     RowBox[{"HubbardDQMCStep", "[", 
      RowBox[{
       SubscriptBox["exp", "\[Tau]k"], ",", 
       SubscriptBox["invexp", "\[Tau]k"], ",", 
       RowBox[{"N", "[", 
        RowBox[{
         SubscriptBox["\[Lambda]", "val"], ",", 
         RowBox[{"4", "MachinePrecision"}]}], "]"}], ",", 
       SubscriptBox["s", "0"], ",", 
       SubscriptBox["site", "order"], ",", 
       SubscriptBox["rnd", "list"]}], "]"}]}], "]"}]}], ";"}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"compare", " ", "updated", " ", 
    RowBox[{"Green", "'"}], "s", " ", "functions", " ", "with", " ", 
    "recomputed", " ", 
    RowBox[{"Green", "'"}], "s", " ", "functions"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{"Norm", "[", 
   RowBox[{"Flatten", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
       SubscriptBox["G", 
        RowBox[{"u", ",", "1"}]], ",", 
       SubscriptBox["G", 
        RowBox[{"d", ",", "1"}]]}], "}"}], "-", 
     RowBox[{"InitializeUpDownGreensFunctions", "[", 
      RowBox[{
       SubscriptBox["exp", "\[Tau]k"], ",", 
       RowBox[{
        RowBox[{"N", "[", 
         RowBox[{
          SubscriptBox["\[Lambda]", "val"], ",", 
          RowBox[{"4", "MachinePrecision"}]}], "]"}], 
        SubscriptBox["s", "1"]}]}], "]"}]}], "]"}], "]"}]}]], "Input"],

Cell[BoxData["0"], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   SubscriptBox["s", "1"], "\[LeftDoubleBracket]", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{"1", ",", "2", ",", 
      RowBox[{"-", "1"}]}], "}"}], ",", 
    RowBox[{"{", 
     RowBox[{"1", ",", "2", ",", 
      RowBox[{"-", "1"}]}], "}"}]}], "\[RightDoubleBracket]"}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"1", ",", 
     RowBox[{"-", "1"}], ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "1"}], ",", "1", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", 
     RowBox[{"-", "1"}], ",", "1"}], "}"}]}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Flatten", "[", 
  RowBox[{
   RowBox[{
    RowBox[{"Transpose", "[", 
     SubscriptBox["s", "1"], "]"}], "/.", 
    RowBox[{"{", 
     RowBox[{"1", "\[Rule]", "0"}], "}"}]}], "/.", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "1"}], "\[Rule]", "1"}], "}"}]}], "]"}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
  "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", 
   ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", 
   ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", 
   ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", 
   "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", 
   ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", 
   "0", ",", "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", 
   ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", 
   ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", 
   ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", 
   "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", 
   ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", 
   "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", 
   ",", "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", 
   "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", 
   ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", 
   "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "1", 
   ",", "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", 
   "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", 
   ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", 
   "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", 
   ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", 
   "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", 
   ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", 
   "1", ",", "0", ",", "0", ",", "0", ",", "0", ",", "0", ",", "0", ",", "0", 
   ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", 
   "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", 
   ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", 
   "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", 
   ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", 
   "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", 
   ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", 
   "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", 
   ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", 
   "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", 
   ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", 
   "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", 
   ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", ",", "0", ",", 
   "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", 
   ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", 
   "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", 
   ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", 
   "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", 
   ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", ",", 
   "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", 
   ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", 
   "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", 
   ",", "0"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["G", 
      RowBox[{"u", ",", "1"}]], "\[LeftDoubleBracket]", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{"1", ",", "2", ",", 
        RowBox[{"-", "1"}]}], "}"}], ",", 
      RowBox[{"{", 
       RowBox[{"1", ",", "2", ",", 
        RowBox[{"-", "1"}]}], "}"}]}], "\[RightDoubleBracket]"}], "//", 
    "MatrixForm"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{
     SubscriptBox["G", 
      RowBox[{"d", ",", "1"}]], "\[LeftDoubleBracket]", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{"1", ",", "2", ",", 
        RowBox[{"-", "1"}]}], "}"}], ",", 
      RowBox[{"{", 
       RowBox[{"1", ",", "2", ",", 
        RowBox[{"-", "1"}]}], "}"}]}], "\[RightDoubleBracket]"}], "//", 
    "MatrixForm"}]}]}]], "Input"],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"0.869584831000293924513350315488674857968384989664025812192839259500673\
964969415265121627889207385`63.81835908076401", 
      RowBox[{
      "-", "0.1036658044602740983010797465641109867859450376677572439433027721\
12104492174849004087566174830865`63.81835908076401"}], 
      "0.002485330263204113781556631390688992937907877833173897712881760725932\
52182195358359463530493801`63.81835908076401"},
     {
      RowBox[{
      "-", "0.2456743536496151757757937224889116846036764401207361412666485562\
56331124549670997994652192968953`63.81835908076401"}], 
      "0.679520191197914170912088900965594788885333232234209774471689257439710\
707715361330771898472477614`63.81835908076401", 
      RowBox[{
      "-", "0.0927277433021352928671983882892348540926135026458856857786324395\
59447508196630631779379799081687`63.81835908076401"}]},
     {
      RowBox[{
      "-", "0.0249280638546788767367264242075386836289428002504972359645685468\
59687769987191992217466991504545`63.81835908076402"}], 
      "0.019622135746464034092096399628813170705319259210498083208079764352159\
297260834697108992079805663`63.81835908076402", 
      "0.947848280739576672484332966768146623669428563055222734173331881471682\
579460827376530882188639542`63.81835908076401"}
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
   MatrixForm[BoxForm`e$]]]], "Output"],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"0.130415168999706075486649684511325142031615010335974187807160740499325\
126748546308172594160200565`63.81835908076401", 
      RowBox[{
      "-", "0.2456743536496151757757937224889116846036764401207361412666485562\
56329056354726504259359638478342`63.81835908076401"}], 
      "0.024928063854678876736726424207538683628942800250497235964568546859687\
541428049050794403290977695`63.81835908076402"},
     {
      RowBox[{
      "-", "0.1036658044602740983010797465641109867859450376677572439433027721\
12108701244820930068683294624057`63.81835908076401"}], 
      "0.320479808802085829087911099034405211114666767765790225528310742560298\
871944809512103479551994129`63.81835908076401", 
      "0.019622135746464034092096399628813170705319259210498083208079764352158\
239439811689673030068874795`63.81835908076402"},
     {
      RowBox[{
      "-", "0.0024853302632041137815566313906889929379078778331738977128817607\
25934911428305663960825485153801`63.81835908076401"}], 
      RowBox[{
      "-", "0.0927277433021352928671983882892348540926135026458856857786324395\
59442071156499063396103929489591`63.81835908076401"}], 
      "0.052151719260423327515667033231853376330571436944777265826668118528316\
820191389854369777005597934`63.81835908076401"}
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
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"largest", " ", "entries"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"Max", "[", 
    RowBox[{"Flatten", "[", 
     RowBox[{"Abs", "[", 
      RowBox[{"N", "[", 
       SubscriptBox["G", 
        RowBox[{"u", ",", "1"}]], "]"}], "]"}], "]"}], "]"}], 
   "\[IndentingNewLine]", 
   RowBox[{"Max", "[", 
    RowBox[{"Flatten", "[", 
     RowBox[{"Abs", "[", 
      RowBox[{"N", "[", 
       SubscriptBox["G", 
        RowBox[{"d", ",", "1"}]], "]"}], "]"}], "]"}], "]"}]}]}]], "Input"],

Cell[BoxData["1.8021593072080877`"], "Output"],

Cell[BoxData["1.8021593072080877`"], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "determinants", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["detG", 
     RowBox[{"u", ",", "1"}]], "=", 
    RowBox[{"Det", "[", 
     SubscriptBox["G", 
      RowBox[{"u", ",", "1"}]], "]"}]}], "\[IndentingNewLine]", 
   RowBox[{
    SubscriptBox["detG", 
     RowBox[{"d", ",", "1"}]], "=", 
    RowBox[{"Det", "[", 
     SubscriptBox["G", 
      RowBox[{"d", ",", "1"}]], "]"}]}]}]}]], "Input"],

Cell[BoxData["3.\
955411089731725005533059260659054402767568290068045199888736471543807746827489\
5570901067073847999`63.81835908076401*^-38"], "Output"],

Cell[BoxData["3.\
548952760849338968980347852770132875881474775435757238410508418568203667955546\
03480399768888598`63.81835908076401*^-40"], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"save", " ", "updated", " ", 
    RowBox[{"Green", "'"}], "s", " ", "functions", " ", "as", " ", 
    "reference", " ", "to", " ", "disk"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{"Table", "[", 
     RowBox[{
      RowBox[{"Export", "[", 
       RowBox[{
        RowBox[{
         SubscriptBox["fn", "base"], "<>", "\"\<_G\>\"", "<>", 
         RowBox[{"ToString", "[", "\[Sigma]", "]"}], "<>", "\"\<1.dat\>\""}], 
        ",", 
        RowBox[{"Flatten", "[", 
         RowBox[{"Transpose", "[", 
          RowBox[{"N", "[", 
           SubscriptBox["G", 
            RowBox[{"\[Sigma]", ",", "1"}]], "]"}], "]"}], "]"}], ",", 
        "\"\<Real64\>\""}], "]"}], ",", 
      RowBox[{"{", 
       RowBox[{"\[Sigma]", ",", 
        RowBox[{"{", 
         RowBox[{"u", ",", "d"}], "}"}]}], "}"}]}], "]"}], ";"}], 
   "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"Table", "[", 
     RowBox[{
      RowBox[{"Export", "[", 
       RowBox[{
        RowBox[{
         SubscriptBox["fn", "base"], "<>", "\"\<_detG\>\"", "<>", 
         RowBox[{"ToString", "[", "\[Sigma]", "]"}], "<>", "\"\<1.dat\>\""}], 
        ",", 
        RowBox[{"N", "[", 
         SubscriptBox["detG", 
          RowBox[{"\[Sigma]", ",", "1"}]], "]"}], ",", "\"\<Real64\>\""}], 
       "]"}], ",", 
      RowBox[{"{", 
       RowBox[{"\[Sigma]", ",", 
        RowBox[{"{", 
         RowBox[{"u", ",", "d"}], "}"}]}], "}"}]}], "]"}], ";"}]}]}]], "Input"]
}, Open  ]]
},
WindowSize->{1475, 969},
WindowMargins->{{Automatic, 168}, {Automatic, 75}},
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
Cell[1486, 35, 40, 0, 43, "Subsection"],
Cell[CellGroupData[{
Cell[1551, 39, 326, 10, 72, "Input"],
Cell[1880, 51, 31, 0, 31, "Output"]
}, Open  ]],
Cell[1926, 54, 219, 7, 52, "Input"],
Cell[2148, 63, 261, 8, 67, "Input"],
Cell[CellGroupData[{
Cell[2434, 75, 337, 11, 52, "Input"],
Cell[2774, 88, 46, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[2857, 93, 295, 9, 69, "Input"],
Cell[3155, 104, 29, 0, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[3233, 110, 41, 0, 43, "Subsection"],
Cell[3277, 112, 320, 10, 72, "Input"],
Cell[CellGroupData[{
Cell[3622, 126, 947, 29, 72, "Input"],
Cell[4572, 157, 1469, 50, 52, "Output"],
Cell[6044, 209, 74, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[6155, 216, 295, 8, 52, "Input"],
Cell[6453, 226, 29, 0, 31, "Output"]
}, Open  ]],
Cell[6497, 229, 1075, 30, 52, "Input"],
Cell[CellGroupData[{
Cell[7597, 263, 333, 10, 72, "Input"],
Cell[7933, 275, 28, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[7998, 280, 355, 10, 52, "Input"],
Cell[8356, 292, 28, 0, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[8433, 298, 45, 0, 43, "Subsection"],
Cell[8481, 300, 1023, 32, 72, "Input"],
Cell[CellGroupData[{
Cell[9529, 336, 442, 13, 52, "Input"],
Cell[9974, 351, 4754, 143, 129, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[14765, 499, 511, 15, 52, "Input"],
Cell[15279, 516, 787, 18, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[16103, 539, 371, 10, 52, "Input"],
Cell[16477, 551, 28, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[16542, 556, 874, 28, 72, "Input"],
Cell[17419, 586, 49, 0, 31, "Output"],
Cell[17471, 588, 50, 0, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[17570, 594, 56, 0, 43, "Subsection"],
Cell[CellGroupData[{
Cell[17651, 598, 556, 18, 72, "Input"],
Cell[18210, 618, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[18322, 625, 388, 11, 52, "Input"],
Cell[18713, 638, 407, 15, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[19157, 658, 305, 10, 31, "Input"],
Cell[19465, 670, 4100, 54, 152, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[23614, 730, 31, 0, 43, "Subsection"],
Cell[23648, 732, 625, 17, 52, "Input"],
Cell[24276, 751, 843, 22, 52, "Input"],
Cell[25122, 775, 513, 13, 72, "Input"],
Cell[25638, 790, 9337, 215, 532, "Input"],
Cell[34978, 1007, 242, 7, 31, "Input"],
Cell[CellGroupData[{
Cell[35245, 1018, 765, 22, 72, "Input"],
Cell[36013, 1042, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[36125, 1049, 393, 11, 52, "Input"],
Cell[36521, 1062, 460, 14, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[37018, 1081, 665, 19, 72, "Input"],
Cell[37686, 1102, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[37798, 1109, 380, 11, 52, "Input"],
Cell[38181, 1122, 28, 0, 31, "Output"]
}, Open  ]],
Cell[38224, 1125, 835, 25, 31, "Input"],
Cell[CellGroupData[{
Cell[39084, 1154, 867, 25, 52, "Input"],
Cell[39954, 1181, 28, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[40019, 1186, 388, 11, 52, "Input"],
Cell[40410, 1199, 323, 11, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[40770, 1215, 305, 10, 31, "Input"],
Cell[41078, 1227, 4100, 54, 152, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[45215, 1286, 877, 27, 72, "Input"],
Cell[46095, 1315, 1813, 39, 71, "Output"],
Cell[47911, 1356, 1814, 39, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[49762, 1400, 559, 17, 72, "Input"],
Cell[50324, 1419, 46, 0, 31, "Output"],
Cell[50373, 1421, 46, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[50456, 1426, 476, 15, 72, "Input"],
Cell[50935, 1443, 153, 2, 31, "Output"],
Cell[51091, 1447, 151, 2, 31, "Output"]
}, Open  ]],
Cell[51257, 1452, 1513, 43, 72, "Input"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature WvDLkdB6qnMrRAKTVcanUzwO *)
