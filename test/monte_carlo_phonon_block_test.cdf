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
NotebookDataLength[     43101,       1433]
NotebookOptionsPosition[     40088,       1308]
NotebookOutlinePosition[     40431,       1323]
CellTagsIndexPosition[     40388,       1320]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["General parameters", "Subsection"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"Coulomb", " ", "coupling", " ", "constant"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["U", "val"], "=", "4"}], ";"}]}]], "Input"],

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

Cell[BoxData["0.7369045906209687`"], "Output"]
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
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"phonon", " ", "frequency"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["\[CapitalOmega]", "val"], "=", 
     RowBox[{"6", "/", "5"}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"N", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData["1.2`"], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"electron", "-", 
    RowBox[{"phonon", " ", "interaction", " ", "strength"}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["g", "val"], "=", 
     RowBox[{"13", "/", "20"}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"N", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData["0.65`"], "Output"]
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
    SubscriptBox["exp", "\[Tau]k"], "=", 
    RowBox[{"FullSimplify", "[", 
     RowBox[{"MatrixExp", "[", 
      RowBox[{
       RowBox[{"-", 
        SubscriptBox["\[Tau]", "val"]}], " ", 
       RowBox[{"(", 
        RowBox[{"-", 
         SubscriptBox["latt", "neigh"]}], ")"}]}], "]"}], "]"}]}], 
   ";"}]}]], "Input"],

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
  RowBox[{"(*", " ", 
   RowBox[{"compare", " ", "with", " ", "numerical", " ", "evaluation"}], " ",
    "*)"}], "\[IndentingNewLine]", 
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
         SubscriptBox["latt", "neigh"]}], ")"}]}], "]"}], "]"}]}], 
   "]"}]}]], "Input"],

Cell[BoxData["9.59646170837076`*^-16"], "Output"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Hubbard-Stratonovich field", "Subsection"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"Hubbard", "-", 
    RowBox[{
    "Stratonovich", " ", "field", " ", "remains", " ", "constant", " ", 
     "during", " ", "phonon", " ", "block", " ", "updates"}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["s", "val"], "=", 
     RowBox[{
      RowBox[{
       RowBox[{"Transpose", "[", 
        RowBox[{"Partition", "[", 
         RowBox[{
          RowBox[{"Import", "[", 
           RowBox[{
            RowBox[{
             RowBox[{"NotebookDirectory", "[", "]"}], "<>", 
             RowBox[{"FileBaseName", "[", 
              RowBox[{"NotebookFileName", "[", "]"}], "]"}], "<>", 
             "\"\<_s.dat\>\""}], ",", "\"\<Integer8\>\""}], "]"}], ",", 
          SubscriptBox["n", "sites"]}], "]"}], "]"}], "/.", 
       RowBox[{"{", 
        RowBox[{"1", "\[Rule]", 
         RowBox[{"-", "1"}]}], "}"}]}], "/.", 
      RowBox[{"{", 
       RowBox[{"0", "\[Rule]", "1"}], "}"}]}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"Dimensions", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"24", ",", "16"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   SubscriptBox["s", "val"], "\[LeftDoubleBracket]", 
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
     RowBox[{"-", "1"}], ",", 
     RowBox[{"-", "1"}]}], "}"}]}], "}"}]], "Output"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Phonon block update", "Subsection"],

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
   RowBox[{"expK_", ",", "\[Lambda]s_", ",", "\[Tau]gX_"}], "]"}], ":=", 
  RowBox[{"{", "\n", "\t", 
   RowBox[{
    RowBox[{"InitializeGreensFunction", "[", 
     RowBox[{"expK", ",", 
      RowBox[{"Exp", "[", 
       RowBox[{
        RowBox[{"-", "\[Lambda]s"}], "-", "\[Tau]gX"}], "]"}]}], "]"}], ",", 
    "\n", "\t", 
    RowBox[{"InitializeGreensFunction", "[", 
     RowBox[{"expK", ",", 
      RowBox[{"Exp", "[", 
       RowBox[{
        RowBox[{"+", "\[Lambda]s"}], "-", "\[Tau]gX"}], "]"}]}], "]"}]}], 
   "}"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"PhononBlockUpdate", "[", 
   RowBox[{
   "expK_", ",", "\[Tau]_", ",", "\[Lambda]_", ",", "\[CapitalOmega]_", ",", 
    "g_", ",", "s_", ",", "X0_", ",", "i_Integer", ",", "\[CapitalDelta]x_"}],
    "]"}], ":=", 
  RowBox[{"Module", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
     "nsites", ",", "L", ",", "Gu0", ",", "Gd0", ",", "Gu1", ",", "Gd1", ",", 
      "\[CapitalDelta]Eph", ",", "X1"}], "}"}], ",", "\[IndentingNewLine]", 
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
       RowBox[{"Gu0", ",", "Gd0"}], "}"}], "=", 
      RowBox[{"InitializeUpDownGreensFunctions", "[", 
       RowBox[{"expK", ",", 
        RowBox[{"\[Lambda]", " ", "s"}], ",", 
        RowBox[{"\[Tau]", " ", "g", " ", "X0"}]}], "]"}]}], ";", 
     "\[IndentingNewLine]", 
     RowBox[{"Print", "[", 
      RowBox[{"\"\<Det[Gu0]Det[Gd0]: \>\"", ",", 
       RowBox[{
        RowBox[{"Det", "[", "Gu0", "]"}], 
        RowBox[{"Det", "[", "Gd0", "]"}]}]}], "]"}], ";", 
     "\[IndentingNewLine]", 
     RowBox[{"(*", " ", 
      RowBox[{
      "calculate", " ", "change", " ", "in", " ", "phonon", " ", "energy"}], 
      " ", "*)"}], "\[IndentingNewLine]", 
     RowBox[{"\[CapitalDelta]Eph", "=", 
      RowBox[{
       FractionBox["1", "2"], 
       SuperscriptBox["\[CapitalOmega]", "2"], 
       RowBox[{"Total", "[", 
        RowBox[{
         SuperscriptBox["\[CapitalDelta]x", "2"], "+", 
         RowBox[{"2", " ", "\[CapitalDelta]x", " ", 
          RowBox[{"X0", "\[LeftDoubleBracket]", 
           RowBox[{"i", ",", ";;"}], "\[RightDoubleBracket]"}]}]}], "]"}]}]}],
      ";", "\[IndentingNewLine]", 
     RowBox[{"X1", "=", "X0"}], ";", 
     RowBox[{
      RowBox[{"X1", "\[LeftDoubleBracket]", 
       RowBox[{"i", ",", ";;"}], "\[RightDoubleBracket]"}], "+=", 
      "\[CapitalDelta]x"}], ";", "\[IndentingNewLine]", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{"Gu1", ",", "Gd1"}], "}"}], "=", 
      RowBox[{"InitializeUpDownGreensFunctions", "[", 
       RowBox[{"expK", ",", 
        RowBox[{"\[Lambda]", " ", "s"}], ",", 
        RowBox[{"\[Tau]", " ", "g", " ", "X1"}]}], "]"}]}], ";", 
     "\[IndentingNewLine]", 
     RowBox[{"Print", "[", 
      RowBox[{"\"\<Det[Gu1]Det[Gd1]: \>\"", ",", 
       RowBox[{
        RowBox[{"Det", "[", "Gu1", "]"}], 
        RowBox[{"Det", "[", "Gd1", "]"}]}]}], "]"}], ";", 
     "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{
        FractionBox[
         RowBox[{
          RowBox[{"Det", "[", "Gu0", "]"}], 
          RowBox[{"Det", "[", "Gd0", "]"}]}], 
         RowBox[{
          RowBox[{"Det", "[", "Gu1", "]"}], 
          RowBox[{"Det", "[", "Gd1", "]"}]}]], 
        RowBox[{"Exp", "[", 
         RowBox[{
          RowBox[{"-", "\[Tau]"}], " ", "\[CapitalDelta]Eph"}], "]"}]}], ",", 
       "X1", ",", "Gu1", ",", "Gd1"}], "}"}]}]}], "]"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"initial", " ", "phonon", " ", "field"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["X", "0"], "=", 
     RowBox[{"Transpose", "[", 
      RowBox[{"Partition", "[", 
       RowBox[{
        RowBox[{"SetPrecision", "[", 
         RowBox[{
          RowBox[{"Import", "[", 
           RowBox[{
            RowBox[{
             RowBox[{"NotebookDirectory", "[", "]"}], "<>", 
             RowBox[{"FileBaseName", "[", 
              RowBox[{"NotebookFileName", "[", "]"}], "]"}], "<>", 
             "\"\<_X0.dat\>\""}], ",", "\"\<Real64\>\""}], "]"}], ",", 
          RowBox[{"4", "MachinePrecision"}]}], "]"}], ",", 
        SubscriptBox["n", "sites"]}], "]"}], "]"}]}], ";"}], 
   "\[IndentingNewLine]", 
   RowBox[{"Dimensions", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"24", ",", "16"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{"N", "[", 
   RowBox[{
    SubscriptBox["X", "0"], "\[LeftDoubleBracket]", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"1", ",", "2", ",", 
       RowBox[{"-", "1"}]}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"1", ",", "2", ",", 
       RowBox[{"-", "1"}]}], "}"}]}], "\[RightDoubleBracket]"}], 
   "]"}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "0.33784826443736615`"}], ",", 
     RowBox[{"-", "0.10214103727985746`"}], ",", 
     RowBox[{"-", "0.7958095962716176`"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "0.6065976708421408`"}], ",", 
     RowBox[{"-", "0.5097989029159002`"}], ",", 
     RowBox[{"-", "0.9241493051683802`"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "0.4280929630257604`"}], ",", "0.07925917911592362`", ",", 
     RowBox[{"-", "0.3791353696242037`"}]}], "}"}]}], "}"}]], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"initial", " ", 
    RowBox[{"Green", "'"}], "s", " ", "function", " ", "matrices"}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"u", ",", "0"}]], ",", 
      SubscriptBox["G", 
       RowBox[{"d", ",", "0"}]]}], "}"}], "=", 
    RowBox[{"InitializeUpDownGreensFunctions", "[", 
     RowBox[{
      SubscriptBox["exp", "\[Tau]k"], ",", 
      RowBox[{
       SubscriptBox["\[Lambda]", "val"], 
       SubscriptBox["s", "val"]}], ",", 
      RowBox[{
       SubscriptBox["\[Tau]", "val"], 
       SubscriptBox["g", "val"], 
       SubscriptBox["X", "0"]}]}], "]"}]}], ";"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"N", "[", 
    RowBox[{
     SubscriptBox["G", 
      RowBox[{"u", ",", "0"}]], "\[LeftDoubleBracket]", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{"1", ",", "2", ",", 
        RowBox[{"-", "1"}]}], "}"}], ",", 
      RowBox[{"{", 
       RowBox[{"1", ",", "2", ",", 
        RowBox[{"-", "1"}]}], "}"}]}], "\[RightDoubleBracket]"}], "]"}], 
   "\[IndentingNewLine]", 
   RowBox[{"N", "[", 
    RowBox[{
     SubscriptBox["G", 
      RowBox[{"d", ",", "0"}]], "\[LeftDoubleBracket]", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{"1", ",", "2", ",", 
        RowBox[{"-", "1"}]}], "}"}], ",", 
      RowBox[{"{", 
       RowBox[{"1", ",", "2", ",", 
        RowBox[{"-", "1"}]}], "}"}]}], "\[RightDoubleBracket]"}], 
    "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"0.5078994792710819`", ",", 
     RowBox[{"-", "0.5238539833337109`"}], ",", "0.1024363020759369`"}], 
    "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "0.2978093500754256`"}], ",", "0.5521626972723901`", ",", 
     "0.06284183194252611`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
    "0.031867680277268426`", ",", "0.15626929718613192`", ",", 
     "0.5136164528922846`"}], "}"}]}], "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"0.39717198365861106`", ",", 
     RowBox[{"-", "0.2865281832170811`"}], ",", 
     RowBox[{"-", "0.006994119545551892`"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "0.41371397623344486`"}], ",", "0.37677442066641226`", ",", 
     "0.08705190192597781`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "0.1001790779335859`"}], ",", "0.07317579596294306`", ",", 
     "0.3888355965471296`"}], "}"}]}], "}"}]], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["\[CapitalDelta]x", "0"], "=", 
   RowBox[{"SetPrecision", "[", 
    RowBox[{"0.90005446944002498", ",", 
     RowBox[{"4", "MachinePrecision"}]}], "]"}]}], ";"}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
   "perform", " ", "a", " ", "phonon", " ", "block", " ", "update", " ", 
    "step"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      SubscriptBox["p", "val"], ",", 
      SubscriptBox["X", "1"], ",", 
      SubscriptBox["G", 
       RowBox[{"u", ",", "1"}]], ",", 
      SubscriptBox["G", 
       RowBox[{"d", ",", "1"}]]}], "}"}], "=", 
    RowBox[{"PhononBlockUpdate", "[", 
     RowBox[{
      SubscriptBox["exp", "\[Tau]k"], ",", 
      SubscriptBox["\[Tau]", "val"], ",", 
      SubscriptBox["\[Lambda]", "val"], ",", 
      SubscriptBox["\[CapitalOmega]", "val"], ",", 
      SubscriptBox["g", "val"], ",", 
      SubscriptBox["s", "val"], ",", 
      SubscriptBox["X", "0"], ",", "16", ",", 
      SubscriptBox["\[CapitalDelta]x", "0"]}], "]"}]}], ";"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 InterpretationBox[
  RowBox[{"\<\"Det[Gu0]Det[Gd0]: \"\>", "\[InvisibleSpace]", 
   "1.400382061609800841128287888811449974394044061017193333517336391413925862\
0433535209255949535`63.203585744828366*^-82"}],
  SequenceForm[
  "Det[Gu0]Det[Gd0]: ", 
   1.4003820616098008411282878888114499743940440610171933335173363914139258620\
433535209255949535`63.203585744828366*^-82],
  Editable->False]], "Print"],

Cell[BoxData[
 InterpretationBox[
  RowBox[{"\<\"Det[Gu1]Det[Gd1]: \"\>", "\[InvisibleSpace]", 
   "4.471127538211942453333437541576928396021996902505968431023687696367430754\
7576277904401152864`63.19390085434241*^-82"}],
  SequenceForm[
  "Det[Gu1]Det[Gd1]: ", 
   4.4711275382119424533334375415769283960219969025059684310236876963674307547\
576277904401152864`63.19390085434241*^-82],
  Editable->False]], "Print"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"acceptance", " ", "probability"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  SubscriptBox["p", "val"]}]], "Input"],

Cell[BoxData["0.\
218590421799713732923128000910264844414710104823480170638353358238736261482891\
324592508219118624`62.65885308370692"], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"N", "[", 
    RowBox[{
     SubscriptBox["X", "1"], "\[LeftDoubleBracket]", 
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
     {
      RowBox[{"-", "0.33784826443736615`"}], 
      RowBox[{"-", "0.10214103727985746`"}], 
      RowBox[{"-", "0.7958095962716176`"}]},
     {
      RowBox[{"-", "0.6065976708421408`"}], 
      RowBox[{"-", "0.5097989029159002`"}], 
      RowBox[{"-", "0.9241493051683802`"}]},
     {
      RowBox[{"-", "0.4280929630257604`"}], "0.07925917911592362`", 
      RowBox[{"-", "0.3791353696242037`"}]}
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
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{"N", "[", 
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"u", ",", "1"}]], "\[LeftDoubleBracket]", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{"1", ",", "2", ",", 
         RowBox[{"-", "1"}]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{"1", ",", "2", ",", 
         RowBox[{"-", "1"}]}], "}"}]}], "\[RightDoubleBracket]"}], "]"}], "//",
     "MatrixForm"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"N", "[", 
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"d", ",", "1"}]], "\[LeftDoubleBracket]", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{"1", ",", "2", ",", 
         RowBox[{"-", "1"}]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{"1", ",", "2", ",", 
         RowBox[{"-", "1"}]}], "}"}]}], "\[RightDoubleBracket]"}], "]"}], "//",
     "MatrixForm"}]}]}]], "Input"],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"0.5079168178793175`", 
      RowBox[{"-", "0.5241859402230448`"}], "0.1026147319081147`"},
     {
      RowBox[{"-", "0.2979606268958537`"}], "0.5521722639478813`", 
      "0.06268709259205553`"},
     {"0.03189991461787961`", "0.1562722362753481`", "0.5136488290457388`"}
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
     {"0.3972103957163743`", 
      RowBox[{"-", "0.2865579740212162`"}], 
      RowBox[{"-", "0.006969135655977277`"}]},
     {
      RowBox[{"-", "0.4137949885600493`"}], "0.3768775883346961`", 
      "0.08690629147758558`"},
     {
      RowBox[{"-", "0.1005326998905071`"}], "0.07344812883936919`", 
      "0.3885956701406351`"}
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

Cell[BoxData["0.9579408370230565`"], "Output"],

Cell[BoxData["0.9373403266556257`"], "Output"]
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
138442990021162709089779408595020188338881015116038053851864980189333292247162\
6410630224135372`63.49302384123082*^-41"], "Output"],

Cell[BoxData["1.\
424632390146361515555486630637116736597251520735139753200834568937583006272575\
5865636661477572`63.49684626949826*^-41"], "Output"]
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
   RowBox[{
   "save", " ", "updated", " ", "phonon", " ", "field", " ", "and", " ", 
    "corresponding", " ", "entrywise", " ", "exponential", " ", "to", " ", 
    "disk"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{"Export", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["fn", "base"], "<>", "\"\<_X1.dat\>\""}], ",", 
      RowBox[{"Flatten", "[", 
       RowBox[{"Transpose", "[", 
        RowBox[{"N", "[", 
         SubscriptBox["X", "1"], "]"}], "]"}], "]"}], ",", "\"\<Real64\>\""}],
      "]"}], ";"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"Export", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["fn", "base"], "<>", "\"\<_expX1.dat\>\""}], ",", 
      RowBox[{"Flatten", "[", 
       RowBox[{"Transpose", "[", 
        RowBox[{"N", "[", 
         RowBox[{"Exp", "[", 
          RowBox[{
           RowBox[{"-", 
            SubscriptBox["\[Tau]", "val"]}], 
           SubscriptBox["g", "val"], 
           SubscriptBox["X", "1"]}], "]"}], "]"}], "]"}], "]"}], ",", 
      "\"\<Real64\>\""}], "]"}], ";"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"save", " ", "updated", " ", 
    RowBox[{"Green", "'"}], "s", " ", "functions", " ", "to", " ", "disk"}], 
   " ", "*)"}], "\[IndentingNewLine]", 
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
WindowSize->{1511, 997},
WindowMargins->{{253, Automatic}, {41, Automatic}},
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
Cell[1529, 37, 227, 7, 52, "Input"],
Cell[1759, 46, 219, 7, 52, "Input"],
Cell[1981, 55, 261, 8, 67, "Input"],
Cell[CellGroupData[{
Cell[2267, 67, 337, 11, 52, "Input"],
Cell[2607, 80, 46, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[2690, 85, 295, 9, 69, "Input"],
Cell[2988, 96, 29, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[3054, 101, 320, 9, 72, "Input"],
Cell[3377, 112, 31, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[3445, 117, 362, 11, 72, "Input"],
Cell[3810, 130, 32, 0, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[3891, 136, 41, 0, 43, "Subsection"],
Cell[3935, 138, 320, 10, 72, "Input"],
Cell[CellGroupData[{
Cell[4280, 152, 947, 29, 72, "Input"],
Cell[5230, 183, 1469, 50, 52, "Output"],
Cell[6702, 235, 74, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[6813, 242, 295, 8, 52, "Input"],
Cell[7111, 252, 29, 0, 31, "Output"]
}, Open  ]],
Cell[7155, 255, 1075, 30, 52, "Input"],
Cell[CellGroupData[{
Cell[8255, 289, 333, 10, 72, "Input"],
Cell[8591, 301, 28, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[8656, 306, 355, 10, 52, "Input"],
Cell[9014, 318, 28, 0, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[9091, 324, 45, 0, 43, "Subsection"],
Cell[9139, 326, 650, 20, 52, "Input"],
Cell[CellGroupData[{
Cell[9814, 350, 442, 13, 52, "Input"],
Cell[10259, 365, 4754, 143, 129, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[15050, 513, 511, 15, 52, "Input"],
Cell[15564, 530, 787, 18, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[16388, 553, 504, 16, 52, "Input"],
Cell[16895, 571, 49, 0, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[16993, 577, 48, 0, 43, "Subsection"],
Cell[CellGroupData[{
Cell[17066, 581, 1090, 30, 72, "Input"],
Cell[18159, 613, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[18271, 620, 390, 11, 52, "Input"],
Cell[18664, 633, 344, 12, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[19057, 651, 41, 0, 43, "Subsection"],
Cell[19101, 653, 625, 17, 52, "Input"],
Cell[19729, 672, 843, 22, 52, "Input"],
Cell[20575, 696, 612, 17, 72, "Input"],
Cell[21190, 715, 3220, 83, 273, "Input"],
Cell[CellGroupData[{
Cell[24435, 802, 870, 24, 72, "Input"],
Cell[25308, 828, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[25420, 835, 429, 13, 52, "Input"],
Cell[25852, 850, 593, 16, 31, "Output"]
}, Open  ]],
Cell[26460, 869, 735, 23, 52, "Input"],
Cell[CellGroupData[{
Cell[27220, 896, 866, 27, 72, "Input"],
Cell[28089, 925, 490, 14, 31, "Output"],
Cell[28582, 941, 525, 14, 31, "Output"]
}, Open  ]],
Cell[29122, 958, 224, 6, 31, "Input"],
Cell[CellGroupData[{
Cell[29371, 968, 884, 25, 52, "Input"],
Cell[CellGroupData[{
Cell[30280, 997, 419, 9, 23, "Print"],
Cell[30702, 1008, 417, 9, 23, "Print"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[31168, 1023, 172, 5, 52, "Input"],
Cell[31343, 1030, 147, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[31527, 1037, 472, 14, 52, "Input"],
Cell[32002, 1053, 958, 26, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[32997, 1084, 957, 29, 72, "Input"],
Cell[33957, 1115, 832, 21, 71, "Output"],
Cell[34792, 1138, 885, 24, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[35714, 1167, 559, 17, 72, "Input"],
Cell[36276, 1186, 46, 0, 31, "Output"],
Cell[36325, 1188, 46, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[36408, 1193, 476, 15, 72, "Input"],
Cell[36887, 1210, 150, 2, 31, "Output"],
Cell[37040, 1214, 150, 2, 31, "Output"]
}, Open  ]],
Cell[37205, 1219, 242, 7, 31, "Input"],
Cell[37450, 1228, 1136, 32, 72, "Input"],
Cell[38589, 1262, 1483, 43, 72, "Input"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature duD77DbgrQwtnCKAnbohIjQa *)
