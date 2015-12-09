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
NotebookDataLength[     56693,       1765]
NotebookOptionsPosition[     53827,       1649]
NotebookOutlinePosition[     54170,       1664]
CellTagsIndexPosition[     54127,       1661]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["General parameters", "Subsection"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
   "number", " ", "of", " ", "orbitals", " ", "per", " ", "unit", " ", 
    "cell"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["n", "orb"], "=", "2"}], ";"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"Coulomb", " ", "coupling", " ", "constants"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["U", "val"], "=", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"9", "/", "2"}], ",", 
       RowBox[{"11", "/", "3"}]}], "}"}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"N", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"4.5`", ",", "3.6666666666666665`"}], "}"}]], "Output"]
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

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["fn", "base"], "=", 
   RowBox[{
    RowBox[{"NotebookDirectory", "[", "]"}], "<>", 
    RowBox[{"FileBaseName", "[", 
     RowBox[{"NotebookFileName", "[", "]"}], "]"}]}]}], ";"}]], "Input"]
}, Open  ]],

Cell[CellGroupData[{

Cell["Rectangular lattice and kinetic operator", "Subsection"],

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
     SubscriptBox["n", "y"], "=", "5"}], ";"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["n", "sites"], "=", 
   RowBox[{
    SubscriptBox["n", "x"], 
    SubscriptBox["n", "y"]}]}], ";"}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{
    "combined", " ", "lattice", " ", "and", " ", "orbital", " ", "indices"}], 
    ",", " ", 
    RowBox[{
     RowBox[{"in", " ", "column"}], "-", 
     RowBox[{"major", " ", 
      RowBox[{"(", "Fortran", ")"}], " ", "order"}]}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["orb", "ind"], "=", 
    RowBox[{"Flatten", "[", 
     RowBox[{
      RowBox[{"Transpose", "[", 
       RowBox[{
        RowBox[{"Outer", "[", 
         RowBox[{
          RowBox[{
           RowBox[{"{", "##", "}"}], "&"}], ",", 
          RowBox[{"Range", "[", 
           RowBox[{"0", ",", 
            RowBox[{
             SubscriptBox["n", "orb"], "-", "1"}]}], "]"}], ",", 
          RowBox[{"Range", "[", 
           RowBox[{"0", ",", 
            RowBox[{
             SubscriptBox["n", "x"], "-", "1"}]}], "]"}], ",", 
          RowBox[{"Range", "[", 
           RowBox[{"0", ",", 
            RowBox[{
             SubscriptBox["n", "y"], "-", "1"}]}], "]"}]}], "]"}], ",", 
        RowBox[{"{", 
         RowBox[{"1", ",", "3", ",", "2"}], "}"}]}], "]"}], ",", "2"}], 
     "]"}]}], "\[IndentingNewLine]", 
   RowBox[{"Dimensions", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"0", ",", "0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "1", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "2", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "3", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "1", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "2", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "3", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "1", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "2", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "3", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "1", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "2", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "3", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "1", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "2", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "3", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "1", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "2", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "3", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "0", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "1", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "2", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "3", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "0", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "1", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "2", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "3", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "0", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "1", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "2", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "3", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "0", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "1", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "2", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "3", ",", "4"}], "}"}]}], "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"40", ",", "3"}], "}"}]], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"hopping", " ", "parameters"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["t", 
      RowBox[{"{", 
       RowBox[{"0", ",", "0"}], "}"}]], "=", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{"0", ",", 
         RowBox[{"-", 
          FractionBox["1", "3"]}]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{"0", ",", "0"}], "}"}]}], "}"}]}], ";"}], 
   "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{
     SubscriptBox["t", 
      RowBox[{"{", 
       RowBox[{"1", ",", "0"}], "}"}]], "=", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{
         FractionBox["4", "9"], ",", 
         FractionBox["3", "16"]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         FractionBox["3", "7"], ",", 
         RowBox[{"-", 
          FractionBox["2", "19"]}]}], "}"}]}], "}"}]}], ";"}], 
   "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{
     SubscriptBox["t", 
      RowBox[{"{", 
       RowBox[{"0", ",", "1"}], "}"}]], "=", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{
         RowBox[{"-", 
          FractionBox["1", "10"]}], ",", 
         FractionBox["2", "17"]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         FractionBox["2", "9"], ",", 
         RowBox[{"-", 
          FractionBox["3", "7"]}]}], "}"}]}], "}"}]}], ";"}], 
   "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{
     SubscriptBox["t", 
      RowBox[{"{", 
       RowBox[{"1", ",", "1"}], "}"}]], "=", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{
         FractionBox["1", "13"], ",", 
         FractionBox["4", "5"]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         FractionBox["3", "19"], ",", 
         FractionBox["1", "6"]}], "}"}]}], "}"}]}], ";"}], 
   "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{
     SubscriptBox["t", 
      RowBox[{"{", 
       RowBox[{"1", ",", 
        RowBox[{"-", "1"}]}], "}"}]], "=", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{
         FractionBox["1", "5"], ",", 
         RowBox[{"-", 
          FractionBox["1", "9"]}]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"-", 
          FractionBox["3", "11"]}], ",", 
         FractionBox["1", "7"]}], "}"}]}], "}"}]}], ";"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"chemical", " ", "potential"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[Mu]", "val"], "=", 
    FractionBox["5", "7"]}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"site", " ", "energies", " ", "for", " ", "each", " ", "orbital"}],
    " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[Epsilon]", "val"], "=", 
    RowBox[{"{", 
     RowBox[{
      RowBox[{
       RowBox[{"-", "1"}], "/", "5"}], ",", 
      RowBox[{"2", "/", "11"}]}], "}"}]}], ";"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"kinetic", " ", "operator"}], ";", " ", 
    RowBox[{"factor", " ", 
     RowBox[{"(", 
      RowBox[{"-", "1"}], ")"}], " ", "from", " ", "negative", " ", "sign", 
     " ", "in", " ", "Hamiltonian"}]}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["K", "val"], "=", 
     RowBox[{"Module", "[", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{"\[CapitalDelta]l", ",", "i", ",", "j"}], "}"}], ",", 
       RowBox[{
        RowBox[{"-", 
         RowBox[{"Outer", "[", 
          RowBox[{
           RowBox[{
            RowBox[{"(", 
             RowBox[{
              RowBox[{"\[CapitalDelta]l", "=", 
               RowBox[{"{", 
                RowBox[{
                 RowBox[{"Mod", "[", 
                  RowBox[{
                   RowBox[{
                    RowBox[{"(", 
                    RowBox[{"#2", "-", "#1"}], ")"}], "\[LeftDoubleBracket]", 
                    "2", "\[RightDoubleBracket]"}], ",", 
                   SubscriptBox["n", "x"], ",", 
                   RowBox[{"-", "1"}]}], "]"}], ",", 
                 RowBox[{"Mod", "[", 
                  RowBox[{
                   RowBox[{
                    RowBox[{"(", 
                    RowBox[{"#2", "-", "#1"}], ")"}], "\[LeftDoubleBracket]", 
                    "3", "\[RightDoubleBracket]"}], ",", 
                   SubscriptBox["n", "y"], ",", 
                   RowBox[{"-", "1"}]}], "]"}]}], "}"}]}], ";", 
              RowBox[{"i", "=", 
               RowBox[{
                RowBox[{"First", "[", "#1", "]"}], "+", "1"}]}], ";", 
              RowBox[{"j", "=", 
               RowBox[{
                RowBox[{"First", "[", "#2", "]"}], "+", "1"}]}], ";", 
              RowBox[{"If", "[", 
               RowBox[{
                RowBox[{
                 RowBox[{
                  RowBox[{"First", "[", "\[CapitalDelta]l", "]"}], "<", "0"}],
                  "\[Or]", 
                 RowBox[{"(", 
                  RowBox[{
                   RowBox[{
                    RowBox[{"First", "[", "\[CapitalDelta]l", "]"}], 
                    "\[Equal]", "0"}], "\[And]", 
                   RowBox[{
                    RowBox[{"Last", "[", "\[CapitalDelta]l", "]"}], "<", 
                    "0"}]}], ")"}]}], ",", 
                RowBox[{
                 RowBox[{"\[CapitalDelta]l", "=", 
                  RowBox[{"-", "\[CapitalDelta]l"}]}], ";", 
                 RowBox[{
                  RowBox[{"{", 
                   RowBox[{"i", ",", "j"}], "}"}], "=", 
                  RowBox[{"{", 
                   RowBox[{"j", ",", "i"}], "}"}]}]}]}], "]"}], ";", 
              RowBox[{"If", "[", 
               RowBox[{
                RowBox[{
                 RowBox[{"Norm", "[", "\[CapitalDelta]l", "]"}], "\[Equal]", 
                 "0"}], ",", 
                RowBox[{
                 RowBox[{
                  SubscriptBox["t", "\[CapitalDelta]l"], 
                  "\[LeftDoubleBracket]", 
                  RowBox[{"i", ",", "j"}], "\[RightDoubleBracket]"}], "+", 
                 RowBox[{
                  SubscriptBox["t", "\[CapitalDelta]l"], 
                  "\[LeftDoubleBracket]", 
                  RowBox[{"j", ",", "i"}], "\[RightDoubleBracket]"}]}], ",", 
                RowBox[{"If", "[", 
                 RowBox[{
                  RowBox[{
                   RowBox[{"Norm", "[", "\[CapitalDelta]l", "]"}], 
                   "\[LessEqual]", 
                   SqrtBox["2"]}], ",", 
                  RowBox[{
                   SubscriptBox["t", "\[CapitalDelta]l"], 
                   "\[LeftDoubleBracket]", 
                   RowBox[{"i", ",", "j"}], "\[RightDoubleBracket]"}], ",", 
                  "0"}], "]"}]}], "]"}]}], ")"}], "&"}], ",", 
           SubscriptBox["orb", "ind"], ",", 
           SubscriptBox["orb", "ind"], ",", "1"}], "]"}]}], "-", 
        RowBox[{"DiagonalMatrix", "[", 
         RowBox[{
          RowBox[{
           RowBox[{
            SubscriptBox["\[Mu]", "val"], "-", 
            RowBox[{
             SubscriptBox["\[Epsilon]", "val"], "\[LeftDoubleBracket]", 
             RowBox[{
              RowBox[{"First", "[", "#", "]"}], "+", "1"}], 
             "\[RightDoubleBracket]"}]}], "&"}], "/@", 
          SubscriptBox["orb", "ind"]}], "]"}]}]}], "]"}]}], ";"}], 
   "\[IndentingNewLine]", 
   RowBox[{"Dimensions", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"40", ",", "40"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["K", "val"], "\[LeftDoubleBracket]", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"1", ",", "2", ",", "3", ",", 
       RowBox[{"-", "1"}]}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"1", ",", "2", ",", "3", ",", 
       RowBox[{"-", "1"}]}], "}"}]}], "\[RightDoubleBracket]"}], "//", 
   "MatrixForm"}]}]], "Input"],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {
      RowBox[{"-", 
       FractionBox["32", "35"]}], 
      RowBox[{"-", 
       FractionBox["4", "9"]}], "0", 
      RowBox[{"-", 
       FractionBox["3", "19"]}]},
     {
      RowBox[{"-", 
       FractionBox["4", "9"]}], 
      RowBox[{"-", 
       FractionBox["32", "35"]}], 
      RowBox[{"-", 
       FractionBox["4", "9"]}], "0"},
     {"0", 
      RowBox[{"-", 
       FractionBox["4", "9"]}], 
      RowBox[{"-", 
       FractionBox["32", "35"]}], 
      FractionBox["1", "9"]},
     {
      RowBox[{"-", 
       FractionBox["3", "19"]}], "0", 
      FractionBox["1", "9"], 
      RowBox[{"-", 
       FractionBox["41", "77"]}]}
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
   RowBox[{"check", ":", " ", "symmetric"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["K", "val"], ";"}], "\[IndentingNewLine]", 
   RowBox[{"Norm", "[", 
    RowBox[{"%", "-", 
     RowBox[{"Transpose", "[", "%", "]"}]}], "]"}]}]}]], "Input"],

Cell[BoxData["0"], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   SuperscriptBox["\[ExponentialE]", 
    RowBox[{
     RowBox[{"-", "\[Tau]"}], " ", "K"}]], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["exp", "\[Tau]k"], "=", 
    RowBox[{"MatrixExp", "[", 
     RowBox[{"N", "[", 
      RowBox[{
       RowBox[{
        RowBox[{"-", 
         SubscriptBox["\[Tau]", "val"]}], 
        SubscriptBox["K", "val"]}], ",", 
       RowBox[{"4", "MachinePrecision"}]}], "]"}], "]"}]}], ";"}]}]], "Input"],

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
     {"1.13536924667298`", "0.061151949420940634`", "0.019875961813431216`"},
     {"0.061151949420940634`", "1.13536924667298`", "0.0013882490710468614`"},
     {"0.019875961813431216`", "0.0013882490710468614`", "1.0827259296636316`"}
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

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   SuperscriptBox["\[ExponentialE]", 
    RowBox[{"\[Tau]", " ", "K"}]], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["invexp", "\[Tau]k"], "=", 
    RowBox[{"MatrixExp", "[", 
     RowBox[{"N", "[", 
      RowBox[{
       RowBox[{
        SubscriptBox["\[Tau]", "val"], 
        SubscriptBox["K", "val"]}], ",", 
       RowBox[{"4", "MachinePrecision"}]}], "]"}], "]"}]}], ";"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "check", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{"Norm", "[", 
   RowBox[{
    RowBox[{
     SubscriptBox["exp", "\[Tau]k"], ".", 
     SubscriptBox["invexp", "\[Tau]k"]}], "-", 
    RowBox[{"IdentityMatrix", "[", 
     RowBox[{
      SubscriptBox["n", "orb"], 
      SubscriptBox["n", "sites"]}], "]"}]}], "]"}]}]], "Input"],

Cell[BoxData["0``63.47147700905429"], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
   "save", " ", "hopping", " ", "parameters", " ", "to", " ", "disk"}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{"Export", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["fn", "base"], "<>", "\"\<_taa.dat\>\""}], ",", 
      RowBox[{"Flatten", "[", 
       RowBox[{"Transpose", "[", 
        SubscriptBox["t", 
         RowBox[{"{", 
          RowBox[{"0", ",", "0"}], "}"}]], "]"}], "]"}], ",", 
      "\"\<Real64\>\""}], "]"}], ";"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"Export", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["fn", "base"], "<>", "\"\<_tab.dat\>\""}], ",", 
      RowBox[{"Flatten", "[", 
       RowBox[{"Transpose", "[", 
        SubscriptBox["t", 
         RowBox[{"{", 
          RowBox[{"1", ",", "0"}], "}"}]], "]"}], "]"}], ",", 
      "\"\<Real64\>\""}], "]"}], ";"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"Export", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["fn", "base"], "<>", "\"\<_tac.dat\>\""}], ",", 
      RowBox[{"Flatten", "[", 
       RowBox[{"Transpose", "[", 
        SubscriptBox["t", 
         RowBox[{"{", 
          RowBox[{"0", ",", "1"}], "}"}]], "]"}], "]"}], ",", 
      "\"\<Real64\>\""}], "]"}], ";"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"Export", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["fn", "base"], "<>", "\"\<_tad.dat\>\""}], ",", 
      RowBox[{"Flatten", "[", 
       RowBox[{"Transpose", "[", 
        SubscriptBox["t", 
         RowBox[{"{", 
          RowBox[{"1", ",", "1"}], "}"}]], "]"}], "]"}], ",", 
      "\"\<Real64\>\""}], "]"}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"(*", " ", 
    RowBox[{
     RowBox[{"hopping", " ", "b"}], " ", "\[Rule]", " ", 
     RowBox[{
     "c", " ", "is", " ", "actually", " ", "transposed", " ", "version", " ", 
      "of", " ", 
      SubscriptBox["t", 
       RowBox[{"{", 
        RowBox[{"1", ",", 
         RowBox[{"-", "1"}]}], "}"}]]}]}], " ", "*)"}], 
   RowBox[{
    RowBox[{"Export", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["fn", "base"], "<>", "\"\<_tbc.dat\>\""}], ",", 
      RowBox[{"Flatten", "[", 
       SubscriptBox["t", 
        RowBox[{"{", 
         RowBox[{"1", ",", 
          RowBox[{"-", "1"}]}], "}"}]], "]"}], ",", "\"\<Real64\>\""}], "]"}],
     ";"}]}]}]], "Input"]
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
         RowBox[{
          SubscriptBox["n", "orb"], 
          SubscriptBox["n", "sites"]}], ",", 
         SubscriptBox["L", "val"]}], "}"}]}], "]"}]}], "-", "1"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"Dimensions", "[", 
  SubscriptBox["s", "0"], "]"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"40", ",", "16"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["s", "0"], "\[LeftDoubleBracket]", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"1", ",", "2", ",", "3", ",", 
       RowBox[{"-", "1"}]}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"1", ",", "2", ",", "3", ",", 
       RowBox[{"-", "1"}]}], "}"}]}], "\[RightDoubleBracket]"}], "//", 
   "MatrixForm"}]}]], "Input"],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {
      RowBox[{"-", "1"}], "1", "1", 
      RowBox[{"-", "1"}]},
     {
      RowBox[{"-", "1"}], 
      RowBox[{"-", "1"}], "1", 
      RowBox[{"-", "1"}]},
     {"1", "1", 
      RowBox[{"-", "1"}], "1"},
     {"1", 
      RowBox[{"-", "1"}], "1", 
      RowBox[{"-", "1"}]}
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

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"save", " ", "initial", " ", "Hubbard"}], "-", 
    RowBox[{"Stratonovich", " ", "field", " ", "to", " ", "disk"}]}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"Export", "[", 
    RowBox[{
     RowBox[{
      SubscriptBox["fn", "base"], "<>", "\"\<_HS0.dat\>\""}], ",", 
     RowBox[{
      RowBox[{
       RowBox[{"Flatten", "[", 
        RowBox[{"Transpose", "[", 
         SubscriptBox["s", "0"], "]"}], "]"}], "/.", 
       RowBox[{"{", 
        RowBox[{"1", "\[Rule]", "0"}], "}"}]}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"-", "1"}], "\[Rule]", "1"}], "}"}]}], ",", 
     "\"\<Integer8\>\""}], "]"}], ";"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"\[Lambda]", " ", "parameters", " ", "for", " ", "Hubbard"}], "-", 
    RowBox[{"Stratonovich", " ", "field"}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["\[Lambda]", "val"], "=", 
     RowBox[{"ArcCosh", "[", 
      RowBox[{"Exp", "[", 
       RowBox[{
        SubscriptBox["\[Tau]", "val"], 
        RowBox[{
         SubscriptBox["U", "val"], "/", "2"}]}], "]"}], "]"}]}], ";"}], 
   "\[IndentingNewLine]", 
   RowBox[{"N", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"0.7856003600934582`", ",", "0.7031327390810828`"}], 
  "}"}]], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["\[Lambda]", "orb"], "=", 
   RowBox[{
    RowBox[{
     RowBox[{
      SubscriptBox["\[Lambda]", "val"], "\[LeftDoubleBracket]", 
      RowBox[{
       RowBox[{"First", "[", "#", "]"}], "+", "1"}], 
      "\[RightDoubleBracket]"}], "&"}], "/@", 
    SubscriptBox["orb", "ind"]}]}], ";"}]], "Input"]
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
   "expK_", ",", "invexpK_", ",", "\[Lambda]_List", ",", "s0_List", ",", 
    "orbcellorder_List", ",", "rnd_List"}], "]"}], ":=", 
  RowBox[{"Module", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"s", "=", "s0"}], ",", "Gu", ",", "Gd", ",", "norbcell", ",", 
      "L", ",", "du", ",", "dd", ",", "cu", ",", "cd", ",", "l", ",", "i", 
      ",", "j"}], "}"}], ",", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{
      RowBox[{"{", 
       RowBox[{"norbcell", ",", "L"}], "}"}], "=", 
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
        RowBox[{
         RowBox[{"DiagonalMatrix", "[", "\[Lambda]", "]"}], ".", "s"}]}], 
       "]"}]}], ";", "\[IndentingNewLine]", 
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
         RowBox[{
         "iterate", " ", "over", " ", "lattice", " ", "cells", " ", "and", 
          " ", "orbitals"}], " ", "*)"}], "\[IndentingNewLine]", 
        RowBox[{"For", "[", 
         RowBox[{
          RowBox[{"j", "=", "1"}], ",", 
          RowBox[{"j", "\[LessEqual]", "norbcell"}], ",", 
          RowBox[{"j", "++"}], ",", "\[IndentingNewLine]", 
          RowBox[{
           RowBox[{"i", "=", 
            RowBox[{"orbcellorder", "\[LeftDoubleBracket]", 
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
                  RowBox[{"+", "2"}], 
                  RowBox[{
                  "\[Lambda]", "\[LeftDoubleBracket]", "i", 
                   "\[RightDoubleBracket]"}], 
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
                  RowBox[{"-", "2"}], 
                  RowBox[{
                  "\[Lambda]", "\[LeftDoubleBracket]", "i", 
                   "\[RightDoubleBracket]"}], 
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
                    RowBox[{"+", "2"}], 
                    RowBox[{
                    "\[Lambda]", "\[LeftDoubleBracket]", "i", 
                    "\[RightDoubleBracket]"}], 
                    RowBox[{"s", "\[LeftDoubleBracket]", 
                    RowBox[{"i", ",", "l"}], "\[RightDoubleBracket]"}]}], 
                   "]"}], "-", "1"}], ")"}], 
                RowBox[{"(", 
                 RowBox[{
                  RowBox[{"UnitVector", "[", 
                   RowBox[{"norbcell", ",", "i"}], "]"}], "-", 
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
                    RowBox[{"-", "2"}], 
                    RowBox[{
                    "\[Lambda]", "\[LeftDoubleBracket]", "i", 
                    "\[RightDoubleBracket]"}], 
                    RowBox[{"s", "\[LeftDoubleBracket]", 
                    RowBox[{"i", ",", "l"}], "\[RightDoubleBracket]"}]}], 
                   "]"}], "-", "1"}], ")"}], 
                RowBox[{"(", 
                 RowBox[{
                  RowBox[{"UnitVector", "[", 
                   RowBox[{"norbcell", ",", "i"}], "]"}], "-", 
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
       RowBox[{
        SubscriptBox["n", "orb"], 
        SubscriptBox["n", "sites"]}]}], "]"}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"Dimensions", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"16", ",", "40"}], "}"}]], "Output"]
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
    "0.95463862700127`", ",", "0.30392075394759765`", ",", 
     "0.8844326011924695`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
    "0.6005706127827714`", ",", "0.8561780715123695`", ",", 
     "0.825407882544888`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
    "0.04729107837330126`", ",", "0.2834439599668278`", ",", 
     "0.6485657532654547`"}], "}"}]}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"pre", "-", 
    RowBox[{"determined", " ", "update", " ", "ordering"}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["orbcell", "order"], "=", 
     RowBox[{"Partition", "[", 
      RowBox[{
       RowBox[{
        RowBox[{"Import", "[", 
         RowBox[{
          RowBox[{
           SubscriptBox["fn", "base"], "<>", "\"\<_orbcellorder.dat\>\""}], 
          ",", "\"\<Integer32\>\""}], "]"}], "+", "1"}], ",", 
       RowBox[{
        SubscriptBox["n", "orb"], 
        SubscriptBox["n", "sites"]}]}], "]"}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"Dimensions", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"16", ",", "40"}], "}"}]], "Output"]
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
        RowBox[{
         SubscriptBox["n", "orb"], 
         SubscriptBox["n", "sites"]}], "]"}]}], "]"}], "&"}], "/@", 
    SubscriptBox["orbcell", "order"]}], "]"}]}]], "Input"],

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
         SubscriptBox["\[Lambda]", "orb"], ",", 
         RowBox[{"4", "MachinePrecision"}]}], "]"}], ",", 
       SubscriptBox["s", "0"], ",", 
       SubscriptBox["orbcell", "order"], ",", 
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
        RowBox[{"DiagonalMatrix", "[", 
         RowBox[{"N", "[", 
          RowBox[{
           SubscriptBox["\[Lambda]", "orb"], ",", 
           RowBox[{"4", "MachinePrecision"}]}], "]"}], "]"}], ".", 
        SubscriptBox["s", "1"]}]}], "]"}]}], "]"}], "]"}]}]], "Input"],

Cell[BoxData["0"], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["s", "1"], "\[LeftDoubleBracket]", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"1", ",", "2", ",", "3", ",", 
       RowBox[{"-", "1"}]}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"1", ",", "2", ",", "3", ",", 
       RowBox[{"-", "1"}]}], "}"}]}], "\[RightDoubleBracket]"}], "//", 
   "MatrixForm"}]}]], "Input"],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"1", 
      RowBox[{"-", "1"}], 
      RowBox[{"-", "1"}], 
      RowBox[{"-", "1"}]},
     {"1", 
      RowBox[{"-", "1"}], 
      RowBox[{"-", "1"}], 
      RowBox[{"-", "1"}]},
     {
      RowBox[{"-", "1"}], 
      RowBox[{"-", "1"}], 
      RowBox[{"-", "1"}], 
      RowBox[{"-", "1"}]},
     {
      RowBox[{"-", "1"}], 
      RowBox[{"-", "1"}], 
      RowBox[{"-", "1"}], 
      RowBox[{"-", "1"}]}
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

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"save", " ", "final", " ", "Hubbard"}], "-", 
    RowBox[{"Stratonovich", " ", "field", " ", "to", " ", "disk"}]}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"Export", "[", 
    RowBox[{
     RowBox[{
      SubscriptBox["fn", "base"], "<>", "\"\<_HS1.dat\>\""}], ",", 
     RowBox[{
      RowBox[{
       RowBox[{"Flatten", "[", 
        RowBox[{"Transpose", "[", 
         SubscriptBox["s", "1"], "]"}], "]"}], "/.", 
       RowBox[{"{", 
        RowBox[{"1", "\[Rule]", "0"}], "}"}]}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"-", "1"}], "\[Rule]", "1"}], "}"}]}], ",", 
     "\"\<Integer8\>\""}], "]"}], ";"}]}]], "Input"],

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
     {"0.010515716802566114651037804219461530822167372173673801004648583202026\
49853511786046567227122333`63.81835908076401", 
      RowBox[{
      "-", "0.0064038748123094914319003509649727047852914405985774687416482845\
40951075444053010214724206925246`63.81835908076401"}], 
      RowBox[{
      "-", "0.0084116287136778602207038269586555866829261240883001029461541563\
75881631341977512976899762780485`63.81835908076401"}]},
     {
      RowBox[{
      "-", "0.0137513963961349524475472249878568076615952942188384717579551414\
82454915966685134923167698785506`63.81835908076401"}], 
      "0.067515010655889973649189523842986194402852453573421523971361012374546\
66303821196490591456536222`63.81835908076401", 
      RowBox[{
      "-", "0.0010522403538976021456491858590636854533694672318649857050367585\
44088226923861085557936879003182`63.81835908076402"}]},
     {
      RowBox[{
      "-", "0.0014418457471853347445312901261859348218628049026794455067889258\
11995022111756496595492830289908`63.81835908076401"}], 
      RowBox[{
      "-", "0.0000418921815958411183608580072087076331231441358823038967498693\
5941795852257393276800235497139`63.81835908076401"}], 
      "0.010601324314361003042343288211957959344799674143316397689803053535209\
216025565844472875519647765`63.81835908076401"}
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
     {"0.939457384668491111238653116267780296012619696739650077372890868562371\
868287110702434435136158925`63.81835908076401", 
      RowBox[{
      "-", "0.0698056777086087478555293789916801177702108374093611074844585216\
77258848112408091567732401098985`63.81835908076401"}], 
      "0.006740125424451534816199145844557583101915967261026597298499803239553\
492835090418745637421226996`63.81835908076401"},
     {
      RowBox[{
      "-", "0.0406359703252602574877335579699195516036347614098507632460554319\
95711140899176257982703168630444`63.81835908076401"}], 
      "0.795534355346729299924923751892424808229763357888462857695583311266324\
379570151228880840005757153`63.81835908076401", 
      RowBox[{
      "-", "0.0018454347599997979521679115785216839734551257327329473081946266\
50890489772840415672368573489788`63.81835908076401"}]},
     {"0.005409872615239048280151766200268404443151562613959706276681450730313\
870770067654951727529218535`63.81835908076401", 
      RowBox[{
      "-", "0.0126839086310488797124197780826691432979627916967950020579389883\
28613011287326427977861202461315`63.81835908076401"}], 
      "0.978698695181088324543742404531592165955050112641962871114066977374380\
028733459560536363770068009`63.81835908076401"}
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

Cell[BoxData["0.9811413140287052`"], "Output"],

Cell[BoxData["0.9814977131811229`"], "Output"]
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

Cell[BoxData["1.\
453670651197544573738083756036978339739825845676523929635482250048303444071882\
71017211`63.81835908076401*^-86"], "Output"],

Cell[BoxData["7.\
510646143973916956800775469233351706543130664251842693789989680232728758784029\
01418`63.81835908076401*^-52"], "Output"]
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
WindowMargins->{{Automatic, 168}, {66, Automatic}},
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
Cell[1529, 37, 264, 8, 52, "Input"],
Cell[CellGroupData[{
Cell[1818, 49, 411, 13, 72, "Input"],
Cell[2232, 64, 94, 2, 31, "Output"]
}, Open  ]],
Cell[2341, 69, 219, 7, 52, "Input"],
Cell[2563, 78, 261, 8, 67, "Input"],
Cell[CellGroupData[{
Cell[2849, 90, 295, 9, 69, "Input"],
Cell[3147, 101, 29, 0, 31, "Output"]
}, Open  ]],
Cell[3191, 104, 242, 7, 31, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[3470, 116, 62, 0, 43, "Subsection"],
Cell[3535, 118, 320, 10, 72, "Input"],
Cell[3858, 130, 161, 6, 31, "Input"],
Cell[CellGroupData[{
Cell[4044, 140, 1259, 38, 72, "Input"],
Cell[5306, 180, 2813, 82, 72, "Output"],
Cell[8122, 264, 74, 2, 31, "Output"]
}, Open  ]],
Cell[8211, 269, 2400, 87, 211, "Input"],
Cell[10614, 358, 236, 7, 67, "Input"],
Cell[10853, 367, 387, 12, 52, "Input"],
Cell[CellGroupData[{
Cell[11265, 383, 4470, 110, 188, "Input"],
Cell[15738, 495, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[15850, 502, 456, 13, 52, "Input"],
Cell[16309, 517, 1199, 41, 126, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[17545, 563, 328, 10, 72, "Input"],
Cell[17876, 575, 28, 0, 31, "Output"]
}, Open  ]],
Cell[17919, 578, 513, 16, 52, "Input"],
Cell[CellGroupData[{
Cell[18457, 598, 480, 14, 52, "Input"],
Cell[18940, 614, 789, 18, 71, "Output"]
}, Open  ]],
Cell[19744, 635, 470, 14, 52, "Input"],
Cell[CellGroupData[{
Cell[20239, 653, 377, 11, 52, "Input"],
Cell[20619, 666, 47, 0, 31, "Output"]
}, Open  ]],
Cell[20681, 669, 2386, 71, 152, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[23104, 745, 56, 0, 43, "Subsection"],
Cell[CellGroupData[{
Cell[23185, 749, 614, 20, 72, "Input"],
Cell[23802, 771, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[23914, 778, 454, 13, 52, "Input"],
Cell[24371, 793, 835, 27, 86, "Output"]
}, Open  ]],
Cell[25221, 823, 737, 22, 52, "Input"],
Cell[CellGroupData[{
Cell[25983, 849, 581, 18, 72, "Input"],
Cell[26567, 869, 112, 3, 31, "Output"]
}, Open  ]],
Cell[26694, 875, 350, 11, 31, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[27081, 891, 31, 0, 43, "Subsection"],
Cell[27115, 893, 625, 17, 52, "Input"],
Cell[27743, 912, 843, 22, 52, "Input"],
Cell[28589, 936, 513, 13, 72, "Input"],
Cell[29105, 951, 9929, 229, 532, "Input"],
Cell[CellGroupData[{
Cell[39059, 1184, 819, 24, 72, "Input"],
Cell[39881, 1210, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[39993, 1217, 393, 11, 52, "Input"],
Cell[40389, 1230, 454, 14, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[40880, 1249, 711, 21, 72, "Input"],
Cell[41594, 1272, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[41706, 1279, 439, 13, 52, "Input"],
Cell[42148, 1294, 28, 0, 31, "Output"]
}, Open  ]],
Cell[42191, 1297, 838, 25, 31, "Input"],
Cell[CellGroupData[{
Cell[43054, 1326, 923, 26, 52, "Input"],
Cell[43980, 1354, 28, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[44045, 1359, 454, 13, 52, "Input"],
Cell[44502, 1374, 967, 33, 86, "Output"]
}, Open  ]],
Cell[45484, 1410, 735, 22, 52, "Input"],
Cell[CellGroupData[{
Cell[46244, 1436, 877, 27, 72, "Input"],
Cell[47124, 1465, 1855, 41, 71, "Output"],
Cell[48982, 1508, 1807, 38, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[50826, 1551, 559, 17, 72, "Input"],
Cell[51388, 1570, 46, 0, 31, "Output"],
Cell[51437, 1572, 46, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[51520, 1577, 476, 15, 72, "Input"],
Cell[51999, 1594, 142, 2, 31, "Output"],
Cell[52144, 1598, 139, 2, 31, "Output"]
}, Open  ]],
Cell[52298, 1603, 1513, 43, 72, "Input"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature Ex0ABso#toXZZA1FRUogRM2W *)
