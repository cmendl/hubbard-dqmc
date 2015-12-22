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
NotebookDataLength[     50470,       1654]
NotebookOptionsPosition[     47297,       1528]
NotebookOutlinePosition[     47640,       1543]
CellTagsIndexPosition[     47597,       1540]
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
     SubscriptBox["n", "y"], "=", "6"}], ";"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"total", " ", "number", " ", "of", " ", "unit", " ", "cells"}], 
   " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["n", "cells"], "=", 
    RowBox[{
     SubscriptBox["n", "x"], 
     SubscriptBox["n", "y"]}]}], ";"}]}]], "Input"],

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
    RowBox[{"0", ",", "0", ",", "5"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "1", ",", "5"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "2", ",", "5"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "3", ",", "5"}], "}"}], ",", 
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
    RowBox[{"1", ",", "3", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "0", ",", "5"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "1", ",", "5"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "2", ",", "5"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "3", ",", "5"}], "}"}]}], "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"48", ",", "3"}], "}"}]], "Output"]
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
         FractionBox["2", "9"]}], "}"}], ",", 
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
         FractionBox["4", "3"], ",", 
         RowBox[{"-", 
          FractionBox["3", "7"]}]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         FractionBox["2", "9"], ",", 
         FractionBox["15", "16"]}], "}"}]}], "}"}]}], ";"}], 
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
         FractionBox["1", "7"], ",", 
         RowBox[{"-", 
          FractionBox["3", "11"]}]}], "}"}]}], "}"}]}], ";"}], 
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
         RowBox[{"-", 
          FractionBox["1", "9"]}], ",", 
         FractionBox["4", "5"]}], "}"}], ",", 
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
       RowBox[{"1", ",", 
        RowBox[{"-", "1"}]}], "}"}]], "=", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{
         FractionBox["3", "19"], ",", 
         FractionBox["1", "6"]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         FractionBox["1", "5"], ",", 
         FractionBox["1", "13"]}], "}"}]}], "}"}]}], ";"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"chemical", " ", "potential"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[Mu]", "val"], "=", 
    FractionBox["2", "7"]}], ";"}]}]], "Input"],

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
      RowBox[{"1", "/", "5"}], ",", 
      RowBox[{
       RowBox[{"-", "3"}], "/", "11"}]}], "}"}]}], ";"}]}]], "Input"],

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
  RowBox[{"48", ",", "48"}], "}"}]], "Output"]
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
       FractionBox["3", "35"]}], 
      RowBox[{"-", 
       FractionBox["4", "3"]}], "0", 
      RowBox[{"-", 
       FractionBox["3", "7"]}]},
     {
      RowBox[{"-", 
       FractionBox["4", "3"]}], 
      RowBox[{"-", 
       FractionBox["3", "35"]}], 
      RowBox[{"-", 
       FractionBox["4", "3"]}], "0"},
     {"0", 
      RowBox[{"-", 
       FractionBox["4", "3"]}], 
      RowBox[{"-", 
       FractionBox["3", "35"]}], 
      RowBox[{"-", 
       FractionBox["1", "6"]}]},
     {
      RowBox[{"-", 
       FractionBox["3", "7"]}], "0", 
      RowBox[{"-", 
       FractionBox["1", "6"]}], 
      RowBox[{"-", 
       FractionBox["43", "77"]}]}
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
     {"1.0502800096504994`", "0.1746622508398583`", "0.05932628151426754`"},
     {"0.1746622508398583`", "1.0502800096504994`", "0.010530829456782832`"},
     {"0.05932628151426754`", "0.010530829456782832`", "1.099246504517032`"}
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
      SubscriptBox["n", "cells"]}], "]"}]}], "]"}]}]], "Input"],

Cell[BoxData["0``63.38967201398729"], "Output"]
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
             SubscriptBox["fn", "base"], "<>", "\"\<_s.dat\>\""}], ",", 
            "\"\<Integer8\>\""}], "]"}], ",", 
          RowBox[{
           SubscriptBox["n", "orb"], 
           SubscriptBox["n", "cells"]}]}], "]"}], "]"}], "/.", 
       RowBox[{"{", 
        RowBox[{"1", "\[Rule]", 
         RowBox[{"-", "1"}]}], "}"}]}], "/.", 
      RowBox[{"{", 
       RowBox[{"0", "\[Rule]", "1"}], "}"}]}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"Dimensions", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"48", ",", "16"}], "}"}]], "Output"]
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
    RowBox[{
     RowBox[{"-", "1"}], ",", "1", ",", 
     RowBox[{"-", "1"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", 
     RowBox[{"-", "1"}], ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "1"}], ",", "1", ",", 
     RowBox[{"-", "1"}]}], "}"}]}], "}"}]], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"Coulomb", " ", "coupling", " ", "constants"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["U", "val"], "=", 
    RowBox[{"{", 
     RowBox[{
      RowBox[{"33", "/", "8"}], ",", 
      RowBox[{"17", "/", "5"}]}], "}"}]}], ";"}]}]], "Input"],

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
  RowBox[{"0.74928560110877`", ",", "0.6752355951860223`"}], "}"}]], "Output"]
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

Cell["Phonon block update", "Subsection"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"phonon", " ", "frequencies"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["\[CapitalOmega]", "val"], "=", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"6", "/", "5"}], ",", 
       RowBox[{"5", "/", "6"}]}], "}"}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"N", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"1.2`", ",", "0.8333333333333334`"}], "}"}]], "Output"]
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
     RowBox[{"{", 
      RowBox[{
       RowBox[{"13", "/", "20"}], ",", 
       RowBox[{"19", "/", "17"}]}], "}"}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"N", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"0.65`", ",", "1.1176470588235294`"}], "}"}]], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["\[CapitalOmega]", "orb"], "=", 
   RowBox[{
    RowBox[{
     RowBox[{
      SubscriptBox["\[CapitalOmega]", "val"], "\[LeftDoubleBracket]", 
      RowBox[{
       RowBox[{"First", "[", "#", "]"}], "+", "1"}], 
      "\[RightDoubleBracket]"}], "&"}], "/@", 
    SubscriptBox["orb", "ind"]}]}], ";"}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["g", "orb"], "=", 
   RowBox[{
    RowBox[{
     RowBox[{
      SubscriptBox["g", "val"], "\[LeftDoubleBracket]", 
      RowBox[{
       RowBox[{"First", "[", "#", "]"}], "+", "1"}], 
      "\[RightDoubleBracket]"}], "&"}], "/@", 
    SubscriptBox["orb", "ind"]}]}], ";"}]], "Input"],

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
   "expK_", ",", "\[Tau]_", ",", "\[Lambda]_List", ",", 
    "\[CapitalOmega]_List", ",", "g_List", ",", "s_List", ",", "X0_List", ",",
     "i_Integer", ",", "\[CapitalDelta]x_"}], "]"}], ":=", 
  RowBox[{"Module", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
     "Gu0", ",", "Gd0", ",", "Gu1", ",", "Gd1", ",", "\[CapitalDelta]Eph", 
      ",", "X1"}], "}"}], ",", "\[IndentingNewLine]", 
    RowBox[{"(*", " ", 
     RowBox[{"compute", " ", "the", " ", "initial", " ", 
      RowBox[{"Green", "'"}], "s", " ", "functions"}], " ", "*)"}], 
    "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{
      RowBox[{"{", 
       RowBox[{"Gu0", ",", "Gd0"}], "}"}], "=", 
      RowBox[{"InitializeUpDownGreensFunctions", "[", 
       RowBox[{"expK", ",", 
        RowBox[{
         RowBox[{"DiagonalMatrix", "[", "\[Lambda]", "]"}], ".", "s"}], ",", 
        RowBox[{"\[Tau]", " ", 
         RowBox[{
          RowBox[{"DiagonalMatrix", "[", "g", "]"}], ".", "X0"}]}]}], "]"}]}],
      ";", "\[IndentingNewLine]", 
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
       SuperscriptBox[
        RowBox[{
        "\[CapitalOmega]", "\[LeftDoubleBracket]", "i", 
         "\[RightDoubleBracket]"}], "2"], 
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
        RowBox[{
         RowBox[{"DiagonalMatrix", "[", "\[Lambda]", "]"}], ".", "s"}], ",", 
        RowBox[{"\[Tau]", " ", 
         RowBox[{
          RowBox[{"DiagonalMatrix", "[", "g", "]"}], ".", "X1"}]}]}], "]"}]}],
      ";", "\[IndentingNewLine]", 
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
             SubscriptBox["fn", "base"], "<>", "\"\<_X0.dat\>\""}], ",", 
            "\"\<Real64\>\""}], "]"}], ",", 
          RowBox[{"4", "MachinePrecision"}]}], "]"}], ",", 
        RowBox[{
         SubscriptBox["n", "orb"], 
         SubscriptBox["n", "cells"]}]}], "]"}], "]"}]}], ";"}], 
   "\[IndentingNewLine]", 
   RowBox[{"Dimensions", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"48", ",", "16"}], "}"}]], "Output"]
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
     RowBox[{"-", "0.9639052773577816`"}], ",", "0.3002895359822837`", ",", 
     RowBox[{"-", "0.34487387435709316`"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "0.20870090118062867`"}], ",", 
     RowBox[{"-", "0.3979167113586759`"}], ",", 
     RowBox[{"-", "0.1039610880383317`"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "0.45808557712004916`"}], ",", 
     RowBox[{"-", "0.8918462880433915`"}], ",", "0.19648963966598743`"}], 
    "}"}]}], "}"}]], "Output"]
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
       RowBox[{"DiagonalMatrix", "[", 
        SubscriptBox["\[Lambda]", "orb"], "]"}], ".", 
       SubscriptBox["s", "val"]}], ",", 
      RowBox[{
       SubscriptBox["\[Tau]", "val"], 
       RowBox[{
        RowBox[{"DiagonalMatrix", "[", 
         SubscriptBox["g", "orb"], "]"}], ".", 
        SubscriptBox["X", "0"]}]}]}], "]"}]}], ";"}]}]], "Input"],

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
    RowBox[{"0.18977881827483797`", ",", 
     RowBox[{"-", "0.3454091251126187`"}], ",", 
     RowBox[{"-", "0.023393529890992908`"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "0.06704790116655956`"}], ",", "0.8007998734640727`", ",", 
     "0.01186238751070249`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
    "0.010695966744322733`", ",", "0.10609869057607492`", ",", 
     "0.07049015337219558`"}], "}"}]}], "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"0.6969430469184574`", ",", 
     RowBox[{"-", "0.07849132495641707`"}], ",", 
     RowBox[{"-", "0.05154765199519996`"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "0.4722162516509891`"}], ",", "0.2120226795821074`", ",", 
     RowBox[{"-", "0.019973757174955323`"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "0.02883642574193392`"}], ",", 
     RowBox[{"-", "0.004426062084514871`"}], ",", "0.8542813701844659`"}], 
    "}"}]}], "}"}]], "Output"]
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
      SubscriptBox["\[Lambda]", "orb"], ",", 
      SubscriptBox["\[CapitalOmega]", "orb"], ",", 
      SubscriptBox["g", "orb"], ",", 
      SubscriptBox["s", "val"], ",", 
      SubscriptBox["X", "0"], ",", "40", ",", 
      SubscriptBox["\[CapitalDelta]x", "0"]}], "]"}]}], ";"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 InterpretationBox[
  RowBox[{"\<\"Det[Gu0]Det[Gd0]: \"\>", "\[InvisibleSpace]", 
   "1.247422081197379117478462858521105617349755951431552862826703705495665832\
974049425065372`58.14125038571763*^-163"}],
  SequenceForm[
  "Det[Gu0]Det[Gd0]: ", 
   1.2474220811973791174784628585211056173497559514315528628267037054956658329\
74049425065372`58.14125038571763*^-163],
  Editable->False]], "Print"],

Cell[BoxData[
 InterpretationBox[
  RowBox[{"\<\"Det[Gu1]Det[Gd1]: \"\>", "\[InvisibleSpace]", 
   "1.010603138626472207118386940451385651739231840578531344329460697409365016\
5645937046527258`58.86794471936846*^-162"}],
  SequenceForm[
  "Det[Gu1]Det[Gd1]: ", 
   1.0106031386264722071183869404513856517392318405785313443294606974093650165\
645937046527258`58.86794471936846*^-162],
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
225366592611184175940576220069746604706967807792530115278198812497060046199436\
412560030552909023`58.06656578073824"], "Output"]
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
      RowBox[{"-", "0.9639052773577816`"}], "0.3002895359822837`", 
      RowBox[{"-", "0.34487387435709316`"}]},
     {
      RowBox[{"-", "0.20870090118062867`"}], 
      RowBox[{"-", "0.3979167113586759`"}], 
      RowBox[{"-", "0.1039610880383317`"}]},
     {
      RowBox[{"-", "0.45808557712004916`"}], 
      RowBox[{"-", "0.8918462880433915`"}], "0.19648963966598743`"}
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
     {"0.18977382352891237`", 
      RowBox[{"-", "0.3453819047873784`"}], 
      RowBox[{"-", "0.02338989588043585`"}]},
     {
      RowBox[{"-", "0.0670351322532808`"}], "0.8008021894425651`", 
      "0.011862084657657845`"},
     {"0.010710078909873253`", "0.10608316236374463`", "0.07048488514713107`"}
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
     {"0.6969079507403365`", 
      RowBox[{"-", "0.07842843555379299`"}], 
      RowBox[{"-", "0.051377625426537656`"}]},
     {
      RowBox[{"-", "0.47220484334945567`"}], "0.21194181005066046`", 
      RowBox[{"-", "0.020132200534340183`"}]},
     {
      RowBox[{"-", "0.028859875514228975`"}], 
      RowBox[{"-", "0.004444119521230867`"}], "0.8543165818311285`"}
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

Cell[BoxData["1.689998921283661`"], "Output"],

Cell[BoxData["1.420446595275789`"], "Output"]
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
492097373712942444690687564684640089122550323969536837391615069707970681698725\
69815194162303`59.727465997573596*^-80"], "Output"],

Cell[BoxData["6.\
773037446689436799886068751793013049530975413782798474565479268324884873583941\
999381086916`58.9325335277861*^-83"], "Output"]
}, Open  ]],

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
           RowBox[{
            RowBox[{"DiagonalMatrix", "[", 
             SubscriptBox["g", "orb"], "]"}], ".", 
            SubscriptBox["X", "1"]}]}], "]"}], "]"}], "]"}], "]"}], ",", 
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
WindowMargins->{{Automatic, 254}, {Automatic, 51}},
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
Cell[1796, 47, 219, 7, 52, "Input"],
Cell[2018, 56, 261, 8, 67, "Input"],
Cell[CellGroupData[{
Cell[2304, 68, 295, 9, 69, "Input"],
Cell[2602, 79, 29, 0, 31, "Output"]
}, Open  ]],
Cell[2646, 82, 242, 7, 31, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[2925, 94, 62, 0, 43, "Subsection"],
Cell[2990, 96, 320, 10, 72, "Input"],
Cell[3313, 108, 317, 10, 52, "Input"],
Cell[CellGroupData[{
Cell[3655, 122, 1259, 38, 72, "Input"],
Cell[4917, 162, 3365, 98, 72, "Output"],
Cell[8285, 262, 74, 2, 31, "Output"]
}, Open  ]],
Cell[8374, 267, 2375, 86, 211, "Input"],
Cell[10752, 355, 236, 7, 67, "Input"],
Cell[10991, 364, 387, 12, 52, "Input"],
Cell[CellGroupData[{
Cell[11403, 380, 4470, 110, 188, "Input"],
Cell[15876, 492, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[15988, 499, 456, 13, 52, "Input"],
Cell[16447, 514, 1240, 43, 126, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[17724, 562, 328, 10, 72, "Input"],
Cell[18055, 574, 28, 0, 31, "Output"]
}, Open  ]],
Cell[18098, 577, 513, 16, 52, "Input"],
Cell[CellGroupData[{
Cell[18636, 597, 480, 14, 52, "Input"],
Cell[19119, 613, 784, 18, 71, "Output"]
}, Open  ]],
Cell[19918, 634, 470, 14, 52, "Input"],
Cell[CellGroupData[{
Cell[20413, 652, 377, 11, 52, "Input"],
Cell[20793, 665, 47, 0, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[20889, 671, 48, 0, 43, "Subsection"],
Cell[CellGroupData[{
Cell[20962, 675, 1025, 30, 72, "Input"],
Cell[21990, 707, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[22102, 714, 390, 11, 52, "Input"],
Cell[22495, 727, 365, 13, 31, "Output"]
}, Open  ]],
Cell[22875, 743, 335, 11, 52, "Input"],
Cell[CellGroupData[{
Cell[23235, 758, 581, 18, 72, "Input"],
Cell[23819, 778, 107, 2, 31, "Output"]
}, Open  ]],
Cell[23941, 783, 350, 11, 31, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[24328, 799, 41, 0, 43, "Subsection"],
Cell[CellGroupData[{
Cell[24394, 803, 405, 12, 72, "Input"],
Cell[24802, 817, 94, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[24933, 824, 447, 14, 72, "Input"],
Cell[25383, 840, 95, 2, 31, "Output"]
}, Open  ]],
Cell[25493, 845, 362, 11, 31, "Input"],
Cell[25858, 858, 334, 11, 31, "Input"],
Cell[26195, 871, 625, 17, 52, "Input"],
Cell[26823, 890, 843, 22, 52, "Input"],
Cell[27669, 914, 612, 17, 72, "Input"],
Cell[28284, 933, 3378, 88, 253, "Input"],
Cell[CellGroupData[{
Cell[31687, 1025, 801, 24, 72, "Input"],
Cell[32491, 1051, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[32603, 1058, 429, 13, 52, "Input"],
Cell[33035, 1073, 578, 16, 31, "Output"]
}, Open  ]],
Cell[33628, 1092, 860, 26, 52, "Input"],
Cell[CellGroupData[{
Cell[34513, 1122, 866, 27, 72, "Input"],
Cell[35382, 1151, 511, 14, 31, "Output"],
Cell[35896, 1167, 560, 15, 31, "Output"]
}, Open  ]],
Cell[36471, 1185, 224, 6, 31, "Input"],
Cell[CellGroupData[{
Cell[36720, 1195, 884, 25, 52, "Input"],
Cell[CellGroupData[{
Cell[37629, 1224, 411, 9, 23, "Print"],
Cell[38043, 1235, 413, 9, 23, "Print"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[38505, 1250, 172, 5, 52, "Input"],
Cell[38680, 1257, 147, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[38864, 1264, 472, 14, 52, "Input"],
Cell[39339, 1280, 937, 25, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[40313, 1310, 957, 29, 72, "Input"],
Cell[41273, 1341, 860, 22, 71, "Output"],
Cell[42136, 1365, 922, 24, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[43095, 1394, 559, 17, 72, "Input"],
Cell[43657, 1413, 45, 0, 31, "Output"],
Cell[43705, 1415, 45, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[43787, 1420, 476, 15, 72, "Input"],
Cell[44266, 1437, 149, 2, 31, "Output"],
Cell[44418, 1441, 145, 2, 31, "Output"]
}, Open  ]],
Cell[44578, 1446, 1217, 34, 72, "Input"],
Cell[45798, 1482, 1483, 43, 72, "Input"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature sxpy9j6wTCZ1vDgfhVe5Srb3 *)
