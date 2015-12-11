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
NotebookDataLength[     74580,       2272]
NotebookOptionsPosition[     70891,       2130]
NotebookOutlinePosition[     71234,       2145]
CellTagsIndexPosition[     71191,       2142]
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
     SubscriptBox["n", "y"], "=", "5"}], ";"}]}]}]], "Input"],

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
         FractionBox["4", "3"]}], "}"}], ",", 
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
         RowBox[{"-", 
          FractionBox["3", "7"]}]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         FractionBox["2", "9"], ",", 
         FractionBox["3", "16"]}], "}"}]}], "}"}]}], ";"}], 
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
         FractionBox["1", "13"], ",", 
         FractionBox["1", "5"]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         FractionBox["3", "7"], ",", 
         RowBox[{"-", 
          FractionBox["2", "19"]}]}], "}"}]}], "}"}]}], ";"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"chemical", " ", "potential"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[Mu]", "val"], "=", 
    FractionBox["6", "7"]}], ";"}]}]], "Input"],

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
      RowBox[{"3", "/", "11"}]}], "}"}]}], ";"}]}]], "Input"],

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
       FractionBox["37", "35"]}], 
      RowBox[{"-", 
       FractionBox["4", "9"]}], "0", 
      RowBox[{"-", 
       FractionBox["3", "19"]}]},
     {
      RowBox[{"-", 
       FractionBox["4", "9"]}], 
      RowBox[{"-", 
       FractionBox["37", "35"]}], 
      RowBox[{"-", 
       FractionBox["4", "9"]}], "0"},
     {"0", 
      RowBox[{"-", 
       FractionBox["4", "9"]}], 
      RowBox[{"-", 
       FractionBox["37", "35"]}], 
      RowBox[{"-", 
       FractionBox["1", "5"]}]},
     {
      RowBox[{"-", 
       FractionBox["3", "19"]}], "0", 
      RowBox[{"-", 
       FractionBox["1", "5"]}], 
      RowBox[{"-", 
       FractionBox["45", "77"]}]}
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
     {"1.1710577459066396`", "0.06393066948421308`", "0.02342921265687421`"},
     {"0.06393066948421308`", "1.1710577459066396`", "0.0022983055520301348`"},
     {"0.02342921265687421`", "0.0022983055520301348`", "1.103650818488647`"}
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

Cell[BoxData["0``63.429583732914985"], "Output"]
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

Cell["Hubbard-Stratonovich field", "Subsection"],

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
          SubscriptBox["n", "cells"]}], ",", 
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
       RowBox[{"10", "/", "3"}], ",", 
       RowBox[{"11", "/", "5"}]}], "}"}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"N", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"3.3333333333333335`", ",", "2.2`"}], "}"}]], "Output"]
}, Open  ]],

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
  RowBox[{"0.6681262480464829`", ",", "0.5365005605374655`"}], 
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

Cell["Phonon field", "Subsection"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"SeedRandom", "[", "42", "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   SubscriptBox["X", "0"], "=", 
   RowBox[{"RandomReal", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
       RowBox[{"-", "4"}], ",", "4"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{
        SubscriptBox["n", "orb"], 
        SubscriptBox["n", "cells"]}], ",", 
       SubscriptBox["L", "val"]}], "}"}], ",", 
     RowBox[{"WorkingPrecision", "\[Rule]", 
      RowBox[{"4", "MachinePrecision"}]}]}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"Dimensions", "[", "%", "]"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"40", ",", "16"}], "}"}]], "Output"]
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
     RowBox[{"-", "0.5927577302722995`"}], ",", "3.163307780879373`", ",", 
     "1.2375801082319826`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
    "2.4395359773435485`", ",", "2.4936758164548483`", ",", 
     "1.782206609362699`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "1.993183461948001`"}], ",", "3.181277987312156`", ",", 
     RowBox[{"-", "0.6825705230247102`"}]}], "}"}]}], "}"}]], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"terms", " ", "in", " ", "the", " ", "phonon", " ", 
    RowBox[{"(", "lattice", ")"}], " ", "energy", " ", "depending", " ", "on",
     " ", 
    SubscriptBox["X", 
     RowBox[{"i", ",", "l"}]]}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"Eph", "[", "Xil_", "]"}], ":=", 
   RowBox[{
    RowBox[{
     FractionBox["1", "2"], 
     SuperscriptBox["\[CapitalOmega]", "2"], 
     SuperscriptBox["Xil", "2"]}], "+", 
    RowBox[{
     FractionBox["1", "2"], 
     SuperscriptBox[
      RowBox[{"(", 
       FractionBox[
        RowBox[{
         SubscriptBox["X", 
          RowBox[{"i", ",", 
           RowBox[{"l", "+", "1"}]}]], "-", "Xil"}], "\[CapitalDelta]\[Tau]"],
        ")"}], "2"]}], "+", 
    RowBox[{
     FractionBox["1", "2"], 
     SuperscriptBox[
      RowBox[{"(", 
       FractionBox[
        RowBox[{"Xil", "-", 
         SubscriptBox["X", 
          RowBox[{"i", ",", 
           RowBox[{"l", "-", "1"}]}]]}], "\[CapitalDelta]\[Tau]"], ")"}], 
      "2"]}]}]}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{
    "change", " ", "in", " ", "phonon", " ", "energy", " ", "when", " ", 
     "updating", " ", 
     SubscriptBox["X", 
      RowBox[{"i", ",", "l"}]]}], " ", "\[Rule]", " ", 
    RowBox[{
     SubscriptBox["X", 
      RowBox[{"i", ",", "l"}]], "+", "\[CapitalDelta]X"}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"Eph", "[", 
     RowBox[{
      SubscriptBox["X", 
       RowBox[{"i", ",", "l"}]], "+", "\[CapitalDelta]X"}], "]"}], "-", 
    RowBox[{"Eph", "[", 
     SubscriptBox["X", 
      RowBox[{"i", ",", "l"}]], "]"}], "-", 
    RowBox[{"\[CapitalDelta]X", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        FractionBox["1", "2"], 
        SuperscriptBox["\[CapitalOmega]", "2"], 
        RowBox[{"(", 
         RowBox[{"\[CapitalDelta]X", "+", 
          RowBox[{"2", " ", 
           SubscriptBox["X", 
            RowBox[{"i", ",", "l"}]]}]}], ")"}]}], "+", 
       FractionBox[
        RowBox[{"\[CapitalDelta]X", "-", 
         RowBox[{"(", 
          RowBox[{
           SubscriptBox["X", 
            RowBox[{"i", ",", 
             RowBox[{"l", "+", "1"}]}]], "-", 
           RowBox[{"2", " ", 
            SubscriptBox["X", 
             RowBox[{"i", ",", "l"}]]}], "+", 
           SubscriptBox["X", 
            RowBox[{"i", ",", 
             RowBox[{"l", "-", "1"}]}]]}], ")"}]}], 
        SuperscriptBox["\[CapitalDelta]\[Tau]", "2"]]}], ")"}]}]}], 
   "]"}]}]], "Input"],

Cell[BoxData["0"], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"phonon", " ", "frequencies"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[CapitalOmega]", "val"], "=", 
    RowBox[{"{", 
     RowBox[{
      RowBox[{"13", "/", "10"}], ",", 
      RowBox[{"7", "/", "8"}]}], "}"}]}], ";"}]}]], "Input"],

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
  RowBox[{"(*", " ", 
   RowBox[{"electron", "-", 
    RowBox[{"phonon", " ", "interaction", " ", "strength"}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["g", "val"], "=", 
    RowBox[{"{", 
     RowBox[{
      RowBox[{"5", "/", "11"}], ",", 
      RowBox[{"7", "/", "10"}]}], "}"}]}], ";"}]}]], "Input"],

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
   "box", " ", "width", " ", "for", " ", "updates", " ", "of", " ", "the", 
    " ", "phonon", " ", "field", " ", "variables"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["box", "width"], "=", "12"}], ";"}]}]], "Input"]
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
  RowBox[{"HubbardDQMCPhononStep", "[", 
   RowBox[{
   "expK_", ",", "invexpK_", ",", "\[Tau]_", ",", "\[Lambda]_List", ",", 
    "\[CapitalOmega]_List", ",", "g_List", ",", "boxwidth_", ",", "s0_List", 
    ",", "X0_List", ",", "orbcellorderHS_List", ",", "orbcellorderPh_List", 
    ",", "rnd_List"}], "]"}], ":=", 
  RowBox[{"Module", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"s", "=", "s0"}], ",", 
      RowBox[{"X", "=", "X0"}], ",", "Gu", ",", "Gd", ",", "norbcell", ",", 
      "L", ",", "\[CapitalDelta]x", ",", "\[CapitalDelta]Eph", ",", "du", ",",
       "dd", ",", "cu", ",", "cd", ",", "l", ",", "i", ",", "j", ",", 
      RowBox[{"n", "=", "1"}]}], "}"}], ",", "\[IndentingNewLine]", 
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
         RowBox[{"DiagonalMatrix", "[", "\[Lambda]", "]"}], ".", "s"}], ",", 
        RowBox[{"\[Tau]", " ", 
         RowBox[{
          RowBox[{"DiagonalMatrix", "[", "g", "]"}], ".", "X0"}]}]}], "]"}]}],
      ";", "\[IndentingNewLine]", 
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
             RowBox[{
              RowBox[{"-", "\[Lambda]"}], " ", 
              RowBox[{"s", "\[LeftDoubleBracket]", 
               RowBox[{";;", ",", "l"}], "\[RightDoubleBracket]"}]}], "-", 
             RowBox[{"\[Tau]", " ", "g", " ", 
              RowBox[{"X", "\[LeftDoubleBracket]", 
               RowBox[{";;", ",", "l"}], "\[RightDoubleBracket]"}]}]}], "]"}],
            "]"}], ".", "expK", ".", "Gu", ".", "invexpK", ".", 
          RowBox[{"DiagonalMatrix", "[", 
           RowBox[{"Exp", "[", 
            RowBox[{
             RowBox[{
              RowBox[{"+", "\[Lambda]"}], " ", 
              RowBox[{"s", "\[LeftDoubleBracket]", 
               RowBox[{";;", ",", "l"}], "\[RightDoubleBracket]"}]}], "+", 
             RowBox[{"\[Tau]", " ", "g", " ", 
              RowBox[{"X", "\[LeftDoubleBracket]", 
               RowBox[{";;", ",", "l"}], "\[RightDoubleBracket]"}]}]}], "]"}],
            "]"}]}]}], ";", "\[IndentingNewLine]", 
        RowBox[{"Gd", "=", 
         RowBox[{
          RowBox[{"DiagonalMatrix", "[", 
           RowBox[{"Exp", "[", 
            RowBox[{
             RowBox[{
              RowBox[{"+", "\[Lambda]"}], " ", 
              RowBox[{"s", "\[LeftDoubleBracket]", 
               RowBox[{";;", ",", "l"}], "\[RightDoubleBracket]"}]}], "-", 
             RowBox[{"\[Tau]", " ", "g", " ", 
              RowBox[{"X", "\[LeftDoubleBracket]", 
               RowBox[{";;", ",", "l"}], "\[RightDoubleBracket]"}]}]}], "]"}],
            "]"}], ".", "expK", ".", "Gd", ".", "invexpK", ".", 
          RowBox[{"DiagonalMatrix", "[", 
           RowBox[{"Exp", "[", 
            RowBox[{
             RowBox[{
              RowBox[{"-", "\[Lambda]"}], " ", 
              RowBox[{"s", "\[LeftDoubleBracket]", 
               RowBox[{";;", ",", "l"}], "\[RightDoubleBracket]"}]}], "+", 
             RowBox[{"\[Tau]", " ", "g", " ", 
              RowBox[{"X", "\[LeftDoubleBracket]", 
               RowBox[{";;", ",", "l"}], "\[RightDoubleBracket]"}]}]}], "]"}],
            "]"}]}]}], ";", "\[IndentingNewLine]", 
        RowBox[{"(*", " ", 
         RowBox[{
          RowBox[{"iterate", " ", "over", " ", "lattice", " ", "sites"}], ",",
           " ", 
          RowBox[{
           RowBox[{"updating", " ", "the", " ", "Hubbard"}], "-", 
           RowBox[{"Stratonovich", " ", "field"}]}]}], " ", "*)"}], 
        "\[IndentingNewLine]", 
        RowBox[{"For", "[", 
         RowBox[{
          RowBox[{"j", "=", "1"}], ",", 
          RowBox[{"j", "\[LessEqual]", "norbcell"}], ",", 
          RowBox[{"j", "++"}], ",", "\[IndentingNewLine]", 
          RowBox[{
           RowBox[{"i", "=", 
            RowBox[{"orbcellorderHS", "\[LeftDoubleBracket]", 
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
               RowBox[{"n", "++"}], "\[RightDoubleBracket]"}], "<", 
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
            "]"}]}]}], "]"}], ";", "\[IndentingNewLine]", 
        RowBox[{"(*", " ", 
         RowBox[{
          RowBox[{"iterate", " ", "over", " ", "lattice", " ", "sites"}], ",",
           " ", 
          RowBox[{"updating", " ", "the", " ", "phonon", " ", "field"}]}], 
         " ", "*)"}], "\[IndentingNewLine]", 
        RowBox[{"For", "[", 
         RowBox[{
          RowBox[{"j", "=", "1"}], ",", 
          RowBox[{"j", "\[LessEqual]", "norbcell"}], ",", 
          RowBox[{"j", "++"}], ",", "\[IndentingNewLine]", 
          RowBox[{
           RowBox[{"i", "=", 
            RowBox[{"orbcellorderPh", "\[LeftDoubleBracket]", 
             RowBox[{"l", ",", "j"}], "\[RightDoubleBracket]"}]}], ";", 
           "\[IndentingNewLine]", 
           RowBox[{"\[CapitalDelta]x", "=", 
            RowBox[{
             RowBox[{"(", 
              RowBox[{
               RowBox[{"rnd", "\[LeftDoubleBracket]", 
                RowBox[{"n", "++"}], "\[RightDoubleBracket]"}], "-", 
               RowBox[{"1", "/", "2"}]}], ")"}], "boxwidth"}]}], ";", 
           "\[IndentingNewLine]", 
           RowBox[{"\[CapitalDelta]Eph", "=", 
            RowBox[{"\[CapitalDelta]x", 
             RowBox[{"(", 
              RowBox[{
               RowBox[{
                FractionBox["1", "2"], 
                SuperscriptBox[
                 RowBox[{
                 "\[CapitalOmega]", "\[LeftDoubleBracket]", "i", 
                  "\[RightDoubleBracket]"}], "2"], 
                RowBox[{"(", 
                 RowBox[{"\[CapitalDelta]x", "+", 
                  RowBox[{"2", " ", 
                   RowBox[{"X", "\[LeftDoubleBracket]", 
                    RowBox[{"i", ",", "l"}], "\[RightDoubleBracket]"}]}]}], 
                 ")"}]}], "+", 
               FractionBox[
                RowBox[{"\[CapitalDelta]x", "-", 
                 RowBox[{"(", 
                  RowBox[{
                   RowBox[{"X", "\[LeftDoubleBracket]", 
                    RowBox[{"i", ",", 
                    RowBox[{"Mod", "[", 
                    RowBox[{
                    RowBox[{"l", "+", "1"}], ",", "L", ",", "1"}], "]"}]}], 
                    "\[RightDoubleBracket]"}], "-", 
                   RowBox[{"2", 
                    RowBox[{"X", "\[LeftDoubleBracket]", 
                    RowBox[{"i", ",", "l"}], "\[RightDoubleBracket]"}]}], "+", 
                   RowBox[{"X", "\[LeftDoubleBracket]", 
                    RowBox[{"i", ",", 
                    RowBox[{"Mod", "[", 
                    RowBox[{
                    RowBox[{"l", "-", "1"}], ",", "L", ",", "1"}], "]"}]}], 
                    "\[RightDoubleBracket]"}]}], ")"}]}], 
                SuperscriptBox["\[Tau]", "2"]]}], ")"}]}]}], ";", 
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
                  RowBox[{"-", "\[Tau]"}], " ", 
                  RowBox[{
                  "g", "\[LeftDoubleBracket]", "i", "\[RightDoubleBracket]"}],
                   "\[CapitalDelta]x"}], "]"}], "-", "1"}], ")"}]}]}]}], ";", 
           "\[IndentingNewLine]", 
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
                  RowBox[{"-", "\[Tau]"}], " ", 
                  RowBox[{
                  "g", "\[LeftDoubleBracket]", "i", "\[RightDoubleBracket]"}],
                   "\[CapitalDelta]x"}], "]"}], "-", "1"}], ")"}]}]}]}], ";", 
           "\[IndentingNewLine]", 
           RowBox[{"If", "[", 
            RowBox[{"(*", 
             RowBox[{"RandomReal", "[", "]"}], "*)"}], 
            RowBox[{
             RowBox[{
              RowBox[{"rnd", "\[LeftDoubleBracket]", 
               RowBox[{"n", "++"}], "\[RightDoubleBracket]"}], "<", 
              RowBox[{"du", " ", "dd", " ", 
               RowBox[{"Exp", "[", 
                RowBox[{
                 RowBox[{"-", "\[Tau]"}], " ", "\[CapitalDelta]Eph"}], 
                "]"}]}]}], ",", "\[IndentingNewLine]", 
             RowBox[{
              RowBox[{"cu", "=", 
               RowBox[{
                RowBox[{"(", 
                 RowBox[{
                  RowBox[{"Exp", "[", 
                   RowBox[{
                    RowBox[{"-", "\[Tau]"}], " ", 
                    RowBox[{
                    "g", "\[LeftDoubleBracket]", "i", 
                    "\[RightDoubleBracket]"}], "\[CapitalDelta]x"}], "]"}], 
                  "-", "1"}], ")"}], 
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
                    RowBox[{"-", "\[Tau]"}], " ", 
                    RowBox[{
                    "g", "\[LeftDoubleBracket]", "i", 
                    "\[RightDoubleBracket]"}], "\[CapitalDelta]x"}], "]"}], 
                  "-", "1"}], ")"}], 
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
               RowBox[{
               "actually", " ", "update", " ", "phonon", " ", "field"}], " ", 
               "*)"}], "\[IndentingNewLine]", 
              RowBox[{
               RowBox[{"X", "\[LeftDoubleBracket]", 
                RowBox[{"i", ",", "l"}], "\[RightDoubleBracket]"}], "+=", 
               "\[CapitalDelta]x"}]}]}], "]"}]}]}], "]"}]}]}], "]"}], ";", 
     "\[IndentingNewLine]", 
     RowBox[{"(*", " ", 
      RowBox[{
       RowBox[{
        RowBox[{"return", " ", "updated", " ", "Hubbard"}], "-", 
        RowBox[{"Stratonovich", " ", "field"}]}], ",", " ", 
       RowBox[{
        RowBox[{"phonon", " ", "field", " ", "and", " ", "spin"}], "-", 
        RowBox[{"up", " ", "and", " ", "spin"}], "-", 
        RowBox[{"down", " ", 
         RowBox[{"Green", "'"}], "s", " ", "functions"}]}]}], " ", "*)"}], 
     "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"s", ",", "X", ",", "Gu", ",", "Gd"}], "}"}]}]}], 
   "]"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{
     RowBox[{"load", " ", "pre"}], "-", 
     RowBox[{
     "computed", " ", "random", " ", "numbers", " ", "from", " ", "disk"}]}], 
    ",", " ", 
    RowBox[{
    "to", " ", "enable", " ", "same", " ", "execution", " ", "path", " ", 
     "as", " ", "C", " ", "implementation"}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["rnd", "list"], "=", 
     RowBox[{"Import", "[", 
      RowBox[{
       RowBox[{
        SubscriptBox["fn", "base"], "<>", "\"\<_rnd.dat\>\""}], ",", 
       "\"\<Real64\>\""}], "]"}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"Length", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData["1920"], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   SubscriptBox["rnd", "list"], "\[LeftDoubleBracket]", 
   RowBox[{"{", 
    RowBox[{"1", ",", "2", ",", 
     RowBox[{"-", "1"}]}], "}"}], "\[RightDoubleBracket]"}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
  "0.95463862700127`", ",", "0.30392075394759765`", ",", 
   "0.7910658402987042`"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"pre", "-", 
    RowBox[{
    "determined", " ", "update", " ", "ordering", " ", "for", " ", "the", " ",
      "Hubbard"}], "-", 
    RowBox[{"Stratonovich", " ", "field"}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["order", "HS"], "=", 
     RowBox[{"Partition", "[", 
      RowBox[{
       RowBox[{
        RowBox[{"Import", "[", 
         RowBox[{
          RowBox[{
           SubscriptBox["fn", "base"], "<>", "\"\<_orbcellorderHS.dat\>\""}], 
          ",", "\"\<Integer32\>\""}], "]"}], "+", "1"}], ",", 
       RowBox[{
        SubscriptBox["n", "orb"], 
        SubscriptBox["n", "cells"]}]}], "]"}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"Dimensions", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"16", ",", "40"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"pre", "-", 
    RowBox[{
    "determined", " ", "update", " ", "ordering", " ", "for", " ", "the", " ",
      "phonon", " ", "field"}]}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["order", "Ph"], "=", 
     RowBox[{"Partition", "[", 
      RowBox[{
       RowBox[{
        RowBox[{"Import", "[", 
         RowBox[{
          RowBox[{
           SubscriptBox["fn", "base"], "<>", "\"\<_orbcellorderPh.dat\>\""}], 
          ",", "\"\<Integer32\>\""}], "]"}], "+", "1"}], ",", 
       RowBox[{
        SubscriptBox["n", "orb"], 
        SubscriptBox["n", "cells"]}]}], "]"}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"Dimensions", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"16", ",", "40"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "check", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"Total", "[", 
    RowBox[{
     RowBox[{
      RowBox[{"Norm", "[", 
       RowBox[{
        RowBox[{"Sort", "[", "#", "]"}], "-", 
        RowBox[{"Range", "[", 
         RowBox[{
          SubscriptBox["n", "orb"], 
          SubscriptBox["n", "cells"]}], "]"}]}], "]"}], "&"}], "/@", 
     SubscriptBox["order", "HS"]}], "]"}], "\[IndentingNewLine]", 
   RowBox[{"Total", "[", 
    RowBox[{
     RowBox[{
      RowBox[{"Norm", "[", 
       RowBox[{
        RowBox[{"Sort", "[", "#", "]"}], "-", 
        RowBox[{"Range", "[", 
         RowBox[{
          SubscriptBox["n", "orb"], 
          SubscriptBox["n", "cells"]}], "]"}]}], "]"}], "&"}], "/@", 
     SubscriptBox["order", "Ph"]}], "]"}]}]}]], "Input"],

Cell[BoxData["0"], "Output"],

Cell[BoxData["0"], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"perform", " ", "simulation", " ", "step"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      SubscriptBox["s", "1"], ",", 
      SubscriptBox["X", "1"], ",", 
      SubscriptBox["G", 
       RowBox[{"u", ",", "1"}]], ",", 
      SubscriptBox["G", 
       RowBox[{"d", ",", "1"}]]}], "}"}], "=", 
    RowBox[{"Block", "[", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{"$MinPrecision", "=", 
        RowBox[{"4", "MachinePrecision"}]}], "}"}], ",", 
      RowBox[{"HubbardDQMCPhononStep", "[", 
       RowBox[{
        SubscriptBox["exp", "\[Tau]k"], ",", 
        SubscriptBox["invexp", "\[Tau]k"], ",", 
        SubscriptBox["\[Tau]", "val"], ",", 
        SubscriptBox["\[Lambda]", "orb"], ",", 
        SubscriptBox["\[CapitalOmega]", "orb"], ",", 
        SubscriptBox["g", "orb"], ",", 
        SubscriptBox["box", "width"], ",", 
        SubscriptBox["s", "0"], ",", 
        SubscriptBox["X", "0"], ",", 
        SubscriptBox["order", "HS"], ",", 
        SubscriptBox["order", "Ph"], ",", 
        RowBox[{"SetPrecision", "[", 
         RowBox[{
          SubscriptBox["rnd", "list"], ",", 
          RowBox[{"4", "MachinePrecision"}]}], "]"}]}], "]"}]}], "]"}]}], 
   ";"}]}]], "Input"],

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
         SubscriptBox["\[Lambda]", "orb"], "]"}], ".", 
        SubscriptBox["s", "1"]}], ",", 
       RowBox[{
        SubscriptBox["\[Tau]", "val"], 
        RowBox[{
         RowBox[{"DiagonalMatrix", "[", 
          SubscriptBox["g", "orb"], "]"}], ".", 
         SubscriptBox["X", "1"]}]}]}], "]"}]}], "]"}], "]"}]}]], "Input"],

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
      RowBox[{"-", "1"}], "1"},
     {"1", "1", 
      RowBox[{"-", "1"}], "1"},
     {
      RowBox[{"-", "1"}], 
      RowBox[{"-", "1"}], "1", 
      RowBox[{"-", "1"}]},
     {
      RowBox[{"-", "1"}], "1", 
      RowBox[{"-", "1"}], "1"}
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
     {"0.250449230644419984002563167707258524800136453301150810284188381066722\
527592690909760560710756019`63.81835908076401", 
      RowBox[{
      "-", "0.0641579264163802492603079720713269814569444240891724482964722170\
93925073948831240901780332723971`63.81835908076401"}], 
      RowBox[{
      "-", "0.0363605177984747026319369475712399192927475379751446169218810853\
54277408282995627915608336904378`63.81835908076401"}]},
     {
      RowBox[{
      "-", "0.0404738651756093869531365629609581887043702310449340037355089458\
828151186191901378291325796556`63.81835908076402"}], 
      "0.543114690616723605467999347152907745031535156491706562068853200407623\
312579318908264221967858942`63.81835908076401", 
      "0.085188425819877205231956442756951655086996053045776486497063524588864\
110187020436657983620352804`63.81835908076402"},
     {
      RowBox[{
      "-", "0.0369301559898929439234043257046731569450019255617620877426734681\
90286322751824251607302301190424`63.81835908076402"}], 
      "0.019852686567528416009439307634695730916373399885568222502774447032302\
08022475309644324528441495`63.81835908076401", 
      "0.349853497289340056641843149939762453583523771025017715969958275530947\
337711208044663307953766996`63.81835908076401"}
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
     {"0.541925299485414038075643274819166023084354089628569347715394503926191\
10116109761525839974251791`63.81835908076401", 
      RowBox[{
      "-", "0.1549688576316925866378607977077381941354352283379355908721215142\
22445682257050734469746109787463`63.81835908076401"}], 
      "0.027014143717429278950527596838542774387541821437378905682717838979830\
030406203967317855205192376`63.81835908076401"},
     {
      RowBox[{
      "-", "0.1232008340464675267275259776614790049901514805197132166951987617\
94494261115966508607013332903932`63.81835908076401"}], 
      "0.293353727112416043815609697314956643432414910976071088664657601816190\
814360597906186410195848852`63.81835908076401", 
      RowBox[{
      "-", "0.0295190780760339393237566807994992203323734635115348321407821336\
52187574028791731983391636691597`63.81835908076401"}]},
     {
      RowBox[{
      "-", "0.0102203950171383735341747913431900152650393247443915078580651344\
99773593122434183995227200452906`63.81835908076401"}], 
      RowBox[{
      "-", "0.0407529682557452627427520444327313765585433184254693105233426622\
7885545179669894565890355775619`63.81835908076401"}], 
      "0.290245545873426656880288650650860019712391603881882893725612262942187\
546307732348401611098063979`63.81835908076402"}
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

Cell[BoxData["1.6844349640025118`"], "Output"],

Cell[BoxData["1.3903239150995421`"], "Output"]
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

Cell[BoxData["7.\
504484897851563772452349009081128945379310057172555634923250835821229143584769\
412663608`63.81835908076401*^-67"], "Output"],

Cell[BoxData["1.\
944675105248578249015308196313931342325270095862120041350346839221084338190797\
25855238043021`63.81835908076401*^-62"], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{
    "save", " ", "initial", " ", "and", " ", "updated", " ", "Hubbard"}], "-", 
    RowBox[{"Stratonovich", " ", "field", " ", "to", " ", "disk"}]}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"Table", "[", 
    RowBox[{
     RowBox[{"Export", "[", 
      RowBox[{
       RowBox[{
        SubscriptBox["fn", "base"], "<>", "\"\<_HS\>\"", "<>", 
        RowBox[{"ToString", "[", "i", "]"}], "<>", "\"\<.dat\>\""}], ",", 
       RowBox[{
        RowBox[{
         RowBox[{"Flatten", "[", 
          RowBox[{"Transpose", "[", 
           SubscriptBox["s", "i"], "]"}], "]"}], "/.", 
         RowBox[{"{", 
          RowBox[{"1", "\[Rule]", "0"}], "}"}]}], "/.", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{"-", "1"}], "\[Rule]", "1"}], "}"}]}], ",", 
       "\"\<Integer8\>\""}], "]"}], ",", 
     RowBox[{"{", 
      RowBox[{"i", ",", "0", ",", "1"}], "}"}]}], "]"}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
   "save", " ", "initial", " ", "and", " ", "updated", " ", "phonon", " ", 
    "fields", " ", "to", " ", "disk"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"Table", "[", 
    RowBox[{
     RowBox[{"Export", "[", 
      RowBox[{
       RowBox[{
        SubscriptBox["fn", "base"], "<>", "\"\<_X\>\"", "<>", 
        RowBox[{"ToString", "[", "i", "]"}], "<>", "\"\<.dat\>\""}], ",", 
       RowBox[{"Flatten", "[", 
        RowBox[{"Transpose", "[", 
         RowBox[{"N", "[", 
          SubscriptBox["X", "i"], "]"}], "]"}], "]"}], ",", 
       "\"\<Real64\>\""}], "]"}], ",", 
     RowBox[{"{", 
      RowBox[{"i", ",", "0", ",", "1"}], "}"}]}], "]"}], ";"}]}]], "Input"],

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
WindowSize->{1508, 978},
WindowMargins->{{Automatic, 211}, {48, Automatic}},
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
Cell[4917, 162, 2813, 82, 72, "Output"],
Cell[7733, 246, 74, 2, 31, "Output"]
}, Open  ]],
Cell[7822, 251, 2374, 86, 211, "Input"],
Cell[10199, 339, 236, 7, 67, "Input"],
Cell[10438, 348, 387, 12, 52, "Input"],
Cell[CellGroupData[{
Cell[10850, 364, 4470, 110, 188, "Input"],
Cell[15323, 476, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[15435, 483, 456, 13, 52, "Input"],
Cell[15894, 498, 1245, 43, 126, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[17176, 546, 328, 10, 72, "Input"],
Cell[17507, 558, 28, 0, 31, "Output"]
}, Open  ]],
Cell[17550, 561, 513, 16, 52, "Input"],
Cell[CellGroupData[{
Cell[18088, 581, 480, 14, 52, "Input"],
Cell[18571, 597, 788, 18, 71, "Output"]
}, Open  ]],
Cell[19374, 618, 470, 14, 52, "Input"],
Cell[CellGroupData[{
Cell[19869, 636, 377, 11, 52, "Input"],
Cell[20249, 649, 48, 0, 31, "Output"]
}, Open  ]],
Cell[20312, 652, 2386, 71, 152, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[22735, 728, 48, 0, 43, "Subsection"],
Cell[CellGroupData[{
Cell[22808, 732, 614, 20, 72, "Input"],
Cell[23425, 754, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[23537, 761, 454, 13, 52, "Input"],
Cell[23994, 776, 835, 27, 86, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[24866, 808, 412, 13, 72, "Input"],
Cell[25281, 823, 94, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[25412, 830, 581, 18, 72, "Input"],
Cell[25996, 850, 112, 3, 31, "Output"]
}, Open  ]],
Cell[26123, 856, 350, 11, 31, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[26510, 872, 34, 0, 43, "Subsection"],
Cell[CellGroupData[{
Cell[26569, 876, 642, 20, 72, "Input"],
Cell[27214, 898, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[27326, 905, 429, 13, 52, "Input"],
Cell[27758, 920, 498, 14, 31, "Output"]
}, Open  ]],
Cell[28271, 937, 1068, 34, 69, "Input"],
Cell[CellGroupData[{
Cell[29364, 975, 1530, 47, 68, "Input"],
Cell[30897, 1024, 28, 0, 31, "Output"]
}, Open  ]],
Cell[30940, 1027, 330, 10, 52, "Input"],
Cell[31273, 1039, 362, 11, 31, "Input"],
Cell[31638, 1052, 368, 12, 52, "Input"],
Cell[32009, 1066, 334, 11, 31, "Input"],
Cell[32346, 1079, 315, 9, 52, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[32698, 1093, 31, 0, 43, "Subsection"],
Cell[32732, 1095, 625, 17, 52, "Input"],
Cell[33360, 1114, 843, 22, 52, "Input"],
Cell[34206, 1138, 612, 17, 72, "Input"],
Cell[34821, 1157, 19037, 437, 888, "Input"],
Cell[CellGroupData[{
Cell[53883, 1598, 720, 22, 72, "Input"],
Cell[54606, 1622, 31, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[54674, 1627, 283, 7, 52, "Input"],
Cell[54960, 1636, 143, 4, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[55140, 1645, 809, 24, 72, "Input"],
Cell[55952, 1671, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[56064, 1678, 769, 22, 72, "Input"],
Cell[56836, 1702, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[56948, 1709, 825, 25, 72, "Input"],
Cell[57776, 1736, 28, 0, 31, "Output"],
Cell[57807, 1738, 28, 0, 31, "Output"]
}, Open  ]],
Cell[57850, 1741, 1309, 37, 72, "Input"],
Cell[CellGroupData[{
Cell[59184, 1782, 1019, 29, 52, "Input"],
Cell[60206, 1813, 28, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[60271, 1818, 454, 13, 52, "Input"],
Cell[60728, 1833, 835, 27, 86, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[61600, 1865, 877, 27, 72, "Input"],
Cell[62480, 1894, 1811, 39, 71, "Output"],
Cell[64294, 1935, 1834, 40, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[66165, 1980, 559, 17, 72, "Input"],
Cell[66727, 1999, 46, 0, 31, "Output"],
Cell[66776, 2001, 46, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[66859, 2006, 476, 15, 72, "Input"],
Cell[67338, 2023, 143, 2, 31, "Output"],
Cell[67484, 2027, 148, 2, 31, "Output"]
}, Open  ]],
Cell[67647, 2032, 990, 28, 52, "Input"],
Cell[68640, 2062, 749, 20, 52, "Input"],
Cell[69392, 2084, 1483, 43, 72, "Input"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature nw0ckPi2Fno2bBKBpBtRuL87 *)
