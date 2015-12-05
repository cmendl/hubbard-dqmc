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
NotebookDataLength[     38000,       1198]
NotebookOptionsPosition[     35882,       1105]
NotebookOutlinePosition[     36225,       1120]
CellTagsIndexPosition[     36182,       1117]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["General parameters", "Subsection"],

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
  RowBox[{"(*", " ", 
   RowBox[{
   "number", " ", "of", " ", "orbitals", " ", "per", " ", "unit", " ", 
    "cell"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["n", "orb"], "=", "2"}], ";"}]}]], "Input"],

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
         RowBox[{"-", 
          FractionBox["2", "3"]}]}], "}"}], ",", 
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
         FractionBox["1", "5"], ",", 
         RowBox[{"-", 
          FractionBox["1", "9"]}]}], "}"}], ",", 
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
         RowBox[{"-", 
          FractionBox["1", "10"]}], ",", 
         FractionBox["2", "17"]}], "}"}]}], "}"}]}], ";"}], 
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
      RowBox[{"1", "/", "5"}], ",", 
      RowBox[{
       RowBox[{"-", "2"}], "/", "11"}]}], "}"}]}], ";"}]}]], "Input"],

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
       FractionBox["23", "35"]}], 
      RowBox[{"-", 
       FractionBox["4", "9"]}], "0", 
      FractionBox["1", "10"]},
     {
      RowBox[{"-", 
       FractionBox["4", "9"]}], 
      RowBox[{"-", 
       FractionBox["23", "35"]}], 
      RowBox[{"-", 
       FractionBox["4", "9"]}], "0"},
     {"0", 
      RowBox[{"-", 
       FractionBox["4", "9"]}], 
      RowBox[{"-", 
       FractionBox["23", "35"]}], 
      RowBox[{"-", 
       FractionBox["1", "6"]}]},
     {
      FractionBox["1", "10"], "0", 
      RowBox[{"-", 
       FractionBox["1", "6"]}], 
      RowBox[{"-", 
       FractionBox["80", "77"]}]}
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
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Hubbard-Stratonovich field", "Subsection"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"SeedRandom", "[", "42", "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
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
         RowBox[{
          SubscriptBox["n", "orb"], 
          SubscriptBox["n", "sites"]}], ",", 
         SubscriptBox["L", "val"]}], "}"}]}], "]"}]}], "-", "1"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"Dimensions", "[", 
  SubscriptBox["s", "val"], "]"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"48", ",", "16"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["s", "val"], "\[LeftDoubleBracket]", 
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
    RowBox[{"save", " ", "Hubbard"}], "-", 
    RowBox[{"Stratonovich", " ", "field", " ", "to", " ", "disk"}]}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"Export", "[", 
    RowBox[{
     RowBox[{
      SubscriptBox["fn", "base"], "<>", "\"\<_HS.dat\>\""}], ",", 
     RowBox[{
      RowBox[{
       RowBox[{"Flatten", "[", 
        RowBox[{"Transpose", "[", 
         SubscriptBox["s", "val"], "]"}], "]"}], "/.", 
       RowBox[{"{", 
        RowBox[{"1", "\[Rule]", "0"}], "}"}]}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"-", "1"}], "\[Rule]", "1"}], "}"}]}], ",", 
     "\"\<Integer8\>\""}], "]"}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{
     RowBox[{"\[Lambda]", " ", "parameter", " ", "for", " ", "Hubbard"}], "-", 
     RowBox[{"Stratonovich", " ", "field"}]}], ",", " ", 
    RowBox[{"per", " ", "orbital"}]}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[Lambda]", "val"], "=", 
    RowBox[{"{", 
     RowBox[{
      RowBox[{"3", "/", "4"}], ",", 
      RowBox[{"5", "/", "23"}]}], "}"}]}], ";"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[{
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
    SubscriptBox["orb", "ind"]}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{"Dimensions", "[", "%", "]"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", "48", "}"}]], "Output"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Phonon field", "Subsection"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"SeedRandom", "[", "42", "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   SubscriptBox["X", "val"], "=", 
   RowBox[{"RandomReal", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
       RowBox[{"-", "4"}], ",", "4"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{
        SubscriptBox["n", "orb"], 
        SubscriptBox["n", "sites"]}], ",", 
       SubscriptBox["L", "val"]}], "}"}], ",", 
     RowBox[{"WorkingPrecision", "\[Rule]", 
      RowBox[{"4", "MachinePrecision"}]}]}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"Dimensions", "[", "%", "]"}]}], "Input"],

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
    SubscriptBox["X", "val"], "\[LeftDoubleBracket]", 
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
    RowBox[{"2.9975697663870666`", ",", 
     RowBox[{"-", "2.960671810704432`"}], ",", 
     RowBox[{"-", "3.757827926459985`"}]}], "}"}]}], "}"}]], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"save", " ", "phonon", " ", "field", " ", "to", " ", "disk"}], " ",
    "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"Export", "[", 
    RowBox[{
     RowBox[{
      RowBox[{"NotebookDirectory", "[", "]"}], "<>", 
      RowBox[{"FileBaseName", "[", 
       RowBox[{"NotebookFileName", "[", "]"}], "]"}], "<>", 
      "\"\<_X.dat\>\""}], ",", 
     RowBox[{"Flatten", "[", 
      RowBox[{"Transpose", "[", 
       RowBox[{"N", "[", 
        SubscriptBox["X", "val"], "]"}], "]"}], "]"}], ",", 
     "\"\<Real64\>\""}], "]"}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"electron", "-", 
     RowBox[{"phonon", " ", "interaction", " ", "strength"}]}], ",", " ", 
    RowBox[{"per", " ", "orbital"}]}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["g", "val"], "=", 
    RowBox[{"{", 
     RowBox[{
      RowBox[{"1", "/", "7"}], ",", 
      RowBox[{"3", "/", "13"}]}], "}"}]}], ";"}]}]], "Input"],

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
    SubscriptBox["orb", "ind"]}]}], ";"}]], "Input"]
}, Open  ]],

Cell[CellGroupData[{

Cell["Construct time flow map", "Subsection"],

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

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{
   SubscriptBox["A", "val"], "=", 
   RowBox[{"Block", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"$MinPrecision", "=", 
       RowBox[{"4", "MachinePrecision"}]}], "}"}], ",", 
     RowBox[{"HubbardTimeFlowMap", "[", 
      RowBox[{
       RowBox[{"MatrixExp", "[", 
        RowBox[{"N", "[", 
         RowBox[{
          RowBox[{
           RowBox[{"-", 
            SubscriptBox["\[Tau]", "val"]}], 
           SubscriptBox["K", "val"]}], ",", 
          RowBox[{"4", "MachinePrecision"}]}], "]"}], "]"}], ",", 
       RowBox[{"Exp", "[", 
        RowBox[{
         RowBox[{"-", 
          RowBox[{
           RowBox[{"DiagonalMatrix", "[", 
            SubscriptBox["\[Lambda]", "orb"], "]"}], ".", 
           SubscriptBox["s", "val"]}]}], "-", 
         RowBox[{
          SubscriptBox["\[Tau]", "val"], " ", 
          RowBox[{
           RowBox[{"DiagonalMatrix", "[", 
            SubscriptBox["g", "orb"], "]"}], ".", 
           SubscriptBox["X", "val"]}]}]}], "]"}]}], "]"}]}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"Dimensions", "[", "%", "]"}]}], "Input"],

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
     {"3628.\
074325885854828899007417768079083441400934494961673761681202441869224981589449\
749126589362783783`63.81835908076401", 
      "13162.63419616803781583896875512894091067843901751522682664603903280466\
1631410372685150850374696416692`63.81835908076401", 
      "580.4799537736649635322789953728926872967697555612931790979747667276242\
09733605888092510448895076657`63.81835908076401"},
     {"2079.\
228043950868561044581773475250174953258361264566230166246690032320569762112932\
906344353554258367`63.81835908076401", 
      "11170.83210275412884301819236482555394459395146687240541038206937563064\
1210733816733639580517217528603`63.81835908076402", 
      "398.1453953613406037280421854058173091964329659420822714314001276257149\
42779976788761810756644603969`63.81835908076401"},
     {
      RowBox[{
      "-", "28.881652411342642487299404951084572978638976947323622444259660093\
259367985477482421003367238918561`63.81835908076401"}], 
      RowBox[{
      "-", "7.2348192897514409291638605516961000871988919095757535088760762373\
44360017726363730230891480624052`63.81835908076401"}], 
      "200.9939070120678848511093270720001767500140485184303824241982602008794\
70701657899469562205519663975`63.81835908076401"}
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
  RowBox[{"Eigenvalues", "[", 
   SubscriptBox["A", "val"], "]"}], "]"}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
  "28086.395021820965`", ",", "12826.135036739914`", ",", 
   "5817.694489608892`", ",", "4003.4489156684417`", ",", 
   "3523.064809665367`", ",", "2725.540116449287`", ",", 
   "2003.1072875019781`", ",", "1212.425913363389`", ",", "802.242293624396`",
    ",", "732.7233203325856`", ",", "637.9773279211716`", ",", 
   "374.67027051142173`", ",", "310.2692912412569`", ",", 
   "231.54284896525428`", ",", "196.78010553838195`", ",", 
   "156.53600373706848`", ",", "107.26566000801793`", ",", 
   RowBox[{"87.63001568728562`", "\[VeryThinSpace]", "+", 
    RowBox[{"7.308275261804366`", " ", "\[ImaginaryI]"}]}], ",", 
   RowBox[{"87.63001568728562`", "\[VeryThinSpace]", "-", 
    RowBox[{"7.308275261804366`", " ", "\[ImaginaryI]"}]}], ",", 
   "69.84536245447849`", ",", "50.172951033560544`", ",", "39.6799008835994`",
    ",", "28.677739727840788`", ",", "14.327007216605779`", ",", 
   "9.519279259556468`", ",", "5.594078255971631`", ",", 
   RowBox[{"4.068063412458823`", "\[VeryThinSpace]", "+", 
    RowBox[{"0.020910762919016588`", " ", "\[ImaginaryI]"}]}], ",", 
   RowBox[{"4.068063412458823`", "\[VeryThinSpace]", "-", 
    RowBox[{"0.020910762919016588`", " ", "\[ImaginaryI]"}]}], ",", 
   "2.6587340021110104`", ",", "2.3273819236951367`", ",", 
   "1.6173917033187666`", ",", "1.2777228629997148`", ",", 
   "0.9157578886557188`", ",", "0.8788919669770905`", ",", 
   "0.5596613156115796`", ",", "0.3880314368071604`", ",", 
   "0.34003160367165686`", ",", "0.24995890234415213`", ",", 
   "0.2300030277164786`", ",", "0.17969312143891633`", ",", 
   "0.12617900916697128`", ",", "0.08200689941245107`", ",", 
   "0.0427712523442815`", ",", "0.026558153490960076`", ",", 
   "0.013833729869148924`", ",", "0.01020217093157794`", ",", 
   "0.003878952285642207`", ",", "0.0014790860445305085`"}], "}"}]], "Output"]
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
  "62274.0045926152`", ",", "26988.730144604608`", ",", "7644.580410292427`", 
   ",", "6806.368953830373`", ",", "5139.393608885351`", ",", 
   "3515.531500349713`", ",", "1588.733813272009`", ",", 
   "1439.6855814572218`", ",", "1030.8820610204114`", ",", 
   "766.2380009400157`", ",", "631.6593866519243`", ",", "419.5746221926564`",
    ",", "330.85436294608417`", ",", "237.2244518737809`", ",", 
   "194.99256031069356`", ",", "176.33788076439308`", ",", 
   "118.64077343179618`", ",", "100.52156015274798`", ",", 
   "75.91235073193563`", ",", "72.48776068719313`", ",", 
   "46.347125380903854`", ",", "40.74997486623602`", ",", 
   "24.724802779293068`", ",", "14.62358715917631`", ",", 
   "8.884941964894383`", ",", "5.2551985876946645`", ",", 
   "4.030001417174745`", ",", "3.6325225644469654`", ",", 
   "2.6031020153078055`", ",", "2.1953118456401586`", ",", 
   "1.537437223327407`", ",", "1.287500913825893`", ",", "0.854962190508695`",
    ",", "0.7956086091770265`", ",", "0.520458458605726`", ",", 
   "0.3870966515075309`", ",", "0.3124371662975475`", ",", 
   "0.25544313654493217`", ",", "0.2230761912738197`", ",", 
   "0.1620601119374107`", ",", "0.12568450866184297`", ",", 
   "0.06203059058465274`", ",", "0.036631409727902936`", ",", 
   "0.020550117871558984`", ",", "0.010895016733680204`", ",", 
   "0.00753343936198815`", ",", "0.0021428061602469166`", ",", 
   "0.0006575563827147119`"}], "}"}]], "Output"],

Cell[BoxData["9.470519369839874`*^7"], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
   "smallest", " ", "singular", " ", "value", " ", "now", " ", "of", " ", 
    "order", " ", "1"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"N", "[", 
    RowBox[{"SingularValueList", "[", 
     RowBox[{
      RowBox[{"IdentityMatrix", "[", 
       RowBox[{
        SubscriptBox["n", "orb"], 
        SubscriptBox["n", "sites"]}], "]"}], "+", 
      SubscriptBox["A", "val"]}], "]"}], "]"}], "\[IndentingNewLine]", 
   RowBox[{"(*", " ", 
    RowBox[{"condition", " ", "number"}], " ", "*)"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"Max", "[", "%", "]"}], "/", 
    RowBox[{"Min", "[", "%", "]"}]}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
  "62274.518750861076`", ",", "26989.1406517933`", ",", "7645.096375349083`", 
   ",", "6807.155265131755`", ",", "5140.193313097351`", ",", 
   "3516.2669961450424`", ",", "1589.1410914136459`", ",", 
   "1440.4139171817367`", ",", "1031.6285684324264`", ",", 
   "767.1425870363299`", ",", "632.4367993277245`", ",", "420.2080938491361`",
    ",", "331.6576891560326`", ",", "238.05751607791757`", ",", 
   "195.732820513518`", ",", "177.11450841192865`", ",", 
   "119.45155888535635`", ",", "101.35864967464224`", ",", 
   "76.7563876955038`", ",", "73.25235194717499`", ",", "47.0354460741086`", 
   ",", "41.50747689738799`", ",", "25.44659215119376`", ",", 
   "15.402663713745996`", ",", "9.741087110848738`", ",", 
   "6.108986934287546`", ",", "4.906674993193011`", ",", "4.515016404828249`",
    ",", "3.4976219530907486`", ",", "3.046382890425379`", ",", 
   "2.4112225618883074`", ",", "2.089327106853559`", ",", "1.67618559008973`",
    ",", "1.6228370469812705`", ",", "1.4147699593611416`", ",", 
   "1.2533372991367968`", ",", "1.1531048132428832`", ",", 
   "1.1244294904184615`", ",", "1.1171877539263035`", ",", 
   "1.0439969924283485`", ",", "0.9806680416629223`", ",", 
   "0.9464728815576519`", ",", "0.9363172089085301`", ",", 
   "0.8254018679184089`", ",", "0.7911164658063126`", ",", 
   "0.7492492365077089`", ",", "0.7202727411639754`", ",", 
   "0.6893011109257003`"}], "}"}]], "Output"],

Cell[BoxData["90344.4340416472`"], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"save", " ", "as", " ", "reference", " ", "to", " ", "disk"}], " ",
    "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"Export", "[", 
    RowBox[{
     RowBox[{
      SubscriptBox["fn", "base"], "<>", "\"\<_A.dat\>\""}], ",", 
     RowBox[{"Flatten", "[", 
      RowBox[{"Transpose", "[", 
       SubscriptBox["A", "val"], "]"}], "]"}], ",", "\"\<Real64\>\""}], "]"}],
    ";"}]}]], "Input"]
}, Open  ]]
},
WindowSize->{1509, 925},
WindowMargins->{{Automatic, 147}, {43, Automatic}},
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
Cell[1529, 37, 219, 7, 52, "Input"],
Cell[1751, 46, 261, 8, 67, "Input"],
Cell[CellGroupData[{
Cell[2037, 58, 295, 9, 69, "Input"],
Cell[2335, 69, 29, 0, 31, "Output"]
}, Open  ]],
Cell[2379, 72, 264, 8, 52, "Input"],
Cell[2646, 82, 242, 7, 31, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[2925, 94, 41, 0, 43, "Subsection"],
Cell[2969, 96, 320, 10, 72, "Input"],
Cell[3292, 108, 161, 6, 31, "Input"],
Cell[CellGroupData[{
Cell[3478, 118, 1259, 38, 72, "Input"],
Cell[4740, 158, 3365, 98, 72, "Output"],
Cell[8108, 258, 74, 2, 31, "Output"]
}, Open  ]],
Cell[8197, 263, 2400, 87, 211, "Input"],
Cell[10600, 352, 236, 7, 67, "Input"],
Cell[10839, 361, 387, 12, 52, "Input"],
Cell[CellGroupData[{
Cell[11251, 377, 4470, 110, 188, "Input"],
Cell[15724, 489, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[15836, 496, 456, 13, 52, "Input"],
Cell[16295, 511, 1199, 41, 126, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[17531, 557, 328, 10, 72, "Input"],
Cell[17862, 569, 28, 0, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[17939, 575, 48, 0, 43, "Subsection"],
Cell[CellGroupData[{
Cell[18012, 579, 618, 20, 72, "Input"],
Cell[18633, 601, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[18745, 608, 456, 13, 52, "Input"],
Cell[19204, 623, 835, 27, 86, "Output"]
}, Open  ]],
Cell[20054, 653, 722, 22, 52, "Input"],
Cell[20779, 677, 482, 14, 52, "Input"],
Cell[CellGroupData[{
Cell[21286, 695, 416, 12, 52, "Input"],
Cell[21705, 709, 51, 1, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[21805, 716, 34, 0, 43, "Subsection"],
Cell[CellGroupData[{
Cell[21864, 720, 644, 20, 72, "Input"],
Cell[22511, 742, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[22623, 749, 431, 13, 52, "Input"],
Cell[23057, 764, 498, 14, 31, "Output"]
}, Open  ]],
Cell[23570, 781, 609, 17, 52, "Input"],
Cell[24182, 800, 428, 13, 52, "Input"],
Cell[24613, 815, 334, 11, 31, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[24984, 831, 45, 0, 43, "Subsection"],
Cell[25032, 833, 625, 17, 52, "Input"],
Cell[CellGroupData[{
Cell[25682, 854, 1128, 33, 52, "Input"],
Cell[26813, 889, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[26925, 896, 436, 13, 52, "Input"],
Cell[27364, 911, 1788, 38, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[29189, 954, 117, 3, 31, "Input"],
Cell[29309, 959, 1875, 31, 72, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[31221, 995, 392, 10, 72, "Input"],
Cell[31616, 1007, 1483, 25, 72, "Output"],
Cell[33102, 1034, 48, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[33187, 1039, 700, 19, 92, "Input"],
Cell[33890, 1060, 1458, 24, 72, "Output"],
Cell[35351, 1086, 44, 0, 31, "Output"]
}, Open  ]],
Cell[35410, 1089, 456, 13, 52, "Input"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature axp@grIicBsHiBKgLWL4LeNe *)
