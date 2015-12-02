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
NotebookDataLength[     40393,       1233]
NotebookOptionsPosition[     38650,       1153]
NotebookOutlinePosition[     38994,       1168]
CellTagsIndexPosition[     38951,       1165]
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
    SubscriptBox["n", "orb"], "=", "3"}], ";"}]}]], "Input"],

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
    RowBox[{"1", ",", "3", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "1", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "2", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "3", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "0", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "1", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "2", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "3", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "0", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "1", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "2", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "3", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "0", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "1", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "2", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "3", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "0", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "1", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "2", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "3", ",", "4"}], "}"}]}], "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"60", ",", "3"}], "}"}]], "Output"]
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
          FractionBox["2", "3"]}], ",", 
         FractionBox["1", "2"]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{"0", ",", "0", ",", 
         FractionBox["2", "11"]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{"0", ",", "0", ",", "0"}], "}"}]}], "}"}]}], ";"}], 
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
         FractionBox["3", "16"], ",", 
         RowBox[{"-", 
          FractionBox["1", "9"]}]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         FractionBox["3", "7"], ",", 
         RowBox[{"-", 
          FractionBox["2", "19"]}], ",", 
         RowBox[{"-", 
          FractionBox["1", "15"]}]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"-", 
          FractionBox["3", "8"]}], ",", "0", ",", 
         FractionBox["1", "17"]}], "}"}]}], "}"}]}], ";"}], 
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
          FractionBox["2", "7"]}], ",", 
         FractionBox["1", "13"]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         FractionBox["2", "9"], ",", 
         FractionBox["1", "13"], ",", 
         FractionBox["3", "7"]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"-", 
          FractionBox["1", "3"]}], ",", 
         FractionBox["1", "5"], ",", 
         FractionBox["2", "9"]}], "}"}]}], "}"}]}], ";"}], 
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
         FractionBox["4", "5"], ",", 
         FractionBox["1", "10"], ",", 
         RowBox[{"-", 
          FractionBox["3", "11"]}]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"-", 
          FractionBox["1", "10"]}], ",", 
         FractionBox["2", "17"], ",", 
         FractionBox["3", "20"]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"-", 
          FractionBox["1", "5"]}], ",", 
         FractionBox["1", "4"], ",", 
         RowBox[{"-", 
          FractionBox["1", "6"]}]}], "}"}]}], "}"}]}], ";"}], 
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
         FractionBox["1", "6"], ",", 
         RowBox[{"-", 
          FractionBox["3", "10"]}]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         FractionBox["1", "7"], ",", 
         RowBox[{"-", 
          FractionBox["1", "4"]}], ",", 
         FractionBox["2", "9"]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"-", 
          FractionBox["2", "5"]}], ",", 
         RowBox[{"-", 
          FractionBox["5", "17"]}], ",", 
         FractionBox["3", "16"]}], "}"}]}], "}"}]}], ";"}]}]}]], "Input"],

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
       RowBox[{"-", "1"}], "/", "3"}], ",", 
      RowBox[{"2", "/", "11"}], ",", 
      RowBox[{"4", "/", "17"}]}], "}"}]}], ";"}]}]], "Input"],

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
  RowBox[{"60", ",", "60"}], "}"}]], "Output"]
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
       FractionBox["22", "21"]}], 
      RowBox[{"-", 
       FractionBox["4", "9"]}], "0", 
      FractionBox["1", "5"]},
     {
      RowBox[{"-", 
       FractionBox["4", "9"]}], 
      RowBox[{"-", 
       FractionBox["22", "21"]}], 
      RowBox[{"-", 
       FractionBox["4", "9"]}], "0"},
     {"0", 
      RowBox[{"-", 
       FractionBox["4", "9"]}], 
      RowBox[{"-", 
       FractionBox["22", "21"]}], 
      FractionBox["3", "10"]},
     {
      FractionBox["1", "5"], "0", 
      FractionBox["3", "10"], 
      RowBox[{"-", 
       FractionBox["57", "119"]}]}
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
          SubscriptBox["n", "x"], 
          SubscriptBox["n", "y"]}], ",", 
         SubscriptBox["L", "val"]}], "}"}]}], "]"}]}], "-", "1"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"Dimensions", "[", 
  SubscriptBox["s", "val"], "]"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"60", ",", "16"}], "}"}]], "Output"]
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
     {
      RowBox[{"-", "1"}], "1", 
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
    RowBox[{"\[Lambda]", " ", "parameters", " ", "for", " ", "Hubbard"}], "-", 
    RowBox[{"Stratonovich", " ", "field"}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[Lambda]", "val"], "=", 
    RowBox[{"{", 
     RowBox[{
      RowBox[{"3", "/", "4"}], ",", 
      RowBox[{"2", "/", "5"}], ",", 
      RowBox[{"1", "/", "9"}]}], "}"}]}], ";"}]}]], "Input"],

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
    SubscriptBox["orb", "ind"]}]}], ";"}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{
   SubscriptBox["exp\[Lambda]s", "val"], "=", 
   RowBox[{"Exp", "[", 
    RowBox[{"-", 
     RowBox[{
      RowBox[{"DiagonalMatrix", "[", 
       SubscriptBox["\[Lambda]", "orb"], "]"}], ".", 
      SubscriptBox["s", "val"]}]}], "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{"Dimensions", "[", "%", "]"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"60", ",", "16"}], "}"}]], "Output"]
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
     "\"\<Integer8\>\""}], "]"}], ";"}]}]], "Input"]
}, Open  ]],

Cell[CellGroupData[{

Cell["Construct time flow map", "Subsection"],

Cell[BoxData[
 RowBox[{
  RowBox[{"HubbardTimeFlowMap", "[", 
   RowBox[{"exp\[Lambda]s_", ",", "exp\[Tau]k_"}], "]"}], ":=", 
  RowBox[{"Fold", "[", 
   RowBox[{
    RowBox[{
     RowBox[{
      RowBox[{"DiagonalMatrix", "[", "#2", "]"}], ".", "exp\[Tau]k", ".", 
      "#1"}], "&"}], ",", 
    RowBox[{"IdentityMatrix", "[", 
     RowBox[{"Length", "[", "exp\[Tau]k", "]"}], "]"}], ",", 
    RowBox[{"Transpose", "[", "exp\[Lambda]s", "]"}]}], "]"}]}]], "Input"],

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
       SubscriptBox["exp\[Lambda]s", "val"], ",", 
       RowBox[{"MatrixExp", "[", 
        RowBox[{"N", "[", 
         RowBox[{
          RowBox[{
           RowBox[{"-", 
            SubscriptBox["\[Tau]", "val"]}], 
           SubscriptBox["K", "val"]}], ",", 
          RowBox[{"4", "MachinePrecision"}]}], "]"}], "]"}]}], "]"}]}], 
    "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{"Dimensions", "[", "%", "]"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"60", ",", "60"}], "}"}]], "Output"]
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
     {"116530.\
670933538274088168967014538282376614176304352662454038277162781384638092813463\
63286748726543488`63.81835908076407", 
      "270382.4761346598780163427179071977816063616177286623116073486432904512\
08823031287436448204186425189186`63.81835908076401", 
      RowBox[{
      "-", "42342.\
258244295829355952563352669806295163752173408729735961926602394414391657427008\
352628407385518201`63.81835908076402"}]},
     {"103624.\
964439016825858172002474069236638973541112255947810836205133964644813543208315\
387423583534029876`63.81835908076402", 
      "246515.5997826520128534544489114293542817778845424996261968161392202143\
98668224637481567535174741454106`63.81835908076403", 
      RowBox[{
      "-", "38432.\
144420059128685041454495292140576964248255909611348653017490552114447775082916\
247034491014698617`63.81835908076401"}]},
     {
      RowBox[{
      "-", "13692.\
704593694147879084922739872922610426939519664463997840163161696815707627191430\
164432918568713298`63.81835908076402"}], 
      RowBox[{
      "-", "31521.\
204946665377976576934464315244871625775862000551356588804310492742113492341722\
00635746349689277`63.81835908076402"}], 
      "5186.991495035959290738989417917003755420540596345442966476605685424403\
889049771644219739485305333795`63.81835908076402"}
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
  "1.1943105805115812`*^6", ",", "31411.377515753768`", ",", 
   "26100.203206045888`", ",", "9367.454639122223`", ",", 
   "3940.110810220406`", ",", "2134.094464722419`", ",", 
   "1777.2633857443395`", ",", "1023.9477183044422`", ",", 
   RowBox[{"744.8472784504262`", "\[VeryThinSpace]", "+", 
    RowBox[{"96.51818503399792`", " ", "\[ImaginaryI]"}]}], ",", 
   RowBox[{"744.8472784504262`", "\[VeryThinSpace]", "-", 
    RowBox[{"96.51818503399792`", " ", "\[ImaginaryI]"}]}], ",", 
   "388.3032863981481`", ",", "347.4316624655997`", ",", 
   "283.27742269427114`", ",", "125.81628985381734`", ",", 
   "102.38894647808476`", ",", 
   RowBox[{"72.01973478554825`", "\[VeryThinSpace]", "+", 
    RowBox[{"0.7760727226661182`", " ", "\[ImaginaryI]"}]}], ",", 
   RowBox[{"72.01973478554825`", "\[VeryThinSpace]", "-", 
    RowBox[{"0.7760727226661182`", " ", "\[ImaginaryI]"}]}], ",", 
   "64.74169638784707`", ",", "54.512300055573704`", ",", 
   RowBox[{"42.90506656004126`", "\[VeryThinSpace]", "+", 
    RowBox[{"1.8799973518578286`", " ", "\[ImaginaryI]"}]}], ",", 
   RowBox[{"42.90506656004126`", "\[VeryThinSpace]", "-", 
    RowBox[{"1.8799973518578286`", " ", "\[ImaginaryI]"}]}], ",", 
   "32.961903795110544`", ",", "29.994718320676206`", ",", 
   "25.271833969035473`", ",", "18.72591062749691`", ",", 
   "17.895907294574044`", ",", "14.798031149327233`", ",", 
   "8.562498062417284`", ",", "7.030838706456286`", ",", "6.23856886884477`", 
   ",", "5.107412910534743`", ",", "4.4529201798857345`", ",", 
   "3.6447498916065353`", ",", "2.590635453675862`", ",", 
   "2.4244583081046516`", ",", "1.9632256223056448`", ",", 
   "1.6748281404080496`", ",", "1.3622310940231457`", ",", 
   "0.8425716908477633`", ",", "0.7645457195071264`", ",", 
   "0.5615166346831996`", ",", "0.42258655075662904`", ",", 
   "0.34602592727376263`", ",", "0.29373712806131114`", ",", 
   "0.2259670458469978`", ",", "0.19309609419982002`", ",", 
   "0.17421850863076685`", ",", "0.14542715798510275`", ",", 
   "0.11288206369163428`", ",", "0.1030142076000632`", ",", 
   "0.07735280025658517`", ",", 
   RowBox[{"0.06496086232919163`", "\[VeryThinSpace]", "+", 
    RowBox[{"0.0024917694371372734`", " ", "\[ImaginaryI]"}]}], ",", 
   RowBox[{"0.06496086232919163`", "\[VeryThinSpace]", "-", 
    RowBox[{"0.0024917694371372734`", " ", "\[ImaginaryI]"}]}], ",", 
   "0.04679608123810363`", ",", "0.03647274793979899`", ",", 
   "0.02358759260015041`", ",", "0.00825117522248877`", ",", 
   "0.0029698528060933586`", ",", "0.0022705543571599025`", ",", 
   "0.0001224325411897875`"}], "}"}]], "Output"]
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
  "1.7112770995274305`*^6", ",", "48168.9859579927`", ",", 
   "34216.83934399199`", ",", "12948.976284175547`", ",", "5104.43662982209`",
    ",", "2722.172567806904`", ",", "2233.897269917039`", ",", 
   "1468.914803929534`", ",", "930.7505291153878`", ",", "771.7257290020098`",
    ",", "473.2706627812935`", ",", "432.5439833404914`", ",", 
   "311.5908472227572`", ",", "145.37087821247775`", ",", 
   "104.59460396080915`", ",", "80.78306218383769`", ",", 
   "71.89646083603733`", ",", "67.24301476085978`", ",", 
   "48.974877059183996`", ",", "47.59944470455417`", ",", 
   "40.461552391639216`", ",", "36.29992505305329`", ",", 
   "30.154106615264762`", ",", "24.322404698351836`", ",", 
   "20.271665928258344`", ",", "18.704526129114853`", ",", 
   "14.864598136915859`", ",", "8.298560476748019`", ",", 
   "7.115219851444566`", ",", "6.139347453928914`", ",", "5.011183159807001`",
    ",", "4.349369000251751`", ",", "3.500861023747595`", ",", 
   "2.559153890333854`", ",", "2.323725591194609`", ",", "1.7839236464582`", 
   ",", "1.5346290689829964`", ",", "1.3178448435953833`", ",", 
   "0.7496563290366226`", ",", "0.7008921265192422`", ",", 
   "0.5149182044690876`", ",", "0.3715376982450651`", ",", 
   "0.3072889118701183`", ",", "0.2640020380249577`", ",", 
   "0.21099647527886162`", ",", "0.17817097103877413`", ",", 
   "0.1676137595511057`", ",", "0.12993558481591808`", ",", 
   "0.106664336014076`", ",", "0.08830550629391005`", ",", 
   "0.06278846890999484`", ",", "0.059823439588242454`", ",", 
   "0.05115137194753147`", ",", "0.03914272510597865`", ",", 
   "0.03186939504965228`", ",", "0.02018161037509554`", ",", 
   "0.007041733113372758`", ",", "0.00246578449699011`", ",", 
   "0.0015062127653343034`", ",", "0.00008447543780462951`"}], 
  "}"}]], "Output"],

Cell[BoxData["2.0257688435840782`*^10"], "Output"]
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
        SubscriptBox["n", "x"], 
        SubscriptBox["n", "y"]}], "]"}], "+", 
      SubscriptBox["A", "val"]}], "]"}], "]"}], "\[IndentingNewLine]", 
   RowBox[{"(*", " ", 
    RowBox[{"condition", " ", "number"}], " ", "*)"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"Max", "[", "%", "]"}], "/", 
    RowBox[{"Min", "[", "%", "]"}]}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
  "1.711277800181291`*^6", ",", "48169.64319210188`", ",", 
   "34217.540073231976`", ",", "12949.575699434208`", ",", 
   "5105.185193408555`", ",", "2722.890066982741`", ",", "2234.50181507089`", 
   ",", "1469.5736629799308`", ",", "931.612290769839`", ",", 
   "772.4414183438049`", ",", "473.9703254996999`", ",", "433.2190424855862`",
    ",", "312.29222866105135`", ",", "146.11127373561732`", ",", 
   "105.37191031413106`", ",", "81.59577863030917`", ",", 
   "72.68947261243353`", ",", "68.08344091597714`", ",", "49.79415477909124`",
    ",", "48.42123239182893`", ",", "41.37221082672775`", ",", 
   "37.1432353122868`", ",", "31.024003890709736`", ",", 
   "25.190898524415964`", ",", "21.158337261077484`", ",", 
   "19.566434528089946`", ",", "15.73210626109161`", ",", 
   "9.188054797027776`", ",", "7.989124909494025`", ",", "7.008259871230386`",
    ",", "5.912698998955182`", ",", "5.2185435036031445`", ",", 
   "4.298070594226643`", ",", "3.4469804898079914`", ",", "3.19099637744037`",
    ",", "2.7087276500517596`", ",", "2.383391541826035`", ",", 
   "2.197993844384455`", ",", "1.651119189867769`", ",", 
   "1.5987251646665275`", ",", "1.4186883637199772`", ",", 
   "1.2537407253770583`", ",", "1.2113068759749137`", ",", 
   "1.1599808888202436`", ",", "1.1381116093303838`", ",", 
   "1.1016095347022712`", ",", "1.069236022339998`", ",", 
   "1.0408630057130537`", ",", "1.0289259201399503`", ",", 
   "1.0023833266035567`", ",", "0.948659354658241`", ",", 
   "0.9360924116473238`", ",", "0.9266061479607107`", ",", 
   "0.9044587523326654`", ",", "0.8945337334209434`", ",", 
   "0.8642115979006003`", ",", "0.8614905621985885`", ",", 
   "0.8281288764416195`", ",", "0.7762002121299282`", ",", 
   "0.7008757692970699`"}], "}"}]], "Output"],

Cell[BoxData["2.4416278535318533`*^6"], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"save", " ", "time", " ", 
    RowBox[{"flow", " ", "'"}], 
    RowBox[{"A", "'"}], " ", "matrix", " ", "as", " ", "reference", " ", "to",
     " ", "disk"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"Export", "[", 
    RowBox[{
     RowBox[{
      RowBox[{"NotebookDirectory", "[", "]"}], "<>", 
      RowBox[{"FileBaseName", "[", 
       RowBox[{"NotebookFileName", "[", "]"}], "]"}], "<>", 
      "\"\<_A.dat\>\""}], ",", 
     RowBox[{"Flatten", "[", 
      RowBox[{"Transpose", "[", 
       SubscriptBox["A", "val"], "]"}], "]"}], ",", "\"\<Real64\>\""}], "]"}],
    ";"}]}]], "Input"]
}, Open  ]]
},
WindowSize->{1509, 867},
WindowMargins->{{Automatic, 217}, {Automatic, 117}},
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
Cell[2925, 94, 62, 0, 43, "Subsection"],
Cell[2990, 96, 320, 10, 72, "Input"],
Cell[CellGroupData[{
Cell[3335, 110, 1259, 38, 72, "Input"],
Cell[4597, 150, 4193, 122, 92, "Output"],
Cell[8793, 274, 74, 2, 31, "Output"]
}, Open  ]],
Cell[8882, 279, 3716, 127, 211, "Input"],
Cell[12601, 408, 236, 7, 67, "Input"],
Cell[12840, 417, 425, 13, 52, "Input"],
Cell[CellGroupData[{
Cell[13290, 434, 4470, 110, 188, "Input"],
Cell[17763, 546, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[17875, 553, 456, 13, 52, "Input"],
Cell[18334, 568, 1154, 39, 126, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[19525, 612, 328, 10, 72, "Input"],
Cell[19856, 624, 28, 0, 31, "Output"]
}, Open  ]],
Cell[19899, 627, 2386, 71, 152, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[22322, 703, 48, 0, 43, "Subsection"],
Cell[CellGroupData[{
Cell[22395, 707, 649, 21, 72, "Input"],
Cell[23047, 730, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[23159, 737, 456, 13, 52, "Input"],
Cell[23618, 752, 857, 28, 86, "Output"]
}, Open  ]],
Cell[24490, 783, 457, 14, 52, "Input"],
Cell[24950, 799, 350, 11, 31, "Input"],
Cell[CellGroupData[{
Cell[25325, 814, 357, 10, 52, "Input"],
Cell[25685, 826, 75, 2, 31, "Output"]
}, Open  ]],
Cell[25775, 831, 722, 22, 52, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[26534, 858, 45, 0, 43, "Subsection"],
Cell[26582, 860, 464, 12, 31, "Input"],
Cell[CellGroupData[{
Cell[27071, 876, 711, 21, 52, "Input"],
Cell[27785, 899, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[27897, 906, 436, 13, 52, "Input"],
Cell[28336, 921, 1856, 44, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[30229, 970, 117, 3, 31, "Input"],
Cell[30349, 975, 2643, 45, 96, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[33029, 1025, 392, 10, 72, "Input"],
Cell[33424, 1037, 1841, 31, 96, "Output"],
Cell[35268, 1070, 50, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[35355, 1075, 729, 20, 92, "Input"],
Cell[36087, 1097, 1815, 30, 96, "Output"],
Cell[37905, 1129, 49, 0, 31, "Output"]
}, Open  ]],
Cell[37969, 1132, 665, 18, 52, "Input"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature 3wT#QezIOAyenC1yLvec5aMt *)
