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
NotebookDataLength[     46689,       1413]
NotebookOptionsPosition[     44469,       1318]
NotebookOutlinePosition[     44812,       1333]
CellTagsIndexPosition[     44769,       1330]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["General parameters", "Subsection"],

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
   RowBox[{"inverse", " ", "temperature"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[Beta]", "val"], "=", "2"}], ";"}]}]], "Input"],

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
    RowBox[{"1", ",", "3", ",", "5"}], "}"}], ",", 
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
    RowBox[{"2", ",", "3", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "0", ",", "5"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "1", ",", "5"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "2", ",", "5"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "3", ",", "5"}], "}"}]}], "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"72", ",", "3"}], "}"}]], "Output"]
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
          FractionBox["1", "2"]}], ",", 
         FractionBox["2", "3"]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{"0", ",", "0", ",", 
         FractionBox["3", "11"]}], "}"}], ",", 
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
         FractionBox["3", "19"], ",", 
         FractionBox["1", "6"], ",", 
         RowBox[{"-", 
          FractionBox["3", "10"]}]}], "}"}], ",", 
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
          FractionBox["3", "8"]}], ",", 
         FractionBox["1", "17"], ",", 
         FractionBox["4", "9"]}], "}"}]}], "}"}]}], ";"}], 
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
         FractionBox["4", "9"], ",", 
         FractionBox["3", "16"], ",", 
         RowBox[{"-", 
          FractionBox["1", "9"]}]}], "}"}]}], "}"}]}], ";"}], 
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
          FractionBox["2", "5"]}], ",", 
         RowBox[{"-", 
          FractionBox["5", "17"]}], ",", 
         FractionBox["3", "16"]}], "}"}], ",", 
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
         RowBox[{"-", 
          FractionBox["1", "3"]}], ",", 
         FractionBox["1", "5"], ",", 
         FractionBox["2", "9"]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         FractionBox["1", "7"], ",", 
         RowBox[{"-", 
          FractionBox["1", "4"]}], ",", 
         FractionBox["2", "9"]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"-", 
          FractionBox["1", "10"]}], ",", 
         FractionBox["2", "17"], ",", 
         FractionBox["3", "20"]}], "}"}]}], "}"}]}], ";"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"chemical", " ", "potential"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[Mu]", "val"], "=", 
    FractionBox["5", "6"]}], ";"}]}]], "Input"],

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
       RowBox[{"-", "1"}], "/", "13"}], ",", 
      RowBox[{"5", "/", "17"}], ",", 
      RowBox[{"4", "/", "11"}]}], "}"}]}], ";"}]}]], "Input"],

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
  RowBox[{"72", ",", "72"}], "}"}]], "Output"]
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
       FractionBox["71", "78"]}], 
      RowBox[{"-", 
       FractionBox["3", "19"]}], "0", 
      FractionBox["1", "5"]},
     {
      RowBox[{"-", 
       FractionBox["3", "19"]}], 
      RowBox[{"-", 
       FractionBox["71", "78"]}], 
      RowBox[{"-", 
       FractionBox["3", "19"]}], "0"},
     {"0", 
      RowBox[{"-", 
       FractionBox["3", "19"]}], 
      RowBox[{"-", 
       FractionBox["71", "78"]}], 
      RowBox[{"-", 
       FractionBox["2", "9"]}]},
     {
      FractionBox["1", "5"], "0", 
      RowBox[{"-", 
       FractionBox["2", "9"]}], 
      RowBox[{"-", 
       FractionBox["31", "66"]}]}
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
  RowBox[{"Block", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"\[Beta]", "=", 
       SubscriptBox["\[Beta]", "val"]}], ",", "L"}], "}"}], ",", 
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
            RowBox[{
             SubscriptBox["n", "orb"], 
             SubscriptBox["n", "sites"]}], ",", "L"}], "}"}]}], "]"}]}], "-", 
       "1"}]}]}]}], "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{"Dimensions", "[", 
  SubscriptBox["s", "val"], "]"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"72", ",", "16"}], "}"}]], "Output"]
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
     {"1", "1", 
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
      RowBox[{"1", "/", "2"}], ",", 
      RowBox[{"2", "/", "9"}], ",", 
      RowBox[{"1", "/", "5"}]}], "}"}]}], ";"}]}]], "Input"],

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
  RowBox[{"72", ",", "16"}], "}"}]], "Output"]
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

Cell["Calculate time flow map", "Subsection"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
   "time", " ", "flow", " ", "map", " ", "generated", " ", "by", " ", "the", 
    " ", "Hubbard", " ", "Hamiltonian"}], " ", "*)"}], "\[IndentingNewLine]", 
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
     RowBox[{"Transpose", "[", "exp\[Lambda]s", "]"}]}], "]"}]}]}]], "Input"],

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
  RowBox[{"72", ",", "72"}], "}"}]], "Output"]
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
     {"2249.\
373553558472690342848597565385830551129512626044173689752622282752541562281203\
002528007518032017`63.81835908076401", 
      RowBox[{
      "-", "848.\
546361405019417064747334994733248748920377840096257420749151565203938268046446\
033792447221046818`63.81835908076401"}], 
      RowBox[{
      "-", "13.919291849747459384229755773320393656659202648417982767459035886\
502718904144137781077493804643878`63.81835908076401"}]},
     {
      RowBox[{
      "-", "234.\
154676061233502159176107055707747086497675443917022429391911270529190046350839\
389559153103221441`63.81835908076401"}], 
      "4944.084206718673278082658226180485328702619951571350803630937879319425\
706617989493359632822725195461`63.81835908076402", 
      "131.2557369646846333668781311245156427853470678285856835623911693858129\
50148602728922314108291015996`63.81835908076401"},
     {
      RowBox[{
      "-", "2.3992145663402266556731683881987477966755386490034185775378535196\
64209344141888929967147076334915`63.81835908076401"}], 
      RowBox[{
      "-", "345.\
589486747499845393902340875365130770337779092524856739937525031634301716459016\
092066095204070444`63.81835908076402"}], 
      "20.05423137036258677721216530807804702502146513003048961428199424325371\
8618098267924889505606600082`63.81835908076402"}
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
  "29644.196351225113`", ",", "15377.674429782011`", ",", 
   "8979.098364227886`", ",", "5206.974280161469`", ",", 
   "2286.7573500471804`", ",", "1288.0687973561455`", ",", "914.3710889696`", 
   ",", "752.0490955952713`", ",", "628.7330293479923`", ",", 
   "496.3253202802233`", ",", "355.5343937985266`", ",", 
   "282.49743218252206`", ",", "224.5161340896359`", ",", 
   "186.78665852132528`", ",", "168.1185750403378`", ",", 
   "128.40954567095855`", ",", "118.70261052182846`", ",", 
   "84.62568406604254`", ",", "77.2653368555953`", ",", "60.975750020960184`",
    ",", "55.399422067919836`", ",", "43.88469755326174`", ",", 
   "39.32231062574385`", ",", "32.34473915357883`", ",", 
   "27.762291703789245`", ",", "21.759342347891206`", ",", 
   "18.25420915048884`", ",", "16.148406029411326`", ",", 
   "14.262291595456222`", ",", "12.587075904234707`", ",", 
   "8.90132439999925`", ",", "8.329041454519224`", ",", "7.375639764767439`", 
   ",", "7.099531808714575`", ",", "6.275137175700959`", ",", 
   "5.397461687157978`", ",", "4.1710162192410865`", ",", 
   "3.7257333947317313`", ",", "3.245274175553451`", ",", 
   "2.974104208490202`", ",", "2.5561804003546102`", ",", 
   "2.3344070251446376`", ",", "1.9709495342418388`", ",", 
   "1.701660730946939`", ",", "1.5714908651011457`", ",", 
   "1.4406014261978428`", ",", "1.1959482228032157`", ",", 
   "1.0888688893670804`", ",", "0.8404700785704189`", ",", 
   "0.6359956807267838`", ",", "0.5803193358805867`", ",", 
   "0.538867257480588`", ",", "0.4909878711151644`", ",", 
   "0.32949690982568675`", ",", "0.24597366511267088`", ",", 
   "0.21073888492776233`", ",", "0.1705014217377063`", ",", 
   "0.1569031773341827`", ",", "0.12974079566565808`", ",", 
   "0.07291905598071416`", ",", "0.07185786053440571`", ",", 
   "0.06158605385913982`", ",", "0.04604593768230283`", ",", 
   "0.03612153750862489`", ",", "0.027699827881565497`", ",", 
   "0.022970739332393582`", ",", "0.015984220306399694`", ",", 
   "0.010294885791289812`", ",", "0.0067383287976807565`", ",", 
   "0.005233893824390784`", ",", "0.004109758934259306`", ",", 
   "0.0014544949755053024`"}], "}"}]], "Output"],

Cell[BoxData["2.038109230382628`*^7"], "Output"]
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
  "29644.90346816242`", ",", "15378.401101356401`", ",", "8979.880843009973`",
    ",", "5207.685171549286`", ",", "2287.640504615514`", ",", 
   "1288.9449836787103`", ",", "915.2058792414092`", ",", 
   "752.8747567287788`", ",", "629.5555364081731`", ",", "497.180279604545`", 
   ",", "356.3585215009262`", ",", "283.4348316561722`", ",", 
   "225.36564503928173`", ",", "187.6113180633596`", ",", 
   "168.97058210808373`", ",", "129.32128977441843`", ",", 
   "119.59892488284402`", ",", "85.5410580530765`", ",", "78.1941505749232`", 
   ",", "61.84944298586887`", ",", "56.32562806632867`", ",", 
   "44.75543124121582`", ",", "40.23416493339816`", ",", 
   "33.245843886047695`", ",", "28.649106121334746`", ",", 
   "22.68443809399725`", ",", "19.15689910359412`", ",", 
   "17.049629490597983`", ",", "15.147367401967596`", ",", 
   "13.527930061283824`", ",", "9.84895895311257`", ",", "9.240315692930876`",
    ",", "8.317131298454004`", ",", "8.038075885438335`", ",", 
   "7.210948507752773`", ",", "6.336658980173299`", ",", "5.108448941645463`",
    ",", "4.664742299552786`", ",", "4.183576234627444`", ",", 
   "3.919767049178549`", ",", "3.455684096677288`", ",", 
   "3.2588731226857655`", ",", "2.931169706003079`", ",", 
   "2.6519854176708293`", ",", "2.483197478687333`", ",", 
   "2.3578586628950933`", ",", "2.127641659294479`", ",", 
   "2.003577076533303`", ",", "1.7752387472352806`", ",", 
   "1.5744510934311386`", ",", "1.5008850454727698`", ",", 
   "1.443737427808705`", ",", "1.425619594632849`", ",", "1.286770986951555`",
    ",", "1.1995568266775223`", ",", "1.135758810571556`", ",", 
   "1.12970769635475`", ",", "1.0997444229683218`", ",", "1.07705339051581`", 
   ",", "1.0329264281888355`", ",", "1.0222847009865292`", ",", 
   "1.0186412760731396`", ",", "1.0012426707742352`", ",", 
   "0.9857024242277319`", ",", "0.9744621024733517`", ",", 
   "0.9506053457719483`", ",", "0.9377542831136102`", ",", 
   "0.917861162371252`", ",", "0.8956343288447202`", ",", 
   "0.8897713471037944`", ",", "0.865192623700201`", ",", 
   "0.8303933876945967`"}], "}"}]], "Output"],

Cell[BoxData["35699.83083615939`"], "Output"]
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
      RowBox[{
       SubscriptBox["n", "orb"], 
       SubscriptBox["n", "sites"]}], "]"}], "+", 
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
     {"0.143411766138543981963879486416646360691582985784013800032910193399886\
287865452987187606677999312`63.81835908076401", 
      RowBox[{
      "-", "0.0340423165292098604759124138620908142442103812505696454162404167\
92848306656813603627561307762428`63.81835908076401"}], 
      "0.081852214999803877593607343812988385463965266816463679966247957895591\
031949993068838812131954244`63.81835908076401"},
     {
      RowBox[{
      "-", "0.0197412889401725193311779115374679900062176995237165023480196961\
22613631485291045716084568224988`63.81835908076401"}], 
      "0.107436840534006204409950900117912985854017550851632262757261496680885\
048298900162433383644793086`63.81835908076401", 
      RowBox[{
      "-", "0.0020742298835716584658138817686556788935136860723115454670614318\
42338475008892167175854023142992`63.81835908076401"}]},
     {"0.030683923539780046563866042507725126472882974897946264562226839170631\
527609856505220976833305254`63.81835908076401", 
      "0.009920619842925400259256349072660452876240033493174025671882938626310\
113805568947541608560973176`63.81835908076401", 
      "0.510824327265275801269033053220028562324408867067671627244357068427746\
414081591348418597191764974`63.81835908076401"}
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
  "1.2042485101865736`", ",", "1.1558119806006475`", ",", 
   "1.1238842465033292`", ",", "1.1165271001725685`", ",", 
   "1.0894893922916904`", ",", "1.0663774274426305`", ",", 
   "1.0519612628392492`", ",", "1.0262071736415697`", ",", 
   "1.014504961559235`", ",", "0.9987588715398292`", ",", 
   "0.9816998618542127`", ",", "0.9782010813963822`", ",", 
   "0.9681231622211761`", ",", "0.9284590799357604`", ",", 
   "0.9093021788652481`", ",", "0.8851847280732172`", ",", 
   "0.8804686265183035`", ",", "0.8336412062859534`", ",", 
   "0.7771390637032203`", ",", "0.701449393488126`", ",", 
   "0.6926467242161848`", ",", "0.666273545076869`", ",", 
   "0.6351419895938081`", ",", "0.5633045141434518`", ",", 
   "0.49910732744569714`", ",", "0.47000395749517226`", ",", 
   "0.42411363146430137`", ",", "0.4027065944544288`", ",", 
   "0.3770759798816218`", ",", "0.34116073114155937`", ",", 
   "0.3068545360169962`", ",", "0.2893783031155888`", ",", 
   "0.2551172014698083`", ",", "0.23902994565343497`", ",", 
   "0.21437411453487393`", ",", "0.19575413426328458`", ",", 
   "0.15781186949288084`", ",", "0.13867801148834452`", ",", 
   "0.12440788246495482`", ",", "0.12023376379616381`", ",", 
   "0.10822141074303668`", ",", "0.10153357372699473`", ",", 
   "0.0739211391151366`", ",", "0.06601807254441475`", ",", 
   "0.05865230095184473`", ",", "0.0522005150516445`", ",", 
   "0.04408308444125048`", ",", "0.034905102999193`", ",", 
   "0.0300789477153164`", ",", "0.024854498699186506`", ",", 
   "0.02234365689854169`", ",", "0.01775390766743705`", ",", 
   "0.01616829435680571`", ",", "0.012788680389101883`", ",", 
   "0.011690292623917747`", ",", "0.008361279175206415`", ",", 
   "0.007732678832266131`", ",", "0.005918189944805539`", ",", 
   "0.005330168831617518`", ",", "0.004437233544738808`", ",", 
   "0.0035281478784974295`", ",", "0.0028061627256397766`", ",", 
   "0.002011342848906629`", ",", "0.0015884222156242762`", ",", 
   "0.0013282421691822607`", ",", "0.001092650323475713`", ",", 
   "0.0007758283035059824`", ",", "0.0004371316201048254`", ",", 
   "0.00019202389680989494`", ",", "0.00011136005226376805`", ",", 
   "0.00006502626595633523`", ",", "0.00003373261110713228`"}], 
  "}"}]], "Output"]
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

Cell[BoxData["0.6832326232780012`"], "Output"],

Cell[BoxData[
 RowBox[{"-", "0.3814985531670959`"}]], "Output"],

Cell[BoxData["5.543778581274894`*^-6"], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "determinant", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   SubscriptBox["detG", "val"], "=", 
   RowBox[{"Det", "[", 
    SubscriptBox["G", "val"], "]"}]}]}]], "Input"],

Cell[BoxData["4.\
674757420982416594614260119798431067677823300093947428264846623717470027220210\
727487661606`63.81835908076401*^-83"], "Output"]
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
WindowSize->{1458, 940},
WindowMargins->{{Automatic, 316}, {72, Automatic}},
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
Cell[1529, 37, 228, 7, 67, "Input"],
Cell[1760, 46, 219, 7, 52, "Input"],
Cell[1982, 55, 264, 8, 52, "Input"],
Cell[2249, 65, 242, 7, 31, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[2528, 77, 62, 0, 43, "Subsection"],
Cell[2593, 79, 320, 10, 72, "Input"],
Cell[2916, 91, 161, 6, 31, "Input"],
Cell[CellGroupData[{
Cell[3102, 101, 1259, 38, 72, "Input"],
Cell[4364, 141, 5021, 146, 112, "Output"],
Cell[9388, 289, 74, 2, 31, "Output"]
}, Open  ]],
Cell[9477, 294, 3744, 128, 211, "Input"],
Cell[13224, 424, 236, 7, 67, "Input"],
Cell[13463, 433, 426, 13, 52, "Input"],
Cell[CellGroupData[{
Cell[13914, 450, 4470, 110, 188, "Input"],
Cell[18387, 562, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[18499, 569, 456, 13, 52, "Input"],
Cell[18958, 584, 1201, 41, 126, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[20196, 630, 328, 10, 72, "Input"],
Cell[20527, 642, 28, 0, 31, "Output"]
}, Open  ]],
Cell[20570, 645, 2386, 71, 152, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[22993, 721, 48, 0, 43, "Subsection"],
Cell[CellGroupData[{
Cell[23066, 725, 922, 29, 90, "Input"],
Cell[23991, 756, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[24103, 763, 456, 13, 52, "Input"],
Cell[24562, 778, 813, 26, 86, "Output"]
}, Open  ]],
Cell[25390, 807, 457, 14, 52, "Input"],
Cell[25850, 823, 350, 11, 31, "Input"],
Cell[CellGroupData[{
Cell[26225, 838, 357, 10, 52, "Input"],
Cell[26585, 850, 75, 2, 31, "Output"]
}, Open  ]],
Cell[26675, 855, 722, 22, 52, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[27434, 882, 45, 0, 43, "Subsection"],
Cell[27482, 884, 679, 17, 52, "Input"],
Cell[CellGroupData[{
Cell[28186, 905, 711, 21, 52, "Input"],
Cell[28900, 928, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[29012, 935, 436, 13, 52, "Input"],
Cell[29451, 950, 1860, 44, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[31348, 999, 392, 10, 72, "Input"],
Cell[31743, 1011, 2205, 37, 112, "Output"],
Cell[33951, 1050, 48, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[34036, 1055, 783, 21, 92, "Input"],
Cell[34822, 1078, 2152, 35, 92, "Output"],
Cell[36977, 1115, 45, 0, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[37071, 1121, 72, 0, 43, "Subsection"],
Cell[37146, 1123, 305, 10, 31, "Input"],
Cell[CellGroupData[{
Cell[37476, 1137, 436, 13, 52, "Input"],
Cell[37915, 1152, 1785, 37, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[39737, 1194, 123, 3, 31, "Input"],
Cell[39863, 1199, 2284, 39, 112, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[42184, 1243, 651, 18, 92, "Input"],
Cell[42838, 1263, 46, 0, 31, "Output"],
Cell[42887, 1265, 63, 1, 31, "Output"],
Cell[42953, 1268, 49, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[43039, 1273, 220, 6, 52, "Input"],
Cell[43262, 1281, 146, 2, 31, "Output"]
}, Open  ]],
Cell[43423, 1286, 1030, 29, 72, "Input"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature VwTbdnRY7cOw7DgiDluNOsmK *)
