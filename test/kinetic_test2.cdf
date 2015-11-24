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
NotebookDataLength[     26877,        860]
NotebookOptionsPosition[     26310,        821]
NotebookOutlinePosition[     26654,        836]
CellTagsIndexPosition[     26611,        833]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{
Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"lattice", " ", "dimensions"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["n", "x"], "=", "5"}], ";"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{
     SubscriptBox["n", "y"], "=", "4"}], ";"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"number", " ", "of", " ", "orbitals"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["n", "orb"], "=", "3"}], ";"}]}]], "Input"],

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
    RowBox[{"0", ",", "4", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "1", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "2", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "3", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "4", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "1", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "2", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "3", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "4", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "1", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "2", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "3", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "4", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "1", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "2", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "3", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "4", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "0", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "1", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "2", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "3", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "4", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "0", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "1", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "2", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "3", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "4", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "0", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "1", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "2", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "3", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "4", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "1", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "2", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "3", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "4", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "0", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "1", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "2", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "3", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "4", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "0", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "1", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "2", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "3", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "4", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "0", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "1", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "2", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "3", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "4", ",", "3"}], "}"}]}], "}"}]], "Output"],

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
         FractionBox["1", "2"], ",", 
         RowBox[{"-", 
          FractionBox["2", "3"]}]}], "}"}], ",", 
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
         FractionBox["3", "19"], ",", 
         FractionBox["2", "9"], ",", 
         RowBox[{"-", 
          FractionBox["3", "10"]}]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         FractionBox["1", "7"], ",", 
         FractionBox["1", "6"], ",", 
         RowBox[{"-", 
          FractionBox["1", "4"]}]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"-", 
          FractionBox["2", "5"]}], ",", 
         FractionBox["3", "16"], ",", 
         RowBox[{"-", 
          FractionBox["5", "17"]}]}], "}"}]}], "}"}]}], ";"}], 
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
       RowBox[{"1", ",", "1"}], "}"}]], "=", 
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
       RowBox[{"1", ",", 
        RowBox[{"-", "1"}]}], "}"}]], "=", 
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
         FractionBox["1", "17"]}], "}"}]}], "}"}]}], ";"}]}]}]], "Input"],

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
       RowBox[{"-", "2"}], "/", "11"}], ",", 
      RowBox[{"1", "/", "3"}], ",", 
      RowBox[{"4", "/", "13"}]}], "}"}]}], ";"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"factor", " ", 
    RowBox[{"(", 
     RowBox[{"-", "1"}], ")"}], " ", "from", " ", "negative", " ", "sign", 
    " ", "in", " ", "Hamiltonian"}], " ", "*)"}], "\[IndentingNewLine]", 
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
       FractionBox["69", "77"]}], 
      RowBox[{"-", 
       FractionBox["3", "19"]}], "0", 
      FractionBox["1", "3"]},
     {
      RowBox[{"-", 
       FractionBox["3", "19"]}], 
      RowBox[{"-", 
       FractionBox["69", "77"]}], 
      RowBox[{"-", 
       FractionBox["3", "19"]}], "0"},
     {"0", 
      RowBox[{"-", 
       FractionBox["3", "19"]}], 
      RowBox[{"-", 
       FractionBox["69", "77"]}], "0"},
     {
      FractionBox["1", "3"], "0", "0", 
      RowBox[{"-", 
       FractionBox["37", "91"]}]}
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
   RowBox[{"time", " ", "step"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[CapitalDelta]t", "val"], "=", 
    FractionBox["1", "9"]}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   SuperscriptBox["\[ExponentialE]", 
    RowBox[{
     RowBox[{"-", "\[CapitalDelta]t"}], " ", "K"}]], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["exp", "K"], "=", 
    RowBox[{"MatrixExp", "[", 
     RowBox[{"N", "[", 
      RowBox[{
       RowBox[{"-", 
        SubscriptBox["\[CapitalDelta]t", "val"]}], 
       SubscriptBox["K", "val"]}], "]"}], "]"}]}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   SuperscriptBox["\[ExponentialE]", 
    RowBox[{"\[CapitalDelta]t", " ", "K"}]], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["invexp", "K"], "=", 
    RowBox[{"MatrixExp", "[", 
     RowBox[{"N", "[", 
      RowBox[{
       SubscriptBox["\[CapitalDelta]t", "val"], 
       SubscriptBox["K", "val"]}], "]"}], "]"}]}], ";"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["exp", "K"], "\[LeftDoubleBracket]", 
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
     {"1.1302165541842928`", "0.03287893893911974`", 
      RowBox[{"-", "0.0430041330966574`"}]},
     {"0.03287893893911974`", "1.1302165541842923`", 
      RowBox[{"-", "0.0007629924827896992`"}]},
     {
      RowBox[{"-", "0.043004133096657425`"}], 
      RowBox[{"-", "0.0007629924827896997`"}], "1.0585704138355398`"}
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
    SubscriptBox["invexp", "K"], "\[LeftDoubleBracket]", 
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
     {"0.9246484506748494`", 
      RowBox[{"-", "0.006267074898421067`"}], "0.033062397792559844`"},
     {
      RowBox[{"-", "0.006267074898421059`"}], "0.92464845067485`", 
      RowBox[{"-", "0.00002771031972013051`"}]},
     {"0.03306239779255987`", 
      RowBox[{"-", "0.000027710319720130336`"}], "0.9665913423233302`"}
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
   RowBox[{"check", ":", " ", 
    RowBox[{"inverse", " ", "matrix"}]}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{"Norm", "[", 
   RowBox[{
    RowBox[{
     SubscriptBox["exp", "K"], ".", 
     SubscriptBox["invexp", "K"]}], "-", 
    RowBox[{"IdentityMatrix", "[", 
     RowBox[{"Length", "[", 
      SubscriptBox["exp", "K"], "]"}], "]"}]}], "]"}]}]], "Input"],

Cell[BoxData["1.7849803102991325`*^-15"], "Output"]
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
     ";"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
   "save", " ", "matrix", " ", "exponentials", " ", "as", " ", "reference", 
    " ", "to", " ", "disk"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{"Export", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["fn", "base"], "<>", "\"\<_expK.dat\>\""}], ",", 
      RowBox[{"Flatten", "[", 
       RowBox[{"Transpose", "[", 
        SubscriptBox["exp", "K"], "]"}], "]"}], ",", "\"\<Real64\>\""}], 
     "]"}], ";"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"Export", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["fn", "base"], "<>", "\"\<_invexpK.dat\>\""}], ",", 
      RowBox[{"Flatten", "[", 
       RowBox[{"Transpose", "[", 
        SubscriptBox["invexp", "K"], "]"}], "]"}], ",", "\"\<Real64\>\""}], 
     "]"}], ";"}]}]}]], "Input"]
},
WindowSize->{1454, 838},
WindowMargins->{{Automatic, 323}, {110, Automatic}},
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
Cell[1464, 33, 320, 10, 72, "Input"],
Cell[1787, 45, 220, 7, 52, "Input"],
Cell[CellGroupData[{
Cell[2032, 56, 1259, 38, 72, "Input"],
Cell[3294, 96, 4193, 122, 92, "Output"],
Cell[7490, 220, 74, 2, 31, "Output"]
}, Open  ]],
Cell[7579, 225, 3716, 127, 211, "Input"],
Cell[11298, 354, 236, 7, 67, "Input"],
Cell[11537, 363, 425, 13, 52, "Input"],
Cell[CellGroupData[{
Cell[11987, 380, 4399, 108, 188, "Input"],
Cell[16389, 490, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[16501, 497, 456, 13, 52, "Input"],
Cell[16960, 512, 1105, 37, 126, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[18102, 554, 328, 10, 72, "Input"],
Cell[18433, 566, 28, 0, 31, "Output"]
}, Open  ]],
Cell[18476, 569, 238, 7, 67, "Input"],
Cell[18717, 578, 461, 15, 52, "Input"],
Cell[19181, 595, 419, 13, 52, "Input"],
Cell[CellGroupData[{
Cell[19625, 612, 436, 13, 52, "Input"],
Cell[20064, 627, 877, 22, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[20978, 654, 439, 13, 52, "Input"],
Cell[21420, 669, 881, 22, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[22338, 696, 415, 12, 52, "Input"],
Cell[22756, 710, 51, 0, 31, "Output"]
}, Open  ]],
Cell[22822, 713, 242, 7, 31, "Input"],
Cell[23067, 722, 2386, 71, 152, "Input"],
Cell[25456, 795, 850, 24, 72, "Input"]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature 0w0#Q3K@#uW56BKmuen1Kow9 *)
