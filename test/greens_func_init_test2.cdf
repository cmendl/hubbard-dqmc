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
NotebookDataLength[     43972,       1234]
NotebookOptionsPosition[     41591,       1131]
NotebookOutlinePosition[     41935,       1146]
CellTagsIndexPosition[     41892,       1143]
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
   RowBox[{"time", " ", "step"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[Tau]", "val"], "=", 
    FractionBox["1", "8"]}], ";"}]}]], "Input"],

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
   RowBox[{
    RowBox[{"-", 
     SubscriptBox["latt", "neigh"]}], "-", 
    RowBox[{
     SubscriptBox["\[Mu]", "val"], 
     RowBox[{"IdentityMatrix", "[", 
      SubscriptBox["n", "sites"], "]"}]}]}]}], ";"}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"symbolic", " ", 
    SuperscriptBox["\[ExponentialE]", 
     RowBox[{
      RowBox[{"-", "\[Tau]"}], " ", "k"}]]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["exp", "\[Tau]k"], "=", 
    RowBox[{"FullSimplify", "[", 
     RowBox[{"MatrixExp", "[", 
      RowBox[{
       RowBox[{"-", 
        SubscriptBox["\[Tau]", "val"]}], 
       SubscriptBox["k", "val"]}], "]"}], "]"}]}], ";"}]}]], "Input"],

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
      FractionBox[
       RowBox[{
        SuperscriptBox[
         RowBox[{"(", 
          RowBox[{"1", "+", 
           SuperscriptBox["\[ExponentialE]", 
            RowBox[{"1", "/", "4"}]]}], ")"}], "2"], " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", 
          RowBox[{"2", " ", 
           SuperscriptBox["\[ExponentialE]", 
            RowBox[{"1", "/", "8"}]]}], "+", 
          RowBox[{"2", " ", 
           SuperscriptBox["\[ExponentialE]", 
            RowBox[{"3", "/", "8"}]]}], "+", 
          SqrtBox["\[ExponentialE]"]}], ")"}]}], 
       RowBox[{"24", " ", 
        SuperscriptBox["\[ExponentialE]", 
         RowBox[{"71", "/", "136"}]]}]], 
      FractionBox[
       RowBox[{
        RowBox[{"-", "1"}], "-", 
        RowBox[{"2", " ", 
         SuperscriptBox["\[ExponentialE]", 
          RowBox[{"1", "/", "8"}]]}], "-", 
        RowBox[{"2", " ", 
         SuperscriptBox["\[ExponentialE]", 
          RowBox[{"3", "/", "8"}]]}], "+", 
        RowBox[{"2", " ", 
         SuperscriptBox["\[ExponentialE]", 
          RowBox[{"5", "/", "8"}]]}], "+", 
        RowBox[{"2", " ", 
         SuperscriptBox["\[ExponentialE]", 
          RowBox[{"7", "/", "8"}]]}], "+", "\[ExponentialE]"}], 
       RowBox[{"24", " ", 
        SuperscriptBox["\[ExponentialE]", 
         RowBox[{"71", "/", "136"}]]}]], 
      FractionBox[
       RowBox[{
        SuperscriptBox[
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "+", 
           SuperscriptBox["\[ExponentialE]", 
            RowBox[{"1", "/", "4"}]]}], ")"}], "2"], " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", 
          SuperscriptBox["\[ExponentialE]", 
           RowBox[{"1", "/", "4"}]]}], ")"}], " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", 
          SuperscriptBox["\[ExponentialE]", 
           RowBox[{"1", "/", "8"}]], "+", 
          SuperscriptBox["\[ExponentialE]", 
           RowBox[{"1", "/", "4"}]]}], ")"}]}], 
       RowBox[{"24", " ", 
        SuperscriptBox["\[ExponentialE]", 
         RowBox[{"71", "/", "136"}]]}]]},
     {
      FractionBox[
       RowBox[{
        RowBox[{"-", "1"}], "-", 
        RowBox[{"2", " ", 
         SuperscriptBox["\[ExponentialE]", 
          RowBox[{"1", "/", "8"}]]}], "-", 
        RowBox[{"2", " ", 
         SuperscriptBox["\[ExponentialE]", 
          RowBox[{"3", "/", "8"}]]}], "+", 
        RowBox[{"2", " ", 
         SuperscriptBox["\[ExponentialE]", 
          RowBox[{"5", "/", "8"}]]}], "+", 
        RowBox[{"2", " ", 
         SuperscriptBox["\[ExponentialE]", 
          RowBox[{"7", "/", "8"}]]}], "+", "\[ExponentialE]"}], 
       RowBox[{"24", " ", 
        SuperscriptBox["\[ExponentialE]", 
         RowBox[{"71", "/", "136"}]]}]], 
      FractionBox[
       RowBox[{
        SuperscriptBox[
         RowBox[{"(", 
          RowBox[{"1", "+", 
           SuperscriptBox["\[ExponentialE]", 
            RowBox[{"1", "/", "4"}]]}], ")"}], "2"], " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", 
          RowBox[{"2", " ", 
           SuperscriptBox["\[ExponentialE]", 
            RowBox[{"1", "/", "8"}]]}], "+", 
          RowBox[{"2", " ", 
           SuperscriptBox["\[ExponentialE]", 
            RowBox[{"3", "/", "8"}]]}], "+", 
          SqrtBox["\[ExponentialE]"]}], ")"}]}], 
       RowBox[{"24", " ", 
        SuperscriptBox["\[ExponentialE]", 
         RowBox[{"71", "/", "136"}]]}]], 
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
        SuperscriptBox["\[ExponentialE]", 
         RowBox[{"71", "/", "136"}]]}]]},
     {
      FractionBox[
       RowBox[{
        SuperscriptBox[
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "+", 
           SuperscriptBox["\[ExponentialE]", 
            RowBox[{"1", "/", "4"}]]}], ")"}], "2"], " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", 
          SuperscriptBox["\[ExponentialE]", 
           RowBox[{"1", "/", "4"}]]}], ")"}], " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", 
          SuperscriptBox["\[ExponentialE]", 
           RowBox[{"1", "/", "8"}]], "+", 
          SuperscriptBox["\[ExponentialE]", 
           RowBox[{"1", "/", "4"}]]}], ")"}]}], 
       RowBox[{"24", " ", 
        SuperscriptBox["\[ExponentialE]", 
         RowBox[{"71", "/", "136"}]]}]], 
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
        SuperscriptBox["\[ExponentialE]", 
         RowBox[{"71", "/", "136"}]]}]], 
      FractionBox[
       RowBox[{
        SuperscriptBox[
         RowBox[{"(", 
          RowBox[{"1", "+", 
           SuperscriptBox["\[ExponentialE]", 
            RowBox[{"1", "/", "4"}]]}], ")"}], "2"], " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", 
          RowBox[{"2", " ", 
           SuperscriptBox["\[ExponentialE]", 
            RowBox[{"1", "/", "8"}]]}], "+", 
          RowBox[{"2", " ", 
           SuperscriptBox["\[ExponentialE]", 
            RowBox[{"3", "/", "8"}]]}], "+", 
          SqrtBox["\[ExponentialE]"]}], ")"}]}], 
       RowBox[{"24", " ", 
        SuperscriptBox["\[ExponentialE]", 
         RowBox[{"71", "/", "136"}]]}]]}
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
     {"1.009131490858961`", "0.1254885300705579`", "0.015564813186582259`"},
     {"0.1254885300705579`", "1.009131490858961`", "0.0019355312417656271`"},
     {"0.015564813186582259`", "0.0019355312417656271`", "1.009131490858961`"}
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
       SubscriptBox["k", "val"]}], "]"}], "]"}]}], "]"}]}]], "Input"],

Cell[BoxData["9.660797011037558`*^-16"], "Output"]
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
     RowBox[{"N", "[", 
      RowBox[{
       SubscriptBox["exp", "\[Tau]k"], ",", 
       RowBox[{"6", "MachinePrecision"}]}], "]"}]}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
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
     {"1.421542046212384685561292583225173572770467647961668121311899137988756\
09389529713356799147031142061788789136739082324373590958`94.22238864282609*^\
10", "8.8333680084578944239233362171637990573411525315659884765774868202672270\
048731930014655331995715711190299774031083886511044035`94.22238864282606*^9", 
      "1.708362340139251092483837320261995101687187263264741861278603529951865\
31657806452297552359781970420725712709835174207174137958`94.2223886428261*^\
10"},
     {"5.008460210565839468128739933714955455201387101749541194541850304102853\
3893403623007285642106778303794284887167683717726190711`94.22238864282605*^9",
       "3.11492720587853005451049441108167952291869063808116745064885814626258\
98978219308676937814124737917444567599638019895978948353`94.22238864282603*^\
9", "6.05240894320384300339582948971723743313493072343938073023333575937060867\
77639076956649021094420298217984721006356460815111392`94.2223886428261*^9"},
     {"4.304322135476647830065088497562125577210702107735104940142134270165736\
246547732562763518337506715744575281324168856860054466`94.22238864282606*^9", 
      "2.679233241071462999836542790122305207588196619251302616731522709434637\
194284939722931311040944728805560111922909903649144508`94.2223886428261*^9", 
      "5.295476275970609990721751697904455183420519997511328765413354874343436\
6024475270664122940057779537216036506431298376923600311`94.22238864282605*^9"}
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
  "1.3541106044926009`*^11", ",", "2.458359081039403`*^9", ",", 
   "4.6413000667705595`*^8", ",", "1.883262362055223`*^7", ",", 
   "1.2662520234656553`*^6", ",", "57405.00927188725`", ",", 
   "41481.02806470855`", ",", "9416.86844321774`", ",", "4837.700253605513`", 
   ",", "803.3948997938639`", ",", "253.53367072371103`", ",", 
   "8.292639705347293`", ",", "1.1406879772355831`", ",", 
   "0.10810471238762913`", ",", "0.0449400447104714`", ",", 
   "0.011087183907893695`", ",", "0.0009355293018800415`", ",", 
   "0.00013214568313003145`", ",", "0.00002216117309065215`", ",", 
   "4.04359758153946`*^-6", ",", "1.1006417242063872`*^-6", ",", 
   "7.707246723956613`*^-9", ",", "1.888792246975202`*^-9", ",", 
   "1.037217196211033`*^-11"}], "}"}]], "Output"],

Cell[BoxData["1.3055227096496119`*^22"], "Output"]
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
  "1.3541106045004163`*^11", ",", "2.4583590816612844`*^9", ",", 
   "4.6413000726669616`*^8", ",", "1.8832624189158835`*^7", ",", 
   "1.266252475024012`*^6", ",", "57405.22810876017`", ",", 
   "41481.367870834496`", ",", "9417.120799548033`", ",", 
   "4838.137017113882`", ",", "803.558923668733`", ",", "254.02245045695975`",
    ",", "8.912921422575774`", ",", "1.8934311456322133`", ",", 
   "1.0047859312459433`", ",", "0.9930448160290607`", ",", 
   "0.9792263227710335`", ",", "0.9669536200190407`", ",", 
   "0.9439711271100101`", ",", "0.903849382392107`", ",", 
   "0.8775969880013166`", ",", "0.8427302679429521`", ",", 
   "0.7824801940999317`", ",", "0.5063060796212253`", ",", 
   "0.4403818721344298`"}], "}"}]], "Output"],

Cell[BoxData["3.074855461096419`*^11"], "Output"]
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
     {"0.110454801252729635124453827097237832300413844301692314753970089655009\
0567802177752669455790065780395194945222451923`94.22238864282596", 
      RowBox[{
      "-", "0.2943851828752893153368483449079830490308365314826317575072946329\
592971718124302866087754637852247069857366151581274`94.22238864282596"}], 
      "0.109812836059017577066617449449904390084038730141099574946366480947130\
3325100582092031267875503509622568902703089735`94.22238864282596"},
     {
      RowBox[{
      "-", "0.1568532776771268862503084155770420801712975602402123606099206326\
516883820315434990214295716205126854963945865170711`94.22238864282596"}], 
      "0.534829119217626123562513039353826125180350245364197268400762009649169\
81250016682521998111995137897132988980840237`94.22238864282596", 
      "0.171491338063785310986467145599210189863356228586772329590387235977351\
3357653162057696793015007355417909595881092824`94.22238864282596"},
     {"0.016769601176051841402560695506870303083172646218401022029498864971654\
6395189853898623840599582558181462143793584763`94.22238864282596", 
      "0.080131547253846174268546740402007489819856469787457274784022283916032\
0852226557444002457065583367961881591737166698`94.22238864282596", 
      "0.504705621223911738133429773185043721026016763397913303105719744858123\
0883940206550705240506462269407030416030975428`94.22238864282596"}
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
  "2.2707565031076995`", ",", "1.9750898522650846`", ",", 
   "1.277987618779637`", ",", "1.1866192992461666`", ",", 
   "1.1394751961004905`", ",", "1.1063790267283504`", ",", 
   "1.059354434983116`", ",", "1.0341757653074495`", ",", 
   "1.0212143778673972`", ",", "1.0070038973656308`", ",", 
   "0.9952368647916788`", ",", "0.5281417295299121`", ",", 
   "0.11219665837814675`", ",", "0.00393665992199156`", ",", 
   "0.0012444638103630218`", ",", "0.00020669112851965795`", ",", 
   "0.00010618956911416006`", ",", "0.000024107208882643887`", ",", 
   "0.000017420016136951776`", ",", "7.897319213382284`*^-7", ",", 
   "5.309934451809742`*^-8", ",", "2.15456872932886`*^-9", ",", 
   "4.06775400493662`*^-10", ",", "7.384921118529594`*^-12"}], 
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

Cell[BoxData["0.9523782319980738`"], "Output"],

Cell[BoxData[
 RowBox[{"-", "1.312611538052992`"}]], "Output"],

Cell[BoxData["0.00004819595202010882`"], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "determinant", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   SubscriptBox["detG", "val"], "=", 
   RowBox[{"Det", "[", 
    SubscriptBox["G", "val"], "]"}]}]}]], "Input"],

Cell[BoxData["6.\
982066207297022680766427358445101294991056327999658560921181914033903427601644\
578539565151949485784360882438769758765707876056`94.22238864282596*^-66"], \
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
Cell[4750, 154, 1075, 30, 52, "Input"],
Cell[CellGroupData[{
Cell[5850, 188, 333, 10, 72, "Input"],
Cell[6186, 200, 28, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[6251, 205, 355, 10, 52, "Input"],
Cell[6609, 217, 28, 0, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[6686, 223, 45, 0, 43, "Subsection"],
Cell[6734, 225, 228, 7, 67, "Input"],
Cell[6965, 234, 258, 8, 67, "Input"],
Cell[7226, 244, 293, 10, 31, "Input"],
Cell[7522, 256, 492, 16, 52, "Input"],
Cell[CellGroupData[{
Cell[8039, 276, 442, 13, 52, "Input"],
Cell[8484, 291, 6548, 187, 140, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[15069, 483, 511, 15, 52, "Input"],
Cell[15583, 500, 786, 18, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[16406, 523, 441, 13, 52, "Input"],
Cell[16850, 538, 50, 0, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[16949, 544, 48, 0, 43, "Subsection"],
Cell[CellGroupData[{
Cell[17022, 548, 823, 26, 90, "Input"],
Cell[17848, 576, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[17960, 583, 567, 15, 92, "Input"],
Cell[18530, 600, 521, 14, 31, "Output"],
Cell[19054, 616, 521, 14, 31, "Output"],
Cell[19578, 632, 559, 16, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[20174, 653, 302, 10, 31, "Input"],
Cell[20479, 665, 8144, 105, 292, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[28672, 776, 45, 0, 43, "Subsection"],
Cell[28720, 778, 733, 19, 52, "Input"],
Cell[29456, 799, 124, 4, 31, "Input"],
Cell[CellGroupData[{
Cell[29605, 807, 448, 14, 52, "Input"],
Cell[30056, 823, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[30168, 830, 436, 13, 52, "Input"],
Cell[30607, 845, 1982, 34, 80, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[32626, 884, 392, 10, 72, "Input"],
Cell[33021, 896, 809, 14, 55, "Output"],
Cell[33833, 912, 50, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[33920, 917, 729, 19, 92, "Input"],
Cell[34652, 938, 780, 14, 55, "Output"],
Cell[35435, 954, 49, 0, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[35533, 960, 72, 0, 43, "Subsection"],
Cell[35608, 962, 253, 8, 31, "Input"],
Cell[CellGroupData[{
Cell[35886, 974, 436, 13, 52, "Input"],
Cell[36325, 989, 1932, 36, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[38294, 1030, 123, 3, 31, "Input"],
Cell[38420, 1035, 811, 15, 55, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[39268, 1055, 651, 18, 92, "Input"],
Cell[39922, 1075, 46, 0, 31, "Output"],
Cell[39971, 1077, 62, 1, 31, "Output"],
Cell[40036, 1080, 50, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[40123, 1085, 220, 6, 52, "Input"],
Cell[40346, 1093, 184, 3, 31, "Output"]
}, Open  ]],
Cell[40545, 1099, 1030, 29, 72, "Input"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature XwDp8ioheYlUsB1oVo3N5pAg *)
