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
NotebookDataLength[     28731,        871]
NotebookOptionsPosition[     26831,        783]
NotebookOutlinePosition[     27175,        798]
CellTagsIndexPosition[     27132,        795]
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
    RowBox[{"\[Lambda]", " ", "parameter", " ", "for", " ", "Hubbard"}], "-", 
    RowBox[{"Stratonovich", " ", "field"}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[Lambda]", "val"], "=", 
    RowBox[{"3", "/", "4"}]}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"electron", "-", 
    RowBox[{"phonon", " ", "interaction", " ", "strength"}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["g", "val"], "=", 
    RowBox[{"7", "/", "10"}]}], ";"}]}]], "Input"]
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
  SubscriptBox["n", "sites"], "=", 
  RowBox[{"Length", "[", 
   SubscriptBox["latt", "sites"], "]"}]}]], "Input"],

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
         SubscriptBox["n", "sites"], ",", 
         SubscriptBox["L", "val"]}], "}"}]}], "]"}]}], "-", "1"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"Dimensions", "[", 
  SubscriptBox["s", "val"], "]"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"24", ",", "16"}], "}"}]], "Output"]
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
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}]}], "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1"}], "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", 
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

Cell["Phonon field", "Subsection"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"SeedRandom", "[", "42", "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   SubscriptBox["X", "val"], "=", 
   RowBox[{
    FractionBox["1", "8"], 
    RowBox[{"RandomReal", "[", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{
        RowBox[{"-", "32"}], ",", "32"}], "}"}], ",", 
      RowBox[{"{", 
       RowBox[{
        SubscriptBox["n", "sites"], ",", 
        SubscriptBox["L", "val"]}], "}"}], ",", 
      RowBox[{"WorkingPrecision", "\[Rule]", 
       RowBox[{"4", "MachinePrecision"}]}]}], "]"}]}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"Dimensions", "[", "%", "]"}]}], "Input"],

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
    RowBox[{
    "0.5933178240671829`", ",", "2.54002679332862`", ",", 
     "3.547967343509361`"}], "}"}]}], "}"}]], "Output"]
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
     "\"\<Real64\>\""}], "]"}], ";"}]}]], "Input"]
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
           RowBox[{"(", 
            RowBox[{"-", 
             SubscriptBox["latt", "neigh"]}], ")"}]}], ",", 
          RowBox[{"4", "MachinePrecision"}]}], "]"}], "]"}], ",", 
       RowBox[{"Exp", "[", 
        RowBox[{
         RowBox[{
          RowBox[{"-", 
           SubscriptBox["\[Lambda]", "val"]}], 
          SubscriptBox["s", "val"]}], "-", 
         RowBox[{
          SubscriptBox["\[Tau]", "val"], " ", 
          SubscriptBox["g", "val"], 
          SubscriptBox["X", "val"]}]}], "]"}]}], "]"}]}], "]"}]}], 
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
     {"89162.\
980987753465854297890330575311996294258472534018476283679398533737785112680602\
760000095107856099`63.818359080764054", 
      "218749.0268921007472277650602814117026108464917318730358557744145731972\
15844643088314424418975100192469`63.81835908076401", 
      "86542.54180333742608755973703443240802262254139781661315145845511014688\
1838661808343305436243236382034`63.818359080764026"},
     {"53006.\
354968884248270494848693491595425412336303086797657209876978649893819195584927\
786994657752645427`63.818359080764075", 
      "132077.6533850949892412603237515221295335785785513915366473887099200729\
27440615040120844954339440955797`63.81835908076402", 
      "55713.48411160466074440862923657296597708682304039631248208927115241387\
2500301737609738652845832600426`63.81835908076401"},
     {"19875.\
763956333254319931694191113812868405810013894282266370226537150340545261198851\
559516918588653035`63.818359080764026", 
      "48802.36648513267187851368730155380735562895112901918079661736922016818\
8715836082996015899169647381099`63.81835908076401", 
      "26106.95198982414389132266217533342941252071338654663320948625374466720\
3712359550633592893867490680792`63.81835908076401"}
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
  RowBox[{"735620.6892761089`", ",", 
   RowBox[{"61882.56010396028`", "\[VeryThinSpace]", "+", 
    RowBox[{"19745.514003939716`", " ", "\[ImaginaryI]"}]}], ",", 
   RowBox[{"61882.56010396028`", "\[VeryThinSpace]", "-", 
    RowBox[{"19745.514003939716`", " ", "\[ImaginaryI]"}]}], ",", 
   "5275.964270669338`", ",", "1202.4967955957043`", ",", 
   "871.6137049903214`", ",", "357.9600834432995`", ",", "137.2022522551753`",
    ",", "107.66932104184818`", ",", "38.0156652985503`", ",", 
   "16.932114085947195`", ",", 
   RowBox[{"9.16855098609305`", "\[VeryThinSpace]", "+", 
    RowBox[{"2.0221377006875065`", " ", "\[ImaginaryI]"}]}], ",", 
   RowBox[{"9.16855098609305`", "\[VeryThinSpace]", "-", 
    RowBox[{"2.0221377006875065`", " ", "\[ImaginaryI]"}]}], ",", 
   "1.3418288925521435`", ",", "0.2551336596546565`", ",", 
   "0.1889391270539019`", ",", "0.0932874689237606`", ",", 
   RowBox[{"0.04444015680567801`", "\[VeryThinSpace]", "+", 
    RowBox[{"0.0064920084667040925`", " ", "\[ImaginaryI]"}]}], ",", 
   RowBox[{"0.04444015680567801`", "\[VeryThinSpace]", "-", 
    RowBox[{"0.0064920084667040925`", " ", "\[ImaginaryI]"}]}], ",", 
   "0.004504402085512505`", ",", "0.0019306475026127294`", ",", 
   "0.00034432879643612406`", ",", "0.0001242658498279786`", ",", 
   "0.000012149764515330047`"}], "}"}]], "Output"]
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
  "1.4149734751049972`*^6", ",", "144711.86935971194`", ",", 
   "108128.4694784551`", ",", "7662.259717482019`", ",", "1906.502348152517`",
    ",", "1021.3693046373086`", ",", "518.4627467667169`", ",", 
   "190.62553356177415`", ",", "116.80705359880444`", ",", "56.466783998658`",
    ",", "15.764815485995015`", ",", "7.382167766255946`", ",", 
   "6.18661441375155`", ",", "1.1501954356072281`", ",", 
   "0.3699215989023006`", ",", "0.09880138563548278`", ",", 
   "0.07176551851156257`", ",", "0.048299740465296226`", ",", 
   "0.033168256618267354`", ",", "0.0025734848497138696`", ",", 
   "0.0015031631523402175`", ",", "0.000205666965447385`", ",", 
   "0.00006857538759696212`", ",", "6.9313733629869975`*^-6"}], 
  "}"}]], "Output"],

Cell[BoxData["2.041404208090776`*^11"], "Output"]
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
  "1.4149740011826477`*^6", ",", "144712.43113509758`", ",", 
   "108128.80632137397`", ",", "7662.714914146222`", ",", 
   "1906.9235197760102`", ",", "1021.8994471127148`", ",", 
   "519.0162505175642`", ",", "191.03703039534008`", ",", 
   "117.41183657724292`", ",", "56.95510315331998`", ",", 
   "16.525094296933148`", ",", "7.88414574744706`", ",", "6.867408668802429`",
    ",", "1.8481803622832536`", ",", "1.201111461803738`", ",", 
   "1.0144055935693068`", ",", "0.9674351003111376`", ",", 
   "0.9158794864838261`", ",", "0.8863630818439282`", ",", 
   "0.8192474626183507`", ",", "0.7051003257261932`", ",", 
   "0.6339246768761118`", ",", "0.516138908104978`", ",", 
   "0.4464252229150884`"}], "}"}]], "Output"],

Cell[BoxData["3.169565536514904`*^6"], "Output"]
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
      RowBox[{"NotebookDirectory", "[", "]"}], "<>", 
      RowBox[{"FileBaseName", "[", 
       RowBox[{"NotebookFileName", "[", "]"}], "]"}], "<>", 
      "\"\<_A.dat\>\""}], ",", 
     RowBox[{"Flatten", "[", 
      RowBox[{"Transpose", "[", 
       RowBox[{"N", "[", 
        SubscriptBox["A", "val"], "]"}], "]"}], "]"}], ",", 
     "\"\<Real64\>\""}], "]"}], ";"}]}]], "Input"]
}, Open  ]]
},
WindowSize->{1509, 867},
WindowMargins->{{Automatic, 239}, {108, Automatic}},
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
Cell[2379, 72, 339, 10, 52, "Input"],
Cell[2721, 84, 287, 9, 52, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[3045, 98, 41, 0, 43, "Subsection"],
Cell[3089, 100, 320, 10, 72, "Input"],
Cell[CellGroupData[{
Cell[3434, 114, 947, 29, 72, "Input"],
Cell[4384, 145, 1469, 50, 52, "Output"],
Cell[5856, 197, 74, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[5967, 204, 138, 4, 31, "Input"],
Cell[6108, 210, 29, 0, 31, "Output"]
}, Open  ]],
Cell[6152, 213, 1075, 30, 52, "Input"],
Cell[CellGroupData[{
Cell[7252, 247, 333, 10, 72, "Input"],
Cell[7588, 259, 28, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[7653, 264, 355, 10, 52, "Input"],
Cell[8011, 276, 28, 0, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[8088, 282, 48, 0, 43, "Subsection"],
Cell[CellGroupData[{
Cell[8161, 286, 560, 18, 72, "Input"],
Cell[8724, 306, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[8836, 313, 567, 15, 92, "Input"],
Cell[9406, 330, 483, 12, 31, "Output"],
Cell[9892, 344, 506, 14, 31, "Output"],
Cell[10401, 360, 502, 13, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[10940, 378, 302, 10, 31, "Input"],
Cell[11245, 390, 4100, 54, 152, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[15394, 450, 34, 0, 43, "Subsection"],
Cell[CellGroupData[{
Cell[15453, 454, 645, 20, 88, "Input"],
Cell[16101, 476, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[16213, 483, 431, 13, 52, "Input"],
Cell[16647, 498, 466, 14, 31, "Output"]
}, Open  ]],
Cell[17128, 515, 609, 17, 52, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[17774, 537, 45, 0, 43, "Subsection"],
Cell[17822, 539, 625, 17, 52, "Input"],
Cell[CellGroupData[{
Cell[18472, 560, 1059, 32, 52, "Input"],
Cell[19534, 594, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[19646, 601, 436, 13, 52, "Input"],
Cell[20085, 616, 1760, 36, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[21882, 657, 117, 3, 31, "Input"],
Cell[22002, 662, 1367, 23, 52, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[23406, 690, 392, 10, 72, "Input"],
Cell[23801, 702, 786, 14, 55, "Output"],
Cell[24590, 718, 49, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[24676, 723, 694, 19, 92, "Input"],
Cell[25373, 744, 767, 14, 55, "Output"],
Cell[26143, 760, 48, 0, 31, "Output"]
}, Open  ]],
Cell[26206, 763, 609, 17, 52, "Input"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature 5vp4L2@8wFm2QCgx74etI@WQ *)
