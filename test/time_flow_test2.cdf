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
NotebookDataLength[     28326,        797]
NotebookOptionsPosition[     26917,        726]
NotebookOutlinePosition[     27260,        741]
CellTagsIndexPosition[     27217,        738]
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
    SubscriptBox["\[Beta]", "val"], "=", 
    RowBox[{"2", "-", 
     FractionBox["1", "11"]}]}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"time", " ", "step"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[Tau]", "val"], "=", 
    FractionBox["1", "11"]}], ";"}]}]], "Input"],

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

Cell[BoxData["21"], "Output"]
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
    RowBox[{"3", "/", "4"}]}], ";"}]}]], "Input"]
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
     SubscriptBox["n", "x"], "=", "5"}], ";"}], "\[IndentingNewLine]", 
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
    RowBox[{"4", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"3", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"4", ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"3", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"4", ",", "2"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"3", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"4", ",", "3"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"3", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"4", ",", "4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "5"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", "5"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2", ",", "5"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"3", ",", "5"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"4", ",", "5"}], "}"}]}], "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"30", ",", "2"}], "}"}]], "Output"]
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
         RowBox[{
          SubscriptBox["n", "x"], 
          SubscriptBox["n", "y"]}], ",", 
         SubscriptBox["L", "val"]}], "}"}]}], "]"}]}], "-", "1"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"Dimensions", "[", 
  SubscriptBox["s", "val"], "]"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"30", ",", "21"}], "}"}]], "Output"]
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
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}]}], "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}]}], "}"}]], "Output"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1", ",", "1", ",", 
   RowBox[{"-", "1"}], ",", "1", ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", 
   RowBox[{"-", "1"}], ",", "1"}], "}"}]], "Output"]
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
  "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", 
   ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", ",", 
   "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", 
   ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", 
   "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", 
   ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", 
   "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", 
   ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", 
   "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", 
   ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", 
   ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", 
   "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", 
   ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", 
   "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", 
   ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", 
   "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", 
   ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", 
   "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", 
   ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", 
   "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", 
   ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", 
   "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", 
   ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", 
   "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", 
   ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", 
   "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", 
   ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", 
   "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", 
   ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", 
   "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", 
   ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", 
   ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", 
   "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", 
   ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", 
   "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", 
   ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", ",", "0", ",", 
   "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", 
   ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", 
   ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", 
   "0", ",", "0", ",", "0", ",", "0", ",", "0", ",", "0", ",", "0", ",", "0", 
   ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", 
   "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", 
   ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", 
   "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", 
   ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", 
   ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", 
   ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", 
   ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", 
   ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", 
   ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", 
   ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", 
   ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", 
   "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", 
   ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", 
   ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", 
   "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", 
   ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", 
   "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", 
   ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", 
   "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", 
   ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", 
   "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", 
   ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", ",", "0", ",", "0", 
   ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", 
   "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", 
   ",", "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", 
   ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", 
   "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", 
   ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", 
   "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", 
   ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0"}], 
  "}"}]], "Output"]
}, Open  ]]
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
       RowBox[{"2", "MachinePrecision"}]}], "}"}], ",", 
     RowBox[{"HubbardTimeFlowMap", "[", 
      RowBox[{
       RowBox[{"Exp", "[", 
        RowBox[{
         RowBox[{"-", 
          SubscriptBox["\[Lambda]", "val"]}], 
         SubscriptBox["s", "val"]}], "]"}], ",", 
       RowBox[{"MatrixExp", "[", 
        RowBox[{"N", "[", 
         RowBox[{
          RowBox[{
           RowBox[{"-", 
            SubscriptBox["\[Tau]", "val"]}], 
           RowBox[{"(", 
            RowBox[{"-", 
             SubscriptBox["latt", "neigh"]}], ")"}]}], ",", 
          RowBox[{"2", "MachinePrecision"}]}], "]"}], "]"}]}], "]"}]}], 
    "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{"Dimensions", "[", "%", "]"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"30", ",", "30"}], "}"}]], "Output"]
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
     {"130293.2910032659931862537578856038692378708540820301477482002703`31.\
90917954038202", 
      "334115.6980469842609229231428613519097584507911823223546685792964`31.\
909179540382016", 
      "150924.5913694387526562543988641304465115208063619999925593831477`31.\
90917954038201"},
     {"23651.6967132148898897739516537575081334278510968410473512397339`31.\
90917954038201", 
      "62229.6823200303734725002668100814286600267714914493150143531248`31.\
909179540382016", 
      "26961.1977973862409578654483084465549460533425767673281392760927`31.\
90917954038201"},
     {"8635.4411476779681452583724249652217880108073581359963226379701`31.\
909179540382006", 
      "21461.7114798753883783212240730240806011960216266097718025836594`31.\
90917954038201", 
      "11337.0455810316393014549914863421512369249570534653113072634631`31.\
90917954038201"}
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
  "903991.1127378811`", ",", "88967.4697270439`", ",", "51561.698238619945`", 
   ",", "25356.756284476745`", ",", "8614.94393181194`", ",", 
   "2271.2096757196505`", ",", "1133.8912362444498`", ",", 
   "755.2927527712854`", ",", "583.9048348921339`", ",", 
   "181.06704227276686`", ",", "61.004243190727095`", ",", 
   "20.11963002248016`", ",", "13.724686754926532`", ",", 
   "10.00142207871787`", ",", "3.3353028063712973`", ",", 
   "1.5341051904445109`", ",", "0.8776440647157964`", ",", 
   RowBox[{"0.4037721958508881`", "\[VeryThinSpace]", "+", 
    RowBox[{"0.026402808443626212`", " ", "\[ImaginaryI]"}]}], ",", 
   RowBox[{"0.4037721958508881`", "\[VeryThinSpace]", "-", 
    RowBox[{"0.026402808443626212`", " ", "\[ImaginaryI]"}]}], ",", 
   "0.16734964782610778`", ",", "0.07566201754468066`", ",", 
   "0.05727413410813141`", ",", "0.011663945873079867`", ",", 
   "0.006057878542715313`", ",", "0.0011066853450720326`", ",", 
   "0.0007039116844469646`", ",", "0.00015507422976824962`", ",", 
   "0.0000841072375772166`", ",", "8.345454230543897`*^-6", ",", 
   "5.306220318636983`*^-6"}], "}"}]], "Output"]
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
  "2.199822813518927`*^6", ",", "206505.85539426643`", ",", 
   "87091.90676552047`", ",", "40961.33697110747`", ",", 
   "12898.004702301327`", ",", "3223.9704112826166`", ",", 
   "1371.104337875546`", ",", "1043.941215951833`", ",", "759.3232352788617`",
    ",", "259.34943754686344`", ",", "63.40749465479545`", ",", 
   "40.12484811455317`", ",", "12.992963426464495`", ",", 
   "10.322001642744393`", ",", "3.6705751762532492`", ",", 
   "1.7611071903827522`", ",", "0.9203793252693295`", ",", 
   "0.3939604893910697`", ",", "0.2414087055193538`", ",", 
   "0.08494811037474824`", ",", "0.05457582869723254`", ",", 
   "0.03903963429417093`", ",", "0.007638837069163889`", ",", 
   "0.004530905277508469`", ",", "0.0008798108045062188`", ",", 
   "0.00037813249680478025`", ",", "0.00010636411796144558`", ",", 
   "0.0000351421068286678`", ",", "5.858262558159103`*^-6", ",", 
   "3.1790056191919035`*^-6"}], "}"}]], "Output"],

Cell[BoxData["6.919845627948645`*^11"], "Output"]
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
  "2.1998232366909743`*^6", ",", "206506.2308722211`", ",", 
   "87092.37156522639`", ",", "40961.93944488191`", ",", 
   "12898.471313997172`", ",", "3224.674285157619`", ",", 
   "1371.4704047633677`", ",", "1044.4209247168608`", ",", 
   "759.9032646227647`", ",", "259.89413381251717`", ",", 
   "64.06050574845946`", ",", "40.31219951078935`", ",", 
   "13.384276805853524`", ",", "11.000121047673296`", ",", 
   "4.233781731418773`", ",", "2.4088686126496324`", ",", 
   "1.592159166749276`", ",", "1.2316431097424234`", ",", 
   "1.0319840333785282`", ",", "1.0238242247551887`", ",", 
   "0.9590037406731784`", ",", "0.9294316249767453`", ",", 
   "0.9076877243123068`", ",", "0.8811717473628687`", ",", 
   "0.7893999476621284`", ",", "0.7142815531869668`", ",", 
   "0.6339986260692161`", ",", "0.572095706210354`", ",", 
   "0.3669412310822608`", ",", "0.2681498838905997`"}], "}"}]], "Output"],

Cell[BoxData["8.203707586121732`*^6"], "Output"]
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
       RowBox[{"NotebookFileName", "[", "]"}], "]"}], "<>", "\"\<.dat\>\""}], 
     ",", 
     RowBox[{"Flatten", "[", 
      RowBox[{"Transpose", "[", 
       RowBox[{"N", "[", 
        SubscriptBox["A", "val"], "]"}], "]"}], "]"}], ",", 
     "\"\<Real64\>\""}], "]"}], ";"}]}]], "Input"]
}, Open  ]]
},
WindowSize->{1350, 867},
WindowMargins->{{Automatic, 244}, {90, Automatic}},
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
Cell[1529, 37, 269, 9, 67, "Input"],
Cell[1801, 48, 229, 7, 67, "Input"],
Cell[CellGroupData[{
Cell[2055, 59, 295, 9, 69, "Input"],
Cell[2353, 70, 29, 0, 31, "Output"]
}, Open  ]],
Cell[2397, 73, 339, 10, 52, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[2773, 88, 41, 0, 43, "Subsection"],
Cell[2817, 90, 320, 10, 72, "Input"],
Cell[CellGroupData[{
Cell[3162, 104, 947, 29, 72, "Input"],
Cell[4112, 135, 1823, 62, 52, "Output"],
Cell[5938, 199, 74, 2, 31, "Output"]
}, Open  ]],
Cell[6027, 204, 1075, 30, 52, "Input"],
Cell[CellGroupData[{
Cell[7127, 238, 333, 10, 72, "Input"],
Cell[7463, 250, 28, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[7528, 255, 355, 10, 52, "Input"],
Cell[7886, 267, 28, 0, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[7963, 273, 48, 0, 43, "Subsection"],
Cell[CellGroupData[{
Cell[8036, 277, 612, 20, 72, "Input"],
Cell[8651, 299, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[8763, 306, 567, 15, 92, "Input"],
Cell[9333, 323, 581, 14, 31, "Output"],
Cell[9917, 339, 638, 17, 31, "Output"],
Cell[10558, 358, 733, 22, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[11328, 385, 302, 10, 31, "Input"],
Cell[11633, 397, 6691, 87, 252, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[18373, 490, 45, 0, 43, "Subsection"],
Cell[18421, 492, 464, 12, 31, "Input"],
Cell[CellGroupData[{
Cell[18910, 508, 892, 27, 52, "Input"],
Cell[19805, 537, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[19917, 544, 436, 13, 52, "Input"],
Cell[20356, 559, 1411, 33, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[21804, 597, 117, 3, 31, "Input"],
Cell[21924, 602, 1167, 20, 76, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[23128, 627, 392, 10, 72, "Input"],
Cell[23523, 639, 975, 17, 55, "Output"],
Cell[24501, 658, 49, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[24587, 663, 694, 19, 92, "Input"],
Cell[25284, 684, 945, 17, 55, "Output"],
Cell[26232, 703, 48, 0, 31, "Output"]
}, Open  ]],
Cell[26295, 706, 606, 17, 52, "Input"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature 4wDNMNEy7G0gLCwkiEysC48b *)
