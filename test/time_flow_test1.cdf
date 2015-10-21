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
NotebookDataLength[     24325,        727]
NotebookOptionsPosition[     22917,        656]
NotebookOutlinePosition[     23261,        671]
CellTagsIndexPosition[     23218,        668]
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
     {"62469.3035198537177778187469633115118678841583223966406830599848`31.\
90917954038203", 
      "163384.7772674144051182502135892017375114840063881828187686812209`31.\
90917954038201", 
      "80594.5149852625648018275253348419216992011205118460928600693882`31.\
90917954038201"},
     {"38895.8694613195893045154806402674350410469451055559312852949717`31.\
90917954038203", 
      "104728.0689226216077803175439270544005555125049907257136446144454`31.\
909179540382016", 
      "54386.8018285640573429594193892028894354279883148831064907971135`31.\
90917954038201"},
     {"18002.3010503149263486275748685662658909455612055797973989570565`31.\
90917954038201", 
      "49050.1594548777643145286964235107564689081040562435333533256288`31.\
90917954038201", 
      "29406.8207223453788235336962133874432466636739576087809081624233`31.\
909179540382006"}
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
  "462477.43356236775`", ",", "45505.83342751887`", ",", "33523.68315378342`",
    ",", "4575.600495988783`", ",", "1250.1053775355217`", ",", 
   "766.5069427755532`", ",", "287.4809714385927`", ",", "109.5182021977702`",
    ",", "77.26349764778679`", ",", "33.3461333373011`", ",", 
   "12.52972787256439`", ",", "8.983571880833011`", ",", "5.401174574048856`",
    ",", "1.2038566540779791`", ",", "0.3923606291870442`", ",", 
   "0.18787057560448653`", ",", 
   RowBox[{"0.077078228191939`", "\[VeryThinSpace]", "+", 
    RowBox[{"0.03097660441444115`", " ", "\[ImaginaryI]"}]}], ",", 
   RowBox[{"0.077078228191939`", "\[VeryThinSpace]", "-", 
    RowBox[{"0.03097660441444115`", " ", "\[ImaginaryI]"}]}], ",", 
   "0.029339898687228805`", ",", "0.0042357700507420345`", ",", 
   "0.001329026411395223`", ",", "0.0003388337252811161`", ",", 
   "0.00014398554364431542`", ",", "0.000014196588521782339`"}], 
  "}"}]], "Output"]
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
  "911568.0992787442`", ",", "79036.65823573529`", ",", "56423.33204385526`", 
   ",", "6684.893089403027`", ",", "1497.0834974282977`", ",", 
   "1038.265549380045`", ",", "381.811489885118`", ",", "172.1099196692093`", 
   ",", "91.4782837746991`", ",", "38.123860866148675`", ",", 
   "11.714360601163769`", ",", "8.032459205917194`", ",", 
   "3.7186938517401513`", ",", "0.9894977317946198`", ",", 
   "0.532580990743791`", ",", "0.0975491304002477`", ",", 
   "0.0908262436515306`", ",", "0.054775418889483195`", ",", 
   "0.029273299922101162`", ",", "0.0029218794802449506`", ",", 
   "0.0010766328042708239`", ",", "0.00020550427151156197`", ",", 
   "0.00007883990958102735`", ",", "8.267844439091696`*^-6"}], 
  "}"}]], "Output"],

Cell[BoxData["1.1025462634114206`*^11"], "Output"]
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
  "911568.6201054314`", ",", "79037.21486284667`", ",", "56423.78162585993`", 
   ",", "6685.326096309581`", ",", "1497.6007584334502`", ",", 
   "1038.647060912178`", ",", "382.4760600841942`", ",", "172.7241591612146`",
    ",", "92.05777667469071`", ",", "38.582006002997645`", ",", 
   "12.42232315455745`", ",", "8.740491790992069`", ",", 
   "4.2890363376276355`", ",", "1.7223361026864494`", ",", 
   "1.3058218510991617`", ",", "1.0196787289728164`", ",", 
   "0.9664311861997479`", ",", "0.9456394425894382`", ",", 
   "0.8690546220105698`", ",", "0.8279588436945453`", ",", 
   "0.7232170774781302`", ",", "0.661022081725821`", ",", 
   "0.5860882169540593`", ",", "0.5159214764282678`"}], "}"}]], "Output"],

Cell[BoxData["1.7668747314343934`*^6"], "Output"]
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
Cell[2379, 72, 339, 10, 52, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[2755, 87, 41, 0, 43, "Subsection"],
Cell[2799, 89, 320, 10, 72, "Input"],
Cell[CellGroupData[{
Cell[3144, 103, 947, 29, 72, "Input"],
Cell[4094, 134, 1469, 50, 52, "Output"],
Cell[5566, 186, 74, 2, 31, "Output"]
}, Open  ]],
Cell[5655, 191, 1075, 30, 52, "Input"],
Cell[CellGroupData[{
Cell[6755, 225, 333, 10, 72, "Input"],
Cell[7091, 237, 28, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[7156, 242, 355, 10, 52, "Input"],
Cell[7514, 254, 28, 0, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[7591, 260, 48, 0, 43, "Subsection"],
Cell[CellGroupData[{
Cell[7664, 264, 612, 20, 72, "Input"],
Cell[8279, 286, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[8391, 293, 567, 15, 92, "Input"],
Cell[8961, 310, 483, 12, 31, "Output"],
Cell[9447, 324, 506, 14, 31, "Output"],
Cell[9956, 340, 502, 13, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[10495, 358, 302, 10, 31, "Input"],
Cell[10800, 370, 4100, 54, 152, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[14949, 430, 45, 0, 43, "Subsection"],
Cell[14997, 432, 464, 12, 31, "Input"],
Cell[CellGroupData[{
Cell[15486, 448, 892, 27, 52, "Input"],
Cell[16381, 477, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[16493, 484, 436, 13, 52, "Input"],
Cell[16932, 499, 1410, 33, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[18379, 537, 117, 3, 31, "Input"],
Cell[18499, 542, 973, 17, 52, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[19509, 564, 392, 10, 72, "Input"],
Cell[19904, 576, 780, 14, 55, "Output"],
Cell[20687, 592, 50, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[20774, 597, 694, 19, 92, "Input"],
Cell[21471, 618, 757, 13, 52, "Output"],
Cell[22231, 633, 49, 0, 31, "Output"]
}, Open  ]],
Cell[22295, 636, 606, 17, 52, "Input"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature 0wpSTUSXkysT4C1rkWLzBhdr *)
