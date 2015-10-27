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
NotebookDataLength[     69921,       2031]
NotebookOptionsPosition[     65895,       1874]
NotebookOutlinePosition[     66238,       1889]
CellTagsIndexPosition[     66195,       1886]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["General parameters", "Subsection"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"Coulomb", " ", "coupling", " ", "constant"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["U", "val"], "=", 
     RowBox[{"9", "/", "2"}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"N", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData["4.5`"], "Output"]
}, Open  ]],

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

Cell[BoxData[{
 RowBox[{
  RowBox[{
   SubscriptBox["\[Lambda]", "val"], "=", 
   RowBox[{"ArcCosh", "[", 
    RowBox[{"Exp", "[", 
     RowBox[{
      SubscriptBox["\[Tau]", "val"], 
      RowBox[{
       SubscriptBox["U", "val"], "/", "2"}]}], "]"}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"N", "[", "%", "]"}]}], "Input"],

Cell[BoxData["0.7856003600934582`"], "Output"]
}, Open  ]],

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
   RowBox[{"phonon", " ", "frequency"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[CapitalOmega]", "val"], "=", 
    RowBox[{"13", "/", "10"}]}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"electron", "-", 
    RowBox[{"phonon", " ", "interaction", " ", "strength"}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["g", "val"], "=", 
    RowBox[{"7", "/", "10"}]}], ";"}]}]], "Input"],

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
   RowBox[{"nearest", " ", 
    RowBox[{"(", "direct", ")"}], " ", "neighbors"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["neigh", "nearest"], "=", 
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
    SubscriptBox["neigh", "nearest"], ";"}], "\[IndentingNewLine]", 
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
      SubscriptBox["neigh", "nearest"]}], ")"}], "-", "4"}], "]"}]}]], "Input"],

Cell[BoxData["0"], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"next", "-", 
    RowBox[{"nearest", " ", 
     RowBox[{"(", "diagonal", ")"}], " ", "neighbors"}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["neigh", "nextnearest"], "=", 
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
               RowBox[{"-", "1"}]}], "]"}]}], "}"}], "]"}], "\[Equal]", 
          SqrtBox["2"]}], ",", "1", ",", "0"}], "]"}], "&"}], ",", 
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
    SubscriptBox["neigh", "nextnearest"], ";"}], "\[IndentingNewLine]", 
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
    RowBox[{
     RowBox[{"every", " ", "site", " ", "has", " ", "4", " ", "next"}], "-", 
     RowBox[{"nearest", " ", "neighbors"}]}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{"Norm", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{"Total", "/@", 
      SubscriptBox["neigh", "nextnearest"]}], ")"}], "-", "4"}], 
   "]"}]}]], "Input"],

Cell[BoxData["0"], "Output"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Kinetic energy operator", "Subsection"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"t", " ", "hopping", " ", "parameter"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["t", "val"], "=", "1"}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    SuperscriptBox["t", "\[Prime]",
     MultilineFunction->None], " ", "hopping", " ", "parameter"}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["tp", "val"], "=", 
    RowBox[{"-", 
     FractionBox["1", "12"]}]}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"chemical", " ", "potential"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[Mu]", "val"], "=", 
    RowBox[{"-", 
     FractionBox["2", "17"]}]}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["k", "val"], "=", 
   RowBox[{"-", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{
       SubscriptBox["t", "val"], 
       SubscriptBox["neigh", "nearest"]}], "+", 
      RowBox[{
       SubscriptBox["tp", "val"], 
       SubscriptBox["neigh", "nextnearest"]}], "+", 
      RowBox[{
       SubscriptBox["\[Mu]", "val"], 
       RowBox[{"IdentityMatrix", "[", 
        SubscriptBox["n", "sites"], "]"}]}]}], ")"}]}]}], ";"}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   SuperscriptBox["\[ExponentialE]", 
    RowBox[{
     RowBox[{"-", "\[Tau]"}], " ", "k"}]], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["exp", "\[Tau]k"], "=", 
    RowBox[{"MatrixExp", "[", 
     RowBox[{"N", "[", 
      RowBox[{
       RowBox[{
        RowBox[{"-", 
         SubscriptBox["\[Tau]", "val"]}], 
        SubscriptBox["k", "val"]}], ",", 
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
     {"1.0161512421630685`", "0.12380186065442117`", "0.004851180264200935`"},
     {"0.12380186065442117`", "1.0161512421630685`", 
      RowBox[{"-", "0.0006630341916546836`"}]},
     {"0.004851180264200935`", 
      RowBox[{"-", "0.0006630341916546836`"}], "1.0161512421630685`"}
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
    RowBox[{"\[Tau]", " ", "k"}]], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["invexp", "\[Tau]k"], "=", 
    RowBox[{"MatrixExp", "[", 
     RowBox[{"N", "[", 
      RowBox[{
       RowBox[{
        SubscriptBox["\[Tau]", "val"], 
        SubscriptBox["k", "val"]}], ",", 
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
     SubscriptBox["n", "sites"], "]"}]}], "]"}]}]], "Input"],

Cell[BoxData["0``63.41384730152976"], "Output"]
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
   SubscriptBox["s", "0"], "=", 
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
  SubscriptBox["s", "0"], "]"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"24", ",", "16"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   SubscriptBox["s", "0"], "\[LeftDoubleBracket]", 
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
    RowBox[{
     RowBox[{"-", "1"}], ",", 
     RowBox[{"-", "1"}], ",", 
     RowBox[{"-", "1"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "1"}], ",", "1", ",", 
     RowBox[{"-", "1"}]}], "}"}]}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Flatten", "[", 
  RowBox[{
   RowBox[{
    RowBox[{"Transpose", "[", 
     SubscriptBox["s", "0"], "]"}], "/.", 
    RowBox[{"{", 
     RowBox[{"1", "\[Rule]", "0"}], "}"}]}], "/.", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "1"}], "\[Rule]", "1"}], "}"}]}], "]"}]], "Input"],

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
   SubscriptBox["X", "0"], "=", 
   RowBox[{"RandomReal", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
       RowBox[{"-", "4"}], ",", "4"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
       SubscriptBox["n", "sites"], ",", 
       SubscriptBox["L", "val"]}], "}"}], ",", 
     RowBox[{"WorkingPrecision", "\[Rule]", 
      RowBox[{"4", "MachinePrecision"}]}]}], "]"}]}], 
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
    "0.5933178240671829`", ",", "2.54002679332862`", ",", 
     "3.547967343509361`"}], "}"}]}], "}"}]], "Output"]
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
}, Open  ]]
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
   "expK_", ",", "invexpK_", ",", "\[Tau]_", ",", "\[Lambda]_", ",", 
    "\[CapitalOmega]_", ",", "g_", ",", "boxwidth_", ",", "s0_", ",", "X0_", 
    ",", "siteorderHS_List", ",", "siteorderPh_List", ",", "rnd_List"}], 
   "]"}], ":=", 
  RowBox[{"Module", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{
      RowBox[{"s", "=", "s0"}], ",", 
      RowBox[{"X", "=", "X0"}], ",", "Gu", ",", "Gd", ",", "nsites", ",", "L",
       ",", "\[CapitalDelta]x", ",", "\[CapitalDelta]Eph", ",", "du", ",", 
      "dd", ",", "cu", ",", "cd", ",", "l", ",", "i", ",", "j", ",", 
      RowBox[{"n", "=", "1"}]}], "}"}], ",", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{
      RowBox[{"{", 
       RowBox[{"nsites", ",", "L"}], "}"}], "=", 
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
        RowBox[{"\[Lambda]", " ", "s"}], ",", 
        RowBox[{"\[Tau]", " ", "g", " ", "X0"}]}], "]"}]}], ";", 
     "\[IndentingNewLine]", 
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
          RowBox[{"j", "\[LessEqual]", "nsites"}], ",", 
          RowBox[{"j", "++"}], ",", "\[IndentingNewLine]", 
          RowBox[{
           RowBox[{"i", "=", 
            RowBox[{"siteorderHS", "\[LeftDoubleBracket]", 
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
                  RowBox[{"+", "2"}], "\[Lambda]", " ", 
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
                  RowBox[{"-", "2"}], "\[Lambda]", " ", 
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
                    RowBox[{"+", "2"}], "\[Lambda]", " ", 
                    RowBox[{"s", "\[LeftDoubleBracket]", 
                    RowBox[{"i", ",", "l"}], "\[RightDoubleBracket]"}]}], 
                   "]"}], "-", "1"}], ")"}], 
                RowBox[{"(", 
                 RowBox[{
                  RowBox[{"UnitVector", "[", 
                   RowBox[{"nsites", ",", "i"}], "]"}], "-", 
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
                    RowBox[{"-", "2"}], "\[Lambda]", " ", 
                    RowBox[{"s", "\[LeftDoubleBracket]", 
                    RowBox[{"i", ",", "l"}], "\[RightDoubleBracket]"}]}], 
                   "]"}], "-", "1"}], ")"}], 
                RowBox[{"(", 
                 RowBox[{
                  RowBox[{"UnitVector", "[", 
                   RowBox[{"nsites", ",", "i"}], "]"}], "-", 
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
          RowBox[{"j", "\[LessEqual]", "nsites"}], ",", 
          RowBox[{"j", "++"}], ",", "\[IndentingNewLine]", 
          RowBox[{
           RowBox[{"i", "=", 
            RowBox[{"siteorderPh", "\[LeftDoubleBracket]", 
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
                SuperscriptBox["\[CapitalOmega]", "2"], 
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
                  RowBox[{"-", "\[Tau]"}], " ", "g", " ", 
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
                  RowBox[{"-", "\[Tau]"}], " ", "g", " ", 
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
                    RowBox[{"-", "\[Tau]"}], " ", "g", " ", 
                    "\[CapitalDelta]x"}], "]"}], "-", "1"}], ")"}], 
                RowBox[{"(", 
                 RowBox[{
                  RowBox[{"UnitVector", "[", 
                   RowBox[{"nsites", ",", "i"}], "]"}], "-", 
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
                    RowBox[{"-", "\[Tau]"}], " ", "g", " ", 
                    "\[CapitalDelta]x"}], "]"}], "-", "1"}], ")"}], 
                RowBox[{"(", 
                 RowBox[{
                  RowBox[{"UnitVector", "[", 
                   RowBox[{"nsites", ",", "i"}], "]"}], "-", 
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

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["fn", "base"], "=", 
   RowBox[{
    RowBox[{"NotebookDirectory", "[", "]"}], "<>", 
    RowBox[{"FileBaseName", "[", 
     RowBox[{"NotebookFileName", "[", "]"}], "]"}]}]}], ";"}]], "Input"],

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

Cell[BoxData["1152"], "Output"]
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
  "0.8915577426527812`", ",", "0.27670149969199476`", ",", 
   "0.3825412917734885`"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"pre", "-", 
    RowBox[{
    "determined", " ", "site", " ", "update", " ", "ordering", " ", "for", 
     " ", "the", " ", "Hubbard"}], "-", 
    RowBox[{"Stratonovich", " ", "field"}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["site", 
      RowBox[{"order", ",", "HS"}]], "=", 
     RowBox[{"Partition", "[", 
      RowBox[{
       RowBox[{
        RowBox[{"Import", "[", 
         RowBox[{
          RowBox[{
           SubscriptBox["fn", "base"], "<>", "\"\<_siteorderHS.dat\>\""}], 
          ",", "\"\<Integer32\>\""}], "]"}], "+", "1"}], ",", 
       SubscriptBox["n", "sites"]}], "]"}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"Dimensions", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"16", ",", "24"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"pre", "-", 
    RowBox[{
    "determined", " ", "site", " ", "update", " ", "ordering", " ", "for", 
     " ", "the", " ", "phonon", " ", "field"}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["site", 
      RowBox[{"order", ",", "ph"}]], "=", 
     RowBox[{"Partition", "[", 
      RowBox[{
       RowBox[{
        RowBox[{"Import", "[", 
         RowBox[{
          RowBox[{
           SubscriptBox["fn", "base"], "<>", "\"\<_siteorderPh.dat\>\""}], 
          ",", "\"\<Integer32\>\""}], "]"}], "+", "1"}], ",", 
       SubscriptBox["n", "sites"]}], "]"}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"Dimensions", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"16", ",", "24"}], "}"}]], "Output"]
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
         SubscriptBox["n", "sites"], "]"}]}], "]"}], "&"}], "/@", 
     SubscriptBox["site", 
      RowBox[{"order", ",", "HS"}]]}], "]"}], "\[IndentingNewLine]", 
   RowBox[{"Total", "[", 
    RowBox[{
     RowBox[{
      RowBox[{"Norm", "[", 
       RowBox[{
        RowBox[{"Sort", "[", "#", "]"}], "-", 
        RowBox[{"Range", "[", 
         SubscriptBox["n", "sites"], "]"}]}], "]"}], "&"}], "/@", 
     SubscriptBox["site", 
      RowBox[{"order", ",", "ph"}]]}], "]"}]}]}]], "Input"],

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
        SubscriptBox["\[Lambda]", "val"], ",", 
        SubscriptBox["\[CapitalOmega]", "val"], ",", 
        SubscriptBox["g", "val"], ",", 
        SubscriptBox["box", "width"], ",", 
        SubscriptBox["s", "0"], ",", 
        SubscriptBox["X", "0"], ",", 
        SubscriptBox["site", 
         RowBox[{"order", ",", "HS"}]], ",", 
        SubscriptBox["site", 
         RowBox[{"order", ",", "ph"}]], ",", 
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
        SubscriptBox["\[Lambda]", "val"], 
        SubscriptBox["s", "1"]}], ",", 
       RowBox[{
        SubscriptBox["\[Tau]", "val"], 
        SubscriptBox["g", "val"], 
        SubscriptBox["X", "1"]}]}], "]"}]}], "]"}], "]"}]}]], "Input"],

Cell[BoxData["0"], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   SubscriptBox["s", "1"], "\[LeftDoubleBracket]", 
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
    RowBox[{"1", ",", 
     RowBox[{"-", "1"}], ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "1"}], ",", "1", ",", 
     RowBox[{"-", "1"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1", ",", 
     RowBox[{"-", "1"}], ",", "1"}], "}"}]}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Flatten", "[", 
  RowBox[{
   RowBox[{
    RowBox[{"Transpose", "[", 
     SubscriptBox["s", "1"], "]"}], "/.", 
    RowBox[{"{", 
     RowBox[{"1", "\[Rule]", "0"}], "}"}]}], "/.", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "1"}], "\[Rule]", "1"}], "}"}]}], "]"}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
  "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", 
   ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", 
   ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", 
   ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", 
   "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", 
   ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", 
   "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", 
   ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", 
   ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", 
   ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", 
   ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "0", ",", 
   "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", 
   ",", "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", ",", "0", 
   ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", 
   "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", 
   ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", 
   "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", 
   ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", 
   "0", ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "0", 
   ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", 
   ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", 
   "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "1", 
   ",", "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", 
   "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", 
   ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", 
   "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", 
   ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", 
   "0", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", 
   ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", 
   "1", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", 
   ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", 
   "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "1", ",", "1", 
   ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", 
   "1", ",", "0", ",", "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", 
   ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", 
   "1", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", ",", "0", ",", "1", 
   ",", "1", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", "1", ",", 
   "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "0", ",", "1", 
   ",", "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "0", ",", "1", ",", 
   "0", ",", "0", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", "1", 
   ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", ",", "1", ",", 
   "0", ",", "1", ",", "0", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", 
   ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", "1", ",", 
   "0", ",", "1", ",", "1", ",", "1", ",", "1", ",", "0", ",", "0", ",", "1", 
   ",", "0"}], "}"}]], "Output"]
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
     {"0.785980262508502862277444160086626375534654186845494926142925181387583\
130028148378816463493582552`63.81835908076402", 
      RowBox[{
      "-", "0.0317580321568469849313885121870861307978723019584785914864916566\
07868633933670997939018444804179`63.81835908076401"}], 
      RowBox[{
      "-", "0.0039260514226327877186614596112500777022905597139785738010247819\
00607029189281304369274770027774`63.81835908076401"}]},
     {
      RowBox[{
      "-", "0.2178573613069781430213529291248626392724291896258907744962217725\
30423376146400170535302902038561`63.81835908076401"}], 
      "0.044814854512678494091541854911307867756520871370568277225801804137943\
451407022617631853821513569`63.81835908076401", 
      "0.035780449892028899176803442423886651121188395807268896067877426330301\
882396383117398224683036161`63.81835908076401"},
     {"0.246581507608330555022805614420229748196315101398447769212398737437142\
006914883622409679692118596`63.81835908076401", 
      "0.039346793219748989738213051714641791355172554924957467460361915230169\
956052027281481073861887591`63.81835908076401", 
      "0.763620600385500426125495150976470193641007031028795057078860950767335\
139534665522350768331937784`63.81835908076401"}
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
     {"0.427460099183264692780428358354238264994398708126290219378912285635272\
308353146095170365016493961`63.81835908076401", 
      RowBox[{
      "-", "0.4107727386530902364767169595502313678540279913728716099941835058\
00708479049452466703795083005187`63.81835908076401"}], 
      RowBox[{
      "-", "0.6247040888425440709850337941422302271970333431922055525968448592\
84711551915509242688498536148618`63.81835908076401"}]},
     {
      RowBox[{
      "-", "0.0241283414393269899683082665059891607670245347103316338536423049\
93555349249644438844142895719266`63.81835908076401"}], 
      "0.958969430909613625102108093758354237012119476688342278726440818701015\
346555831770635150280438313`63.81835908076401", 
      "0.023783765819311922626956530360775707318519705063019064702522158056575\
894501275693047532055279603`63.81835908076401"},
     {"0.020549103342282534852965196225022120369291636236043756421509583126925\
756154631006204271845342556`63.81835908076401", 
      "0.042556943589495381205763427679902928747745724100060891747919259750744\
698750215824056772778672393`63.81835908076401", 
      "0.332013713509182934406361914583342666475143269783730089550854612853645\
238693041639844595084443945`63.81835908076401"}
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

Cell[BoxData["1.5713283030369798`"], "Output"],

Cell[BoxData["0.988123029850554`"], "Output"]
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

Cell[BoxData["5.\
911885300376497560042768746055837626100722686350539287592540201336906907701036\
410887026607`63.81835908076401*^-45"], "Output"],

Cell[BoxData["1.\
410496866741598978881996121708435680965368795693677675302783447112844571414757\
33712965`63.81835908076401*^-29"], "Output"]
}, Open  ]],

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
WindowMargins->{{Automatic, 157}, {Automatic, 70}},
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
Cell[CellGroupData[{
Cell[1551, 39, 326, 10, 72, "Input"],
Cell[1880, 51, 31, 0, 31, "Output"]
}, Open  ]],
Cell[1926, 54, 219, 7, 52, "Input"],
Cell[2148, 63, 261, 8, 67, "Input"],
Cell[CellGroupData[{
Cell[2434, 75, 337, 11, 52, "Input"],
Cell[2774, 88, 46, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[2857, 93, 295, 9, 69, "Input"],
Cell[3155, 104, 29, 0, 31, "Output"]
}, Open  ]],
Cell[3199, 107, 248, 7, 52, "Input"],
Cell[3450, 116, 287, 9, 52, "Input"],
Cell[3740, 127, 315, 9, 52, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[4092, 141, 41, 0, 43, "Subsection"],
Cell[4136, 143, 320, 10, 72, "Input"],
Cell[CellGroupData[{
Cell[4481, 157, 947, 29, 72, "Input"],
Cell[5431, 188, 1469, 50, 52, "Output"],
Cell[6903, 240, 74, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[7014, 247, 295, 8, 52, "Input"],
Cell[7312, 257, 29, 0, 31, "Output"]
}, Open  ]],
Cell[7356, 260, 1122, 32, 52, "Input"],
Cell[CellGroupData[{
Cell[8503, 296, 336, 10, 72, "Input"],
Cell[8842, 308, 28, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[8907, 313, 358, 10, 52, "Input"],
Cell[9268, 325, 28, 0, 31, "Output"]
}, Open  ]],
Cell[9311, 328, 1167, 33, 59, "Input"],
Cell[CellGroupData[{
Cell[10503, 365, 340, 10, 72, "Input"],
Cell[10846, 377, 28, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[10911, 382, 426, 13, 52, "Input"],
Cell[11340, 397, 28, 0, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[11417, 403, 45, 0, 43, "Subsection"],
Cell[11465, 405, 221, 7, 52, "Input"],
Cell[11689, 414, 331, 11, 67, "Input"],
Cell[12023, 427, 258, 8, 67, "Input"],
Cell[12284, 437, 486, 16, 31, "Input"],
Cell[12773, 455, 513, 16, 52, "Input"],
Cell[CellGroupData[{
Cell[13311, 475, 480, 14, 52, "Input"],
Cell[13794, 491, 835, 20, 71, "Output"]
}, Open  ]],
Cell[14644, 514, 470, 14, 52, "Input"],
Cell[CellGroupData[{
Cell[15139, 532, 327, 9, 52, "Input"],
Cell[15469, 543, 47, 0, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[15565, 549, 48, 0, 43, "Subsection"],
Cell[CellGroupData[{
Cell[15638, 553, 556, 18, 72, "Input"],
Cell[16197, 573, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[16309, 580, 388, 11, 52, "Input"],
Cell[16700, 593, 407, 15, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[17144, 613, 305, 10, 31, "Input"],
Cell[17452, 625, 4100, 54, 152, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[21601, 685, 34, 0, 43, "Subsection"],
Cell[CellGroupData[{
Cell[21660, 689, 588, 18, 72, "Input"],
Cell[22251, 709, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[22363, 716, 429, 13, 52, "Input"],
Cell[22795, 731, 466, 14, 31, "Output"]
}, Open  ]],
Cell[23276, 748, 1068, 34, 69, "Input"],
Cell[CellGroupData[{
Cell[24369, 786, 1530, 47, 68, "Input"],
Cell[25902, 835, 28, 0, 31, "Output"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[25979, 841, 31, 0, 43, "Subsection"],
Cell[26013, 843, 625, 17, 52, "Input"],
Cell[26641, 862, 843, 22, 52, "Input"],
Cell[27487, 886, 612, 17, 72, "Input"],
Cell[28102, 905, 17840, 409, 888, "Input"],
Cell[45945, 1316, 242, 7, 31, "Input"],
Cell[CellGroupData[{
Cell[46212, 1327, 720, 22, 72, "Input"],
Cell[46935, 1351, 31, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[47003, 1356, 283, 7, 52, "Input"],
Cell[47289, 1365, 145, 4, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[47471, 1374, 795, 23, 72, "Input"],
Cell[48269, 1399, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[48381, 1406, 758, 22, 72, "Input"],
Cell[49142, 1430, 75, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[49254, 1437, 769, 23, 72, "Input"],
Cell[50026, 1462, 28, 0, 31, "Output"],
Cell[50057, 1464, 28, 0, 31, "Output"]
}, Open  ]],
Cell[50100, 1467, 1375, 39, 72, "Input"],
Cell[CellGroupData[{
Cell[51500, 1510, 891, 26, 52, "Input"],
Cell[52394, 1538, 28, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[52459, 1543, 388, 11, 52, "Input"],
Cell[52850, 1556, 344, 12, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[53231, 1573, 305, 10, 31, "Input"],
Cell[53539, 1585, 4100, 54, 152, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[57676, 1644, 877, 27, 72, "Input"],
Cell[58556, 1673, 1785, 37, 71, "Output"],
Cell[60344, 1712, 1785, 37, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[62166, 1754, 559, 17, 72, "Input"],
Cell[62728, 1773, 46, 0, 31, "Output"],
Cell[62777, 1775, 45, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[62859, 1780, 476, 15, 72, "Input"],
Cell[63338, 1797, 146, 2, 31, "Output"],
Cell[63487, 1801, 142, 2, 31, "Output"]
}, Open  ]],
Cell[63644, 1806, 749, 20, 52, "Input"],
Cell[64396, 1828, 1483, 43, 72, "Input"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature 2wTkUWCYLsersAwme#tmCVhG *)
