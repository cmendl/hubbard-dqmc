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
NotebookDataLength[      7429,        255]
NotebookOptionsPosition[      7679,        242]
NotebookOutlinePosition[      8022,        257]
CellTagsIndexPosition[      7979,        254]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{
Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"Eq", " ", 
    RowBox[{"(", "15", ")"}]}], " ", "*)"}], "\n", 
  RowBox[{
   RowBox[{"GreenSpinFlipUpdate", "[", 
    RowBox[{"i_", ",", "\[Delta]_", ",", "G_"}], "]"}], ":=", 
   RowBox[{"Module", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"c", "=", 
       RowBox[{"\[Delta]", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"UnitVector", "[", 
           RowBox[{
            RowBox[{"Length", "[", "G", "]"}], ",", "i"}], "]"}], "-", 
          RowBox[{
          "G", "\[LeftDoubleBracket]", "i", "\[RightDoubleBracket]"}]}], 
         ")"}]}]}], "}"}], ",", 
     RowBox[{"G", "-", 
      RowBox[{"KroneckerProduct", "[", 
       RowBox[{
        RowBox[{
         RowBox[{"G", "\[LeftDoubleBracket]", 
          RowBox[{";;", ",", "i"}], "\[RightDoubleBracket]"}], "/", 
         RowBox[{"(", 
          RowBox[{"1", "+", 
           RowBox[{
           "c", "\[LeftDoubleBracket]", "i", "\[RightDoubleBracket]"}]}], 
          ")"}]}], ",", "c"}], "]"}]}]}], "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["\[Delta]", "val"], "=", 
   RowBox[{"1", "/", "7"}]}], ";"}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["n", "sites"], "=", "24"}], ";"}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"initial", " ", 
    RowBox[{"Green", "'"}], "s", " ", "function", " ", "matrix"}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{"SeedRandom", "[", "42", "]"}], ";"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{
     SubscriptBox["G", "0"], "=", 
     RowBox[{"Table", "[", 
      RowBox[{
       FractionBox[
        RowBox[{"RandomInteger", "[", 
         RowBox[{"{", 
          RowBox[{
           RowBox[{"-", "12"}], ",", "12"}], "}"}], "]"}], 
        RowBox[{"RandomInteger", "[", 
         RowBox[{"{", 
          RowBox[{"1", ",", "12"}], "}"}], "]"}]], ",", 
       RowBox[{"{", 
        RowBox[{"i", ",", 
         SubscriptBox["n", "sites"]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{"j", ",", 
         SubscriptBox["n", "sites"]}], "}"}]}], "]"}]}], ";"}]}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["G", "0"], "\[LeftDoubleBracket]", 
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
      FractionBox["1", "11"], 
      RowBox[{"-", "10"}], 
      RowBox[{"-", "1"}]},
     {
      RowBox[{"-", 
       FractionBox["1", "3"]}], 
      RowBox[{"-", 
       FractionBox["12", "7"]}], "2"},
     {
      FractionBox["4", "3"], 
      FractionBox["2", "3"], 
      FractionBox["4", "3"]}
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
   RowBox[{"apply", " ", "Eq", " ", 
    RowBox[{"(", "15", ")"}]}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["G", "1"], "=", 
    RowBox[{"GreenSpinFlipUpdate", "[", 
     RowBox[{"12", ",", 
      SubscriptBox["\[Delta]", "val"], ",", 
      SubscriptBox["G", "0"]}], "]"}]}], ";"}]}]], "Input"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["G", "1"], "\[LeftDoubleBracket]", 
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
      FractionBox["38", "55"], 
      RowBox[{"-", 
       FractionBox["552", "55"]}], 
      RowBox[{"-", 
       FractionBox["7", "5"]}]},
     {
      RowBox[{"-", 
       FractionBox["31", "120"]}], 
      RowBox[{"-", 
       FractionBox["2647", "1540"]}], 
      FractionBox["39", "20"]},
     {
      FractionBox["25", "12"], 
      FractionBox["41", "66"], 
      FractionBox["5", "6"]}
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
   RowBox[{"save", " ", "initial", " ", "and", " ", "updated", " ", 
    RowBox[{"Green", "'"}], "s", " ", "function", " ", "to", " ", "disk"}], 
   " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"Table", "[", 
    RowBox[{
     RowBox[{"Export", "[", 
      RowBox[{
       RowBox[{
        RowBox[{"NotebookDirectory", "[", "]"}], "<>", 
        RowBox[{"FileBaseName", "[", 
         RowBox[{"NotebookFileName", "[", "]"}], "]"}], "<>", "\"\<_G\>\"", "<>", 
        RowBox[{"ToString", "[", "i", "]"}], "<>", "\"\<.dat\>\""}], ",", 
       RowBox[{"N", "[", 
        RowBox[{"Flatten", "[", 
         RowBox[{"Transpose", "[", 
          SubscriptBox["G", "i"], "]"}], "]"}], "]"}], ",", 
       "\"\<Real64\>\""}], "]"}], ",", 
     RowBox[{"{", 
      RowBox[{"i", ",", "0", ",", "1"}], "}"}]}], "]"}], ";"}]}]], "Input"]
},
WindowSize->{1351, 778},
WindowMargins->{{Automatic, 385}, {87, Automatic}},
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
Cell[1464, 33, 1065, 31, 52, "Input"],
Cell[2532, 66, 123, 4, 31, "Input"],
Cell[2658, 72, 95, 3, 31, "Input"],
Cell[2756, 77, 888, 27, 89, "Input"],
Cell[CellGroupData[{
Cell[3669, 108, 434, 13, 52, "Input"],
Cell[4106, 123, 860, 28, 101, "Output"]
}, Open  ]],
Cell[4981, 154, 379, 11, 52, "Input"],
Cell[CellGroupData[{
Cell[5385, 169, 434, 13, 52, "Input"],
Cell[5822, 184, 954, 31, 101, "Output"]
}, Open  ]],
Cell[6791, 218, 884, 22, 52, "Input"]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature #v0zqCFiDph@HDwk1Yvp9gkh *)
