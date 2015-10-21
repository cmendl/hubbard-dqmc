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
NotebookDataLength[      6798,        222]
NotebookOptionsPosition[      7125,        210]
NotebookOutlinePosition[      7469,        225]
CellTagsIndexPosition[      7426,        222]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell[BoxData[{
 RowBox[{
  RowBox[{
   SubscriptBox["A", "mat"], "=", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"2", ",", 
       RowBox[{
        RowBox[{"-", "1"}], "/", "2"}], ",", "1", ",", 
       RowBox[{"5", "/", "2"}], ",", 
       RowBox[{
        RowBox[{"-", "3"}], "/", "2"}]}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{
        RowBox[{"-", "1"}], "/", "2"}], ",", "3", ",", 
       RowBox[{"1", "/", "2"}], ",", 
       RowBox[{"5", "/", "2"}], ",", 
       RowBox[{"3", "/", "2"}]}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"1", ",", 
       RowBox[{"1", "/", "2"}], ",", 
       RowBox[{"-", "3"}], ",", 
       RowBox[{"-", "1"}], ",", 
       RowBox[{"-", "1"}]}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"5", "/", "2"}], ",", 
       RowBox[{"5", "/", "2"}], ",", 
       RowBox[{"-", "1"}], ",", 
       RowBox[{
        RowBox[{"-", "5"}], "/", "2"}], ",", "2"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{
        RowBox[{"-", "3"}], "/", "2"}], ",", 
       RowBox[{"3", "/", "2"}], ",", 
       RowBox[{"-", "1"}], ",", "2", ",", "1"}], "}"}]}], "}"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"%", "//", "MatrixForm"}]}], "Input"],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"2", 
      RowBox[{"-", 
       FractionBox["1", "2"]}], "1", 
      FractionBox["5", "2"], 
      RowBox[{"-", 
       FractionBox["3", "2"]}]},
     {
      RowBox[{"-", 
       FractionBox["1", "2"]}], "3", 
      FractionBox["1", "2"], 
      FractionBox["5", "2"], 
      FractionBox["3", "2"]},
     {"1", 
      FractionBox["1", "2"], 
      RowBox[{"-", "3"}], 
      RowBox[{"-", "1"}], 
      RowBox[{"-", "1"}]},
     {
      FractionBox["5", "2"], 
      FractionBox["5", "2"], 
      RowBox[{"-", "1"}], 
      RowBox[{"-", 
       FractionBox["5", "2"]}], "2"},
     {
      RowBox[{"-", 
       FractionBox["3", "2"]}], 
      FractionBox["3", "2"], 
      RowBox[{"-", "1"}], "2", "1"}
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
  RowBox[{"(*", " ", "symmetric", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{"Norm", "[", 
   RowBox[{
    SubscriptBox["A", "mat"], "-", 
    RowBox[{"Transpose", "[", 
     SubscriptBox["A", "mat"], "]"}]}], "]"}]}]], "Input"],

Cell[BoxData["0"], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{"MatrixExp", "[", 
  RowBox[{"N", "[", 
   SubscriptBox["A", "mat"], "]"}], "]"}], "\[IndentingNewLine]", 
 RowBox[{"%", "//", "MatrixForm"}]}], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"25.955160669119373`", ",", 
     RowBox[{"-", "9.770400042675018`"}], ",", "4.871553760371625`", ",", 
     "2.7134175545724095`", ",", 
     RowBox[{"-", "15.334376171601765`"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "9.770400042675028`"}], ",", "105.8959791472407`", ",", 
     RowBox[{"-", "8.733225951221918`"}], ",", "50.05805070380273`", ",", 
     "66.66275673512335`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"4.871553760371625`", ",", 
     RowBox[{"-", "8.733225951221916`"}], ",", "1.6299718135475225`", ",", 
     RowBox[{"-", "3.2373011106885214`"}], ",", 
     RowBox[{"-", "7.6198789670106155`"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2.713417554572406`", ",", "50.05805070380271`", ",", 
     RowBox[{"-", "3.2373011106885223`"}], ",", "26.201262114829625`", ",", 
     "29.683741571682877`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"-", "15.334376171601768`"}], ",", "66.66275673512334`", ",", 
     RowBox[{"-", "7.6198789670106155`"}], ",", "29.683741571682877`", ",", 
     "47.30452009463436`"}], "}"}]}], "}"}]], "Output"],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"25.955160669119373`", 
      RowBox[{"-", "9.770400042675018`"}], "4.871553760371625`", 
      "2.7134175545724095`", 
      RowBox[{"-", "15.334376171601765`"}]},
     {
      RowBox[{"-", "9.770400042675028`"}], "105.8959791472407`", 
      RowBox[{"-", "8.733225951221918`"}], "50.05805070380273`", 
      "66.66275673512335`"},
     {"4.871553760371625`", 
      RowBox[{"-", "8.733225951221916`"}], "1.6299718135475225`", 
      RowBox[{"-", "3.2373011106885214`"}], 
      RowBox[{"-", "7.6198789670106155`"}]},
     {"2.713417554572406`", "50.05805070380271`", 
      RowBox[{"-", "3.2373011106885223`"}], "26.201262114829625`", 
      "29.683741571682877`"},
     {
      RowBox[{"-", "15.334376171601768`"}], "66.66275673512334`", 
      RowBox[{"-", "7.6198789670106155`"}], "29.683741571682877`", 
      "47.30452009463436`"}
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
}, Open  ]]
},
WindowSize->{1350, 770},
WindowMargins->{{163, Automatic}, {155, Automatic}},
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
Cell[1486, 35, 1250, 40, 52, "Input"],
Cell[2739, 77, 1261, 43, 151, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[4037, 125, 256, 7, 52, "Input"],
Cell[4296, 134, 28, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[4361, 139, 177, 4, 52, "Input"],
Cell[4541, 145, 1169, 26, 52, "Output"],
Cell[5713, 173, 1396, 34, 101, "Output"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature buDDe1S59n29hDKbN#90Gztu *)
