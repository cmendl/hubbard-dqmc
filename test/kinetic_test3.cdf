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
NotebookDataLength[     19364,        602]
NotebookOptionsPosition[     18814,        563]
NotebookOutlinePosition[     19158,        578]
CellTagsIndexPosition[     19115,        575]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{
Cell["Hexagonal lattice with nearest neighbor hopping", "Text"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"lattice", " ", "dimensions", " ", 
    RowBox[{"(", 
     RowBox[{
     "tilted", " ", "rectangular", " ", "lattice", " ", "forms", " ", 
      "triangular", " ", "lattice"}], ")"}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
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
    RowBox[{"hexagonal", " ", "lattice", " ", "sites"}], ",", " ", 
    RowBox[{
     RowBox[{"in", " ", "column"}], "-", 
     RowBox[{"major", " ", 
      RowBox[{"(", "Fortran", ")"}], " ", "order"}]}]}], " ", "*)"}], 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["latt", "sites"], "=", 
     RowBox[{"Flatten", "[", 
      RowBox[{
       RowBox[{"Transpose", "[", 
        RowBox[{
         RowBox[{"Table", "[", 
          RowBox[{
           RowBox[{
            RowBox[{"i", 
             RowBox[{"{", 
              RowBox[{"1", ",", "0"}], "}"}]}], "+", 
            RowBox[{"j", 
             RowBox[{"{", 
              RowBox[{
               FractionBox["1", "2"], ",", 
               FractionBox[
                SqrtBox["3"], "2"]}], "}"}]}], "+", 
            RowBox[{"\[Sigma]", 
             RowBox[{"{", 
              RowBox[{
               FractionBox["1", "2"], ",", 
               FractionBox["1", 
                RowBox[{"2", " ", 
                 SqrtBox["3"]}]]}], "}"}]}]}], ",", 
           RowBox[{"{", 
            RowBox[{"\[Sigma]", ",", "0", ",", "1"}], "}"}], ",", 
           RowBox[{"{", 
            RowBox[{"i", ",", "0", ",", 
             RowBox[{
              SubscriptBox["n", "x"], "-", "1"}]}], "}"}], ",", 
           RowBox[{"{", 
            RowBox[{"j", ",", "0", ",", 
             RowBox[{
              SubscriptBox["n", "y"], "-", "1"}]}], "}"}]}], "]"}], ",", 
         RowBox[{"{", 
          RowBox[{"1", ",", "3", ",", "2"}], "}"}]}], "]"}], ",", "2"}], 
      "]"}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{"Dimensions", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"60", ",", "2"}], "}"}]], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"HexMod", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{"x_", ",", "y_"}], "}"}], ",", 
    RowBox[{"{", 
     RowBox[{"nx_", ",", "ny_"}], "}"}]}], "]"}], ":=", 
  RowBox[{"{", 
   RowBox[{
    RowBox[{"Mod", "[", 
     RowBox[{"x", ",", "nx", ",", 
      RowBox[{"-", "1"}]}], "]"}], ",", 
    RowBox[{"Mod", "[", 
     RowBox[{"y", ",", 
      RowBox[{
       FractionBox[
        SqrtBox["3"], "2"], "ny"}], ",", 
      RowBox[{"-", 
       FractionBox[
        SqrtBox["3"], "2"]}]}], "]"}]}], "}"}]}]], "Input"],

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
           RowBox[{"FullSimplify", "[", 
            RowBox[{"HexMod", "[", 
             RowBox[{
              RowBox[{"#1", "-", "#2"}], ",", 
              RowBox[{"{", 
               RowBox[{
                SubscriptBox["n", "x"], ",", 
                SubscriptBox["n", "y"]}], "}"}]}], "]"}], "]"}], "]"}], 
          "\[Equal]", 
          FractionBox["1", 
           SqrtBox["3"]]}], ",", "1", ",", "0"}], "]"}], "&"}], ",", 
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
   RowBox[{
   "every", " ", "site", " ", "should", " ", "have", " ", "3", " ", "nearest",
     " ", "neighbors"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{"Norm", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{"Total", "/@", 
      SubscriptBox["neigh", "nearest"]}], ")"}], "-", "3"}], "]"}]}]], "Input"],

Cell[BoxData["0"], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
   "visualize", " ", "lattice", " ", "and", " ", "some", " ", "nearest", " ", 
    "neighbors"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{"Show", "[", 
   RowBox[{
    RowBox[{"ListPlot", "[", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{
        SubscriptBox["latt", "sites"], "\[LeftDoubleBracket]", 
        RowBox[{"1", ";;", 
         RowBox[{
          RowBox[{"Length", "[", 
           SubscriptBox["latt", "sites"], "]"}], "/", "2"}]}], 
        "\[RightDoubleBracket]"}], ",", 
       RowBox[{
        SubscriptBox["latt", "sites"], "\[LeftDoubleBracket]", 
        RowBox[{
         RowBox[{
          RowBox[{
           RowBox[{"Length", "[", 
            SubscriptBox["latt", "sites"], "]"}], "/", "2"}], "+", "1"}], ";;", 
         RowBox[{"-", "1"}]}], "\[RightDoubleBracket]"}]}], "}"}], "]"}], ",", 
    RowBox[{"Graphics", "[", 
     RowBox[{"{", 
      RowBox[{"Dashed", ",", 
       RowBox[{"Table", "[", 
        RowBox[{
         RowBox[{"If", "[", 
          RowBox[{
           RowBox[{
            RowBox[{
             SubscriptBox["neigh", "nearest"], "\[LeftDoubleBracket]", 
             RowBox[{"i", ",", "j"}], "\[RightDoubleBracket]"}], "\[Equal]", 
            "1"}], ",", 
           RowBox[{"Line", "[", 
            RowBox[{"{", 
             RowBox[{
              RowBox[{
               SubscriptBox["latt", "sites"], "\[LeftDoubleBracket]", "i", 
               "\[RightDoubleBracket]"}], ",", 
              RowBox[{
               SubscriptBox["latt", "sites"], "\[LeftDoubleBracket]", "j", 
               "\[RightDoubleBracket]"}]}], "}"}], "]"}], ",", 
           RowBox[{"{", "}"}]}], "]"}], ",", 
         RowBox[{"{", 
          RowBox[{"j", ",", 
           RowBox[{"Length", "[", 
            SubscriptBox["latt", "sites"], "]"}]}], "}"}], ",", 
         RowBox[{"{", 
          RowBox[{"i", ",", 
           RowBox[{"{", 
            RowBox[{"2", ",", "17", ",", "55"}], "}"}]}], "}"}]}], "]"}]}], 
      "}"}], "]"}]}], "]"}]}]], "Input"],

Cell[BoxData[
 GraphicsBox[{{{}, {
     {RGBColor[0.368417, 0.506779, 0.709798], PointSize[
      0.016666666666666666`], AbsoluteThickness[1.6], PointBox[CompressedData["
1:eJxTTMoPSmViYGCQA2IQjR18sEcTcEDlcqDxBdD4D+xX+US8qNr2GmrODzQ+
iwMqnweNL4TG/wDV/xvmLgdUPgcaXwCNL4LG/2HfYOmkd+cEiwPMPah8HjS+
EBpfDI0Pcw+3A6p7YHwBNL4IGl8Cjc/i8Frf/CJPkCDcPah8ITS+GBpfCs4H
AFyybLo=
       "]]}, 
     {RGBColor[0.880722, 0.611041, 0.142051], PointSize[
      0.016666666666666666`], AbsoluteThickness[1.6], PointBox[CompressedData["
1:eJxTTMoPSmViYGCQA2IQDQEP7GWMJ7gur7xkD+H/QOOzOKDyedD4Qmj8D1D9
n6B8BgdUPgcaXwCNL4LG/2G/U+deTJEWgwPMPah8HjS+EBpfDI3P4PDY/ouY
oAS7A8w9qHwBNL4IGl8Cjc/iwBfMdXE9Kx/cPah8ITS+GBpfCo0PCx8huHtQ
+SJofAk0vgycDwAnSmSN
       "]]}}, {}}, 
   {Dashing[{
     Small, Small}], {{}, {}, {}}, {{}, {}, {}}, {{}, {}, {}}, {{}, {}, {}}, \
{{}, {}, {}}, {{}, {}, {}}, {{}, {}, {}}, {{}, {}, {}}, {{}, {}, {}}, {{}, \
{}, {}}, {{}, {}, {}}, {{}, {}, {}}, {{}, {}, {}}, {{}, {}, {}}, {{}, {}, \
{}}, {{}, {}, {}}, {{}, {}, {}}, {{}, {}, {}}, {{}, {}, {}}, {{}, {}, {}}, \
{{}, {}, LineBox[
      NCache[{{Rational[13, 2], Rational[1, 2] 3^Rational[-1, 2] + 
         2 3^Rational[1, 2]}, {2, 2 3^Rational[1, 2]}}, {{6.5, 
        3.7527767497325675`}, {
        2, 3.4641016151377544`}}]]}, {{}, {}, {}}, {{}, {}, {}}, {{}, {}, \
{}}, {{}, {}, 
     LineBox[NCache[{{
        Rational[13, 2], Rational[1, 2] 3^Rational[-1, 2] + 
         2 3^Rational[1, 2]}, {6, 2 3^Rational[1, 2]}}, {{6.5, 
        3.7527767497325675`}, {
        6, 3.4641016151377544`}}]]}, {{}, {}, {}}, {{}, {}, {}}, {{}, {}, \
{}}, {{}, {}, {}}, {{}, {}, 
     LineBox[NCache[{{
        Rational[13, 2], Rational[1, 2] 3^Rational[-1, 2] + 
         2 3^Rational[1, 2]}, {
        Rational[13, 2], Rational[5, 2] 3^Rational[1, 2]}}, {{6.5, 
       3.7527767497325675`}, {6.5, 4.330127018922193}}]]}, {
     LineBox[NCache[{{1, 0}, {
        Rational[1, 2], Rational[1, 2] 3^Rational[-1, 2]}}, {{1, 0}, {0.5, 
        0.2886751345948129}}]], {}, {}}, {
     LineBox[
      NCache[{{1, 0}, {Rational[3, 2], Rational[1, 2] 3^Rational[-1, 2]}}, {{
        1, 0}, {1.5, 
        0.2886751345948129}}]], {}, {}}, {{}, {}, {}}, {{}, {}, {}}, {{}, {}, \
{}}, {{}, {}, {}}, {{}, {}, {}}, {{}, {}, {}}, {{}, {}, {}}, {{}, {}, {}}, \
{{}, {}, {}}, {{}, 
     LineBox[NCache[{{Rational[5, 2], Rational[3, 2] 3^Rational[1, 2]}, {
        Rational[5, 2], Rational[1, 2] 3^Rational[-1, 2] + 
         3^Rational[1, 2]}}, {{2.5, 2.598076211353316}, {2.5, 
       2.0207259421636903`}}]], {}}, {{}, {}, {}}, {{}, {}, {}}, {{}, {}, \
{}}, {{}, 
     LineBox[NCache[{{Rational[5, 2], Rational[3, 2] 3^Rational[1, 2]}, {
        2, Rational[1, 2] 3^Rational[-1, 2] + 
         Rational[3, 2] 3^Rational[1, 2]}}, {{2.5, 2.598076211353316}, {
        2, 2.886751345948129}}]], {}}, {{}, 
     LineBox[NCache[{{Rational[5, 2], Rational[3, 2] 3^Rational[1, 2]}, {
        3, Rational[1, 2] 3^Rational[-1, 2] + 
         Rational[3, 2] 3^Rational[1, 2]}}, {{2.5, 2.598076211353316}, {
        3, 2.886751345948129}}]], {}}, {{}, {}, {}}, {{}, {}, {}}, {{}, {}, \
{}}, {{}, {}, {}}, {{}, {}, {}}, {{}, {}, {}}, {{}, {}, {}}, {{}, {}, {}}, \
{{}, {}, {}}, {{}, {}, {}}, {{}, {}, {}}, {
     LineBox[NCache[{{1, 0}, {
        6, Rational[1, 2] 3^Rational[-1, 2] + 
         Rational[5, 2] 3^Rational[1, 2]}}, {{1, 0}, {
        6, 4.618802153517006}}]], {}, {}}, {{}, {}, {}}}},
  AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
  Axes->{True, True},
  AxesLabel->{None, None},
  AxesOrigin->{0, 0},
  DisplayFunction->Identity,
  Frame->{{False, False}, {False, False}},
  FrameLabel->{{None, None}, {None, None}},
  FrameTicks->{{Automatic, Automatic}, {Automatic, Automatic}},
  GridLines->{None, None},
  GridLinesStyle->Directive[
    GrayLevel[0.5, 0.4]],
  Method->{},
  PlotRange->{{0, 7.}, {0, 4.618802153517006}},
  PlotRangeClipping->True,
  PlotRangePadding->{{
     Scaled[0.02], 
     Scaled[0.02]}, {
     Scaled[0.02], 
     Scaled[0.05]}},
  Ticks->{Automatic, Automatic}]], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"t", " ", "hopping", " ", "parameter"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["t", "val"], "=", 
    FractionBox["4", "5"]}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"chemical", " ", "potential"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[Mu]", "val"], "=", 
    FractionBox["2", "9"]}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"time", " ", "step"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[CapitalDelta]t", "val"], "=", 
    FractionBox["1", "7"]}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"factor", " ", 
    RowBox[{"(", 
     RowBox[{"-", "1"}], ")"}], " ", "from", " ", "negative", " ", "sign", 
    " ", "in", " ", "Hamiltonian"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["K", "val"], "=", 
    RowBox[{"-", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        SubscriptBox["t", "val"], 
        SubscriptBox["neigh", "nearest"]}], "+", 
       RowBox[{
        SubscriptBox["\[Mu]", "val"], 
        RowBox[{"IdentityMatrix", "[", 
         RowBox[{"2", 
          SubscriptBox["n", "x"], 
          SubscriptBox["n", "y"]}], "]"}]}]}], ")"}]}]}], ";"}]}]], "Input"],

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
     {"1.0525894484146352`", "0.006800150592028771`", 
      "5.055417907590732`*^-7"},
     {"0.006800150592028771`", "1.0525894484146352`", 
      "0.0005166460945511059`"},
     {"5.055417907590743`*^-7", "0.0005166460945511058`", 
      "1.0525894484146352`"}
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
     {"0.9878357955216168`", "0.006381815987098333`", 
      RowBox[{"-", "4.744416522472982`*^-7"}]},
     {"0.006381815987098333`", "0.9878357955216168`", 
      RowBox[{"-", "0.0004848628366764579`"}]},
     {
      RowBox[{"-", "4.7444165224729933`*^-7"}], 
      RowBox[{"-", "0.0004848628366764578`"}], "0.9878357955216168`"}
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

Cell[BoxData["2.0061041888474687`*^-15"], "Output"]
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
   RowBox[{"save", " ", "as", " ", "reference", " ", "to", " ", "disk"}], " ",
    "*)"}], "\[IndentingNewLine]", 
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
WindowSize->{1350, 770},
WindowMargins->{{Automatic, 323}, {Automatic, 177}},
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
Cell[1464, 33, 63, 0, 30, "Text"],
Cell[1530, 35, 477, 15, 72, "Input"],
Cell[CellGroupData[{
Cell[2032, 54, 1686, 49, 102, "Input"],
Cell[3721, 105, 74, 2, 31, "Output"]
}, Open  ]],
Cell[3810, 110, 554, 20, 51, "Input"],
Cell[4367, 132, 942, 28, 73, "Input"],
Cell[CellGroupData[{
Cell[5334, 164, 336, 10, 72, "Input"],
Cell[5673, 176, 28, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[5738, 181, 361, 10, 52, "Input"],
Cell[6102, 193, 28, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[6167, 198, 2061, 56, 72, "Input"],
Cell[8231, 256, 4051, 85, 249, "Output"]
}, Open  ]],
Cell[12297, 344, 244, 8, 67, "Input"],
Cell[12544, 354, 236, 7, 67, "Input"],
Cell[12783, 363, 238, 7, 67, "Input"],
Cell[13024, 372, 682, 21, 52, "Input"],
Cell[13709, 395, 461, 15, 52, "Input"],
Cell[14173, 412, 419, 13, 52, "Input"],
Cell[CellGroupData[{
Cell[14617, 429, 436, 13, 52, "Input"],
Cell[15056, 444, 816, 21, 77, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[15909, 470, 439, 13, 52, "Input"],
Cell[16351, 485, 884, 22, 77, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[17272, 512, 415, 12, 52, "Input"],
Cell[17690, 526, 51, 0, 31, "Output"]
}, Open  ]],
Cell[17756, 529, 242, 7, 31, "Input"],
Cell[18001, 538, 809, 23, 72, "Input"]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature kvpWe#D4WKV0RDgl82Eg6ypk *)
