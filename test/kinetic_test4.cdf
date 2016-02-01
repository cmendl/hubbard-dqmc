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
NotebookDataLength[     22989,        670]
NotebookOptionsPosition[     22439,        631]
NotebookOutlinePosition[     22782,        646]
CellTagsIndexPosition[     22739,        643]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{
Cell["Kagome lattice with nearest neighbor hopping", "Text"],

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
    RowBox[{"Kagome", " ", "lattice", " ", "sites"}], ",", " ", 
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
            RowBox[{
             RowBox[{"{", 
              RowBox[{
               RowBox[{"{", 
                RowBox[{"0", ",", "0"}], "}"}], ",", 
               RowBox[{"{", 
                RowBox[{
                 FractionBox["1", "4"], ",", 
                 FractionBox[
                  SqrtBox["3"], "4"]}], "}"}], ",", 
               RowBox[{"{", 
                RowBox[{
                 RowBox[{"-", 
                  FractionBox["1", "4"]}], ",", 
                 FractionBox[
                  SqrtBox["3"], "4"]}], "}"}]}], "}"}], 
             "\[LeftDoubleBracket]", "\[Sigma]", "\[RightDoubleBracket]"}]}], 
           ",", 
           RowBox[{"{", 
            RowBox[{"\[Sigma]", ",", "3"}], "}"}], ",", 
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
  RowBox[{"90", ",", "2"}], "}"}]], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"TriMod", "[", 
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
            RowBox[{"TriMod", "[", 
             RowBox[{
              RowBox[{"#1", "-", "#2"}], ",", 
              RowBox[{"{", 
               RowBox[{
                SubscriptBox["n", "x"], ",", 
                SubscriptBox["n", "y"]}], "}"}]}], "]"}], "]"}], "]"}], 
          "\[Equal]", 
          FractionBox["1", "2"]}], ",", "1", ",", "0"}], "]"}], "&"}], ",", 
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
   "every", " ", "site", " ", "should", " ", "have", " ", "4", " ", "nearest",
     " ", "neighbors"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{"Norm", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{"Total", "/@", 
      SubscriptBox["neigh", "nearest"]}], ")"}], "-", "4"}], "]"}]}]], "Input"],

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
     RowBox[{
      RowBox[{"{", 
       RowBox[{
        RowBox[{
         SubscriptBox["latt", "sites"], "\[LeftDoubleBracket]", 
         RowBox[{"1", ";;", 
          RowBox[{
           FractionBox["1", "3"], 
           RowBox[{"Length", "[", 
            SubscriptBox["latt", "sites"], "]"}]}]}], 
         "\[RightDoubleBracket]"}], ",", 
        RowBox[{
         SubscriptBox["latt", "sites"], "\[LeftDoubleBracket]", 
         RowBox[{
          RowBox[{
           RowBox[{
            FractionBox["1", "3"], 
            RowBox[{"Length", "[", 
             SubscriptBox["latt", "sites"], "]"}]}], "+", "1"}], ";;", 
          RowBox[{
           FractionBox["2", "3"], 
           RowBox[{"Length", "[", 
            SubscriptBox["latt", "sites"], "]"}]}]}], 
         "\[RightDoubleBracket]"}], ",", 
        RowBox[{
         SubscriptBox["latt", "sites"], "\[LeftDoubleBracket]", 
         RowBox[{
          RowBox[{
           RowBox[{
            FractionBox["2", "3"], 
            RowBox[{"Length", "[", 
             SubscriptBox["latt", "sites"], "]"}]}], "+", "1"}], ";;", 
          RowBox[{"-", "1"}]}], "\[RightDoubleBracket]"}]}], "}"}], ",", 
      RowBox[{"PlotRange", "\[Rule]", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"{", 
          RowBox[{
           RowBox[{
            RowBox[{"-", "1"}], "/", "4"}], ",", "7"}], "}"}], ",", 
         RowBox[{"{", 
          RowBox[{
           RowBox[{
            RowBox[{"-", "1"}], "/", "4"}], ",", "7"}], "}"}]}], "}"}]}], ",", 
      RowBox[{"AspectRatio", "\[Rule]", "1"}]}], "]"}], ",", 
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
            RowBox[{"2", ",", "17", ",", "46", ",", "72"}], "}"}]}], "}"}]}], 
        "]"}]}], "}"}], "]"}]}], "]"}]}]], "Input"],

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
1:eJxTTMoPSmViYGCQA2IQDQEX7Ff5RLyo2nbbHsL/gsZnckDlc6HxBdH4L+wb
LJ307pz4AuX/QeOzOaDy+dD4wmj8L/av9c0v8gQxOsDcg8rnQuMLovFF0fh/
7EWdfbtX23M4wNyDyudD4wuj8cXR+EwODmHJrl66/HD3oPIF0fiiaHxJND6b
w9Yte+pLeIXh7kHlC6PxxdH40nA+AP3qb9E=
       "]]}, 
     {RGBColor[0.560181, 0.691569, 0.194885], PointSize[
      0.016666666666666666`], AbsoluteThickness[1.6], PointBox[CompressedData["
1:eJxTTMoPSmViYGCQA2IQDQEX9q/yiXhRte22PYT/wh6V/weNz+aAyudD41+w
b7B00rtz4guU/wWNz+SAyudC4wui8V/Yv9Y3v8gTxOgAcw8qn80Blc+HxhdG
43+xF3X27V5tz+EAcw8qnwuNL4jGF0Xj/7F3CEt29dLlh7sHlc+HxhdG44uj
8Zkctm7ZU1/CKwx3DypfEI0visaXhPMBwEFyvQ==
       "]]}}, {}}, 
   {Dashing[{
     Small, Small}], {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, \
{{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, \
{}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, \
{}, LineBox[
      NCache[{{Rational[7, 4], Rational[5, 4] 3^Rational[1, 2]}, {
        2, 3^Rational[1, 2]}}, {{1.75, 2.1650635094610964`}, {
        2, 1.7320508075688772`}}]]}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, \
{}, {}, {}}, {{}, {}, 
     LineBox[NCache[{{Rational[7, 4], Rational[7, 4] 3^Rational[1, 2]}, {
        Rational[3, 2], Rational[3, 2] 3^Rational[1, 2]}}, {{1.75, 
       3.031088913245535}, {1.5, 2.598076211353316}}]], 
     LineBox[NCache[{{Rational[7, 4], Rational[5, 4] 3^Rational[1, 2]}, {
        Rational[3, 2], Rational[3, 2] 3^Rational[1, 2]}}, {{1.75, 
       2.1650635094610964`}, {1.5, 
       2.598076211353316}}]]}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, \
{}, {}}, {{}, {}, {}, {}}, {{}, {}, 
     LineBox[NCache[{{Rational[7, 4], Rational[7, 4] 3^Rational[1, 2]}, {
        2, 2 3^Rational[1, 2]}}, {{1.75, 3.031088913245535}, {
        2, 3.4641016151377544`}}]], {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, \
{{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, \
{}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {
     LineBox[NCache[{{1, 0}, {
        Rational[5, 4], Rational[1, 4] 3^Rational[1, 2]}}, {{1, 0}, {1.25, 
        0.4330127018922193}}]], {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, \
{}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, \
{{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, 
     LineBox[NCache[{{Rational[7, 4], Rational[5, 4] 3^Rational[1, 2]}, {
        Rational[5, 4], Rational[5, 4] 3^Rational[1, 2]}}, {{1.75, 
       2.1650635094610964`}, {1.25, 2.1650635094610964`}}]]}, {{}, 
     LineBox[NCache[{{Rational[5, 2], Rational[3, 2] 3^Rational[1, 2]}, {
        Rational[9, 4], Rational[5, 4] 3^Rational[1, 2]}}, {{2.5, 
       2.598076211353316}, {2.25, 2.1650635094610964`}}]], {}, 
     LineBox[NCache[{{Rational[7, 4], Rational[5, 4] 3^Rational[1, 2]}, {
        Rational[9, 4], Rational[5, 4] 3^Rational[1, 2]}}, {{1.75, 
       2.1650635094610964`}, {2.25, 
       2.1650635094610964`}}]]}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, \
{}, {}}, {{}, {}, {}, {}}, {{}, 
     LineBox[NCache[{{Rational[5, 2], Rational[3, 2] 3^Rational[1, 2]}, {
        Rational[11, 4], Rational[7, 4] 3^Rational[1, 2]}}, {{2.5, 
       2.598076211353316}, {2.75, 
       3.031088913245535}}]], {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, \
{{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, \
{}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, \
{}, {}}, {
     LineBox[NCache[{{1, 0}, {
        Rational[23, 4], Rational[11, 4] 3^Rational[1, 2]}}, {{1, 0}, {5.75, 
        4.763139720814412}}]], {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, \
{}}, {LineBox[
      NCache[{{1, 0}, {Rational[3, 4], Rational[1, 4] 3^Rational[1, 2]}}, {{1,
         0}, {0.75, 
        0.4330127018922193}}]], {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, \
{}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, \
{{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, 
     LineBox[NCache[{{Rational[5, 2], Rational[3, 2] 3^Rational[1, 2]}, {
        Rational[11, 4], Rational[5, 4] 3^Rational[1, 2]}}, {{2.5, 
       2.598076211353316}, {2.75, 
       2.1650635094610964`}}]], {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, \
{{}, {}, LineBox[
      NCache[{{Rational[7, 4], Rational[7, 4] 3^Rational[1, 2]}, {
        Rational[5, 4], Rational[7, 4] 3^Rational[1, 2]}}, {{1.75, 
       3.031088913245535}, {1.25, 3.031088913245535}}]], {}}, {{}, 
     LineBox[NCache[{{Rational[5, 2], Rational[3, 2] 3^Rational[1, 2]}, {
        Rational[9, 4], Rational[7, 4] 3^Rational[1, 2]}}, {{2.5, 
       2.598076211353316}, {2.25, 3.031088913245535}}]], 
     LineBox[NCache[{{Rational[7, 4], Rational[7, 4] 3^Rational[1, 2]}, {
        Rational[9, 4], Rational[7, 4] 3^Rational[1, 2]}}, {{1.75, 
       3.031088913245535}, {2.25, 
       3.031088913245535}}]], {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, \
{}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, \
{}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, {}}, {{}, {}, {}, \
{}}, {{}, {}, {}, {}}, {
     LineBox[NCache[{{1, 0}, {
        Rational[25, 4], Rational[11, 4] 3^Rational[1, 2]}}, {{1, 0}, {6.25, 
        4.763139720814412}}]], {}, {}, {}}}},
  AspectRatio->1,
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
  PlotRange->NCache[{{
      Rational[-1, 4], 7}, {
      Rational[-1, 4], 7}}, {{-0.25, 7}, {-0.25, 7}}],
  PlotRangeClipping->True,
  PlotRangePadding->{{0, 0}, {0, 0}},
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
    FractionBox["6", "5"]}], ";"}]}]], "Input"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"chemical", " ", "potential"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["\[Mu]", "val"], "=", 
    FractionBox["2", "7"]}], ";"}]}]], "Input"],

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
         RowBox[{"3", 
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
     {"1.0709935634758616`", "0.0009097475306832663`", 
      "0.0004701180512032681`"},
     {"0.0009097475306832649`", "1.0709935634758605`", "0.14989380414453857`"},
     {"0.0004701180512032676`", "0.14989380414453854`", "1.0709935634758603`"}
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
     {"1.002004932045422`", 
      RowBox[{"-", "0.0006989928628410956`"}], 
      RowBox[{"-", "0.00033793868524501045`"}]},
     {
      RowBox[{"-", "0.0006989928628410949`"}], "1.002004932045422`", 
      RowBox[{"-", "0.12306461132589853`"}]},
     {
      RowBox[{"-", "0.00033793868524501007`"}], 
      RowBox[{"-", "0.12306461132589853`"}], "1.0020049320454234`"}
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

Cell[BoxData["2.5705268929063077`*^-15"], "Output"]
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
WindowSize->{1406, 900},
WindowMargins->{{Automatic, 272}, {67, Automatic}},
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
Cell[1464, 33, 60, 0, 30, "Text"],
Cell[1527, 35, 477, 15, 72, "Input"],
Cell[CellGroupData[{
Cell[2029, 54, 2079, 60, 95, "Input"],
Cell[4111, 116, 74, 2, 31, "Output"]
}, Open  ]],
Cell[4200, 121, 554, 20, 51, "Input"],
Cell[4757, 143, 921, 27, 67, "Input"],
Cell[CellGroupData[{
Cell[5703, 174, 336, 10, 72, "Input"],
Cell[6042, 186, 28, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[6107, 191, 361, 10, 52, "Input"],
Cell[6471, 203, 28, 0, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[6536, 208, 3030, 84, 123, "Input"],
Cell[9569, 294, 6313, 115, 382, "Output"]
}, Open  ]],
Cell[15897, 412, 244, 8, 67, "Input"],
Cell[16144, 422, 236, 7, 67, "Input"],
Cell[16383, 431, 238, 7, 67, "Input"],
Cell[16624, 440, 682, 21, 52, "Input"],
Cell[17309, 463, 461, 15, 52, "Input"],
Cell[17773, 480, 419, 13, 52, "Input"],
Cell[CellGroupData[{
Cell[18217, 497, 436, 13, 52, "Input"],
Cell[18656, 512, 800, 19, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[19493, 536, 439, 13, 52, "Input"],
Cell[19935, 551, 925, 24, 71, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[20897, 580, 415, 12, 52, "Input"],
Cell[21315, 594, 51, 0, 31, "Output"]
}, Open  ]],
Cell[21381, 597, 242, 7, 31, "Input"],
Cell[21626, 606, 809, 23, 72, "Input"]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature Jx0nm@#GTW0StAgNNICWUM4a *)
