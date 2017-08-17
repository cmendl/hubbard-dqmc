(* Content-type: application/vnd.wolfram.cdf.text *)

(*** Wolfram CDF File ***)
(* http://www.wolfram.com/cdf *)

(* CreatedBy='Mathematica 11.1' *)

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
NotebookDataLength[     74744,       2335]
NotebookOptionsPosition[     66453,       2118]
NotebookOutlinePosition[     66798,       2133]
CellTagsIndexPosition[     66755,       2130]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["\<\
Correlation functions in terms of the single-particle \
Green\[CloseCurlyQuote]s function\
\>", "Subsection",ExpressionUUID->"04451415-e7ea-42cb-b75b-d7fcf4206ed1"],

Cell[CellGroupData[{

Cell["General definitions", "Subsubsection",ExpressionUUID->"5dcc4ab5-6a1d-4e92-a2a7-0f59b76cf7e5"],

Cell[BoxData[
 RowBox[{
  RowBox[{"ExpandNCM", "[", 
   RowBox[{"x_", "+", "y_"}], "]"}], ":=", 
  RowBox[{
   RowBox[{"ExpandNCM", "[", "x", "]"}], "+", 
   RowBox[{"ExpandNCM", "[", "y", "]"}]}]}]], "Input",ExpressionUUID->\
"44b438e9-d7d0-46f1-a26a-0f094a593ed0"],

Cell[BoxData[
 RowBox[{
  RowBox[{"ExpandNCM", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{"h", ":", "NonCommutativeMultiply"}], ")"}], "[", 
    RowBox[{"x___", ",", "y_Plus", ",", "z___"}], "]"}], "]"}], ":=", 
  RowBox[{"Distribute", "[", 
   RowBox[{
    RowBox[{"h", "[", 
     RowBox[{"x", ",", "y", ",", "z"}], "]"}], ",", "Plus", ",", "h", ",", 
    "Plus", ",", 
    RowBox[{
     RowBox[{"ExpandNCM", "[", 
      RowBox[{"h", "[", "##", "]"}], "]"}], "&"}]}], "]"}]}]], "Input",Express\
ionUUID->"32c5a07f-a275-45ce-ae5d-f6eaf5d04b98"],

Cell[BoxData[
 RowBox[{
  RowBox[{"ExpandNCM", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{"h", ":", "NonCommutativeMultiply"}], ")"}], "[", 
    RowBox[{"x___", ",", 
     RowBox[{
      RowBox[{"(", 
       RowBox[{"y1_", "?", "NumericQ"}], ")"}], "y2_"}], ",", "z___"}], "]"}],
    "]"}], ":=", 
  RowBox[{"y1", " ", 
   RowBox[{"ExpandNCM", "[", 
    RowBox[{"h", "[", 
     RowBox[{"x", ",", "y2", ",", "z"}], "]"}], "]"}]}]}]], "Input",Expression\
UUID->"8c16f1ce-61f3-46c8-93a9-c21c9861d11c"],

Cell[BoxData[
 RowBox[{
  RowBox[{"ExpandNCM", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{"h", ":", "NonCommutativeMultiply"}], ")"}], "[", "]"}], "]"}], ":=",
   "1"}]], "Input",ExpressionUUID->"ab9cc2de-a827-45ec-80cb-d921aba76197"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"definition", " ", "of", " ", "the", " ", "single"}], "-", 
    RowBox[{"particle", " ", 
     RowBox[{"Green", "'"}], "s", " ", "function", " ", 
     SubscriptBox["G", 
      RowBox[{"i", ",", "j", ",", "\[Sigma]"}]]}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{"ExpandNCM", "[", 
    RowBox[{
     RowBox[{
      SubscriptBox["c", 
       RowBox[{"i_", ",", "\[Sigma]_"}]], "[", "t_", "]"}], "**", 
     RowBox[{
      SubscriptBox["cd", 
       RowBox[{"j_", ",", "\[Tau]_"}]], "[", "s_", "]"}]}], "]"}], ":=", 
   RowBox[{
    RowBox[{"KroneckerDelta", "[", 
     RowBox[{"\[Sigma]", "-", "\[Tau]"}], "]"}], 
    RowBox[{
     SubscriptBox["G", 
      RowBox[{"i", ",", "j", ",", "\[Sigma]"}]], "[", 
     RowBox[{"t", ",", "s"}], "]"}]}]}]}]], "Input",ExpressionUUID->"4e8b3e40-\
cbcf-4172-a6c4-0429a6b91e2a"],

Cell[BoxData[
 RowBox[{
  RowBox[{"ExpandNCM", "[", 
   RowBox[{
    RowBox[{
     SubscriptBox["cd", 
      RowBox[{"j_", ",", "\[Tau]_"}]], "[", "s_", "]"}], "**", 
    RowBox[{
     SubscriptBox["c", 
      RowBox[{"i_", ",", "\[Sigma]_"}]], "[", "t_", "]"}]}], "]"}], ":=", 
  RowBox[{
   RowBox[{"KroneckerDelta", "[", 
    RowBox[{"\[Sigma]", "-", "\[Tau]"}], "]"}], 
   RowBox[{"(", 
    RowBox[{
     RowBox[{
      RowBox[{"KroneckerDelta", "[", 
       RowBox[{"i", "-", "j"}], "]"}], 
      RowBox[{"KroneckerDelta", "[", 
       RowBox[{"t", "-", "s"}], "]"}]}], "-", 
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"i", ",", "j", ",", "\[Sigma]"}]], "[", 
      RowBox[{"t", ",", "s"}], "]"}]}], ")"}]}]}]], "Input",ExpressionUUID->\
"f2c886cd-b96f-42b7-9563-9c0d8b657200"],

Cell[BoxData[
 RowBox[{"(*", " ", 
  RowBox[{"apply", " ", 
   RowBox[{"Wick", "'"}], "s", " ", "theorem", " ", "to", " ", "four", " ", 
   "fermionic", " ", "operators"}], " ", "*)"}]], "Input",ExpressionUUID->\
"839e042b-6efe-411d-9900-fc2d74ebbaaa"],

Cell[BoxData[
 RowBox[{
  RowBox[{"ExpandNCM", "[", 
   RowBox[{
    RowBox[{
     SubscriptBox["c", 
      RowBox[{"i1_", ",", "\[Sigma]1_"}]], "[", "t1_", "]"}], "**", 
    RowBox[{
     SubscriptBox["cd", 
      RowBox[{"j1_", ",", "\[Tau]1_"}]], "[", "s1_", "]"}], "**", 
    RowBox[{
     SubscriptBox["c", 
      RowBox[{"i2_", ",", "\[Sigma]2_"}]], "[", "t2_", "]"}], "**", 
    RowBox[{
     SubscriptBox["cd", 
      RowBox[{"j2_", ",", "\[Tau]2_"}]], "[", "s2_", "]"}]}], "]"}], ":=", 
  RowBox[{
   RowBox[{
    RowBox[{"ExpandNCM", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["c", 
        RowBox[{"i1", ",", "\[Sigma]1"}]], "[", "t1", "]"}], "**", 
      RowBox[{
       SubscriptBox["cd", 
        RowBox[{"j1", ",", "\[Tau]1"}]], "[", "s1", "]"}]}], "]"}], 
    RowBox[{"ExpandNCM", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["c", 
        RowBox[{"i2", ",", "\[Sigma]2"}]], "[", "t2", "]"}], "**", 
      RowBox[{
       SubscriptBox["cd", 
        RowBox[{"j2", ",", "\[Tau]2"}]], "[", "s2", "]"}]}], "]"}]}], "+", 
   RowBox[{
    RowBox[{"ExpandNCM", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["c", 
        RowBox[{"i1", ",", "\[Sigma]1"}]], "[", "t1", "]"}], "**", 
      RowBox[{
       SubscriptBox["cd", 
        RowBox[{"j2", ",", "\[Tau]2"}]], "[", "s2", "]"}]}], "]"}], 
    RowBox[{"ExpandNCM", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["cd", 
        RowBox[{"j1", ",", "\[Tau]1"}]], "[", "s1", "]"}], "**", 
      RowBox[{
       SubscriptBox["c", 
        RowBox[{"i2", ",", "\[Sigma]2"}]], "[", "t2", "]"}]}], 
     "]"}]}]}]}]], "Input",ExpressionUUID->"32bb2983-417a-41c4-86b8-\
0da088df5fd4"],

Cell[BoxData[
 RowBox[{
  RowBox[{"ExpandNCM", "[", 
   RowBox[{
    RowBox[{
     SubscriptBox["cd", 
      RowBox[{"j1_", ",", "\[Tau]1_"}]], "[", "s1_", "]"}], "**", 
    RowBox[{
     SubscriptBox["c", 
      RowBox[{"i1_", ",", "\[Sigma]1_"}]], "[", "t1_", "]"}], "**", 
    RowBox[{
     SubscriptBox["cd", 
      RowBox[{"j2_", ",", "\[Tau]2_"}]], "[", "s2_", "]"}], "**", 
    RowBox[{
     SubscriptBox["c", 
      RowBox[{"i2_", ",", "\[Sigma]2_"}]], "[", "t2_", "]"}]}], "]"}], ":=", 
  RowBox[{
   RowBox[{
    RowBox[{"ExpandNCM", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["cd", 
        RowBox[{"j1", ",", "\[Tau]1"}]], "[", "s1", "]"}], "**", 
      RowBox[{
       SubscriptBox["c", 
        RowBox[{"i1", ",", "\[Sigma]1"}]], "[", "t1", "]"}]}], "]"}], 
    RowBox[{"ExpandNCM", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["cd", 
        RowBox[{"j2", ",", "\[Tau]2"}]], "[", "s2", "]"}], "**", 
      RowBox[{
       SubscriptBox["c", 
        RowBox[{"i2", ",", "\[Sigma]2"}]], "[", "t2", "]"}]}], "]"}]}], "+", 
   RowBox[{
    RowBox[{"ExpandNCM", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["cd", 
        RowBox[{"j1", ",", "\[Tau]1"}]], "[", "s1", "]"}], "**", 
      RowBox[{
       SubscriptBox["c", 
        RowBox[{"i2", ",", "\[Sigma]2"}]], "[", "t2", "]"}]}], "]"}], 
    RowBox[{"ExpandNCM", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["c", 
        RowBox[{"i1", ",", "\[Sigma]1"}]], "[", "t1", "]"}], "**", 
      RowBox[{
       SubscriptBox["cd", 
        RowBox[{"j2", ",", "\[Tau]2"}]], "[", "s2", "]"}]}], 
     "]"}]}]}]}]], "Input",ExpressionUUID->"49e89b23-d9f3-47d1-922b-\
f0056862b500"],

Cell[BoxData[
 RowBox[{
  RowBox[{"ExpandNCM", "[", 
   RowBox[{
    RowBox[{
     SubscriptBox["c", 
      RowBox[{"i1_", ",", "\[Sigma]1_"}]], "[", "t1_", "]"}], "**", 
    RowBox[{
     SubscriptBox["c", 
      RowBox[{"i2_", ",", "\[Sigma]2_"}]], "[", "t2_", "]"}], "**", 
    RowBox[{
     SubscriptBox["cd", 
      RowBox[{"j1_", ",", "\[Tau]1_"}]], "[", "s1_", "]"}], "**", 
    RowBox[{
     SubscriptBox["cd", 
      RowBox[{"j2_", ",", "\[Tau]2_"}]], "[", "s2_", "]"}]}], "]"}], ":=", 
  RowBox[{
   RowBox[{
    RowBox[{"KroneckerDelta", "[", 
     RowBox[{"\[Sigma]2", "-", "\[Tau]1"}], "]"}], 
    RowBox[{"KroneckerDelta", "[", 
     RowBox[{"i2", "-", "j1"}], "]"}], 
    RowBox[{"KroneckerDelta", "[", 
     RowBox[{"t2", "-", "s1"}], "]"}], 
    RowBox[{"ExpandNCM", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["c", 
        RowBox[{"i1", ",", "\[Sigma]1"}]], "[", "t1", "]"}], "**", 
      RowBox[{
       SubscriptBox["cd", 
        RowBox[{"j2", ",", "\[Tau]2"}]], "[", "s2", "]"}]}], "]"}]}], "-", 
   RowBox[{"ExpandNCM", "[", 
    RowBox[{
     RowBox[{
      SubscriptBox["c", 
       RowBox[{"i1", ",", "\[Sigma]1"}]], "[", "t1", "]"}], "**", 
     RowBox[{
      SubscriptBox["cd", 
       RowBox[{"j1", ",", "\[Tau]1"}]], "[", "s1", "]"}], "**", 
     RowBox[{
      SubscriptBox["c", 
       RowBox[{"i2", ",", "\[Sigma]2"}]], "[", "t2", "]"}], "**", 
     RowBox[{
      SubscriptBox["cd", 
       RowBox[{"j2", ",", "\[Tau]2"}]], "[", "s2", "]"}]}], 
    "]"}]}]}]], "Input",ExpressionUUID->"81e44988-eec1-41f7-8e84-\
bf43e7153d78"]
}, Open  ]],

Cell[CellGroupData[{

Cell["Total particle number (charge) correlations", "Subsubsection",ExpressionUUID->"2bc99d75-87bf-4dc0-92e1-936596065360"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "example", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{"ExpandNCM", "[", 
   RowBox[{
    RowBox[{
     SubscriptBox["c", 
      RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}], "**", 
    RowBox[{
     SubscriptBox["cd", 
      RowBox[{"j", ",", "\[Sigma]"}]], "[", "s", "]"}]}], "]"}]}]], "Input",Ex\
pressionUUID->"d93ee796-f4bb-4e00-b8e0-1ed35e963262"],

Cell[BoxData[
 RowBox[{
  SubscriptBox["G", 
   RowBox[{"i", ",", "j", ",", "\[Sigma]"}]], "[", 
  RowBox[{"t", ",", "s"}], "]"}]], "Output",ExpressionUUID->"17e97fc6-a0e2-\
4e6b-b970-fdce8498b8df"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"example", ":", " ", 
    RowBox[{"number", " ", "operator"}]}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{"ExpandNCM", "[", 
   RowBox[{
    RowBox[{
     SubscriptBox["cd", 
      RowBox[{"i", ",", "\[Sigma]"}]], "[", "s", "]"}], "**", 
    RowBox[{
     SubscriptBox["c", 
      RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}]}], "]"}]}]], "Input",Ex\
pressionUUID->"db3e7010-bcb1-49c6-ac24-2d09c0cf0926"],

Cell[BoxData[
 RowBox[{
  TemplateBox[{RowBox[{
      RowBox[{"s", "-", "t"}]}]},
   "KroneckerDeltaSeq"], "-", 
  RowBox[{
   SubscriptBox["G", 
    RowBox[{"i", ",", "i", ",", "\[Sigma]"}]], "[", 
   RowBox[{"t", ",", "s"}], "]"}]}]], "Output",ExpressionUUID->"fbbccbfc-e479-\
4763-b03e-54d673cb9ff7"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["\[Rho]", "corr"], "=", 
   RowBox[{
    RowBox[{
     RowBox[{"(", 
      RowBox[{"2", "-", 
       RowBox[{
        SubscriptBox["G", 
         RowBox[{"i", ",", "i", ",", "up"}]], "[", 
        RowBox[{"t", ",", "t"}], "]"}], "-", 
       RowBox[{
        SubscriptBox["G", 
         RowBox[{"i", ",", "i", ",", "dn"}]], "[", 
        RowBox[{"t", ",", "t"}], "]"}]}], ")"}], 
     RowBox[{"(", 
      RowBox[{"2", "-", 
       RowBox[{
        SubscriptBox["G", 
         RowBox[{"j", ",", "j", ",", "up"}]], "[", 
        RowBox[{"s", ",", "s"}], "]"}], "-", 
       RowBox[{
        SubscriptBox["G", 
         RowBox[{"j", ",", "j", ",", "dn"}]], "[", 
        RowBox[{"s", ",", "s"}], "]"}]}], ")"}]}], "-", 
    RowBox[{
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"i", ",", "j", ",", "up"}]], "[", 
      RowBox[{"t", ",", "s"}], "]"}], 
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"j", ",", "i", ",", "up"}]], "[", 
      RowBox[{"s", ",", "t"}], "]"}]}], "-", 
    RowBox[{
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"i", ",", "j", ",", "dn"}]], "[", 
      RowBox[{"t", ",", "s"}], "]"}], 
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"j", ",", "i", ",", "dn"}]], "[", 
      RowBox[{"s", ",", "t"}], "]"}]}]}]}], ";"}]], "Input",ExpressionUUID->\
"3c77d700-7922-4693-b81a-894372d33ec0"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"check", " ", "case", " ", "i"}], "\[NotEqual]", "j"}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"ExpandNCM", "[", 
       RowBox[{
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            SubscriptBox["n", 
             RowBox[{"i", ",", "up"}]], "[", "t", "]"}], "+", 
           RowBox[{
            SubscriptBox["n", 
             RowBox[{"i", ",", "dn"}]], "[", "t", "]"}]}], ")"}], "**", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            SubscriptBox["n", 
             RowBox[{"j", ",", "up"}]], "[", "s", "]"}], "+", 
           RowBox[{
            SubscriptBox["n", 
             RowBox[{"j", ",", "dn"}]], "[", "s", "]"}]}], ")"}]}], "/.", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"i_", ",", "\[Sigma]_"}]], "[", "t_", "]"}], 
          "\[RuleDelayed]", 
          RowBox[{
           RowBox[{
            SubscriptBox["cd", 
             RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}], "**", 
           RowBox[{
            SubscriptBox["c", 
             RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}]}]}], "}"}]}], 
       "]"}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{
         RowBox[{"KroneckerDelta", "[", 
          RowBox[{"up", "-", "dn"}], "]"}], "\[Rule]", "0"}], ",", 
        RowBox[{
         RowBox[{"KroneckerDelta", "[", 
          RowBox[{"i", "-", "j"}], "]"}], "\[Rule]", "0"}]}], "}"}]}], ")"}], 
    "-", 
    SubscriptBox["\[Rho]", "corr"]}], "]"}]}]], "Input",ExpressionUUID->\
"b59d2b5c-6b6c-4d1d-ac06-3adfb6add2f2"],

Cell[BoxData["0"], "Output",ExpressionUUID->"9464890a-5ac0-4120-8db5-9248dd0dd1fd"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"check", " ", "case", " ", "t"}], "\[NotEqual]", "s"}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"ExpandNCM", "[", 
       RowBox[{
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            SubscriptBox["n", 
             RowBox[{"i", ",", "up"}]], "[", "t", "]"}], "+", 
           RowBox[{
            SubscriptBox["n", 
             RowBox[{"i", ",", "dn"}]], "[", "t", "]"}]}], ")"}], "**", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            SubscriptBox["n", 
             RowBox[{"j", ",", "up"}]], "[", "s", "]"}], "+", 
           RowBox[{
            SubscriptBox["n", 
             RowBox[{"j", ",", "dn"}]], "[", "s", "]"}]}], ")"}]}], "/.", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"i_", ",", "\[Sigma]_"}]], "[", "t_", "]"}], 
          "\[RuleDelayed]", 
          RowBox[{
           RowBox[{
            SubscriptBox["cd", 
             RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}], "**", 
           RowBox[{
            SubscriptBox["c", 
             RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}]}]}], "}"}]}], 
       "]"}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{
         RowBox[{"KroneckerDelta", "[", 
          RowBox[{"up", "-", "dn"}], "]"}], "\[Rule]", "0"}], ",", 
        RowBox[{
         RowBox[{"KroneckerDelta", "[", 
          RowBox[{"t", "-", "s"}], "]"}], "\[Rule]", "0"}]}], "}"}]}], ")"}], 
    "-", 
    SubscriptBox["\[Rho]", "corr"]}], "]"}]}]], "Input",ExpressionUUID->\
"da735f2c-61e9-4d51-8120-6c49737fb482"],

Cell[BoxData["0"], "Output",ExpressionUUID->"a5242d33-11d0-401a-9dd8-ba226de14a07"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"special", " ", "case", " ", "i"}], "\[Equal]", 
    RowBox[{"j", " ", "and", " ", "t"}], "\[Equal]", "s"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{"Expand", "[", 
   RowBox[{"FullSimplify", "[", 
    RowBox[{
     RowBox[{"ExpandNCM", "[", 
      RowBox[{
       RowBox[{
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"i", ",", "up"}]], "[", "t", "]"}], "+", 
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"i", ",", "dn"}]], "[", "t", "]"}]}], ")"}], "**", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"i", ",", "up"}]], "[", "t", "]"}], "+", 
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"i", ",", "dn"}]], "[", "t", "]"}]}], ")"}]}], "/.", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{
          SubscriptBox["n", 
           RowBox[{"i_", ",", "\[Sigma]_"}]], "[", "t_", "]"}], 
         "\[RuleDelayed]", 
         RowBox[{
          RowBox[{
           SubscriptBox["cd", 
            RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}], "**", 
          RowBox[{
           SubscriptBox["c", 
            RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}]}]}], "}"}]}], 
      "]"}], "/.", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"KroneckerDelta", "[", 
        RowBox[{"up", "-", "dn"}], "]"}], "\[Rule]", "0"}], "}"}]}], "]"}], 
   "]"}]}]], "Input",ExpressionUUID->"43f954af-94a5-4d11-bfa8-7c93a41c4461"],

Cell[BoxData[
 RowBox[{"4", "-", 
  RowBox[{"3", " ", 
   RowBox[{
    SubscriptBox["G", 
     RowBox[{"i", ",", "i", ",", "dn"}]], "[", 
    RowBox[{"t", ",", "t"}], "]"}]}], "-", 
  RowBox[{"3", " ", 
   RowBox[{
    SubscriptBox["G", 
     RowBox[{"i", ",", "i", ",", "up"}]], "[", 
    RowBox[{"t", ",", "t"}], "]"}]}], "+", 
  RowBox[{"2", " ", 
   RowBox[{
    SubscriptBox["G", 
     RowBox[{"i", ",", "i", ",", "dn"}]], "[", 
    RowBox[{"t", ",", "t"}], "]"}], " ", 
   RowBox[{
    SubscriptBox["G", 
     RowBox[{"i", ",", "i", ",", "up"}]], "[", 
    RowBox[{"t", ",", "t"}], "]"}]}]}]], "Output",ExpressionUUID->"e20bea6a-\
5b56-4e97-a7fb-c6058909be6d"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Spin-up and spin-down correlations", "Subsubsection",ExpressionUUID->"d74a3517-f66a-45ff-adb4-c09bccbbfef3"],

Cell[BoxData[{
 RowBox[{
  RowBox[{
   SubscriptBox["up", "corr"], "=", 
   RowBox[{
    RowBox[{
     RowBox[{"(", 
      RowBox[{"1", "-", 
       RowBox[{
        SubscriptBox["G", 
         RowBox[{"i", ",", "i", ",", "up"}]], "[", 
        RowBox[{"t", ",", "t"}], "]"}]}], ")"}], 
     RowBox[{"(", 
      RowBox[{"1", "-", 
       RowBox[{
        SubscriptBox["G", 
         RowBox[{"j", ",", "j", ",", "up"}]], "[", 
        RowBox[{"s", ",", "s"}], "]"}]}], ")"}]}], "-", 
    RowBox[{
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"i", ",", "j", ",", "up"}]], "[", 
      RowBox[{"t", ",", "s"}], "]"}], 
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"j", ",", "i", ",", "up"}]], "[", 
      RowBox[{"s", ",", "t"}], "]"}]}]}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   SubscriptBox["dn", "corr"], "=", 
   RowBox[{
    RowBox[{
     RowBox[{"(", 
      RowBox[{"1", "-", 
       RowBox[{
        SubscriptBox["G", 
         RowBox[{"i", ",", "i", ",", "dn"}]], "[", 
        RowBox[{"t", ",", "t"}], "]"}]}], ")"}], 
     RowBox[{"(", 
      RowBox[{"1", "-", 
       RowBox[{
        SubscriptBox["G", 
         RowBox[{"j", ",", "j", ",", "dn"}]], "[", 
        RowBox[{"s", ",", "s"}], "]"}]}], ")"}]}], "-", 
    RowBox[{
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"i", ",", "j", ",", "dn"}]], "[", 
      RowBox[{"t", ",", "s"}], "]"}], 
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"j", ",", "i", ",", "dn"}]], "[", 
      RowBox[{"s", ",", "t"}], "]"}]}]}]}], ";"}]}], "Input",ExpressionUUID->\
"3ee05fa0-4305-4fc5-931e-821c80e3f187"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"check", " ", "case", " ", "i"}], "\[NotEqual]", "j"}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"ExpandNCM", "[", 
       RowBox[{
        RowBox[{
         RowBox[{
          SubscriptBox["n", 
           RowBox[{"i", ",", "up"}]], "[", "t", "]"}], "**", 
         RowBox[{
          SubscriptBox["n", 
           RowBox[{"j", ",", "up"}]], "[", "s", "]"}]}], "/.", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"i_", ",", "\[Sigma]_"}]], "[", "t_", "]"}], 
          "\[RuleDelayed]", 
          RowBox[{
           RowBox[{
            SubscriptBox["cd", 
             RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}], "**", 
           RowBox[{
            SubscriptBox["c", 
             RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}]}]}], "}"}]}], 
       "]"}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"KroneckerDelta", "[", 
         RowBox[{"i", "-", "j"}], "]"}], "\[Rule]", "0"}], "}"}]}], ")"}], 
    "-", 
    SubscriptBox["up", "corr"]}], "]"}]}]], "Input",ExpressionUUID->"c70e2f82-\
1fbb-4cce-9135-8c74f8cc0564"],

Cell[BoxData["0"], "Output",ExpressionUUID->"5648fb8c-6dac-4027-a29a-86f5dff0fe98"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"check", " ", "case", " ", "t"}], "\[NotEqual]", "s"}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"ExpandNCM", "[", 
       RowBox[{
        RowBox[{
         RowBox[{
          SubscriptBox["n", 
           RowBox[{"i", ",", "up"}]], "[", "t", "]"}], "**", 
         RowBox[{
          SubscriptBox["n", 
           RowBox[{"j", ",", "up"}]], "[", "s", "]"}]}], "/.", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"i_", ",", "\[Sigma]_"}]], "[", "t_", "]"}], 
          "\[RuleDelayed]", 
          RowBox[{
           RowBox[{
            SubscriptBox["cd", 
             RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}], "**", 
           RowBox[{
            SubscriptBox["c", 
             RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}]}]}], "}"}]}], 
       "]"}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"KroneckerDelta", "[", 
         RowBox[{"t", "-", "s"}], "]"}], "\[Rule]", "0"}], "}"}]}], ")"}], 
    "-", 
    SubscriptBox["up", "corr"]}], "]"}]}]], "Input",ExpressionUUID->"2399a700-\
f5ad-4d65-9feb-29a87d4f00ec"],

Cell[BoxData["0"], "Output",ExpressionUUID->"43ca4a33-be5a-4562-a6d5-9dc3a7cf9bfc"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"special", " ", "case", " ", "i"}], "\[Equal]", 
    RowBox[{"j", " ", "and", " ", "t"}], "\[Equal]", "s"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{"ExpandNCM", "[", 
    RowBox[{
     RowBox[{
      RowBox[{
       SubscriptBox["n", 
        RowBox[{"i", ",", "up"}]], "[", "t", "]"}], "**", 
      RowBox[{
       SubscriptBox["n", 
        RowBox[{"i", ",", "up"}]], "[", "t", "]"}]}], "/.", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{
        SubscriptBox["n", 
         RowBox[{"i_", ",", "\[Sigma]_"}]], "[", "t_", "]"}], 
       "\[RuleDelayed]", 
       RowBox[{
        RowBox[{
         SubscriptBox["cd", 
          RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}], "**", 
        RowBox[{
         SubscriptBox["c", 
          RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}]}]}], "}"}]}], 
    "]"}], "]"}]}]], "Input",ExpressionUUID->"50008626-57b9-4111-81cb-\
17670cdd7b38"],

Cell[BoxData[
 RowBox[{"1", "-", 
  RowBox[{
   SubscriptBox["G", 
    RowBox[{"i", ",", "i", ",", "up"}]], "[", 
   RowBox[{"t", ",", "t"}], "]"}]}]], "Output",ExpressionUUID->"4a25734f-b56c-\
4176-ba9f-9a229f6fc457"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"check", " ", "case", " ", "i"}], "\[NotEqual]", "j"}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"ExpandNCM", "[", 
       RowBox[{
        RowBox[{
         RowBox[{
          SubscriptBox["n", 
           RowBox[{"i", ",", "dn"}]], "[", "t", "]"}], "**", 
         RowBox[{
          SubscriptBox["n", 
           RowBox[{"j", ",", "dn"}]], "[", "s", "]"}]}], "/.", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"i_", ",", "\[Sigma]_"}]], "[", "t_", "]"}], 
          "\[RuleDelayed]", 
          RowBox[{
           RowBox[{
            SubscriptBox["cd", 
             RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}], "**", 
           RowBox[{
            SubscriptBox["c", 
             RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}]}]}], "}"}]}], 
       "]"}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"KroneckerDelta", "[", 
         RowBox[{"i", "-", "j"}], "]"}], "\[Rule]", "0"}], "}"}]}], ")"}], 
    "-", 
    SubscriptBox["dn", "corr"]}], "]"}]}]], "Input",ExpressionUUID->"6b1ef567-\
a761-442d-9828-9e5d8c0f76b6"],

Cell[BoxData["0"], "Output",ExpressionUUID->"f2890f77-661b-46a6-bdaf-b0b444a9f68d"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"check", " ", "case", " ", "t"}], "\[NotEqual]", "s"}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"ExpandNCM", "[", 
       RowBox[{
        RowBox[{
         RowBox[{
          SubscriptBox["n", 
           RowBox[{"i", ",", "dn"}]], "[", "t", "]"}], "**", 
         RowBox[{
          SubscriptBox["n", 
           RowBox[{"j", ",", "dn"}]], "[", "s", "]"}]}], "/.", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"i_", ",", "\[Sigma]_"}]], "[", "t_", "]"}], 
          "\[RuleDelayed]", 
          RowBox[{
           RowBox[{
            SubscriptBox["cd", 
             RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}], "**", 
           RowBox[{
            SubscriptBox["c", 
             RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}]}]}], "}"}]}], 
       "]"}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"KroneckerDelta", "[", 
         RowBox[{"t", "-", "s"}], "]"}], "\[Rule]", "0"}], "}"}]}], ")"}], 
    "-", 
    SubscriptBox["dn", "corr"]}], "]"}]}]], "Input",ExpressionUUID->"5db9708f-\
4837-480a-86ba-485eca0efbef"],

Cell[BoxData["0"], "Output",ExpressionUUID->"4c3ce712-410a-4360-82d4-88eaf248bdaa"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"special", " ", "case", " ", "i"}], "\[Equal]", 
    RowBox[{"j", " ", "and", " ", "t"}], "\[Equal]", "s"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{"ExpandNCM", "[", 
    RowBox[{
     RowBox[{
      RowBox[{
       SubscriptBox["n", 
        RowBox[{"i", ",", "dn"}]], "[", "t", "]"}], "**", 
      RowBox[{
       SubscriptBox["n", 
        RowBox[{"i", ",", "dn"}]], "[", "t", "]"}]}], "/.", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{
        SubscriptBox["n", 
         RowBox[{"i_", ",", "\[Sigma]_"}]], "[", "t_", "]"}], 
       "\[RuleDelayed]", 
       RowBox[{
        RowBox[{
         SubscriptBox["cd", 
          RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}], "**", 
        RowBox[{
         SubscriptBox["c", 
          RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}]}]}], "}"}]}], 
    "]"}], "]"}]}]], "Input",ExpressionUUID->"481af523-bce2-4b53-9ad3-\
686258ce9e29"],

Cell[BoxData[
 RowBox[{"1", "-", 
  RowBox[{
   SubscriptBox["G", 
    RowBox[{"i", ",", "i", ",", "dn"}]], "[", 
   RowBox[{"t", ",", "t"}], "]"}]}]], "Output",ExpressionUUID->"1559858c-45a9-\
4d47-b55d-7e92e1b4cbec"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"cross", " ", "terms"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["updown", "cross"], "=", 
    RowBox[{
     RowBox[{
      RowBox[{"(", 
       RowBox[{"1", "-", 
        RowBox[{
         SubscriptBox["G", 
          RowBox[{"i", ",", "i", ",", "up"}]], "[", 
         RowBox[{"t", ",", "t"}], "]"}]}], ")"}], 
      RowBox[{"(", 
       RowBox[{"1", "-", 
        RowBox[{
         SubscriptBox["G", 
          RowBox[{"j", ",", "j", ",", "dn"}]], "[", 
         RowBox[{"s", ",", "s"}], "]"}]}], ")"}]}], "+", 
     RowBox[{
      RowBox[{"(", 
       RowBox[{"1", "-", 
        RowBox[{
         SubscriptBox["G", 
          RowBox[{"i", ",", "i", ",", "dn"}]], "[", 
         RowBox[{"t", ",", "t"}], "]"}]}], ")"}], 
      RowBox[{"(", 
       RowBox[{"1", "-", 
        RowBox[{
         SubscriptBox["G", 
          RowBox[{"j", ",", "j", ",", "up"}]], "[", 
         RowBox[{"s", ",", "s"}], "]"}]}], ")"}]}]}]}], ";"}]}]], "Input",Expr\
essionUUID->"b48889d7-ae23-45aa-b3fc-fd0495a1f5c5"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "check", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"ExpandNCM", "[", 
       RowBox[{
        RowBox[{
         RowBox[{
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"i", ",", "up"}]], "[", "t", "]"}], "**", 
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"j", ",", "dn"}]], "[", "s", "]"}]}], "+", 
         RowBox[{
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"i", ",", "dn"}]], "[", "t", "]"}], "**", 
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"j", ",", "up"}]], "[", "s", "]"}]}]}], "/.", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"i_", ",", "\[Sigma]_"}]], "[", "t_", "]"}], 
          "\[RuleDelayed]", 
          RowBox[{
           RowBox[{
            SubscriptBox["cd", 
             RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}], "**", 
           RowBox[{
            SubscriptBox["c", 
             RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}]}]}], "}"}]}], 
       "]"}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"KroneckerDelta", "[", 
         RowBox[{"up", "-", "dn"}], "]"}], "\[Rule]", "0"}], "}"}]}], ")"}], 
    "-", 
    SubscriptBox["updown", "cross"]}], "]"}]}]], "Input",ExpressionUUID->\
"192838cf-098e-4a98-9a42-41716ed15600"],

Cell[BoxData["0"], "Output",ExpressionUUID->"dd7e5078-2863-4518-83ea-5299689dc2b8"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
   "sum", " ", "of", " ", "three", " ", "terms", " ", "gives", " ", "charge", 
    " ", "correlations"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      SubscriptBox["up", "corr"], "+", 
      SubscriptBox["dn", "corr"], "+", 
      SubscriptBox["updown", "cross"]}], ")"}], "-", 
    SubscriptBox["\[Rho]", "corr"]}], "]"}]}]], "Input",ExpressionUUID->\
"35174cfd-8dd2-44fd-9cc1-067797823c87"],

Cell[BoxData["0"], "Output",ExpressionUUID->"180f4cb4-24cf-4e69-bcd6-5fd05f2f7e12"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Spin difference correlations", "Subsubsection",ExpressionUUID->"7409eea4-b238-40c2-9be8-e0319f0504ce"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["spin", 
    RowBox[{"corr", ",", "z"}]], "=", 
   RowBox[{
    RowBox[{
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        SubscriptBox["G", 
         RowBox[{"i", ",", "i", ",", "up"}]], "[", 
        RowBox[{"t", ",", "t"}], "]"}], "-", 
       RowBox[{
        SubscriptBox["G", 
         RowBox[{"i", ",", "i", ",", "dn"}]], "[", 
        RowBox[{"t", ",", "t"}], "]"}]}], ")"}], 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        SubscriptBox["G", 
         RowBox[{"j", ",", "j", ",", "up"}]], "[", 
        RowBox[{"s", ",", "s"}], "]"}], "-", 
       RowBox[{
        SubscriptBox["G", 
         RowBox[{"j", ",", "j", ",", "dn"}]], "[", 
        RowBox[{"s", ",", "s"}], "]"}]}], ")"}]}], "-", 
    RowBox[{
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"i", ",", "j", ",", "up"}]], "[", 
      RowBox[{"t", ",", "s"}], "]"}], 
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"j", ",", "i", ",", "up"}]], "[", 
      RowBox[{"s", ",", "t"}], "]"}]}], "-", 
    RowBox[{
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"i", ",", "j", ",", "dn"}]], "[", 
      RowBox[{"t", ",", "s"}], "]"}], 
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"j", ",", "i", ",", "dn"}]], "[", 
      RowBox[{"s", ",", "t"}], "]"}]}]}]}], ";"}]], "Input",ExpressionUUID->\
"15327b87-0323-4dae-9ff9-59fe45e1038b"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"check", " ", "case", " ", "i"}], "\[NotEqual]", "j"}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"ExpandNCM", "[", 
       RowBox[{
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            SubscriptBox["n", 
             RowBox[{"i", ",", "up"}]], "[", "t", "]"}], "-", 
           RowBox[{
            SubscriptBox["n", 
             RowBox[{"i", ",", "dn"}]], "[", "t", "]"}]}], ")"}], "**", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            SubscriptBox["n", 
             RowBox[{"j", ",", "up"}]], "[", "s", "]"}], "-", 
           RowBox[{
            SubscriptBox["n", 
             RowBox[{"j", ",", "dn"}]], "[", "s", "]"}]}], ")"}]}], "/.", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"i_", ",", "\[Sigma]_"}]], "[", "t_", "]"}], 
          "\[RuleDelayed]", 
          RowBox[{
           RowBox[{
            SubscriptBox["cd", 
             RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}], "**", 
           RowBox[{
            SubscriptBox["c", 
             RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}]}]}], "}"}]}], 
       "]"}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{
         RowBox[{"KroneckerDelta", "[", 
          RowBox[{"up", "-", "dn"}], "]"}], "\[Rule]", "0"}], ",", 
        RowBox[{
         RowBox[{"KroneckerDelta", "[", 
          RowBox[{"i", "-", "j"}], "]"}], "\[Rule]", "0"}]}], "}"}]}], ")"}], 
    "-", 
    SubscriptBox["spin", 
     RowBox[{"corr", ",", "z"}]]}], "]"}]}]], "Input",ExpressionUUID->\
"2dee7dab-31d6-4033-9f82-8571f1e1f9b2"],

Cell[BoxData["0"], "Output",ExpressionUUID->"5985ad9f-b7b3-42db-90f4-8672a6a0427d"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"check", " ", "case", " ", "t"}], "\[NotEqual]", "s"}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"ExpandNCM", "[", 
       RowBox[{
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            SubscriptBox["n", 
             RowBox[{"i", ",", "up"}]], "[", "t", "]"}], "-", 
           RowBox[{
            SubscriptBox["n", 
             RowBox[{"i", ",", "dn"}]], "[", "t", "]"}]}], ")"}], "**", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            SubscriptBox["n", 
             RowBox[{"j", ",", "up"}]], "[", "s", "]"}], "-", 
           RowBox[{
            SubscriptBox["n", 
             RowBox[{"j", ",", "dn"}]], "[", "s", "]"}]}], ")"}]}], "/.", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"i_", ",", "\[Sigma]_"}]], "[", "t_", "]"}], 
          "\[RuleDelayed]", 
          RowBox[{
           RowBox[{
            SubscriptBox["cd", 
             RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}], "**", 
           RowBox[{
            SubscriptBox["c", 
             RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}]}]}], "}"}]}], 
       "]"}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{
         RowBox[{"KroneckerDelta", "[", 
          RowBox[{"up", "-", "dn"}], "]"}], "\[Rule]", "0"}], ",", 
        RowBox[{
         RowBox[{"KroneckerDelta", "[", 
          RowBox[{"t", "-", "s"}], "]"}], "\[Rule]", "0"}]}], "}"}]}], ")"}], 
    "-", 
    SubscriptBox["spin", 
     RowBox[{"corr", ",", "z"}]]}], "]"}]}]], "Input",ExpressionUUID->\
"e9b5e707-25dd-493e-801a-14b416c7bdd6"],

Cell[BoxData["0"], "Output",ExpressionUUID->"3f0c8c17-e3f7-4ced-8c9f-f1e006a6a834"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"special", " ", "case", " ", "i"}], "\[Equal]", 
    RowBox[{"j", " ", "and", " ", "t"}], "\[Equal]", "s"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{"Expand", "[", 
   RowBox[{"FullSimplify", "[", 
    RowBox[{
     RowBox[{"ExpandNCM", "[", 
      RowBox[{
       RowBox[{
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"i", ",", "up"}]], "[", "t", "]"}], "-", 
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"i", ",", "dn"}]], "[", "t", "]"}]}], ")"}], "**", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"i", ",", "up"}]], "[", "t", "]"}], "-", 
          RowBox[{
           SubscriptBox["n", 
            RowBox[{"i", ",", "dn"}]], "[", "t", "]"}]}], ")"}]}], "/.", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{
          SubscriptBox["n", 
           RowBox[{"i_", ",", "\[Sigma]_"}]], "[", "t_", "]"}], 
         "\[RuleDelayed]", 
         RowBox[{
          RowBox[{
           SubscriptBox["cd", 
            RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}], "**", 
          RowBox[{
           SubscriptBox["c", 
            RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}]}]}], "}"}]}], 
      "]"}], "/.", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"KroneckerDelta", "[", 
        RowBox[{"up", "-", "dn"}], "]"}], "\[Rule]", "0"}], "}"}]}], "]"}], 
   "]"}]}]], "Input",ExpressionUUID->"fda87cf2-4d3f-47f8-85ae-946b1c1d5553"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["G", 
    RowBox[{"i", ",", "i", ",", "dn"}]], "[", 
   RowBox[{"t", ",", "t"}], "]"}], "+", 
  RowBox[{
   SubscriptBox["G", 
    RowBox[{"i", ",", "i", ",", "up"}]], "[", 
   RowBox[{"t", ",", "t"}], "]"}], "-", 
  RowBox[{"2", " ", 
   RowBox[{
    SubscriptBox["G", 
     RowBox[{"i", ",", "i", ",", "dn"}]], "[", 
    RowBox[{"t", ",", "t"}], "]"}], " ", 
   RowBox[{
    SubscriptBox["G", 
     RowBox[{"i", ",", "i", ",", "up"}]], "[", 
    RowBox[{"t", ",", "t"}], "]"}]}]}]], "Output",ExpressionUUID->"7ce4942e-\
2c31-4f50-ab6b-f9fefcfbef13"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["spin", 
    RowBox[{"corr", ",", "x"}]], "=", 
   RowBox[{
    RowBox[{
     RowBox[{"-", 
      RowBox[{
       SubscriptBox["G", 
        RowBox[{"i", ",", "j", ",", "up"}]], "[", 
       RowBox[{"t", ",", "s"}], "]"}]}], 
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"j", ",", "i", ",", "dn"}]], "[", 
      RowBox[{"s", ",", "t"}], "]"}]}], "-", 
    RowBox[{
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"i", ",", "j", ",", "dn"}]], "[", 
      RowBox[{"t", ",", "s"}], "]"}], 
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"j", ",", "i", ",", "up"}]], "[", 
      RowBox[{"s", ",", "t"}], "]"}]}]}]}], ";"}]], "Input",ExpressionUUID->\
"d0938527-7f5e-470f-96da-4b2162531716"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"check", " ", "case", " ", "i"}], "\[NotEqual]", "j"}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"ExpandNCM", "[", 
       RowBox[{
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            SubscriptBox["x", 
             RowBox[{"i", ",", "1"}]], "[", "t", "]"}], "+", 
           RowBox[{
            SubscriptBox["x", 
             RowBox[{"i", ",", 
              RowBox[{"-", "1"}]}]], "[", "t", "]"}]}], ")"}], "**", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            SubscriptBox["x", 
             RowBox[{"j", ",", "1"}]], "[", "s", "]"}], "+", 
           RowBox[{
            SubscriptBox["x", 
             RowBox[{"j", ",", 
              RowBox[{"-", "1"}]}]], "[", "s", "]"}]}], ")"}]}], "/.", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{
           RowBox[{
            SubscriptBox["x", 
             RowBox[{"i_", ",", "1"}]], "[", "t_", "]"}], "\[RuleDelayed]", 
           RowBox[{
            RowBox[{
             SubscriptBox["cd", 
              RowBox[{"i", ",", "dn"}]], "[", "t", "]"}], "**", 
            RowBox[{
             SubscriptBox["c", 
              RowBox[{"i", ",", "up"}]], "[", "t", "]"}]}]}], ",", 
          RowBox[{
           RowBox[{
            SubscriptBox["x", 
             RowBox[{"i_", ",", 
              RowBox[{"-", "1"}]}]], "[", "t_", "]"}], "\[RuleDelayed]", 
           RowBox[{
            RowBox[{
             SubscriptBox["cd", 
              RowBox[{"i", ",", "up"}]], "[", "t", "]"}], "**", 
            RowBox[{
             SubscriptBox["c", 
              RowBox[{"i", ",", "dn"}]], "[", "t", "]"}]}]}]}], "}"}]}], 
       "]"}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{
         RowBox[{"KroneckerDelta", "[", 
          RowBox[{"up", "-", "dn"}], "]"}], "\[Rule]", "0"}], ",", 
        RowBox[{
         RowBox[{"KroneckerDelta", "[", 
          RowBox[{"i", "-", "j"}], "]"}], "\[Rule]", "0"}]}], "}"}]}], ")"}], 
    "-", 
    SubscriptBox["spin", 
     RowBox[{"corr", ",", "x"}]]}], "]"}]}]], "Input",ExpressionUUID->\
"791e684e-2697-41d9-ab96-127eccd848d6"],

Cell[BoxData["0"], "Output",ExpressionUUID->"f0ab60d8-4b69-4f75-9140-3cc0e0312589"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"check", " ", "case", " ", "t"}], "\[NotEqual]", "s"}], " ", 
   "*)"}], "\[IndentingNewLine]", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"ExpandNCM", "[", 
       RowBox[{
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            SubscriptBox["x", 
             RowBox[{"i", ",", "1"}]], "[", "t", "]"}], "+", 
           RowBox[{
            SubscriptBox["x", 
             RowBox[{"i", ",", 
              RowBox[{"-", "1"}]}]], "[", "t", "]"}]}], ")"}], "**", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            SubscriptBox["x", 
             RowBox[{"j", ",", "1"}]], "[", "s", "]"}], "+", 
           RowBox[{
            SubscriptBox["x", 
             RowBox[{"j", ",", 
              RowBox[{"-", "1"}]}]], "[", "s", "]"}]}], ")"}]}], "/.", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{
           RowBox[{
            SubscriptBox["x", 
             RowBox[{"i_", ",", "1"}]], "[", "t_", "]"}], "\[RuleDelayed]", 
           RowBox[{
            RowBox[{
             SubscriptBox["cd", 
              RowBox[{"i", ",", "dn"}]], "[", "t", "]"}], "**", 
            RowBox[{
             SubscriptBox["c", 
              RowBox[{"i", ",", "up"}]], "[", "t", "]"}]}]}], ",", 
          RowBox[{
           RowBox[{
            SubscriptBox["x", 
             RowBox[{"i_", ",", 
              RowBox[{"-", "1"}]}]], "[", "t_", "]"}], "\[RuleDelayed]", 
           RowBox[{
            RowBox[{
             SubscriptBox["cd", 
              RowBox[{"i", ",", "up"}]], "[", "t", "]"}], "**", 
            RowBox[{
             SubscriptBox["c", 
              RowBox[{"i", ",", "dn"}]], "[", "t", "]"}]}]}]}], "}"}]}], 
       "]"}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{
         RowBox[{"KroneckerDelta", "[", 
          RowBox[{"up", "-", "dn"}], "]"}], "\[Rule]", "0"}], ",", 
        RowBox[{
         RowBox[{"KroneckerDelta", "[", 
          RowBox[{"t", "-", "s"}], "]"}], "\[Rule]", "0"}]}], "}"}]}], ")"}], 
    "-", 
    SubscriptBox["spin", 
     RowBox[{"corr", ",", "x"}]]}], "]"}]}]], "Input",ExpressionUUID->\
"5fb9c248-d3c9-4f31-a187-992d5741ce83"],

Cell[BoxData["0"], "Output",ExpressionUUID->"2b2688df-2c87-4ab1-b5b8-6f07467894c3"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{"special", " ", "case", " ", "i"}], "\[Equal]", 
    RowBox[{"j", " ", "and", " ", "t"}], "\[Equal]", "s"}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{"Expand", "[", 
   RowBox[{"FullSimplify", "[", 
    RowBox[{
     RowBox[{"ExpandNCM", "[", 
      RowBox[{
       RowBox[{
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           SubscriptBox["x", 
            RowBox[{"i", ",", "1"}]], "[", "t", "]"}], "+", 
          RowBox[{
           SubscriptBox["x", 
            RowBox[{"i", ",", 
             RowBox[{"-", "1"}]}]], "[", "t", "]"}]}], ")"}], "**", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           SubscriptBox["x", 
            RowBox[{"i", ",", "1"}]], "[", "t", "]"}], "+", 
          RowBox[{
           SubscriptBox["x", 
            RowBox[{"i", ",", 
             RowBox[{"-", "1"}]}]], "[", "t", "]"}]}], ")"}]}], "/.", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{
          RowBox[{
           SubscriptBox["x", 
            RowBox[{"i_", ",", "1"}]], "[", "t_", "]"}], "\[RuleDelayed]", 
          RowBox[{
           RowBox[{
            SubscriptBox["cd", 
             RowBox[{"i", ",", "dn"}]], "[", "t", "]"}], "**", 
           RowBox[{
            SubscriptBox["c", 
             RowBox[{"i", ",", "up"}]], "[", "t", "]"}]}]}], ",", 
         RowBox[{
          RowBox[{
           SubscriptBox["x", 
            RowBox[{"i_", ",", 
             RowBox[{"-", "1"}]}]], "[", "t_", "]"}], "\[RuleDelayed]", 
          RowBox[{
           RowBox[{
            SubscriptBox["cd", 
             RowBox[{"i", ",", "up"}]], "[", "t", "]"}], "**", 
           RowBox[{
            SubscriptBox["c", 
             RowBox[{"i", ",", "dn"}]], "[", "t", "]"}]}]}]}], "}"}]}], "]"}],
      "/.", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"KroneckerDelta", "[", 
        RowBox[{"up", "-", "dn"}], "]"}], "\[Rule]", "0"}], "}"}]}], "]"}], 
   "]"}]}]], "Input",ExpressionUUID->"ad3a100b-01b0-415c-b08a-3528accb96d6"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["G", 
    RowBox[{"i", ",", "i", ",", "dn"}]], "[", 
   RowBox[{"t", ",", "t"}], "]"}], "+", 
  RowBox[{
   SubscriptBox["G", 
    RowBox[{"i", ",", "i", ",", "up"}]], "[", 
   RowBox[{"t", ",", "t"}], "]"}], "-", 
  RowBox[{"2", " ", 
   RowBox[{
    SubscriptBox["G", 
     RowBox[{"i", ",", "i", ",", "dn"}]], "[", 
    RowBox[{"t", ",", "t"}], "]"}], " ", 
   RowBox[{
    SubscriptBox["G", 
     RowBox[{"i", ",", "i", ",", "up"}]], "[", 
    RowBox[{"t", ",", "t"}], "]"}]}]}]], "Output",ExpressionUUID->"cf87e8e9-\
c082-4fa1-b161-7eeb4db3eac9"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
   "xx", " ", "and", " ", "yy", " ", "correlations", " ", "seem", " ", "to", 
    " ", "agree"}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"ExpandNCM", "[", 
       RowBox[{
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            SubscriptBox["y", 
             RowBox[{"i", ",", "1"}]], "[", "t", "]"}], "+", 
           RowBox[{
            SubscriptBox["y", 
             RowBox[{"i", ",", 
              RowBox[{"-", "1"}]}]], "[", "t", "]"}]}], ")"}], "**", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            SubscriptBox["y", 
             RowBox[{"j", ",", "1"}]], "[", "s", "]"}], "+", 
           RowBox[{
            SubscriptBox["y", 
             RowBox[{"j", ",", 
              RowBox[{"-", "1"}]}]], "[", "s", "]"}]}], ")"}]}], "/.", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{
           RowBox[{
            SubscriptBox["y", 
             RowBox[{"i_", ",", "1"}]], "[", "t_", "]"}], "\[RuleDelayed]", 
           RowBox[{"\[ImaginaryI]", " ", 
            RowBox[{
             RowBox[{
              SubscriptBox["cd", 
               RowBox[{"i", ",", "dn"}]], "[", "t", "]"}], "**", 
             RowBox[{
              SubscriptBox["c", 
               RowBox[{"i", ",", "up"}]], "[", "t", "]"}]}]}]}], ",", 
          RowBox[{
           RowBox[{
            SubscriptBox["y", 
             RowBox[{"i_", ",", 
              RowBox[{"-", "1"}]}]], "[", "t_", "]"}], "\[RuleDelayed]", 
           RowBox[{
            RowBox[{"-", "\[ImaginaryI]"}], " ", 
            RowBox[{
             RowBox[{
              SubscriptBox["cd", 
               RowBox[{"i", ",", "up"}]], "[", "t", "]"}], "**", 
             RowBox[{
              SubscriptBox["c", 
               RowBox[{"i", ",", "dn"}]], "[", "t", "]"}]}]}]}]}], "}"}]}], 
       "]"}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{
         RowBox[{"KroneckerDelta", "[", 
          RowBox[{"up", "-", "dn"}], "]"}], "\[Rule]", "0"}], ",", 
        RowBox[{
         RowBox[{"KroneckerDelta", "[", 
          RowBox[{"i", "-", "j"}], "]"}], "\[Rule]", "0"}]}], "}"}]}], ")"}], 
    "-", 
    SubscriptBox["spin", 
     RowBox[{"corr", ",", "x"}]]}], "]"}]}]], "Input",ExpressionUUID->\
"fbcc240c-d686-4565-b66d-606779d11e87"],

Cell[BoxData["0"], "Output",ExpressionUUID->"cfecf615-62f4-4d36-a62a-c4e176ea7d3f"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Superconducting susceptibilities", "Subsubsection",ExpressionUUID->"b2761392-2125-4a92-968f-1ad236858cbc"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"s", "-", 
    RowBox[{"wave", " ", "pairing"}]}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["sc", "sw"], "=", 
    RowBox[{
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"i", ",", "j", ",", "up"}]], "[", 
      RowBox[{"t", ",", "s"}], "]"}], 
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"i", ",", "j", ",", "dn"}]], "[", 
      RowBox[{"t", ",", "s"}], "]"}]}]}], ";"}]}]], "Input",ExpressionUUID->\
"890beb2a-0213-4cd3-a0f8-cd93143911d2"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "check", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"ExpandNCM", "[", 
       RowBox[{
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           SubscriptBox["c", 
            RowBox[{"i", ",", "dn"}]], "[", "t", "]"}], "**", 
          RowBox[{
           SubscriptBox["c", 
            RowBox[{"i", ",", "up"}]], "[", "t", "]"}]}], ")"}], "**", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           SubscriptBox["cd", 
            RowBox[{"j", ",", "up"}]], "[", "s", "]"}], "**", 
          RowBox[{
           SubscriptBox["cd", 
            RowBox[{"j", ",", "dn"}]], "[", "s", "]"}]}], ")"}]}], "]"}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"KroneckerDelta", "[", 
         RowBox[{"up", "-", "dn"}], "]"}], "\[Rule]", "0"}], "}"}]}], ")"}], 
    "-", 
    SubscriptBox["sc", "sw"]}], "]"}]}]], "Input",ExpressionUUID->"1bf83b97-\
5b9f-48c5-999a-f2a36d30193f"],

Cell[BoxData["0"], "Output",ExpressionUUID->"55302f1b-3d56-487e-858a-d86fa024568f"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"d", "-", 
    RowBox[{"wave", " ", "pairing"}]}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["sc", "dw"], "=", 
    RowBox[{
     FractionBox["1", "4"], " ", 
     RowBox[{
      SubscriptBox["G", 
       RowBox[{"i", ",", "j", ",", "up"}]], "[", 
      RowBox[{"t", ",", "s"}], "]"}], " ", 
     RowBox[{
      SubscriptBox["G", 
       RowBox[{
        RowBox[{"i", "+", "\[Delta]"}], ",", 
        RowBox[{"j", "+", "\[Epsilon]"}], ",", "dn"}]], "[", 
      RowBox[{"t", ",", "s"}], "]"}]}]}], ";"}]}]], "Input",ExpressionUUID->\
"658ce7f9-745e-48a8-ad94-7bf5ac15e484"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "check", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{
       FractionBox["1", "4"], 
       RowBox[{"ExpandNCM", "[", 
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            SubscriptBox["c", 
             RowBox[{
              RowBox[{"i", "+", "\[Delta]"}], ",", "dn"}]], "[", "t", "]"}], "**", 
           RowBox[{
            SubscriptBox["c", 
             RowBox[{"i", ",", "up"}]], "[", "t", "]"}]}], ")"}], "**", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            SubscriptBox["cd", 
             RowBox[{"j", ",", "up"}]], "[", "s", "]"}], "**", 
           RowBox[{
            SubscriptBox["cd", 
             RowBox[{
              RowBox[{"j", "+", "\[Epsilon]"}], ",", "dn"}]], "[", "s", 
            "]"}]}], ")"}]}], "]"}]}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"KroneckerDelta", "[", 
         RowBox[{"up", "-", "dn"}], "]"}], "\[Rule]", "0"}], "}"}]}], ")"}], 
    "-", 
    SubscriptBox["sc", "dw"]}], "]"}]}]], "Input",ExpressionUUID->"348c8a86-\
0fc5-45c0-99c6-4afa87a8169f"],

Cell[BoxData["0"], "Output",ExpressionUUID->"cd22a849-6aad-4d9d-968f-10e7db1beb42"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Current correlations", "Subsubsection",ExpressionUUID->"9b971066-7c7c-4378-b39e-f08f836f4af5"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["j", 
    RowBox[{"corr", ",", "diag"}]], "=", 
   RowBox[{
    RowBox[{"-", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        SubscriptBox["G", 
         RowBox[{"i\[Delta]", ",", "i", ",", "up"}]], "[", 
        RowBox[{"t", ",", "t"}], "]"}], "-", 
       RowBox[{
        SubscriptBox["G", 
         RowBox[{"i", ",", "i\[Delta]", ",", "up"}]], "[", 
        RowBox[{"t", ",", "t"}], "]"}], "+", 
       RowBox[{
        SubscriptBox["G", 
         RowBox[{"i\[Delta]", ",", "i", ",", "dn"}]], "[", 
        RowBox[{"t", ",", "t"}], "]"}], "-", 
       RowBox[{
        SubscriptBox["G", 
         RowBox[{"i", ",", "i\[Delta]", ",", "dn"}]], "[", 
        RowBox[{"t", ",", "t"}], "]"}]}], ")"}]}], " ", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{
       SubscriptBox["G", 
        RowBox[{"j\[Epsilon]", ",", "j", ",", "up"}]], "[", 
       RowBox[{"s", ",", "s"}], "]"}], "-", 
      RowBox[{
       SubscriptBox["G", 
        RowBox[{"j", ",", "j\[Epsilon]", ",", "up"}]], "[", 
       RowBox[{"s", ",", "s"}], "]"}], "+", 
      RowBox[{
       SubscriptBox["G", 
        RowBox[{"j\[Epsilon]", ",", "j", ",", "dn"}]], "[", 
       RowBox[{"s", ",", "s"}], "]"}], "-", 
      RowBox[{
       SubscriptBox["G", 
        RowBox[{"j", ",", "j\[Epsilon]", ",", "dn"}]], "[", 
       RowBox[{"s", ",", "s"}], "]"}]}], ")"}]}]}], ";"}]], "Input",Expression\
UUID->"68783b49-a8ab-4e30-a244-128a73712eed"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["j", 
    RowBox[{"corr", ",", "offdiag"}]], "=", 
   RowBox[{"Sum", "[", 
    RowBox[{
     RowBox[{
      RowBox[{
       RowBox[{
        SubscriptBox["G", 
         RowBox[{"i\[Delta]", ",", "j\[Epsilon]", ",", "\[Sigma]"}]], "[", 
        RowBox[{"t", ",", "s"}], "]"}], 
       RowBox[{"(", 
        RowBox[{
         RowBox[{
          TemplateBox[{RowBox[{"t", "-", "s"}]},
           "KroneckerDeltaSeq"], 
          TemplateBox[{RowBox[{
              RowBox[{"i", "-", "j"}]}]},
           "KroneckerDeltaSeq"]}], "-", 
         RowBox[{
          SubscriptBox["G", 
           RowBox[{"j", ",", "i", ",", "\[Sigma]"}]], "[", 
          RowBox[{"s", ",", "t"}], "]"}]}], ")"}]}], "-", 
      RowBox[{
       RowBox[{
        SubscriptBox["G", 
         RowBox[{"i\[Delta]", ",", "j", ",", "\[Sigma]"}]], "[", 
        RowBox[{"t", ",", "s"}], "]"}], 
       RowBox[{"(", 
        RowBox[{
         RowBox[{
          TemplateBox[{RowBox[{"t", "-", "s"}]},
           "KroneckerDeltaSeq"], 
          TemplateBox[{RowBox[{"i", "-", "j\[Epsilon]"}]},
           "KroneckerDeltaSeq"]}], "-", 
         RowBox[{
          SubscriptBox["G", 
           RowBox[{"j\[Epsilon]", ",", "i", ",", "\[Sigma]"}]], "[", 
          RowBox[{"s", ",", "t"}], "]"}]}], ")"}]}], "-", 
      RowBox[{
       RowBox[{
        SubscriptBox["G", 
         RowBox[{"i", ",", "j\[Epsilon]", ",", "\[Sigma]"}]], "[", 
        RowBox[{"t", ",", "s"}], "]"}], 
       RowBox[{"(", 
        RowBox[{
         RowBox[{
          TemplateBox[{RowBox[{"t", "-", "s"}]},
           "KroneckerDeltaSeq"], 
          TemplateBox[{RowBox[{"i\[Delta]", "-", "j"}]},
           "KroneckerDeltaSeq"]}], "-", 
         RowBox[{
          SubscriptBox["G", 
           RowBox[{"j", ",", "i\[Delta]", ",", "\[Sigma]"}]], "[", 
          RowBox[{"s", ",", "t"}], "]"}]}], ")"}]}], "+", 
      RowBox[{
       RowBox[{
        SubscriptBox["G", 
         RowBox[{"i", ",", "j", ",", "\[Sigma]"}]], "[", 
        RowBox[{"t", ",", "s"}], "]"}], 
       RowBox[{"(", 
        RowBox[{
         RowBox[{
          TemplateBox[{RowBox[{"t", "-", "s"}]},
           "KroneckerDeltaSeq"], 
          TemplateBox[{RowBox[{"i\[Delta]", "-", "j\[Epsilon]"}]},
           "KroneckerDeltaSeq"]}], "-", 
         RowBox[{
          SubscriptBox["G", 
           RowBox[{"j\[Epsilon]", ",", "i\[Delta]", ",", "\[Sigma]"}]], "[", 
          RowBox[{"s", ",", "t"}], "]"}]}], ")"}]}]}], ",", 
     RowBox[{"{", 
      RowBox[{"\[Sigma]", ",", 
       RowBox[{"{", 
        RowBox[{"up", ",", "dn"}], "}"}]}], "}"}]}], "]"}]}], ";"}]], "Input",\
ExpressionUUID->"5967b3ee-c567-4558-9a6c-22f63c7ade8f"],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["j", "corr"], "=", 
   RowBox[{
    SubscriptBox["j", 
     RowBox[{"corr", ",", "diag"}]], "+", 
    SubscriptBox["j", 
     RowBox[{"corr", ",", "offdiag"}]]}]}], ";"}]], "Input",ExpressionUUID->\
"2c309710-7858-48d6-8c59-bc23d96a63ab"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "check", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"ExpandNCM", "[", 
       RowBox[{
        RowBox[{
         RowBox[{
          SubscriptBox["J", 
           RowBox[{"i\[Delta]", ",", "i"}]], "[", "t", "]"}], "**", 
         RowBox[{
          SubscriptBox["J", 
           RowBox[{"j\[Epsilon]", ",", "j"}]], "[", "s", "]"}]}], "/.", 
        RowBox[{"{", 
         RowBox[{
          RowBox[{
           SubscriptBox["J", 
            RowBox[{"i_", ",", "j_"}]], "[", "t_", "]"}], "\[RuleDelayed]", 
          RowBox[{"\[ImaginaryI]", " ", 
           RowBox[{"Sum", "[", 
            RowBox[{
             RowBox[{
              RowBox[{
               RowBox[{
                SubscriptBox["cd", 
                 RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}], "**", 
               RowBox[{
                SubscriptBox["c", 
                 RowBox[{"j", ",", "\[Sigma]"}]], "[", "t", "]"}]}], "-", 
              RowBox[{
               RowBox[{
                SubscriptBox["cd", 
                 RowBox[{"j", ",", "\[Sigma]"}]], "[", "t", "]"}], "**", 
               RowBox[{
                SubscriptBox["c", 
                 RowBox[{"i", ",", "\[Sigma]"}]], "[", "t", "]"}]}]}], ",", 
             RowBox[{"{", 
              RowBox[{"\[Sigma]", ",", 
               RowBox[{"{", 
                RowBox[{"up", ",", "dn"}], "}"}]}], "}"}]}], "]"}]}]}], 
         "}"}]}], "]"}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"KroneckerDelta", "[", 
         RowBox[{"up", "-", "dn"}], "]"}], "\[Rule]", "0"}], "}"}]}], ")"}], 
    "-", 
    SubscriptBox["j", "corr"]}], "]"}]}]], "Input",ExpressionUUID->"87a470ef-\
ae01-4cae-b137-265b7004ec15"],

Cell[BoxData["0"], "Output",ExpressionUUID->"c5987598-0702-4c45-b478-4f9dfa0168df"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Raman response", "Subsubsection",ExpressionUUID->"5e41d14f-6a91-40bd-823e-9bc8988d4949"],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"Raman", " ", 
    SubscriptBox["B", 
     RowBox[{"1", "g"}]]}], " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    SubscriptBox["raman", "B1g"], "=", 
    RowBox[{
     FractionBox["1", "16"], 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           SubscriptBox["G", 
            RowBox[{"i", ",", 
             RowBox[{"i", "+", "\[Delta]"}], ",", "up"}]], "[", 
           RowBox[{"t", ",", "t"}], "]"}], "+", 
          RowBox[{
           SubscriptBox["G", 
            RowBox[{"i", ",", 
             RowBox[{"i", "+", "\[Delta]"}], ",", "dn"}]], "[", 
           RowBox[{"t", ",", "t"}], "]"}]}], ")"}], " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           SubscriptBox["G", 
            RowBox[{"j", ",", 
             RowBox[{"j", "+", "\[Epsilon]"}], ",", "up"}]], "[", 
           RowBox[{"s", ",", "s"}], "]"}], "+", 
          RowBox[{
           SubscriptBox["G", 
            RowBox[{"j", ",", 
             RowBox[{"j", "+", "\[Epsilon]"}], ",", "dn"}]], "[", 
           RowBox[{"s", ",", "s"}], "]"}]}], ")"}]}], "+", 
       RowBox[{
        RowBox[{
         SubscriptBox["G", 
          RowBox[{"i", ",", 
           RowBox[{"j", "+", "\[Epsilon]"}], ",", "up"}]], "[", 
         RowBox[{"t", ",", "s"}], "]"}], 
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           TemplateBox[{RowBox[{"t", "-", "s"}]},
            "KroneckerDeltaSeq"], " ", 
           TemplateBox[{RowBox[{"i", "+", "\[Delta]", "-", "j"}]},
            "KroneckerDeltaSeq"]}], "-", 
          RowBox[{
           SubscriptBox["G", 
            RowBox[{"j", ",", 
             RowBox[{"i", "+", "\[Delta]"}], ",", "up"}]], "[", 
           RowBox[{"s", ",", "t"}], "]"}]}], ")"}]}], "+", 
       RowBox[{
        RowBox[{
         SubscriptBox["G", 
          RowBox[{"i", ",", 
           RowBox[{"j", "+", "\[Epsilon]"}], ",", "dn"}]], "[", 
         RowBox[{"t", ",", "s"}], "]"}], 
        RowBox[{"(", 
         RowBox[{
          RowBox[{
           TemplateBox[{RowBox[{"t", "-", "s"}]},
            "KroneckerDeltaSeq"], " ", 
           TemplateBox[{RowBox[{"i", "+", "\[Delta]", "-", "j"}]},
            "KroneckerDeltaSeq"]}], "-", 
          RowBox[{
           SubscriptBox["G", 
            RowBox[{"j", ",", 
             RowBox[{"i", "+", "\[Delta]"}], ",", "dn"}]], "[", 
           RowBox[{"s", ",", "t"}], "]"}]}], ")"}]}]}], ")"}]}]}], 
   ";"}]}]], "Input",ExpressionUUID->"5d64089e-66aa-43ff-83b3-d72b319a8676"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", "check", " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{
       FractionBox["1", "16"], 
       RowBox[{"ExpandNCM", "[", 
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            RowBox[{
             SubscriptBox["cd", 
              RowBox[{
               RowBox[{"i", "+", "\[Delta]"}], ",", "up"}]], "[", "t", "]"}], 
            "**", 
            RowBox[{
             SubscriptBox["c", 
              RowBox[{"i", ",", "up"}]], "[", "t", "]"}]}], "+", 
           RowBox[{
            RowBox[{
             SubscriptBox["cd", 
              RowBox[{
               RowBox[{"i", "+", "\[Delta]"}], ",", "dn"}]], "[", "t", "]"}], 
            "**", 
            RowBox[{
             SubscriptBox["c", 
              RowBox[{"i", ",", "dn"}]], "[", "t", "]"}]}]}], ")"}], "**", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{
            RowBox[{
             SubscriptBox["cd", 
              RowBox[{
               RowBox[{"j", "+", "\[Epsilon]"}], ",", "up"}]], "[", "s", 
             "]"}], "**", 
            RowBox[{
             SubscriptBox["c", 
              RowBox[{"j", ",", "up"}]], "[", "s", "]"}]}], "+", 
           RowBox[{
            RowBox[{
             SubscriptBox["cd", 
              RowBox[{
               RowBox[{"j", "+", "\[Epsilon]"}], ",", "dn"}]], "[", "s", 
             "]"}], "**", 
            RowBox[{
             SubscriptBox["c", 
              RowBox[{"j", ",", "dn"}]], "[", "s", "]"}]}]}], ")"}]}], 
        "]"}]}], "/.", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{
         RowBox[{"KroneckerDelta", "[", 
          RowBox[{"up", "-", "dn"}], "]"}], "\[Rule]", "0"}], ",", 
        RowBox[{
         RowBox[{"KroneckerDelta", "[", "\[Delta]", "]"}], "\[Rule]", "0"}], 
        ",", 
        RowBox[{
         RowBox[{"KroneckerDelta", "[", "\[Epsilon]", "]"}], "\[Rule]", 
         "0"}]}], "}"}]}], ")"}], "-", 
    SubscriptBox["raman", "B1g"]}], "]"}]}]], "Input",ExpressionUUID->\
"4807b7f1-96eb-4bf2-9615-0e2d2d289b5f"],

Cell[BoxData["0"], "Output",ExpressionUUID->"0c5bfebf-6c4f-47c8-b175-cb0c2ec6727a"]
}, Open  ]]
}, Open  ]]
}, Open  ]]
},
WindowSize->{1579, 894},
WindowMargins->{{Automatic, 248}, {Automatic, 79}},
FrontEndVersion->"11.1 for Microsoft Windows (64-bit) (April 18, 2017)",
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
Cell[1486, 35, 174, 3, 43, "Subsection", "ExpressionUUID" -> \
"04451415-e7ea-42cb-b75b-d7fcf4206ed1"],
Cell[CellGroupData[{
Cell[1685, 42, 99, 0, 34, "Subsubsection", "ExpressionUUID" -> \
"5dcc4ab5-6a1d-4e92-a2a7-0f59b76cf7e5"],
Cell[1787, 44, 266, 7, 30, "Input", "ExpressionUUID" -> \
"44b438e9-d7d0-46f1-a26a-0f094a593ed0"],
Cell[2056, 53, 550, 15, 30, "Input", "ExpressionUUID" -> \
"32c5a07f-a275-45ce-ae5d-f6eaf5d04b98"],
Cell[2609, 70, 504, 15, 30, "Input", "ExpressionUUID" -> \
"8c16f1ce-61f3-46c8-93a9-c21c9861d11c"],
Cell[3116, 87, 240, 6, 30, "Input", "ExpressionUUID" -> \
"ab9cc2de-a827-45ec-80cb-d921aba76197"],
Cell[3359, 95, 907, 26, 55, "Input", "ExpressionUUID" -> \
"4e8b3e40-cbcf-4172-a6c4-0429a6b91e2a"],
Cell[4269, 123, 793, 24, 33, "Input", "ExpressionUUID" -> \
"f2c886cd-b96f-42b7-9563-9c0d8b657200"],
Cell[5065, 149, 252, 5, 30, "Input", "ExpressionUUID" -> \
"839e042b-6efe-411d-9900-fc2d74ebbaaa"],
Cell[5320, 156, 1675, 52, 55, "Input", "ExpressionUUID" -> \
"32bb2983-417a-41c4-86b8-0da088df5fd4"],
Cell[6998, 210, 1675, 52, 55, "Input", "ExpressionUUID" -> \
"49e89b23-d9f3-47d1-922b-f0056862b500"],
Cell[8676, 264, 1573, 47, 55, "Input", "ExpressionUUID" -> \
"81e44988-eec1-41f7-8e84-bf43e7153d78"]
}, Open  ]],
Cell[CellGroupData[{
Cell[10286, 316, 123, 0, 34, "Subsubsection", "ExpressionUUID" -> \
"2bc99d75-87bf-4dc0-92e1-936596065360"],
Cell[CellGroupData[{
Cell[10434, 320, 405, 11, 52, "Input", "ExpressionUUID" -> \
"d93ee796-f4bb-4e00-b8e0-1ed35e963262"],
Cell[10842, 333, 198, 5, 33, "Output", "ExpressionUUID" -> \
"17e97fc6-a0e2-4e6b-b970-fdce8498b8df"]
}, Open  ]],
Cell[CellGroupData[{
Cell[11077, 343, 471, 13, 52, "Input", "ExpressionUUID" -> \
"db3e7010-bcb1-49c6-ac24-2d09c0cf0926"],
Cell[11551, 358, 303, 9, 33, "Output", "ExpressionUUID" -> \
"fbbccbfc-e479-4763-b03e-54d673cb9ff7"]
}, Open  ]],
Cell[11869, 370, 1391, 44, 33, "Input", "ExpressionUUID" -> \
"3c77d700-7922-4693-b81a-894372d33ec0"],
Cell[CellGroupData[{
Cell[13285, 418, 1767, 53, 52, "Input", "ExpressionUUID" -> \
"b59d2b5c-6b6c-4d1d-ac06-3adfb6add2f2"],
Cell[15055, 473, 83, 0, 30, "Output", "ExpressionUUID" -> \
"9464890a-5ac0-4120-8db5-9248dd0dd1fd"]
}, Open  ]],
Cell[CellGroupData[{
Cell[15175, 478, 1767, 53, 52, "Input", "ExpressionUUID" -> \
"da735f2c-61e9-4d51-8120-6c49737fb482"],
Cell[16945, 533, 83, 0, 30, "Output", "ExpressionUUID" -> \
"a5242d33-11d0-401a-9dd8-ba226de14a07"]
}, Open  ]],
Cell[CellGroupData[{
Cell[17065, 538, 1586, 47, 52, "Input", "ExpressionUUID" -> \
"43f954af-94a5-4d11-bfa8-7c93a41c4461"],
Cell[18654, 587, 666, 21, 33, "Output", "ExpressionUUID" -> \
"e20bea6a-5b56-4e97-a7fb-c6058909be6d"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[19369, 614, 114, 0, 34, "Subsubsection", "ExpressionUUID" -> \
"d74a3517-f66a-45ff-adb4-c09bccbbfef3"],
Cell[19486, 616, 1604, 53, 55, "Input", "ExpressionUUID" -> \
"3ee05fa0-4305-4fc5-931e-821c80e3f187"],
Cell[CellGroupData[{
Cell[21115, 673, 1274, 39, 52, "Input", "ExpressionUUID" -> \
"c70e2f82-1fbb-4cce-9135-8c74f8cc0564"],
Cell[22392, 714, 83, 0, 30, "Output", "ExpressionUUID" -> \
"5648fb8c-6dac-4027-a29a-86f5dff0fe98"]
}, Open  ]],
Cell[CellGroupData[{
Cell[22512, 719, 1274, 39, 52, "Input", "ExpressionUUID" -> \
"2399a700-f5ad-4d65-9feb-29a87d4f00ec"],
Cell[23789, 760, 83, 0, 30, "Output", "ExpressionUUID" -> \
"43ca4a33-be5a-4562-a6d5-9dc3a7cf9bfc"]
}, Open  ]],
Cell[CellGroupData[{
Cell[23909, 765, 1010, 31, 52, "Input", "ExpressionUUID" -> \
"50008626-57b9-4111-81cb-17670cdd7b38"],
Cell[24922, 798, 218, 6, 33, "Output", "ExpressionUUID" -> \
"4a25734f-b56c-4176-ba9f-9a229f6fc457"]
}, Open  ]],
Cell[CellGroupData[{
Cell[25177, 809, 1274, 39, 52, "Input", "ExpressionUUID" -> \
"6b1ef567-a761-442d-9828-9e5d8c0f76b6"],
Cell[26454, 850, 83, 0, 30, "Output", "ExpressionUUID" -> \
"f2890f77-661b-46a6-bdaf-b0b444a9f68d"]
}, Open  ]],
Cell[CellGroupData[{
Cell[26574, 855, 1274, 39, 52, "Input", "ExpressionUUID" -> \
"5db9708f-4837-480a-86ba-485eca0efbef"],
Cell[27851, 896, 83, 0, 30, "Output", "ExpressionUUID" -> \
"4c3ce712-410a-4360-82d4-88eaf248bdaa"]
}, Open  ]],
Cell[CellGroupData[{
Cell[27971, 901, 1010, 31, 52, "Input", "ExpressionUUID" -> \
"481af523-bce2-4b53-9ad3-686258ce9e29"],
Cell[28984, 934, 218, 6, 33, "Output", "ExpressionUUID" -> \
"1559858c-45a9-4d47-b55d-7e92e1b4cbec"]
}, Open  ]],
Cell[29217, 943, 1093, 34, 52, "Input", "ExpressionUUID" -> \
"b48889d7-ae23-45aa-b3fc-fd0495a1f5c5"],
Cell[CellGroupData[{
Cell[30335, 981, 1474, 44, 52, "Input", "ExpressionUUID" -> \
"192838cf-098e-4a98-9a42-41716ed15600"],
Cell[31812, 1027, 83, 0, 30, "Output", "ExpressionUUID" -> \
"dd7e5078-2863-4518-83ea-5299689dc2b8"]
}, Open  ]],
Cell[CellGroupData[{
Cell[31932, 1032, 524, 14, 51, "Input", "ExpressionUUID" -> \
"35174cfd-8dd2-44fd-9cc1-067797823c87"],
Cell[32459, 1048, 83, 0, 30, "Output", "ExpressionUUID" -> \
"180f4cb4-24cf-4e69-bcd6-5fd05f2f7e12"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[32591, 1054, 108, 0, 34, "Subsubsection", "ExpressionUUID" -> \
"7409eea4-b238-40c2-9be8-e0319f0504ce"],
Cell[32702, 1056, 1394, 45, 34, "Input", "ExpressionUUID" -> \
"15327b87-0323-4dae-9ff9-59fe45e1038b"],
Cell[CellGroupData[{
Cell[34121, 1105, 1791, 54, 54, "Input", "ExpressionUUID" -> \
"2dee7dab-31d6-4033-9f82-8571f1e1f9b2"],
Cell[35915, 1161, 83, 0, 30, "Output", "ExpressionUUID" -> \
"5985ad9f-b7b3-42db-90f4-8672a6a0427d"]
}, Open  ]],
Cell[CellGroupData[{
Cell[36035, 1166, 1791, 54, 54, "Input", "ExpressionUUID" -> \
"e9b5e707-25dd-493e-801a-14b416c7bdd6"],
Cell[37829, 1222, 83, 0, 30, "Output", "ExpressionUUID" -> \
"3f0c8c17-e3f7-4ced-8c9f-f1e006a6a834"]
}, Open  ]],
Cell[CellGroupData[{
Cell[37949, 1227, 1586, 47, 52, "Input", "ExpressionUUID" -> \
"fda87cf2-4d3f-47f8-85ae-946b1c1d5553"],
Cell[39538, 1276, 602, 19, 33, "Output", "ExpressionUUID" -> \
"7ce4942e-2c31-4f50-ab6b-f9fefcfbef13"]
}, Open  ]],
Cell[40155, 1298, 759, 25, 34, "Input", "ExpressionUUID" -> \
"d0938527-7f5e-470f-96da-4b2162531716"],
Cell[CellGroupData[{
Cell[40939, 1327, 2284, 68, 54, "Input", "ExpressionUUID" -> \
"791e684e-2697-41d9-ab96-127eccd848d6"],
Cell[43226, 1397, 83, 0, 30, "Output", "ExpressionUUID" -> \
"f0ab60d8-4b69-4f75-9140-3cc0e0312589"]
}, Open  ]],
Cell[CellGroupData[{
Cell[43346, 1402, 2284, 68, 54, "Input", "ExpressionUUID" -> \
"5fb9c248-d3c9-4f31-a187-992d5741ce83"],
Cell[45633, 1472, 83, 0, 30, "Output", "ExpressionUUID" -> \
"2b2688df-2c87-4ab1-b5b8-6f07467894c3"]
}, Open  ]],
Cell[CellGroupData[{
Cell[45753, 1477, 2064, 61, 52, "Input", "ExpressionUUID" -> \
"ad3a100b-01b0-415c-b08a-3528accb96d6"],
Cell[47820, 1540, 602, 19, 33, "Output", "ExpressionUUID" -> \
"cf87e8e9-c082-4fa1-b161-7eeb4db3eac9"]
}, Open  ]],
Cell[CellGroupData[{
Cell[48459, 1564, 2440, 71, 54, "Input", "ExpressionUUID" -> \
"fbcc240c-d686-4565-b66d-606779d11e87"],
Cell[50902, 1637, 83, 0, 30, "Output", "ExpressionUUID" -> \
"cfecf615-62f4-4d36-a62a-c4e176ea7d3f"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[51034, 1643, 112, 0, 34, "Subsubsection", "ExpressionUUID" -> \
"b2761392-2125-4a92-968f-1ad236858cbc"],
Cell[51149, 1645, 551, 17, 52, "Input", "ExpressionUUID" -> \
"890beb2a-0213-4cd3-a0f8-cd93143911d2"],
Cell[CellGroupData[{
Cell[51725, 1666, 1044, 31, 52, "Input", "ExpressionUUID" -> \
"1bf83b97-5b9f-48c5-999a-f2a36d30193f"],
Cell[52772, 1699, 83, 0, 30, "Output", "ExpressionUUID" -> \
"55302f1b-3d56-487e-858a-d86fa024568f"]
}, Open  ]],
Cell[52870, 1702, 664, 20, 74, "Input", "ExpressionUUID" -> \
"658ce7f9-745e-48a8-ad94-7bf5ac15e484"],
Cell[CellGroupData[{
Cell[53559, 1726, 1209, 36, 75, "Input", "ExpressionUUID" -> \
"348c8a86-0fc5-45c0-99c6-4afa87a8169f"],
Cell[54771, 1764, 83, 0, 30, "Output", "ExpressionUUID" -> \
"cd22a849-6aad-4d9d-968f-10e7db1beb42"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[54903, 1770, 100, 0, 34, "Subsubsection", "ExpressionUUID" -> \
"9b971066-7c7c-4378-b39e-f08f836f4af5"],
Cell[55006, 1772, 1470, 43, 33, "Input", "ExpressionUUID" -> \
"68783b49-a8ab-4e30-a244-128a73712eed"],
Cell[56479, 1817, 2703, 77, 33, "Input", "ExpressionUUID" -> \
"5967b3ee-c567-4558-9a6c-22f63c7ade8f"],
Cell[59185, 1896, 289, 9, 33, "Input", "ExpressionUUID" -> \
"2c309710-7858-48d6-8c59-bc23d96a63ab"],
Cell[CellGroupData[{
Cell[59499, 1909, 1808, 50, 52, "Input", "ExpressionUUID" -> \
"87a470ef-ae01-4cae-b137-265b7004ec15"],
Cell[61310, 1961, 83, 0, 30, "Output", "ExpressionUUID" -> \
"c5987598-0702-4c45-b478-4f9dfa0168df"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[61442, 1967, 94, 0, 34, "Subsubsection", "ExpressionUUID" -> \
"5e41d14f-6a91-40bd-823e-9bc8988d4949"],
Cell[61539, 1969, 2600, 74, 77, "Input", "ExpressionUUID" -> \
"5d64089e-66aa-43ff-83b3-d72b319a8676"],
Cell[CellGroupData[{
Cell[64164, 2047, 2163, 64, 106, "Input", "ExpressionUUID" -> \
"4807b7f1-96eb-4bf2-9615-0e2d2d289b5f"],
Cell[66330, 2113, 83, 0, 30, "Output", "ExpressionUUID" -> \
"0c5bfebf-6c4f-47c8-b175-cb0c2ec6727a"]
}, Open  ]]
}, Open  ]]
}, Open  ]]
}
]
*)

(* NotebookSignature 3v0W4xlzi#J5bAgW4RAd@DLJ *)
