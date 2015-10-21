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
NotebookDataLength[      4543,        171]
NotebookOptionsPosition[      4806,        157]
NotebookOutlinePosition[      5150,        172]
CellTagsIndexPosition[      5107,        169]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell[BoxData[
 RowBox[{"TrigToExp", "[", 
  RowBox[{"ArcCosh", "[", "x", "]"}], "]"}]], "Input"],

Cell[BoxData[
 RowBox[{"Log", "[", 
  RowBox[{"x", "+", 
   RowBox[{
    SqrtBox[
     RowBox[{
      RowBox[{"-", "1"}], "+", "x"}]], " ", 
    SqrtBox[
     RowBox[{"1", "+", "x"}]]}]}], "]"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"FullSimplify", "[", 
  RowBox[{
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"Exp", "[", 
       RowBox[{"2", 
        RowBox[{"ArcCosh", "[", "x", "]"}]}], "]"}], "-", "1"}], ")"}], "-", 
    RowBox[{"2", " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{
        RowBox[{"(", 
         RowBox[{"x", "-", "1"}], ")"}], 
        RowBox[{"(", 
         RowBox[{"x", "+", "1"}], ")"}]}], "+", 
       RowBox[{"x", 
        SqrtBox[
         RowBox[{
          RowBox[{"(", 
           RowBox[{"x", "-", "1"}], ")"}], 
          RowBox[{"(", 
           RowBox[{"x", "+", "1"}], ")"}]}]]}]}], ")"}]}]}], ",", 
   RowBox[{"Assumptions", "\[Rule]", 
    RowBox[{"{", 
     RowBox[{"x", ">", "1"}], "}"}]}]}], "]"}]], "Input"],

Cell[BoxData["0"], "Output"]
}, Open  ]],

Cell[BoxData[{
 RowBox[{
  RowBox[{
   SubscriptBox["U", "val"], "=", "4"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   SubscriptBox["\[Tau]", "val"], "=", 
   RowBox[{"1", "/", "11"}]}], ";"}]}], "Input"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  SubscriptBox["\[Lambda]", "val"], "=", 
  RowBox[{"ArcCosh", "[", 
   RowBox[{"Exp", "[", 
    RowBox[{
     SubscriptBox["U", "val"], 
     RowBox[{
      SubscriptBox["\[Tau]", "val"], "/", "2"}]}], "]"}], 
   "]"}]}], "\[IndentingNewLine]", 
 RowBox[{"N", "[", "%", "]"}]}], "Input"],

Cell[BoxData[
 RowBox[{"ArcCosh", "[", 
  SuperscriptBox["\[ExponentialE]", 
   RowBox[{"2", "/", "11"}]], "]"}]], "Output"],

Cell[BoxData["0.6214513424467174`"], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{"save", " ", "reference", " ", "values", " ", "to", " ", "disk"}], 
   " ", "*)"}], "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{"Block", "[", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{"\[Lambda]", "=", 
        SubscriptBox["\[Lambda]", "val"]}], "}"}], ",", 
      RowBox[{"N", "[", 
       RowBox[{
        RowBox[{"{", 
         RowBox[{
          RowBox[{"Exp", "[", 
           RowBox[{"-", "\[Lambda]"}], "]"}], ",", 
          RowBox[{"Exp", "[", "\[Lambda]", "]"}], ",", 
          RowBox[{"Exp", "[", "\[Lambda]", "]"}], ",", 
          RowBox[{"Exp", "[", 
           RowBox[{"-", "\[Lambda]"}], "]"}], ",", 
          RowBox[{
           RowBox[{"Exp", "[", 
            RowBox[{"2", "\[Lambda]"}], "]"}], "-", "1"}], ",", 
          RowBox[{
           RowBox[{"Exp", "[", 
            RowBox[{
             RowBox[{"-", "2"}], "\[Lambda]"}], "]"}], "-", "1"}]}], "}"}], 
        ",", 
        RowBox[{"2", "MachinePrecision"}]}], "]"}]}], "]"}], ";"}], 
   "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"Export", "[", 
     RowBox[{
      RowBox[{
       RowBox[{"NotebookDirectory", "[", "]"}], "<>", 
       RowBox[{"FileBaseName", "[", 
        RowBox[{"NotebookFileName", "[", "]"}], "]"}], "<>", "\"\<.dat\>\""}],
       ",", 
      RowBox[{"N", "[", "%", "]"}], ",", "\"\<Real64\>\""}], "]"}], 
    ";"}]}]}]], "Input"]
},
WindowSize->{1199, 583},
WindowMargins->{{196, Automatic}, {Automatic, 183}},
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
Cell[1486, 35, 96, 2, 31, "Input"],
Cell[1585, 39, 206, 8, 36, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[1828, 52, 762, 26, 39, "Input"],
Cell[2593, 80, 28, 0, 31, "Output"]
}, Open  ]],
Cell[2636, 83, 217, 7, 52, "Input"],
Cell[CellGroupData[{
Cell[2878, 94, 313, 10, 52, "Input"],
Cell[3194, 106, 124, 3, 33, "Output"],
Cell[3321, 111, 46, 0, 31, "Output"]
}, Open  ]],
Cell[3382, 114, 1420, 41, 72, "Input"]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature Twp0db5c0jf8RBK#iHdURVJv *)
