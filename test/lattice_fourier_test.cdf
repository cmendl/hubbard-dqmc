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
NotebookDataLength[     20608,        610]
NotebookOptionsPosition[     20247,        575]
NotebookOutlinePosition[     20591,        590]
CellTagsIndexPosition[     20548,        587]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["Dimension 4\[Times]6", "Subsection"],

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
   RowBox[{"arbitrary", " ", "\"\<lattice vector\>\""}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["v", "val"], "=", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{
         FractionBox["7", "10"], ",", 
         RowBox[{"-", 
          FractionBox["9", "10"]}], ",", 
         RowBox[{"-", 
          FractionBox["3", "10"]}], ",", 
         RowBox[{"-", 
          FractionBox["1", "10"]}], ",", 
         RowBox[{"-", "1"}], ",", 
         FractionBox["3", "10"]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"-", 
          FractionBox["1", "5"]}], ",", "1", ",", 
         RowBox[{"-", 
          FractionBox["4", "5"]}], ",", 
         RowBox[{"-", 
          FractionBox["4", "5"]}], ",", 
         FractionBox["2", "5"], ",", 
         RowBox[{"-", 
          FractionBox["4", "5"]}]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         FractionBox["1", "5"], ",", 
         FractionBox["3", "10"], ",", 
         FractionBox["1", "2"], ",", 
         RowBox[{"-", 
          FractionBox["2", "5"]}], ",", 
         RowBox[{"-", 
          FractionBox["4", "5"]}], ",", 
         FractionBox["9", "10"]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         FractionBox["1", "10"], ",", 
         FractionBox["3", "10"], ",", 
         RowBox[{"-", 
          FractionBox["3", "5"]}], ",", 
         RowBox[{"-", 
          FractionBox["3", "10"]}], ",", 
         RowBox[{"-", 
          FractionBox["4", "5"]}], ",", 
         FractionBox["3", "10"]}], "}"}]}], "}"}]}], ";"}], 
   "\[IndentingNewLine]", 
   RowBox[{"Dimensions", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"4", ",", "6"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{
    "Fourier", " ", "transform", " ", "without", " ", "normalization"}], ";", 
    " ", 
    RowBox[{
    "conventions", " ", "of", " ", "forward", " ", "and", " ", "backward", 
     " ", "Fourier", " ", "transforms", " ", "differ"}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   SubscriptBox["v", 
    RowBox[{"F", ",", "val"}]], "=", 
   RowBox[{
    SqrtBox[
     RowBox[{
      SubscriptBox["n", "x"], 
      SubscriptBox["n", "y"]}]], 
    RowBox[{"InverseFourier", "[", 
     SubscriptBox["v", "val"], "]"}]}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{
     RowBox[{
      RowBox[{"-", "2.8000000000000003`"}], "+", 
      RowBox[{"0.`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"4.8`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.8660254037844389`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"0.20000000000000007`", "\[VeryThinSpace]", "+", 
      RowBox[{"0.8660254037844384`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{
      RowBox[{"-", "2.4`"}], "+", 
      RowBox[{"0.`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"0.20000000000000007`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.8660254037844384`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"4.8`", "\[VeryThinSpace]", "+", 
      RowBox[{"0.8660254037844389`", " ", "\[ImaginaryI]"}]}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{
      RowBox[{"-", "1.9999999999999998`"}], "+", 
      RowBox[{"0.2000000000000002`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{
      RowBox[{"-", "0.5464101615137754`"}], "+", 
      RowBox[{"1.5392304845413263`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{
      RowBox[{"-", "0.571281292110204`"}], "+", 
      RowBox[{"1.0999999999999999`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"0.9999999999999999`", "\[VeryThinSpace]", "-", 
      RowBox[{"1.6000000000000003`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"4.971281292110204`", "\[VeryThinSpace]", "+", 
      RowBox[{"1.0999999999999999`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"0.14641016151377506`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.5392304845413267`", " ", "\[ImaginaryI]"}]}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"1.6000000000000005`", "\[VeryThinSpace]", "+", 
      RowBox[{"0.`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"0.2000000000000002`", "\[VeryThinSpace]", "+", 
      RowBox[{"0.5196152422706632`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"1.5999999999999996`", "\[VeryThinSpace]", "+", 
      RowBox[{"5.715767664977295`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"0.8000000000000003`", "\[VeryThinSpace]", "+", 
      RowBox[{"0.`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"1.5999999999999996`", "\[VeryThinSpace]", "-", 
      RowBox[{"5.715767664977295`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"0.2000000000000002`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.5196152422706632`", " ", "\[ImaginaryI]"}]}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{
      RowBox[{"-", "1.9999999999999998`"}], "-", 
      RowBox[{"0.2000000000000002`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"0.14641016151377506`", "\[VeryThinSpace]", "+", 
      RowBox[{"0.5392304845413267`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"4.971281292110204`", "\[VeryThinSpace]", "-", 
      RowBox[{"1.0999999999999999`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"0.9999999999999999`", "\[VeryThinSpace]", "+", 
      RowBox[{"1.6000000000000003`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{
      RowBox[{"-", "0.571281292110204`"}], "-", 
      RowBox[{"1.0999999999999999`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{
      RowBox[{"-", "0.5464101615137754`"}], "-", 
      RowBox[{"1.5392304845413263`", " ", "\[ImaginaryI]"}]}]}], "}"}]}], 
  "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
   "complex", " ", "conjugation", " ", "effectively", " ", "flips", " ", 
    "signs", " ", "of", " ", "indices"}], "*)"}], "\[IndentingNewLine]", 
  RowBox[{"Norm", "[", 
   RowBox[{
    RowBox[{"Conjugate", "[", 
     SubscriptBox["v", 
      RowBox[{"F", ",", "val"}]], "]"}], "-", 
    RowBox[{"Table", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["v", 
        RowBox[{"F", ",", "val"}]], "\[LeftDoubleBracket]", 
       RowBox[{
        RowBox[{
         RowBox[{"Mod", "[", 
          RowBox[{
           RowBox[{"-", "i"}], ",", 
           SubscriptBox["n", "x"]}], "]"}], "+", "1"}], ",", 
        RowBox[{
         RowBox[{"Mod", "[", 
          RowBox[{
           RowBox[{"-", "j"}], ",", 
           SubscriptBox["n", "y"]}], "]"}], "+", "1"}]}], 
       "\[RightDoubleBracket]"}], ",", 
      RowBox[{"{", 
       RowBox[{"i", ",", "0", ",", 
        RowBox[{
         SubscriptBox["n", "x"], "-", "1"}]}], "}"}], ",", 
      RowBox[{"{", 
       RowBox[{"j", ",", "0", ",", 
        RowBox[{
         SubscriptBox["n", "y"], "-", "1"}]}], "}"}]}], "]"}]}], 
   "]"}]}]], "Input"],

Cell[BoxData["0.`"], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["fn", "base"], "=", 
   RowBox[{
    RowBox[{"NotebookDirectory", "[", "]"}], "<>", 
    RowBox[{"FileBaseName", "[", 
     RowBox[{"NotebookFileName", "[", "]"}], "]"}], "<>", "\"\<_\>\"", "<>", 
    RowBox[{"ToString", "[", 
     SubscriptBox["n", "x"], "]"}], "<>", "\"\<x\>\"", "<>", 
    RowBox[{"ToString", "[", 
     SubscriptBox["n", "y"], "]"}]}]}], ";"}]], "Input"],

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
       SubscriptBox["fn", "base"], "<>", "\"\<_v.dat\>\""}], ",", 
      RowBox[{"N", "[", 
       RowBox[{"Flatten", "[", 
        RowBox[{"Transpose", "[", 
         SubscriptBox["v", "val"], "]"}], "]"}], "]"}], ",", 
      "\"\<Real64\>\""}], "]"}], ";"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"Export", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["fn", "base"], "<>", "\"\<_vF.dat\>\""}], ",", 
      RowBox[{
       RowBox[{
        RowBox[{"{", 
         RowBox[{
          RowBox[{"Re", "[", "#", "]"}], ",", 
          RowBox[{"Im", "[", "#", "]"}]}], "}"}], "&"}], "/@", 
       RowBox[{"Flatten", "[", 
        RowBox[{"Transpose", "[", 
         SubscriptBox["v", 
          RowBox[{"F", ",", "val"}]], "]"}], "]"}]}], ",", "\"\<Real64\>\""}],
      "]"}], ";"}]}]}]], "Input"]
}, Open  ]],

Cell[CellGroupData[{

Cell["Dimension 5\[Times]6", "Subsection"],

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
   RowBox[{"arbitrary", " ", "\"\<lattice vector\>\""}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{
     SubscriptBox["v", "val"], "=", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{
         FractionBox["7", "10"], ",", 
         RowBox[{"-", 
          FractionBox["3", "10"]}], ",", 
         FractionBox["2", "5"], ",", "0", ",", 
         RowBox[{"-", 
          FractionBox["1", "2"]}], ",", 
         RowBox[{"-", 
          FractionBox["3", "10"]}]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"-", 
          FractionBox["1", "5"]}], ",", 
         RowBox[{"-", "1"}], ",", 
         FractionBox["1", "2"], ",", 
         FractionBox["4", "5"], ",", 
         FractionBox["2", "5"], ",", 
         RowBox[{"-", 
          FractionBox["1", "5"]}]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"-", 
          FractionBox["1", "10"]}], ",", "0", ",", "1", ",", 
         FractionBox["3", "5"], ",", 
         FractionBox["1", "5"], ",", 
         RowBox[{"-", 
          FractionBox["1", "10"]}]}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"-", 
          FractionBox["9", "10"]}], ",", 
         RowBox[{"-", 
          FractionBox["1", "2"]}], ",", 
         FractionBox["2", "5"], ",", "0", ",", 
         RowBox[{"-", 
          FractionBox["1", "5"]}], ",", "0"}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"-", 
          FractionBox["1", "5"]}], ",", 
         FractionBox["1", "10"], ",", 
         FractionBox["2", "5"], ",", 
         FractionBox["1", "5"], ",", 
         RowBox[{"-", 
          FractionBox["3", "10"]}], ",", 
         FractionBox["3", "5"]}], "}"}]}], "}"}]}], ";"}], 
   "\[IndentingNewLine]", 
   RowBox[{"Dimensions", "[", "%", "]"}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"5", ",", "6"}], "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
    RowBox[{
    "Fourier", " ", "transform", " ", "without", " ", "normalization"}], ";", 
    " ", 
    RowBox[{
    "conventions", " ", "of", " ", "forward", " ", "and", " ", "backward", 
     " ", "Fourier", " ", "transforms", " ", "differ"}]}], " ", "*)"}], 
  "\[IndentingNewLine]", 
  RowBox[{
   SubscriptBox["v", 
    RowBox[{"F", ",", "val"}]], "=", 
   RowBox[{
    SqrtBox[
     RowBox[{
      SubscriptBox["n", "x"], 
      SubscriptBox["n", "y"]}]], 
    RowBox[{"InverseFourier", "[", 
     SubscriptBox["v", "val"], "]"}]}]}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{
     RowBox[{"1.5000000000000002`", "\[VeryThinSpace]", "+", 
      RowBox[{"0.`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{
      RowBox[{"-", "4.3`"}], "-", 
      RowBox[{"1.2124355652982142`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"0.5999999999999998`", "\[VeryThinSpace]", "+", 
      RowBox[{"4.156921938165305`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"1.6999999999999997`", "\[VeryThinSpace]", "+", 
      RowBox[{"0.`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"0.5999999999999998`", "\[VeryThinSpace]", "-", 
      RowBox[{"4.156921938165305`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{
      RowBox[{"-", "4.3`"}], "+", 
      RowBox[{"1.2124355652982142`", " ", "\[ImaginaryI]"}]}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"0.016311896062463076`", "\[VeryThinSpace]", "-", 
      RowBox[{"1.170270448271348`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"2.2231032722963233`", "\[VeryThinSpace]", "+", 
      RowBox[{"1.9683534214857474`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"1.354663028285239`", "\[VeryThinSpace]", "-", 
      RowBox[{"1.4248797001452125`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"0.9072949016875156`", "\[VeryThinSpace]", "-", 
      RowBox[{"2.467446886053801`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"2.25607605877714`", "\[VeryThinSpace]", "-", 
      RowBox[{"1.585449640543103`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"1.5550120226411297`", "\[VeryThinSpace]", "+", 
      RowBox[{"1.8583240425238456`", " ", "\[ImaginaryI]"}]}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{
      RowBox[{"-", "0.7663118960624632`"}], "+", 
      RowBox[{"2.956850871772666`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"2.5029868997779534`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.34627847500602327`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"0.6963750097508049`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.2204022594963477`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"1.2427050983124845`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.4735038167780705`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"0.34288590318681555`", "\[VeryThinSpace]", "+", 
      RowBox[{"0.1999753020368742`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"0.2688978052845936`", "\[VeryThinSpace]", "+", 
      RowBox[{"2.4484296556876375`", " ", "\[ImaginaryI]"}]}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{
      RowBox[{"-", "0.7663118960624632`"}], "-", 
      RowBox[{"2.956850871772666`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"0.2688978052845936`", "\[VeryThinSpace]", "-", 
      RowBox[{"2.4484296556876375`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"0.34288590318681555`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.1999753020368742`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"1.2427050983124845`", "\[VeryThinSpace]", "+", 
      RowBox[{"0.4735038167780705`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"0.6963750097508049`", "\[VeryThinSpace]", "+", 
      RowBox[{"0.2204022594963477`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"2.5029868997779534`", "\[VeryThinSpace]", "+", 
      RowBox[{"0.34627847500602327`", " ", "\[ImaginaryI]"}]}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"0.016311896062463076`", "\[VeryThinSpace]", "+", 
      RowBox[{"1.170270448271348`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"1.5550120226411297`", "\[VeryThinSpace]", "-", 
      RowBox[{"1.8583240425238456`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"2.25607605877714`", "\[VeryThinSpace]", "+", 
      RowBox[{"1.585449640543103`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"0.9072949016875156`", "\[VeryThinSpace]", "+", 
      RowBox[{"2.467446886053801`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"1.354663028285239`", "\[VeryThinSpace]", "+", 
      RowBox[{"1.4248797001452125`", " ", "\[ImaginaryI]"}]}], ",", 
     RowBox[{"2.2231032722963233`", "\[VeryThinSpace]", "-", 
      RowBox[{"1.9683534214857474`", " ", "\[ImaginaryI]"}]}]}], "}"}]}], 
  "}"}]], "Output"]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(*", " ", 
   RowBox[{
   "complex", " ", "conjugation", " ", "effectively", " ", "flips", " ", 
    "signs", " ", "of", " ", "indices"}], "*)"}], "\[IndentingNewLine]", 
  RowBox[{"Norm", "[", 
   RowBox[{
    RowBox[{"Conjugate", "[", 
     SubscriptBox["v", 
      RowBox[{"F", ",", "val"}]], "]"}], "-", 
    RowBox[{"Table", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["v", 
        RowBox[{"F", ",", "val"}]], "\[LeftDoubleBracket]", 
       RowBox[{
        RowBox[{
         RowBox[{"Mod", "[", 
          RowBox[{
           RowBox[{"-", "i"}], ",", 
           SubscriptBox["n", "x"]}], "]"}], "+", "1"}], ",", 
        RowBox[{
         RowBox[{"Mod", "[", 
          RowBox[{
           RowBox[{"-", "j"}], ",", 
           SubscriptBox["n", "y"]}], "]"}], "+", "1"}]}], 
       "\[RightDoubleBracket]"}], ",", 
      RowBox[{"{", 
       RowBox[{"i", ",", "0", ",", 
        RowBox[{
         SubscriptBox["n", "x"], "-", "1"}]}], "}"}], ",", 
      RowBox[{"{", 
       RowBox[{"j", ",", "0", ",", 
        RowBox[{
         SubscriptBox["n", "y"], "-", "1"}]}], "}"}]}], "]"}]}], 
   "]"}]}]], "Input"],

Cell[BoxData["0.`"], "Output"]
}, Open  ]],

Cell[BoxData[
 RowBox[{
  RowBox[{
   SubscriptBox["fn", "base"], "=", 
   RowBox[{
    RowBox[{"NotebookDirectory", "[", "]"}], "<>", 
    RowBox[{"FileBaseName", "[", 
     RowBox[{"NotebookFileName", "[", "]"}], "]"}], "<>", "\"\<_\>\"", "<>", 
    RowBox[{"ToString", "[", 
     SubscriptBox["n", "x"], "]"}], "<>", "\"\<x\>\"", "<>", 
    RowBox[{"ToString", "[", 
     SubscriptBox["n", "y"], "]"}]}]}], ";"}]], "Input"],

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
       SubscriptBox["fn", "base"], "<>", "\"\<_v.dat\>\""}], ",", 
      RowBox[{"N", "[", 
       RowBox[{"Flatten", "[", 
        RowBox[{"Transpose", "[", 
         SubscriptBox["v", "val"], "]"}], "]"}], "]"}], ",", 
      "\"\<Real64\>\""}], "]"}], ";"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"Export", "[", 
     RowBox[{
      RowBox[{
       SubscriptBox["fn", "base"], "<>", "\"\<_vF.dat\>\""}], ",", 
      RowBox[{
       RowBox[{
        RowBox[{"{", 
         RowBox[{
          RowBox[{"Re", "[", "#", "]"}], ",", 
          RowBox[{"Im", "[", "#", "]"}]}], "}"}], "&"}], "/@", 
       RowBox[{"Flatten", "[", 
        RowBox[{"Transpose", "[", 
         SubscriptBox["v", 
          RowBox[{"F", ",", "val"}]], "]"}], "]"}]}], ",", "\"\<Real64\>\""}],
      "]"}], ";"}]}]}]], "Input"]
}, Open  ]]
},
WindowSize->{1350, 770},
WindowMargins->{{384, Automatic}, {Automatic, 166}},
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
Cell[1486, 35, 42, 0, 43, "Subsection"],
Cell[1531, 37, 320, 10, 72, "Input"],
Cell[CellGroupData[{
Cell[1876, 51, 1737, 55, 88, "Input"],
Cell[3616, 108, 73, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[3726, 115, 612, 20, 61, "Input"],
Cell[4341, 137, 3261, 67, 92, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[7639, 209, 1163, 36, 52, "Input"],
Cell[8805, 247, 30, 0, 31, "Output"]
}, Open  ]],
Cell[8850, 250, 426, 11, 31, "Input"],
Cell[9279, 263, 1052, 31, 72, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[10368, 299, 42, 0, 43, "Subsection"],
Cell[10413, 301, 320, 10, 72, "Input"],
Cell[CellGroupData[{
Cell[10758, 315, 1897, 59, 88, "Input"],
Cell[12658, 376, 73, 2, 31, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[12768, 383, 612, 20, 61, "Input"],
Cell[13383, 405, 4119, 77, 112, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[17539, 487, 1163, 36, 52, "Input"],
Cell[18705, 525, 30, 0, 31, "Output"]
}, Open  ]],
Cell[18750, 528, 426, 11, 31, "Input"],
Cell[19179, 541, 1052, 31, 72, "Input"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature pupRMouvIFHgNB1#eWtOFhRF *)
