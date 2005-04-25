% Minimum LaTeX version needed
\NeedsTeXFormat{LaTeX2e} 

% Document class provided
\ProvidesClass{jkthesis}

% Rquire package "ifthen"
\RequirePackage{ifthen}

% Declare boolean for language
\newboolean{@german} 
\setboolean{@german}{false}

% Declare language option
\DeclareOption{german}{\setboolean{@german}{true}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}} 
\ProcessOptions 

% Load document class "article"
\LoadClass[12pt,a4paper,pdftex,reqno,twoside]{article}

% Load package "hyperref"
\RequirePackage[a4paper=true, breaklinks=true, pdftex=true, hyperindex=true, colorlinks=true, 
                linkcolor=blue, anchorcolor=blue, citecolor=blue, filecolor=blue, menucolor=blue, 
                pagecolor=blue, urlcolor=blue, bookmarks=true, bookmarksopen=true, bookmarksnumbered=true, 
                pdfstartview=FitH]{hyperref}
\hypersetup
  {
    pdftitle = {}
    pdfsubject = {}
    pdfauthor = {\textcopyright\ Jens Keiner}
    pdfkeywords = {}
    pdfcreator = {}
    pdfproducer = {}
}

% Load some other packages
\RequirePackage{amsmath}
% \RequirePackage[backrefs]{amsrefs}
\RequirePackage{amsthm}
\RequirePackage{amssymb}
\RequirePackage{mathrsfs}
\RequirePackage{graphicx}
\RequirePackage{epsfig}
%\RequirePackage[amsmath,amsthm]{ntheorem}
\RequirePackage{t1enc}
\RequirePackage[latin1]{inputenc}
\RequirePackage[T1]{fontenc}
\RequirePackage{ae}
\RequirePackage{exscale}
%\RequirePackage{jkmath}

% Page setup

% Vertical
\setlength{\voffset}{0mm}
\setlength{\topmargin}{-15mm}
\setlength{\headheight}{10mm} 
\setlength{\headsep}{10mm}
\setlength{\textheight}{237mm}
\setlength{\footskip}{10mm}

% Horizontal
\setlength{\hoffset}{0mm}
\setlength{\oddsidemargin}{5mm}
\setlength{\evensidemargin}{5mm}
\setlength{\textwidth}{150mm}

% Set pagestyle
\pagestyle{headings}

% Define mathematical environments in English
\newcommand{\englishnames}
{
  \newtheorem{theorem}{Theorem}[section]
  \newtheorem{corollary}[theorem]{Corollary}
  \newtheorem{lemma}[theorem]{Lemma}
  \newtheorem{proposition}[theorem]{Proposition}
  \theoremstyle{definition}
  \newtheorem{definition}[theorem]{Definition}
  \newtheorem{conclusion}[theorem]{Conclusion}
  \newtheorem{example}[theorem]{Example}
  \newtheorem{Proof}[theorem]{Proof}
  \theoremstyle{remark}
  \newtheorem{remark}[theorem]{Remark}
  \numberwithin{equation}{section}
  %\newenvironment{proof}{\begin{Proof}}{\qed\end{Proof}}  
}

% Define mathematical environments in German
\newcommand{\germannames}
{
  \renewtheorem{theorem}{Satz}[section]
  \renewtheorem{corollary}[theorem]{Korollar}
  \renewtheorem{lemma}[theorem]{Lemma}
  \renewtheorem{proposition}[theorem]{Proposition}
  \theoremstyle{definition}
  \renewtheorem{definition}[theorem]{Definition}
  \newtheorem{conclusion}[theorem]{Folgerung}
  \renewtheorem{example}[theorem]{Beispiel}
  \theoremstyle{remark}
  \renewtheorem{remark}[theorem]{Bemerkung}
  \numberwithin{equation}{section}
  %\renewenvironment{proof}{\begin{Proof}}{\qed\end{Proof}}  
}

% Switch by language.
\ifthenelse{\boolean{@german}}
{\RequirePackage{ngerman}\germannames}
{\englishnames}

\newcommand{\jkNorm}[1]{\ensuremath{\left\|#1\right\|}}
\newcommand{\jkAbs}[1]{\ensuremath{\left\vert#1\right\vert}}
\newcommand{\jkSet}[1]{\ensuremath{\left\{#1\right\rbrace}}
\newcommand{\jkSetV}[3]{\ensuremath{\left\{#1\ \left#2\ #3\right.\right\rbrace}}
\newcommand{\jkEps}{\ensuremath{\varepsilon}}
\newcommand{\jkTo}{\ensuremath{\longrightarrow}}
%\newcommand\qed{\hspace*{\fill}$\Box$}
\newcommand{\N}{\ensuremath{\mathbb{N}}}
\newcommand{\NZ}{\ensuremath{\mathbb{N}_{0}}}
\newcommand{\Z}{\ensuremath{\mathbb{Z}}}
\newcommand{\Q}{\ensuremath{\mathbb{Q}}}
\newcommand{\R}{\ensuremath{\mathbb{R}}}
\newcommand{\C}{\ensuremath{\mathbb{C}}}
\newcommand{\Rp}{\ensuremath{\mathbb{R}^{+}}}
\newcommand{\Rn}{\ensuremath{\mathbb{R}^n}}
\newcommand{\Rnn}{\ensuremath{\mathbb{R}^{n \times n}}}
\newcommand{\jkS}{\ensuremath{\mathbb{S}}}
\newcommand{\jkSS}{\ensuremath{\mathbb{S}^2}}
\newcommand{\jkL}[1][\R]{\ensuremath{\mathbb{L}^1\left(#1\right)}}
\newcommand{\jkLL}[1][\R]{\ensuremath{\mathbb{L}^2\left(#1\right)}}
\newcommand{\jkInterv}[4]{\ensuremath{\left#1\left.#2,#3\right#4\right.}}
\newcommand{\jkFun}[2]{\ensuremath{#1\left(#2\right)}}
\newcommand{\jkPar}[1]{\ensuremath{\left(#1\right)}}
\newcommand{\jkEncl}[3]{\ensuremath{\left#1#2\right#3}}
\newcommand{\jkO}[1]{\ensuremath{\mathcal{O}\jkPar{#1}}}
\newcommand{\jko}[1]{\ensuremath{o\jkPar{#1}}}
\newcommand{\jkTheta}[1]{\ensuremath{\Theta\jkPar{#1}}}
\newcommand{\jkSX}{\ensuremath{\xi}}
\newcommand{\jkSY}{\ensuremath{\eta}}
\newcommand{\jkSpacer}{\ensuremath{{}^{}}}
\newcommand{\jkScalProd}[2]{\ensuremath{\left<#1,#2\right>}}



