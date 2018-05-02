# MetabolomicWorkflows
---
title: <b> Metabolomic Workflows </b>
subtitle: "Biological Analysis"
author: "Pol Castellano Escuder"
date: "May 2018"
output:
  html_document:
---

# Whole analysis softwares (top-down) 

## Metaboanalyst 

MetaboAnalyst is a web server designed to permit comprehensive metabolomic data analysis, visualization and interpretation [@Xia2015]. It supports a wide range of complex statistical calculations and high quality graphical rendering functions that require significant computational resources.

#### <u> Upload dataset </u>

MetaboAnalyst accepts data from either targeted profiling (concentration tables) or metabolic fingerprinting approaches (spectral bins, peak lists) produced from either NMR, LC-MS, or GC-MS. 

GC/LC-MS spectra saved as open data format (NetCDF, mzDATA, mzXML) can also be processed using the XCMS package [@Tautenhahn2012]. 

#### <u> Processing </u>

**1. Missing values:**

 + Remove features with too many missing values:
 
 By default, MetaboAnalyst remove all features in dataset that have > 50% of missing values. This parameter can be modified by user.

 + Impute values:
 
 MetaboAnalyst offers a variety of methods to deal with missing values [@Armitage2015]. By default, the missing values are treated as the result of **low signal intensity**. 
