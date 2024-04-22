# Notes

## Title

**Seesaw vs. sawtooth??**

Seesaw year effects in index standardization: why they happen and how to avoid them

Data integration across scientific surveys: what causes oscillating year effects and how to avoid them

A cautionary tale about spatiotemporal index standardization with irregular sampling

## Journals

- CJFAS
- Fisheries Research
- ICES JMS
- Ecography?? too fishy?
- MEE??
- Eco. Appl.?

## Random uncategorized thoughts

- partial coverage (e.g., arrowtooth WCHG survey 2014) causes problems
- take a full survey, show how degrading causes an outlying year!
- show self simulation test --- spirit of Quang's inside quillback


## Introduction

- spatiotemporal modelling; standard factor/IID structure; widely used
- biennial, moving this way, moving towards stitching, budget cuts: dropping parts some years
- regression discontinuity design
- early signs this is a problem!
- here we explore when this happens, why this happens, what makes it better/worse, and how we can fix
- goal is to provide guidance on what to do

## Methods

## Results

## Discussion

- can happen to lesser degree with less extreme examples and be harder to detect

## Guidance

- Look for it (do we have guidance when it's not an obvious N/S thing? e.g., plot average latitude against average index?? use judgement: is there any systematic pattern with sampling?)
- Penalize time sufficiently
- Self- and cross-simulation test... can I recapture? have I penalized time too much or too little?
- simulate with new random fields?? and with full vs. the actual survey coverage... do they substantially differ? (maybe chi-squared test idea from Rufener et al.)
  - do with underlying IID and with underlying RW... must be self consistent to pass

## Future worries

- incorporate penalty into likelihood?
- biologically based prior on random walk/AR1 SD?
- how to detect when less obvious changes?
- how do we penalize time just enough?
- to what extent does this happen on much smaller scale
- speeding this up!

## Important points to make

- This isn't *just* a bienniel survey problem - irregular sampling causes it too, joining regions causes it too
  - Should the case studies emphasize this? which region joining one to use??

## Things to make sure to answer

- Does having a ton of data first mean you're OK?? Think IPHC.

## Things others can do

- Norwegian case study from Brian?
- Collect examples of where this happens worldwide for table (which may have to go in supplement with summary in text)
- A region joining case study?

## Case studies

- synoptic stitching - done
- norwegian - done?

## Tables

S1. Examples of where this happens worldwide

## Figures

1. example bad indexes with example biennial sampling illustration?
2. simulation sampling design? (first to drop)
3. subset of example indexes from matrix of options
4. showing getting the spatial pattern wrong and attributing spatial variance to time
5a. Forest plot of correlates
5b. RMSE, Mean SE, coverage forest plot
6. application to real data (as done by Jillian); pick a couple

S1. full matrix of simulation examples (~3 replicates?)
S2. extra examples on real data

