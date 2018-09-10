# Code Functionality
This code segments out fiber tracts from a number of white matter structures including the Arcuate, MdLF, pArc, TPC, Cingulum, and many others. 

### Authors
- Dan Bullock(dnbulloc@indiana.edu)

### Funding 

### References 
Under review

## Running the App 

### On Brainlife.io

Use the compiled and dockerized version of this application found at https://github.com/brain-life/app-wmaSeg

### Running Locally (on your machine)

1.  git clone this repo.
2.  Ensure that you have the requisite dependancies installed (VISTASOFT, MBA, & Freesurfer)
3.  In matlab, pass your whole brain tractogram (path or object, currently known to work with .tck and .fg, trk functionality coming soon) along with the path to the corresponding subject's freesurfer directory (i.e. contains /labels, /surf, /mri, etc) into wma_wrapperDev(wbFG,fsDIR).

### Sample Datasets

If you don't have your own input files, you can download sample datasets from Brainlife.io, or you can use [Brainlife CLI](https://github.com/brain-life/cli).

## Input
-wbFG: either a path to a saved whole brain fiber group / tractogram (fg, fe, or tck; trk soon) or the object itself. 

-fsDIR: path  to THIS SUBJECT'S freesurfer directory (i.e. contains /labels, /surf, /mri, etc).

## Output

-classification:  A structure whose **names** field correspond to the names of the various fiber tracts that were segmented in this function.  The **index** field contains a vector (with as many values as the source whole brain tractogram) that features a number for those streamline indexes which were identified as being part of a segmented tract.  The numeral corresponds to the index of the corresponding name in list of names under the **names** field.

### Dependencies

  - VISTASOFT: https://github.com/vistalab/vistasoft/
  - MBA: https://github.com/francopestilli/mba
  - Freesurfer: https://github.com/freesurfer/freesurfer
