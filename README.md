*Update Nov 12, 2013*: You can now escape the evil (Carbon-based) clutches of Net Station and do everything within Matlab/FieldTrip, including but not limited to: filtering, event segmentation, artifact detection/rejection, ICA component rejection, and analysis. All that is required is writing your own trial segmentation function [as described here](http://fieldtrip.fcdonders.nl/example/making_your_own_trialfun_for_conditional_trial_definition). This obviates some of what is described below, including using the ERP PCA toolkit for ICA. Contact me for more information if you are interested, as it may take me a while to update this documentation.

# High-Level Summary

This toolkit contains [MATLAB](http://www.mathworks.com/) functions and scripts that I've written for importing and analyzing EEG data with [FieldTrip](http://www.ru.nl/neuroimaging/fieldtrip/). The main goal is to provide wrapper functions for quickly and easily running analyses and making plots in the way that I've determined has been fruitful for using FieldTrip.

I've set up the functions to import and process [Net Station](http://www.egi.com/) data; they don't currently deal with other formats, but it should relatively easy to modify mat-mvm to import anything that FieldTrip can import. Also included are scripts for using [PyEPL](http://pyepl.sourceforge.net/)-based experimental data in the context of Net Station and FieldTrip, dealing with Net Station files, interfacing with the [ERP PCA Toolkit](http://sourceforge.net/projects/erppcatoolkit/), running statistics, running ICA blink correction, ANOVAs, and more. There are other scripts as well, many of which are specific to my experiments

The toolkit has been tuned for my workflow, and so I can't guarantee how well they will work for others, but I wanted to put them out there so I can keep track of them in a version control system and so that others can use them. I would enjoy hearing feedback that anyone might have.

## Getting Started and Documentation

See the StartingOut page in the [Wiki](http://code.google.com/p/mat-mvm/w/list) tab to get started.

Some functions may not be documented particularly well, and that's because I haven't gotten around to doing that yet. Hopefully the more important functions are documented in the function header, accessible with the 'help' command in MATLAB.

The Functions page has a brief overview of some of the most useful functions to use after you get your data in FieldTrip format.

If you have any questions, please email me. You can find my email address on the project homepage under the members section.

## Licensing

Some functions, like a few in the mat-mvm/stat/rmaov directory, have been downloaded from the MATLAB File Exchange and modified, and they are included because I use them. The original download URL and copyright information is contained in those files, which are covered under the Simplified BSD License. The Simplified BSD License is GPL-compatible, and thus these files can be included here.

## About me

I am a [graduate student](http://psych.colorado.edu/~mollison/) in [Tim Curran's](http://psych.colorado.edu/~tcurran/) [lab](http://psych.colorado.edu/~tclab/) in the [Department of Psychology and Neuroscience](http://psych.colorado.edu/) at the [University of Colorado Boulder](http://www.colorado.edu/).

# How to start using the scripts

*NB*: As noted on the main page, I now do all my processing and analysis within Matlab/FieldTrip, so many of the instructions below are out of date. See mat-mvm/space/space_ft_seg_tla.m and the accompanying mat-mvm/space/space_trialfun.m for examples of the new workflow.

In order to generalize to other experiments, things may need to change, but may not. To get trial metadata during segmentation (which gets put in the trialinfo field by FieldTrip), you will need to either have .evt files exported from your Net Station recordings, events.mat files created from other behavioral data, and/or other ways of reading behavioral data that might need to get built in to seg2ft or ideally would get put in your trialfun function. An example of event file creation can be found at another one of my projects, [https://github.com/warmlogic/expertTrain expertTrain], an experiment framework that started as an expertise training experiment but has grown. See expertTrain/analysis/space_createEvents.m for an example of event creation.

## Introduction

These tools are essentially wrapper functions that I've written for using the [FieldTrip](http://www.ru.nl/neuroimaging/fieldtrip/) toolbox with data recorded in [Net Station](http://www.egi.com/). You'll need to download [FieldTrip](http://www.ru.nl/neuroimaging/fieldtrip/), and possibly the [ERP PCA Toolkit](http://sourceforge.net/projects/erppcatoolkit/) and [EEGLAB](http://sccn.ucsd.edu/eeglab/). For initial behavioral data processing, I also use some functions contained in the Kahana lab's [http://memory.psych.upenn.edu/Software eeg_toolbox]; this will likely only apply to you if you use [PyEPL](http://pyepl.sourceforge.net/) for your experiments.

A few subfolders in mat-mvm (excluding data, eeg, eptoolkit, ft_general, netstation, stat, utilities, etc.) hold experiment-specific scripts. Given the lack of anything resembling a tutorial, it would behoove you to explore (one of) these directories. Read the details below before starting.

## Initial Setup Information

  * The scripts are mostly set up to deal with Net Station (NS) files and may not work well with other formats, though they might.
    * I've never tried anything besides NS and EEGLAB format.
    * For an example of EEGLAB files, look at the `kahn2` scripts, though note that these are set up for a particular file/directory setup that I did not create.
    * If you're using another format, note that you will need to modify parts of some of the main scripts that check on the file extension.
  * For the format of NS files you should use either `EGIS` or `raw` (a.k.a. `sbin`, for simple binary). `EGIS` loads faster.
  * When using NS files, you will need to know the names given to each segment type in the Segmentation Tool. Do not split the segment types into separate files; instead, you should have all of the segment types in a single `EGIS` or `sbin`/`raw` file, with the subject number at the beginning of each filename.
    * Take note that the files need to be stored in a particular path on your computer or the network. This is defined in your analysis script, the important functions mentioned below, and `mm_ft_setSaveDirs.m`.
    * The scripts expect NS files with individual trials, but you can import average files if you'd like. Note that if you use average files, Net Station will do some funky re-naming of your segments. I believe it adds `Sub001` to the beginning of the segment name, but you should check on this. Error checking code in `seg2ft.m` will tell you the available segment names if you have input the wrong ones.
  * I would eventually like to split the scripts into separate data processing and data visualization/analysis scripts, but I haven't gotten around to that yet.
    * I have separated the data processing steps in my scripts for running on a computer cluster; these scripts are in the experiment folders and are (not intuitively) named with `*_ftprocess_*.m`.

## Actually getting your data in MATLAB

To get started in analyzing your data:

  * The most important functions that get your data in FieldTrip format are `create_ft_struct.m`, `seg2ft.m`, and `process_ft_data.m`. They have extensive help sections that you can view with the 'help' command in MATLAB. The help sections may not explain all of the options (yet), so feel free to peruse the top of the code for some of these.
  * You will need to file your such that the tools know where to find it. There are some assumptions made regarding file locations.
  * Take a look at some of the current experiment scripts such as `mat-mvm/soco/soco_ft_seg_tla.m` or `mat-mvm/tnt/tnt_ft_seg_pow.m`, but note that they may not be tuned for your experiment and all options may not be present.
    * I will give a detailed initial setup here in the future. For now, try to glean them from the example experiment scripts.

## Preprocessing

Here's how I preprocessed my data at one point (it has since changed a bit):

  * Net Station: filter, segment/epoch, eye artifact detection, bad channel replacement/interpolation
  * ERP PCA Toolkit: ICA blink removal, baseline correction
  * Net Station: artifact detection, bad channel replace, average rereference
  * FieldTrip: use these tools to read in NS files and get the data in FieldTrip format

## Analyzing and Plotting

See the Functions section.

## Tips

To automatically add mat-mvm to your path, put this in ~/Documents/MATLAB/startup.m:

    myMatlabDir = fullfile(getenv('HOME'),'Documents','MATLAB');
    addpath(genpath(fullfile(myMatlabDir,'mat-mvm')));

To remove version control (.git, .svn, and CVS) directories from your path (I don't know if they would cause any conflicts, but I do this anyway), put this in ~/Documents/MATLAB/startup.m:

    %% remove version control directories from path (.git, .svn, CVS)
    entries = regexp(path, ['[^',pathsep,']*',pathsep], 'match');
    % find the version control entries
    vc_entries = cell2mat(cellfun(@(x) ~isempty(strfind(x,'.git')) | ~isempty(strfind(x,'.svn')) | ~isempty(strfind(x,'CVS')), entries, 'UniformOutput', false));
    % remove them
    rmpath(sprintf(repmat('%s',1,sum(vc_entries)),entries{vc_entries}));


# Brief overview of some of the more important functions

Be sure to read the StartingOut section above before jumping in too deep.

NB: Functions are usually suffixed with ER for ERP data (Event Related) and TFR for time–frequency data (TF Representation). This is adopted from the typical FieldTrip function naming convention.

## Plotting

mm_ft_simpleplotER: An extremely simple function for plotting ERP data

mm_ft_plotER and mm_ft_plotTFR: The most basic grand average plotting functions

mm_ft_subjplotER and mm_ft_subjplotTFR: Plot individual subjects

mm_ft_contrastER and mm_ft_contrastTFR: Plot condition contrasts

mm_ft_clusterplotER and mm_ft_clusterplotTFR: Plot non-parametric statistical clustering results

## Statistics

mm_ft_ttestER and mm_ft_ttestTFR: Run t-tests, comparing two conditions

mm_ft_rmaov2ER, mm_ft_rmaov2TFR, mm_ft_rmaov33ER, and mm_ft_rmaov33TFR: Run ANOVAs comparing multiple conditions, ROIs, and time windows.

mm_ft_clusterstatER and mm_ft_clusterstatTFR: Run non-parametric statistical clustering

mm_ft_corr_dprimeER: Run correlations between ERP voltage and d'. You could probably input other data besides d'.

## Other

There are many other functions, mostly for bookkeeping and getting your data ready to analyze. Check out some of the experiment scripts for examples.

For example, check out send_gmail, which will send you an email (via your gmail account) when processing is done. !ProTip: use Cmd-Enter or Cmd-Shift-Enter to run sections of a script that has been divided up using cell notation (%%) to start your processing, and then follow that up with an email call to let you know when it's done.

[![githalytics.com alpha](https://cruel-carlota.pagodabox.com/6382a7ac8cfc3532d269c944a6155764 "githalytics.com")](http://githalytics.com/warmlogic/mat-mvm)
