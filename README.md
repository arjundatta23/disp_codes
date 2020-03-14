# dispersion_codes
Code accompanying Datta 2017: "A review of array based methods for measurement of multimode surface wave phase velocity dispersion"

Refer to the above for any abbreviations/acronyms used in this README.

This package consists of three main programs, implementing the three methods discussed in Datta 2017:

1. disp_ParkMethod/implement_park.py: implements the frequency-domain slant stack technique.
2. disp_fkMUSIC/implement_MUSIC.py: implements the fk-MUSIC technique.
3. disp_ucd/implement_ucd.py: implements the UC-diagram technique.

Modules which are used by all or some of these programs are in the directory "modules_common".

Apart from these three main programs, this package consists of a few supplementary programs or scripts which perform ancillary functions related to storage or presentation of results. These supplementary programs are described in sepaate README files within the "disp_ucd" and "disp_fkMUSIC" directories.

This parent README discusses only the three main programs.

**********************************************************************************************
I. BASIC CODE USAGE

This package is based on the ObsPy Python framework (https://github.com/obspy/obspy/wiki) and works with SAC-format input seismograms. All three main programs have the same usage format:

python <name of program> <input dir> <file ext>

Here <input dir> is (the path to) a directory containing all the input seismograms (in SAC format) and <file ext> is the extension (along with the leading dot) of the file names for the seismograms. All input seismograms must have the same <file ext>.

Hence for instance, to run the frequency-domain slant-stack code on the provided example directory, do

python disp_ParkMethod/implement_park.py example_2010181072232 .LHZ

Note that the input directory "example_2010181072232" contains other files apart from seismograms, but the code considers only those files in the directory which end, in this case, with ".LHZ".

The fk-MUSIC and UC-diagram codes are run in exactly the same way.

**********************************************************************************************
II. INTERACTIVE USER PROMPTS

All three programs take a number of user inputs before/after running their respective algorithms. Here I first describe those that are common to all three and then those of each program individually.

Common user inputs/prompts:
(not necessarily consecutive, may be separated by other prompts specific to a particular method)

1. "Enter frequency range of interest (Hz): " <f1 f2> f1 and f2 are lower and upper bounds of frequency in Hz. This defines the frequency range in which calculations are performed.
2. "Use all stations (y/n): " This gives the user the option of using all stations/seismograms in <input dir>, or only a subset of them. Enter 'y' if all are to be used. If you enter 'n', you will face the following prompt:
	2a. "Start and end station numbers: " <x1 x2>. x1 and x2 are integers denoting station numbers, when stations are ordered according to epicentral distance. If number of stations is N, "1" denotes the nearest (to source) station in the profile and "N" denotes the farthest. Each of x1 and x2 can be any integer in the range 1-N, with x1 less than x2. The ordering of stations is done internally by the code but the user needs to know the numbers of stations (in the ordered profile) he/she wishes to include in the analysis. This can be known using the "modules_common/seisarray_data.py" module.

NB: this code is meant to be used with data from roughly linear arrays of stations so that an ordering based on epicentral distance leads to a record section type ordering.

**********************************************************************************************
III. USER CONTROLS IN THE SOURCE CODE

windowing etc.


