# parama
A tool for sterechemical quality assessment of super high-resolution protein structures
********************************************************************************************************
About PARAMA:
********************************************************************************************************
Please visit http://pauling.mbu.iisc.ac.in/parama/about.html for the science behind this tool. It is also available as a webserver at  http://pauling.mbu.iisc.ac.in/parama .
********************************************************************************************************
PRE-REQUISITES
********************************************************************************************************
PARAMA is entirely written in python and requires atleast python 2.6 installed in your system. Other python packages which are used in the tool and might need to be installed additionally by you are:
  1. numpy
  2. matplotlib
  3. cPickle
  4. PIL
 
INSTRUCTIONS FOR USE:
*********************************************************************************************************
After downloading the software from github, unzip it. 
parama_results.py is the main script, which imports and uses other helper scripts.
To run the program, type
  python parama_results.py <pdb_file_name> <chain_id>
All the results generated from the run will be stored in the results folder. A log file will be generated in the current working directory.

PLEASE NOTE:
  1. The PDB file should be present within the parama folder, with a .pdb extension. Please make sure the file is in the standard PDB format.
  2. Chain ID MUST be provided.
  3. During the run, the PDB file might be rewritten, so please ensure you have a copy of the original PDB file elsewhere.
  
**********************************************************************************************************
