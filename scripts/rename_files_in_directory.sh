
# Script by Eleni Nov 22, 2016
#cd /mnt/hgfs/D/sequencing_data/Herring_PopulationStructure/individualData_Trimmed_90bp

#rename renames the named files according to the regular expression perlexpr.
# rename syntax: rename [ -v ] [ -n ] [ -f ] perlexpr [ files ]
#Options

#	-v, --verbose	Verbose: print names of files successfully renamed.
#	-n, --no-act	No Action: show what files would have been renamed.
#	-f, --force	Force: overwrite existing files.
# commands
#	s: substitute one expression for another To substitute one expression for another, the form of perlexpr is: s/expr1/expr2/[gi]
rename "s/oldExtension/newExtension/" *.txt

# \ escapes special character after it, so that it is read literally
# *.txt : do this to all files that end in .txt

rename "s/\.fq\.trimmed/_trimmed\.fq/" *.trimmed ## Substitute the ending .fq.trimmed for :  _trimmed.fq in all files that end with .trimmed in a directory




