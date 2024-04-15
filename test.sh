#!/usr/bin/bash
# QIB Core Bioinformatics 2024

set -euo pipefail

# Check if ./cc.bin exists
if [ ! -f "./cc.bin" ]; then
    echo "Error: ./cc.bin does not exist"
    exit 1
fi

VER=$(./cc.bin -v | cut -f 3 -d " ")

if ! echo "$VER" | perl -ne 'chomp; if ($_=~/(\d+\.\d+)$/) { print "OK\tVersion valid\n" } else { die "Version not passing" } ' ; then
    echo "Error: Output does not end with a version number: $VER"
    false
fi
# Check if test/Matrix.mat.scaled exists
if [ ! -f "test/Matrix.mat.scaled" ]; then
    echo "Error: test/Matrix.mat.scaled does not exist"
    exit 1
fi

# Check if test/MB2.clusters.ext.can.Rhcl exists
if [ ! -f "test/MB2.clusters.ext.can.Rhcl" ]; then
    echo "Error: test/MB2.clusters.ext.can.Rhcl does not exist"
    exit 1
fi

echo -e "OK\tAll required files exist"

# -------- COMMAND 1 --------- #
mkdir -p test/out/
./cc.bin -i test/Matrix.mat.scaled -o test/out/clusters.txt \
		-c test/out/profiles.txt -p MGS --sampleDistMatFile test/out/smpl_dist.mat --sampleDistLog \
		test/out/autocorr.log --sampleMinDist 0.15 --dont_use_mmap -n 44 --progress_stat_file test/out/progress.txt \
		--profile_measure 75Q -b --stop_criteria 100000 --filter_max_top3_sample_contribution 0.7 --max_canopy_dist 0.1 --max_merge_dist 0.1 2> test/out/err.txt > test/out/out.txt

# Check output files
# Check if test/out/smpl_dist.mat exists and is not empty
if [ ! -f "test/out/smpl_dist.mat" ] || [ ! -s "test/out/smpl_dist.mat" ]; then
    echo "Error: test/out/smpl_dist.mat does not exist or is empty"
    exit 1
fi

# Check if test/out/profiles.txt exists and is not empty
if [ ! -f "test/out/profiles.txt" ] || [ ! -s "test/out/profiles.txt" ]; then
    echo "Error: test/out/profiles.txt does not exist or is empty"
    exit 1
fi

# Check if test/out/clusters.txt exists and is not empty
if [ ! -f "test/out/clusters.txt" ] || [ ! -s "test/out/clusters.txt" ]; then
    echo "Error: test/out/clusters.txt does not exist or is empty"
    exit 1
fi

perl -ne 'chomp; next if /^$/; 
@F = split; 
die "Invalid format: $_" unless scalar @F == 2 && $F[0] =~ /^MGS\d+$/ && $F[1] =~ /^\d+(\.\d+)?$/;
END { print "OK\tFile has valid format\n" }' test/out/clusters.txt

# ----------- COMMAND 2 ---------- #

mkdir -p test/out2/
./cc.bin -i test/Matrix.mat.scaled -c test/out/profiles.txt \
	  --referenceMB2 test/MB2.clusters.ext.can.Rhcl \
	  --maxMB2genes 1000 --dont_use_mmap -n 24 -b --stop_criteria 0 \
	  --filter_max_top3_sample_contribution 1 --cag_filter_min_sample_obs 2 \
	  --dont_create_progress_stat_file -o test/out2/MB2.ext.can.pear.test.corr \
	  --max_canopy_dist 0.35 > test/out2/out.txt 2> test/out2/err.txt

perl -ne 'chomp; next if /^$/; @F = split; die "Not a table with 3 columns" unless scalar @F == 3; END { print "OK\tFile is a table with 3 columns\n" }' test/out2/MB2.ext.can.pear.test.corr

echo -e "OK\tTEST TERMINATED"
