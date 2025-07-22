in_dir=$1
d=$2
in_fq=$3
python -m bin_jurgen.postprocessing.cli transport stage --in-dir $in_dir
# test backup
#python -m bin_jurgen.postprocessing.cli transport backup --destination r --in-dir $in_dir --in-fq /home/magi/ALE_TEST/old/
python -m bin_jurgen.postprocessing.cli transport backup --destination $d --in-dir $in_dir --in-fq $in_fq
python -m bin_jurgen.postprocessing.cli export variants
python -m bin_jurgen.postprocessing.cli export refertazione
python -m bin_jurgen.postprocessing.cli export checklist

python -m bin_jurgen.postprocessing.cli clean stage