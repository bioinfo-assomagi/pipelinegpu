in_dir=$1
python -m bin_jurgen.postprocessing.cli transport stage --in-dir $in_dir
# python -m bin_jurgen.postprocessing.cli transport backup --in-dir $in_dir
python -m bin_jurgen.postprocessing.cli export variants
python -m bin_jurgen.postprocessing.cli export refertazione
python -m bin_jurgen.postprocessing.cli export checklist
