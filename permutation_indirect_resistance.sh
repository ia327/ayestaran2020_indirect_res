#!/bin/bash
while [[ $# -gt 0 ]]
do key="$1"
case $key in
  -h)  printf "############   HELP   ###############\nOPTIONS\n"
  printf "\t-h\tThis Help\n"
  printf "\t-C\tForce the generation of combined_data.RData\n"
  printf "\t-A\tRun tissue-specific ANOVA\n"
  printf "\t-g|--gdsc\tUse GDSC dataset\n"
  printf "\t-c|--ctrp\tUse CTRP dataset\n"
  exit 0
  ;;
  -C)
  force_combined_data=1
  shift
  ;;
  -A)
  run_ANOVA=1
  shift
  ;;
  -g|--gdsc)
  dataset='gdsc'
  shift
  ;;
  -c|--ctrp)
  dataset='ctrp'
  shift
  ;;
esac
done


DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

upperdataset=`echo $dataset | tr '[:lower:]' '[:upper:]'`

if [ "$dataset" != "gdsc" ] && [ "$dataset" != "ctrp" ]
then
  printf "Error: choose a valid dataset: GDSC or CTRP\n"
  exit 0
else
  printf "Analysing data from $upperdataset data set\n"
fi

# Run from main directory

if [ "$force_combined_data" == 1 ]
then
  printf "Creating combined data file...\n"
  Rscript --vanilla "$DIR/scripts/create_combined_data.R" "$DIR"
fi

mkdir -p "$DIR/permutation_test"
mkdir -p "$DIR/permutation_test/$dataset"
mkdir -p "$DIR/permutation_test/$dataset/ANOVA_out/"

for PERM_ID in {1..100}
do
  printf "Permutation test #$PERM_ID\n"

  if [ "$run_ANOVA" == 0 ] && [ ! -d "$DIR/results/ANOVA_out" ]
    then printf "Error: results from ANOVA are missing, try using -A\n"
    exit 0
  fi

  if [ "$run_ANOVA" == 1 ]
    then
    printf "Running tissue specific ANOVA...\n"
    Rscript --vanilla "$DIR/scripts/run_all_anovas.R" "$DIR" "$dataset" \
    "$DIR/permutation_test/$dataset/ANOVA_out/P${PERM_ID}_" "permute"
  fi

  mkdir -p "$DIR/permutation_test/$dataset/ANOVA_filtered/"
  printf "Filtering ANOVA output...\n"
  Rscript --vanilla "$DIR/scripts/anova_out_filter.R" "$DIR" "$dataset" \
  "$DIR/permutation_test/$dataset/ANOVA_out/P${PERM_ID}_" \
  "$DIR/permutation_test/$dataset/ANOVA_filtered/P${PERM_ID}_" "permute"

  mkdir -p "$DIR/permutation_test/$dataset/UNRES_detection/"
  printf "Running variance change analysis...\n"
  Rscript --vanilla "$DIR/scripts/variance_analysis.R" "$DIR" "$dataset" \
  "$DIR/permutation_test/$dataset/ANOVA_filtered/P${PERM_ID}_" \
  "$DIR/permutation_test/$dataset/UNRES_detection/P${PERM_ID}_"

  mkdir -p "$DIR/permutation_test/$dataset/putative_markers/"
  printf "Looking for putative markers within found outliers...\n"
  Rscript --vanilla "$DIR/scripts/seq_cn_lookup.R" "$DIR" "$dataset" \
  "$DIR/permutation_test/$dataset/UNRES_detection/P${PERM_ID}_" \
  "$DIR/permutation_test/$dataset/putative_markers/P${PERM_ID}_" "permute" \
  "$DIR/permutation_test/$dataset/ANOVA_out/P${PERM_ID}_"

done

printf "Done\n"
